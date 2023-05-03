#include "PIPSIPMppInterface.hpp"
#include "PIPSIPMppOptions.h"
#include "Residuals.h"
#include "DistributedInputTree.h"
#include "DistributedFactory.hpp"
#include "gmspips_reader.hpp"
#include "Problem.hpp"

#include "mpi.h"

#include <cstdlib>
#include <cstring>
#include <iostream>

static void setParams(ScalerType& scaler_type, bool& stepDiffLp, bool& presolve, bool& printsol, bool& hierarchical, const char* paramname) {
   if (strcmp(paramname, "scale") == 0 || strcmp(paramname, "scaleEqui") == 0)
      scaler_type = ScalerType::EQUILIBRIUM;
   else if (strcmp(paramname, "scaleGeo") == 0)
      scaler_type = ScalerType::GEOMETRIC_MEAN;
   else if (strcmp(paramname, "scaleGeoEqui") == 0)
      scaler_type = ScalerType::GEOMETRIC_MEAN_EQUILIBRIUM;
   else if (strcmp(paramname, "scaleCurtisReid") == 0)
      scaler_type = ScalerType::CURTIS_REID;
   else if (strcmp(paramname, "stepLp") == 0)
      stepDiffLp = true;
   else if (strcmp(paramname, "presolve") == 0)
      presolve = true;
   else if (strcmp(paramname, "printsol") == 0)
      printsol = true;
   else if (strcmp(paramname, "hierarchical") == 0)
      hierarchical = true;
}

struct Result {
    TerminationStatus status;
    double objective;
    double constraint_violation;
};

Result run_pips(PIPSIPMppInterface& pipsIpm, const Problem& problem) {
    // set the problem
    pipsIpm.get_presolved_problem() = problem;

    // solve the continuous problem
    TerminationStatus status = TerminationStatus::DID_NOT_RUN;
    double objective = std::numeric_limits<double>::infinity();
    double constraint_violation = std::numeric_limits<double>::infinity();
    int number_iterations = -1;
    try {
        status = pipsIpm.run();
        objective = pipsIpm.getObjective();
        number_iterations = pipsIpm.n_iterations();
        constraint_violation = pipsIpm.get_residuals().constraint_violation();
    }
    catch (const std::exception&) {
    }
    if (PIPS_MPIgetRank() == 0) {
        std::cout << "Status: " << status << " reached in " << number_iterations << " iterations.\n";
    }
    return {status, objective, constraint_violation};
}

int main(int argc, char** argv) {

   MPI_Init(&argc, &argv);
   const int my_rank = PIPS_MPIgetRank();

   MPI_Barrier(MPI_COMM_WORLD);
   const double t0 = MPI_Wtime();

   ScalerType scaler_type = ScalerType::NONE;
   bool primal_dual_step_length = true;
   bool presolve = false;
   bool printsol = false;
   bool hierarchical = false;

   if ((argc < 3) || (argc > 9)) {
      std::cout << "Usage: " << argv[0] << " numBlocks file_name[.gdx] GDXLibDir [scale] [stepLp] [presolve] [printsol] [hierarchical]\n";
      std::cout << "Expecting files to be of the form \"file_nameXX\", where XX specifies a block from 0 to num_blocks!\n";
      exit(1);
   }

   const size_t numBlocks = atoi(argv[1]);
   const std::string file_name{argv[2]};
   const std::string path_to_gams{argv[3]};

   for (int i = 5; i <= argc; ++i) {
      setParams(scaler_type, primal_dual_step_length, presolve, printsol, hierarchical, argv[i - 1]);
   }

   if (my_rank == 0) {
      std::cout << "reading " << file_name << "\n";
      std::cout << "GAMS located at " << path_to_gams << "\n";
   }

   gmspips_reader gams_reader(file_name, path_to_gams, numBlocks);
   std::unique_ptr<DistributedInputTree> tree{gams_reader.read_problem()};

   MPI_Barrier(MPI_COMM_WORLD);
   if (my_rank == 0)
      std::cout << "Using a total of " << PIPS_MPIgetSize() << " MPI processes.\n";

   if (hierarchical) {
      if (my_rank == 0)
         std::cout << "Using Hierarchical approach\n";
      pipsipmpp_options::activate_hierarchial_approach();
   }

   pipsipmpp_options::set_int_parameter("OUTER_SOLVE", 2);
   if (my_rank == 0)
      std::cout << "Using outer BICGSTAB\n";

   if (my_rank == 0 && pipsipmpp_options::get_int_parameter("INNER_SC_SOLVE") == 2)
      std::cout << "Using inner BICGSTAB\n";

   pipsipmpp_options::set_bool_parameter("GONDZIO_ADAPTIVE_LINESEARCH", !primal_dual_step_length);
   if (primal_dual_step_length && my_rank == 0) {
      std::cout << "Different steplengths in primal and dual direction are used.\n";
   }

   if( my_rank == 0 )
      std::cout << "Creating PIPSIpmInterface ...\n";

   PIPSIPMppInterface pipsIpm(tree.get(), primal_dual_step_length ? InteriorPointMethodType::PRIMAL_DUAL : InteriorPointMethodType::PRIMAL, MPI_COMM_WORLD, scaler_type,
         presolve ? PresolverType::PRESOLVE : PresolverType::NONE);

   if (my_rank == 0) {
      std::cout << "PIPSIPMppInterface created\n";
   }

   // relax the integers

    // get information from the problem
    //const Vector<double>& initial_point = pipsIpm.get_primal_variables();
    const Problem& problem = pipsIpm.get_presolved_problem();
    //assert (problem.variable_integrality_type != nullptr && "The problem has only continuous variables");
    if (my_rank == 0) {
        if (problem.variable_integrality_type == nullptr) {
            std::cout << "All variables are continuous\n";
        }
        else {
            std::cout << "The problem has discrete variables\n";
        }
    }

    const Vector<double>& integrality = *problem.variable_integrality_type;
    std::unique_ptr<DistributedFactory> factory = std::make_unique<DistributedFactory>(tree.get(), MPI_COMM_WORLD);
    const size_t number_variables = factory->tree->nx();
    if (my_rank == 0) {
        std::cout << "The problem has " << number_variables << " variables\n";
        std::cout << "solving...\n";
    }
    // lower bounding: solve the LP relaxation
    // if number_blocks = 5: continuous relaxation: Objective: 29.0771 reached in 24 iterations
    const Result LP_result = run_pips(pipsIpm, problem);
    if (my_rank == 0) {
        std::cout << "\nComputing the lower bound (LP relaxation):\n";
        std::cout << "LP status: " << LP_result.status << "\n";
        std::cout << "LP objective: " << LP_result.objective << "\n";
        std::cout << "LP constraint_violation: " << LP_result.constraint_violation << "\n";
    }

    // generate an initial integer point
    pipsIpm.get_primal_variables().transform_value([&](double current_value, double /*lb*/, double /*ub*/, double variable_integrality) {
        if (0 < variable_integrality) {
            return std::round(current_value);
        }
        else {
            return current_value;
        }
    }, *problem.primal_lower_bounds, *problem.primal_upper_bounds, integrality);

   const TerminationStatus status = pipsIpm.run();
   const double objective = pipsIpm.getObjective();

   const double t1 = MPI_Wtime();

   if (my_rank == 0){
      if (status != TerminationStatus::SUCCESSFUL_TERMINATION) {
         std::cout << "Failed to solve Instance successfully.\n";
      } else {
         std::cout << "Solving finished.\n";
      }

      std::cout << "---Objective value: " << objective << "\n";
      std::cout << "solving took " << t1 - t0 << " seconds\n";
   }

   if (printsol) {
      gams_reader.write_solution(pipsIpm, file_name);
   }

   MPI_Barrier(MPI_COMM_WORLD);
   const double tn = MPI_Wtime();

   if (my_rank == 0)
      std::cout << "---total time (in sec.): " << tn - t0 << "\n";

   MPI_Finalize();

   return 0;
}
