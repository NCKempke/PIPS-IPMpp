#include "PIPSIPMppInterface.hpp"
#include "DistributedInputTree.h"

#include "PIPSIPMppOptions.h"
#include "PreprocessType.h"
#include "InteriorPointMethodType.hpp"

#include "mpi.h"

#include "pipsdef.h"
#include "gmspips_reader.hpp"
#include "gmspipsio.h"

#include <cstdlib>
#include <cstring>
#include <iostream>

static void setParams(ScalerType& scaler_type, bool& stepDiffLp, bool& presolve, bool& printsol, bool& hierarchical, const char* paramname) {
   if (strcmp(paramname, "scale") == 0 || strcmp(paramname, "scaleEqui") == 0)
      scaler_type = ScalerType::SCALER_EQUI_STOCH;
   else if (strcmp(paramname, "scaleGeo") == 0)
      scaler_type = ScalerType::SCALER_GEO_STOCH;
   else if (strcmp(paramname, "scaleGeoEqui") == 0)
      scaler_type = ScalerType::SCALER_GEO_EQUI_STOCH;
   else if (strcmp(paramname, "stepLp") == 0)
      stepDiffLp = true;
   else if (strcmp(paramname, "presolve") == 0)
      presolve = true;
   else if (strcmp(paramname, "printsol") == 0)
      printsol = true;
   else if (strcmp(paramname, "hierarchical") == 0)
      hierarchical = true;
}

void writeSolutionGDX(PIPSIPMppInterface& ipm_interface, const std::string& file_name, const std::string& path_to_gams) {
   std::vector<double> primalSolVec;
   std::vector<double> dualSolEqVec;
   std::vector<double> dualSolIneqVec;
   std::vector<double> dualSolVarBounds;

   std::vector<double> eqValues;
   std::vector<double> ineqValues;

   const double objective = ipm_interface.getObjective();
   primalSolVec = ipm_interface.gatherPrimalSolution();
   dualSolEqVec = ipm_interface.gatherDualSolutionEq();
   dualSolIneqVec = ipm_interface.gatherDualSolutionIneq();
   dualSolVarBounds = ipm_interface.gatherDualSolutionVarBounds();
   eqValues = ipm_interface.gatherEqualityConsValues();
   ineqValues = ipm_interface.gatherInequalityConsValues();

   if (PIPS_MPIgetRank() == 0) {
      const int rc = writeSolution(file_name.c_str(), primalSolVec.size(), dualSolEqVec.size(), dualSolIneqVec.size(), objective, &primalSolVec[0], &dualSolVarBounds[0],
         &eqValues[0], &ineqValues[0], &dualSolEqVec[0], &dualSolIneqVec[0], path_to_gams.c_str());
      if (0 == rc)
         std::cout << "Solution written to " << file_name << "_sol.gdx\n";
      else if (-1 == rc)
         std::cout << "Could not access " << file_name << ".map\n";
      else
         std::cout << "Other error writing solution: rc=" << rc << "\n";
   }
}

int main(int argc, char** argv) {

   MPI_Init(&argc, &argv);
   const int my_rank = PIPS_MPIgetRank();

   MPI_Barrier(MPI_COMM_WORLD);
   const double t0 = MPI_Wtime();

   ScalerType scaler_type = ScalerType::SCALER_NONE;
   bool primal_dual_step_length = false;
   bool presolve = false;
   bool printsol = false;
   bool hierarchical = false;

   if ((argc < 3) || (argc > 9)) {
      std::cout << "Usage: " << argv[0] << " numBlocks file_name[.gdx] GDXLibDir [scale] [stepLp] [presolve] [printsol] [hierarchical_approach]\n";
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

   gmspips_reader reader(file_name, path_to_gams, numBlocks);

   std::unique_ptr<DistributedInputTree> root{reader.read_problem()};

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

   PIPSIPMppInterface pipsIpm(root.get(), primal_dual_step_length ? InteriorPointMethodType::PRIMAL_DUAL : InteriorPointMethodType::PRIMAL, MPI_COMM_WORLD, scaler_type,
         presolve ? PresolverType::PRESOLVER_STOCH : PresolverType::PRESOLVER_NONE);

   if (my_rank == 0) {
      std::cout << "PIPSIPMppInterface created\n";
      std::cout << "solving...\n";
   }

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
      writeSolutionGDX(pipsIpm, file_name, path_to_gams);
   }

   MPI_Barrier(MPI_COMM_WORLD);
   const double tn = MPI_Wtime();

   if (my_rank == 0)
      std::cout << "---total time (in sec.): " << tn - t0 << "\n";

   MPI_Finalize();

   return 0;
}
