#include "../../../Core/Interface/PIPSIPMppInterface.hpp"
#include "../../../Core/Options/PIPSIPMppOptions.h"
#include "DistributedInputTree.h"
#include "gmspips_reader.hpp"

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

int main(int argc, char** argv) {

   MPI_Init(&argc, &argv);
   const int my_rank = PIPS_MPIgetRank();

   MPI_Barrier(MPI_COMM_WORLD);
   const double t0 = MPI_Wtime();

   ScalerType scaler_type = ScalerType::NONE;
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
         presolve ? PresolverType::PRESOLVE : PresolverType::NONE);

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
      reader.write_solution(pipsIpm, file_name);
   }

   MPI_Barrier(MPI_COMM_WORLD);
   const double tn = MPI_Wtime();

   if (my_rank == 0)
      std::cout << "---total time (in sec.): " << tn - t0 << "\n";

   MPI_Finalize();

   return 0;
}
