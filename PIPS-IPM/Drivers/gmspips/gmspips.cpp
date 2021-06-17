#include "PIPSIPMppInterface.hpp"
#include "DistributedInputTree.h"

#include "PIPSIPMppOptions.h"
#include "PreprocessType.h"
#include "MehrotraStrategyType.h"

#include "mpi.h"

#include "pipsdef.h"
#include "gmspips_reader.hpp"
#include "gmspipsio.h"

#include <cstdio>
#include <cstdlib>
#include <cassert>
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
      std::cout << "Usage: " << argv[0]
                << " numBlocks all.gdx|blockstem [GDXLibDir] [scale] [stepLp] [presolve] [printsol] [hierarchical_approach]\n";
      exit(1);
   }

   const int numBlocks = atoi(argv[1]);
   const std::string file_name{argv[2]};
   const std::string path_to_gams{argv[3]};

   for (int i = 5; i <= argc; i++) {
      setParams(scaler_type, primal_dual_step_length, presolve, printsol, hierarchical, argv[i - 1]);
   }

   if (my_rank == 0) {
      std::cout << "reading " << file_name << "\n";
      std::cout << "GAMS located at " << path_to_gams << "\n";
   }

   gmspips_reader reader(file_name, path_to_gams, numBlocks);

   std::unique_ptr<DistributedInputTree> root{reader.read_problem()};

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

   std::vector<double> primalSolVec;
   std::vector<double> dualSolEqVec;
   std::vector<double> dualSolIneqVec;
   std::vector<double> dualSolVarBounds;

   std::vector<double> eqValues;
   std::vector<double> ineqValues;

   pipsipmpp_options::set_bool_parameter("GONDZIO_ADAPTIVE_LINESEARCH", !primal_dual_step_length);
   if (primal_dual_step_length && my_rank == 0) {
      std::cout << "Different steplengths in primal and dual direction are used.\n";
   }

   // create the PIPS-IPM++ interface
   PIPSIPMppInterface pipsIpm(root.get(), primal_dual_step_length ? MehrotraStrategyType::PRIMAL_DUAL : MehrotraStrategyType::PRIMAL, MPI_COMM_WORLD, scaler_type,
         presolve ? PresolverType::PRESOLVER_STOCH : PresolverType::PRESOLVER_NONE);

   if (my_rank == 0) {
      std::cout << "PIPSIPMppInterface created\n";
      std::cout << "solving...\n";
   }

   // run PIPS-IPM++
   pipsIpm.run();
   double objective = pipsIpm.getObjective();

   if (presolve) {
      pipsIpm.postsolveComputedSolution();
   }

   if (printsol) {
      primalSolVec = pipsIpm.gatherPrimalSolution();
      dualSolEqVec = pipsIpm.gatherDualSolutionEq();
      dualSolIneqVec = pipsIpm.gatherDualSolutionIneq();
      dualSolVarBounds = pipsIpm.gatherDualSolutionVarBounds();
      eqValues = pipsIpm.gatherEqualityConsValues();
      ineqValues = pipsIpm.gatherInequalityConsValues();
   }

   if (my_rank == 0)
      std::cout << "solving finished. \n ---Objective value: " << objective << "\n";

   if (printsol && my_rank == 0) {
      const int rc = writeSolution(file_name.c_str(), primalSolVec.size(), dualSolEqVec.size(), dualSolIneqVec.size(), objective, &primalSolVec[0], &dualSolVarBounds[0],
            &eqValues[0], &ineqValues[0], &dualSolEqVec[0], &dualSolIneqVec[0], path_to_gams.c_str());
      if (0 == rc)
         std::cout << "Solution written to " << file_name << "_sol.gdx\n";
      else if (-1 == rc)
         std::cout << "Could not access " << file_name << ".map\n";
      else
         std::cout << "Other error writing solution: rc=" << rc << "\n";
   }

   MPI_Barrier(MPI_COMM_WORLD);
   const double t1 = MPI_Wtime();

   if (my_rank == 0)
      std::cout << "---total time (in sec.): " << t1 - t0 << "\n";

   MPI_Finalize();

   return 0;
}
