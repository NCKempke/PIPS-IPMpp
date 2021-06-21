/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef PIPSIPMPPINTERFACE_H
#define PIPSIPMPPINTERFACE_H

#include <vector>
#include <string>
#include <memory>

#include "TerminationStatus.hpp"
#include "MehrotraStrategyType.h"
#include "pipsdef.h"
#include "PreprocessType.h"

template<typename T>
class Vector;

class DistributedInputTree;
class DistributedFactory;
class PreprocessFactory;
class Residuals;
class Variables;
class DistributedProblem;
class Problem;
class Presolver;
class Postsolver;
class Scaler;
class InteriorPointMethod;

class PIPSIPMppInterface {
public:
   PIPSIPMppInterface(DistributedInputTree* tree, MehrotraStrategyType mehrotra_heuristic, MPI_Comm = MPI_COMM_WORLD, ScalerType scaler_type = ScalerType::NONE,
         PresolverType presolver_type = PresolverType::NONE, const std::string& settings = "PIPSIPMpp.opt");

   ~PIPSIPMppInterface();

   TerminationStatus run();

   double getObjective();

   [[nodiscard]] int n_iterations() const;

   [[nodiscard]] double getFirstStageObjective() const;

   std::vector<double> gatherPrimalSolution();

   std::vector<double> gatherDualSolutionEq();

   std::vector<double> gatherDualSolutionIneq();

   std::vector<double> gatherDualSolutionIneqUpp();

   std::vector<double> gatherDualSolutionIneqLow();

   std::vector<double> gatherDualSolutionVarBounds();

   std::vector<double> gatherDualSolutionVarBoundsUpp();

   std::vector<double> gatherDualSolutionVarBoundsLow();

   std::vector<double> gatherSlacksInequalityUp();

   std::vector<double> gatherSlacksInequalityLow();

   std::vector<double> gatherSlacksVarsUp();

   std::vector<double> gatherSlacksVarsLow();

   std::vector<double> gatherPrimalResidsEQ();

   std::vector<double> gatherPrimalResidsIneqUp();

   std::vector<double> gatherPrimalResidsIneqLow();

   std::vector<double> gatherDualResids();

   [[nodiscard]] std::vector<double> getFirstStagePrimalColSolution() const;

   [[nodiscard]] std::vector<double> getSecondStagePrimalColSolution(int scen) const;

   void allgatherBlocksizes(std::vector<unsigned int>& block_lengths_col,
      std::vector<unsigned int>& block_lengths_A, std::vector<unsigned int>& block_lengths_C) const;

   void postsolveComputedSolution();

   std::vector<double> gatherEqualityConsValues();

   std::vector<double> gatherInequalityConsValues();

   void getVarsUnscaledUnperm();

   void getResidsUnscaledUnperm();
   //more get methods to follow here

private:
   static void printComplementarityResiduals(const Variables& vars) ;

   std::vector<double> gatherFromSolution(std::unique_ptr<Vector<double>> Variables::* member_to_gather);
   std::vector<double> gatherFromResiduals(std::unique_ptr<Vector<double>> Residuals::* member_to_gather);

protected:
   std::unique_ptr<DistributedFactory> factory;
   std::unique_ptr<PreprocessFactory> preprocess_factory;

   std::unique_ptr<Problem> presolved_problem; // possibly presolved problem
   std::unique_ptr<Problem> dataUnpermNotHier; // data after presolve before permutation, scaling and hierarchical data
   std::unique_ptr<Problem> original_problem; // original data
   std::unique_ptr<Variables> variables;
   std::unique_ptr<Variables> unscaleUnpermNotHierVars;
   std::unique_ptr<Variables> postsolved_variables;

   std::unique_ptr<Residuals> residuals;
   std::unique_ptr<Residuals> unscaleUnpermNotHierResids;
   std::unique_ptr<Residuals> postsolvedResids;

   std::unique_ptr<Presolver> presolver;
   std::unique_ptr<Postsolver> postsolver;
   std::unique_ptr<Scaler> scaler;
   std::unique_ptr<InteriorPointMethod> solver;

   MPI_Comm comm = MPI_COMM_NULL;
   const int my_rank = -1;
   TerminationStatus result{TerminationStatus::DID_NOT_RUN};
   bool ran_solver = false;
};

#endif // PIPSIPMPPINTERFACE_H
