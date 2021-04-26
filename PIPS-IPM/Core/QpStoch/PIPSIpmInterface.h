/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef PIPSIPM_INTERFACE
#define PIPSIPM_INTERFACE

#include <algorithm>
#include <functional>
#include <cstdlib>
#include <stdexcept>
#include <memory>

//#include "stochasticInput.hpp"
//#include "sTreeImpl.h"

#include "sTree.h"
#include "DistributedQP.hpp"
#include "DistributedResiduals.hpp"
#include "sVars.h"
#include "StochMonitor.h"


#include "PreprocessFactory.h"
#include "Scaler.h"
#include "Presolver.h"
#include "Postsolver.h"

#include "sTreeCallbacks.h"
#include "pipsport.h"

#include "StochOptions.h"

//#define PRESOLVE_POSTSOLVE_ONLY // will not call solve routine an just presolve and then postsolve the problem - for debugging presolve and postsolve operations

template<class FORMULATION, class IPMSOLVER>
class PIPSIpmInterface {
public:
   PIPSIpmInterface(StochInputTree* in, MPI_Comm = MPI_COMM_WORLD, ScalerType scaler_type = SCALER_NONE,
         PresolverType presolver_type = PRESOLVER_NONE, std::string settings = "PIPSIPMpp.opt");

   ~PIPSIpmInterface() = default;

   void run();

   double getObjective();

   double getFirstStageObjective() const;

   void setPrimalTolerance(double val);

   void setDualTolerance(double val);

   std::vector<double> gatherPrimalSolution();

   std::vector<double> gatherDualSolutionEq();

   std::vector<double> gatherDualSolutionIneq();

   std::vector<double> gatherDualSolutionIneqUpp();

   std::vector<double> gatherDualSolutionIneqLow();

   std::vector<double> gatherDualSolutionVarBounds();

   std::vector<double> gatherDualSolutionVarBoundsUpp();

   std::vector<double> gatherDualSolutionVarBoundsLow();

   std::vector<double> getFirstStagePrimalColSolution() const;

   std::vector<double> getSecondStagePrimalColSolution(int scen) const;

   std::vector<double> getFirstStageDualRowSolution() const;

   std::vector<double> getSecondStageDualRowSolution(int scen) const;

   void postsolveComputedSolution();

private:
   void printComplementarityResiduals(const sVars& vars) const;

   std::vector<double> gatherFromSolution(OoqpVectorHandle sVars::* member_to_gather);

public:
   std::vector<double> gatherEqualityConsValues();

   std::vector<double> gatherInequalityConsValues();

   void getVarsUnscaledUnperm();

   void getResidsUnscaledUnperm();
   //more get methods to follow here

   static bool isDistributed() { return true; }

protected:
   std::unique_ptr<FORMULATION> formulation_factory{};
   std::unique_ptr<PreprocessFactory> preprocess_factory{};

   std::unique_ptr<DistributedQP> presolved_problem{};       // possibly presolved problem
   std::unique_ptr<DistributedQP> dataUnpermNotHier{}; // data after presolve before permutation, scaling and hierarchical data
   std::unique_ptr<DistributedQP> original_problem{};   // original data
   std::unique_ptr<sVars> vars{};
   std::unique_ptr<sVars> unscaleUnpermNotHierVars{};
   std::unique_ptr<sVars> postsolvedVars{};

   std::unique_ptr<DistributedResiduals> residuals{};
   std::unique_ptr<DistributedResiduals> unscaleUnpermNotHierResids{};
   std::unique_ptr<DistributedResiduals> postsolvedResids{};

   std::unique_ptr<Presolver> presolver{};
   std::unique_ptr<Postsolver> postsolver{};
   std::unique_ptr<Scaler> scaler{};
   std::unique_ptr<IPMSOLVER> solver{};

   MPI_Comm comm = MPI_COMM_NULL;
   const int my_rank = -1;
   int result{-1};
   bool ran_solver = false;
};

//----------------------------------------------------------------------
// IMPLEMENTATION
//----------------------------------------------------------------------

template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::PIPSIpmInterface(StochInputTree* in, MPI_Comm comm, ScalerType scaler_type, PresolverType presolver_type,
      std::string settings) : comm(comm), my_rank(PIPS_MPIgetRank()) {
   pips_options::setOptions(settings);
   const bool postsolve = pips_options::getBoolParameter("POSTSOLVE");

   MPI_Barrier(comm);
   const double t0 = MPI_Wtime();

   formulation_factory.reset(new FORMULATION(in, comm));
#ifdef TIMING
   if( my_rank == 0 ) printf("factory created\n");
#endif

   preprocess_factory.reset(new PreprocessFactory());
#ifdef TIMING
   if( my_rank == 0 ) printf("prefactory created\n");
#endif

   // presolving activated?
   if (presolver_type != PRESOLVER_NONE) {
      original_problem.reset(dynamic_cast<DistributedQP*>(formulation_factory->create_problem()));

      MPI_Barrier(comm);
      const double t0_presolve = MPI_Wtime();

      if (postsolve)
         postsolver.reset(preprocess_factory->makePostsolver(original_problem.get()));

      presolver.reset(
            preprocess_factory->makePresolver(dynamic_cast<sFactory*>(formulation_factory.get())->tree, original_problem.get(), presolver_type,
                  postsolver.get()));

      presolved_problem.reset(dynamic_cast<DistributedQP*>(presolver->presolve()));

      formulation_factory->data = presolved_problem.get(); // todo update also sTree* of factory

      MPI_Barrier(comm);
      const double t_presolve = MPI_Wtime();
      if (my_rank == 0)
         std::cout << "---presolve time (in sec.): " << t_presolve - t0_presolve << "\n";
   }
   else {
      presolved_problem.reset(dynamic_cast<DistributedQP*>(formulation_factory->create_problem()));
      assert(presolved_problem);
   }

//  data->writeToStreamDense(std::cout);

#if 0
   ofstream myfile;
   myfile.open ("PipsToMPS_prslv.mps");
   data->writeMPSformat(myfile);
   myfile.close();
#endif

#ifdef TIMING
   if( my_rank == 0 ) printf("data created\n");
#endif

   dataUnpermNotHier.reset(presolved_problem->cloneFull());

   // after identifying the linking structure switch to hierarchical data structure -> will this do anything to the scaler?
   if (pips_options::getBoolParameter("PARDISO_FOR_GLOBAL_SC"))
      presolved_problem->activateLinkStructureExploitation();

   // TODO : save "old" data somewhere?
   if (pips_options::getBoolParameter("HIERARCHICAL")) {
      if (my_rank == 0)
         std::cout << "Using hierarchical approach!\n";

      presolved_problem.reset(dynamic_cast<DistributedQP*>(formulation_factory->switchToHierarchicalData(presolved_problem.release())));

      if (pips_options::getBoolParameter("HIERARCHICAL_PRINT_HIER_DATA"))
         presolved_problem->writeToStreamDense(std::cout);
   }

   vars.reset(dynamic_cast<sVars*>( formulation_factory->makeVariables(presolved_problem.get())));
#ifdef TIMING
   if( my_rank == 0 ) printf("variables created\n");
#endif

   residuals.reset(dynamic_cast<DistributedResiduals*>( formulation_factory->makeResiduals(presolved_problem.get())));
#ifdef TIMING
   if( my_rank == 0 ) printf("resids created\n");
#endif

   scaler.reset(preprocess_factory->makeScaler(presolved_problem.get(), scaler_type));

#ifdef TIMING
   if( my_rank == 0 ) printf("scaler created\n");
#endif

   if (scaler) {
      MPI_Barrier(comm);
      const double t0_scaling = MPI_Wtime();

      scaler->scale();

      MPI_Barrier(comm);
      const double t_scaling = MPI_Wtime();
      if (my_rank == 0)
         std::cout << "---scaling time (in sec.): " << t_scaling - t0_scaling << "\n";
   }

   solver.reset(new IPMSOLVER(*formulation_factory, *presolved_problem, scaler.get()));
   solver->addMonitor(new StochMonitor(formulation_factory.get(), scaler.get()));
#ifdef TIMING
   if( my_rank == 0 ) printf("solver created\n");
   //solver->monitorSelf();
#endif

   MPI_Barrier(comm);
   const double t1 = MPI_Wtime();
   if (my_rank == 0)
      std::cout << "---reading time (in sec.): " << t1 - t0 << "\n";
}


template<typename FORMULATION, typename IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::run() {
   if (my_rank == 0)
      std::cout << "solving ...\n";

   if (my_rank == 0) {
      // TODO : use unlifted data....
      if (!pips_options::getBoolParameter("HIERARCHICAL")) {
         std::cout << "1st stage " << presolved_problem->getLocalnx() << " variables, " << presolved_problem->getLocalmy()
                   << " equality constraints, " << presolved_problem->getLocalmz() << " inequality constraints.\n";

         const int nscens = presolved_problem->children.size();
         if (nscens) {
            std::cout << "2nd stage " << presolved_problem->children[0]->getLocalnx() << " variables, "
                      << presolved_problem->children[0]->getLocalmy() << " equality constraints, " << presolved_problem->children[0]->getLocalmz()
                      << " inequality constraints.\n";

            std::cout << nscens << " scenarios." << "\n";
            std::cout << "Total " << presolved_problem->getLocalnx() + nscens * presolved_problem->children[0]->getLocalnx() << " variables, "
                      << presolved_problem->getLocalmy() + nscens * presolved_problem->children[0]->getLocalmy() << " equality constraints, "
                      << presolved_problem->getLocalmz() + nscens * presolved_problem->children[0]->getLocalmz() << " inequality constraints.\n";
         }
      }
   }
#ifdef TIMING
   double tmElapsed=MPI_Wtime();
#endif

#if defined(PRESOLVE_POSTSOLVE_ONLY) && !defined(NDEBUG)
   const int result = 0;
#else
   //---------------------------------------------
   //result = solver->solve(presolved_problem.get(), vars.get(), resids.get());
   result = solver->solve(*presolved_problem, *vars, *residuals);
   //---------------------------------------------
#endif

   if (result != 0 && my_rank == 0)
      std::cout << "failed to solve instance, result code: " << result << "\n";

   ran_solver = true;

#ifdef TIMING
   if ( 0 != result )
      return;

   tmElapsed = MPI_Wtime()-tmElapsed;

   const double objective = getObjective();

   if( my_rank == 0 ) {
    //std::cout << " " << data->nx << " variables, " << data->my
    // << " equality constraints, " << data->mz << " inequality constraints.\n";

    std::cout << " Iterates: " << solver->iter <<",    Optimal Solution:  "
          << objective << "\n";

    std::cout << "Solve time: " << tmElapsed << " seconds." << endl;

    char *var = getenv("OMP_NUM_THREADS");
    if(var != nullptr) {
      int num_threads;
      sscanf( var, "%d", &num_threads );
      std::cout << "Num threads: " << num_threads << "\n";
    }
  }
#endif

#if !defined(NDEBUG) && defined(PRESOLVE_POSTSOLVE_ONLY)
   postsolveComputedSolution();
#endif
}

template<typename FORMULATION, typename SOLVER>
double PIPSIpmInterface<FORMULATION, SOLVER>::getObjective() {

   if (!ran_solver)
      throw std::logic_error("Must call go() and start solution process before trying to retrieve original solution");

   if (postsolver != nullptr && postsolvedVars == nullptr)
      this->postsolveComputedSolution();

   double obj;
   if (postsolvedVars != nullptr)
      obj = original_problem->objective_value(postsolvedVars.get());
   else {
      obj = presolved_problem->objective_value(vars.get());
      if (scaler)
         obj = scaler->getObjUnscaled(obj);
   }

   return obj;
}


template<typename FORMULATION, typename SOLVER>
double PIPSIpmInterface<FORMULATION, SOLVER>::getFirstStageObjective() const {
   OoqpVector& x = *(dynamic_cast<StochVector&>(*vars->x).first);
   OoqpVector& c = *(dynamic_cast<StochVector&>(*presolved_problem->g).first);
   return c.dotProductWith(x);
}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::getVarsUnscaledUnperm() {
   assert(unscaleUnpermNotHierVars == nullptr);
   assert(dataUnpermNotHier);

   if (!ran_solver)
      throw std::logic_error("Must call go() and start solution process before trying to retrieve unscaled unpermuted solution");
   if (scaler) {
      std::unique_ptr<sVars> unscaled_vars{dynamic_cast<sVars*>(scaler->getVariablesUnscaled(*vars))};
      unscaleUnpermNotHierVars.reset(presolved_problem->getVarsUnperm(*unscaled_vars, *dataUnpermNotHier));
   }
   else
      unscaleUnpermNotHierVars.reset(presolved_problem->getVarsUnperm(*vars, *dataUnpermNotHier));

}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::getResidsUnscaledUnperm() {
   assert(unscaleUnpermNotHierResids == nullptr);
   assert(dataUnpermNotHier);

   if (!ran_solver)
      throw std::logic_error("Must call go() and start solution process before trying to retrieve unscaled unpermuted residuals");
   if (scaler) {
      std::unique_ptr<DistributedResiduals> unscaled_resids{dynamic_cast<DistributedResiduals*>(scaler->getResidualsUnscaled(*residuals))};
      unscaleUnpermNotHierResids.reset(presolved_problem->getResidsUnperm(*unscaled_resids, *dataUnpermNotHier));
   }
   else
      unscaleUnpermNotHierResids.reset(presolved_problem->getResidsUnperm(*residuals, *dataUnpermNotHier));
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherFromSolution(OoqpVectorHandle sVars::* member_to_gather) {
   if (unscaleUnpermNotHierVars == nullptr)
      this->getVarsUnscaledUnperm();

   if (postsolver != nullptr && postsolvedVars == nullptr)
      this->postsolveComputedSolution();

   std::vector<double> vec;
   if (postsolver == nullptr)
      vec = dynamic_cast<const StochVector&>(*(unscaleUnpermNotHierVars.get()->*member_to_gather)).gatherStochVector();
   else
      vec = dynamic_cast<const StochVector&>(*(postsolvedVars.get()->*member_to_gather)).gatherStochVector();

   return vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherPrimalSolution() {
   return gatherFromSolution(&sVars::x);
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionEq() {
   return gatherFromSolution(&sVars::y);
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneq() {
   return gatherFromSolution(&sVars::z);
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneqUpp() {
   return gatherFromSolution(&sVars::pi);
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneqLow() {
   return gatherFromSolution(&sVars::lambda);
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBounds() {
   std::vector<double> duals_varbounds_upp = gatherDualSolutionVarBoundsUpp();
   std::vector<double> duals_varbounds_low = gatherDualSolutionVarBoundsLow();

   assert(duals_varbounds_low.size() == duals_varbounds_upp.size());

   std::vector<double> duals_varbounds;
   duals_varbounds.reserve(duals_varbounds_low.size());

   std::transform(duals_varbounds_low.begin(), duals_varbounds_low.end(), duals_varbounds_upp.begin(), std::back_inserter(duals_varbounds),
         std::minus<double>());

   return duals_varbounds;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBoundsUpp() {
   return gatherFromSolution(&sVars::phi);
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBoundsLow() {
   return gatherFromSolution(&sVars::gamma);
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherEqualityConsValues() {
   if (unscaleUnpermNotHierResids == nullptr)
      this->getResidsUnscaledUnperm();

   if (postsolver != nullptr && postsolvedVars == nullptr)
      this->postsolveComputedSolution();

   StochVector* eq_vals = (postsolvedVars == nullptr) ? dynamic_cast<StochVector*>(unscaleUnpermNotHierResids->rA->cloneFull())
                                                      : dynamic_cast<StochVector*>(postsolvedResids->rA->cloneFull());

   if (original_problem == nullptr || postsolvedVars == nullptr)
      eq_vals->axpy(1.0, *presolved_problem->bA);
   else
      eq_vals->axpy(1.0, *original_problem->bA);

   std::vector<double> eq_vals_vec = eq_vals->gatherStochVector();

   delete eq_vals;

   return eq_vals_vec;
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherInequalityConsValues() {
   if (unscaleUnpermNotHierVars == nullptr)
      this->getVarsUnscaledUnperm();

   if (unscaleUnpermNotHierResids == nullptr)
      this->getResidsUnscaledUnperm();

   if (postsolver != nullptr && postsolvedVars == nullptr)
      this->postsolveComputedSolution();

   StochVector* ineq_vals = (postsolvedVars == nullptr) ? dynamic_cast<StochVector*>(unscaleUnpermNotHierResids->rC->cloneFull())
                                                        : dynamic_cast<StochVector*>(postsolvedResids->rC->cloneFull());

   if (postsolvedVars == nullptr)
      ineq_vals->axpy(1.0, *unscaleUnpermNotHierVars->s);
   else
      ineq_vals->axpy(1.0, *postsolvedVars->s);

   std::vector<double> ineq_vals_vec = ineq_vals->gatherStochVector();

   delete ineq_vals;

   return ineq_vals_vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getFirstStagePrimalColSolution() const {
   SimpleVector const& v = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->x).first);
   return std::vector<double>(&v[0], &v[0] + v.length());
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getSecondStagePrimalColSolution(int scen) const {
   SimpleVector const& v = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->x).children[scen]->first);
   if (!v.length())
      return std::vector<double>(); //this vector is not on this processor
   else
      return std::vector<double>(&v[0], &v[0] + v.length());
}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::printComplementarityResiduals(const sVars& svars) const {
   const int my_rank = PIPS_MPIgetRank();

   /* complementarity residuals before postsolve */
   std::unique_ptr<OoqpVectorBase<double>> t_clone{svars.t->cloneFull()};
   std::unique_ptr<OoqpVectorBase<double>> u_clone{svars.u->cloneFull()};
   std::unique_ptr<OoqpVectorBase<double>> v_clone{svars.v->cloneFull()};
   std::unique_ptr<OoqpVectorBase<double>> w_clone{svars.w->cloneFull()};

   t_clone->componentMult(*svars.lambda);
   t_clone->selectNonZeros(*svars.iclow);

   u_clone->componentMult(*svars.pi);
   u_clone->selectNonZeros(*svars.icupp);

   v_clone->componentMult(*svars.gamma);
   v_clone->selectNonZeros(*svars.ixlow);

   w_clone->componentMult(*svars.phi);
   w_clone->selectNonZeros(*svars.ixupp);

   const double rlambda_infnorm = t_clone->infnorm();
   const double rpi_infnorm = u_clone->infnorm();
   const double rgamma_infnorm = v_clone->infnorm();
   const double rphi_infnorm = w_clone->infnorm();

   if (my_rank == 0) {
      std::cout << " rl norm = " << rlambda_infnorm << "\n";
      std::cout << " rp norm = " << rpi_infnorm << "\n";
      std::cout << " rg norm = " << rgamma_infnorm << "\n";
      std::cout << " rf norm = " << rphi_infnorm << "\n";
      std::cout << "\n";
   }
}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::postsolveComputedSolution() {
   const bool print_residuals = pips_options::getBoolParameter("POSTSOLVE_PRINT_RESIDS");
   const int my_rank = PIPS_MPIgetRank(comm);

   assert(original_problem);
   assert(presolved_problem);

#if !defined(NDEBUG) && defined(PRESOLVE_POSTSOLVE_ONLY) // todo : resids for C also need recomputation.. - s variable
   /* todo: randomize all vectors x since it has not actually been set to anything */
   vars->x->setToConstant(0.1);
   resids->evaluate(data, vars);
#endif

   if (unscaleUnpermNotHierVars == nullptr)
      this->getVarsUnscaledUnperm();

   if (unscaleUnpermNotHierResids == nullptr)
      this->getResidsUnscaledUnperm();

   if (postsolvedVars != nullptr || postsolvedResids != nullptr)
      return;

   if (postsolver == nullptr) {
      assert("no postsolver available" && 0);
      return;
   }

   if (print_residuals) {
      if (my_rank == 0)
         std::cout << "\n" << "Residuals before postsolve:" << "\n";
      residuals->evaluate(*presolved_problem.get(), vars.get(), print_residuals);
      printComplementarityResiduals(*vars);

      MPI_Barrier(comm);
      if (my_rank == 0)
         std::cout << "Residuals after unscaling/permuting:" << "\n";
      unscaleUnpermNotHierResids->evaluate(*dataUnpermNotHier.get(), unscaleUnpermNotHierVars.get(), print_residuals);
      printComplementarityResiduals(*unscaleUnpermNotHierVars);
   }

   MPI_Barrier(comm);
   const double t0_postsolve = MPI_Wtime();


   if (pips_options::getBoolParameter("HIERARCHICAL"))
      formulation_factory->switchToOriginalTree();

   dynamic_cast<sTreeCallbacks*>(formulation_factory->tree)->switchToOriginalData();
   formulation_factory->data = original_problem.get();

   postsolvedVars.reset(dynamic_cast<sVars*>( formulation_factory->makeVariables(original_problem.get())));

   postsolvedResids.reset(dynamic_cast<DistributedResiduals*>( formulation_factory->makeResiduals(original_problem.get())));
   postsolver->postsolve(*unscaleUnpermNotHierVars, *postsolvedVars, result);

   double obj_postsolved = original_problem->objective_value(postsolvedVars.get());

   MPI_Barrier(comm);
   const double t_postsolve = MPI_Wtime();

   if (my_rank == 0) {
      std::cout << "---postsolve time (in sec.): " << t_postsolve - t0_postsolve << "\n";
      std::cout << "Objective value after postsolve: " << obj_postsolved << "\n";
   }

   /* compute residuals for postprocessed solution and check for feasibility */
   if (print_residuals) {
      if (my_rank == 0)
         std::cout << "\n" << "Residuals after postsolve:" << "\n";
      postsolvedResids->evaluate(*original_problem.get(), postsolvedVars.get(), print_residuals);

      printComplementarityResiduals(*postsolvedVars);
   }
}

#endif
