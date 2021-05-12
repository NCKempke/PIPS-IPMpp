/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef PIPSIPMPPINTERFACE_H
#define PIPSIPMPPINTERFACE_H

#include <algorithm>
#include <functional>
#include <cstdlib>
#include <stdexcept>
#include <memory>
#include <MehrotraStrategy.hpp>
#include <InteriorPointMethod.h>
#include "DistributedTree.h"
#include "DistributedQP.hpp"
#include "DistributedResiduals.hpp"
#include "DistributedVariables.h"
#include "Statistics.hpp"
#include "PreprocessFactory.h"
#include "Scaler.h"
#include "Presolver.h"
#include "Postsolver.h"
#include "DistributedTreeCallbacks.h"
#include "pipsport.h"
#include "PIPSIPMppOptions.h"

class PIPSIPMppInterface {
public:
   PIPSIPMppInterface(DistributedInputTree* tree, MehrotraHeuristic mehrotra_heuristic, MPI_Comm = MPI_COMM_WORLD, ScalerType scaler_type = SCALER_NONE,
         PresolverType presolver_type = PRESOLVER_NONE, std::string settings = "PIPSIPMpp.opt");

   ~PIPSIPMppInterface() = default;

   void run();

   double getObjective();

   [[nodiscard]] double getFirstStageObjective() const;

   std::vector<double> gatherPrimalSolution();

   std::vector<double> gatherDualSolutionEq();

   std::vector<double> gatherDualSolutionIneq();

   std::vector<double> gatherDualSolutionIneqUpp();

   std::vector<double> gatherDualSolutionIneqLow();

   std::vector<double> gatherDualSolutionVarBounds();

   std::vector<double> gatherDualSolutionVarBoundsUpp();

   std::vector<double> gatherDualSolutionVarBoundsLow();

   [[nodiscard]] std::vector<double> getFirstStagePrimalColSolution() const;

   [[nodiscard]] std::vector<double> getSecondStagePrimalColSolution(int scen) const;

   void postsolveComputedSolution();

   std::vector<double> gatherEqualityConsValues();

   std::vector<double> gatherInequalityConsValues();

   void getVarsUnscaledUnperm();

   void getResidsUnscaledUnperm();
   //more get methods to follow here

private:
   void printComplementarityResiduals(const DistributedVariables& vars) const;

   std::vector<double> gatherFromSolution(SmartPointer<Vector<double> > DistributedVariables::* member_to_gather);

protected:
   DistributedFactory factory;
   std::unique_ptr<PreprocessFactory> preprocess_factory{};

   std::unique_ptr<DistributedQP> presolved_problem{};       // possibly presolved problem
   std::unique_ptr<DistributedQP> dataUnpermNotHier{}; // data after presolve before permutation, scaling and hierarchical data
   std::unique_ptr<DistributedQP> original_problem{};   // original data
   std::unique_ptr<DistributedVariables> variables{};
   std::unique_ptr<DistributedVariables> unscaleUnpermNotHierVars{};
   std::unique_ptr<DistributedVariables> postsolved_variables{};

   std::unique_ptr<DistributedResiduals> residuals{};
   std::unique_ptr<DistributedResiduals> unscaleUnpermNotHierResids{};
   std::unique_ptr<DistributedResiduals> postsolvedResids{};

   std::unique_ptr<Presolver> presolver{};
   std::unique_ptr<Postsolver> postsolver{};
   std::unique_ptr<Scaler> scaler{};
   std::unique_ptr<InteriorPointMethod> solver{};

   MPI_Comm comm = MPI_COMM_NULL;
   const int my_rank = -1;
   int result{-1};
   bool ran_solver = false;
};

#endif // PIPSIPMPPINTERFACE_H


// implementation
#ifdef PIPSIPMPPINTERFACE_H

PIPSIPMppInterface::PIPSIPMppInterface(DistributedInputTree* tree, MehrotraHeuristic mehrotra_heuristic, MPI_Comm comm, ScalerType
scaler_type, PresolverType presolver_type, std::string settings) : factory(tree, comm), comm(comm), my_rank(PIPS_MPIgetRank()) {
   pipsipmpp_options::set_options(settings);
   const bool postsolve = pipsipmpp_options::get_bool_parameter("POSTSOLVE");

   MPI_Barrier(comm);
   const double t0 = MPI_Wtime();

#ifdef TIMING
   if( my_rank == 0 ) printf("factory created\n");
#endif

   preprocess_factory.reset(new PreprocessFactory());
#ifdef TIMING
   if( my_rank == 0 ) printf("prefactory created\n");
#endif

   // presolving activated?
   if (presolver_type != PRESOLVER_NONE) {
      original_problem.reset(dynamic_cast<DistributedQP*>(factory.make_problem()));

      MPI_Barrier(comm);
      const double t0_presolve = MPI_Wtime();

      if (postsolve)
         postsolver.reset(preprocess_factory->makePostsolver(original_problem.get()));

      presolver.reset(
            preprocess_factory->makePresolver(dynamic_cast<DistributedFactory*>(&factory)->tree, original_problem.get(), presolver_type,
                  postsolver.get()));

      presolved_problem.reset(dynamic_cast<DistributedQP*>(presolver->presolve()));

      factory.problem = presolved_problem.get(); // todo update also sTree* of factory

      MPI_Barrier(comm);
      const double t_presolve = MPI_Wtime();
      if (my_rank == 0)
         std::cout << "---presolve time (in sec.): " << t_presolve - t0_presolve << "\n";
   }
   else {
      presolved_problem.reset(dynamic_cast<DistributedQP*>(factory.make_problem()));
      assert(presolved_problem);
   }

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
   if (pipsipmpp_options::get_bool_parameter("PARDISO_FOR_GLOBAL_SC"))
      presolved_problem->activateLinkStructureExploitation();

   // TODO : save "old" data somewhere?
   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL")) {
      if (my_rank == 0)
         std::cout << "Using hierarchical approach!\n";

      presolved_problem.reset(dynamic_cast<DistributedQP*>(factory.switchToHierarchicalData(presolved_problem.release())));

      if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL_PRINT_HIER_DATA"))
         presolved_problem->writeToStreamDense(std::cout);
   }

   variables.reset(dynamic_cast<DistributedVariables*>( factory.make_variables(*presolved_problem)));
#ifdef TIMING
   if( my_rank == 0 ) printf("variables created\n");
#endif

   residuals.reset(dynamic_cast<DistributedResiduals*>( factory.make_residuals(*presolved_problem)));
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

   solver.reset(new InteriorPointMethod(factory, *presolved_problem, mehrotra_heuristic, scaler.get()));
#ifdef TIMING
   if( my_rank == 0 ) printf("solver created\n");
#endif

   MPI_Barrier(comm);
   const double t1 = MPI_Wtime();
   if (my_rank == 0)
      std::cout << "---reading time (in sec.): " << t1 - t0 << "\n";
}

void PIPSIPMppInterface::run() {
   if (my_rank == 0)
      std::cout << "solving ...\n";

   if (my_rank == 0) {
      // TODO : use unlifted data....
      if (!pipsipmpp_options::get_bool_parameter("HIERARCHICAL")) {
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
   result = solver->solve(*presolved_problem, *variables, *residuals);
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

double PIPSIPMppInterface::getObjective() {

   if (!ran_solver)
      throw std::logic_error("Must call go() and start solution process before trying to retrieve original solution");

   if (postsolver != nullptr && postsolved_variables == nullptr)
      this->postsolveComputedSolution();

   double obj;
   if (postsolved_variables != nullptr)
      obj = original_problem->objective_value(*postsolved_variables.get());
   else {
      obj = presolved_problem->objective_value(*variables.get());
      if (scaler)
         obj = scaler->get_unscaled_objective(obj);
   }

   return obj;
}

double PIPSIPMppInterface::getFirstStageObjective() const {
   Vector<double>& x = *(dynamic_cast<DistributedVector<double>&>(*variables->x).first);
   Vector<double>& c = *(dynamic_cast<DistributedVector<double>&>(*presolved_problem->g).first);
   return c.dotProductWith(x);
}

void PIPSIPMppInterface::getVarsUnscaledUnperm() {
   assert(unscaleUnpermNotHierVars == nullptr);
   assert(dataUnpermNotHier);

   if (!ran_solver)
      throw std::logic_error("Must call go() and start solution process before trying to retrieve unscaled unpermuted solution");
   if (scaler) {
      std::unique_ptr<DistributedVariables> unscaled_vars{dynamic_cast<DistributedVariables*>(scaler->get_unscaled_variables(*variables))};
      unscaleUnpermNotHierVars.reset(presolved_problem->getVarsUnperm(*unscaled_vars, *dataUnpermNotHier));
   }
   else
      unscaleUnpermNotHierVars.reset(presolved_problem->getVarsUnperm(*variables, *dataUnpermNotHier));

}

void PIPSIPMppInterface::getResidsUnscaledUnperm() {
   assert(unscaleUnpermNotHierResids == nullptr);
   assert(dataUnpermNotHier);

   if (!ran_solver)
      throw std::logic_error("Must call go() and start solution process before trying to retrieve unscaled unpermuted residuals");
   if (scaler) {
      std::unique_ptr<DistributedResiduals> unscaled_resids{dynamic_cast<DistributedResiduals*>(scaler->get_unscaled_residuals(*residuals))};
      unscaleUnpermNotHierResids.reset(presolved_problem->getResidsUnperm(*unscaled_resids, *dataUnpermNotHier));
   }
   else
      unscaleUnpermNotHierResids.reset(presolved_problem->getResidsUnperm(*residuals, *dataUnpermNotHier));
}

std::vector<double> PIPSIPMppInterface::gatherFromSolution(SmartPointer<Vector<double> > DistributedVariables::* member_to_gather) {
   if (unscaleUnpermNotHierVars == nullptr)
      this->getVarsUnscaledUnperm();

   if (postsolver != nullptr && postsolved_variables == nullptr)
      this->postsolveComputedSolution();

   std::vector<double> vec;
   if (postsolver == nullptr)
      vec = dynamic_cast<const DistributedVector<double>&>(*(unscaleUnpermNotHierVars.get()->*member_to_gather)).gatherStochVector();
   else
      vec = dynamic_cast<const DistributedVector<double>&>(*(postsolved_variables.get()->*member_to_gather)).gatherStochVector();

   return vec;
}


std::vector<double> PIPSIPMppInterface::gatherPrimalSolution() {
   return gatherFromSolution(&DistributedVariables::x);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionEq() {
   return gatherFromSolution(&DistributedVariables::y);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionIneq() {
   return gatherFromSolution(&DistributedVariables::z);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionIneqUpp() {
   return gatherFromSolution(&DistributedVariables::pi);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionIneqLow() {
   return gatherFromSolution(&DistributedVariables::lambda);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionVarBounds() {
   std::vector<double> duals_varbounds_upp = gatherDualSolutionVarBoundsUpp();
   std::vector<double> duals_varbounds_low = gatherDualSolutionVarBoundsLow();

   assert(duals_varbounds_low.size() == duals_varbounds_upp.size());

   std::vector<double> duals_varbounds;
   duals_varbounds.reserve(duals_varbounds_low.size());

   std::transform(duals_varbounds_low.begin(), duals_varbounds_low.end(), duals_varbounds_upp.begin(), std::back_inserter(duals_varbounds),
         std::minus<double>());

   return duals_varbounds;
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionVarBoundsUpp() {
   return gatherFromSolution(&DistributedVariables::phi);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionVarBoundsLow() {
   return gatherFromSolution(&DistributedVariables::gamma);
}


std::vector<double> PIPSIPMppInterface::gatherEqualityConsValues() {
   if (unscaleUnpermNotHierResids == nullptr)
      this->getResidsUnscaledUnperm();

   if (postsolver != nullptr && postsolved_variables == nullptr)
      this->postsolveComputedSolution();

   DistributedVector<double>* eq_vals = (postsolved_variables == nullptr)
                                        ? dynamic_cast<DistributedVector<double>*>(unscaleUnpermNotHierResids->rA->cloneFull())
                                        : dynamic_cast<DistributedVector<double>*>(postsolvedResids->rA->cloneFull());

   if (original_problem == nullptr || postsolved_variables == nullptr)
      eq_vals->axpy(1.0, *presolved_problem->bA);
   else
      eq_vals->axpy(1.0, *original_problem->bA);

   std::vector<double> eq_vals_vec = eq_vals->gatherStochVector();

   delete eq_vals;

   return eq_vals_vec;
}



std::vector<double> PIPSIPMppInterface::gatherInequalityConsValues() {
   if (unscaleUnpermNotHierVars == nullptr)
      this->getVarsUnscaledUnperm();

   if (unscaleUnpermNotHierResids == nullptr)
      this->getResidsUnscaledUnperm();

   if (postsolver != nullptr && postsolved_variables == nullptr)
      this->postsolveComputedSolution();

   DistributedVector<double>* ineq_vals = (postsolved_variables == nullptr)
                                          ? dynamic_cast<DistributedVector<double>*>(unscaleUnpermNotHierResids->rC->cloneFull())
                                          : dynamic_cast<DistributedVector<double>*>(postsolvedResids->rC->cloneFull());

   if (postsolved_variables == nullptr)
      ineq_vals->axpy(1.0, *unscaleUnpermNotHierVars->s);
   else
      ineq_vals->axpy(1.0, *postsolved_variables->s);

   std::vector<double> ineq_vals_vec = ineq_vals->gatherStochVector();

   delete ineq_vals;

   return ineq_vals_vec;
}


std::vector<double> PIPSIPMppInterface::getFirstStagePrimalColSolution() const {
   auto const& v = *dynamic_cast<SimpleVector<double> const*>(dynamic_cast<DistributedVector<double> const&>(*variables->x).first);
   return std::vector<double>(&v[0], &v[0] + v.length());
}


std::vector<double> PIPSIPMppInterface::getSecondStagePrimalColSolution(int scen) const {
   auto const& v = *dynamic_cast<SimpleVector<double> const*>(dynamic_cast<DistributedVector<double> const&>(*variables->x).children[scen]->first);
   if (!v.length())
      return std::vector<double>(); //this vector is not on this processor
   else
      return std::vector<double>(&v[0], &v[0] + v.length());
}


void PIPSIPMppInterface::printComplementarityResiduals(const DistributedVariables& svars) const {
   const int my_rank = PIPS_MPIgetRank();

   /* complementarity residuals before postsolve */
   std::unique_ptr<Vector<double>> t_clone{svars.t->cloneFull()};
   std::unique_ptr<Vector<double>> u_clone{svars.u->cloneFull()};
   std::unique_ptr<Vector<double>> v_clone{svars.v->cloneFull()};
   std::unique_ptr<Vector<double>> w_clone{svars.w->cloneFull()};

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


void PIPSIPMppInterface::postsolveComputedSolution() {
   const bool print_residuals = pipsipmpp_options::get_bool_parameter("POSTSOLVE_PRINT_RESIDS");
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

   if (postsolved_variables != nullptr || postsolvedResids != nullptr)
      return;

   if (postsolver == nullptr) {
      assert("no postsolver available" && 0);
      return;
   }

   if (print_residuals) {
      if (my_rank == 0)
         std::cout << "\n" << "Residuals before postsolve:" << "\n";
      residuals->evaluate(*presolved_problem.get(), *variables.get(), print_residuals);
      printComplementarityResiduals(*variables);

      MPI_Barrier(comm);
      if (my_rank == 0)
         std::cout << "Residuals after unscaling/permuting:" << "\n";
      unscaleUnpermNotHierResids->evaluate(*dataUnpermNotHier.get(), *unscaleUnpermNotHierVars.get(), print_residuals);
      printComplementarityResiduals(*unscaleUnpermNotHierVars);
   }

   MPI_Barrier(comm);
   const double t0_postsolve = MPI_Wtime();


   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL"))
      factory.switchToOriginalTree();

   dynamic_cast<DistributedTreeCallbacks*>(factory.tree)->switchToOriginalData();
   factory.problem = original_problem.get();

   postsolved_variables.reset(dynamic_cast<DistributedVariables*>( factory.make_variables(*original_problem)));

   postsolvedResids.reset(dynamic_cast<DistributedResiduals*>(factory.make_residuals(*original_problem)));
   postsolver->postsolve(*unscaleUnpermNotHierVars, *postsolved_variables, result);

   double obj_postsolved = original_problem->objective_value(*postsolved_variables.get());

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
      postsolvedResids->evaluate(*original_problem.get(), *postsolved_variables.get(), print_residuals);

      printComplementarityResiduals(*postsolved_variables);
   }
}

#endif