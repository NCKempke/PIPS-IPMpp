//
// Created by charlie on 12.05.21.
//

#include "PIPSIPMppInterface.hpp"
#include "PreprocessType.h"
#include "DistributedFactory.hpp"
#include "DistributedQP.hpp"
#include "DistributedResiduals.hpp"
#include "DistributedVariables.h"
#include "PreprocessFactory.h"
#include "Scaler.hpp"
#include "PIPSIPMppOptions.h"
#include "InteriorPointMethod.hpp"
#include "DistributedTreeCallbacks.h"

#include <functional>
#include <memory>

PIPSIPMppInterface::PIPSIPMppInterface(DistributedInputTree* tree, MehrotraStrategyType mehrotra_heuristic, MPI_Comm comm, ScalerType
scaler_type, PresolverType presolver_type, const std::string& settings) : comm(comm), my_rank(PIPS_MPIgetRank(comm)) {
   factory = std::make_unique<DistributedFactory>(tree, comm);
   pipsipmpp_options::set_options(settings);
   const bool postsolve = pipsipmpp_options::get_bool_parameter("POSTSOLVE");

   MPI_Barrier(comm);
   const double t0 = MPI_Wtime();

#ifdef TIMING
   if( my_rank == 0 ) printf("factory created\n");
#endif

   preprocess_factory = std::make_unique<PreprocessFactory>();
#ifdef TIMING
   if( my_rank == 0 ) printf("prefactory created\n");
#endif

   // presolving activated?
   if (presolver_type != PresolverType::PRESOLVER_NONE) {
      original_problem.reset(dynamic_cast<DistributedQP*>(factory->make_problem()));

      MPI_Barrier(comm);
      const double t0_presolve = MPI_Wtime();

      if (postsolve)
         postsolver.reset(preprocess_factory->make_postsolver(original_problem.get()));

      presolver.reset(preprocess_factory->make_presolver(factory->tree, original_problem.get(), presolver_type, postsolver.get()));

      presolved_problem.reset(dynamic_cast<DistributedQP*>(presolver->presolve()));

      MPI_Barrier(comm);
      const double t_presolve = MPI_Wtime();
      if (my_rank == 0)
         std::cout << "---presolve time (in sec.): " << t_presolve - t0_presolve << "\n";
   }
   else {
      presolved_problem.reset(dynamic_cast<DistributedQP*>(factory->make_problem()));
   }
   assert(presolved_problem);

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

      presolved_problem.reset(dynamic_cast<DistributedQP*>(factory->switchToHierarchicalData(presolved_problem.release())));

      if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL_PRINT_HIER_DATA"))
         presolved_problem->write_to_streamDense(std::cout);
   }

   variables.reset(dynamic_cast<DistributedVariables*>(factory->make_variables(*presolved_problem)));
#ifdef TIMING
   if( my_rank == 0 ) printf("variables created\n");
#endif

   residuals.reset(dynamic_cast<DistributedResiduals*>(factory->make_residuals(*presolved_problem)));
#ifdef TIMING
   if( my_rank == 0 ) printf("resids created\n");
#endif

   scaler.reset(preprocess_factory->make_scaler(*presolved_problem, scaler_type));

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

   solver = std::make_unique<InteriorPointMethod>(*factory, *presolved_problem, mehrotra_heuristic, scaler.get());
#ifdef TIMING
   if( my_rank == 0 ) printf("solver created\n");
#endif

   MPI_Barrier(comm);
   const double t1 = MPI_Wtime();
   if (my_rank == 0)
      std::cout << "---reading time (in sec.): " << t1 - t0 << "\n";
}

PIPSIPMppInterface::~PIPSIPMppInterface() = default;;

TerminationStatus PIPSIPMppInterface::run() {
   if (my_rank == 0)
      std::cout << "solving ...\n";

   if (my_rank == 0) {
      // TODO : use unlifted data....
      if (!pipsipmpp_options::get_bool_parameter("HIERARCHICAL")) {
         std::cout << "1st stage " << presolved_problem->getLocalnx() << " variables, " << presolved_problem->getLocalmy()
                   << " equality constraints, " << presolved_problem->getLocalmz() << " inequality constraints.\n";

         const size_t nscens = presolved_problem->children.size();
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

   if (result != TerminationStatus::SUCCESSFUL_TERMINATION && my_rank == 0)
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

   return result;
}

int PIPSIPMppInterface::n_iterations() const {
   if(!ran_solver)
      throw std::logic_error("Must call run() and start solution process before trying to retrieve the iteration count!");

   return solver->n_iterations();
}

double PIPSIPMppInterface::getObjective() {

   if (!ran_solver)
      throw std::logic_error("Must call run() and start solution process before trying to retrieve original solution");

   if (postsolver != nullptr && postsolved_variables == nullptr)
      this->postsolveComputedSolution();

   double obj;
   if (postsolved_variables != nullptr)
      obj = original_problem->objective_value(*postsolved_variables);
   else {
      obj = presolved_problem->objective_value(*variables);
      if (scaler)
         obj = scaler->get_unscaled_objective(obj);
   }

   return obj;
}

double PIPSIPMppInterface::getFirstStageObjective() const {
   Vector<double>& x = *(dynamic_cast<DistributedVector<double>&>(*variables->primals).first);
   Vector<double>& c = *(dynamic_cast<DistributedVector<double>&>(*presolved_problem->g).first);
   return c.dotProductWith(x);
}

void PIPSIPMppInterface::getVarsUnscaledUnperm() {
   assert(unscaleUnpermNotHierVars == nullptr);
   assert(dataUnpermNotHier);

   if (!ran_solver)
      throw std::logic_error("Must call run() and start solution process before trying to retrieve unscaled unpermuted solution");
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
      throw std::logic_error("Must call run() and start solution process before trying to retrieve unscaled unpermuted residuals");
   if (scaler) {
      std::unique_ptr<DistributedResiduals> unscaled_resids{dynamic_cast<DistributedResiduals*>(scaler->get_unscaled_residuals(*residuals))};
      unscaleUnpermNotHierResids.reset(presolved_problem->getResidsUnperm(*unscaled_resids, *dataUnpermNotHier));
   }
   else
      unscaleUnpermNotHierResids.reset(presolved_problem->getResidsUnperm(*residuals, *dataUnpermNotHier));
}

std::vector<double> PIPSIPMppInterface::gatherFromSolution(std::unique_ptr<Vector<double>> DistributedVariables::* member_to_gather) {
   if (!unscaleUnpermNotHierVars)
      this->getVarsUnscaledUnperm();

   if (postsolver && !postsolved_variables)
      this->postsolveComputedSolution();

   std::vector<double> vec;
   if (!postsolver)
      vec = dynamic_cast<const DistributedVector<double>&>(*(*unscaleUnpermNotHierVars.*member_to_gather)).gatherStochVector();
   else
      vec = dynamic_cast<const DistributedVector<double>&>(*(*postsolved_variables.*member_to_gather)).gatherStochVector();

   return vec;
}

std::vector<double> PIPSIPMppInterface::gatherFromResiduals(std::unique_ptr<Vector<double>> DistributedResiduals::* member_to_gather) {
   if(!unscaleUnpermNotHierResids)
      this->getResidsUnscaledUnperm();

   if(postsolver && !postsolved_variables)
      this->postsolveComputedSolution();

   std::vector<double> vec;
   if(!postsolver)
      vec = dynamic_cast<const DistributedVector<double>&>(*(*unscaleUnpermNotHierResids.*member_to_gather)).gatherStochVector();
   else
      vec = dynamic_cast<const DistributedVector<double>&>(*(*postsolvedResids.*member_to_gather)).gatherStochVector();

   return vec;
}

std::vector<double> PIPSIPMppInterface::gatherPrimalSolution() {
   return gatherFromSolution(&DistributedVariables::primals);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionEq() {
   return gatherFromSolution(&DistributedVariables::equality_duals);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionIneq() {
   return gatherFromSolution(&DistributedVariables::inequality_duals);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionIneqUpp() {
   return gatherFromSolution(&DistributedVariables::slack_upper_bound_gap_dual);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionIneqLow() {
   return gatherFromSolution(&DistributedVariables::slack_lower_bound_gap_dual);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionVarBounds() {
   std::vector<double> duals_varbounds_upp = gatherDualSolutionVarBoundsUpp();
   std::vector<double> duals_varbounds_low = gatherDualSolutionVarBoundsLow();

   assert(duals_varbounds_low.size() == duals_varbounds_upp.size());

   std::vector<double> duals_varbounds;
   duals_varbounds.reserve(duals_varbounds_low.size());

   std::transform(duals_varbounds_low.begin(), duals_varbounds_low.end(), duals_varbounds_upp.begin(), std::back_inserter(duals_varbounds),
         std::minus<>());

   return duals_varbounds;
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionVarBoundsUpp() {
   return gatherFromSolution(&DistributedVariables::primal_upper_bound_gap_dual);
}


std::vector<double> PIPSIPMppInterface::gatherDualSolutionVarBoundsLow() {
   return gatherFromSolution(&DistributedVariables::primal_lower_bound_gap_dual);
}

std::vector<double> PIPSIPMppInterface::gatherEqualityConsValues() {
   if (unscaleUnpermNotHierResids == nullptr)
      this->getResidsUnscaledUnperm();

   if (postsolver != nullptr && postsolved_variables == nullptr)
      this->postsolveComputedSolution();

   DistributedVector<double>* eq_vals = (postsolved_variables == nullptr)
                                        ? dynamic_cast<DistributedVector<double>*>(unscaleUnpermNotHierResids->equality_residuals->cloneFull())
                                        : dynamic_cast<DistributedVector<double>*>(postsolvedResids->equality_residuals->cloneFull());

   if (!original_problem || !postsolved_variables)
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
                                          ? dynamic_cast<DistributedVector<double>*>(unscaleUnpermNotHierResids->inequality_residuals->cloneFull())
                                          : dynamic_cast<DistributedVector<double>*>(postsolvedResids->inequality_residuals->cloneFull());

   if (postsolved_variables == nullptr)
      ineq_vals->axpy(1.0, *unscaleUnpermNotHierVars->slacks);
   else
      ineq_vals->axpy(1.0, *postsolved_variables->slacks);

   std::vector<double> ineq_vals_vec = ineq_vals->gatherStochVector();

   delete ineq_vals;

   return ineq_vals_vec;
}

std::vector<double> PIPSIPMppInterface::gatherSlacksInequalityUp() {
   return gatherFromSolution(&DistributedVariables::slack_upper_bound_gap);
}

std::vector<double> PIPSIPMppInterface::gatherSlacksInequalityLow() {
   return gatherFromSolution(&DistributedVariables::slack_lower_bound_gap);
}

std::vector<double> PIPSIPMppInterface::gatherSlacksVarsUp() {
   return gatherFromSolution(&DistributedVariables::primal_upper_bound_gap);
}

std::vector<double> PIPSIPMppInterface::gatherSlacksVarsLow() {
   return gatherFromSolution(&DistributedVariables::primal_lower_bound_gap);

}

std::vector<double> PIPSIPMppInterface::gatherPrimalResidsEQ() {
   return gatherFromResiduals(&DistributedResiduals::equality_residuals);
}

std::vector<double> PIPSIPMppInterface::gatherPrimalResidsIneqUp() {
   return gatherFromResiduals(&DistributedResiduals::ru);
}

std::vector<double> PIPSIPMppInterface::gatherPrimalResidsIneqLow() {
   return gatherFromResiduals(&DistributedResiduals::rt);
}

std::vector<double> PIPSIPMppInterface::gatherDualResids() {
   return gatherFromResiduals(&DistributedResiduals::lagrangian_gradient);
}

std::vector<double> PIPSIPMppInterface::getFirstStagePrimalColSolution() const {
   auto const& v = dynamic_cast<const SimpleVector<double>&>(*dynamic_cast<DistributedVector<double> const&>(*variables->primals).first);
   return std::vector<double>(&v[0], &v[0] + v.length());
}

std::vector<double> PIPSIPMppInterface::getSecondStagePrimalColSolution(int scen) const {
   auto const& v = dynamic_cast<const SimpleVector<double>&>(*dynamic_cast<DistributedVector<double> const&>(*variables->primals).children[scen]->first);
   if (!v.length())
      return std::vector<double>(); //this vector is not on this processor
   else
      return std::vector<double>(&v[0], &v[0] + v.length());
}

void PIPSIPMppInterface::allgatherBlocksizes(std::vector<unsigned int>& block_lengths_col,
   std::vector<unsigned int>& block_lengths_A, std::vector<unsigned int>& block_lengths_C) const
{
   /// gather col lengths
   const auto& col_vec = presolver ? dynamic_cast<const DistributedVector<double>&>(*original_problem->g) :
      dynamic_cast<const DistributedVector<double>&>(*presolved_problem->g);

   assert( block_lengths_col.size() == col_vec.children.size() + 1);
   assert( !col_vec.last );

   if( my_rank == 0 )
      block_lengths_col[0] = col_vec.first->length();

   for( unsigned int i = 0; i < col_vec.children.size(); ++i )
   {
      if( !col_vec.children[i]->isKindOf(kStochDummy) )
      {
         assert( col_vec.children[i]->first );
         assert( !col_vec.children[i]->last );
         block_lengths_col[i + 1] = col_vec.children[i]->first->length();
      }
   }

   PIPS_MPIsumArrayInPlace(block_lengths_col, MPI_COMM_WORLD);

   /// gather row lengths
   const auto& row_A_vec = presolver ? dynamic_cast<const DistributedVector<double>&>(*original_problem->bA) :
      dynamic_cast<const DistributedVector<double>&>(*presolved_problem->bA);
   const auto& row_C_vec = presolver ? dynamic_cast<const DistributedVector<double>&>(*original_problem->inequality_upper_bound_indicators) :
      dynamic_cast<const DistributedVector<double>&>(*presolved_problem->inequality_upper_bound_indicators);

   assert( block_lengths_A.size() == row_A_vec.children.size() + 2);
   assert( block_lengths_C.size() == row_C_vec.children.size() + 2);
   assert( block_lengths_A.size() == block_lengths_C.size() );

   if( my_rank == 0 )
   {
      assert( row_A_vec.last );
      assert( row_C_vec.last );
      assert( row_A_vec.first );
      assert( row_C_vec.first );

      block_lengths_A[0] = row_A_vec.first->length();
      block_lengths_C[0] = row_C_vec.first->length();
      block_lengths_A[ row_A_vec.children.size() + 1] = row_A_vec.last->length();
      block_lengths_C[ row_C_vec.children.size() + 1] = row_C_vec.last->length();
   }

   for( unsigned int i = 0; i < row_A_vec.children.size(); ++i )
   {
      if (!row_A_vec.children[i]->isKindOf(kStochDummy)) {
         assert(row_A_vec.children[i]->first);
         assert(row_C_vec.children[i]->first);
         assert(!row_A_vec.children[i]->last);
         assert(!row_C_vec.children[i]->last);

         block_lengths_A[i + 1] = row_A_vec.children[i]->first->length();
         block_lengths_C[i + 1] = row_C_vec.children[i]->first->length();
      }
   }

   PIPS_MPIsumArrayInPlace(block_lengths_A, MPI_COMM_WORLD);
   PIPS_MPIsumArrayInPlace(block_lengths_C, MPI_COMM_WORLD);
}

void PIPSIPMppInterface::printComplementarityResiduals(const DistributedVariables& svars) {
   const int my_rank = PIPS_MPIgetRank();

   /* complementarity residuals before postsolve */
   std::unique_ptr<Vector<double>> t_clone{svars.slack_lower_bound_gap->cloneFull()};
   std::unique_ptr<Vector<double>> u_clone{svars.slack_upper_bound_gap->cloneFull()};
   std::unique_ptr<Vector<double>> v_clone{svars.primal_lower_bound_gap->cloneFull()};
   std::unique_ptr<Vector<double>> w_clone{svars.primal_upper_bound_gap->cloneFull()};

   t_clone->componentMult(*svars.slack_lower_bound_gap_dual);
   t_clone->selectNonZeros(*svars.iclow);

   u_clone->componentMult(*svars.slack_upper_bound_gap_dual);
   u_clone->selectNonZeros(*svars.icupp);

   v_clone->componentMult(*svars.primal_lower_bound_gap_dual);
   v_clone->selectNonZeros(*svars.ixlow);

   w_clone->componentMult(*svars.primal_upper_bound_gap_dual);
   w_clone->selectNonZeros(*svars.ixupp);

   const double rlambda_infnorm = t_clone->inf_norm();
   const double rpi_infnorm = u_clone->inf_norm();
   const double rgamma_infnorm = v_clone->inf_norm();
   const double rphi_infnorm = w_clone->inf_norm();

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
         std::cout << "\n" << "Residuals before postsolve:\n";
      residuals->evaluate(*presolved_problem, *variables, print_residuals);
      printComplementarityResiduals(*variables);

      MPI_Barrier(comm);
      if (my_rank == 0)
         std::cout << "Residuals after unscaling/permuting:\n";
      unscaleUnpermNotHierResids->evaluate(*dataUnpermNotHier, *unscaleUnpermNotHierVars, print_residuals);
      printComplementarityResiduals(*unscaleUnpermNotHierVars);
   }

   MPI_Barrier(comm);
   const double t0_postsolve = MPI_Wtime();

   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL"))
      factory->switchToOriginalTree();

   dynamic_cast<DistributedTreeCallbacks*>(factory->tree)->switchToOriginalData();

   postsolved_variables.reset(dynamic_cast<DistributedVariables*>(factory->make_variables(*original_problem)));

   postsolvedResids.reset(dynamic_cast<DistributedResiduals*>(factory->make_residuals(*original_problem)));
   postsolver->postsolve(*unscaleUnpermNotHierVars, *postsolved_variables, result);

   double obj_postsolved = original_problem->objective_value(*postsolved_variables);

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
      postsolvedResids->evaluate(*original_problem, *postsolved_variables, print_residuals);

      printComplementarityResiduals(*postsolved_variables);
   }
}