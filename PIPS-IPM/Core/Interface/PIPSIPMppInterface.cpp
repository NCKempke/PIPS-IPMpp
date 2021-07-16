//
// Created by charlie on 12.05.21.
//

#include "PIPSIPMppInterface.hpp"
#include "PreprocessType.h"
#include "DistributedFactory.hpp"
#include "DistributedProblem.hpp"
#include "DistributedResiduals.hpp"
#include "DistributedVariables.h"
#include "PreprocessFactory.h"
#include "Scaler.hpp"
#include "PIPSIPMppOptions.h"
#include "PIPSIPMppSolver.hpp"
#include "DistributedTreeCallbacks.h"

#include <functional>
#include <memory>

PIPSIPMppInterface::PIPSIPMppInterface(DistributedInputTree* tree, InteriorPointMethodType mehrotra_heuristic, MPI_Comm comm, ScalerType
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
   if (presolver_type != PresolverType::NONE) {
      original_problem = factory->make_problem();

      MPI_Barrier(comm);
      const double t0_presolve = MPI_Wtime();

      if (postsolve)
         postsolver = preprocess_factory->make_postsolver(original_problem.get());

      presolver = preprocess_factory->make_presolver(*factory->tree, original_problem.get(), presolver_type, postsolver.get());

      // this changes the tree and switches it to the one of the presolved data!
      presolved_problem.reset(dynamic_cast<DistributedProblem*>(presolver->presolve()));

      MPI_Barrier(comm);
      const double t_presolve = MPI_Wtime();
      if (my_rank == 0)
         std::cout << "---presolve time (in sec.): " << t_presolve - t0_presolve << "\n";
   }
   else {
      presolved_problem = factory->make_problem();
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

   dataUnpermNotHier = presolved_problem->clone_full();

   // after identifying the linking structure switch to hierarchical data structure
   if (pipsipmpp_options::get_bool_parameter("PARDISO_FOR_GLOBAL_SC"))
      dynamic_cast<DistributedProblem&>(*presolved_problem).activateLinkStructureExploitation();

   // TODO : save "old" data somewhere?
   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL")) {
      if (my_rank == 0)
         std::cout << "Using hierarchical approach!\n";

      presolved_problem.reset(dynamic_cast<DistributedProblem*>(factory->switchToHierarchicalData(dynamic_cast<DistributedProblem*>(presolved_problem.release()))));

      if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL_PRINT_HIER_DATA"))
         presolved_problem->write_to_streamDense(std::cout);
   }

   variables = factory->make_variables(*presolved_problem);
#ifdef TIMING
   if( my_rank == 0 ) printf("variables created\n");
#endif

   residuals = factory->make_residuals(*presolved_problem);
#ifdef TIMING
   if( my_rank == 0 ) printf("resids created\n");
#endif

   scaler = preprocess_factory->make_scaler(*factory, *presolved_problem, scaler_type);

#ifdef TIMING
   if( my_rank == 0 ) printf("scaler created\n");
#endif

   if (scaler) {
      MPI_Barrier(comm);
      const double t0_scaling = MPI_Wtime();

      scaler->scale();

      MPI_Barrier(comm);
      const double t_scaling = MPI_Wtime();
      if (my_rank == 0) {
         std::cout << "---scaling time (in sec.): " << t_scaling - t0_scaling << "\n";
      }
      presolved_problem->print_ranges();
   }

   solver = std::make_unique<PIPSIPMppSolver>(*factory, *presolved_problem, mehrotra_heuristic, scaler.get());
#ifdef TIMING
   if( my_rank == 0 ) printf("solver created\n");
#endif

   MPI_Barrier(comm);
   const double t1 = MPI_Wtime();
   if (my_rank == 0)
      std::cout << "---reading time (in sec.): " << t1 - t0 << "\n";
}

PIPSIPMppInterface::~PIPSIPMppInterface() = default;

TerminationStatus PIPSIPMppInterface::run() {
   if (my_rank == 0)
      std::cout << "solving ...\n";

   if (my_rank == 0) {
      const auto& qp_presolved_problem = dynamic_cast<DistributedProblem&>(*presolved_problem);
      // TODO : use unlifted data....
      if (!pipsipmpp_options::get_bool_parameter("HIERARCHICAL")) {
         std::cout << "1st stage " << qp_presolved_problem.getLocalnx() << " variables, " << qp_presolved_problem.getLocalmy()
                   << " equality constraints, " << qp_presolved_problem.getLocalmz() << " inequality constraints.\n";

         const size_t nscens = qp_presolved_problem.children.size();
         if (nscens) {
            std::cout << "2nd stage " << qp_presolved_problem.children[0]->getLocalnx() << " variables, "
                      << qp_presolved_problem.children[0]->getLocalmy() << " equality constraints, " << qp_presolved_problem.children[0]->getLocalmz()
                      << " inequality constraints.\n";

            std::cout << nscens << " scenarios." << "\n";
            std::cout << "Total " << qp_presolved_problem.getLocalnx() + nscens * qp_presolved_problem.children[0]->getLocalnx() << " variables, "
                      << qp_presolved_problem.getLocalmy() + nscens * qp_presolved_problem.children[0]->getLocalmy() << " equality constraints, "
                      << qp_presolved_problem.getLocalmz() + nscens * qp_presolved_problem.children[0]->getLocalmz() << " inequality constraints.\n";
         }
      }
   }
#ifdef TIMING
   double tmElapsed=MPI_Wtime();
#endif

#if defined(PRESOLVE_POSTSOLVE_ONLY) && !defined(NDEBUG)
   const int result = TerminationStatus::DID_NOT_RUN;
#else
   //---------------------------------------------
   result = solver->solve(*presolved_problem, *variables, *residuals);
   //---------------------------------------------
#endif

   if (result != TerminationStatus::SUCCESSFUL_TERMINATION && my_rank == 0)
      std::cout << "failed to solve instance, result code: " << result << "\n";

   ran_solver = true;

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

   if (postsolver && !postsolved_variables)
      this->postsolveComputedSolution();

   double obj;
   if (postsolved_variables != nullptr)
      obj = original_problem->evaluate_objective(*postsolved_variables);
   else {
      obj = presolved_problem->evaluate_objective(*variables);
      if (scaler)
         obj = scaler->get_unscaled_objective(obj);
   }

   return obj;
}

double PIPSIPMppInterface::getFirstStageObjective() const {
   const auto& x = *(dynamic_cast<DistributedVector<double>&>(*variables->primals).first);
   const auto& c = *(dynamic_cast<DistributedVector<double>&>(*presolved_problem->objective_gradient).first);

   return c.dotProductWith(x);
}

void PIPSIPMppInterface::getVarsUnscaledUnperm() {
   assert(!unscaleUnpermNotHierVars);
   assert(dataUnpermNotHier);

   const auto& qp_presolved_problem = dynamic_cast<DistributedProblem&>(*presolved_problem);

   if (!ran_solver)
      throw std::logic_error("Must call run() and start solution process before trying to retrieve unscaled unpermuted solution");
   if (scaler) {
      std::unique_ptr<DistributedVariables> unscaled_vars = std::make_unique<DistributedVariables>(dynamic_cast<const DistributedVariables&>(*variables));
      scaler->unscale_variables(*unscaled_vars);
      unscaleUnpermNotHierVars.reset(qp_presolved_problem.getVarsUnperm(*unscaled_vars, *dataUnpermNotHier));
   }
   else
      unscaleUnpermNotHierVars.reset(qp_presolved_problem.getVarsUnperm(*variables, *dataUnpermNotHier));

}

void PIPSIPMppInterface::getResidsUnscaledUnperm() {
   assert(!unscaleUnpermNotHierResids);
   assert(dataUnpermNotHier);

   const auto& qp_presolved_problem = dynamic_cast<DistributedProblem&>(*presolved_problem);

   if (!ran_solver)
      throw std::logic_error("Must call run() and start solution process before trying to retrieve unscaled unpermuted residuals");
   if (scaler) {
      std::unique_ptr<DistributedResiduals> unscaled_residuals = std::make_unique<DistributedResiduals>(dynamic_cast<const DistributedResiduals&>(*residuals));
      scaler->unscale_residuals(*unscaled_residuals);
      unscaleUnpermNotHierResids.reset(qp_presolved_problem.getResidsUnperm(*unscaled_residuals, *dataUnpermNotHier));
   }
   else
      unscaleUnpermNotHierResids.reset(qp_presolved_problem.getResidsUnperm(*residuals, *dataUnpermNotHier));
}

std::vector<double> PIPSIPMppInterface::gatherFromSolution(std::unique_ptr<Vector<double>> Variables::* member_to_gather) {
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

std::vector<double> PIPSIPMppInterface::gatherFromResiduals(std::unique_ptr<Vector<double>> Residuals::* member_to_gather) {
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
                                        ? dynamic_cast<DistributedVector<double>*>(unscaleUnpermNotHierResids->equality_residuals->clone_full())
                                        : dynamic_cast<DistributedVector<double>*>(postsolvedResids->equality_residuals->clone_full());

   if (!original_problem || !postsolved_variables)
      eq_vals->add(1.0, *presolved_problem->equality_rhs);
   else
      eq_vals->add(1.0, *original_problem->equality_rhs);

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
                                          ? dynamic_cast<DistributedVector<double>*>(unscaleUnpermNotHierResids->inequality_residuals->clone_full())
                                          : dynamic_cast<DistributedVector<double>*>(postsolvedResids->inequality_residuals->clone_full());

   if (postsolved_variables == nullptr)
      ineq_vals->add(1.0, *unscaleUnpermNotHierVars->slacks);
   else
      ineq_vals->add(1.0, *postsolved_variables->slacks);

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
   auto const& v = dynamic_cast<const DenseVector<double>&>(*dynamic_cast<DistributedVector<double> const&>(*variables->primals).first);
   return std::vector<double>(&v[0], &v[0] + v.length());
}

std::vector<double> PIPSIPMppInterface::getSecondStagePrimalColSolution(int scen) const {
   auto const& v = dynamic_cast<const DenseVector<double>&>(*dynamic_cast<DistributedVector<double> const&>(*variables->primals).children[scen]->first);
   if (!v.length())
      return std::vector<double>(); //this vector is not on this processor
   else
      return std::vector<double>(&v[0], &v[0] + v.length());
}

void PIPSIPMppInterface::allgatherBlocksizes(std::vector<unsigned int>& block_lengths_col,
   std::vector<unsigned int>& block_lengths_A, std::vector<unsigned int>& block_lengths_C) const
{
   /// gather col lengths
   const auto& col_vec = presolver ? dynamic_cast<const DistributedVector<double>&>(*original_problem->objective_gradient) :
      dynamic_cast<const DistributedVector<double>&>(*presolved_problem->objective_gradient);

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
   const auto& row_A_vec = presolver ? dynamic_cast<const DistributedVector<double>&>(*original_problem->equality_rhs) :
      dynamic_cast<const DistributedVector<double>&>(*presolved_problem->equality_rhs);
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

void PIPSIPMppInterface::printComplementarityResiduals(const Variables& svars) {
   const int my_rank = PIPS_MPIgetRank();

   std::unique_ptr<Vector<double>> t_clone{svars.slack_lower_bound_gap->clone_full()};
   std::unique_ptr<Vector<double>> u_clone{svars.slack_upper_bound_gap->clone_full()};
   std::unique_ptr<Vector<double>> v_clone{svars.primal_lower_bound_gap->clone_full()};
   std::unique_ptr<Vector<double>> w_clone{svars.primal_upper_bound_gap->clone_full()};

   t_clone->componentMult(*svars.slack_lower_bound_gap_dual);
   t_clone->selectNonZeros(*svars.inequality_lower_bound_indicators);

   u_clone->componentMult(*svars.slack_upper_bound_gap_dual);
   u_clone->selectNonZeros(*svars.inequality_upper_bound_indicators);

   v_clone->componentMult(*svars.primal_lower_bound_gap_dual);
   v_clone->selectNonZeros(*svars.primal_lower_bound_indicators);

   w_clone->componentMult(*svars.primal_upper_bound_gap_dual);
   w_clone->selectNonZeros(*svars.primal_upper_bound_indicators);

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

   if (!unscaleUnpermNotHierVars)
      this->getVarsUnscaledUnperm();

   if (!unscaleUnpermNotHierResids)
      this->getResidsUnscaledUnperm();

   if (postsolved_variables || postsolvedResids)
      return;

   if (!postsolver) {
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

   dynamic_cast<DistributedTreeCallbacks&>(*factory->tree).switchToOriginalData();

   postsolved_variables = factory->make_variables(*original_problem);

   postsolvedResids = factory->make_residuals(*original_problem);
   postsolver->postsolve(*unscaleUnpermNotHierVars, *postsolved_variables, result);

   double obj_postsolved = original_problem->evaluate_objective(*postsolved_variables);

   MPI_Barrier(comm);
   const double t_postsolve = MPI_Wtime();

   if (my_rank == 0) {
      std::cout << "---postsolve time (in sec.): " << t_postsolve - t0_postsolve << "\n";
      std::cout << "Objective value after postsolve: " << obj_postsolved << "\n";
   }

   /* compute residuals for postprocessed solution and check for feasibility */
   if (print_residuals) {
      if (my_rank == 0)
         std::cout << "\n" << "Residuals after postsolve:\n";
      postsolvedResids->evaluate(*original_problem, *postsolved_variables, print_residuals);

      printComplementarityResiduals(*postsolved_variables);
   }
}