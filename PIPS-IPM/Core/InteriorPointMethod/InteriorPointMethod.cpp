#include "InteriorPointMethod.hpp"
#include "Status.h"
#include "DistributedOptions.h"
#include "DistributedRootLinearSystem.h"
#include "DistributedFactory.h"

extern int gOoqpPrintLevel;
extern double g_iterNumber;

InteriorPointMethod::InteriorPointMethod(DistributedFactory& problem_formulation, Problem& problem, const Scaler* scaler) : Solver(problem_formulation,
      problem, scaler), step_length_type(PRIMAL), // TO CHANGE
      n_linesearch_points(pips_options::get_int_parameter("GONDZIO_STOCH_N_LINESEARCH")),
      dynamic_corrector_schedule(pips_options::get_bool_parameter("GONDZIO_STOCH_USE_DYNAMIC_CORRECTOR_SCHEDULE")),
      additional_correctors_small_comp_pairs(pips_options::get_bool_parameter("GONDZIO_STOCH_ADDITIONAL_CORRECTORS_SMALL_VARS")),
      max_additional_correctors(pips_options::get_int_parameter("GONDZIO_STOCH_ADDITIONAL_CORRECTORS_MAX")),
      first_iter_small_correctors(pips_options::get_int_parameter("GONDZIO_STOCH_FIRST_ITER_SMALL_CORRECTORS")),
      max_alpha_small_correctors(pips_options::get_double_parameter("GONDZIO_STOCH_MAX_ALPHA_SMALL_CORRECTORS")), NumberSmallCorrectors(0),
      push_converged_vars_from_bound(pips_options::get_bool_parameter("GONDZIO_STOCH_PUSH_CONVERGED_VARS_FROM_BOUND")),
      fequency_push_converged_vars_from_bound(pips_options::get_int_parameter("GONDZIO_STOCH_FREQUENCY_PUSH_CONVERGED_VARS")),
      mu_limit_push_converged_vars_from_bound(pips_options::get_double_parameter("GONDZIO_STOCH_MU_LIMIT_PUSH_CONVERGED_VARS")),
      bicgstab_skipped(false), bicgstab_converged(true), bigcstab_norm_res_rel(0.), bicg_iterations(0),
      dynamic_bicg_tol(pips_options::get_bool_parameter("OUTER_BICG_DYNAMIC_TOL")) {
   assert(max_additional_correctors > 0);
   assert(first_iter_small_correctors >= 0);
   assert(0 < max_alpha_small_correctors && max_alpha_small_correctors < 1);
   assert(n_linesearch_points > 0);

   if (pips_options::get_bool_parameter("GONDZIO_STOCH_ADAPTIVE_LINESEARCH")) {
      const int size = PIPS_MPIgetSize();

      if (size > 1)
         this->n_linesearch_points = std::min(unsigned(size) + this->n_linesearch_points, max_linesearch_points);
   }

   // the two StepFactor constants set targets for increase in step length for each corrector
   step_factor0 = 0.3;
   step_factor1 = 1.5;

   temp_step = factory.make_variables(problem);
}

TerminationCode InteriorPointMethod::solve(Problem& problem, Variables& iterate, Residuals& residuals) {
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   double mu, muaff;
   TerminationCode status_code;
   double sigma = 1.;
   double step_length_primal = 1., step_length_dual = 1.;

   DistributedFactory* distributed_factory = dynamic_cast<DistributedFactory*>(&factory);
   g_iterNumber = 0.;

   bool pure_centering_step = false;
   bool numerical_troubles = false;
   bool precond_decreased = true;

   set_problem_norm(problem);

   // initialization of (x,y,z) and factorization routine.
   linear_system = factory.make_linear_system(problem);

   // register as observer for the BiCGStab solves
   registerBiCGStabOvserver(linear_system);
   setBiCGStabTol(-1);

   distributed_factory->iterate_started();
   this->solve_linear_system(iterate, problem, residuals, *step);
   distributed_factory->iterate_ended();

   iteration = 0;
   number_gondzio_corrections = 0;
   mu = iterate.mu();

   while (true) {
      iteration++;

      setBiCGStabTol(iteration);
      bool small_corr = false;

      distributed_factory->iterate_started();

      // evaluate residuals and update algorithm status:
      residuals.evaluate(problem, iterate);

      //  termination test:
      status_code = this->do_status(&problem, &iterate, &residuals, iteration, mu, SUCCESSFUL_TERMINATION);

      if (status_code != NOT_FINISHED)
         break;

      std::cout << "Step 1: " << step_length_primal << ", " << step_length_dual << "\n";
      if (gOoqpPrintLevel >= 10) {
         this->do_monitor_Pd(&problem, &iterate, &residuals, step_length_primal, step_length_dual, sigma, iteration, mu, status_code, 0);
      }

      // *** Predictor step ***
      if (!pure_centering_step) {
         compute_predictor_step(problem, iterate, residuals);
         check_linear_system_solve_numerical_troubles(&residuals, numerical_troubles, small_corr);
      }
      else
         step->setToZero();

      // compute the step length
      if (step_length_type == PRIMAL) {
         double step_length = iterate.stepbound(step);
         step_length_primal = step_length;
         step_length_dual = step_length;
      }
      else {
         iterate.stepbound_pd(step, step_length_primal, step_length_dual);
      }
      std::cout << "Step 2: " << step_length_primal << ", " << step_length_dual << "\n";

      // calculate centering parameter
      muaff = iterate.mustep_pd(step, step_length_primal, step_length_dual);

      assert(!PIPSisZero(mu));
      sigma = pow(muaff / mu, tsig);

      if (gOoqpPrintLevel >= 10) {
         // TODO this varies
         this->do_monitor_Pd(&problem, &iterate, &residuals, step_length_primal, step_length_dual, sigma, iteration, mu, status_code, 2);
      }

      g_iterNumber += 1.;

      // *** Corrector step ***
      compute_corrector_step(problem, iterate, sigma, mu);
      check_linear_system_solve_numerical_troubles(&residuals, numerical_troubles, small_corr);

      // calculate weighted predictor-corrector step
      double weight_primal_candidate, weight_dual_candidate = -1.;
      if (step_length_type == PRIMAL) {
         calculate_alpha_PD_weight_candidate(&iterate, step, corrector_step, step_length_primal, step_length_dual, step_length_primal,
               step_length_dual, weight_primal_candidate, weight_dual_candidate);
      }
      else {
         calculate_alpha_PD_weight_candidate(&iterate, step, corrector_step, step_length_primal, step_length_dual, step_length_primal,
               step_length_dual, weight_primal_candidate, weight_dual_candidate);
      }
      std::cout << "Step 3: " << step_length_primal << ", " << step_length_dual << "\n";

      assert(weight_primal_candidate >= 0. && weight_primal_candidate <= 1.);
      assert(weight_dual_candidate >= 0. && weight_dual_candidate <= 1.);

      step->saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

      // prepare for Gondzio corrector loop: zero out the corrector_residuals structure:
      corrector_residuals->clear_linear_residuals();

      // calculate the target box:
      const double rmin = sigma * mu * beta_min;
      const double rmax = sigma * mu * beta_max;

      number_gondzio_corrections = 0;
      NumberSmallCorrectors = 0;

      // if small_corr_aggr only try small correctors
      // enter the Gondzio correction loop:
      while (number_gondzio_corrections < maximum_correctors && NumberSmallCorrectors < max_additional_correctors &&
             (PIPSisLT(step_length_primal, 1.) || PIPSisLT(step_length_dual, 1.))) {
         if (dynamic_corrector_schedule)
            adjustLimitGondzioCorrectors();
         corrector_step->copy(&iterate);

         // calculate target steplength
         const double step_length_primal_target = std::min(1., step_factor1 * step_length_primal + step_factor0);
         const double step_length_dual_target = std::min(1., step_factor1 * step_length_dual + step_factor0);

         PIPSdebugMessage("corrector loop: %d step_length_primal: %f step_length_dual %f \n", number_gondzio_corrections, step_length_primal,
                  step_length_dual);

         // add a step of this length to corrector_step
         corrector_step->saxpy_pd(step, step_length_primal_target, step_length_dual_target);
         // corrector_step is now x_k + alpha_target * delta_p (a trial point)

         /* compute corrector step */
         compute_gondzio_corrector(&problem, &iterate, rmin, rmax, small_corr);
         const bool was_small_corr = small_corr;
         check_linear_system_solve_numerical_troubles(&residuals, numerical_troubles, small_corr);

         if (numerical_troubles) {
            if (!was_small_corr && small_corr)
               continue;
            else
               // exit corrector loop if small correctors have already been tried or are not allowed
               break;
         }

         // calculate weighted predictor-corrector step
         double step_length_primal_enhanced, step_length_dual_enhanced;
         calculate_alpha_PD_weight_candidate(&iterate, step, corrector_step, step_length_primal_target, step_length_dual_target,
               step_length_primal_enhanced, step_length_dual_enhanced, weight_primal_candidate, weight_dual_candidate);

         // if the enhanced step length is actually 1, make it official and stop correcting
         if (PIPSisEQ(step_length_primal_enhanced, 1.) && PIPSisEQ(step_length_dual_enhanced, 1.)) {
            PIPSdebugMessage("both 1. \n");

            step->saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);
            step_length_primal = step_length_primal_enhanced;
            step_length_dual = step_length_dual_enhanced;
            std::cout << "Step 4: " << step_length_primal << ", " << step_length_dual << "\n";

            if (small_corr)
               NumberSmallCorrectors++;

            number_gondzio_corrections++;

            // exit Gondzio correction loop
            break;
         }
         else if (step_length_primal_enhanced >= (1. + acceptance_tolerance) * step_length_primal &&
                  step_length_dual_enhanced >= (1. + acceptance_tolerance) * step_length_dual) {
            PIPSdebugMessage("both better \n");

            // if enhanced step length is significantly better than the current alpha, make the enhanced step official, but maybe keep correcting
            step->saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);
            step_length_primal = step_length_primal_enhanced;
            step_length_dual = step_length_dual_enhanced;
            std::cout << "Step 5: " << step_length_primal << ", " << step_length_dual << "\n";

            if (small_corr)
               NumberSmallCorrectors++;

            number_gondzio_corrections++;
         }
         else if (step_length_primal_enhanced >= (1. + acceptance_tolerance) * step_length_primal) {
            PIPSdebugMessage("primal better \n");

            step->saxpy_pd(corrector_step, weight_primal_candidate, 0.);

            step_length_primal = step_length_primal_enhanced;
            std::cout << "Step 6: " << step_length_primal << ", " << step_length_dual << "\n";

            if (small_corr)
               NumberSmallCorrectors++;

            number_gondzio_corrections++;
         }
         else if (step_length_dual_enhanced >= (1. + acceptance_tolerance) * step_length_dual) {
            PIPSdebugMessage("dual better \n");

            step->saxpy_pd(corrector_step, 0., weight_dual_candidate);

            step_length_dual = step_length_dual_enhanced;
            std::cout << "Step 7: " << step_length_primal << ", " << step_length_dual << "\n";

            if (small_corr)
               NumberSmallCorrectors++;

            number_gondzio_corrections++;
         }
            /* if not done yet because correctors were not good enough - try a small corrector if enabled */
         else if (additional_correctors_small_comp_pairs && !small_corr && iteration >= first_iter_small_correctors) {
            if (step_length_primal < max_alpha_small_correctors || step_length_dual < max_alpha_small_correctors) {
               // try and center small pairs
               small_corr = true;
               if (my_rank == 0) {
                  std::cout << "Switching to small corrector " << std::endl;
                  std::cout << "Alpha when switching: " << step_length_primal << " " << step_length_dual << std::endl;
               }
            }
            else
               // exit Gondzio correction loop
               break;
         }
         else {
            // exit Gondzio correction loop
            break;
         }
      }

      // We've finally decided on a step direction, now calculate the length using Mehrotra's heuristic.x
      if (step_length_type == PRIMAL) {
         step_length_primal = mehrotra_step_length(&iterate, step);
         step_length_dual = step_length_primal;
      }
      else {
         mehrotra_step_length_PD(&iterate, step, step_length_primal, step_length_dual);
      }

      std::cout << "Step 8: " << step_length_primal << ", " << step_length_dual << "\n";
      assert(step_length_primal != 0);
      assert(step_length_dual != 0);

      // if we encountered numerical troubles while computing the step check enter a probing round
      if (numerical_troubles) {
         if (precond_decreased)
            precond_decreased = decreasePreconditionerImpact(linear_system);

         doProbing_pd(&problem, &iterate, &residuals, step_length_primal, step_length_dual);
         std::cout << "Step 9: " << step_length_primal << ", " << step_length_dual << "\n";
         const double alpha_max = std::max(step_length_primal, step_length_dual);

         if (restartIterateBecauseOfPoorStep(pure_centering_step, precond_decreased, alpha_max))
            continue;
      }

      // actually take the step and calculate the new mu
      iterate.saxpy_pd(step, step_length_primal, step_length_dual);
      mu = iterate.mu();

      pure_centering_step = false;
      numerical_troubles = false;

      distributed_factory->iterate_ended();
   }

   residuals.evaluate(problem, iterate);
   if (gOoqpPrintLevel >= 10) {
      this->do_monitor_Pd(&problem, &iterate, &residuals, step_length_primal, step_length_dual, sigma, iteration, mu, status_code, 1);
   }

   return status_code;
}

void InteriorPointMethod::compute_predictor_step(Problem& problem, Variables& iterate, Residuals& residuals) {
   residuals.set_complementarity_residual(iterate, 0.);
   linear_system->factorize(&problem, &iterate);
   linear_system->solve(&problem, &iterate, &residuals, step);
   step->negate();
}

void InteriorPointMethod::compute_corrector_step(Problem& problem, Variables& iterate, double sigma, double mu) {
   corrector_residuals->clear_linear_residuals();
   // form right hand side of linear system:
   corrector_residuals->set_complementarity_residual(*step, -sigma * mu);

   linear_system->solve(&problem, &iterate, corrector_residuals, corrector_step);
   corrector_step->negate();
}

void InteriorPointMethod::compute_gondzio_corrector(Problem* problem, Variables* iterate, double rmin, double rmax, bool small_corr) {
   // place XZ into the r3 component of corrector_residuals
   corrector_residuals->set_complementarity_residual(*corrector_step, 0.);
   if (small_corr)
      assert(additional_correctors_small_comp_pairs);
   // do the projection operation
   if (small_corr)
      corrector_residuals->project_r3(rmin, std::numeric_limits<double>::infinity());
   else
      corrector_residuals->project_r3(rmin, rmax);

   // solve for corrector direction
   linear_system->solve(problem, iterate, corrector_residuals, corrector_step); // corrector_step is now delta_m
}

void InteriorPointMethod::check_linear_system_solve_numerical_troubles(Residuals* residuals, bool& numerical_troubles, bool& small_corr) const {
   if (!bicgstab_converged && bigcstab_norm_res_rel * 1e2 * dnorm > residuals->residualNorm()) {
      PIPSdebugMessage("Step computation in BiCGStab failed");
      numerical_troubles = true;
      if (additional_correctors_small_comp_pairs && !small_corr) {
         if (PIPS_MPIgetRank() == 0)
            std::cout << "switching to small correctors\n";
         small_corr = true;
      }
   }
}

void InteriorPointMethod::registerBiCGStabOvserver(AbstractLinearSystem* sys) {
   /* every linsys handed to the GondzioStoch should be observable */
   assert(dynamic_cast<Subject*>(sys));
   setSubject(dynamic_cast<Subject*>(sys));
}

void InteriorPointMethod::setBiCGStabTol(int iteration) const {
   if (!dynamic_bicg_tol)
      return;

   assert(iteration >= -1);

   if (iteration == -1)
      pips_options::set_double_parameter("OUTER_BICG_TOL", 1e-10);
   else if (iteration <= 4)
      pips_options::set_double_parameter("OUTER_BICG_TOL", 1e-8);
   else if (iteration <= 8)
      pips_options::set_double_parameter("OUTER_BICG_TOL", 1e-9);
   else
      pips_options::set_double_parameter("OUTER_BICG_TOL", 1e-10);
}

void InteriorPointMethod::adjustLimitGondzioCorrectors() {
   assert(bicg_iterations >= 0);
   if (dynamic_corrector_schedule) {
      if (bicgstab_skipped)
         maximum_correctors = 5;
      else if (bicg_iterations < 2)
         maximum_correctors = 4;
      else if (bicg_iterations <= 15)
         maximum_correctors = 3;
      else if (bicg_iterations < 25)
         maximum_correctors = 2;
      else if (bicg_iterations > 35)
         maximum_correctors = 1;
   }
}

bool InteriorPointMethod::decreasePreconditionerImpact(AbstractLinearSystem* sys) const {
   bool success = false;
   dynamic_cast<DistributedRootLinearSystem*>(sys)->precondSC.decreaseDiagDomBound(success);
   if (!success) {
      if (PIPS_MPIgetRank() == 0)
         std::cout << "Cannot increase precision in preconditioner anymore\n";
   }
   return success;
}

void InteriorPointMethod::computeProbingStep_pd(Variables* probing_step, const Variables* iterate, const Variables* step, double alpha_primal,
      double alpha_dual) const {
   probing_step->copy(iterate);
   probing_step->saxpy_pd(step, alpha_primal, alpha_dual);
}

void InteriorPointMethod::doProbing_pd(Problem* prob, Variables* iterate, Residuals* resid, double& alpha_pri, double& alpha_dual) {
   const double mu_last = iterate->mu();
   const double resids_norm_last = resid->residualNorm();

   computeProbingStep_pd(temp_step, iterate, step, alpha_pri, alpha_dual);

   resid->evaluate(*prob, *temp_step, false);
   const double mu_probing = temp_step->mu();
   const double resids_norm_probing = resid->residualNorm();

   const double factor = computeStepFactorProbing(resids_norm_last, resids_norm_probing, mu_last, mu_probing);

   alpha_pri = factor * alpha_pri;
   alpha_dual = factor * alpha_dual;
}

/* when numerical troubles occurred we only allow controlled steps that worsen the residuals and mu by at most a factor of 10 */
double InteriorPointMethod::computeStepFactorProbing(double resids_norm_last, double resids_norm_probing, double mu_last, double mu_probing) const {
   assert(resids_norm_last > 0.);
   assert(resids_norm_probing > 0.);

   double factor = 1.;
   const double limit_resids = 10. * resids_norm_last;

   if (resids_norm_probing > limit_resids) {
      assert(resids_norm_probing > resids_norm_last);
      assert(limit_resids > resids_norm_last);

      const double resids_diff = resids_norm_probing - resids_norm_last;
      const double resids_max_change = limit_resids - resids_norm_last;

      assert(resids_diff > 0.);
      assert(resids_max_change > 0.);
      assert(resids_max_change < resids_diff);

      factor = std::min(factor, resids_max_change / resids_diff * 0.9995);
   }

   const double mu_limit = 10. * mu_last;
   if (mu_probing > mu_limit) {
      assert(mu_probing > mu_last);
      assert(mu_limit > mu_last);
      const double mu_diff = mu_probing - mu_last;
      const double mu_max_change = mu_limit - mu_last;
      assert(mu_diff > 0.);
      assert(mu_max_change > 0.);
      assert(mu_max_change < mu_diff);

      factor = std::min(factor, mu_max_change / mu_diff * 0.9995);
   }

   assert(0. < factor);
   assert(factor <= 1.);

   return factor;
}

bool InteriorPointMethod::restartIterateBecauseOfPoorStep(bool& pure_centering_step, bool precond_decreased, double alpha_max) const {
   const int my_rank = PIPS_MPIgetRank();

   if (!pure_centering_step && alpha_max < mutol * 1e-2) {
      if (my_rank == 0)
         std::cout << "poor step computed - trying pure centering step\n";
      pure_centering_step = true;
      return true;
   }
   else if (alpha_max < mutol * 1e-2 && pure_centering_step && precond_decreased) {
      if (my_rank == 0)
         std::cout << "poor step computed - refactorization with decreased preconditioning\n";
      return true;
   }

   if (my_rank == 0)
      std::cout << "poor step computed but keeping it\n";
   return false;
}

InteriorPointMethod::~InteriorPointMethod() {
   delete temp_step;
}

void InteriorPointMethod::notifyFromSubject() {
   const Subject& subj = *getSubject();

   bicgstab_skipped = subj.getBoolValue("BICG_SKIPPED");
   if (!bicgstab_skipped)
      bicgstab_converged = subj.getBoolValue("BICG_CONVERGED");
   else
      bicgstab_converged = true;
   bigcstab_norm_res_rel = subj.getDoubleValue("BICG_RELRESNORM");
   bicg_iterations = subj.getIntValue("BICG_NITERATIONS");
   if (!bicgstab_converged)
      PIPSdebugMessage("BiGCStab had troubles converging\n");
}

void InteriorPointMethod::calculate_alpha_PD_weight_candidate(Variables* iterate, Variables* predictor_step, Variables* corrector_step,
      double alpha_primal, double alpha_dual, double& step_length_primal_candidate, double& step_length_dual_candidate,
      double& weight_primal_candidate, double& weight_dual_candidate) {
   assert(alpha_primal > 0. && alpha_primal <= 1.);
   assert(alpha_dual > 0. && alpha_dual <= 1.);

   double alpha_primal_best = -1., alpha_dual_best = -1.;
   double weight_primal_best = -1., weight_dual_best = -1.;
   const double weight_min = alpha_primal * alpha_dual;
   const double weight_interval_length = 1. - weight_min;

   // main loop
   for (unsigned int i = 0; i <= n_linesearch_points; i++) {
      double weight_current = weight_min + (weight_interval_length / (n_linesearch_points)) * i;

      weight_current = std::min(weight_current, 1.);
      assert(weight_current > 0.);

      temp_step->copy(predictor_step);
      temp_step->saxpy(corrector_step, weight_current);

      double alpha_primal_current = 1., alpha_dual_current = 1.;

      // compute the step length
      if (step_length_type == PRIMAL) {
         double step_length = iterate->stepbound(step);
         alpha_primal_current = step_length;
         alpha_dual_current = step_length;
      }
      else {
         iterate->stepbound_pd(temp_step, alpha_primal_current, alpha_dual_current);
      }

      assert(alpha_primal_current > 0. && alpha_primal_current <= 1.);
      assert(alpha_dual_current > 0. && alpha_dual_current <= 1.);

      // find out the maximal alpha
      if (alpha_primal_current > alpha_primal_best) {
         alpha_primal_best = alpha_primal_current;
         weight_primal_best = weight_current;
      }
      if (alpha_dual_current > alpha_dual_best) {
         alpha_dual_best = alpha_dual_current;
         weight_dual_best = weight_current;
      }
   }

   assert(alpha_primal_best >= 0. && weight_primal_best >= 0.);
   assert(alpha_dual_best >= 0. && weight_dual_best >= 0.);

   weight_primal_candidate = weight_primal_best;
   weight_dual_candidate = weight_dual_best;

   step_length_primal_candidate = alpha_primal_best;
   step_length_dual_candidate = alpha_dual_best;
}