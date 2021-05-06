//
// Created by charlie on 26.04.21.
//

#include <cassert>
#include <Observer.h>
#include <DistributedOptions.h>
#include <Monitor.h>
#include "Status.h"
#include "MehrotraHeuristic.hpp"
#include "Problem.h"
#include "Variables.h"
#include "Residuals.h"
#include "DistributedFactory.h"
#include "DistributedRootLinearSystem.h"

extern int gOoqpPrintLevel;
extern double g_iterNumber;

int gLackOfAccuracy = 0;
const unsigned int max_linesearch_points = 50;

MehrotraHeuristic::MehrotraHeuristic(DistributedFactory& factory, Problem& problem, const Scaler* scaler = nullptr) :
      scaler(scaler), corrector_step(factory.make_variables(problem)), corrector_residuals(factory.make_residuals(problem)),
      n_linesearch_points(pips_options::get_int_parameter("GONDZIO_STOCH_N_LINESEARCH")),
      temp_step(factory.make_variables(problem)),
      bicgstab_skipped(false), bicgstab_converged(true), bigcstab_norm_res_rel(0.), bicg_iterations(0),
      dynamic_corrector_schedule(pips_options::get_bool_parameter("GONDZIO_STOCH_USE_DYNAMIC_CORRECTOR_SCHEDULE")),
      additional_correctors_small_comp_pairs(pips_options::get_bool_parameter("GONDZIO_STOCH_ADDITIONAL_CORRECTORS_SMALL_VARS")),
      max_additional_correctors(pips_options::get_int_parameter("GONDZIO_STOCH_ADDITIONAL_CORRECTORS_MAX")),
      first_iter_small_correctors(pips_options::get_int_parameter("GONDZIO_STOCH_FIRST_ITER_SMALL_CORRECTORS")),
      max_alpha_small_correctors(pips_options::get_double_parameter("GONDZIO_STOCH_MAX_ALPHA_SMALL_CORRECTORS")),
      NumberSmallCorrectors(0),
      maximum_correctors(qpgen_options::getIntParameter("GONDZIO_MAX_CORRECTORS")),
      number_gondzio_corrections(0),
      step_factor0(0.3), step_factor1(1.5), acceptance_tolerance(0.01), beta_min(0.1), beta_max(10),
      dynamic_bicg_tol(pips_options::get_bool_parameter("OUTER_BICG_DYNAMIC_TOL")),
      tsig(3.), max_iterations(300) {
   assert(max_additional_correctors > 0);
   assert(first_iter_small_correctors >= 0);
   assert(0 < max_alpha_small_correctors && max_alpha_small_correctors < 1);
   assert(n_linesearch_points > 0);

   if (base_options::getBoolParameter("IP_ACCURACY_REDUCED")) {
      artol = 1.e-3;
      mutol = 1.e-5;
   }

   if (pips_options::get_bool_parameter("GONDZIO_STOCH_ADAPTIVE_LINESEARCH")) {
      const int size = PIPS_MPIgetSize();

      if (size > 1)
         this->n_linesearch_points = std::min(unsigned(size) + this->n_linesearch_points, max_linesearch_points);
   }
   
   if (base_options::getBoolParameter("IP_STEPLENGTH_CONSERVATIVE")) {
      steplength_factor = 0.99;
      gamma_f = 0.95;
   }
   gamma_a = 1. / (1. - gamma_f);

   if (base_options::getBoolParameter("IP_PRINT_TIMESTAMP")) {
      print_timestamp = true;
      start_time = MPI_Wtime();
   }

   // allocate space to track the sequence problem_formulation complementarity gaps, residual norms, and merit functions.
   mu_history = new double[max_iterations];
   residual_norm_history = new double[max_iterations];
   phi_history = new double[max_iterations];
   phi_min_history = new double[max_iterations];
}

TerminationCode
MehrotraHeuristic::corrector_predictor(DistributedFactory& factory, Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
      AbstractLinearSystem& linear_system) {
   set_problem_norm(problem);
   // register as observer for the BiCGStab solves
   registerBiCGStabOvserver(&linear_system);
   set_BiCGStab_tolerance(-1);

   int iteration = 0;
   int number_gondzio_corrections = 0;
   double mu = iterate.mu();
   double alpha = 1., sigma = 1.;
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   double mu_affine;
   double alpha_target, alpha_enhanced;
   TerminationCode status_code;

   g_iterNumber = 0.;

   bool pure_centering_step = false;
   bool numerical_troubles = false;
   bool precond_decreased = true;
   while (true) {
      iteration++;

      set_BiCGStab_tolerance(iteration);
      bool small_corr = false;

      factory.iterate_started();

// evaluate residuals and update algorithm status:
      residuals.evaluate(problem, iterate);

//  termination test:
      status_code = this->default_status(&problem, &iterate, &residuals, dnorm, iteration, mu, SUCCESSFUL_TERMINATION);

      if (status_code != NOT_FINISHED)
         break;

      if (gOoqpPrintLevel >= 10) {
         this->do_monitor(&problem, &iterate, &residuals, dnorm, alpha, sigma, iteration, mu, status_code, 0);
      }

// *** Predictor step ***
      if (!pure_centering_step) {
         compute_predictor_step(problem, iterate, residuals, linear_system, step);
         checkLinsysSolveNumericalTroublesAndReact(&residuals, numerical_troubles, small_corr);
      }
      else
         step.setToZero();

      alpha = iterate.stepbound(&step);

// calculate centering parameter
      mu_affine = iterate.mustep_pd(&step, alpha, alpha);

      assert(!PIPSisZero(mu));
      sigma = pow(mu_affine / mu, tsig);

      if (gOoqpPrintLevel >= 10) {
         this->do_monitor(&problem, &iterate, &residuals, dnorm, alpha, sigma, iteration, mu, status_code, 2);
      }

      g_iterNumber += 1.;

      compute_corrector_step(problem, iterate, linear_system, step, sigma, mu);
      checkLinsysSolveNumericalTroublesAndReact(&residuals, numerical_troubles, small_corr);

// calculate weighted predictor-corrector step
      double weight_candidate = -1.;
      const double alpha_predictor = alpha;

      calculate_alpha_weight_candidate(&iterate, &step, corrector_step, alpha_predictor, alpha, weight_candidate);

      assert(weight_candidate >= 0. && weight_candidate <= 1.);

      step.saxpy(corrector_step, weight_candidate);

// prepare for Gondzio corrector loop: zero out the
// corrector_residuals structure:
      corrector_residuals->clear_linear_residuals();

// calculate the target box:
      const double rmin = sigma * mu * beta_min;
      const double rmax = sigma * mu * beta_max;

      number_gondzio_corrections = 0;
      NumberSmallCorrectors = 0;

// if small_corr_aggr only try small correctors

// enter the Gondzio correction loop:
      while (number_gondzio_corrections < maximum_correctors && NumberSmallCorrectors < max_additional_correctors && PIPSisLT(alpha, 1.)) {
         if (dynamic_corrector_schedule)
            adjustLimitGondzioCorrectors();
         corrector_step->copy(&iterate);

// calculate target steplength
         alpha_target = step_factor1 * alpha + step_factor0;
         if (alpha_target > 1.)
            alpha_target = 1.;

// add a step of this length to corrector_step
         corrector_step->saxpy(&step, alpha_target);
// corrector_step is now x_k + alpha_target * delta_p (a trial point)

/* compute corrector step */
         compute_gondzio_corrector(problem, iterate, linear_system, rmin, rmax, small_corr);
         const bool was_small_corr = small_corr;
         checkLinsysSolveNumericalTroublesAndReact(&residuals, numerical_troubles, small_corr);

         if (numerical_troubles) {
            if (!was_small_corr && small_corr)
               continue;
            else
// exit corrector loop if small correctors have already been tried or are not allowed
               break;
         }

// calculate weighted predictor-corrector step
         calculate_alpha_weight_candidate(&iterate, &step, corrector_step, alpha_target, alpha_enhanced, weight_candidate);

// if the enhanced step length is actually 1, make it official
// and stop correcting
         if (PIPSisEQ(alpha_enhanced, 1.)) {
            step.saxpy(corrector_step, weight_candidate);
            alpha = alpha_enhanced;

            if (small_corr)
               NumberSmallCorrectors++;

            number_gondzio_corrections++;

// exit Gondzio correction loop
            break;
         }
         else if (alpha_enhanced >= (1. + acceptance_tolerance) * alpha) {
// if enhanced step length is significantly better than the
// current alpha, make the enhanced step official, but maybe
// keep correcting
            step.saxpy(corrector_step, weight_candidate);
            alpha = alpha_enhanced;

            if (small_corr)
               NumberSmallCorrectors++;

            number_gondzio_corrections++;
         }
/* if not done yet because correctors were not good enough - try a small corrector if enabled */
         else if (additional_correctors_small_comp_pairs && !small_corr && iteration >= first_iter_small_correctors) {
            if (alpha < max_alpha_small_correctors) {
               small_corr = true;
               if (my_rank == 0) {
                  std::cout << "Switching to small corrector " << std::endl;
                  std::cout << "Alpha when switching: " << alpha << std::endl;
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

// We've finally decided on a step direction, now calculate the
// length using Mehrotra's heuristic.x
      alpha = mehrotra_step_length(&iterate, &step);
      assert(alpha != 0);

// alternatively, just use a crude step scaling factor.
// alpha = 0.995 * iterate->stepbound( step );

// if we encountered numerical troubles while computing the step enter a probing round
      if (numerical_troubles) {
         if (precond_decreased)
            precond_decreased = decreasePreconditionerImpact(&linear_system);
         doProbing(&problem, &iterate, &residuals, &step, alpha);
         if (restartIterateBecauseOfPoorStep(pure_centering_step, precond_decreased, alpha))
            continue;
      }

// actually take the step (at last!) and calculate the new mu
      iterate.saxpy(&step, alpha);
      mu = iterate.mu();

      pure_centering_step = false;
      numerical_troubles = false;

      factory.iterate_ended();
   }
   residuals.evaluate(problem, iterate);
   if (gOoqpPrintLevel >= 10) {
      this->do_monitor(&problem, &iterate, &residuals, dnorm, alpha, sigma, iteration, mu, status_code, 1);
   }
   return status_code;
}

TerminationCode
MehrotraHeuristic::corrector_predictor_pd(DistributedFactory& factory, Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
      AbstractLinearSystem& linear_system) {
   set_problem_norm(problem);
   // register as observer for the BiCGStab solves
   registerBiCGStabOvserver(&linear_system);
   set_BiCGStab_tolerance(-1);

   int iteration = 0;
   int number_gondzio_corrections = 0;
   double mu = iterate.mu();
   double alpha_primal = 1., alpha_dual = 1., sigma = 1.;
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   double alpha_target, alpha_enhanced;
   TerminationCode status_code;

   g_iterNumber = 0.;

   bool pure_centering_step = false;
   bool numerical_troubles = false;
   bool precond_decreased = true;

   while (true) {
      iteration++;

      set_BiCGStab_tolerance(iteration);
      bool small_corr = false;

      factory.iterate_started();

      // evaluate residuals and update algorithm status:
      residuals.evaluate(problem, iterate);

      //  termination test:
      status_code = this->default_status(&problem, &iterate, &residuals, dnorm, iteration, mu, SUCCESSFUL_TERMINATION);

      if (status_code != NOT_FINISHED)
         break;

      if (gOoqpPrintLevel >= 10) {
         this->do_monitor_Pd(&problem, &iterate, &residuals, dnorm, alpha_primal, alpha_dual, sigma, iteration, mu, status_code, 0);
      }

      // *** Predictor step ***
      if (!pure_centering_step) {
         compute_predictor_step(problem, iterate, residuals, linear_system, step);
         checkLinsysSolveNumericalTroublesAndReact(&residuals, numerical_troubles, small_corr);
      }
      else
         step.setToZero();

      iterate.stepbound_pd(&step, alpha_primal, alpha_dual);

      // calculate centering parameter
      double mu_affine = iterate.mustep_pd(&step, alpha_primal, alpha_dual);

      assert(!PIPSisZero(mu));
      sigma = pow(mu_affine / mu, tsig);

      if (gOoqpPrintLevel >= 10) {
         this->do_monitor_Pd(&problem, &iterate, &residuals, dnorm, alpha_primal, alpha_dual, sigma, iteration, mu, status_code, 2);
      }

      g_iterNumber += 1.;

      // *** Corrector step ***
      compute_corrector_step(problem, iterate, linear_system, step, sigma, mu);
      checkLinsysSolveNumericalTroublesAndReact(&residuals, numerical_troubles, small_corr);

      // calculate weighted predictor-corrector step
      double weight_primal_candidate, weight_dual_candidate = -1.;

      calculateAlphaPDWeightCandidate(&iterate, &step, corrector_step, alpha_primal, alpha_dual, alpha_primal, alpha_dual, weight_primal_candidate,
            weight_dual_candidate);

      assert(weight_primal_candidate >= 0. && weight_primal_candidate <= 1.);
      assert(weight_dual_candidate >= 0. && weight_dual_candidate <= 1.);

      step.saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

      // prepare for Gondzio corrector loop: zero out the corrector_residuals structure:
      corrector_residuals->clear_linear_residuals();

      // calculate the target box:
      const double rmin = sigma * mu * beta_min;
      const double rmax = sigma * mu * beta_max;

      number_gondzio_corrections = 0;
      NumberSmallCorrectors = 0;

      // enter the Gondzio correction loop:
      while (number_gondzio_corrections < maximum_correctors && NumberSmallCorrectors < max_additional_correctors &&
             (PIPSisLT(alpha_primal, 1.) || PIPSisLT(alpha_dual, 1.))) {
         if (dynamic_corrector_schedule)
            adjustLimitGondzioCorrectors();
         corrector_step->copy(&iterate);

         double alpha_primal_enhanced, alpha_dual_enhanced;
         const double alpha_p_target = std::min(1., step_factor1 * alpha_primal + step_factor0);
         const double alpha_dual_target = std::min(1., step_factor1 * alpha_dual + step_factor0);

         PIPSdebugMessage("corrector loop: %d alpha_primal: %f alpha_dual %f \n", number_gondzio_corrections, alpha_primal, alpha_dual);

         // add a step of this length to corrector_step
         corrector_step->saxpy_pd(&step, alpha_p_target, alpha_dual_target);
         // corrector_step is now x_k + alpha_target * delta_p (a trial point)

         /* compute corrector step */
         compute_gondzio_corrector(problem, iterate, linear_system, rmin, rmax, small_corr);
         const bool was_small_corr = small_corr;
         checkLinsysSolveNumericalTroublesAndReact(&residuals, numerical_troubles, small_corr);

         if (numerical_troubles) {
            if (!was_small_corr && small_corr)
               continue;
            else
               // exit corrector loop if small correctors have already been tried or are not allowed
               break;
         }

         // calculate weighted predictor-corrector step
         calculateAlphaPDWeightCandidate(&iterate, &step, corrector_step, alpha_p_target, alpha_dual_target, alpha_primal_enhanced,
               alpha_dual_enhanced,weight_primal_candidate, weight_dual_candidate);

         // if the enhanced step length is actually 1, make it official
         // and stop correcting
         if (PIPSisEQ(alpha_primal_enhanced, 1.) && PIPSisEQ(alpha_dual_enhanced, 1.)) {
            PIPSdebugMessage("both 1. \n");

            step.saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

            alpha_primal = alpha_primal_enhanced;
            alpha_dual = alpha_dual_enhanced;

            if (small_corr)
               NumberSmallCorrectors++;

            number_gondzio_corrections++;

            // exit Gondzio correction loop
            break;
         }
         else if (alpha_primal_enhanced >= (1. + acceptance_tolerance) * alpha_primal &&
                  alpha_dual_enhanced >= (1. + acceptance_tolerance) * alpha_dual) {
            PIPSdebugMessage("both better \n");

            // if enhanced step length is significantly better than the
            // current alpha, make the enhanced step official, but maybe
            // keep correcting
            step.saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

            alpha_primal = alpha_primal_enhanced;
            alpha_dual = alpha_dual_enhanced;

            if (small_corr)
               NumberSmallCorrectors++;

            number_gondzio_corrections++;
         }
         else if (alpha_primal_enhanced >= (1. + acceptance_tolerance) * alpha_primal) {
            PIPSdebugMessage("primal better \n");

            step.saxpy_pd(corrector_step, weight_primal_candidate, 0.);

            alpha_primal = alpha_primal_enhanced;

            if (small_corr)
               NumberSmallCorrectors++;

            number_gondzio_corrections++;
         }
         else if (alpha_dual_enhanced >= (1. + acceptance_tolerance) * alpha_dual) {
            PIPSdebugMessage("dual better \n");

            step.saxpy_pd(corrector_step, 0., weight_dual_candidate);

            alpha_dual = alpha_dual_enhanced;

            if (small_corr)
               NumberSmallCorrectors++;

            number_gondzio_corrections++;
         }
            /* if not done yet because correctors were not good enough - try a small corrector if enabled */
         else if (additional_correctors_small_comp_pairs && !small_corr && iteration >= first_iter_small_correctors) {
            if (alpha_primal < max_alpha_small_correctors || alpha_dual < max_alpha_small_correctors) {
               // try and center small pairs
               small_corr = true;
               if (my_rank == 0) {
                  std::cout << "Switching to small corrector " << std::endl;
                  std::cout << "Alpha when switching: " << alpha_primal << " " << alpha_dual << std::endl;
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

      // We've finally decided on a step direction, now calculate the
      // length using Mehrotra's heuristic.x
      mehrotra_step_length_PD(&iterate, &step, alpha_primal, alpha_dual);

      // if we encountered numerical troubles while computing the step check enter a probing round
      if (numerical_troubles) {
         if (precond_decreased)
            precond_decreased = decreasePreconditionerImpact(&linear_system);

         doProbing_pd(&problem, &iterate, &residuals, &step, alpha_primal, alpha_dual);
         const double alpha_max = std::max(alpha_primal, alpha_dual);

         if (restartIterateBecauseOfPoorStep(pure_centering_step, precond_decreased, alpha_max))
            continue;
      }

      // actually take the step and calculate the new mu
      iterate.saxpy_pd(&step, alpha_primal, alpha_dual);
      mu = iterate.mu();

      pure_centering_step = false;
      numerical_troubles = false;

      factory.iterate_ended();
   }
   residuals.evaluate(problem, iterate);
   if (gOoqpPrintLevel >= 10) {
      this->do_monitor_Pd(&problem, &iterate, &residuals, dnorm, alpha_primal, alpha_dual, sigma, iteration, mu, status_code, 1);
   }
   return status_code;
}

void MehrotraHeuristic::compute_predictor_step(Problem& problem, Variables& iterate, Residuals& residuals, AbstractLinearSystem& linear_system,
      Variables& step) {
   residuals.set_complementarity_residual(iterate, 0.);
   linear_system.factorize(&problem, &iterate);
   linear_system.solve(&problem, &iterate, &residuals, &step);
   step.negate();
}

void
MehrotraHeuristic::compute_corrector_step(Problem& problem, Variables& iterate, AbstractLinearSystem& linear_system, Variables& step, double sigma,
      double mu) {
   corrector_residuals->clear_linear_residuals();
   // form right hand side of linear system:
   corrector_residuals->set_complementarity_residual(step, -sigma * mu);

   linear_system.solve(&problem, &iterate, corrector_residuals, corrector_step);
   corrector_step->negate();
}

void MehrotraHeuristic::compute_gondzio_corrector(Problem& problem, Variables& iterate, AbstractLinearSystem& linear_system, double rmin, double rmax,
      bool small_corr) {
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
   linear_system.solve(&problem, &iterate, corrector_residuals, corrector_step); // corrector_step is now delta_m
}

void
MehrotraHeuristic::calculate_alpha_weight_candidate(Variables* iterate, Variables* predictor_step, Variables* corrector_step, double alpha_predictor,
      double& alpha_candidate, double& weight_candidate) {
   assert(corrector_step);
   assert(predictor_step);
   assert(temp_step);

   assert(alpha_predictor > 0. && alpha_predictor <= 1.);

   double alpha_best = -1.;
   double weight_best = -1.;
   const double weight_min = alpha_predictor * alpha_predictor;
   const double weight_interval_length = 1. - weight_min;

   // main loop
   for (unsigned int n = 0; n <= n_linesearch_points; n++) {
      double weight_curr = weight_min + (weight_interval_length / (n_linesearch_points)) * n;

      weight_curr = std::min(weight_curr, 1.);

      assert(weight_curr > 0. && weight_curr <= 1.);

      temp_step->copy(predictor_step);
      temp_step->saxpy(corrector_step, weight_curr);

      const double alpha_curr = iterate->stepbound(temp_step);
      assert(alpha_curr > 0. && alpha_curr <= 1.);

      if (alpha_curr > alpha_best) {
         alpha_best = alpha_curr;
         weight_best = weight_curr;
      }
   }

   assert(alpha_best >= 0. && weight_best >= 0.);

   weight_candidate = weight_best;
   alpha_candidate = alpha_best;
}

void MehrotraHeuristic::calculateAlphaPDWeightCandidate(Variables* iterate, Variables* predictor_step, Variables* corrector_step, double alpha_primal,
      double alpha_dual, double& alpha_primal_candidate, double& alpha_dual_candidate, double& weight_primal_candidate,
      double& weight_dual_candidate) {
   assert(alpha_primal > 0. && alpha_primal <= 1.);
   assert(alpha_dual > 0. && alpha_dual <= 1.);

   double alpha_primal_best = -1., alpha_dual_best = -1.;
   double weight_primal_best = -1., weight_dual_best = -1.;
   const double weight_min = alpha_primal * alpha_dual;
   const double weight_intervallength = 1. - weight_min;

   // main loop
   for (unsigned int n = 0; n <= n_linesearch_points; n++) {
      double weight_curr = weight_min + (weight_intervallength / (n_linesearch_points)) * n;

      weight_curr = std::min(weight_curr, 1.);

      assert(weight_curr > 0. && weight_curr <= 1.);

      temp_step->copy(predictor_step);
      temp_step->saxpy(corrector_step, weight_curr);

      double alpha_primal_curr = 1., alpha_dual_curr = 1.;
      iterate->stepbound_pd(temp_step, alpha_primal_curr, alpha_dual_curr);
      assert(alpha_primal_curr > 0. && alpha_primal_curr <= 1.);
      assert(alpha_dual_curr > 0. && alpha_dual_curr <= 1.);

      if (alpha_primal_curr > alpha_primal_best) {
         alpha_primal_best = alpha_primal_curr;
         weight_primal_best = weight_curr;
      }
      if (alpha_dual_curr > alpha_dual_best) {
         alpha_dual_best = alpha_dual_curr;
         weight_dual_best = weight_curr;
      }
   }

   assert(alpha_primal_best >= 0. && weight_primal_best >= 0.);
   assert(alpha_dual_best >= 0. && weight_dual_best >= 0.);

   weight_primal_candidate = weight_primal_best;
   weight_dual_candidate = weight_dual_best;

   alpha_primal_candidate = alpha_primal_best;
   alpha_dual_candidate = alpha_dual_best;
}

void MehrotraHeuristic::doProbing(Problem* problem, Variables* iterate, Residuals* residuals, Variables* step, double& alpha) {
   const double mu_last = iterate->mu();
   const double resids_norm_last = residuals->residualNorm();

   computeProbingStep(temp_step, iterate, step, alpha);

   residuals->evaluate(*problem, *temp_step, false);
   const double mu_probing = temp_step->mu();
   const double resids_norm_probing = residuals->residualNorm();

   const double factor = compute_step_factor_probing(resids_norm_last, resids_norm_probing, mu_last, mu_probing);

   alpha = factor * alpha;
}

void MehrotraHeuristic::doProbing_pd(Problem* problem, Variables* iterate, Residuals* residuals, Variables* step, double& alpha_primal, double& alpha_dual) {
   const double mu_last = iterate->mu();
   const double resids_norm_last = residuals->residualNorm();

   computeProbingStep_pd(temp_step, iterate, step, alpha_primal, alpha_dual);

   residuals->evaluate(*problem, *temp_step, false);
   const double mu_probing = temp_step->mu();
   const double resids_norm_probing = residuals->residualNorm();

   const double factor = compute_step_factor_probing(resids_norm_last, resids_norm_probing, mu_last, mu_probing);

   alpha_primal = factor * alpha_primal;
   alpha_dual = factor * alpha_dual;
}

bool MehrotraHeuristic::restartIterateBecauseOfPoorStep(bool& pure_centering_step, bool precond_decreased, double alpha_max) const {
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

void MehrotraHeuristic::computeProbingStep(Variables* probing_step, const Variables* iterate, const Variables* step, double alpha) const {
   probing_step->copy(iterate);
   probing_step->saxpy(step, alpha);
}

void MehrotraHeuristic::computeProbingStep_pd(Variables* probing_step, const Variables* iterate, const Variables* step, double alpha_primal, double
alpha_dual) const {
   probing_step->copy(iterate);
   probing_step->saxpy_pd(step, alpha_primal, alpha_dual);
}

/* when numerical troubles occurred we only allow controlled steps that worsen the residuals and mu by at most a factor of 10 */
double MehrotraHeuristic::compute_step_factor_probing(double resids_norm_last, double resids_norm_probing, double mu_last, double mu_probing) const {
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

bool MehrotraHeuristic::decreasePreconditionerImpact(AbstractLinearSystem* sys) const {
   bool success = false;
   dynamic_cast<DistributedRootLinearSystem*>(sys)->precondSC.decreaseDiagDomBound(success);
   if (!success) {
      if (PIPS_MPIgetRank() == 0)
         std::cout << "Cannot increase precision in preconditioner anymore\n";
   }
   return success;
}

void MehrotraHeuristic::adjustLimitGondzioCorrectors() {
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

void MehrotraHeuristic::set_BiCGStab_tolerance(int iteration) const {
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

void MehrotraHeuristic::notifyFromSubject() {
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

void MehrotraHeuristic::checkLinsysSolveNumericalTroublesAndReact(Residuals* residuals, bool& numerical_troubles, bool& small_corr) const {
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

void MehrotraHeuristic::add_monitor(Monitor* m) {
   // Push the monitor onto the list
   m->nextMonitor = itsMonitors;
   itsMonitors = m;
}

void
MehrotraHeuristic::do_monitor(const Problem* problem, const Variables* iterate, const Residuals* residuals, double dnorm, double alpha, double sigma,
      int i, double mu, int stop_code, int level) {
   Monitor* m = itsMonitors;

   while (m) {
      m->doIt(problem, iterate, residuals, dnorm, alpha, sigma, i, mu, stop_code, level);
      m = m->nextMonitor;
   }
}

void MehrotraHeuristic::do_monitor_Pd(const Problem* problem, const Variables* iterate, const Residuals* residuals, double dnorm, double alpha_primal,
      double alpha_dual, double sigma, int i, double mu, int stop_code, int level) {
   Monitor* m = itsMonitors;

   while (m) {
      m->doItPd(problem, iterate, residuals, dnorm, alpha_primal, alpha_dual, sigma, i, mu, stop_code, level);
      m = m->nextMonitor;
   }
}

TerminationCode
MehrotraHeuristic::default_status(const Problem* data, const Variables* iterate /* iterate */, const Residuals* residuals, double dnorm_orig,
      int iteration, double mu, TerminationCode level) {
   const int myrank = PIPS_MPIgetRank();
   TerminationCode stop_code = NOT_FINISHED;

   const std::pair<double, double> gap_norm = compute_unscaled_gap_and_residual_norm(*residuals);
   const double gap = gap_norm.first;
   const double rnorm = gap_norm.second;

   int index = std::min(std::max(0, iteration - 1), max_iterations - 1);

   // store the historical record
   mu_history[index] = mu;
   residual_norm_history[index] = rnorm;
   double phi = (rnorm + gap) / dnorm_orig;
   phi_history[index] = phi;

   if (index > 0) {
      phi_min_history[index] = phi_min_history[index - 1];
      if (phi < phi_min_history[index])
         phi_min_history[index] = phi;
   }
   else
      phi_min_history[index] = phi;

   if (iteration >= max_iterations)
      stop_code = MAX_ITS_EXCEEDED;
   else if (mu <= mutol && rnorm <= artol * dnorm_orig)
      stop_code = SUCCESSFUL_TERMINATION;

   if (myrank == 0) {
      std::cout << "mu/mutol: " << mu << "  " << mutol << "  ....   rnorm/limit: " << rnorm << " " << artol * dnorm_orig << std::endl;

      if (print_timestamp) {
         const double timestamp = MPI_Wtime() - start_time;
         std::cout << "time stamp: " << timestamp << std::endl;
      }
   }

   if (stop_code != NOT_FINISHED)
      return stop_code;

   // check infeasibility condition
   if (index >= 10 && phi >= 1.e-8 && phi >= 1.e4 * phi_min_history[index]) {
#ifdef TIMING
      if( myrank == 0 )
         std::cout << "possible INFEASIBLITY detected, phi: " << phi << std::endl;
#endif
      stop_code = INFEASIBLE;
   }

   if (stop_code != NOT_FINISHED)
      return stop_code;

   // check for unknown status: slow convergence first
   if (index >= 350 && phi_min_history[index] >= 0.5 * phi_min_history[index - 30]) {
      stop_code = UNKNOWN;
      printf("dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm_orig, artol);
   }

   if (index >= 350 && rnorm > artol * dnorm_orig &&
       residual_norm_history[index] * mu_history[0] >= 1.e8 * mu_history[index] * residual_norm_history[0]) {
      stop_code = UNKNOWN;
      printf("dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm_orig, artol);
   }

   if (mu * dnorm_orig < 1e5 * rnorm) {
      gLackOfAccuracy = 1;
   }
   else {
      gLackOfAccuracy = 1;
   }
   return stop_code;
}

double MehrotraHeuristic::mehrotra_step_length(Variables* iterate, Variables* step) {
   double primalValue = -std::numeric_limits<double>::max();
   double primalStep = -std::numeric_limits<double>::max();
   double dualValue = -std::numeric_limits<double>::max();
   double dualStep = -std::numeric_limits<double>::max();


#ifdef TIMING
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

   int firstOrSecond = -1;
   const double maximum_step_length = iterate->findBlocking(step, primalValue, primalStep, dualValue, dualStep, firstOrSecond);
   const double mu_full = iterate->mustep_pd(step, maximum_step_length, maximum_step_length) / gamma_a;

   double step_length = 1.;
   switch (firstOrSecond) {
      case 0:
         step_length = 1; // No constraints were blocking
         break;
      case 1:
         step_length = (-primalValue + mu_full / (dualValue + maximum_step_length * dualStep)) / primalStep;
#ifdef TIMING
         if( myrank == 0 )
            std::cout << "(primal) original alpha " << alpha << std::endl;
#endif
         break;
      case 2:
         step_length = (-dualValue + mu_full / (primalValue + maximum_step_length * primalStep)) / dualStep;
#ifdef TIMING
         if( myrank == 0 )
            std::cout << "(dual) original alpha " << alpha << std::endl;
#endif
         break;
      default:
         std::cout << "Can't get here: firstOrSecond=" << firstOrSecond << std::endl;
         assert(0 && "Can't get here");
         break;
   }
   // safeguard against numerical troubles in the above computations
   step_length = std::min(maximum_step_length, step_length);

   // make it at least gamma_f * maxStep
   if (step_length < gamma_f * maximum_step_length)
      step_length = gamma_f * maximum_step_length;

   // back off just a touch (or a bit more)
   step_length *= steplength_factor;

   assert(step_length < 1.);

   return step_length;
}

void MehrotraHeuristic::mehrotra_step_length_PD(Variables* iterate, Variables* step, double& alpha_primal, double& alpha_dual) {
   double primalValue_p = -std::numeric_limits<double>::max();
   double primalStep_p = -std::numeric_limits<double>::max();
   double dualValue_p = -std::numeric_limits<double>::max();
   double dualStep_p = -std::numeric_limits<double>::max();
   double maxAlpha_p;

   double primalValue_d = -std::numeric_limits<double>::max();
   double primalStep_d = -std::numeric_limits<double>::max();
   double dualValue_d = -std::numeric_limits<double>::max();
   double dualStep_d = -std::numeric_limits<double>::max();
   double maxAlpha_d;

   bool primalBlocking, dualBlocking;

   iterate->findBlocking_pd(step, primalValue_p, primalStep_p, dualValue_p, dualStep_p, primalValue_d, primalStep_d, dualValue_d, dualStep_d,
         maxAlpha_p, maxAlpha_d, primalBlocking, dualBlocking);

   const double mufull = iterate->mustep_pd(step, maxAlpha_p, maxAlpha_d) / gamma_a;

   // No primal constraints were blocking?
   if (!primalBlocking) {
      alpha_primal = 1.;
   }
   else {
      const double dualValueEstim_p = dualValue_p + maxAlpha_d * dualStep_p;

      if (PIPSisEQ(dualValueEstim_p, 0.)) {
         alpha_primal = 0.; // to be corrected below
      }
      else {
         alpha_primal = (-primalValue_p + mufull / (dualValueEstim_p)) / primalStep_p;
      }
   }

   // No dual constraints were blocking?
   if (!dualBlocking) {
      alpha_dual = 1.;
   }
   else {
      const double primValueEstim_d = primalValue_d + maxAlpha_p * primalStep_d;

      if (PIPSisEQ(primValueEstim_d, 0.)) {
         alpha_dual = 0.; // to be corrected below
      }
      else {
         alpha_dual = (-dualValue_d + mufull / (primValueEstim_d)) / dualStep_d;
      }
   }

   assert(alpha_primal <= 1.);
   assert(alpha_dual <= 1.);

   // safeguard against numerical troubles in the above computations
   alpha_primal = std::min(alpha_primal, maxAlpha_p);
   alpha_dual = std::min(alpha_dual, maxAlpha_d);

   // make it at least gamma_f * maxAlpha and no bigger than 1
   if (alpha_primal < gamma_f * maxAlpha_p)
      alpha_primal = gamma_f * maxAlpha_p;
   if (alpha_dual < gamma_f * maxAlpha_d)
      alpha_dual = gamma_f * maxAlpha_d;

   alpha_primal *= steplength_factor;
   alpha_dual *= steplength_factor;

   assert(alpha_primal < 1. && alpha_dual < 1.);
   assert(alpha_primal >= 0 && alpha_dual >= 0);
}

void MehrotraHeuristic::set_problem_norm(const Problem& problem) {
   dnorm = problem.datanorm();

   if (scaler)
      dnorm_orig = scaler->getDnormOrig();
   else
      dnorm_orig = dnorm;
}

std::pair<double, double> MehrotraHeuristic::compute_unscaled_gap_and_residual_norm(const Residuals& residuals) {
   if (!scaler)
      return std::make_pair(std::fabs(residuals.duality_gap()), residuals.residualNorm());
   else {
      if (!residuals_unscaled)
         residuals_unscaled.reset(scaler->get_unscaled_residuals(residuals));
      else {
         residuals_unscaled->copy(residuals);
         scaler->unscaleResiduals(*residuals_unscaled);
      }

      return std::make_pair(std::fabs(residuals_unscaled->duality_gap()), residuals_unscaled->residualNorm());
   }
}

void MehrotraHeuristic::default_monitor(const Problem* problem /* problem */, const Variables* iterate /* iterate */, const Residuals* residuals,
      double alpha, double sigma, int i, double mu, int status_code, int level) const {
   switch (level) {
      case 0 :
      case 1: {

         const Residuals* residuals_unscaled = residuals;
         if (scaler)
            residuals_unscaled = scaler->get_unscaled_residuals(*residuals);

         const double gap = residuals_unscaled->duality_gap();
         const double rnorm = residuals_unscaled->residualNorm();

         if (scaler)
            delete residuals_unscaled;

         std::cout << std::endl << "Duality Gap: " << gap << std::endl;

         if (i > 1) {
            std::cout << " Number of Corrections = " << number_gondzio_corrections << " alpha = " << alpha << std::endl;
         }
         std::cout << " *** Iteration " << i << " *** " << std::endl;
         std::cout << " mu = " << mu << " relative residual norm = " << rnorm / dnorm_orig << std::endl;

         if (level == 1) {
            // Termination has been detected by the status check; print
            // appropriate message
            if (status_code == SUCCESSFUL_TERMINATION)
               std::cout << std::endl << " *** SUCCESSFUL TERMINATION ***" << std::endl;
            else if (status_code == MAX_ITS_EXCEEDED)
               std::cout << std::endl << " *** MAXIMUM ITERATIONS REACHED *** " << std::endl;
            else if (status_code == INFEASIBLE)
               std::cout << std::endl << " *** TERMINATION: PROBABLY INFEASIBLE *** " << std::endl;
            else if (status_code == UNKNOWN)
               std::cout << std::endl << " *** TERMINATION: STATUS UNKNOWN *** " << std::endl;
         }
      }
         break;
      case 2:
         std::cout << " *** sigma = " << sigma << std::endl;
         break;
   }
}

void MehrotraHeuristic::registerBiCGStabOvserver(AbstractLinearSystem* linear_system) {
   /* every linsys handed to the GondzioStoch should be observable */
   assert(dynamic_cast<Subject*>(linear_system));
   setSubject(dynamic_cast<Subject*>(linear_system));
}

MehrotraHeuristic::~MehrotraHeuristic() {
   delete[] mu_history;
   delete[] residual_norm_history;
   delete[] phi_history;
   delete[] phi_min_history;

   Monitor* m = itsMonitors;
   while (m) {
      Monitor* n = m->nextMonitor;
      delete m;
      m = n;
   }
   delete temp_step;
}