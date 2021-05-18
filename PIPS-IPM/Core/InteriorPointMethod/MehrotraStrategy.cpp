//
// Created by charlie on 26.04.21.
//

#include <cassert>
#include "MehrotraStrategy.hpp"
#include "PIPSIPMppOptions.h"
#include "Problem.h"
#include "Variables.h"
#include "Residuals.h"
#include "DistributedFactory.h"
#include "DistributedRootLinearSystem.h"

extern int print_level;
extern double g_iterNumber;

int gLackOfAccuracy = 0;
const unsigned int max_linesearch_points = 50;

PrimalMehrotraStrategy::PrimalMehrotraStrategy(DistributedFactory& factory, Problem& problem, const Scaler* scaler) : MehrotraStrategy(factory,
      problem, scaler) {
}

PrimalDualMehrotraStrategy::PrimalDualMehrotraStrategy(DistributedFactory& factory, Problem& problem, const Scaler* scaler) : MehrotraStrategy(
      factory, problem, scaler) {
}

std::unique_ptr<MehrotraStrategy>
MehrotraFactory::create(DistributedFactory& factory, Problem& problem, MehrotraHeuristic mehrotra_heuristic, const Scaler* scaler) {
   if (mehrotra_heuristic == PRIMAL) {
      return std::make_unique<PrimalMehrotraStrategy>(factory, problem, scaler);
   }
   else {
      return std::make_unique<PrimalDualMehrotraStrategy>(factory, problem, scaler);
   }
}


MehrotraStrategy::MehrotraStrategy(DistributedFactory& factory, Problem& problem, const Scaler* scaler) : scaler(scaler),
      corrector_step(factory.make_variables(problem)), corrector_residuals(factory.make_residuals(problem)),
      n_linesearch_points(pipsipmpp_options::get_int_parameter("GONDZIO_STOCH_N_LINESEARCH")), temp_step(factory.make_variables(problem)),
      statistics(factory, scaler), bicgstab_skipped(false), bicgstab_converged(true), bigcstab_norm_res_rel(0.), bicg_iterations(0),
      dynamic_corrector_schedule(pipsipmpp_options::get_bool_parameter("GONDZIO_STOCH_USE_DYNAMIC_CORRECTOR_SCHEDULE")),
      additional_correctors_small_comp_pairs(pipsipmpp_options::get_bool_parameter("GONDZIO_STOCH_ADDITIONAL_CORRECTORS_SMALL_VARS")),
      max_additional_correctors(pipsipmpp_options::get_int_parameter("GONDZIO_STOCH_ADDITIONAL_CORRECTORS_MAX")),
      first_iter_small_correctors(pipsipmpp_options::get_int_parameter("GONDZIO_STOCH_FIRST_ITER_SMALL_CORRECTORS")),
      max_alpha_small_correctors(pipsipmpp_options::get_double_parameter("GONDZIO_STOCH_MAX_ALPHA_SMALL_CORRECTORS")), NumberSmallCorrectors(0),
      maximum_correctors(options::getIntParameter("GONDZIO_MAX_CORRECTORS")), number_gondzio_corrections(0), step_factor0(0.3), step_factor1(1.5),
      acceptance_tolerance(0.01), beta_min(0.1), beta_max(10), dynamic_bicg_tol(pipsipmpp_options::get_bool_parameter("OUTER_BICG_DYNAMIC_TOL")),
      tsig(3.), max_iterations(300) {
   assert(max_additional_correctors > 0);
   assert(first_iter_small_correctors >= 0);
   assert(0 < max_alpha_small_correctors && max_alpha_small_correctors < 1);
   assert(n_linesearch_points > 0);

   if (abstract_options::get_bool_parameter("IP_ACCURACY_REDUCED")) {
      artol = 1.e-3;
      mutol = 1.e-5;
   }

   if (pipsipmpp_options::get_bool_parameter("GONDZIO_STOCH_ADAPTIVE_LINESEARCH")) {
      const int size = PIPS_MPIgetSize();

      if (size > 1)
         this->n_linesearch_points = std::min(unsigned(size) + this->n_linesearch_points, max_linesearch_points);
   }

   if (abstract_options::get_bool_parameter("IP_STEPLENGTH_CONSERVATIVE")) {
      steplength_factor = 0.99;
      gamma_f = 0.95;
   }
   gamma_a = 1. / (1. - gamma_f);

   if (abstract_options::get_bool_parameter("IP_PRINT_TIMESTAMP")) {
      print_timestamp = true;
      start_time = MPI_Wtime();
   }

   // allocate space to track the sequence problem_formulation complementarity gaps, residual norms, and merit functions.
   mu_history = new double[max_iterations];
   residual_norm_history = new double[max_iterations];
   phi_history = new double[max_iterations];
   phi_min_history = new double[max_iterations];
}

void PrimalMehrotraStrategy::fraction_to_boundary_rule(Variables& iterate, Variables& step) {
   this->primal_step_length = iterate.stepbound(&step);
   return;
}

double PrimalMehrotraStrategy::compute_centering_parameter(Variables& iterate, Variables& step) {
   double mu = iterate.mu();
   double mu_affine = iterate.mustep_pd(&step, this->primal_step_length, this->primal_step_length);

   assert(!PIPSisZero(mu));
   double sigma = pow(mu_affine / mu, tsig);
   return sigma;
}

void PrimalMehrotraStrategy::take_step(Variables& iterate, Variables& step) {
   iterate.saxpy(&step, this->primal_step_length);
   return;
}

TerminationStatus
PrimalMehrotraStrategy::corrector_predictor(DistributedFactory& factory, Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
      AbstractLinearSystem& linear_system) {
   this->set_problem_norm(problem);
   this->set_BiCGStab_tolerance(-1);

   int iteration = 0;
   double sigma = 1.;
   this->primal_step_length = 1.;

   TerminationStatus status_code;

   g_iterNumber = 0.;

   bool pure_centering_step = false;
   bool numerical_troubles = false;
   bool precond_decreased = true;
   bool termination = false;
   while (!termination) {
      iteration++;
      set_BiCGStab_tolerance(iteration);
      bool small_corr = false;

      factory.iterate_started();

      // evaluate residuals and update algorithm status:
      residuals.evaluate(problem, iterate);

      double mu = iterate.mu();
      assert(!PIPSisZero(mu));

      //  termination test
      status_code = this->compute_status(&problem, &iterate, &residuals, iteration, mu);
      if (status_code == NOT_FINISHED) {
         if (print_level >= 10) {
            this->print_statistics(&problem, &iterate, &residuals, dnorm, sigma, iteration, mu, status_code, 0);
         }

         // predictor step
         if (!pure_centering_step) {
            compute_predictor_step(problem, iterate, residuals, linear_system, step);
            check_numerical_troubles(&residuals, numerical_troubles, small_corr);
         }
         else {
            step.setToZero();
         }

         this->fraction_to_boundary_rule(iterate, step);

         // calculate centering parameter
         sigma = this->compute_centering_parameter(iterate, step);

         if (print_level >= 10) {
            this->print_statistics(&problem, &iterate, &residuals, dnorm, sigma, iteration, mu, status_code, 2);
         }
         g_iterNumber += 1.;

         compute_corrector_step(problem, iterate, linear_system, step, sigma, mu);
         check_numerical_troubles(&residuals, numerical_troubles, small_corr);

         // calculate weighted predictor-corrector step
         std::tie(alpha_candidate, weight_candidate) = calculate_alpha_weight_candidate(&iterate, &step, corrector_step, this->primal_step_length);
         assert(weight_candidate >= 0. && weight_candidate <= 1.);

         step.saxpy(corrector_step, weight_candidate);

         // prepare for Gondzio corrector loop: zero out the corrector_residuals structure:
         corrector_residuals->clear_linear_residuals();

         // Gondzio correction loop:
         this->gondzio_correction_loop(problem, iterate, residuals, step, linear_system, iteration, sigma, mu, small_corr, numerical_troubles);

         // We've finally decided on a step direction, now calculate the length using Mehrotra's heuristic
         this->primal_step_length = mehrotra_step_length_primal(&iterate, &step);
         assert(this->primal_step_length != 0);

         // if we encountered numerical troubles while computing the step, enter a probing round
         if (numerical_troubles) {
            if (precond_decreased)
               precond_decreased = decrease_preconditioner_impact(&linear_system);
            do_probing(&problem, &iterate, &residuals, &step, this->primal_step_length);
            if (is_poor_step(pure_centering_step, precond_decreased, this->primal_step_length))
               continue;
         }
         pure_centering_step = false;
         numerical_troubles = false;

         // take the step (at last!)
         this->take_step(iterate, step);
         factory.iterate_ended();
      }
      else {
         termination = true;
      }
   }
   residuals.evaluate(problem, iterate);
   double mu = iterate.mu();
   if (print_level >= 10) {
      this->print_statistics(&problem, &iterate, &residuals, dnorm, sigma, iteration, mu, status_code, 1);
   }
   return status_code;
}

void PrimalDualMehrotraStrategy::fraction_to_boundary_rule(Variables& iterate, Variables& step) {
   std::tie(this->primal_step_length, this->dual_step_length) = iterate.stepbound_pd(&step);
   return;
}

double PrimalDualMehrotraStrategy::compute_centering_parameter(Variables& iterate, Variables& step) {
   double mu = iterate.mu();
   double mu_affine = iterate.mustep_pd(&step, this->primal_step_length, this->dual_step_length);

   assert(!PIPSisZero(mu));
   double sigma = pow(mu_affine / mu, tsig);
   return sigma;
}

void PrimalDualMehrotraStrategy::take_step(Variables& iterate, Variables& step) {
   iterate.saxpy_pd(&step, this->primal_step_length, this->dual_step_length);
   return;
}

TerminationStatus
PrimalDualMehrotraStrategy::corrector_predictor(DistributedFactory& factory, Problem& problem, Variables& iterate, Residuals& residuals,
      Variables& step, AbstractLinearSystem& linear_system) {
   this->set_problem_norm(problem);
   this->set_BiCGStab_tolerance(-1);

   int iteration = 0;
   double mu = iterate.mu();
   double sigma = 1.;
   this->primal_step_length = 1.;
   this->dual_step_length = 1.;
   TerminationStatus status_code;

   g_iterNumber = 0.;

   bool pure_centering_step = false;
   bool numerical_troubles = false;
   bool precond_decreased = true;

   bool termination = false;
   while (!termination) {
      iteration++;

      set_BiCGStab_tolerance(iteration);
      bool small_corr = false;

      factory.iterate_started();

      // evaluate residuals and update algorithm status:
      residuals.evaluate(problem, iterate);

      //  termination test:
      status_code = this->compute_status(&problem, &iterate, &residuals, iteration, mu);

      if (status_code == NOT_FINISHED) {
         if (print_level >= 10) {
            this->print_statistics(&problem, &iterate, &residuals, dnorm, sigma, iteration, mu, status_code, 0);
         }

         // *** Predictor step ***
         if (!pure_centering_step) {
            compute_predictor_step(problem, iterate, residuals, linear_system, step);
            check_numerical_troubles(&residuals, numerical_troubles, small_corr);
         }
         else
            step.setToZero();

         this->fraction_to_boundary_rule(iterate, step);

         // calculate centering parameter
         sigma = this->compute_centering_parameter(iterate, step);

         if (print_level >= 10) {
            this->print_statistics(&problem, &iterate, &residuals, dnorm, sigma, iteration, mu, status_code, 2);
         }

         g_iterNumber += 1.;

         // *** Corrector step ***
         compute_corrector_step(problem, iterate, linear_system, step, sigma, mu);
         check_numerical_troubles(&residuals, numerical_troubles, small_corr);

         // calculate weighted predictor-corrector step
         std::tie(alpha_primal_candidate, alpha_dual_candidate, weight_primal_candidate, weight_dual_candidate) =
               calculate_alpha_pd_weight_candidate(&iterate, &step, corrector_step, this->primal_step_length, this->dual_step_length);

         assert(weight_primal_candidate >= 0. && weight_primal_candidate <= 1.);
         assert(weight_dual_candidate >= 0. && weight_dual_candidate <= 1.);

         step.saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

         // prepare for Gondzio corrector loop: zero out the corrector_residuals structure:
         corrector_residuals->clear_linear_residuals();

         this->gondzio_correction_loop(problem, iterate, residuals, step, linear_system, iteration, sigma, mu, small_corr, numerical_troubles);

         // We've finally decided on a step direction, now calculate the length using Mehrotra's heuristic.x
         std::tie(this->primal_step_length, this->dual_step_length) = mehrotra_step_length_primal_dual(&iterate, &step);

         // if we encountered numerical troubles while computing the step check enter a probing round
         if (numerical_troubles) {
            if (precond_decreased)
               precond_decreased = decrease_preconditioner_impact(&linear_system);

            do_probing(&problem, &iterate, &residuals, &step, this->primal_step_length, this->dual_step_length);
            const double alpha_max = std::max(this->primal_step_length, this->dual_step_length);

            if (is_poor_step(pure_centering_step, precond_decreased, alpha_max))
               continue;
         }

         // actually take the step and calculate the new mu
         this->take_step(iterate, step);
         mu = iterate.mu();

         pure_centering_step = false;
         numerical_troubles = false;

         factory.iterate_ended();
      }
      else {
         termination = true;
      }
   }
   residuals.evaluate(problem, iterate);
   if (print_level >= 10) {
      this->print_statistics(&problem, &iterate, &residuals, dnorm, sigma, iteration, mu, status_code, 1);
   }
   return status_code;
}

void PrimalDualMehrotraStrategy::gondzio_correction_loop(Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
      AbstractLinearSystem& linear_system, int iteration, double sigma, double mu, bool& small_corr, bool& numerical_troubles) {
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);
   // calculate the target box:
   const double rmin = sigma * mu * beta_min;
   const double rmax = sigma * mu * beta_max;

   number_gondzio_corrections = 0;
   NumberSmallCorrectors = 0;

   // enter the Gondzio correction loop:
   while (number_gondzio_corrections < maximum_correctors && NumberSmallCorrectors < max_additional_correctors &&
          (PIPSisLT(this->alpha_primal_candidate, 1.) || PIPSisLT(this->alpha_dual_candidate, 1.))) {
      if (dynamic_corrector_schedule)
         adjust_limit_gondzio_correctors();
      corrector_step->copy(&iterate);

      const double alpha_p_target = std::min(1., step_factor1 * this->alpha_primal_candidate + step_factor0);
      const double alpha_dual_target = std::min(1., step_factor1 * this->alpha_dual_candidate + step_factor0);

      PIPSdebugMessage("corrector loop: %d alpha_primal: %f alpha_dual %f \n", number_gondzio_corrections, this->alpha_primal_candidate,
            this->alpha_dual_candidate);

      // add a step of this length to corrector_step
      corrector_step->saxpy_pd(&step, alpha_p_target, alpha_dual_target);
      // corrector_step is now x_k + alpha_target * delta_p (a trial point)

      /* compute corrector step */
      compute_gondzio_corrector(problem, iterate, linear_system, rmin, rmax, small_corr);
      const bool was_small_corr = small_corr;
      check_numerical_troubles(&residuals, numerical_troubles, small_corr);

      if (numerical_troubles) {
         if (!was_small_corr && small_corr)
            continue;
         else
            // exit corrector loop if small correctors have already been tried or are not allowed
            break;
      }

      // calculate weighted predictor-corrector step
      std::tie(alpha_primal_enhanced, alpha_dual_enhanced, weight_primal_candidate, weight_dual_candidate) = calculate_alpha_pd_weight_candidate(
            &iterate, &step, corrector_step, alpha_p_target, alpha_dual_target);

      // if the enhanced step length is actually 1, make it official
      // and stop correcting
      if (PIPSisEQ(alpha_primal_enhanced, 1.) && PIPSisEQ(alpha_dual_enhanced, 1.)) {
         PIPSdebugMessage("both 1. \n");

         step.saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

         this->alpha_primal_candidate = alpha_primal_enhanced;
         this->alpha_dual_candidate = alpha_dual_enhanced;

         if (small_corr)
            NumberSmallCorrectors++;

         number_gondzio_corrections++;

         // exit Gondzio correction loop
         break;
      }
      else if (alpha_primal_enhanced >= (1. + acceptance_tolerance) * this->alpha_primal_candidate &&
               alpha_dual_enhanced >= (1. + acceptance_tolerance) * this->alpha_dual_candidate) {
         PIPSdebugMessage("both better \n");

         // if enhanced step length is significantly better than the
         // current alpha, make the enhanced step official, but maybe
         // keep correcting
         step.saxpy_pd(corrector_step, weight_primal_candidate, weight_dual_candidate);

         this->alpha_primal_candidate = alpha_primal_enhanced;
         this->alpha_dual_candidate = alpha_dual_enhanced;

         if (small_corr)
            NumberSmallCorrectors++;

         number_gondzio_corrections++;
      }
      else if (alpha_primal_enhanced >= (1. + acceptance_tolerance) * this->alpha_primal_candidate) {
         PIPSdebugMessage("primal better \n");

         step.saxpy_pd(corrector_step, weight_primal_candidate, 0.);

         this->alpha_primal_candidate = alpha_primal_enhanced;

         if (small_corr)
            NumberSmallCorrectors++;

         number_gondzio_corrections++;
      }
      else if (alpha_dual_enhanced >= (1. + acceptance_tolerance) * this->alpha_dual_candidate) {
         PIPSdebugMessage("dual better \n");

         step.saxpy_pd(corrector_step, 0., weight_dual_candidate);

         this->alpha_dual_candidate = alpha_dual_enhanced;

         if (small_corr)
            NumberSmallCorrectors++;

         number_gondzio_corrections++;
      }
         /* if not done yet because correctors were not good enough - try a small corrector if enabled */
      else if (additional_correctors_small_comp_pairs && !small_corr && iteration >= first_iter_small_correctors) {
         if (this->alpha_primal_candidate < max_alpha_small_correctors || this->alpha_dual_candidate < max_alpha_small_correctors) {
            // try and center small pairs
            small_corr = true;
            if (my_rank == 0) {
               std::cout << "Switching to small corrector " << std::endl;
               std::cout << "Alpha when switching: " << this->alpha_primal_candidate << " " << this->alpha_dual_candidate << std::endl;
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
}

void PrimalMehrotraStrategy::gondzio_correction_loop(Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
      AbstractLinearSystem& linear_system, int iteration, double sigma, double mu, bool& small_corr, bool& numerical_troubles) {
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);
   int number_gondzio_corrections = 0;
   NumberSmallCorrectors = 0;
   // calculate the target box:
   const double rmin = sigma * mu * beta_min;
   const double rmax = sigma * mu * beta_max;

   while (number_gondzio_corrections < maximum_correctors && NumberSmallCorrectors < max_additional_correctors && PIPSisLT(this->alpha_candidate, 1.)) {
      if (dynamic_corrector_schedule)
         adjust_limit_gondzio_correctors();
      corrector_step->copy(&iterate);

      // calculate target step length
      double alpha_target = std::min(step_factor1 * this->alpha_candidate + step_factor0, 1.);

      // add a step of this length to corrector_step
      corrector_step->saxpy(&step, alpha_target);
      // corrector_step (a trial point) is now x_k + alpha_target * delta_p

      /* compute corrector step */
      compute_gondzio_corrector(problem, iterate, linear_system, rmin, rmax, small_corr);
      const bool was_small_corr = small_corr;
      check_numerical_troubles(&residuals, numerical_troubles, small_corr);

      if (numerical_troubles) {
         if (!was_small_corr && small_corr)
            continue;
         else {
            // exit corrector loop if small correctors have already been tried or are not allowed
            break;
         }
      }

      // calculate weighted predictor-corrector step
      auto[alpha_enhanced, weight_candidate] = calculate_alpha_weight_candidate(&iterate, &step, corrector_step, alpha_target);

      // if the enhanced step length is actually 1, make it official
      // and stop correcting
      if (PIPSisEQ(alpha_enhanced, 1.)) {
         step.saxpy(corrector_step, weight_candidate);
         this->alpha_candidate = alpha_enhanced;

         if (small_corr)
            NumberSmallCorrectors++;

         number_gondzio_corrections++;

         // exit Gondzio correction loop
         break;
      }
      else if (alpha_enhanced >= (1. + acceptance_tolerance) * this->alpha_candidate) {
         // if enhanced step length is significantly better than the
         // current alpha, make the enhanced step official, but maybe
         // keep correcting
         step.saxpy(corrector_step, weight_candidate);
         this->alpha_candidate = alpha_enhanced;

         if (small_corr)
            NumberSmallCorrectors++;

         number_gondzio_corrections++;
      }
         /* if not done yet because correctors were not good enough - try a small corrector if enabled */
      else if (additional_correctors_small_comp_pairs && !small_corr && iteration >= first_iter_small_correctors) {
         if (this->alpha_candidate < max_alpha_small_correctors) {
            small_corr = true;
            if (my_rank == 0) {
               std::cout << "Switching to small corrector " << std::endl;
               std::cout << "Alpha when switching: " << this->alpha_candidate << std::endl;
            }
         }
         else {
            // exit Gondzio correction loop
            break;
         }
      }
      else {
         // exit Gondzio correction loop
         break;
      }
   }
}

void MehrotraStrategy::compute_predictor_step(Problem& problem, Variables& iterate, Residuals& residuals, AbstractLinearSystem& linear_system,
      Variables& step) {
   residuals.set_complementarity_residual(iterate, 0.);
   linear_system.factorize(&problem, &iterate);

   linear_system.solve(&problem, &iterate, &residuals, &step);
   step.negate();
}

void
MehrotraStrategy::compute_corrector_step(Problem& problem, Variables& iterate, AbstractLinearSystem& linear_system, Variables& step, double sigma,
      double mu) {
   corrector_residuals->clear_linear_residuals();
   // form right hand side of linear system:
   corrector_residuals->set_complementarity_residual(step, -sigma * mu);

   linear_system.solve(&problem, &iterate, corrector_residuals, corrector_step);
   corrector_step->negate();
}

void MehrotraStrategy::compute_gondzio_corrector(Problem& problem, Variables& iterate, AbstractLinearSystem& linear_system, double rmin, double rmax,
      bool small_corr) {
   // place XZ into the r3 component of corrector_residuals
   corrector_residuals->set_complementarity_residual(*corrector_step, 0.);
   if (small_corr)
      assert(additional_correctors_small_comp_pairs);
   // do the projection operation
   corrector_residuals->project_r3(rmin, small_corr ? std::numeric_limits<double>::infinity() : rmax);

   // solve for corrector direction
   linear_system.solve(&problem, &iterate, corrector_residuals, corrector_step); // corrector_step is now delta_m
}

std::pair<double, double>
MehrotraStrategy::calculate_alpha_weight_candidate(Variables* iterate, Variables* predictor_step, Variables* corrector_step, double alpha_predictor) {
   assert(corrector_step);
   assert(predictor_step);
   assert(temp_step);

   assert(alpha_predictor > 0. && alpha_predictor <= 1.);

   double alpha_candidate = -1.;
   double weight_candidate = -1.;
   const double weight_min = alpha_predictor * alpha_predictor;
   const double weight_interval_length = 1. - weight_min;

   // main loop
   for (unsigned int n = 0; n <= n_linesearch_points; n++) {
      double weight_curr = std::min(1., weight_min + (weight_interval_length / (n_linesearch_points)) * n);
      assert(weight_curr > 0. && weight_curr <= 1.);

      temp_step->copy(predictor_step);
      temp_step->saxpy(corrector_step, weight_curr);

      const double alpha_curr = iterate->stepbound(temp_step);
      assert(alpha_curr > 0. && alpha_curr <= 1.);

      if (alpha_curr > alpha_candidate) {
         alpha_candidate = alpha_curr;
         weight_candidate = weight_curr;
      }
   }
   assert(alpha_candidate >= 0. && weight_candidate >= 0.);
   return std::make_pair(alpha_candidate, weight_candidate);
}

std::tuple<double, double, double, double>
MehrotraStrategy::calculate_alpha_pd_weight_candidate(Variables* iterate, Variables* predictor_step, Variables* corrector_step, double alpha_primal,
      double alpha_dual) {
   assert(alpha_primal > 0. && alpha_primal <= 1.);
   assert(alpha_dual > 0. && alpha_dual <= 1.);

   double alpha_primal_best = -1., alpha_dual_best = -1.;
   double weight_primal_best = -1., weight_dual_best = -1.;
   const double weight_min = alpha_primal * alpha_dual;
   const double weight_interval_length = 1. - weight_min;

   // main loop
   for (unsigned int n = 0; n <= n_linesearch_points; n++) {
      double weight_curr = std::min(1., weight_min + (weight_interval_length / (n_linesearch_points)) * n);
      assert(weight_curr > 0. && weight_curr <= 1.);

      temp_step->copy(predictor_step);
      temp_step->saxpy(corrector_step, weight_curr);

      auto[alpha_primal_curr, alpha_dual_curr] = iterate->stepbound_pd(temp_step);
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
   return std::make_tuple(alpha_primal_best, alpha_dual_best, weight_primal_best, weight_dual_best);
}

void MehrotraStrategy::do_probing(Problem* problem, Variables* iterate, Residuals* residuals, Variables* step, double& alpha) {
   const double mu_last = iterate->mu();
   const double resids_norm_last = residuals->residual_norm();

   compute_probing_step(temp_step, iterate, step, alpha);

   residuals->evaluate(*problem, *temp_step, false);
   const double mu_probing = temp_step->mu();
   const double resids_norm_probing = residuals->residual_norm();

   const double factor = compute_step_factor_probing(resids_norm_last, resids_norm_probing, mu_last, mu_probing);

   alpha = factor * alpha;
}

void
MehrotraStrategy::do_probing(Problem* problem, Variables* iterate, Residuals* residuals, Variables* step, double& alpha_primal, double& alpha_dual) {
   const double mu_last = iterate->mu();
   const double resids_norm_last = residuals->residual_norm();

   this->compute_probing_step(temp_step, iterate, step, alpha_primal, alpha_dual);

   residuals->evaluate(*problem, *temp_step, false);
   const double mu_probing = temp_step->mu();
   const double resids_norm_probing = residuals->residual_norm();

   const double factor = this->compute_step_factor_probing(resids_norm_last, resids_norm_probing, mu_last, mu_probing);

   alpha_primal = factor * alpha_primal;
   alpha_dual = factor * alpha_dual;
}

bool MehrotraStrategy::is_poor_step(bool& pure_centering_step, bool precond_decreased, double alpha_max) const {
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

void MehrotraStrategy::compute_probing_step(Variables* probing_step, const Variables* iterate, const Variables* step, double alpha) const {
   probing_step->copy(iterate);
   probing_step->saxpy(step, alpha);
}

void MehrotraStrategy::compute_probing_step(Variables* probing_step, const Variables* iterate, const Variables* step, double alpha_primal,
      double alpha_dual) const {
   probing_step->copy(iterate);
   probing_step->saxpy_pd(step, alpha_primal, alpha_dual);
}

/* when numerical troubles occurred we only allow controlled steps that worsen the residuals and mu by at most a factor of 10 */
double MehrotraStrategy::compute_step_factor_probing(double resids_norm_last, double resids_norm_probing, double mu_last, double mu_probing) const {
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

bool MehrotraStrategy::decrease_preconditioner_impact(AbstractLinearSystem* sys) const {
   bool success = false;
   dynamic_cast<DistributedRootLinearSystem*>(sys)->precondSC.decreaseDiagDomBound(success);
   if (!success) {
      if (PIPS_MPIgetRank() == 0)
         std::cout << "Cannot increase precision in preconditioner any more\n";
   }
   return success;
}

void MehrotraStrategy::adjust_limit_gondzio_correctors() {
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

void MehrotraStrategy::set_BiCGStab_tolerance(int iteration) const {
   if (!dynamic_bicg_tol)
      return;

   assert(iteration >= -1);

   if (iteration == -1)
      pipsipmpp_options::set_double_parameter("OUTER_BICG_TOL", 1e-10);
   else if (iteration <= 4)
      pipsipmpp_options::set_double_parameter("OUTER_BICG_TOL", 1e-8);
   else if (iteration <= 8)
      pipsipmpp_options::set_double_parameter("OUTER_BICG_TOL", 1e-9);
   else
      pipsipmpp_options::set_double_parameter("OUTER_BICG_TOL", 1e-10);
}

void MehrotraStrategy::check_numerical_troubles(Residuals* residuals, bool& numerical_troubles, bool& small_corr) const {
   if (!bicgstab_converged && bigcstab_norm_res_rel * 1e2 * dnorm > residuals->residual_norm()) {
      PIPSdebugMessage("Step computation in BiCGStab failed");
      numerical_troubles = true;
      if (additional_correctors_small_comp_pairs && !small_corr) {
         if (PIPS_MPIgetRank() == 0) {
            std::cout << "switching to small correctors\n";
         }
         small_corr = true;
      }
   }
}

void
PrimalMehrotraStrategy::print_statistics(const Problem* problem, const Variables* iterate, const Residuals* residuals, double dnorm, double sigma,
      int i, double mu, int stop_code, int level) {
   statistics.print(problem, iterate, residuals, dnorm, this->primal_step_length, sigma, i, mu, stop_code, level);
}

void
PrimalDualMehrotraStrategy::print_statistics(const Problem* problem, const Variables* iterate, const Residuals* residuals, double dnorm, double sigma,
      int i, double mu, int stop_code, int level) {
   statistics.print(problem, iterate, residuals, dnorm, this->primal_step_length, this->dual_step_length, sigma, i, mu, stop_code, level);
}

void
MehrotraStrategy::print_statistics(const Problem* problem, const Variables* iterate, const Residuals* residuals, double dnorm, double alpha_primal,
      double alpha_dual, double sigma, int i, double mu, int stop_code, int level) {
   statistics.print(problem, iterate, residuals, dnorm, alpha_primal, alpha_dual, sigma, i, mu, stop_code, level);
}

TerminationStatus
MehrotraStrategy::compute_status(const Problem* data, const Variables* iterate /* iterate */, const Residuals* residuals, int iteration, double mu) {
   const int myrank = PIPS_MPIgetRank();
   TerminationStatus status = NOT_FINISHED;

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
      status = MAX_ITS_EXCEEDED;
   else if (mu <= mutol && rnorm <= artol * dnorm_orig)
      status = SUCCESSFUL_TERMINATION;

   if (myrank == 0) {
      std::cout << "mu/mutol: " << mu << "  " << mutol << "  ....   rnorm/limit: " << rnorm << " " << artol * dnorm_orig << std::endl;

      if (print_timestamp) {
         const double timestamp = MPI_Wtime() - start_time;
         std::cout << "time stamp: " << timestamp << std::endl;
      }
   }

   if (status != NOT_FINISHED)
      return status;

   // check infeasibility condition
   if (index >= 10 && phi >= 1.e-8 && phi >= 1.e4 * phi_min_history[index]) {
#ifdef TIMING
      if( myrank == 0 )
         std::cout << "possible INFEASIBLITY detected, phi: " << phi << std::endl;
#endif
      status = INFEASIBLE;
   }

   if (status != NOT_FINISHED)
      return status;

   // check for unknown status: slow convergence first
   if (index >= 350 && phi_min_history[index] >= 0.5 * phi_min_history[index - 30]) {
      status = UNKNOWN;
      printf("dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm_orig, artol);
   }

   if (index >= 350 && rnorm > artol * dnorm_orig &&
       residual_norm_history[index] * mu_history[0] >= 1.e8 * mu_history[index] * residual_norm_history[0]) {
      status = UNKNOWN;
      printf("dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm_orig, artol);
   }

   if (mu * dnorm_orig < 1e5 * rnorm) {
      gLackOfAccuracy = 1;
   }
   else {
      gLackOfAccuracy = 1;
   }
   return status;
}

double MehrotraStrategy::mehrotra_step_length_primal(Variables* iterate, Variables* step) {
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

std::pair<double, double> MehrotraStrategy::mehrotra_step_length_primal_dual(Variables* iterate, Variables* step) {
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

   double alpha_primal = 1., alpha_dual = 1.;
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
   return std::make_pair(alpha_primal, alpha_dual);
}

void MehrotraStrategy::set_problem_norm(const Problem& problem) {
   dnorm = problem.datanorm();

   if (scaler)
      dnorm_orig = scaler->getDnormOrig();
   else
      dnorm_orig = dnorm;
}

std::pair<double, double> MehrotraStrategy::compute_unscaled_gap_and_residual_norm(const Residuals& residuals) {
   if (!scaler)
      return std::make_pair(std::fabs(residuals.duality_gap()), residuals.residual_norm());
   else {
      if (!residuals_unscaled)
         residuals_unscaled.reset(scaler->get_unscaled_residuals(residuals));
      else {
         residuals_unscaled->copy(residuals);
         scaler->unscaleResiduals(*residuals_unscaled);
      }

      return std::make_pair(std::fabs(residuals_unscaled->duality_gap()), residuals_unscaled->residual_norm());
   }
}

void MehrotraStrategy::default_monitor(const Problem* problem /* problem */, const Variables* iterate /* iterate */, const Residuals* residuals,
      double alpha, double sigma, int i, double mu, int status_code, int level) const {
   switch (level) {
      case 0 :
      case 1: {

         const Residuals* residuals_unscaled = residuals;
         if (scaler)
            residuals_unscaled = scaler->get_unscaled_residuals(*residuals);

         const double gap = residuals_unscaled->duality_gap();
         const double rnorm = residuals_unscaled->residual_norm();

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

MehrotraStrategy::~MehrotraStrategy() {
   delete[] mu_history;
   delete[] residual_norm_history;
   delete[] phi_history;
   delete[] phi_min_history;
   delete temp_step;
}
