//
// Created by charlie on 26.04.21.
//

#include <cassert>
#include "InteriorPointMethod.hpp"
#include "PIPSIPMppOptions.h"
#include "Problem.hpp"
#include "Variables.h"
#include "Residuals.h"
#include "DistributedFactory.hpp"
#include "DistributedRootLinearSystem.h"

extern int print_level;
extern double g_iterNumber;

const unsigned int max_linesearch_points = 50;

InteriorPointMethod::InteriorPointMethod(DistributedFactory& factory, Problem& problem, double dnorm, const Scaler* scaler) :
      dnorm(dnorm),
      corrector_step(factory.make_variables(problem)), corrector_residuals(factory.make_residuals(problem)),
      n_linesearch_points(pipsipmpp_options::get_int_parameter("GONDZIO_STOCH_N_LINESEARCH")), temp_step(factory.make_variables(problem)),
      statistics(factory, scaler), bicgstab_skipped(false), bicgstab_converged(true), bigcstab_norm_res_rel(0.),
      bicg_iterations(0), dynamic_corrector_schedule(pipsipmpp_options::get_bool_parameter("GONDZIO_STOCH_USE_DYNAMIC_CORRECTOR_SCHEDULE")),
      additional_correctors_small_comp_pairs(pipsipmpp_options::get_bool_parameter("GONDZIO_STOCH_ADDITIONAL_CORRECTORS_SMALL_VARS")),
      max_additional_correctors(pipsipmpp_options::get_int_parameter("GONDZIO_STOCH_ADDITIONAL_CORRECTORS_MAX")),
      first_iter_small_correctors(pipsipmpp_options::get_int_parameter("GONDZIO_STOCH_FIRST_ITER_SMALL_CORRECTORS")),
      max_alpha_small_correctors(pipsipmpp_options::get_double_parameter("GONDZIO_STOCH_MAX_ALPHA_SMALL_CORRECTORS")), number_small_correctors(0),
      maximum_correctors(options::getIntParameter("GONDZIO_MAX_CORRECTORS")), number_gondzio_corrections(0), step_factor0(0.3), step_factor1(1.5),
      acceptance_tolerance(0.01), beta_min(0.1), beta_max(10), dynamic_bicg_tol(pipsipmpp_options::get_bool_parameter("OUTER_BICG_DYNAMIC_TOL")),
      tsig(3.), pure_centering_step(false), numerical_troubles(false), precond_decreased(true) {
   assert(max_additional_correctors > 0);
   assert(first_iter_small_correctors >= 0);
   assert(0 < max_alpha_small_correctors && max_alpha_small_correctors < 1);
   assert(n_linesearch_points > 0);

   if (abstract_options::get_bool_parameter("IP_ACCURACY_REDUCED")) {
      mutol = 1e-5;
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

   g_iterNumber = 0.;
   this->sigma = 1.;
}

PrimalInteriorPointMethod::PrimalInteriorPointMethod(DistributedFactory& factory, Problem& problem, double dnorm, const Scaler* scaler) :
   InteriorPointMethod(factory, problem, dnorm, scaler), primal_step_length(1.) {
   this->set_BiCGStab_tolerance(-1);
}

PrimalDualInteriorPointMethod::PrimalDualInteriorPointMethod(DistributedFactory& factory, Problem& problem, double dnorm, const Scaler* scaler) :
   InteriorPointMethod(factory, problem, dnorm, scaler), primal_step_length(1.), dual_step_length(1.) {
   this->set_BiCGStab_tolerance(-1);
}

bool InteriorPointMethod::compute_predictor_step(Problem& problem, Variables& current_iterate, Residuals& residuals, Variables& step,
      AbstractLinearSystem& linear_system, int iteration) {
   set_BiCGStab_tolerance(iteration);

   if (print_level >= 10) {
      double mu = current_iterate.mu();
      this->print_statistics(&problem, &current_iterate, &residuals, dnorm, sigma, iteration, mu, TerminationStatus::NOT_FINISHED, 0);
   }

   bool small_corr = false;

   // compute the affine predictor step
   if (!pure_centering_step) {
      residuals.set_complementarity_residual(current_iterate, 0.);
      linear_system.factorize(current_iterate);

      linear_system.solve(current_iterate, residuals, step);
      step.negate();
      check_numerical_troubles(&residuals, numerical_troubles, small_corr);
   }
   else {
      step.set_to_zero();
   }

   // compute fraction-to-boundary rule
   this->project_to_bounds(current_iterate, step);
   return small_corr;
}

void PrimalInteriorPointMethod::project_to_bounds(const Variables& iterate, const Variables& step) {
   this->primal_step_length = iterate.stepbound(step);
}

void PrimalDualInteriorPointMethod::project_to_bounds(const Variables& iterate, const Variables& step) {
   std::tie(this->primal_step_length, this->dual_step_length) = iterate.stepbound_pd(step);
}

void PrimalInteriorPointMethod::compute_corrector_step(Problem& problem, Variables& current_iterate, Residuals& residuals, Variables& step,
      AbstractLinearSystem& linear_system, int iteration, bool small_corr) {
   // calculate centering parameter
   sigma = this->compute_centering_parameter(current_iterate, step);
   double mu = current_iterate.mu();

   if (print_level >= 10) {
      this->print_statistics(&problem, &current_iterate, &residuals, dnorm, sigma, iteration, mu, TerminationStatus::NOT_FINISHED, 2);
   }
   g_iterNumber += 1.;

   // form right hand side of linear system
   corrector_residuals->clear_linear_residuals();
   corrector_residuals->set_complementarity_residual(step, -sigma * mu);

   // solve the linear system
   linear_system.solve(current_iterate, *corrector_residuals, *corrector_step);
   corrector_step->negate();
   check_numerical_troubles(&residuals, numerical_troubles, small_corr);

   // calculate weighted predictor-corrector step
   std::tie(alpha_candidate, weight_candidate) = calculate_alpha_weight_candidate(current_iterate, step, *corrector_step, this->primal_step_length);
   assert(weight_candidate >= 0. && weight_candidate <= 1.);

   step.saxpy(*corrector_step, weight_candidate);

   // prepare for Gondzio corrector loop: zero out the corrector_residuals structure:
   corrector_residuals->clear_linear_residuals();

   // Gondzio correction loop:
   this->gondzio_correction_loop(problem, current_iterate, residuals, step, linear_system, iteration, sigma, mu, small_corr, numerical_troubles);

   // We've finally decided on a step direction, now calculate the length using Mehrotra's heuristic
   this->mehrotra_step_length(current_iterate, step);

   // if we encountered numerical troubles while computing the step, enter a probing round
   if (numerical_troubles) {
      if (precond_decreased)
         precond_decreased = decrease_preconditioner_impact(&linear_system);
      do_probing(problem, current_iterate, residuals, step);
      if (is_poor_step(pure_centering_step, precond_decreased))
         return;
   }
   pure_centering_step = false;
   numerical_troubles = false;

   const double step_inf_norm = step.inf_norm();
   if (PIPS_MPIgetRank() == 0) {
      auto[primal_step_length, dual_step_length] = this->get_step_lengths();
      std::cout << "Direction has length: " << step_inf_norm << "\n";
      std::cout << "Step length: " << primal_step_length << "\n";
   }
}

double PrimalInteriorPointMethod::compute_centering_parameter(Variables& iterate, const Variables& step) {
   double mu = iterate.mu();
   double mu_affine = iterate.mustep_pd(step, this->primal_step_length, this->primal_step_length);

   assert(!PIPSisZero(mu));
   return pow(mu_affine / mu, tsig);
}

double PrimalDualInteriorPointMethod::compute_centering_parameter(Variables& iterate, const Variables& step) {
   double mu = iterate.mu();
   double mu_affine = iterate.mustep_pd(step, this->primal_step_length, this->dual_step_length);

   assert(!PIPSisZero(mu));
   return pow(mu_affine / mu, tsig);
}

void PrimalInteriorPointMethod::take_step(Variables& iterate, Variables& step, double step_length) {
   iterate.saxpy(step, this->primal_step_length*step_length);
}

void PrimalDualInteriorPointMethod::take_step(Variables& iterate, Variables& step, double step_length) {
   iterate.saxpy_pd(step, this->primal_step_length, this->dual_step_length*step_length);
}

void PrimalDualInteriorPointMethod::compute_corrector_step(Problem& problem, Variables& current_iterate, Residuals& residuals, Variables& step,
      AbstractLinearSystem& linear_system, int iteration, bool small_corr) {
   // calculate centering parameter
   sigma = this->compute_centering_parameter(current_iterate, step);
   double mu = current_iterate.mu();

   if (print_level >= 10) {
      this->print_statistics(&problem, &current_iterate, &residuals, dnorm, sigma, iteration, mu, TerminationStatus::NOT_FINISHED, 2);
   }

   g_iterNumber += 1.;

   // *** Corrector step ***
   corrector_residuals->clear_linear_residuals();
   // form right hand side of linear system:
   corrector_residuals->set_complementarity_residual(step, -sigma * mu);

   linear_system.solve(current_iterate, *corrector_residuals, *corrector_step);
   corrector_step->negate();
   check_numerical_troubles(&residuals, numerical_troubles, small_corr);

   // calculate weighted predictor-corrector step
   std::tie(alpha_primal_candidate, alpha_dual_candidate, weight_primal_candidate, weight_dual_candidate) = calculate_alpha_pd_weight_candidate(
         current_iterate, step, *corrector_step, this->primal_step_length, this->dual_step_length);

   assert(weight_primal_candidate >= 0. && weight_primal_candidate <= 1.);
   assert(weight_dual_candidate >= 0. && weight_dual_candidate <= 1.);

   step.saxpy_pd(*corrector_step, weight_primal_candidate, weight_dual_candidate);

   // prepare for Gondzio corrector loop: zero out the corrector_residuals structure:
   corrector_residuals->clear_linear_residuals();

   this->gondzio_correction_loop(problem, current_iterate, residuals, step, linear_system, iteration, sigma, mu, small_corr,
         numerical_troubles);

   // We've finally decided on a step direction, now calculate the length using Mehrotra's heuristic
   this->mehrotra_step_length(current_iterate, step);

   // if we encountered numerical troubles while computing the step check enter a probing round
   if (numerical_troubles) {
      if (precond_decreased)
         precond_decreased = decrease_preconditioner_impact(&linear_system);

      do_probing(problem, current_iterate, residuals, step);
      if (is_poor_step(pure_centering_step, precond_decreased))
         return;
   }

   const double step_inf_norm = step.inf_norm();
   if (PIPS_MPIgetRank() == 0) {
      auto[initial_primal_step_length, initial_dual_step_length] = this->get_step_lengths();
      std::cout << "Direction has length: " << step_inf_norm << "\n";
      std::cout << "Step lengths: " << initial_primal_step_length << " (primal), " << initial_dual_step_length << " (dual)\n";
   }

   pure_centering_step = false;
   numerical_troubles = false;
}

void PrimalDualInteriorPointMethod::gondzio_correction_loop(Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
      AbstractLinearSystem& linear_system, int iteration, double sigma, double mu, bool& small_corr, bool& numerical_troubles) {
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);
   // calculate the target box:
   const double rmin = sigma * mu * beta_min;
   const double rmax = sigma * mu * beta_max;

   number_gondzio_corrections = 0;
   number_small_correctors = 0;

   // enter the Gondzio correction loop:
   while (number_gondzio_corrections < maximum_correctors && number_small_correctors < max_additional_correctors &&
          (PIPSisLT(this->alpha_primal_candidate, 1.) || PIPSisLT(this->alpha_dual_candidate, 1.))) {
      if (dynamic_corrector_schedule)
         adjust_limit_gondzio_correctors();
      corrector_step->copy(iterate);

      const double alpha_p_target = std::min(1., step_factor1 * this->alpha_primal_candidate + step_factor0);
      const double alpha_dual_target = std::min(1., step_factor1 * this->alpha_dual_candidate + step_factor0);

      PIPSdebugMessage("corrector loop: %d alpha_primal: %f alpha_dual %f \n", number_gondzio_corrections, this->alpha_primal_candidate,
               this->alpha_dual_candidate);

      // add a step of this length to corrector_step
      corrector_step->saxpy_pd(step, alpha_p_target, alpha_dual_target);
      // corrector_step is now x_k + alpha_target * delta_p (a trial point)

      /* compute corrector step */
      compute_gondzio_corrector(iterate, linear_system, rmin, rmax, small_corr);
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
            iterate, step, *corrector_step, alpha_p_target, alpha_dual_target);

      // if the enhanced step length is actually 1, make it official
      // and stop correcting
      if (PIPSisEQ(alpha_primal_enhanced, 1.) && PIPSisEQ(alpha_dual_enhanced, 1.)) {
         PIPSdebugMessage("both 1. \n");

         step.saxpy_pd(*corrector_step, weight_primal_candidate, weight_dual_candidate);

         this->alpha_primal_candidate = alpha_primal_enhanced;
         this->alpha_dual_candidate = alpha_dual_enhanced;

         if (small_corr)
            number_small_correctors++;

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
         step.saxpy_pd(*corrector_step, weight_primal_candidate, weight_dual_candidate);

         this->alpha_primal_candidate = alpha_primal_enhanced;
         this->alpha_dual_candidate = alpha_dual_enhanced;

         if (small_corr)
            number_small_correctors++;

         number_gondzio_corrections++;
      }
      else if (alpha_primal_enhanced >= (1. + acceptance_tolerance) * this->alpha_primal_candidate) {
         PIPSdebugMessage("primal better \n");

         step.saxpy_pd(*corrector_step, weight_primal_candidate, 0.);

         this->alpha_primal_candidate = alpha_primal_enhanced;

         if (small_corr)
            number_small_correctors++;

         number_gondzio_corrections++;
      }
      else if (alpha_dual_enhanced >= (1. + acceptance_tolerance) * this->alpha_dual_candidate) {
         PIPSdebugMessage("dual better \n");

         step.saxpy_pd(*corrector_step, 0., weight_dual_candidate);

         this->alpha_dual_candidate = alpha_dual_enhanced;

         if (small_corr)
            number_small_correctors++;

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

void PrimalInteriorPointMethod::gondzio_correction_loop(Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
      AbstractLinearSystem& linear_system, int iteration, double sigma, double mu, bool& small_corr, bool& numerical_troubles) {
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);
   int number_gondzio_corrections = 0;
   number_small_correctors = 0;
   // calculate the target box:
   const double rmin = sigma * mu * beta_min;
   const double rmax = sigma * mu * beta_max;

   while (number_gondzio_corrections < maximum_correctors && number_small_correctors < max_additional_correctors &&
          PIPSisLT(this->alpha_candidate, 1.)) {
      if (dynamic_corrector_schedule)
         adjust_limit_gondzio_correctors();
      corrector_step->copy(iterate);

      // calculate target step length
      double alpha_target = std::min(step_factor1 * this->alpha_candidate + step_factor0, 1.);

      // add a step of this length to corrector_step
      corrector_step->saxpy(step, alpha_target);
      // corrector_step (a trial point) is now x_k + alpha_target * delta_p

      /* compute corrector step */
      compute_gondzio_corrector(iterate, linear_system, rmin, rmax, small_corr);
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
      auto[alpha_enhanced, weight_candidate] = calculate_alpha_weight_candidate(iterate, step, *corrector_step, alpha_target);

      // if the enhanced step length is actually 1, make it official
      // and stop correcting
      if (PIPSisEQ(alpha_enhanced, 1.)) {
         step.saxpy(*corrector_step, weight_candidate);
         this->alpha_candidate = alpha_enhanced;

         if (small_corr)
            number_small_correctors++;

         number_gondzio_corrections++;

         // exit Gondzio correction loop
         break;
      }
      else if (alpha_enhanced >= (1. + acceptance_tolerance) * this->alpha_candidate) {
         // if enhanced step length is significantly better than the
         // current alpha, make the enhanced step official, but maybe
         // keep correcting
         step.saxpy(*corrector_step, weight_candidate);
         this->alpha_candidate = alpha_enhanced;

         if (small_corr)
            number_small_correctors++;

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

void InteriorPointMethod::compute_gondzio_corrector(Variables& iterate, AbstractLinearSystem& linear_system, double rmin, double rmax, bool small_corr) {
   // place XZ into the r3 component of corrector_residuals
   corrector_residuals->set_complementarity_residual(*corrector_step, 0.);
   if (small_corr)
      assert(additional_correctors_small_comp_pairs);
   // do the projection operation
   corrector_residuals->project_r3(rmin, small_corr ? std::numeric_limits<double>::infinity() : rmax);

   // solve for corrector direction
   linear_system.solve(iterate, *corrector_residuals, *corrector_step); // corrector_step is now delta_m
}

std::pair<double, double>
InteriorPointMethod::calculate_alpha_weight_candidate(Variables& iterate, Variables& predictor_step, Variables& corrector_step, double alpha_predictor) {
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

      const double alpha_curr = iterate.stepbound(*temp_step);
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
InteriorPointMethod::calculate_alpha_pd_weight_candidate(Variables& iterate, Variables& predictor_step, Variables& corrector_step, double alpha_primal,
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

      auto[alpha_primal_curr, alpha_dual_curr] = iterate.stepbound_pd(*temp_step);
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

double InteriorPointMethod::compute_probing_factor(Problem& problem, Variables& iterate, Residuals& residuals, Variables& step) {
   const double resids_norm_last = residuals.get_residual_norm();
   this->compute_probing_step(*temp_step, iterate, step);
   residuals.evaluate(problem, *temp_step, false);
   const double resids_norm_probing = residuals.get_residual_norm();
   const double mu_last = iterate.mu();
   const double mu_probing = temp_step->mu();
   return compute_step_factor_probing(resids_norm_last, resids_norm_probing, mu_last, mu_probing);
}

void PrimalInteriorPointMethod::do_probing(Problem& problem, Variables& iterate, Residuals& residuals, Variables& step) {
   const double factor = InteriorPointMethod::compute_probing_factor(problem, iterate, residuals, step);
   this->primal_step_length = factor * this->primal_step_length;
}

void PrimalDualInteriorPointMethod::do_probing(Problem& problem, Variables& iterate, Residuals& residuals, Variables& step) {
   const double factor = InteriorPointMethod::compute_probing_factor(problem, iterate, residuals, step);
   this->primal_step_length = factor * this->primal_step_length;
   this->dual_step_length = factor * this->dual_step_length;
}

bool InteriorPointMethod::is_poor_step(bool& pure_centering_step, bool precond_decreased, double alpha_max) const {
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

bool PrimalInteriorPointMethod::is_poor_step(bool& pure_centering_step, bool precond_decreased) const {
   return InteriorPointMethod::is_poor_step(pure_centering_step, precond_decreased, this->primal_step_length);
}

bool PrimalDualInteriorPointMethod::is_poor_step(bool& pure_centering_step, bool precond_decreased) const {
   const double alpha_max = std::max(this->primal_step_length, this->dual_step_length);
   return InteriorPointMethod::is_poor_step(pure_centering_step, precond_decreased, alpha_max);
}

void PrimalInteriorPointMethod::compute_probing_step(Variables& probing_step, const Variables& iterate, const Variables& step) const {
   probing_step.copy(iterate);
   probing_step.saxpy(step, this->primal_step_length);
}

void PrimalDualInteriorPointMethod::compute_probing_step(Variables& probing_step, const Variables& iterate, const Variables& step) const {
   probing_step.copy(iterate);
   probing_step.saxpy_pd(step, this->primal_step_length, this->dual_step_length);
}

/* when numerical troubles occurred we only allow controlled steps that worsen the residuals and mu by at most a factor of 10 */
double InteriorPointMethod::compute_step_factor_probing(double resids_norm_last, double resids_norm_probing, double mu_last, double mu_probing) {
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

bool InteriorPointMethod::decrease_preconditioner_impact(AbstractLinearSystem* sys) {
   bool success = false;
   dynamic_cast<DistributedRootLinearSystem*>(sys)->precondSC.decreaseDiagDomBound(success);
   if (!success) {
      if (PIPS_MPIgetRank() == 0)
         std::cout << "Cannot increase precision in preconditioner any more\n";
   }
   return success;
}

void InteriorPointMethod::adjust_limit_gondzio_correctors() {
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

void InteriorPointMethod::set_BiCGStab_tolerance(int iteration) const {
   if (!dynamic_bicg_tol)
      return;

   assert(iteration >= -1);

   if (iteration == -1)
      pipsipmpp_options::set_double_parameter("OUTER_BICG_TOL", 1e-10);
   else if (iteration <= 3)
      pipsipmpp_options::set_double_parameter("OUTER_BICG_TOL", 1e-8);
   else if (iteration <= 7)
      pipsipmpp_options::set_double_parameter("OUTER_BICG_TOL", 1e-9);
   else
      pipsipmpp_options::set_double_parameter("OUTER_BICG_TOL", 1e-10);
}

void InteriorPointMethod::check_numerical_troubles(Residuals* residuals, bool& numerical_troubles, bool& small_corr) const {
   if (!bicgstab_converged && bigcstab_norm_res_rel * 1e2 * dnorm > residuals->get_residual_norm()) {
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
PrimalInteriorPointMethod::print_statistics(const Problem* problem, const Variables* iterate, const Residuals* residuals, double dnorm, double sigma,
      int i, double mu, TerminationStatus stop_code, int level) {
   statistics.print(problem, iterate, residuals, dnorm, this->primal_step_length, sigma, i, mu, stop_code, level);
}

void
PrimalDualInteriorPointMethod::print_statistics(const Problem* problem, const Variables* iterate, const Residuals* residuals, double dnorm, double sigma,
      int i, double mu, TerminationStatus stop_code, int level) {
   statistics.print(problem, iterate, residuals, dnorm, this->primal_step_length, this->dual_step_length, sigma, i, mu, stop_code, level);
}

void PrimalInteriorPointMethod::mehrotra_step_length(Variables& iterate, Variables& step) {
   double primalValue = -std::numeric_limits<double>::max();
   double primalStep = -std::numeric_limits<double>::max();
   double dualValue = -std::numeric_limits<double>::max();
   double dualStep = -std::numeric_limits<double>::max();


#ifdef TIMING
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

   int firstOrSecond = -1;
   const double maximum_step_length = iterate.find_blocking(step, primalValue, primalStep, dualValue, dualStep, firstOrSecond);
   const double mu_full = iterate.mustep_pd(step, maximum_step_length, maximum_step_length) / gamma_a;

   this->primal_step_length = 1.;
   switch (firstOrSecond) {
      case 0:
         this->primal_step_length = 1; // No constraints were blocking
         break;
      case 1:
         this->primal_step_length = (-primalValue + mu_full / (dualValue + maximum_step_length * dualStep)) / primalStep;
#ifdef TIMING
         if( myrank == 0 )
            std::cout << "(primal) original alpha " << alpha << std::endl;
#endif
         break;
      case 2:
         this->primal_step_length = (-dualValue + mu_full / (primalValue + maximum_step_length * primalStep)) / dualStep;
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
   this->primal_step_length = std::min(this->primal_step_length, maximum_step_length);
   this->primal_step_length = std::max(this->primal_step_length, gamma_f * maximum_step_length);

   // back off just a touch (or a bit more)
   this->primal_step_length *= steplength_factor;

   assert(0. < this->primal_step_length && this->primal_step_length < 1.);
}

void PrimalDualInteriorPointMethod::mehrotra_step_length(Variables& iterate, Variables& step) {
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

   iterate.find_blocking(step, primalValue_p, primalStep_p, dualValue_p, dualStep_p, primalValue_d, primalStep_d, dualValue_d, dualStep_d, maxAlpha_p,
         maxAlpha_d, primalBlocking, dualBlocking);

   const double mufull = iterate.mustep_pd(step, maxAlpha_p, maxAlpha_d) / gamma_a;

   this->primal_step_length = 1.;
   this->dual_step_length = 1.;
   // No primal constraints were blocking?
   if (!primalBlocking) {
      this->primal_step_length = 1.;
   }
   else {
      const double dualValueEstim_p = dualValue_p + maxAlpha_d * dualStep_p;

      if (PIPSisEQ(dualValueEstim_p, 0.)) {
         this->primal_step_length = 0.; // to be corrected below
      }
      else {
         this->primal_step_length = (-primalValue_p + mufull / (dualValueEstim_p)) / primalStep_p;
      }
   }

   // No dual constraints were blocking?
   if (!dualBlocking) {
      this->dual_step_length = 1.;
   }
   else {
      const double primValueEstim_d = primalValue_d + maxAlpha_p * primalStep_d;

      if (PIPSisEQ(primValueEstim_d, 0.)) {
         this->dual_step_length = 0.; // to be corrected below
      }
      else {
         this->dual_step_length = (-dualValue_d + mufull / (primValueEstim_d)) / dualStep_d;
      }
   }

   assert(this->primal_step_length <= 1.);
   assert(this->dual_step_length <= 1.);

   // safeguard against numerical troubles in the above computations
   this->primal_step_length = std::min(this->primal_step_length, maxAlpha_p);
   this->dual_step_length = std::min(this->dual_step_length, maxAlpha_d);

   // make it at least gamma_f * maxAlpha and no bigger than 1
   if (this->primal_step_length < gamma_f * maxAlpha_p)
      this->primal_step_length = gamma_f * maxAlpha_p;
   if (this->dual_step_length < gamma_f * maxAlpha_d)
      this->dual_step_length = gamma_f * maxAlpha_d;

   this->primal_step_length *= steplength_factor;
   this->dual_step_length *= steplength_factor;

   assert(this->primal_step_length < 1. && this->dual_step_length < 1.);
   assert(this->primal_step_length >= 0 && this->dual_step_length >= 0);
}


void InteriorPointMethod::notify_from_subject() {
   const Subject& subject = *getSubject();

   bicgstab_skipped = subject.getBoolValue("BICG_SKIPPED");
   if (!bicgstab_skipped)
      bicgstab_converged = subject.getBoolValue("BICG_CONVERGED");
   else
      bicgstab_converged = true;
   bigcstab_norm_res_rel = subject.getDoubleValue("BICG_RELRESNORM");
   bicg_iterations = subject.getIntValue("BICG_NITERATIONS");
   if (!bicgstab_converged)
      PIPSdebugMessage("BiGCStab had troubles converging\n");
}

void InteriorPointMethod::register_observer(AbstractLinearSystem* linear_system) {
   /* every linsys handed to the GondzioStoch should be observable */
   assert(dynamic_cast<Subject*>(linear_system));
   set_subject(dynamic_cast<Subject*>(linear_system));
}

std::pair<double, double> PrimalInteriorPointMethod::get_step_lengths() const {
   return std::make_pair(this->primal_step_length, this->primal_step_length);
}

std::pair<double, double> PrimalDualInteriorPointMethod::get_step_lengths() const {
   return std::make_pair(this->primal_step_length, this->dual_step_length);
}

std::unique_ptr<InteriorPointMethod>
MehrotraFactory::create(DistributedFactory& factory, Problem& problem, double dnorm, InteriorPointMethodType interior_point_method_type, const Scaler* scaler) {
   if (interior_point_method_type == InteriorPointMethodType::PRIMAL) {
      return std::make_unique<PrimalInteriorPointMethod>(factory, problem, dnorm, scaler);
   }
   else {
      return std::make_unique<PrimalDualInteriorPointMethod>(factory, problem, dnorm, scaler);
   }
}
