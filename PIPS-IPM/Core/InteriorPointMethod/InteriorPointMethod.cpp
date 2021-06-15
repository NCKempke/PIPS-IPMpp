#include <AbstractOptions.h>
#include "InteriorPointMethod.hpp"
#include "Variables.h"
#include "Residuals.h"
#include "AbstractLinearSystem.h"
#include "Problem.hpp"
#include "DistributedFactory.hpp"

int gLackOfAccuracy = 0;
extern int print_level;

InteriorPointMethod::InteriorPointMethod(DistributedFactory& factory, Problem& problem, MehrotraStrategyType mehrotra_strategy_type, const Scaler* scaler)
: Solver(factory, problem), scaler(scaler), max_iterations(300), dnorm(problem.datanorm()), dnorm_orig(scaler ? scaler->getDnormOrig() : dnorm),
mehrotra_strategy(MehrotraFactory::create(factory, problem, dnorm, mehrotra_strategy_type, scaler)),
      // allocate space to track the sequence problem_formulation complementarity gaps, residual norms, and merit functions.
      mu_history(max_iterations), residual_norm_history(max_iterations), phi_history(max_iterations), phi_min_history(max_iterations)
{
   if (abstract_options::get_bool_parameter("IP_PRINT_TIMESTAMP")) {
      print_timestamp = true;
      start_time = MPI_Wtime();
   }

   if (abstract_options::get_bool_parameter("IP_ACCURACY_REDUCED")) {
      artol = 1.e-3;
      mutol = 1.e-5;
   }
}

TerminationStatus InteriorPointMethod::solve(Problem& problem, Variables& iterate, Residuals& residuals) {
   // initialization of (x,y,z) and factorization routine.
   std::unique_ptr<AbstractLinearSystem> linear_system = factory.make_linear_system(problem);

   // solve the augmented linear system
   this->factory.iterate_started();
   this->solve_linear_system(iterate, problem, residuals, *step, *linear_system);
   this->factory.iterate_ended();
   // register the linear system to the step computation strategy
   this->mehrotra_strategy->register_observer(linear_system.get());

   TerminationStatus status;
   bool termination = false;
   int iteration = 0;
   while (!termination) {
      if (iteration >= max_iterations - 1) {
         status = TerminationStatus::MAX_ITS_EXCEEDED;
         termination = true;
      }
      else {
         residuals.evaluate(problem, iterate);
         double mu = iterate.mu();
         assert(!PIPSisZero(mu));
         if (PIPS_MPIgetRank() == 0 && this->verbose) {
            std::cout << "Norm of x: " << iterate.primals->inf_norm() << "\n";
            std::cout << "Norm of s: " << iterate.slacks->inf_norm() << "\n";
            std::cout << "Norm of eq. duals: " << iterate.equality_duals->inf_norm() << "\n";
            std::cout << "Norm of ineq. duals: " << iterate.inequality_duals->inf_norm() << "\n";
         }

         // termination test
         const auto [duality_gap, residual_norm] = this->compute_unscaled_gap_and_residual_norm(residuals);
         this->update_history(duality_gap, residual_norm, iteration, mu);
         status = this->compute_status(duality_gap, residual_norm, iteration, mu);

         if (status == TerminationStatus::NOT_FINISHED) {
            // run Gondzio's multiple corrector scheme
            this->factory.iterate_started();
            this->mehrotra_strategy->corrector_predictor(problem, iterate, residuals, *step, *linear_system, iteration);
            this->factory.iterate_ended();
            iteration++;
         }
         else {
            termination = true;
            residuals.evaluate(problem, iterate);
            if (print_level >= 10) {
               double mu = iterate.mu();
               this->mehrotra_strategy->print_statistics(&problem, &iterate, &residuals, dnorm, this->mehrotra_strategy->sigma, iteration, mu, status, 1);
            }
         }
      }
   }
   return status;
}

double InteriorPointMethod::barrier_directional_derivative(Problem& problem, Variables& iterate, Variables& direction) {
   double mu = iterate.mu();
   double result = problem.g->dotProductWith(*direction.primals);
   if (0 < problem.number_primal_lower_bounds) {
      result -= direction.primals->special_operation(*iterate.primals, *problem.primal_lower_bounds, *problem.primal_lower_bound_indicators, mu);
   }
   if (0 < problem.number_primal_upper_bounds) {
      result -= direction.primals->special_operation(*iterate.primals, *problem.primal_upper_bounds, *problem.primal_upper_bound_indicators, mu);
   }
   return result;
}

double InteriorPointMethod::predicted_reduction(Problem& problem, Variables& iterate, Variables& direction, double step_length) {
   // scale it with the step length and return the (positive) predicted reduction
   return -step_length * InteriorPointMethod::barrier_directional_derivative(problem, iterate, direction);
}

std::pair<double, double> InteriorPointMethod::compute_unscaled_gap_and_residual_norm(const Residuals& residuals) {
   if (!scaler)
      return std::make_pair(std::fabs(residuals.get_duality_gap()), residuals.get_residual_norm());
   else {
      if (!residuals_unscaled)
         residuals_unscaled.reset(scaler->get_unscaled_residuals(residuals));
      else {
         residuals_unscaled->copy(residuals);
         scaler->unscale_residuals(*residuals_unscaled);
      }

      return std::make_pair(std::fabs(residuals_unscaled->get_duality_gap()), residuals_unscaled->get_residual_norm());
   }
}

void InteriorPointMethod::update_history(double duality_gap, double residual_norm, int iteration, double mu) {
   // store the historical record
   this->mu_history[iteration] = mu;
   this->residual_norm_history[iteration] = residual_norm;
   double phi = (residual_norm + duality_gap) / dnorm_orig;
   this->phi_history[iteration] = phi;

   if (iteration > 0) {
      this->phi_min_history[iteration] = std::min(phi, this->phi_min_history[iteration - 1]);
   }
   else {
      this->phi_min_history[iteration] = phi;
   }
   return;
}

TerminationStatus InteriorPointMethod::compute_status(double duality_gap, double residual_norm, int iteration, double mu) {
   const int myrank = PIPS_MPIgetRank();
   TerminationStatus status = TerminationStatus::NOT_FINISHED;
   double phi = (residual_norm + duality_gap) / dnorm_orig;

   if (mu <= mutol && residual_norm <= artol * dnorm_orig)
      status = TerminationStatus::SUCCESSFUL_TERMINATION;

   if (myrank == 0) {
      std::cout << "mu/mutol: " << mu << "  " << mutol << "  ....   rnorm/limit: " << residual_norm << " " << artol * dnorm_orig << std::endl;

      if (print_timestamp) {
         const double timestamp = MPI_Wtime() - start_time;
         std::cout << "time stamp: " << timestamp << std::endl;
      }
   }

   if (status != TerminationStatus::NOT_FINISHED)
      return status;

   // check infeasibility condition
   if (iteration >= 10 && phi >= 1.e-8 && phi >= 1.e4 * phi_min_history[iteration]) {
#ifdef TIMING
      if( myrank == 0 )
         std::cout << "possible INFEASIBLITY detected, phi: " << phi << std::endl;
#endif
      status = TerminationStatus::INFEASIBLE;
   }

   if (status != TerminationStatus::NOT_FINISHED)
      return status;

   // check for unknown status: slow convergence first
   if (iteration >= 350 && phi_min_history[iteration] >= 0.5 * phi_min_history[iteration - 30]) {
      status = TerminationStatus::UNKNOWN;
      printf("dnorm=%g rnorm=%g artol=%g\n", residual_norm, dnorm_orig, artol);
   }

   if (iteration >= 350 && residual_norm > artol * dnorm_orig &&
       residual_norm_history[iteration] * mu_history[0] >= 1.e8 * mu_history[iteration] * residual_norm_history[0]) {
      status = TerminationStatus::UNKNOWN;
      printf("dnorm=%g rnorm=%g artol=%g\n", residual_norm, dnorm_orig, artol);
   }

   if (mu * dnorm_orig < 1e5 * residual_norm) {
      gLackOfAccuracy = 1;
   }
   else {
      gLackOfAccuracy = 1;
   }
   return status;
}