#include "Solver.h"
#include "OoqpMonitor.h"
#include "Status.h"
#include "Problem.h"
#include "Variables.h"
#include "Residuals.h"
#include "AbstractLinearSystem.h"
#include "Options.h"
#include "ProblemFactory.h"
#include "QpGenOptions.h"
#include <iostream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <limits>
#include "mpi.h"

double g_iterNumber = 0.0;
int gOoqpPrintLevel = 1000;
int gLackOfAccuracy = 0;
int onSafeSolver = 0;

int gOuterBiCGFails = 0;
int gOuterBiCGIter = 0;
double gOuterBiCGIterAvg = 0.0;
int gInnerBiCGIter = 0;
int gInnerBiCGFails = 0;

// gmu is needed by MA57!
double gmu;
// double grnorm;
extern int gOoqpPrintLevel;

Solver::Solver(ProblemFactory& problem_formulation, Problem& problem, const Scaler* scaler) : scaler(scaler), factory(problem_formulation),
      step(factory.make_variables(problem)), corrector_step(factory.make_variables(problem)), corrector_residuals(factory.make_residuals(problem)) {
   if (base_options::getBoolParameter("IP_STEPLENGTH_CONSERVATIVE")) {
      steplength_factor = 0.99;
      gamma_f = 0.95;
   }

   gamma_a = 1.0 / (1.0 - gamma_f);

   if (base_options::getBoolParameter("IP_ACCURACY_REDUCED")) {
      artol = 1.e-3;
      mutol = 1.e-5;
   }

   if (base_options::getBoolParameter("IP_PRINT_TIMESTAMP")) {
      print_timestamp = true;
      start_time = MPI_Wtime();
   }
   max_iterations = 300;
   print_level = 0; // has no meaning right now
   tsig = 3.0;     // the usual value for the centering exponent (tau)

   number_gondzio_corrections = 0;

   maximum_correctors = qpgen_options::getIntParameter("GONDZIO_MAX_CORRECTORS");

   // the two StepFactor constants set targets for increase in step
   // length for each corrector
   step_factor0 = 0.08;
   step_factor1 = 1.08;

   // accept the enhanced step if it produces a small improvement in the step length
   acceptance_tolerance = 0.01;

   //define the Gondzio correction box
   beta_min = 0.1;
   beta_max = 10.0;

   // allocate space to track the sequence problem_formulation complementarity gaps,
   // residual norms, and merit functions.
   mu_history = new double[max_iterations];
   residual_norm_history = new double[max_iterations];
   phi_history = new double[max_iterations];
   phi_min_history = new double[max_iterations];

   // Use the defaultStatus method
   status = 0;
}

void Solver::solve_linear_system(Variables& iterate, Problem& problem, Residuals& residuals, Variables& step) {
   double problem_norm = std::sqrt(dnorm);
   iterate.push_to_interior(problem_norm, problem_norm);

   residuals.evaluate(problem, iterate);
   residuals.set_complementarity_residual(iterate, 0.0);

   linear_system->factorize(&problem, &iterate);
   linear_system->solve(&problem, &iterate, &residuals, &step);

   step.negate();

   // Take the full affine scaling step
   iterate.saxpy(&step, 1.0);
   double shift = 1.e3 + 2 * iterate.violation();
   iterate.shiftBoundVariables(shift, shift);
   return;
}

double Solver::mehrotra_step_length(Variables* iterate, Variables* step) {
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

   double step_length = 1.0;
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

   assert(step_length < 1.0);

   return step_length;
}

void Solver::mehrotra_step_length_PD(Variables* iterate, Variables* step, double& alpha_primal, double& alpha_dual) {
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
      alpha_primal = 1.0;
   }
   else {
      const double dualValueEstim_p = dualValue_p + maxAlpha_d * dualStep_p;

      if (PIPSisEQ(dualValueEstim_p, 0.0)) {
         alpha_primal = 0.0; // to be corrected below
      }
      else {
         alpha_primal = (-primalValue_p + mufull / (dualValueEstim_p)) / primalStep_p;
      }
   }

   // No dual constraints were blocking?
   if (!dualBlocking) {
      alpha_dual = 1.0;
   }
   else {
      const double primValueEstim_d = primalValue_d + maxAlpha_p * primalStep_d;

      if (PIPSisEQ(primValueEstim_d, 0.0)) {
         alpha_dual = 0.0; // to be corrected below
      }
      else {
         alpha_dual = (-dualValue_d + mufull / (primValueEstim_d)) / dualStep_d;
      }
   }

   assert(alpha_primal <= 1.0);
   assert(alpha_dual <= 1.0);

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

   assert(alpha_primal < 1.0 && alpha_dual < 1.0);
   assert(alpha_primal >= 0 && alpha_dual >= 0);
}


void Solver::do_monitor(const Problem* data, const Variables* iterate, const Residuals* residuals, double alpha, double sigma, int i, double mu,
      int stop_code, int level) {
   OoqpMonitor* m = itsMonitors;

   while (m) {
      m->doIt(this, data, iterate, residuals, alpha, sigma, i, mu, stop_code, level);
      m = m->nextMonitor;
   }
}

void
Solver::do_monitor_Pd(const Problem* data, const Variables* iterate, const Residuals* residuals, double alpha_primal, double alpha_dual, double sigma,
      int i, double mu, int stop_code, int level) {
   OoqpMonitor* m = itsMonitors;

   while (m) {
      m->doItPd(this, data, iterate, residuals, alpha_primal, alpha_dual, sigma, i, mu, stop_code, level);
      m = m->nextMonitor;
   }
}


TerminationCode
Solver::do_status(const Problem* problem, const Variables* iterate, const Residuals* residuals, int i, double mu, TerminationCode level) {
   if (status) {
      return status->doIt(this, problem, iterate, residuals, i, mu, level);
   }
   else {
      return this->default_status(problem, iterate, residuals, i, mu, level);
   }
}


void Solver::monitorSelf() {
   this->add_monitor(new OoqpSelfMonitor);
}

void Solver::add_monitor(OoqpMonitor* m) {
   // Push the monitor onto the list
   m->nextMonitor = itsMonitors;
   itsMonitors = m;
}

std::pair<double, double> Solver::compute_unscaled_gap_and_residual_norm(const Residuals& residuals) {
   if (!scaler)
      return std::make_pair(std::fabs(residuals.duality_gap()), residuals.residualNorm());
   else {
      if (!residuals_unscaled)
         residuals_unscaled.reset(scaler->getResidualsUnscaled(residuals));
      else {
         residuals_unscaled->copy(residuals);
         scaler->unscaleResiduals(*residuals_unscaled);
      }

      return std::make_pair(std::fabs(residuals_unscaled->duality_gap()), residuals_unscaled->residualNorm());
   }
}


TerminationCode
Solver::default_status(const Problem* data, const Variables* iterate /* iterate */, const Residuals* residuals, int iteration, double mu,
      TerminationCode /* level */level) {
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
      printf("hehe dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm_orig, artol);
   }

   if (index >= 350 && rnorm > artol * dnorm_orig &&
       residual_norm_history[index] * mu_history[0] >= 1.e8 * mu_history[index] * residual_norm_history[0]) {
      stop_code = UNKNOWN;
      printf("dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm_orig, artol);
   }

   if (mu * dnorm_orig < 1.0e5 * rnorm) {
      gLackOfAccuracy = 1;
   }
   else {
      gLackOfAccuracy = 1;
   }
   return stop_code;
}

void Solver::set_problem_norm(const Problem& problem) {
   dnorm = problem.datanorm();

   if (scaler)
      dnorm_orig = scaler->getDnormOrig();
   else
      dnorm_orig = dnorm;
}

void Solver::default_monitor(const Problem* problem /* problem */, const Variables* iterate /* iterate */, const Residuals* residuals, double alpha,
      double sigma, int i, double mu, int status_code, int level) const {
   switch (level) {
      case 0 :
      case 1: {

         const Residuals* residuals_unscaled = residuals;
         if (scaler)
            residuals_unscaled = scaler->getResidualsUnscaled(*residuals);

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

Solver::~Solver() {
   OoqpMonitor* m = itsMonitors;
   while (m) {
      OoqpMonitor* n = m->nextMonitor;
      delete m;
      m = n;
   }

   delete corrector_residuals;
   delete corrector_step;
   delete step;
   delete linear_system;

   delete[] mu_history;
   delete[] residual_norm_history;
   delete[] phi_history;
   delete[] phi_min_history;
}
