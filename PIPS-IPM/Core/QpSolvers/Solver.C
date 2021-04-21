#include "Solver.h"
#include "OoqpMonitor.h"
#include "Status.h"
#include "Problem.h"
#include "Variables.h"
#include "Residuals.h"
#include "LinearSystem.h"
#include "OoqpStartStrategy.h"
#include "Options.h"
#include "ProblemFormulation.h"
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

Solver::Solver(ProblemFormulation& problem_formulation, Problem& problem, const Scaler* scaler) : scaler(scaler), factory(problem_formulation) {
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
      printTimeStamp = true;
      startTime = MPI_Wtime();
   }
   step = factory.makeVariables(&problem);
   corrector_step = factory.makeVariables(&problem);
   corrector_residuals = factory.makeResiduals(&problem);

   maxit = 300;
   printlevel = 0; // has no meaning right now
   tsig = 3.0;     // the usual value for the centering exponent (tau)

   NumberGondzioCorrections = 0;

   maximum_correctors = qpgen_options::getIntParameter("GONDZIO_MAX_CORRECTORS");

   // the two StepFactor constants set targets for increase in step
   // length for each corrector
   StepFactor0 = 0.08;
   StepFactor1 = 1.08;

   // accept the enhanced step if it produces a small improvement in the step length
   AcceptTol = 0.01;

   //define the Gondzio correction box
   beta_min = 0.1;
   beta_max = 10.0;

   // allocate space to track the sequence problem_formulation complementarity gaps,
   // residual norms, and merit functions.
   mu_history = new double[maxit];
   rnorm_history = new double[maxit];
   phi_history = new double[maxit];
   phi_min_history = new double[maxit];

   // Use the defaultStatus method
   status = 0;
}

void Solver::start(ProblemFormulation* formulation, Variables* iterate, Problem* prob, Residuals* resid, Variables* step) {
   if (startStrategy) {
      startStrategy->doIt(this, formulation, iterate, prob, resid, step);
   }
   else {
      this->defaultStart(formulation, iterate, prob, resid, step);
      //this->dumbstart( formulation, iterate, prob, resid, step );
   }
}

void Solver::defaultStart(ProblemFormulation* /* formulation */, Variables* iterate, Problem* prob, Residuals* resid, Variables* step) {
   double sdatanorm = std::sqrt(dnorm);
   double a = sdatanorm;
   double b = sdatanorm;
   iterate->interiorPoint(a, b);

   resid->evaluate(*prob, iterate);
   resid->set_r3_xz_alpha(iterate, 0.0);

   linear_system->factorize(prob, iterate);
   linear_system->solve(prob, iterate, resid, step);

   step->negate();

   // Take the full affine scaling step
   iterate->saxpy(step, 1.0);
   // resid->evaluate(*prob, iterate); // Calc the resids if debugging.
   double shift = 1.e3 + 2 * iterate->violation();
   iterate->shiftBoundVariables(shift, shift);
}

void Solver::dumbstart(ProblemFormulation* /* formulation */, Variables* iterate, Problem* /* prob */, Residuals* /*resid*/, Variables* /*step*/  ) {
   const double a = 1.e3;
   const double b = 1.e5;

   const double bigstart = a * dnorm + b;
   iterate->interiorPoint(bigstart, bigstart);
}

double Solver::finalStepLength(Variables* iterate, Variables* step) {
   double primalValue = -std::numeric_limits<double>::max();
   double primalStep = -std::numeric_limits<double>::max();
   double dualValue = -std::numeric_limits<double>::max();
   double dualStep = -std::numeric_limits<double>::max();
   int firstOrSecond = -1;


#ifdef TIMING
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

   const double maxAlpha = iterate->findBlocking(step, primalValue, primalStep, dualValue, dualStep, firstOrSecond);

   const double mufull = iterate->mustep_pd(step, maxAlpha, maxAlpha) / gamma_a;

   double alpha = 1.0;
   switch (firstOrSecond) {
      case 0:
         alpha = 1; // No constraints were blocking
         break;
      case 1:
         alpha = (-primalValue + mufull / (dualValue + maxAlpha * dualStep)) / primalStep;
#ifdef TIMING
         if( myrank == 0 )
            std::cout << "(primal) original alpha " << alpha << std::endl;
#endif
         break;
      case 2:
         alpha = (-dualValue + mufull / (primalValue + maxAlpha * primalStep)) / dualStep;
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
   alpha = std::min(maxAlpha, alpha);

   // make it at least gamma_f * maxStep
   if (alpha < gamma_f * maxAlpha)
      alpha = gamma_f * maxAlpha;

   // back off just a touch (or a bit more)
   alpha *= steplength_factor;

   assert(alpha < 1.0);

   return alpha;
}

void Solver::finalStepLength_PD(Variables* iterate, Variables* step, double& alpha_primal, double& alpha_dual) {
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


void
Solver::doMonitor(const Problem* data, const Variables* vars, const Residuals* resids, double alpha, double sigma, int i, double mu, int stop_code,
      int level) {
   OoqpMonitor* m = itsMonitors;

   while (m) {
      m->doIt(this, data, vars, resids, alpha, sigma, i, mu, stop_code, level);
      m = m->nextMonitor;
   }
}

void
Solver::doMonitorPd(const Problem* data, const Variables* vars, const Residuals* resids, double alpha_primal, double alpha_dual, double sigma, int i,
      double mu, int stop_code, int level) {
   OoqpMonitor* m = itsMonitors;

   while (m) {
      m->doItPd(this, data, vars, resids, alpha_primal, alpha_dual, sigma, i, mu, stop_code, level);
      m = m->nextMonitor;
   }
}


TerminationCode Solver::doStatus(const Problem* problem, const Variables* vars, const Residuals* resids, int i, double mu, TerminationCode level) {
   if (status) {
      return status->doIt(this, problem, vars, resids, i, mu, level);
   }
   else {
      return this->defaultStatus(problem, vars, resids, i, mu, level);
   }
}


void Solver::monitorSelf() {
   this->addMonitor(new OoqpSelfMonitor);
}

void Solver::addMonitor(OoqpMonitor* m) {
   // Push the monitor onto the list
   m->nextMonitor = itsMonitors;
   itsMonitors = m;
}

std::pair<double, double> Solver::computeUnscaledGapAndResidualNorm(const Residuals& residuals) {
   if (!scaler)
      return std::make_pair(std::fabs(residuals.dualityGap()), residuals.residualNorm());
   else {
      if (!residuals_unscaled)
         residuals_unscaled.reset(scaler->getResidualsUnscaled(residuals));
      else {
         residuals_unscaled->copyFrom(residuals);
         scaler->unscaleResiduals(*residuals_unscaled);
      }

      return std::make_pair(std::fabs(residuals_unscaled->dualityGap()), residuals_unscaled->residualNorm());
   }
}


TerminationCode
Solver::defaultStatus(const Problem*, const Variables* /* vars */, const Residuals* resids, int iterate, double mu, TerminationCode /* level */) {
   const int myrank = PIPS_MPIgetRank();
   TerminationCode stop_code = NOT_FINISHED;
   int idx;

   const std::pair<double, double> gap_norm = computeUnscaledGapAndResidualNorm(*resids);
   const double gap = gap_norm.first;
   const double rnorm = gap_norm.second;

   idx = iterate - 1;
   if (idx < 0)
      idx = 0;
   if (idx >= maxit)
      idx = maxit - 1;

   // store the historical record
   mu_history[idx] = mu;
   rnorm_history[idx] = rnorm;
   phi = (rnorm + gap) / dnorm_orig;
   phi_history[idx] = phi;

   if (idx > 0) {
      phi_min_history[idx] = phi_min_history[idx - 1];
      if (phi < phi_min_history[idx])
         phi_min_history[idx] = phi;
   }
   else
      phi_min_history[idx] = phi;

   if (iterate >= maxit)
      stop_code = MAX_ITS_EXCEEDED;
   else if (mu <= mutol && rnorm <= artol * dnorm_orig)
      stop_code = SUCCESSFUL_TERMINATION;

   if (myrank == 0) {
      std::cout << "mu/mutol: " << mu << "  " << mutol << "  ....   rnorm/limit: " << rnorm << " " << artol * dnorm_orig << std::endl;

      if (printTimeStamp) {
         const double timestamp = MPI_Wtime() - startTime;
         std::cout << "time stamp: " << timestamp << std::endl;
      }
   }

   if (stop_code != NOT_FINISHED)
      return stop_code;

   // check infeasibility condition
   if (idx >= 10 && phi >= 1.e-8 && phi >= 1.e4 * phi_min_history[idx]) {
#ifdef TIMING
      if( myrank == 0 )
         std::cout << "possible INFEASIBLITY detected, phi: " << phi << std::endl;
#endif
      stop_code = INFEASIBLE;
   }

   if (stop_code != NOT_FINISHED)
      return stop_code;

   // check for unknown status: slow convergence first
   if (idx >= 350 && phi_min_history[idx] >= 0.5 * phi_min_history[idx - 30]) {
      stop_code = UNKNOWN;
      printf("hehe dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm_orig, artol);
   }

   if (idx >= 350 && rnorm > artol * dnorm_orig && rnorm_history[idx] * mu_history[0] >= 1.e8 * mu_history[idx] * rnorm_history[0]) {
      stop_code = UNKNOWN;
      printf("dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm_orig, artol);
   }

   //if(dnorm * mu < 50 * rnorm || mu < 1e-5) {
   if (mu * dnorm_orig < 1.0e5 * rnorm) {
      //if(!onSafeSolver) {
      gLackOfAccuracy = 1;
      //cout << "Lack of accuracy detected ---->" << mu << ":" << rnorm/dnorm << endl;
   }
   else {
      //if(dnorm_orig * mu > 1e7 * rnorm && mu > 1.0e4)
      //gLackOfAccuracy=-1;
      //else
      gLackOfAccuracy = 1;
   }
   //onSafeSolver=1;
   //}
   //gLackOfAccuracy=-1; //disable iter refin in sLinsysRootAug
   return stop_code;
}

void Solver::setDnorm(const Problem& data) {
   dnorm = data.datanorm();

   if (scaler)
      dnorm_orig = scaler->getDnormOrig();
   else
      dnorm_orig = dnorm;
}

void Solver::defaultMonitor(const Problem* /* problem */, const Variables* /* vars */, const Residuals* resids, double alpha, double sigma, int i,
      double mu, int status_code, int level) const {
   switch (level) {
      case 0 :
      case 1: {

         const Residuals* resids_unscaled = resids;
         if (scaler)
            resids_unscaled = scaler->getResidualsUnscaled(*resids);

         const double gap = resids_unscaled->dualityGap();
         const double rnorm = resids_unscaled->residualNorm();

         if (scaler)
            delete resids_unscaled;

         std::cout << std::endl << "Duality Gap: " << gap << std::endl;

         if (i > 1) {
            std::cout << " Number of Corrections = " << NumberGondzioCorrections << " alpha = " << alpha << std::endl;
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
   delete[] rnorm_history;
   delete[] phi_history;
   delete[] phi_min_history;
}
