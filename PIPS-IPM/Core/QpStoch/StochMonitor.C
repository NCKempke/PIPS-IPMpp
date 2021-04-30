#include "StochMonitor.h"
#include "Status.h"

#include "Residuals.h"
#include "Solver.h"
#include "Problem.h"
#include "Variables.h"
#include <iostream>
#include <cstdio>

StochMonitor::StochMonitor(Scaler* scaler) : scaler{scaler}, mpiComm{MPI_COMM_WORLD}, myRank{PIPS_MPIgetRank(mpiComm)}, myGlobRank{myRank} {
}

StochMonitor::StochMonitor(DistributedFactory* qp, Scaler* scaler) : scaler{scaler}, mpiComm{qp->tree->getCommWorkers()},
      myRank{PIPS_MPIgetRank(mpiComm)}, myGlobRank{PIPS_MPIgetRank()} {
}

void StochMonitor::doIt(const Solver* solver, const Problem* data, const Variables* vars, const Residuals* resids, double alpha, double sigma, int i,
      double mu, int status_code, int level) {
   StochMonitor::doItStoch(solver, data, vars, resids, alpha, -1.0, sigma, i, mu, status_code, level);
}

void StochMonitor::doItPd(const Solver* solver, const Problem* data, const Variables* vars, const Residuals* resids, double alpha_primal,
      double alpha_dual, double sigma, int i, double mu, int status_code, int level) {
   StochMonitor::doItStoch(solver, data, vars, resids, alpha_primal, alpha_dual, sigma, i, mu, status_code, level);
}

void StochMonitor::doItStoch(const Solver* solver, const Problem* problem, const Variables* vars, const Residuals* resids, double alpha_primal,
      double alpha_dual, double, int i, double mu, int status_code, int level) const {
   double objective = problem->objective_value(*vars);

   const Residuals* resids_unscaled = resids;
   if (scaler) {
      objective = scaler->getObjUnscaled(objective);
      resids_unscaled = scaler->getResidualsUnscaled(*resids);
   }

   const double dnorm = solver->dataNormOrig();
   const double rnorm = resids_unscaled->residualNorm();
   const double gap = resids_unscaled->duality_gap();

   if (scaler)
      delete resids_unscaled;

   // log only on the first proc
   if (myRank > 0)
      return;

   switch (level) {
      case 0:
      case 1: {
         std::cout << " --- Iteration " << i << " --- (rank " << myGlobRank << ")" << "\n";
         if (i == 1)
            printf(" mu = %16.12e  rel.res.norm=%16.12e  datanorm=%16.12e\n", mu, rnorm / dnorm, dnorm);
         else
            printf(" mu = %16.12e  rel.res.norm=%16.12e\n", mu, rnorm / dnorm);
         //cout << " mu = " << mu << " relative residual norm = "
         //cout << resids->residualNorm() / dnorm << "\n";
         std::cout << " Duality Gap:  " << gap << "\n";
         if (i > 1) {
            if (alpha_dual != -1.0) {
               std::cout << " alpha primal = " << alpha_primal << "\n";
               std::cout << " alpha dual = " << alpha_dual << "\n";
            }
            else
               std::cout << " alpha = " << alpha_primal << "\n";
         }
         std::cout << " Objective: " << objective << "\n";
         std::cout << "\n";
         if (level == 1) {
            // Termination has been detected by the status check; print
            // appropriate message
            switch (status_code) {
               case SUCCESSFUL_TERMINATION:
                  std::cout << "\n *** SUCCESSFUL TERMINATION ***\n";
                  break;
               case MAX_ITS_EXCEEDED:
                  std::cout << "\n *** MAXIMUM ITERATIONS REACHED ***\n";
                  break;
               case INFEASIBLE:
                  std::cout << "\n *** TERMINATION: PROBABLY INFEASIBLE ***\n";
                  break;
               case UNKNOWN:
                  std::cout << "\n *** TERMINATION: STATUS UNKNOWN ***\n";
                  break;
            }
         }
      }
         break; // end case 0: case 1:
   } // end switch(level)
}
