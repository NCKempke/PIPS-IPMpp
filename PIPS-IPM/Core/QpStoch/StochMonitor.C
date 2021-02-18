#include "StochMonitor.h"
#include "Status.h"

#include "Residuals.h"
#include "Solver.h"
#include "QpGenData.h"
#include "QpGenVars.h"
#include "sTree.h"
#include <iostream>
#include <cstdio>
#include "pipsport.h"

StochMonitor::StochMonitor(QpGenStoch* qp_, Scaler* scaler)
  : qp(qp_), scaler(scaler)
{
  mpiComm=MPI_COMM_WORLD; //default for old version
  MPI_Comm_rank(mpiComm, &myRank);
  myGlobRank = myRank;
}

StochMonitor::StochMonitor(sFactory* qp_, Scaler* scaler)
  : qp(nullptr), scaler(scaler)
{
  mpiComm = qp_->tree->getCommWorkers();
  MPI_Comm_rank(mpiComm, &myRank);
  MPI_Comm_rank(MPI_COMM_WORLD, &myGlobRank);
}

void StochMonitor::doIt( const Solver * solver, const Data * data, const Variables * vars,
			 const Residuals * resids,
			 double alpha, double sigma,
			 int i, double mu,
			 int status_code,
			 int level )
{
   StochMonitor::doItStoch(solver, data, vars, resids, alpha, -1.0, sigma, i, mu, status_code, level);
}

void StochMonitor::doItPd( const Solver * solver, const Data * data, const Variables * vars,
              const Residuals * resids,
              double alpha_primal, double alpha_dual, double sigma,
              int i, double mu,
              int status_code,
              int level )
{
   StochMonitor::doItStoch(solver, data, vars, resids, alpha_primal, alpha_dual, sigma, i, mu, status_code, level);
}

void
StochMonitor::doItStoch(const Solver *solver, const Data *data,
      const Variables *vars, const Residuals *resids, double alpha_primal,
      double alpha_dual, double, int i, double mu, int status_code,
      int level) const
{
   const double objective_scaled = dynamic_cast<const QpGenData*>(data)->objectiveValue(
         dynamic_cast<const QpGenVars*>(vars));
   const double gap_scaled = resids->dualityGap();
   const Residuals *resids_unscaled = resids;

   if( scaler )
      resids_unscaled = scaler->getResidualsUnscaled(*resids);

   const double dnorm = solver->dataNormOrig();
   const double rnorm = resids_unscaled->residualNorm();

   const double gap_unscaled = resids_unscaled->dualityGap();
   const double objective_unscaled = scaler ? scaler->getObjUnscaled(objective_scaled) : objective_scaled;

   if( scaler )
      delete resids_unscaled;

   // log only on the first proc
   if( myRank > 0 )
      return;

   switch( level )
      {
      case 0:
      case 1:
         {
            std::cout << " --- Iteration " << i << " --- (rank " << myGlobRank
                  << ")" << "\n";
            if( i == 1 )
               printf(" mu = %16.12e  rel.res.norm=%16.12e  datanorm=%16.12e\n",
                     mu, rnorm / dnorm, dnorm);
            else
               printf(" mu = %16.12e  rel.res.norm=%16.12e\n", mu,
                     rnorm / dnorm);
            //cout << " mu = " << mu << " relative residual norm = "
            //cout << resids->residualNorm() / dnorm << "\n";
            std::cout << " Duality Gap (unscaled/scaled):  " << gap_unscaled << "/" << gap_scaled << "\n";
            std::cout << " Relaitve Duality Gap (unscaled/scaled):  " << gap_unscaled / objective_unscaled << "/" << gap_scaled / objective_scaled << "\n";
            if( i > 1 )
            {
               if( alpha_dual != -1.0 )
               {
                  std::cout << " alpha primal = " << alpha_primal << "\n";
                  std::cout << " alpha dual = " << alpha_dual << "\n";
               }
               else
                  std::cout << " alpha = " << alpha_primal << "\n";
            }
            std::cout << " Objective (unscaled/scaled): " << objective_unscaled << "/" << objective_scaled << "\n";
            std::cout << "\n";
            if( level == 1 )
            {
               // Termination has been detected by the status check; print
               // appropriate message
               switch( status_code )
                  {
                  case SUCCESSFUL_TERMINATION:
                     std::cout << "\n" << " *** SUCCESSFUL TERMINATION ***"
                           << "\n";
                     break;
                  case MAX_ITS_EXCEEDED:
                     std::cout << "\n" << " *** MAXIMUM ITERATIONS REACHED *** "
                           << "\n";
                     break;
                  case INFEASIBLE:
                     std::cout << "\n"
                           << " *** TERMINATION: PROBABLY INFEASIBLE *** "
                           << "\n";
                     break;
                  case UNKNOWN:
                     std::cout << "\n"
                           << " *** TERMINATION: STATUS UNKNOWN *** " << "\n";
                     break;
                  } // end switch(statusCode)
            }
         }
         break; // end case 0: case 1:
      } // end switch(level)
}
