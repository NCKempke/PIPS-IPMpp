#include <iostream>
#include <cstdio>
#include "Statistics.hpp"
#include "Residuals.h"
#include "Problem.h"
#include "Variables.h"
#include "TerminationStatus.h"

Statistics::Statistics(Scaler* scaler) : scaler{scaler}, mpi_comm{MPI_COMM_WORLD}, rank{PIPS_MPIgetRank(mpi_comm)}, global_rank{rank} {
}

Statistics::Statistics(const DistributedFactory& factory, const Scaler* scaler) : scaler{scaler}, mpi_comm{factory.tree->getCommWorkers()},
      rank{PIPS_MPIgetRank(mpi_comm)}, global_rank{PIPS_MPIgetRank()} {
}

void Statistics::print(const Problem* problem, const Variables* variables, const Residuals* residuals, double dnorm, double alpha, double sigma, int i,
      double mu, int status_code, int level) {
   Statistics::print(problem, variables, residuals, dnorm, alpha, -1.0, sigma, i, mu, status_code, level);
}

void
Statistics::print(const Problem* problem, const Variables* variables, const Residuals* residuals, double dnorm, double alpha_primal, double alpha_dual,
      double sigma, int i, double mu, int status_code, int level) const {
   double objective = problem->objective_value(*variables);

   const Residuals* unscaled_residuals = residuals;
   if (scaler) {
      objective = scaler->get_unscaled_objective(objective);
      unscaled_residuals = scaler->get_unscaled_residuals(*residuals);
   }

   const double residual_norm = unscaled_residuals->residual_norm();
   const double duality_gap = unscaled_residuals->duality_gap();

   if (scaler)
      delete unscaled_residuals;

   // log only on the first proc
   if (rank > 0)
      return;

   switch (level) {
      case 0:
      case 1: {
         std::cout << " --- Iteration " << i << " --- (rank " << global_rank << ")" << "\n";
         if (i == 1)
            printf(" mu = %16.12e  rel.res.norm=%16.12e  datanorm=%16.12e\n", mu, residual_norm / dnorm, dnorm);
         else
            printf(" mu = %16.12e  rel.res.norm=%16.12e\n", mu, residual_norm / dnorm);
         //cout << " mu = " << mu << " relative residual norm = "
         //cout << resids->residualNorm() / dnorm << "\n";
         std::cout << " Duality Gap:  " << duality_gap << "\n";
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
            // Termination has been detected by the status check; print appropriate message
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
