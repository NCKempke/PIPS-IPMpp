#include <iostream>
#include <cstdio>
#include "Statistics.hpp"
#include "Residuals.h"
#include "Problem.hpp"
#include "Variables.h"
#include "TerminationStatus.hpp"

Statistics::Statistics(Scaler* scaler) : scaler{scaler}, mpi_comm{MPI_COMM_WORLD}, rank{PIPS_MPIgetRank(mpi_comm)}, global_rank{rank} {
}

Statistics::Statistics(const DistributedFactory& factory, const Scaler* scaler) : scaler{scaler}, mpi_comm{factory.tree->getCommWorkers()},
      rank{PIPS_MPIgetRank(mpi_comm)}, global_rank{PIPS_MPIgetRank()} {
}

void Statistics::print(const Problem& problem, const Variables& variables, const Residuals& residuals, double dnorm, double alpha, int i,
      double mu, TerminationStatus status_code, int level) const {
   Statistics::print(problem, variables, residuals, dnorm, alpha, -1.0, i, mu, status_code, level);
}

void Statistics::print(const Problem& problem, const Variables& variables, const Residuals& residuals, double dnorm, double alpha_primal, double
alpha_dual, int i, double mu, TerminationStatus status_code, int level) const {
   double objective = problem.evaluate_objective(variables);

   double residual_norm;
   double duality_gap;

   if (scaler) {
      std::unique_ptr<Residuals> unscaled_residuals = residuals.clone_full();
      objective = scaler->get_unscaled_objective(objective);
      scaler->unscale_residuals(*unscaled_residuals);

      residual_norm = unscaled_residuals->get_residual_norm();
      duality_gap = unscaled_residuals->get_duality_gap();
   } else {
      residual_norm = residuals.get_residual_norm();
      duality_gap = residuals.get_duality_gap();
   }

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
               case TerminationStatus::SUCCESSFUL_TERMINATION:
                  std::cout << "\n *** SUCCESSFUL TERMINATION ***\n";
                  break;
               case TerminationStatus::MAX_ITS_EXCEEDED:
                  std::cout << "\n *** MAXIMUM ITERATIONS REACHED ***\n";
                  break;
               case TerminationStatus::INFEASIBLE:
                  std::cout << "\n *** TERMINATION: PROBABLY INFEASIBLE ***\n";
                  break;
               case TerminationStatus::UNKNOWN:
                  std::cout << "\n *** TERMINATION: STATUS UNKNOWN ***\n";
                  break;
               case TerminationStatus::DID_NOT_RUN:
                  std::cout << "\n *** TERMINATION: DID NOT RUN ***\n";
                  break;
               case TerminationStatus::NOT_FINISHED:
                  std::cout << "\n *** TERMINATION: NOT FINISHED YET ***\n";
                  break;
            }
         }
      }
         break; // end case 0: case 1:
   } // end switch(level)
}
