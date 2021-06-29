#ifndef MONITOR_H
#define MONITOR_H

#include "TerminationStatus.hpp"
#include "DistributedFactory.hpp"
#include "Scaler.hpp"

class Problem;

class Variables;

class Residuals;

class Statistics {
public:
   Statistics() = default;
   virtual ~Statistics() = default;
   explicit Statistics(Scaler* scaler = nullptr);
   explicit Statistics(const DistributedFactory& factory, const Scaler* scaler = nullptr);

   void print(const Problem* problem, const Variables* variables, const Residuals* residuals, double dnorm, double alpha, int i, double mu,
      TerminationStatus status_code, int level) const;

   void print(const Problem* problem, const Variables* variables, const Residuals* residuals, double dnorm, double alpha_primal, double alpha_dual,
         int i, double mu, TerminationStatus status_code, int level) const;

protected:
   const Scaler* scaler{};

   MPI_Comm mpi_comm;
   int rank{-1};
   int global_rank{-1};
};

#endif
