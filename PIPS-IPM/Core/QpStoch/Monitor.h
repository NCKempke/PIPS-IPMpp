#ifndef MONITOR_H
#define MONITOR_H

#include "DistributedFactory.h"
#include "Scaler.h"
#include "pipsport.h"

class Solver;

class Problem;

class Variables;

class Residuals;

class Monitor {
public:
   Monitor* nextMonitor{};

   Monitor() = default;
   virtual ~Monitor() = default;
   Monitor(Scaler* scaler = nullptr);
   Monitor(DistributedFactory* qp, Scaler* scaler = nullptr);

   void doIt(const Solver* solver, const Problem* data, const Variables* vars, const Residuals* resids, double alpha, double sigma, int i, double mu,
         int status_code, int level);

   void doItPd(const Solver* solver, const Problem* data, const Variables* vars, const Residuals* resids, double alpha_primal, double alpha_dual,
         double sigma, int i, double mu, int status_code, int level);

private:
   void
   doItStoch(const Solver* solver, const Problem* problem, const Variables* vars, const Residuals* resids, double alpha_primal, double alpha_dual,
         double, int i, double mu, int status_code, int level) const;

protected:
   Scaler* scaler{};

   MPI_Comm mpiComm;
   int myRank{-1};
   int myGlobRank{-1};
};

#endif
