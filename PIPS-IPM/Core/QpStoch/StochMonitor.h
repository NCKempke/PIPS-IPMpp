#ifndef STOCH_MONITOR
#define STOCH_MONITOR

#include "OoqpMonitor.h"
#include "DistributedFactory.h"
#include "Scaler.h"
#include "pipsport.h"

class Solver;

class Problem;

class Variables;

class Residuals;

class StochMonitor : public OoqpMonitor {
private:
   void
   doItStoch(const Solver* solver, const Problem* problem, const Variables* vars, const Residuals* resids, double alpha_primal, double alpha_dual,
         double, int i, double mu, int status_code, int level) const;

public:
   StochMonitor(Scaler* scaler = nullptr);
   StochMonitor(DistributedFactory* qp, Scaler* scaler = nullptr);

   void doIt(const Solver* solver, const Problem* data, const Variables* vars, const Residuals* resids, double alpha, double sigma, int i, double mu,
         int status_code, int level) override;

   void doItPd(const Solver* solver, const Problem* data, const Variables* vars, const Residuals* resids, double alpha_primal, double alpha_dual,
         double sigma, int i, double mu, int status_code, int level) override;

protected:
   Scaler* scaler{};

   MPI_Comm mpiComm;
   int myRank{-1};
   int myGlobRank{-1};
};

#endif
