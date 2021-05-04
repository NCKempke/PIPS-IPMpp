/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "Status.h"
#include "Solver.h"

const char* TerminationStrings[] = {"SUCCESSFUL_TERMINATION", "NOT_FINISHED", "MAX_ITS_EXCEEDED", "INFEASIBLE", "UNKNOWN"};


Status::~Status() {
}


CStatus::CStatus(StatusCFunc doItC_, void* ctx_) {
   doItC = doItC_;
   ctx = ctx_;
}

TerminationCode
CStatus::doIt(const Solver* solver, const Problem* qpdata, const Variables* vars, const Residuals* resids, int i, double mu, TerminationCode level) {
   StatusData data;

   data.solver = (void*) solver;
   data.data = (void*) qpdata;
   data.vars = (void*) vars;
   data.resids = (void*) resids;
   data.i = i;
   data.mu = mu;
   data.level = level;
   data.ctx = ctx;
   data.dnorm = solver->dataNorm();

   return doItC((void*) &data);
}
