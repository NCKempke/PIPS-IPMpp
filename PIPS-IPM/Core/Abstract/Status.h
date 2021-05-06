/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef STATUS_H
#define STATUS_H

#include "pipsport.h"

class Solver;

class Problem;

class Variables;

class Residuals;

enum TerminationCode {
   SUCCESSFUL_TERMINATION = 0, NOT_FINISHED, MAX_ITS_EXCEEDED, INFEASIBLE, UNKNOWN
};

extern const char* TerminationStrings[];

/**
 * Class for representing termination conditions of a QP solver.
 *
 * @ingroup QpSolvers
 */
class Status {
public:
   virtual TerminationCode
   doIt(const Solver* solver, const Problem* data, const Variables* vars, const Residuals* resids, int i, double mu, TerminationCode level) = 0;

   virtual ~Status();
};

extern "C" {
typedef TerminationCode (* StatusCFunc)(void* data);
}

#endif










