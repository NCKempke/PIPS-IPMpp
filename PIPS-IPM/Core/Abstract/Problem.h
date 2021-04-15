/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef PROBLEM_H
#define PROBLEM_H

#include "QpGenVars.h"
#include "OoqpVector.h"

class Variables;

/**
 * Stores the data of the the QP in a format appropriate for the problem 
 * structure.
 *
 * @ingroup AbstractProblemFormulation
 */
class Problem
{
public:
   virtual ~Problem() = default;

   virtual double objective_value(const QpGenVars* x) const = 0;
   virtual void objective_gradient(const QpGenVars* vars, OoqpVector& gradient) = 0;

    /** compute the norm of the problem data */
    virtual double datanorm() const = 0;

    /** print the problem data */
    virtual void print() = 0;
};

#endif
