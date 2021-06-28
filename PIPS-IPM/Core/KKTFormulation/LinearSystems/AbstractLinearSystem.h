/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef ABSTRACTLINEARSYSTEM_H
#define ABSTRACTLINEARSYSTEM_H

class Problem;

class Variables;

class Residuals;

/** Implements the main solver for linear systems that arise in
 * primal-dual interior-point methods for QP
 *
 * @ingroup AbstractProblemFormulation 
 */
class AbstractLinearSystem {
public:
   /** factorizes the matrix, stores data related to the factorization to prepare for later calls to "solve" */
   virtual void factorize(const Variables& iterate) = 0;

   /** assuming the "factor" call was successful, supplies the right-hand side and solves the system. */
   virtual void solve(const Variables& iterate, const Residuals& residuals, Variables& step) = 0;

   virtual ~AbstractLinearSystem() = default;
};


#endif


