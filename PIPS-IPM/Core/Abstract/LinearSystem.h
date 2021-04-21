/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

class Problem;

class Variables;

class Residuals;

/** Implements the main solver for linear systems that arise in
 * primal-dual interior-point methods for QP
 *
 * @ingroup AbstractProblemFormulation 
 */
class LinearSystem {
public:
   /** factorizes the matrix, stores data related to the factorization
    * to prepare for later calls to "solve"
    */
   virtual void factorize(Problem* problem, Variables* iterate) = 0;

   /** assuming the "factor" call was successful, supplies the right-hand side and solves the system.
    */
   virtual void solve(Problem* problem, Variables* iterate, Residuals* residuals, Variables* step) = 0;

   virtual ~LinearSystem() = default;
};


#endif


