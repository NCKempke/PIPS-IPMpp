/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef VARIABLES_H
#define VARIABLES_H

#include <cassert>
#include "OoqpVectorHandle.h"

class QpGen;

class LinearAlgebraPackage;

class Problem;

/** Indicates what type is the blocking variable in the step length
 * determination. If tblock, then the blocking variable is one of the
 * slack variables t for a general lower bound, and so on. Special
 * value no_block is for the case in which a step length of 1 can be
 * taken without hitting the bound.  */

enum {
   no_block = 0, t_block = 1, lambda_block = 2, u_block = 3, pi_block = 4, v_block = 5, gamma_block = 6, w_block = 7, phi_block = 8
};

/**
 * @file Variables.h
 * @ingroup AbstractProblemFormulation
 */

class LinearSystem;

class Residuals;

class Solver;

class MpsReader;

#ifdef TESTING
class VariablesTester;
#endif

/**
 * Holds the variables used by the interior point solver. In terms of
 * in our abstract problem formulation, these variables are the
 * vectors x, y, z and s.
 *
 * @ingroup AbstractProblemFormulation */
class Variables {
public:
#ifdef TESTING
   friend VariablesTester;
#endif

   /** number of complementary primal-dual variables. */
   long long nComplementaryVariables;

   long long nx, nxupp, nxlow;
   long long my;
   long long mz, mcupp, mclow;

   OoqpVectorHandle ixlow;
   OoqpVectorHandle ixupp;
   OoqpVectorHandle icupp;
   OoqpVectorHandle iclow;
   Variables() {};

   OoqpVectorHandle x;
   OoqpVectorHandle s;
   OoqpVectorHandle y;
   OoqpVectorHandle z;

   OoqpVectorHandle v;
   OoqpVectorHandle gamma;

   OoqpVectorHandle w;
   OoqpVectorHandle phi;

   OoqpVectorHandle t;
   OoqpVectorHandle lambda;

   OoqpVectorHandle u;
   OoqpVectorHandle pi;

   /** constructor in which the data and variable pointers are set to
       point to the given arguments */
   Variables(OoqpVector* x_in, OoqpVector* s_in, OoqpVector* y_in, OoqpVector* z_in, OoqpVector* v_in, OoqpVector* gamma_in, OoqpVector* w_in,
         OoqpVector* phi_in, OoqpVector* t_in, OoqpVector* lambda_in, OoqpVector* u_in, OoqpVector* pi_in, OoqpVector* ixlow_in, OoqpVector* ixupp_in,
         OoqpVector* iclow_in, OoqpVector* icupp_in);

   /** constructor that creates variables objects of specified
       dimensions. */
   Variables(LinearAlgebraPackage* la, long long nx_, long long my_, long long mz_, OoqpVector* ixlow, OoqpVector* ixupp, OoqpVector* iclow,
         OoqpVector* icupp);

   Variables(const Variables& vars);

   double getAverageDistanceToBoundForConvergedVars(const Problem&, double tol) const;

   void pushSlacksFromBound(double tol, double amount);

   /** computes mu = (t'lambda +u'pi + v'gamma + w'phi)/(mclow+mcupp+nxlow+nxupp) */
   double mu();

   double mustep_pd(const Variables* step, double alpha_primal, double alpha_dual);

   void saxpy(const Variables* b, double alpha);
   void saxpy_pd(const Variables* b, double alpha_primal, double alpha_dual);

   void negate();

   /** calculate the largest alpha in (0,1] such that the nonnegative
    * variables stay nonnegative in the given search direction. In the
    * general QP problem formulation, this is the largest value of
    * alpha such that (t,u,v,w,lambda,pi,phi,gamma) + alpha *
    * (b->t,b->u,b->v,b->w,b->lambda,b->pi,b->phi,b->gamma) >= 0.
    *
    * @see findBlocking */
   double stepbound(const Variables* b);

   /** calculate the largest alpha_primal and alpha_dual in (0,1] such that the nonnegative
    * variables stay nonnegative in the given search direction b. In the
    * abstract problem formulation, this is the largest value of alphas
    * such that (s,z) + alpha_primal * (b->s,0) + alpha_dual * (0,b->z) >= 0.
    *
    * @see stepbound
    */
   void stepbound_pd(const Variables* b, double& alpha_primal, double& alpha_dual);

   /** Performs the same function as stepbound, and supplies additional
    * information about which component of the nonnegative variables is
    * responsible for restricting alpha. In terms of the abstract
    * formulation, the components have the following meanings.
    *
    * @param primalValue the value of the blocking component of the
    * primal variables (u,t,v,w).
    *
    * @param primalStep the corresponding value of the blocking
    * component of the primal step variables (b->u,b->t,b->v,b->w).
    *
    * @param dualValue the value of the blocking component of the dual
    * variables (lambda,pi,phi,gamma).
    *
    * @param dualStep the corresponding value of the blocking component
    * of the dual step variables (b->lambda,b->pi,b->phi,b->gamma).
    *
    * @param firstOrSecond  1 if the primal step is blocking, 2 if the dual
    * step is block, 0 if no step is blocking.
    *
    * @see stepbound
    * */
   double findBlocking(const Variables* step, double& primalValue, double& primalStep, double& dualValue, double& dualStep, int& firstOrSecond);

   void findBlocking_pd(const Variables* step, double& primalValue, double& primalStep, double& dualValue, double& dualStep, double& primalValue_d,
         double& primalStep_d, double& dualValue_d, double& dualStep_d, double& alphaPrimal, double& alphaDual, bool& primalBlocking,
         bool& dualBlocking);

   /** sets components of (u,t,v,w) to alpha and of (lambda,pi,phi,gamma) to beta */
   void push_to_interior(double alpha, double beta);

   /** add alpha to components of (u,t,v,w) and beta to components of
       (lambda,pi,phi,gamma) */
   void shiftBoundVariables(double alpha, double beta);

   /** check whether this is an interior point. Useful as a sanity check. */
   int isInteriorPoint();

   double violation();

   void print();
   void printSolution(MpsReader* reader, Problem* problem, int& iErr);

   void unscaleSolution(Problem* problem);
   void unscaleBounds(Problem* problem);

   int validNonZeroPattern();

   void copy(const Variables* b);

   double onenorm() const;
   double infnorm() const;

   void setToZero();

   void printNorms() const;

   void setNotIndicatedBoundsTo(Problem& data, double value);

   virtual ~Variables() = default;
};

#endif


