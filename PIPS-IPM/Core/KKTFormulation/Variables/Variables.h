/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef VARIABLES_H
#define VARIABLES_H

#include <cassert>
#include <memory>
#include "Vector.hpp"

class Problem;

class MpsReader;

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

/**
 * Holds the variables used by the interior point solver. In terms of
 * in our abstract problem formulation, these variables are the
 * vectors x, y, z and s.
 *
 * @ingroup AbstractProblemFormulation */
class Variables {
public:

   /** number of complementary primal-dual variables. */
   long long number_complementarity_pairs;

   long long nx, nxupp, nxlow;
   long long my;
   long long mz, mcupp, mclow;

   std::shared_ptr<Vector<double>> ixlow;
   std::shared_ptr<Vector<double>> ixupp;
   std::shared_ptr<Vector<double>> iclow;
   std::shared_ptr<Vector<double>> icupp;

   std::unique_ptr<Vector<double>> primals;
   std::unique_ptr<Vector<double>> slacks;
   std::unique_ptr<Vector<double>> equality_duals;
   std::unique_ptr<Vector<double>> inequality_duals;

   std::unique_ptr<Vector<double>> primal_lower_bound_gap;
   std::unique_ptr<Vector<double>> primal_lower_bound_gap_dual;

   std::unique_ptr<Vector<double>> primal_upper_bound_gap;
   std::unique_ptr<Vector<double>> primal_upper_bound_gap_dual;

   std::unique_ptr<Vector<double>> slack_lower_bound_gap;
   std::unique_ptr<Vector<double>> slack_lower_bound_gap_dual;

   std::unique_ptr<Vector<double>> slack_upper_bound_gap;
   std::unique_ptr<Vector<double>> slack_upper_bound_gap_dual;

   /** constructor in which the data and variable pointers are set to
       point to the given arguments */
   Variables(std::unique_ptr<Vector<double>> x_in, std::unique_ptr<Vector<double>> s_in, std::unique_ptr<Vector<double>> y_in, std::unique_ptr<Vector<double>> z_in, std::unique_ptr<Vector<double>> v_in,
      std::unique_ptr<Vector<double>> gamma_in, std::unique_ptr<Vector<double>> w_in, std::unique_ptr<Vector<double>> phi_in, std::unique_ptr<Vector<double>> t_in, std::unique_ptr<Vector<double>> lambda_in, std::unique_ptr<Vector<double>> u_in,
      std::unique_ptr<Vector<double>> pi_in, std::shared_ptr<Vector<double>> ixlow_in, std::shared_ptr<Vector<double>> ixupp_in, std::shared_ptr<Vector<double>> iclow_in,
      std::shared_ptr<Vector<double>> icupp_in);

   Variables(const Variables& vars);

   /** computes mu = (t'lambda +u'pi + v'gamma + w'phi)/(mclow+mcupp+nxlow+nxupp) */
   [[nodiscard]] double mu() const;

   double mustep_pd(const Variables& step, double alpha_primal, double alpha_dual);

   void saxpy(const Variables& b, double alpha);
   void saxpy_pd(const Variables& b, double alpha_primal, double alpha_dual);

   void negate();

   /** calculate the largest alpha in (0,1] such that the nonnegative
    * variables stay nonnegative in the given search direction. In the
    * general QP problem formulation, this is the largest value of
    * alpha such that (t,u,v,w,lambda,pi,phi,gamma) + alpha *
    * (b->t,b->u,b->v,b->w,b->lambda,b->pi,b->phi,b->gamma) >= 0.
    *
    * @see find_blocking */
   double stepbound(const Variables& iterate);

   /** calculate the largest alpha_primal and alpha_dual in (0,1] such that the nonnegative
    * variables stay nonnegative in the given search direction b. In the
    * abstract problem formulation, this is the largest value of alphas
    * such that (s,z) + alpha_primal * (b->s,0) + alpha_dual * (0,b->z) >= 0.
    *
    * @see stepbound
    */
   std::pair<double, double> stepbound_pd(const Variables& iterate);

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
   double find_blocking(const Variables& step_in, double& primalValue, double& primalStep, double& dualValue, double& dualStep, int& firstOrSecond) const;

   void find_blocking(const Variables& step_in, double& primalValue, double& primalStep, double& dualValue, double& dualStep, double& primalValue_d,
         double& primalStep_d, double& dualValue_d, double& dualStep_d, double& alphaPrimal, double& alphaDual, bool& primalBlocking,
         bool& dualBlocking) const;

   double get_average_distance_to_bound_for_converged_vars(const Problem&, double tol) const;

   void push_slacks_from_bound(double tol, double amount);

   /** sets components of (u,t,v,w) to alpha and of (lambda,pi,phi,gamma) to beta */
   void push_to_interior(double alpha, double beta);

   /** add alpha to components of (u,t,v,w) and beta to components of (lambda,pi,phi,gamma) */
   void shift_bound_variables(double alpha, double beta);

   [[nodiscard]] double violation() const;

   void print() const;
   void print_solution(MpsReader* reader, Problem* problem, int& iErr);

   void unscale_solution(Problem* problem);
   void unscale_bounds(Problem* problem);

   int valid_non_zero_pattern();

   void copy(const Variables& b);

   [[nodiscard]] double one_norm() const;
   [[nodiscard]] double inf_norm() const;

   void set_to_zero();

   virtual ~Variables() = default;
};

#endif


