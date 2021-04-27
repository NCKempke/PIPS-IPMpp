/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef OPTIMIZATIONFACTORY
#define OPTIMIZATIONFACTORY

#include <iostream>
#include "OoqpVector_fwd.h"

/**
 *  @defgroup AbstractProblemFormulation
 *
 *  Abstract base classes for defining a problem formulation.
 *
 * A quadratic program QP takes the form
 * @code
 * minimize    c'* x + (1/2) * x' * Q * x
 * subject to  A x  = b
 *             C x >= d
 * @endcode
 *
 * However, for many (possibly most) QP's, the matrices in the
 * formulation have structure that may be exploited to solve the
 * problem more efficiently. The AbstractProblemFormulation module
 * contains abstract base classes upon which these specialized
 * formulations are based. The optimality conditions of the simple QP
 * defined above as are follows:
 *
 * @code
 * rQ  = c + Q * x - A' * y - C' * z = 0
 * rA  = A * x - b                   = 0
 * rC  = C * x - s - d               = 0
 * r3  = S * z                       = 0
 * s, z >= 0
 * @endcode
 *
 * Where rQ, rA, rC and r3 newly defined quantities known as residual
 * vectors and x, y, z and s are variables of used in the solution of
 * the QPs.  
 * @{
 */
class Problem;

class Residuals;

class LinearSystem;

class Variables;

class LinearAlgebraPackage;

/**
 * Creates a compatible set of components representing a problem formulation
 * specialized by structure.
 */
class ProblemFormulation {
public:
   virtual void join_right_hand_side(OoqpVector& rhs_in, const OoqpVector& rhs1_in, const OoqpVector& rhs2_in, const OoqpVector& rhs3_in) const = 0;

   virtual void separate_variables(OoqpVector& x_in, OoqpVector& y_in, OoqpVector& z_in, const OoqpVector& vars_in) const = 0;

   /** create x shaped vector using LinearAlgebraPackage */
   virtual OoqpVector* make_primal_vector() const;
   /** create dual A shaped vector using LinearAlgebraPackage */
   virtual OoqpVector* make_equalities_dual_vector() const;
   /** create dual C shaped vector using LinearAlgebraPackage */
   virtual OoqpVector* make_inequalities_dual_vector() const;
   /** create a rhs vector for the augmented system */
   virtual OoqpVector* make_right_hand_side() const;

   /** create the Residuals class for the relevant formulation */
   virtual Residuals* make_residuals(Problem& problem) = 0;

   /** creates the LinearSystem class for the relevant formulation */
   virtual LinearSystem* make_linear_system(Problem& problem) = 0;

   /** creates the Variables class for the relevant formulation */
   virtual Variables* make_variables(Problem& problem) = 0;

   virtual ~ProblemFormulation() = default;

protected:
   LinearAlgebraPackage* la{};
   /** number of elements in x */
   long long nx{0};

   /** number of rows in A and b including linking rows (sFactory..) */
   long long my{0};

   /** number of rows in C including linking rows */
   long long mz{0};

   ProblemFormulation() = default;
   ProblemFormulation(int nx_, int my_, int mz_);
};

//@}
#endif


