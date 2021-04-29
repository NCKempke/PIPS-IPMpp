#ifndef RESIDUALS_H
#define RESIDUALS_H

#include <OoqpVectorHandle.h>
#include <iostream>
#include <limits>

/**
 * @file Residuals.h
 * @ingroup AbstractProblemFormulation.
 */

class Problem;

class Variables;

/**
 * Represents the residuals of a QP when solved by an interior point
 * QP solver. In terms of our abstract QP formulation, these residuals
 * are rQ, rA, rC and r3.
 * @ingroup AbstractProblemFormulation
 */
class Residuals {
protected:
   double mResidualNorm{std::numeric_limits<double>::infinity()};
   double mDualityGap{std::numeric_limits<double>::infinity()};
   double primal_objective{std::numeric_limits<double>::infinity()};
   double dual_objective{std::numeric_limits<double>::infinity()};

   long long nx{-1};
   long long my{-1};
   long long mz{-1};

   long long nxupp{-1};
   OoqpVectorHandle ixupp;

   long long nxlow{-1};
   OoqpVectorHandle ixlow;

   long long mcupp{-1};
   OoqpVectorHandle icupp;

   long long mclow{-1};
   OoqpVectorHandle iclow;

   Residuals() = default;

public:
   int m{-1}, n{-1};

   OoqpVectorHandle lagrangian_gradient;
   OoqpVectorHandle rA;
   OoqpVectorHandle rC;
   OoqpVectorHandle rz;
   OoqpVectorHandle rv;
   OoqpVectorHandle rw;
   OoqpVectorHandle rt;
   OoqpVectorHandle ru;
   OoqpVectorHandle rgamma;
   OoqpVectorHandle rphi;
   OoqpVectorHandle rlambda;
   OoqpVectorHandle rpi;

   Residuals(const Residuals& residuals);

   /** The norm of the residuals, ommiting the complementariy conditions */
   double residualNorm() const { return mResidualNorm; }

   /** A quantity that measures progress toward feasibility. IN terms of the abstract problem formulation, this quantity is defined as
    *  @code
    * x' * Q * x - b' * y + c' * x - d' * z
    *  @endcode
    */
   double duality_gap() const { return mDualityGap; };

   double primalObjective() const { return primal_objective; };

   double dualObjective() const { return dual_objective; };

   /** calculate residuals, their norms, and duality/complementarity gap, given a problem and variable set.  */
   void evaluate(Problem& problem, Variables& iterate_in, bool print_resids = false);

   /** Modify the "complementarity" component of the residuals, by
   * adding the pairwise products of the complementary variables plus
   * a constant alpha to this term.  
   */
   void add_to_complementarity_residual(const Variables& variables, double alpha);

   /** Set the "complementarity" component of the residuals to the
    * pairwise products of the complementary variables plus a constant
    * alpha
    */
   void set_complementarity_residual(const Variables& variables, double alpha);

   /** set the complementarity component of the residuals to 0. */
   void clear_complementarity_residual();

   /** set the noncomplementarity components of the residual (the terms
    *  arising from the linear equalities in the KKT conditions) to 0.  */
   void clear_linear_residuals();

   /** perform the projection operation required by Gondzio algorithm:
    * replace each component r3_i of the complementarity component of
    * the residuals by r3p_i - r3_i, where r3p_i is the projection of
    * r3_i onto the box [rmin, rmax]. Then if the resulting value is
    * less than -rmax, replace it by -rmax.
    *
    * @see SimpleVector::gondzioProjection */
   void project_r3(double rmin, double rmax);

   void copy(const Residuals&);

   virtual ~Residuals() = default;

   double recompute_residual_norm();

   int valid_non_zero_pattern();

   void writeToStream(std::ostream& out);

   const long long& getNxupp() { return nxupp; };

   const long long& getNxlow() { return nxlow; };

   const long long& getMcupp() { return mcupp; };

   const long long& getMclow() { return mclow; };

   double constraint_violation();
   double optimality_measure(double mu);
   double feasibility_measure(double mu);
};

#endif



