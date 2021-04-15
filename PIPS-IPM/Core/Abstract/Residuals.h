#ifndef RESIDUALS_H
#define RESIDUALS_H

#include <OoqpVectorHandle.h>
#include <iostream>

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
   double mResidualNorm;
   double mDualityGap;
   double primal_objective;
   double dual_objective;

   long long nx{0};
   long long my{0};
   long long mz{0};

   long long nxupp{0};
   OoqpVectorHandle ixupp;

   long long nxlow{0};
   OoqpVectorHandle ixlow;

   long long mcupp{0};
   OoqpVectorHandle icupp;

   long long mclow{0};
   OoqpVectorHandle iclow;

   Residuals() = default;

public:
   int m, n;

   OoqpVectorHandle rQ;
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
   double dualityGap() const { return mDualityGap; };

   double primalObjective() const { return primal_objective; };

   double dualObjective() const { return dual_objective; };

   /** calculate residuals, their norms, and duality/complementarity gap, given a problem and variable set.  */
   void evaluate(Problem &problem, Variables *iterate_in, bool print_resids = false);

   /** Modify the "complementarity" component of the residuals, by
   * adding the pairwise products of the complementary variables plus
   * a constant alpha to this term.  
   */
   void add_r3_xz_alpha(const Variables *vars, double alpha);

   /** Set the "complementarity" component of the residuals to the
    * pairwise products of the complementary variables plus a constant
    * alpha
    */
   void set_r3_xz_alpha(const Variables *vars, double alpha);

   /** set the complementarity component of the residuals to 0. */
   void clear_r3();

   /** set the noncomplementarity components of the residual (the terms
    *  arising from the linear equalities in the KKT conditions) to 0.  */
   void clear_r1r2();

   /** perform the projection operation required by Gondzio algorithm:
    * replace each component r3_i of the complementarity component of
    * the residuals by r3p_i - r3_i, where r3p_i is the projection of
    * r3_i onto the box [rmin, rmax]. Then if the resulting value is
    * less than -rmax, replace it by -rmax.
    *
    * @see SimpleVector::gondzioProjection */
   void project_r3(double rmin, double rmax);

   void copyFrom(const Residuals &);

   virtual ~Residuals() = default;

   double recomputeResidualNorm();

   int validNonZeroPattern();

   void writeToStream(std::ostream &out);

   const long long &getNxupp() { return nxupp; };

   const long long &getNxlow() { return nxlow; };

   const long long &getMcupp() { return mcupp; };

   const long long &getMclow() { return mclow; };
};

#endif



