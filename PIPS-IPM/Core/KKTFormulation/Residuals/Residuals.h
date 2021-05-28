#ifndef RESIDUALS_H
#define RESIDUALS_H

#include "Vector.hpp"
#include <memory>
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

class Residuals {
protected:
   double residual_norm{std::numeric_limits<double>::infinity()};
   double duality_gap{std::numeric_limits<double>::infinity()};
   double primal_objective{std::numeric_limits<double>::infinity()};
   double dual_objective{std::numeric_limits<double>::infinity()};

   long long nx{-1};
   long long my{-1};
   long long mz{-1};

   std::shared_ptr<Vector<double>> ixupp;
   long long nxupp{-1};

   std::shared_ptr<Vector<double>> ixlow;
   long long nxlow{-1};

   std::shared_ptr<Vector<double>> icupp;
   long long mcupp{-1};

   std::shared_ptr<Vector<double>> iclow;
   long long mclow{-1};

   Residuals() = default;

public:
   std::unique_ptr<Vector<double>> lagrangian_gradient;
   std::unique_ptr<Vector<double>> equality_residuals;
   std::unique_ptr<Vector<double>> inequality_residuals;
   std::unique_ptr<Vector<double>> inequality_dual_residuals;
   std::unique_ptr<Vector<double>> rv;
   std::unique_ptr<Vector<double>> rw;
   std::unique_ptr<Vector<double>> rt;
   std::unique_ptr<Vector<double>> ru;
   std::unique_ptr<Vector<double>> rgamma;
   std::unique_ptr<Vector<double>> rphi;
   std::unique_ptr<Vector<double>> rlambda;
   std::unique_ptr<Vector<double>> rpi;

   Residuals(std::unique_ptr<Vector<double>> rQ_, std::unique_ptr<Vector<double>> rA_,
      std::unique_ptr<Vector<double>> rC_, std::unique_ptr<Vector<double>> rz_, std::unique_ptr<Vector<double>> rt_,
      std::unique_ptr<Vector<double>> rlambda_, std::unique_ptr<Vector<double>> ru_,
      std::unique_ptr<Vector<double>> rpi_,
      std::unique_ptr<Vector<double>> rv_, std::unique_ptr<Vector<double>> rgamma_, std::unique_ptr<Vector<double>> rw_,
      std::unique_ptr<Vector<double>> rphi_, std::shared_ptr<Vector<double>> ixlow_,
      std::shared_ptr<Vector<double>> ixupp_, std::shared_ptr<Vector<double>> iclow_,
      std::shared_ptr<Vector<double>> icupp_);

   Residuals(const Residuals& residuals);

   /** The norm of the residuals, ommiting the complementarity conditions */
   [[nodiscard]] double get_residual_norm() const { return residual_norm; }

   /** A quantity that measures progress toward feasibility. IN terms of the abstract problem formulation, this quantity is defined as
    *  @code
    * x' * Q * x + c' * x - b' * y - d' * z
    *  @endcode
    */
   [[nodiscard]] double get_duality_gap() const { return duality_gap; };

   [[nodiscard]] double get_primal_objective() const { return primal_objective; };

   [[nodiscard]] double get_dual_objective() const { return dual_objective; };

   /** calculate residuals, their norms, and duality/complementarity gap, given a problem and variable set.  */
   void evaluate(Problem& problem, Variables& iterate_in, bool print_residuals = false);

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

   void write_to_stream(std::ostream& out);

   [[nodiscard]] const long long& getNxupp() const { return nxupp; };

   [[nodiscard]] const long long& getNxlow() const { return nxlow; };

   [[nodiscard]] const long long& getMcupp() const { return mcupp; };

   [[nodiscard]] const long long& getMclow() const { return mclow; };

   double constraint_violation();

   double optimality_measure(double mu);

   double feasibility_measure(double mu);
};

#endif



