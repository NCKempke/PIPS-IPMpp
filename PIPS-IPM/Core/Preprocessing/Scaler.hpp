/*
 * QpScaler.h
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

#ifndef SCALER_HPP
#define SCALER_HPP

#include <memory>

template<typename T>
class Vector;

class Problem;

class Variables;

class Residuals;

class ProblemFactory;

class GeneralMatrix;

/**
 * Abstract base class for scalers.
 */
class Scaler {
public:

   explicit Scaler(const ProblemFactory& problem_factory, const Problem& problem, bool bitshifting = false, bool usesides = false);
   virtual ~Scaler() = default;

   /** return norm of unscaled problem */
   [[nodiscard]] double getDnormOrig() const { return dnorm_orig; }
   /** scale */
   virtual void scale() = 0;

   [[nodiscard]] double get_unscaled_objective(double objval) const;

   void unscale_variables(Variables& variables) const;
   void unscale_residuals(Residuals& residuals) const;

   [[nodiscard]] Vector<double>* get_primal_unscaled(const Vector<double>& primal_solution) const;
   [[nodiscard]] Vector<double>* get_dual_eq_unscaled(const Vector<double>& dual_solution) const;
   [[nodiscard]] Vector<double>* get_dual_ineq_unscaled(const Vector<double>& dual_solution) const;
   [[nodiscard]] Vector<double>* get_dual_var_bounds_upp_unscaled(const Vector<double>& dual_solution) const;
   [[nodiscard]] Vector<double>* get_dual_var_bounds_low_unscaled(const Vector<double>& dual_solution) const;

protected:
   static void invertAndRound(bool round, Vector<double>& vector);

   void create_scaling_vectors();
   void applyScaling() const;
   virtual void doObjScaling() const = 0;

   /** get maximum absolute row ratio and write maximum row entries into vectors */
   double
   maxRowRatio(Vector<double>& maxvecA, Vector<double>& maxvecC, Vector<double>& minvecA, Vector<double>& minvecC, const Vector<double>* colScalevec) const;

   /** get maximum absolute column ratio and write maximum column entries into vectors */
   double maxColRatio(Vector<double>& maxvec, Vector<double>& minvec, const Vector<double>* rowScaleVecA, const Vector<double>* rowScaleVecC) const;

   void scaleObjVector(double scaling_factor);

   /** print row col ratios */
   void printRowColRatio() const;

   void setScalingVecsToOne();

   // problem factory to create scaling vectors
   const ProblemFactory& problem_factory;

   // scaling vectors
   std::unique_ptr<Vector<double>> vec_rowscaleQ{};
   std::unique_ptr<Vector<double>> vec_rowscaleA{};
   std::unique_ptr<Vector<double>> vec_rowscaleC{};
   std::unique_ptr<Vector<double>> vec_colscale{};

   // problem data
   std::shared_ptr<GeneralMatrix> A;
   std::shared_ptr<GeneralMatrix> C;
   std::shared_ptr<Vector<double>> obj;
   std::shared_ptr<Vector<double>> bA;
   std::shared_ptr<Vector<double>> bux;
   std::shared_ptr<Vector<double>> blx;
   std::shared_ptr<Vector<double>> rhsC;
   std::shared_ptr<Vector<double>> lhsC;

   // scaling factor for objective
   double factor_objscale;

   const double dnorm_orig;

   const bool do_bitshifting; // only scale by power of two factors?
   const bool with_sides; // consider lhs/rhs?
   const bool scaling_output{false};

   // has scaling been applied
   mutable bool scaling_applied{false};
};

#endif /* SCALER_HPP */
