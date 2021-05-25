/*
 * Scaler.h
 *
 *  Created on: 19.12.2017
 *      Authors: Daniel Rehfeldt, Svenja Uslu
 */

#ifndef PIPS_IPM_CORE_ABSTRACT_SCALER_H_
#define PIPS_IPM_CORE_ABSTRACT_SCALER_H_

template<typename name>
class Vector;

class Problem;

class Variables;

class Residuals;

/**  * @defgroup Preprocessing
 *
 * Interior-point scalers
 * @{
 */

/**
 * Abstract base class for scalers.
 */


class Scaler {
protected:
   Problem* const problem;
   const bool do_bitshifting; // only scale by power of two factors?
   const bool with_sides; // consider lhs/rhs?

   const double dnorm_orig;

public:
   Scaler(Problem* prob, bool bitshifting = false, bool usesides = false);
   virtual ~Scaler() = default;

   /** scale */
   virtual void scale() = 0;

   /** return norm of unscaled problem */
   virtual double getDnormOrig() const { return dnorm_orig; }

   /** unscale given objective value */
   virtual double get_unscaled_objective(double objective_value) const = 0;

   /** compute original variables from given ones */
   virtual Variables* get_unscaled_variables(const Variables& variables) const = 0;

   /** compute original residuals from given ones */
   virtual Residuals* get_unscaled_residuals(const Residuals& residuals) const = 0;

   /** compute original residuals from given ones */
   virtual void unscaleResiduals(Residuals& residuals) const = 0;

   /** compute original variables from given ones */
   virtual void unscaleVariables(Variables& variables) const = 0;

   /** compute original vector from given primal vector */
   virtual Vector<double>* getPrimalUnscaled(const Vector<double>& primal_solution) const = 0;

   /** compute original vector from given dual vector */
   virtual Vector<double>* getDualEqUnscaled(const Vector<double>& dual_solution) const = 0;

   /** compute original vector from given dual vector */
   virtual Vector<double>* getDualIneqUnscaled(const Vector<double>& dual_solution) const = 0;

   /** compute original vector from given dual vector */
   virtual Vector<double>* getDualVarBoundsUppUnscaled(const Vector<double>& dual_solution) const = 0;

   /** compute original vector from given dual vector */
   virtual Vector<double>* getDualVarBoundsLowUnscaled(const Vector<double>& dual_solution) const = 0;
};

//@}


#endif /* PIPS_IPM_CORE_ABSTRACT_SCALER_H_ */
