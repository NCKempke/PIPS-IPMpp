/*
 * Scaler.h
 *
 *  Created on: 19.12.2017
 *      Authors: Daniel Rehfeldt, Svenja Uslu
 */

#ifndef SCALER_H
#define SCALER_H

template<typename name>
class Vector;

class Problem;

class Variables;

class Residuals;

/**
 * Abstract base class for scalers.
 */

class Scaler {
protected:
   const bool do_bitshifting; // only scale by power of two factors?
   const bool with_sides; // consider lhs/rhs?
   const double dnorm_orig;

public:
   Scaler(const Problem& problem, bool bitshifting = false, bool usesides = false);
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
   virtual void unscale_residuals(Residuals& residuals) const = 0;

   /** compute original variables from given ones */
   virtual void unscale_variables(Variables& variables) const = 0;

   /** compute original vector from given primal vector */
   virtual Vector<double>* get_primal_unscaled(const Vector<double>& primal_solution) const = 0;

   /** compute original vector from given dual vector */
   virtual Vector<double>* get_dual_eq_unscaled(const Vector<double>& dual_solution) const = 0;

   /** compute original vector from given dual vector */
   virtual Vector<double>* get_dual_ineq_unscaled(const Vector<double>& dual_solution) const = 0;

   /** compute original vector from given dual vector */
   virtual Vector<double>* get_dual_var_bounds_upp_unscaled(const Vector<double>& dual_solution) const = 0;

   /** compute original vector from given dual vector */
   virtual Vector<double>* get_dual_var_bounds_low_unscaled(const Vector<double>& dual_solution) const = 0;
};

#endif /* SCALER_H */
