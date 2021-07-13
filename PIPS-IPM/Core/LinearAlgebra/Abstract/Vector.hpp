/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */
/* Modified by Cosmin Petra and Miles Lubin */


#ifndef VECTOR_H
#define VECTOR_H

#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <functional>
#include "../../Utilities/pipsdef.h"

/** An abstract class representing the implementation of a OoqpVectorTemplate.
 *
 *  Do not create instances of Vector. Create instance of subclasses
 *  of Vector instead.
 *
 *  @ingroup AbstractLinearAlgebra
 */
template<typename T>
class Vector : public std::enable_shared_from_this<Vector<T>> {
protected:
   int n{0};
public:
   /** Return the length of this vector. */
   int length() const { return n; }

   Vector(int n_);

   Vector() = default;
   virtual ~Vector() = default;

   virtual void pushAwayFromZero(double /*tol*/, double /*amount*/, const Vector<T>* /*select */ ) = 0;

   /** Set all elements of this Vector<double> to zero. */
   virtual void setToZero() = 0;
   /** Check if all elemens in the vector are equal to zero. */
   virtual bool isZero() const = 0;
   /** Set all elements of this Vector<double> to the constant value c */
   virtual void setToConstant(T c) = 0;
   /** Copy the elements of v into this Vector<double> object. */
   virtual void copyFrom(const Vector<T>& v) = 0;

   /** Return the infinity norm of this Vector<double> object. */
   [[nodiscard]] virtual double two_norm() const = 0;
   /** Return the infinity norm of this Vector<double> object. */
   virtual T inf_norm() const = 0;
   /** Return the one norm of this Vector<double> object. */
   virtual T one_norm() const = 0;

   /** Return number of elements in this vector not considered zero */
   [[nodiscard]] virtual int getNnzs() const = 0;

   /** Multiply the components of this Vector<double> by the components of v. */
   virtual void componentMult(const Vector<T>& v) = 0;
   /** Divide the components of this Vector<double> by the components of v. */
   virtual void componentDiv(const Vector<T>& v) = 0;
   /** Compare two vectors with tolerance tol and tell if the are different */
   virtual bool componentEqual(const Vector<T>& v, T tol) const = 0;
   /** Assert that elements in the vector are unequal value */
   virtual bool componentNotEqual(const T val, const T tol = pips_eps0) const = 0;
   /** Multiply the components of this Vector<double> by num. */
   virtual void scalarMult(T num) = 0;

   /** sets entries for which in[i] = zero to val */
   virtual void setNotIndicatedEntriesToVal(const T val, const Vector<T>& ind) = 0;

   /** Write the components of this Vector<double>, one element per line, to
    *  the stream out.
    *
    *  Offsets are used to visualize the StochVector structure (if any)
    */
   virtual void write_to_stream(std::ostream& out, int offset = 0) const = 0;

   /** Scale each element of this Vector<double> by the constant alpha */
   virtual void scale(T alpha) = 0;

   /** this += alpha * x */
   virtual void add(T alpha, const Vector<T>& x) = 0;
   /** this += alpha * x * z */
   virtual void add_product(T alpha, const Vector<T>& x, const Vector<T>& z) = 0;
   /** this += alpha * x / z */
   virtual void add_quotient(T alpha, const Vector<T>& x, const Vector<T>& z) = 0;

   /** Add c to the elements of this Vector<double> object */
   virtual void add_constant(T c) = 0;

   /** Perform the projection needed by Gondzio's multiple corrector method.
    *
    * @see SimpleVector::gondzioProjection
    */
   virtual void gondzioProjection(T rmin, T rmax) = 0;

   /** Return the minimum value in this vector, and the index at
    *  which it occurs. */
   virtual void min(T& m, int& index) const = 0;

   /** Return the maximum value in this vector, and the index at
    *  which it occurs. */
   virtual void max(T& m, int& index) const = 0;

   /** Return the absolute minimum value of this vector */
   virtual void absmin(T& m) const = 0;

   virtual void absminVecUpdate(Vector<T>& absminvec) const = 0;

   virtual void absmaxVecUpdate(Vector<T>& absmaxvec) const = 0;

   /** Return the absolute minimum value of this vector bigger than eps_zero, or -1 if none could be found */
   virtual void absminNonZero(T& m, T zero_eps) const = 0;

   /** Return the dot product of this Vector<double> with v */
   virtual T dotProductWith(const Vector<T>& v) const = 0;

   /** Return the scaled dot product of this (scaled) vector with itself  */
   virtual T dotProductSelf(T scaleFactor) const = 0;

   [[nodiscard]] virtual T scaled_dot_product_self(const Vector<T>& scale) const = 0;

   /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
    */
   virtual T shiftedDotProductWith(T alpha, const Vector<T>& mystep, const Vector<T>& yvec, T beta, const Vector<T>& ystep) const = 0;
   /** Negate all the elements of this Vector<double> object. */
   virtual void negate() = 0;

   /** Invert (1/x) the elements of this Vector<double>. */
   virtual void invert() = 0;

   /** Invert (1/x) the elements of this Vector<double>, but don't divide by zero and replace by zeroReplacementVal instead */
   virtual void safe_invert(T zero_replacement_value = 0.0) = 0;

   /** Take the square root of each element of this Vector<double>. */
   virtual void sqrt() = 0;

   /** Rounds vector entries to nearest power of two values */
   virtual void roundToPow2() = 0;

   /** True if all elements of this Vector<double> are positive. */
   virtual bool all_positive() const = 0;

   virtual void transform(const std::function<T(const T&)>& transformation) = 0;

   virtual T sum_reduce(const std::function<T(const T& a, const T& b)>& reduce) const = 0;

   virtual bool all_of(const std::function<bool(const T&)>& pred) const = 0;

   /** Return the number of non-zero elements in this Vector<double>. */
   virtual long long number_nonzeros() const = 0;

   /** True if this Vector<double> has the same non-zero pattern as select. */
   virtual bool matchesNonZeroPattern(const Vector<T>& select) const = 0;

   /** Set each element of this Vector<double> to zero if the corresponding
    *  element in select is zero.
    */
   virtual void selectNonZeros(const Vector<T>& select) = 0;

   /** Set all negative elements to zero.
    */
   virtual void selectPositive() = 0;

   /** Set all positive elements to zero.
    */
   virtual void selectNegative() = 0;

   /** Add the constant c to some of the elements of this Vector<double>
    *  @param c The constant to be added
    *  @param mask a mask Vector<double>. The constant c is added to an element
    *         of this Vector<double> only if the corresponding element of mask is
    *         non-zero.
    */
   virtual void add_constant(T c, const Vector<T>& mask) = 0;

   /** this += alpha * x / z
    *  @param select only perform the division on elements of x and z if the
    *         corresponding element of select is non-zero. Generally we avoid
    *         performing the division if we know that it will result in
    *         division by zero. The Vector<double> select may be x, z or a third
    *         Vector<double>.
    */
   virtual void add_quotient(T alpha, const Vector<T>& x, const Vector<T>& z, const Vector<T>& select) = 0;

   /** True if selected elements of this Vector<double> are positive
    *  @param select Each element of this Vector<double> must be positive
    *                if the corresponding element of select is non-zero.
    */
   virtual bool are_positive(const Vector<T>& select) const = 0;
   /** Divide selected elements of this Vector<double> by the corresponding
    *  element in div.
    *  @param div If element i of this Vector<double> is selected, then it will
    *             be divided by element i of div.
    *  @param select Perform division on elements of this Vector<double> only if
    *             the corresponding element in select is non-zero
    */
   virtual void divideSome(const Vector<T>& div, const Vector<T>& select) = 0;

   /** True if this Vector<double> identifies itself as having the type kind.
    *
    *  Classes overriding this method must call the inherited version, so that
    *  the class hierarchy is supported.
    */
   virtual bool isKindOf(int kind) const = 0;

   /** Return the largest value of alpha in the interval [0, maxStep] for
    *  which: this + alpha * v >= (1 - fraction)*this). Set firstBlocking to be the
    *  "blocking" index i - the one that limits the step length to alpha.
    */
   virtual T fraction_to_boundary(const Vector<T>& step_in, T fraction = T{1}) const = 0;

   /** Return the largest value of alpha in the interval [0,1] for
    *  which: this + alpha * wstep_vec >= 0 and u_vec + alpha *
    *  ustep_vec >= 0. Also return the components this[i],
    *  wstep_vec[i], u_vec[i], ustep_vec[i], where i is the index of
    *  the "blocking" component - the one that limits the step length
    *  to alpha. Also return first_or_second=1 if the blocking
    *  component is in "this", and first_or_second=2 if the blocking
    *  component is in u_vec.
    */
   virtual T
   find_blocking(const Vector<T>& wstep_vec, const Vector<T>& u_vec, const Vector<T>& ustep_vec, T maxStep, T* w_elt, T* wstep_elt, T* u_elt,
         T* ustep_elt, int& first_or_second) const = 0;

   virtual void
   find_blocking_pd(const Vector<T>& wstep_vec, const Vector<T>& u_vec, const Vector<T>& ustep_vec, T& maxStepPri, T& maxStepDual, T& w_elt_p,
         T& wstep_elt_p, T& u_elt_p, T& ustep_elt_p, T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d, bool& primalBlocking,
         bool& dualBlocking) const = 0;

   /** Copy the elements of this Vector<double> into the C-style array v. */
   virtual void copyIntoArray(T v[]) const = 0;
   /** Copy the elements of the C-style array v into this Vector<double>. */
   virtual void copyFromArray(const T v[]) = 0;
   /** Copy the elements of the C-style char array v into this Vector<double>. */
   virtual void copyFromArray(const char v[]) = 0;

   /** remove entries i for which select[i] == 0 */
   virtual void removeEntries(const Vector<int>&) = 0;

   /** Copy the absolute values of elements of v_in into this Vector<double> object. */
   virtual void copyFromAbs(const Vector<T>& v) = 0;

   /** compute the sum of all entries smaller considered zero with tol and count how many*/
   virtual void getSumCountIfSmall(double, double&, int&, const Vector<T>*) const = 0;

   virtual void pushSmallComplementarityPairs(Vector<T>& /*other_vec*/, const Vector<T>& /*select_in*/, double /*tol_this*/, double /*tol_other*/,
         double /*tol_pairs*/) = 0;

   virtual Vector<T>* clone() const = 0;

   virtual Vector<T>* clone_full() const = 0;

   /** copy [this] = [vx, vy, vz] */
   virtual void jointCopyFrom(const Vector<T>& vx, const Vector<T>& vy, const Vector<T>& vz) = 0;
   /** copy [vx, vy, vz] = [this] */
   virtual void jointCopyTo(Vector<T>& vx, Vector<T>& vy, Vector<T>& vz) const = 0;

   [[nodiscard]] virtual double barrier_directional_derivative(const Vector<T>& x, const Vector<T>& bound, const Vector<T>& bound_indicator) const = 0;
   [[nodiscard]] virtual double barrier_directional_derivative(const Vector<T>& x, double bound, const Vector<T>& bound_indicator) const = 0;

   /*** find and the return the set a[i] b[i] c[i] d[i] where a[i]/b[i] + c[i]/d[i] is min/max ***/
   [[nodiscard]] virtual std::tuple<double, double, double, double> find_abs_nonzero_max_min_pair_a_by_b_plus_c_by_d(const Vector<T>& a,
      const Vector<T>& b, const Vector<T>& select_ab, bool use_ab, const Vector<T>& c, const Vector<T>& d, const Vector<T>& select_cd, bool use_cd, bool find_min) const = 0;
};

// TODO : should be a proper enum class..
enum {
   kSimpleVector = 0, kStochVector, kStochDummy,
};

#endif
