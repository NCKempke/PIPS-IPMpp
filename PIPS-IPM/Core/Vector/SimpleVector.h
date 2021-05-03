/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SIMPLEVECTOR
#define SIMPLEVECTOR

#include "Vector.hpp"
#include "SmartPointer.h"
#include "pipsdef.h"

#include <vector>

//#define RANGECHECKS

/**
 * Simple sequential vectors with element access.
 * @ingroup SparseLinearAlgebra
 * @ingroup DenseLinearAlgebra
 */
template<typename T>
class SimpleVector : public Vector<T> {
protected:
   int preserveVec;
   T* v;
   std::vector<T> test;

public:
   explicit SimpleVector(int nx = 0);
   SimpleVector(T* v, int nx);
   SimpleVector(const SimpleVector<T>& other);
   //@{
   /**
    * Access the individual elements of this vector.
    */
   T& operator[](int i) {
#ifdef RANGECHECKS
      assert( i >= 0 && i < this->n );
#endif
      return v[i];
   }
   const T& operator[](int i) const {
#ifdef RANGECHECKS
      assert( i >= 0 && i < this->n );
#endif
      return v[i];
   }
   //@}

   Vector<T>* clone() const override;
   /* copy vector entries as well */
   Vector<T>* cloneFull() const override;

   virtual ~SimpleVector();

   void copyIntoArray(T v[]) const override;
   void copyFromArray(const T v[]) override;
   void copyFromArray(const char v[]) override;
   bool isZero() const override;
   void setToZero() override;
   void setToConstant(T c) override;
   void randomize(T, T, T*) override { assert(0 && "not implemented here"); };
   void copyFrom(const Vector<T>& v) override;
   void copyFromAbs(const Vector<T>& v) override;
   double twonorm() const override;
   T infnorm() const override;
   T onenorm() const override;
   void min(T& m, int& index) const override;
   void max(T& m, int& index) const override;
   void absminVecUpdate(Vector<T>& absminvec) const override;
   void absmaxVecUpdate(Vector<T>& absmaxvec) const override;
   void absmin(T& m) const override;
   void absminNonZero(T& m, T zero_eps) const override;
   int getNnzs() const override;

   void componentMult(const Vector<T>& v) override;
   void scalarMult(T num) override;
   void printSolutionToStdErr() const;
   void componentDiv(const Vector<T>& v) override;
   bool componentEqual(const Vector<T>& vec, T tol) const override;
   bool componentNotEqual(const T val, const T tol) const override;
   void setNotIndicatedEntriesToVal(const T val, const Vector<T>& ind) override;

   void writeToStream(std::ostream& out, int offset = 0) const override;
   void writefToStream(std::ostream& out, const char format[]) const override;

   void writeMPSformatOnlyRhs(std::ostream& out, const std::string rowName, const Vector<T>* irhs) const override;
   void writeMPSformatBoundsWithVar(std::ostream& out, const std::string varStub, const Vector<T>* ix, bool upperBound) const override;

   void scale(T alpha) override;

   void axpy(T alpha, const Vector<T>& x) override;
   void axzpy(T alpha, const Vector<T>& x, const Vector<T>& z) override;
   void axdzpy(T alpha, const Vector<T>& x, const Vector<T>& z) override;

   void addConstant(T c) override;

/** perform the projection operation required by Gondzio algorithm:
   * replace each component of the vector v by vp_i - v_i, where vp_i
   * is the projection of v_i onto the box [rmin, rmax]. Then if the
   * resulting value is less than -rmax, replace it by -rmax.
   * */
   void gondzioProjection(T rmin, T rmax) override;
   T dotProductWith(const Vector<T>& v) const override;
   T dotProductSelf(T scaleFactor) const override;

   T shiftedDotProductWith(T alpha, const Vector<T>& mystep, const Vector<T>& yvec, T beta, const Vector<T>& ystep) const override;
   void negate() override;
   void invert() override;
   void invertSave(T zeroReplacementVal = 0.0) override;
   void applySqrt() override;
   void roundToPow2() override;
   bool allPositive() const override;
   bool allOf(const std::function<bool(const T&)>& pred) const override;

   long long numberOfNonzeros() const override;

   bool matchesNonZeroPattern(const Vector<T>& select) const override;
   void selectNonZeros(const Vector<T>& select) override;
   void selectPositive() override;
   void selectNegative() override;
   void addSomeConstants(T c, const Vector<T>& select) override;
   void writefSomeToStream(std::ostream& out, const char format[], const Vector<T>& select) const override;
   void axdzpy(T alpha, const Vector<T>& x, const Vector<T>& z, const Vector<T>& select) override;

   bool isKindOf(int kind) const override;

   bool somePositive(const Vector<T>& select) const override;
   void divideSome(const Vector<T>& div, const Vector<T>& select) override;

   T stepbound(const Vector<T>& v, T maxStep) const override;
   T find_blocking(const Vector<T>& wstep_vec, const Vector<T>& u_vec, const Vector<T>& ustep_vec, T maxStep, T* w_elt, T* wstep_elt, T* u_elt,
         T* ustep_elt, int& first_or_second) const override;

   void find_blocking_pd(const Vector<T>& wstep_vec, const Vector<T>& u_vec, const Vector<T>& ustep_vec, T& maxStepPri, T& maxStepDual, T& w_elt_p,
         T& wstep_elt_p, T& u_elt_p, T& ustep_elt_p, T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d, bool& primalBlocking,
         bool& dualBlocking) const override;

   void removeEntries(const Vector<int>& select) override;

   void permuteEntries(const std::vector<unsigned int>& permvec);

   void getSumCountIfSmall(double tol, double& sum_small, int& n_close, const Vector<T>* select) const override;

   void pushAwayFromZero(double tol, double amount, const Vector<T>* select) override;

   void pushSmallComplementarityPairs(Vector<T>& other_vec, const Vector<T>& select_in, double tol_this, double tol_other, double tol_pairs) override;

   /** Returns a pointer to the elements of this vector. */
   T* elements() const { return v; };

   void appendToFront(unsigned int n_to_add, const T& value);
   void appendToFront(const SimpleVector<T>& other);

   void appendToBack(unsigned int n_to_add, const T& value);
   void appendToBack(const SimpleVector<T>& other);

   void jointCopyFrom(const Vector<T>& vx, const Vector<T>& vy, const Vector<T>& vz) override;
   void jointCopyTo(Vector<T>& vx, Vector<T>& vy, Vector<T>& vz) const override;

   virtual SimpleVector<T>* shaveBorder(int n_shave, bool shave_top);
};

#endif
