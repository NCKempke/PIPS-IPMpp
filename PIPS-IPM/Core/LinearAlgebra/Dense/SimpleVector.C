/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "VectorUtilities.h"
#include "SimpleVector.h"
#include "OoqpBlas.h"
#include <pipsdef.h>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <algorithm>
#include <functional>
#include <memory>

template<typename T>
SimpleVector<T>::SimpleVector(int n_) : Vector<T>(n_), test(n_) {
   preserveVec = 0;
   v = new T[this->n];

   std::uninitialized_fill(v, v + this->n, T{});
}

template<typename T>
long long SimpleVector<T>::number_nonzeros() const {
   long long i, count = 0;
   for (i = 0; i < this->n; i++) {
      if (v[i] != (T) 0.0)
         count++;
   }
   return count;
}

template<typename T>
void SimpleVector<T>::min(T& m, int& index) const {
   if (this->n == 0) {
      m = std::numeric_limits<T>::max();
      return;
   }
   index = 0;
   m = v[0];
   for (int i = 0; i < this->n; i++) {
      if (v[i] < m) {
         m = v[i];
         index = i;
      }
   }
}

template<typename T>
void SimpleVector<T>::absminVecUpdate(Vector<T>& absminvec) const {
   SimpleVector<T>& absminvecSimple = dynamic_cast<SimpleVector<T>&>(absminvec);
   assert(absminvecSimple.length() == this->n);
   T* const absminvecArr = absminvecSimple.elements();

   for (int i = 0; i < this->n; i++) {
      const T abs = fabs(v[i]);
      if (abs < absminvecArr[i] && abs > pips_eps)
         absminvecArr[i] = abs;
   }
}

template<typename T>
void SimpleVector<T>::absmaxVecUpdate(Vector<T>& absmaxvec) const {
   SimpleVector<T>& absmaxvecSimple = dynamic_cast<SimpleVector<T>&>(absmaxvec);
   assert(absmaxvecSimple.length() == this->n);
   T* const absmaxvecArr = absmaxvecSimple.elements();

   for (int i = 0; i < this->n; i++) {
      const T abs = std::abs(v[i]);
      if (abs > absmaxvecArr[i])
         absmaxvecArr[i] = abs;
   }
}

template<typename T>
void SimpleVector<T>::absmin(T& m) const {
   m = std::numeric_limits<T>::infinity();

   if (this->n == 0)
      return;

   for (int i = 0; i < this->n; ++i)
      m = std::min(m, std::abs(v[i]));
}

/** Compute the min absolute value that is larger than zero_eps.
 * If there is no such value, return inf */
template<typename T>
void SimpleVector<T>::absminNonZero(T& m, T zero_eps) const {
   assert(zero_eps >= 0.0);

   m = std::numeric_limits<T>::infinity();

   if (this->n == 0)
      return;

   for (int i = 0; i < this->n; i++) {
      if (std::abs(v[i]) < m && std::abs(v[i]) > zero_eps)
         m = std::abs(v[i]);
   }
}

template<typename T>
int SimpleVector<T>::getNnzs() const {
   int non_zeros = 0;
   for (int i = 0; i < this->n; i++) {
      if (!PIPSisZero(v[i]))
         non_zeros++;
   }

   return non_zeros;
}

template<typename T>
void SimpleVector<T>::max(T& m, int& index) const {
   m = -std::numeric_limits<T>::max();
   index = -1;

   for (int i = 0; i < this->n; i++) {
      if (v[i] > m) {
         m = v[i];
         index = i;
      }
   }
}

template<typename T>
bool SimpleVector<T>::isKindOf(int kind) const {
   return (kind == kSimpleVector);
}

template<typename T>
void SimpleVector<T>::copyIntoArray(T w[]) const {
   std::copy(v, v + this->n, w);
}

template<typename T>
void SimpleVector<T>::copyFromArray(const T w[]) {
   std::copy(w, w + this->n, v);
}

template<typename T>
void SimpleVector<T>::copyFromArray(const char w[]) {
   for (int i = 0; i < this->n; i++) {
      this->v[i] = w[i];
   }
}

template<typename T>
void SimpleVector<T>::pushAwayFromZero(double tol, double amount, const Vector<T>* select) {
   assert(0 < amount);
   assert(0 < tol);

   const SimpleVector<T>* selects = dynamic_cast<const SimpleVector<T>*>(select);

   if (selects)
      assert(this->n == selects->n);

   for (int i = 0; i < this->n; ++i) {
      if (selects && (*selects)[i] == 0)
         continue;

      if (PIPSisZero(v[i], tol))
         v[i] += amount;
   }
}

template<typename T>
SimpleVector<T>::SimpleVector(const SimpleVector<T>& other)
      : SimpleVector<T>(other.n) {
   std::copy(other.v, other.v + this->n, v);
}

template<typename T>
SimpleVector<T>::SimpleVector(T* v_, int n_)
      : Vector<T>(n_) {
   preserveVec = 1;
   v = v_;
}

template<typename T>
SimpleVector<T>::~SimpleVector() {
   if (!preserveVec) {
      delete[] v;
   }
}

template<typename T>
Vector<T>* SimpleVector<T>::cloneFull() const {
   SimpleVector<T>* clone = new SimpleVector<T>(this->n);
   clone->copyFromArray(v);

   return clone;
}

template<typename T>
Vector<T>* SimpleVector<T>::clone() const {
   return new SimpleVector<T>(this->n);
}

template<typename T>
bool SimpleVector<T>::isZero() const {
   bool is_zero = true;

   for (int i = 0; i < this->n; ++i)
      is_zero = (is_zero && v[i] == 0.0);

   return is_zero;
}

template<typename T>
void SimpleVector<T>::setToZero() {
   setToConstant(0.0);
}

template<typename T>
void SimpleVector<T>::setToConstant(T c) {
   std::fill(v, v + this->n, c);
}

template<typename T>
void SimpleVector<T>::copyFrom(const Vector<T>& vec) {
   assert(vec.length() == this->n);

   vec.copyIntoArray(this->v);
}

template<typename T>
void SimpleVector<T>::copyFromAbs(const Vector<T>& vec) {
   const auto& vecSimple = dynamic_cast<const SimpleVector<T>&>(vec);
   assert(vec.length() == this->n);
   T* const vecArr = vecSimple.elements();

   for (int i = 0; i < this->n; i++)
      v[i] = std::abs(vecArr[i]);
}

template<typename T>
T SimpleVector<T>::inf_norm() const {

   if (this->n == 0)
      return -std::numeric_limits<T>::max();

   T temp, norm = 0;
   for (int i = 0; i < this->n; i++) {
      temp = std::abs(v[i]);
      // Subtle reversal of the logic to handle NaNs
      if (!(temp <= norm))
         norm = temp;
   }

   return norm;
}

template<>
double SimpleVector<double>::inf_norm() const {
   if (this->n == 0)
      return -std::numeric_limits<double>::max();

   const int one = 1;
   return std::fabs(v[idamax_(&this->n, v, &one) - 1]);
}

template<typename T>
T SimpleVector<T>::one_norm() const {
   T temp, norm = 0;
   for (int i = 0; i < this->n; i++) {
      temp = std::abs(v[i]);
      norm += temp;
   }
   return norm;
}

template<typename T>
double SimpleVector<T>::two_norm() const {
   T temp = dotProductWith(*this);
   return std::sqrt(temp);
}

template<typename T>
void SimpleVector<T>::componentMult(const Vector<T>& vec) {
   assert(this->n == vec.length());
   const auto& sv = dynamic_cast<const SimpleVector<T>&>(vec);
   const T* y = sv.v;
   for (int i = 0; i < this->n; i++)
      v[i] *= y[i];
}

template<typename T>
bool SimpleVector<T>::componentEqual(const Vector<T>& vec, T tol) const {
   assert(this->n == vec.length());
   const SimpleVector<T>& sv = dynamic_cast<const SimpleVector<T>&>(vec);

   for (int i = 0; i < this->n; ++i) {
      /* two comparisons - a numerical one and one for stuff like infinity/nan/max/min */
      if (!PIPSisRelEQ(v[i], sv[i], tol) && v[i] != sv[i]) {
//         std::cout << v[i] << " != " << sv[i] << "\n";
         return false;
      }
   }
   return true;
}

template<typename T>
bool SimpleVector<T>::componentNotEqual(const T val, const T tol) const {
   for (int i = 0; i < this->n; ++i) {
      /* two comparisons - a numerical one and one for stuff like infinity/nan/max/min */
      if (PIPSisRelEQ(v[i], val, tol) || v[i] == val) {
         return false;
      }
   }
   return true;
}

template<typename T>
void SimpleVector<T>::setNotIndicatedEntriesToVal(const T val, const Vector<T>& ind) {
   const SimpleVector<T>& ind_vec = dynamic_cast<const SimpleVector<T>&>(ind);
   assert(ind_vec.length() == this->length());

   for (int i = 0; i < ind_vec.length(); ++i) {
      if (ind_vec[i] == 0)
         this->v[i] = val;
   }
}


template<typename T>
void SimpleVector<T>::scalarMult(T num) {
   for (int i = 0; i < this->n; i++)
      v[i] *= num;
}

// Print first 10 entries of solution vector to stderr.
// Useful for debugging purposes...
template<typename T>
void SimpleVector<T>::printSolutionToStdErr() const {
   for (int i = 0; i < std::min(10, this->n); i++) {
      std::cerr << v[i] << "\n";
   }
   std::cerr << "******" << "\n";
}

template<typename T>
void SimpleVector<T>::componentDiv(const Vector<T>& vec) {
   assert(this->n == vec.length());
   T* pv = v, * lv = v + this->n;

   const auto& sv = dynamic_cast<const SimpleVector<T>&>(vec);
   const T* y = sv.v;

   for (; pv < lv; pv++, y++)
      *pv /= *y;
}

template<typename T>
void SimpleVector<T>::writeToStream(std::ostream& out, int offset) const {
   for (int i = 0; i < this->n; i++) {
      for (int j = 0; j < offset; ++j)
         out << "\t";
      out << v[i] << "\n";
   }
}

template<typename T>
void SimpleVector<T>::writefToStream(std::ostream& out, const char format[]) const {
   std::unique_ptr<SimpleVector<T> > empty = std::make_unique<SimpleVector<T>>(0);
   this->writefSomeToStream(out, format, *empty);
}

template<typename T>
void SimpleVector<T>::writefSomeToStream(std::ostream& out, const char format[], const Vector<T>& select) const {
   const auto& sselect = dynamic_cast<const SimpleVector<T>&>(select);
   T* s = 0;
   if (select.length() > 0) {
      s = sselect.v;
   }
   for (int i = 0; i < this->n; i++) {
      if (!s || s[i] != 0.0) {
         int j = 0;
         char c;
         while ((c = format[j]) != 0) {
            if (c != '%') {
               out << c;
            }
            else {
               // Brain-dead variable substitution, but good enough for this
               // simple case
               if (0 == strncmp("{value}", &format[j + 1], 7)) {
                  out << v[i];
                  j += 7;
               }
               else if (0 == strncmp("{index}", &format[j + 1], 7)) {
                  out << i;
                  j += 7;
               }
               else {
                  out << c;
               }
            }
            j++;
         }
         out << "\n";
      }
   }
}

template<>
void SimpleVector<double>::scale(double alpha) {
   if (this->n == 0)
      return;

   const int one = 1;
   dscal_(&this->n, &alpha, v, &one);
}

// generic implementation without boost 
template<typename T>
void SimpleVector<T>::scale(T) {
   assert(0 && "not implemented here");
   // std::transform( this->v, this->v + this->n, this->v, [alpha](T a)->T { return alpha * a; } );
}

template<>
void SimpleVector<double>::axpy(double alpha, const Vector<double>& vec) {
   assert(this->n == vec.length());
   if (this->n == 0)
      return;

   const auto& sv = dynamic_cast<const SimpleVector<double>&>(vec);
   const int one = 1;
   daxpy_(&this->n, &alpha, sv.v, &one, v, &one);
}

template<typename T>
void SimpleVector<T>::axpy(T alpha, const Vector<T>& vec) {
   assert(this->n == vec.length());
   const auto& sv = dynamic_cast<const SimpleVector<T>&>(vec);

#ifndef PRE_CPP11
   std::transform(this->v, this->v + this->n, sv.v, this->v, [alpha](T a, T b) -> T { return a + alpha * b; });
#else
   for(int i = 0; i < this->n; ++i)
     this->v[i] += alpha*sv[i];
#endif
}

template<typename T>
void SimpleVector<T>::addConstant(T c) {
   for (int i = 0; i < this->n; i++)
      v[i] += c;
}

template<typename T>
void SimpleVector<T>::gondzioProjection(T rmin, T rmax) {
   for (int i = 0; i < this->n; i++) {
      if (v[i] < rmin) {
         v[i] = rmin - v[i];
      }
      else if (v[i] > rmax) {
         v[i] = rmax - v[i];
      }
      else {
         v[i] = 0.0;
      }

      if (v[i] < -rmax)
         v[i] = -rmax;
   }
}

template<typename T>
void SimpleVector<T>::axzpy(T alpha, const Vector<T>& xvec, const Vector<T>& zvec) {
   assert(this->n == xvec.length() && this->n == zvec.length());

   const auto& sxvec = dynamic_cast<const SimpleVector<T>&>(xvec);
   const auto& szvec = dynamic_cast<const SimpleVector<T>&>(zvec);

   T* x = sxvec.v;
   T* z = szvec.v;
   T* lx = x + this->n;
   T* w = v;

   if (alpha == 1.0) {
      while (x < lx) {
         *w += *x * *z;
         w++;
         x++;
         z++;
      }
   }
   else if (alpha == -1) {
      while (x < lx) {
         *w -= *x * *z;
         w++;
         x++;
         z++;
      }
   }
   else {
      while (x < lx) {
         *w += alpha * *x * *z;
         w++;
         x++;
         z++;
      }
   }
}

template<typename T>
void SimpleVector<T>::axdzpy(T alpha, const Vector<T>& xvec, const Vector<T>& zvec) {
   const SimpleVector<T>& sxvec = dynamic_cast<const SimpleVector<T>&>(xvec);
   T* x = sxvec.v;
   const auto& szvec = dynamic_cast<const SimpleVector<T>&>(zvec);
   T* z = szvec.v;

   assert(this->n == xvec.length() && this->n == zvec.length());

   for (int i = 0; i < this->n; i++) {
      //if(x[i] > 0 && z[i] > 0)
      assert(z[i] != 0);
      v[i] += alpha * x[i] / z[i];
   }
}

template<typename T>
void SimpleVector<T>::axdzpy(T alpha, const Vector<T>& xvec, const Vector<T>& zvec, const Vector<T>& select) {
   assert(this->n == xvec.length() && this->n == zvec.length());

   const auto& sxvec = dynamic_cast<const SimpleVector<T>&>(xvec);
   T* x = sxvec.v;
   const auto& szvec = dynamic_cast<const SimpleVector<T>&>(zvec);
   T* z = szvec.v;
   const auto& sselect = dynamic_cast<const SimpleVector<T>&>(select);
   T* s = sselect.v;
   if (alpha == 1.0) {
      for (int i = 0; i < this->n; i++) {
         if (0.0 != s[i])
            v[i] += x[i] / z[i];
      }
   }
   else if (alpha == -1.0) {
      for (int i = 0; i < this->n; i++) {
         if (0.0 != s[i])
            v[i] -= x[i] / z[i];
      }
   }
   else {
      for (int i = 0; i < this->n; i++) {
         if (0.0 != s[i])
            v[i] += alpha * x[i] / z[i];
      }
   }
}

template<>
double SimpleVector<double>::dotProductWith(const Vector<double>& vec) const {
   assert(this->n == vec.length());
   if (this->n == 0)
      return 0.0;

   const SimpleVector<double>& svec = dynamic_cast<const SimpleVector<double>&>(vec);

   const int incx = 1;
   return ddot_(&this->n, v, &incx, svec.v, &incx);
}

template<typename T>
T SimpleVector<T>::dotProductWith(const Vector<T>& vec) const {
   assert(this->n == vec.length());
   const SimpleVector<T>& svec = dynamic_cast<const SimpleVector<T>&>(vec);

   T* vvec = svec.v;

   T dot1 = 0.0;
   T dot2 = 0.0;

   const int size = 8196;
   int kmax = this->n / size;
   int i = 0;

   int k, imax;
   for (k = 0; k < kmax; k++) {
      imax = (k + 1) * size;
      for (; i < imax; i++) {
         dot1 += v[i] * vvec[i];
      }
      dot2 += dot1;
      dot1 = 0;
   }
   for (; i < this->n; i++) {
      dot1 += v[i] * vvec[i];
   }

   return dot2 + dot1;
}

template<typename T>
T SimpleVector<T>::dotProductSelf(T scaleFactor) const {
   assert(scaleFactor >= 0.0);
   return scaleFactor*scaleFactor*dotProductWith(*this);
}

template<typename T>
T SimpleVector<T>::shiftedDotProductWith(T alpha, const Vector<T>& mystep, const Vector<T>& yvec, T beta, const Vector<T>& ystep) const {
   assert(this->n == mystep.length() && this->n == yvec.length() && this->n == ystep.length());

   const SimpleVector<T>& syvec = dynamic_cast<const SimpleVector<T>&>(yvec);
   T* y = syvec.v;

   const SimpleVector<T>& smystep = dynamic_cast<const SimpleVector<T>&>(mystep);
   T* p = smystep.v;

   const SimpleVector<T>& systep = dynamic_cast<const SimpleVector<T>&>(ystep);
   T* q = systep.v;

   T dot1 = 0.0;
   T dot2 = 0.0;

   const int size = 8196;
   int kmax = this->n / size;
   int i = 0;

   int k, imax;
   for (k = 0; k < kmax; k++) {
      imax = (k + 1) * 8196;
      for (; i < imax; i++) {
         dot1 += (v[i] + alpha * p[i]) * (y[i] + beta * q[i]);
      }
      dot2 += dot1;
      dot1 = 0;
   }
   for (; i < this->n; i++) {
      dot1 += (v[i] + alpha * p[i]) * (y[i] + beta * q[i]);
   }

   return dot2 + dot1;
}

template<typename T>
void SimpleVector<T>::negate() {
   for (int i = 0; i < this->n; i++)
      v[i] = -v[i];
}

template<typename T>
void SimpleVector<T>::invert() {
   for (int i = 0; i < this->n; i++) {
      assert(v[i] != 0.0);
      v[i] = 1.0 / v[i];
   }
}

template<typename T>
void SimpleVector<T>::safe_invert(T zero_replacement_value) {
   for (int i = 0; i < this->n; i++) {
      if (v[i] != 0.0)
         v[i] = 1 / v[i];
      else
         v[i] = zero_replacement_value;
   }
}

template<typename T>
void SimpleVector<T>::sqrt() {
   for (int i = 0; i < this->n; i++) {
      assert(v[i] >= 0.0);
      v[i] = std::sqrt(v[i]);
   }
}

template<typename T>
void SimpleVector<T>::roundToPow2() {
   for (int i = 0; i < this->n; i++) {
      int exp;
#if 0
      const double mantissa = std::frexp(v[i], &exp);

      if( mantissa >= 0.75 )
         v[i] = std::ldexp(0.5, exp + 1);
      else
         v[i] = std::ldexp(0.5, exp);
#else
      (void) std::frexp(v[i], &exp);
      v[i] = std::ldexp(0.5, exp);
#endif
   }
}

template<typename T>
bool SimpleVector<T>::all_positive() const {
   for (int i = 0; i < this->n; i++) {
      if (v[i] <= 0)
         return false;
   }
   return true;
}

template<typename T>
bool SimpleVector<T>::all_of(const std::function<bool(const T&)>& pred) const {
   return std::all_of(v, v + this->n, pred);
}


template<typename T>
T SimpleVector<T>::stepbound(const Vector<T>& pvec, T maxStep) const {
   assert(this->n == pvec.length());

   const SimpleVector<T>& spvec = dynamic_cast<const SimpleVector<T>&>(pvec);
   T* p = spvec.v;
   T* w = v;
   T bound = maxStep;

   for (int i = 0; i < this->n; i++) {
      T temp = p[i];
      if (w[i] >= 0 && temp < 0) {
         temp = -w[i] / temp;
         if (temp < bound) {
            bound = temp;
         }
      }
   }
   return bound;
}

template<typename T>
T SimpleVector<T>::find_blocking(const Vector<T>& wstep_vec, const Vector<T>& u_vec, const Vector<T>& ustep_vec, T maxStep, T* w_elt, T* wstep_elt,
      T* u_elt, T* ustep_elt, int& first_or_second) const {
   T* w = v;
   const SimpleVector<T>& swstep = dynamic_cast<const SimpleVector<T>&>(wstep_vec);
   T* wstep = swstep.v;

   const SimpleVector<T>& su_vec = dynamic_cast<const SimpleVector<T>&>(u_vec);
   T* u = su_vec.v;

   const SimpleVector<T>& sustep_vec = dynamic_cast<const SimpleVector<T>&>(ustep_vec);
   T* ustep = sustep_vec.v;

   return ::find_blocking(w, this->n, 1, wstep, 1, u, 1, ustep, 1, maxStep, w_elt, wstep_elt, u_elt, ustep_elt, first_or_second);
}

template<typename T>
void SimpleVector<T>::find_blocking_pd(const Vector<T>& wstep_vec, const Vector<T>& u_vec, const Vector<T>& ustep_vec, T& maxStepPri, T& maxStepDual,
      T& w_elt_p, T& wstep_elt_p, T& u_elt_p, T& ustep_elt_p, T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d, bool& primalBlocking,
      bool& dualBlocking) const {
   const T* w = v;
   const SimpleVector<T>& swstep = dynamic_cast<const SimpleVector<T>&>(wstep_vec);
   const T* wstep = swstep.v;

   const SimpleVector<T>& su_vec = dynamic_cast<const SimpleVector<T>&>(u_vec);
   const T* u = su_vec.v;

   const SimpleVector<T>& sustep_vec = dynamic_cast<const SimpleVector<T>&>(ustep_vec);
   const T* ustep = sustep_vec.v;

   ::find_blocking_pd(w, this->n, wstep, u, ustep, maxStepPri, maxStepDual, w_elt_p, wstep_elt_p, u_elt_p, ustep_elt_p, w_elt_d, wstep_elt_d, u_elt_d,
         ustep_elt_d, primalBlocking, dualBlocking);
}

template<typename T>
bool SimpleVector<T>::matchesNonZeroPattern(const Vector<T>& select) const {
   const SimpleVector<T>& sselect = dynamic_cast<const SimpleVector<T>&>(select);
   T* map = sselect.v;

   T* lmap = map + this->n;
   assert(this->n == select.length());

   T* w = v;
   while (map < lmap) {
      if (*map == 0.0 && *w != 0.0)
         return false;
      map++;
      w++;
   }

   return true;
}

template<typename T>
void SimpleVector<T>::selectNonZeros(const Vector<T>& select) {
   const SimpleVector<T>& sselect = dynamic_cast<const SimpleVector<T>&>(select);
   T* map = sselect.v;

   assert(this->n == select.length());
   for (int i = 0; i < this->n; i++) {
      if (map[i] == 0.)
         v[i] = 0.0;
   }
}

template<typename T>
void SimpleVector<T>::selectPositive() {
   // TODO : std::transform and lambda for C++11
   for (int i = 0; i < this->n; i++) {
      if (v[i] <= 0.0)
         v[i] = 0.0;
   }
}

template<typename T>
void SimpleVector<T>::selectNegative() {
   // TODO : std::transform and lambda for C++11
   for (int i = 0; i < this->n; i++) {
      if (v[i] >= 0.0)
         v[i] = 0.0;
   }
}

template<typename T>
void SimpleVector<T>::add_constant(T c, const Vector<T>& select) {
   const SimpleVector<T>& sselect = dynamic_cast<const SimpleVector<T>&>(select);
   T* map = sselect.v;

   assert(this->n == select.length());
   for (int i = 0; i < this->n; i++) {
      if (map[i])
         v[i] += c;
   }
}

template<typename T>
bool SimpleVector<T>::somePositive(const Vector<T>& select) const {
   const SimpleVector<T>& sselect = dynamic_cast<const SimpleVector<T>&>(select);
   T* map = sselect.v;

   assert(this->n == select.length());
   for (int i = 0; i < this->n; i++) {
      if (0.0 != map[i] && v[i] <= 0) {
         std::cout << "Element " << i << " is nonpositive: " << v[i] << "\n";
         return false;
      }
   }
   return true;
}

template<typename T>
void SimpleVector<T>::divideSome(const Vector<T>& div, const Vector<T>& select) {
   if (this->n == 0)
      return;

   const SimpleVector<T>& sselect = dynamic_cast<const SimpleVector<T>&>(select);
   T* map = sselect.v;

   const SimpleVector<T>& sdiv = dynamic_cast<const SimpleVector<T>&>(div);
   T* q = sdiv.v;
   assert(this->n == div.length() && this->n == select.length());

   for (int i = 0; i < this->n; i++) {
      if (0.0 != map[i]) {
         assert(!PIPSisZero(map[i]));
         assert(!PIPSisZero(q[i]));

         v[i] /= q[i];
      }
   }
}

template<typename T>
void SimpleVector<T>::removeEntries(const Vector<int>& select) {
   const SimpleVector<int>& selectSimple = dynamic_cast<const SimpleVector<int>&>(select);
   const int* const selectArr = selectSimple.elements();
   assert(selectArr);

   assert(this->n == selectSimple.length());

   int nNew = 0;

   for (int i = 0; i < this->n; i++)
      if (selectArr[i] != 0)
         v[nNew++] = v[i];

   this->n = nNew;
}

template<typename T>
void SimpleVector<T>::permuteEntries(const std::vector<unsigned int>& permvec) {
   if (this->n == 0)
      return;

   assert(this->n > 0);
   assert(permvec.size() == size_t(this->n));

   T* buffer = new T[this->n];

   for (size_t i = 0; i < permvec.size(); i++) {
      assert(permvec[i] < unsigned(this->n));
      buffer[i] = v[permvec[i]];
   }

   std::copy(buffer, buffer + this->n, v);

   delete[] buffer;
}

template<typename T>
void SimpleVector<T>::appendToFront(unsigned int n_to_add, const T& value) {
   assert(!preserveVec);

   const int new_len = this->n + n_to_add;

   T* new_v = new T[new_len];

   std::uninitialized_fill(new_v, new_v + n_to_add, value);
   std::uninitialized_copy(this->v, this->v + this->n, new_v + n_to_add);

   delete[] this->v;

   this->n = new_len;
   this->v = new_v;
}

template<typename T>
void SimpleVector<T>::appendToFront(const SimpleVector<T>& other) {
   assert(!preserveVec);

   const int new_len = this->n + other.n;

   T* new_v = new T[new_len];

   std::uninitialized_copy(other.v, other.v + other.n, new_v);
   std::uninitialized_copy(this->v, this->v + this->n, new_v + other.n);

   delete[] this->v;

   this->n = new_len;
   this->v = new_v;
}

template<typename T>
void SimpleVector<T>::jointCopyFrom(const Vector<T>& vx, const Vector<T>& vy, const Vector<T>& vz) {
   assert(this->length() == vx.length() + vy.length() + vz.length());

   const SimpleVector<T>& x = dynamic_cast<const SimpleVector<T>&>(vx);
   const SimpleVector<T>& y = dynamic_cast<const SimpleVector<T>&>(vy);
   const SimpleVector<T>& z = dynamic_cast<const SimpleVector<T>&>(vz);

   std::copy(x.v, x.v + x.length(), v);
   std::copy(y.v, y.v + y.length(), v + x.length());
   std::copy(z.v, z.v + z.length(), v + x.length() + y.length());
}

template<typename T>
void SimpleVector<T>::jointCopyTo(Vector<T>& vx, Vector<T>& vy, Vector<T>& vz) const {
   assert(this->length() == vx.length() + vy.length() + vz.length());

   SimpleVector<T>& x = dynamic_cast<SimpleVector<T>&>(vx);
   SimpleVector<T>& y = dynamic_cast<SimpleVector<T>&>(vy);
   SimpleVector<T>& z = dynamic_cast<SimpleVector<T>&>(vz);

   std::copy(v, v + x.length(), x.v);
   std::copy(v + x.length(), v + x.length() + y.length(), y.v);
   std::copy(v + x.length() + y.length(), v + x.length() + y.length() + z.length(), z.v);
}

template<typename T>
void SimpleVector<T>::appendToBack(unsigned int n_to_add, const T& value) {
   assert(!preserveVec);

   const int new_len = this->n + n_to_add;

   T* new_v = new T[new_len];

   std::uninitialized_copy(this->v, this->v + this->n, new_v);
   std::uninitialized_fill(new_v + this->n, new_v + new_len, value);

   delete[] this->v;

   this->n = new_len;
   this->v = new_v;
}

template<typename T>
void SimpleVector<T>::appendToBack(const SimpleVector<T>& other) {
   assert(!preserveVec);

   const int new_len = this->n + other.n;

   T* new_v = new T[new_len];

   std::uninitialized_copy(this->v, this->v + this->n, new_v);
   std::uninitialized_copy(other.v, other.v + other.n, new_v + this->n);

   delete[] this->v;

   this->n = new_len;
   this->v = new_v;
}

template<typename T>
SimpleVector<T>* SimpleVector<T>::shaveBorder(int n_shave, bool shave_top) {
   // TODO : adjust for n_shave == this->n
   assert(n_shave <= this->n);
   assert(0 <= n_shave);
   auto* vec_new = new SimpleVector<T>(n_shave);

   T* vec_shaved = new T[this->n - n_shave];

   if (shave_top) {
      std::copy(v, v + n_shave, vec_new->v);
      std::copy(v + n_shave, v + this->n, vec_shaved);
   }
   else {
      std::copy(v + this->n - n_shave, v + this->n, vec_new->v);
      std::copy(v, v + this->n - n_shave, vec_shaved);
   }
   delete[] v;
   v = vec_shaved;
   this->n -= n_shave;

   return vec_new;
}


template<typename T>
void SimpleVector<T>::getSumCountIfSmall(double tol, double& sum_small, int& n_close, const Vector<T>* select) const {
   if (this->n == 0)
      return;

   const auto* selects = dynamic_cast<const SimpleVector<T>*>(select);
   if (selects)
      assert(this->n == selects->n);

   for (int i = 0; i < this->n; ++i) {
      const bool small = PIPSisZero(v[i], tol) && ((selects && PIPSisEQ((*selects)[i], 1.0)) || selects == nullptr);

      if (small) {
         sum_small += v[i];
         ++n_close;
      }
   }
}

template<typename T>
void SimpleVector<T>::pushSmallComplementarityPairs(Vector<T>& other_vec_in, const Vector<T>& select_in, double tol_this, double tol_other,
      double tol_pairs) {
   assert(tol_other > 0);
   assert(tol_this > 0);
   assert(tol_pairs > 0);
   assert(this->n == other_vec_in.length());
   assert(this->n == select_in.length());

   SimpleVector<T>& other_vec = dynamic_cast<SimpleVector<T>&>(other_vec_in);
   const SimpleVector<T>& select = dynamic_cast<const SimpleVector<T>&>(select_in);

   for (int i = 0; i < this->n; ++i) {
      if (PIPSisEQ(1.0, select[i])) {
         T& x = this->v[i];
         T& y = other_vec[i];

         const T complementary_product = x * y;
         assert(complementary_product >= 0.0);
         if (complementary_product < tol_pairs) {
            if (x < y)
               x = tol_pairs / y;
            else
               y = tol_pairs / x;
         }

         if (x < tol_this)
            x += tol_this;
         if (y < tol_other)
            y += tol_other;
      }
   }
}


template
class SimpleVector<int>;

template
class SimpleVector<double>;
