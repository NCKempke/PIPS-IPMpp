/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include <cmath>
#include <cstring>
#include <cassert>
#include "SparseStorage.h"
#include "Vector.hpp"
#include "SimpleVector.h"
#include "pipsdef.h"
#include <limits>
#include <fstream>
#include <string>
#include <algorithm>

int SparseStorage::instances = 0;

SparseStorage::SparseStorage(int m_, int n_, int len_) : m{m_}, n{n_}, len{len_} {
   assert(m_ >= 0);
   assert(len_ >= 0);
   if (n_ <= 0)
      assert(len_ == 0);
   neverDeleteElts = 0;

   jcolM = new int[len];
   krowM = new int[m + 1];
   M = new double[len];

   std::fill(krowM, krowM + m + 1, 0);
   std::fill(jcolM, jcolM + len, 0);
   std::fill(M, M + len, 0);

   SparseStorage::instances++;
}

SparseStorage::SparseStorage(int m_, int n_, int len_, int* krowM_, int* jcolM_, double* M_, int deleteElts) : neverDeleteElts{!deleteElts}, m{m_},
      n{n_}, len{len_}, jcolM{jcolM_}, krowM{krowM_}, M{M_} {
   assert(m_ >= 0);
   assert(n_ >= 0);
   assert(len_ >= 0);

   SparseStorage::instances++;
}

SparseStorage::~SparseStorage() {
   if (!neverDeleteElts) {
      delete[] jcolM;
      delete[] krowM;
      delete[] M;
   }

   SparseStorage::instances--;
}


void SparseStorage::copyFrom(int* krowM_, int* jcolM_, double* M_) const {
   memcpy(jcolM_, jcolM, len * sizeof(jcolM[0]));
   memcpy(M_, M, len * sizeof(M[0]));
   memcpy(krowM_, krowM, (m + 1) * sizeof(krowM[0]));
}

std::pair<int, int> SparseStorage::n_rows_columns() const {
   return {m, n};
}

int SparseStorage::n_rows() const {
   return m;
}

int SparseStorage::n_columns() const {
   return n;
}

void SparseStorage::fromGetDiagonal(int idiag, Vector<double>& vec_in) const {
   auto& vec = dynamic_cast<SimpleVector<double>&>(vec_in);
   int extent = vec.length();

   assert(idiag + extent <= m);
   assert(idiag + extent <= n);

   for (int i = idiag; i < idiag + extent; i++) {

      vec[i - idiag] = 0.0;
      for (int k = krowM[i]; k < krowM[i + 1]; k++) {

         const int j = jcolM[k];
         if (i == j) {
            vec[i - idiag] = M[k];
         }
      }
   }
}

void SparseStorage::columnScale(const Vector<double>& scale_in) {
   const auto& scale = dynamic_cast<const SimpleVector<double>&>(scale_in);
   assert(scale.length() == n);

   for (int i = 0; i < m; ++i) {
      for (int k = krowM[i]; k < krowM[i + 1]; ++k) {
         const int j = jcolM[k];
         M[k] = M[k] * scale[j];
      }
   }
}

void SparseStorage::rowScale(const Vector<double>& scale_in) {
   const auto& scale = dynamic_cast<const SimpleVector<double>&>(scale_in);
   assert(scale.length() == m);

   for (int i = 0; i < m; ++i) {
      for (int k = krowM[i]; k < krowM[i + 1]; ++k)
         M[k] = M[k] * scale[i];
   }
}

void SparseStorage::symmetricScale(const Vector<double>& scale_in) {
   const auto& scale = dynamic_cast<const SimpleVector<double>&>(scale_in);

   assert(scale.length() == n);
   assert(scale.length() == m);

   for (int i = 0; i < m; ++i) {
      for (int k = krowM[i]; k < krowM[i + 1]; ++k) {
         const int j = jcolM[k];
         M[k] = M[k] * scale[j] * scale[i];
      }
   }
}

void SparseStorage::scalarMult(double num) {
   for (int i = 0; i < m; ++i) {
      for (int k = krowM[i]; k < krowM[i + 1]; ++k)
         M[k] = M[k] * num;
   }
}

void SparseStorage::getDiagonal(Vector<double>& vec_in) const {
   this->fromGetDiagonal(0, vec_in);
}

void SparseStorage::setToDiagonal(const Vector<double>& vec_in) {
   const auto& vec = dynamic_cast<const SimpleVector<double>&>(vec_in);
   int diagExtent = std::min(m, n);

   assert(diagExtent == vec.length());

   int i;
   for (i = 0; i <= diagExtent; i++) {
      krowM[i] = i; // Initialize to the diagonal matrix of all zeros
   }
   for (; i <= m; i++) {
      krowM[i] = diagExtent;
   }

   const double* v = &vec[0];
   for (i = 0; i < diagExtent; i++) {
      jcolM[i] = i;
      M[i] = v[i];
   }
}

bool SparseStorage::isValid() const {
   assert(krowM && jcolM && M);

   if (m < 0 || n < 0 || len < 0) {
      printf("isValid: negative size parameter \n");
      return false;
   }

   if (krowM[0] != 0 || krowM[m] != len) {
      printf("isValid: krowM broken \n");
      return false;
   }

   for (int i = 0; i < len; i++)
      if (jcolM[i] < 0 || jcolM[i] >= n) {
         printf("isValid: column index out of bounds \n");
         return false;
      }

   for (int i = 0; i < m; i++)
      if (krowM[i] > krowM[i + 1]) {
         printf("isValid: row indices wrongly ordered \n");
         return false;
      }

   return true;
}

bool SparseStorage::isSorted() const {
   assert(isValid());

   for (int i = 0; i < m; i++) {
      for (int j = krowM[i] + 1; j < krowM[i + 1]; j++) {
         const int col = jcolM[j];
         const int prevcol = jcolM[j - 1];

         if (col <= prevcol)
            return false;
      }
   }

   return true;
}

void SparseStorage::fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const {

   assert(row >= 0 && row + rowExtent <= m);
   assert(col >= 0 && col + colExtent <= n);

   for (int i = row; i < row + rowExtent; i++) {

      int jcurrent = col - 1;
      for (int k = krowM[i]; k < krowM[i + 1]; k++) {

         const int j = jcolM[k];
         if (j >= col) {
            // j is big enough to be within range
            if (j < col + colExtent) {
               // j is small enough to be within range
               for (jcurrent++; jcurrent < j; jcurrent++) {
                  A[(i - row) * lda + jcurrent - col] = 0.0;
               }
               jcurrent = j;
               A[(i - row) * lda + j - col] = M[k];
            }
            else { // j is too big.
               // There will be no more interesting elements in this row
               break;
            }
         }
      }
      for (jcurrent++; jcurrent < col + colExtent; jcurrent++) {
         A[(i - row) * lda + jcurrent - col] = 0.0;
      }
   }
}

// used in backsolves
// get a dense block of columns *in column-major format*
// A must be zero'd on input
// allzero is true if there are actually no nonzero entries in this block
// colSparsity contains on exit the sparsity pattern of union of columns
void SparseStorage::fromGetColBlock(int col, double* A, int lda, int colExtent, int* colSparsity, bool& allzero) const {
   for (int i = 0; i < m; i++) {
      for (int k = krowM[i]; k < krowM[i + 1]; k++) {
         const int j = jcolM[k];
         if (j >= col) {
            if (j < col + colExtent) {
               A[i + (j - col) * lda] = M[k];
               allzero = false;
               colSparsity[i] = 1;
            }
            else {
               break;
            }
         }
      }
   }
}

void SparseStorage::fromGetRowsBlock(double* rows_array_dense, size_t row_start, size_t n_rows, size_t array_line_size, size_t array_line_offest,
      int* row_sparsity) const {
   assert(array_line_offest + n <= array_line_size);
   assert(rows_array_dense);
   assert(static_cast<int>(row_start + n_rows) <= m);

   for (size_t i = 0; i < n_rows; ++i) {
      const size_t row = row_start + i;

      if (krowM[row] == krowM[row + 1])
         continue;

      const size_t offset = i * array_line_size + array_line_offest;
      assert(offset >= 0);
      for (int k = krowM[row]; k < krowM[row + 1]; ++k) {
         const int col = jcolM[k];
         const double val = M[k];

         assert(offset + col < array_line_size * n_rows);
         if (row_sparsity)
            row_sparsity[array_line_offest + col] = 1;

         rows_array_dense[offset + col] = val;

      }
   }
}


void SparseStorage::fromGetRowsBlock(const int* rowIndices, int nRows, int arrayLineSize, int arrayLineOffset, double* rowsArrayDense,
      int* rowSparsity) const {
   assert(rowsArrayDense && rowIndices);
   assert(arrayLineSize >= 0 && arrayLineOffset >= 0);

   // todo use OMP?
   for (int i = 0; i < nRows; i++) {
      const int r = rowIndices[i];
      assert(r >= 0 && r < m);

      if (krowM[r] == krowM[r + 1])
         continue;

      const int offset = i * arrayLineSize + arrayLineOffset;

      for (int c = krowM[r]; c < krowM[r + 1]; c++) {
         const int col = jcolM[c];
         assert(offset >= 0);

         if (rowSparsity)
            rowSparsity[arrayLineOffset + col] = 1;

         rowsArrayDense[offset + col] = M[c];
      }
   }
}

void SparseStorage::getLinkVarsNnz(std::vector<int>& vec) const {
   assert(int(vec.size()) == n);

   for (int i = 0; i < len; i++) {
      const int col = jcolM[i];
      assert(col < n);

      vec[col]++;
   }
}

void SparseStorage::atPutSpRow(int row, const double A[], int lenA, const int jcolA[], int& info) {
   int ik;
   int ka = lenA - 1;
   int km_f = krowM[row + 1] - 1;
   int km = km_f;
   int km_s = krowM[row];
   int count = 0;

   assert(row >= 0 && row < m);

   while (ka >= 0) {
      if (km < km_s) {
         // There are no more elements in M. All the rest of A must be
         // inserted.
         count += ka + 1;
         break;
      }
      else if (jcolM[km] == jcolA[ka]) {
         // The element in A will replace an element in M
         km--;
         ka--;
      }
      else if (jcolM[km] > jcolA[ka]) {
         assert(jcolA[ka] >= 0);
         // This element is in M but not in A
         km--;
      }
      else {
         // The element is in A, but not in M and so must be inserted.
         assert(jcolA[ka] < n);
         ka--;
         count++;
      }
   }

   if (count > 0) {
      this->shiftRows(row + 1, count, info);
   }
   else {
      info = 0;
   }

   if (0 == info) {
      ka = lenA - 1;
      km = km_f;
      ik = krowM[row + 1] - 1;

      while (ka >= 0) {
         if (km < km_s) {
            // There are no more elements in M. All the rest of A must be
            // inserted.
            for (; ka >= 0; ka--, ik--) {
               jcolM[ik] = jcolA[ka];
               assert(jcolM[ik] >= 0 && jcolM[ik] < n);
               M[ik] = A[ka];
            }
            break;
         }
         else if (jcolM[km] == jcolA[ka]) {
            // The element in A will replace an element in M
            jcolM[ik] = jcolM[km];
            M[ik] = A[ka];
            km--;
            ka--;
            ik--;
         }
         else if (jcolM[km] > jcolA[ka]) {
            // This element is in M but not in A
            jcolM[ik] = jcolM[km];
            M[ik] = M[km];
            km--;
            ik--;
         }
         else {
            // The element is in A, but not in M.
            jcolM[ik] = jcolA[ka];
            assert(jcolM[ik] >= 0 && jcolM[ik] < n);
            M[ik] = A[ka];

            ka--;
            ik--;
         }
      }
   }
}

void SparseStorage::atPutDense(int row, int col, const double* A, int lda, int rowExtent, int colExtent) {
   int info, count;
   int i, km_f, km, ka, k;

   assert(row >= 0 && row + rowExtent <= m);
   assert(col >= 0 && col + colExtent <= n);

   for (i = row; i < row + rowExtent; i++) {
      // Loop over all rows in range.
      km_f = krowM[i + 1] - 1;
      km = km_f;
      ka = colExtent - 1;
      count = 0;
      while (km >= krowM[i] || ka >= 0) {
         // The current row in M and the current row in A are not
         // both empty
         if (ka < 0) {
            // The current row of A is empty. Insert an element from M
            km--;
         }
         else if (km < krowM[i]) {
            // The current row of M is empty. Insert an element from A.
            if (A[(i - row) * lda + ka] == 0) {
               // The current element of A is zero. Skip it.
               ka--;
            }
            else {
               count++;
               ka--;
            }
         }
         else if (ka + col > jcolM[km]) {
            // The current element in A comes first.
            if (A[(i - row) * lda + ka] == 0) {
               // The current element of A is zero. Skip it.
               ka--;
            }
            else {
               count++;
               ka--;
            }
         }
         else if (ka + col < jcolM[km]) {
            // The current element in M comes first.
            km--;
         }
         else {
            // The current element in M is overwritten by the element in A.
            km--;
            ka--;
         }
      } // End while the current row...are not both empty
      if (count > 0) {
         this->shiftRows(i + 1, count, info);
      }
      else {
         info = 0;
      }
      if (info != 0) {
         std::cout << "bing\n";
         return;
      }

      km = km_f;
      k = krowM[i + 1] - 1;
      ka = colExtent - 1;

      while (km >= krowM[i] || ka >= 0) {
         // The current row in M and the current row in A are not
         // both empty
         if (ka < 0) {
            // The current row of A is empty. Insert an elemnt from M
            M[k] = M[km];
            jcolM[k] = jcolM[km];
            k--;
            km--;
         }
         else if (km < krowM[i]) {
            // The current row of M is empty. Insert an element from A.
            if (A[(i - row) * lda + ka] == 0) {
               // The current element of A is zero. Skip it.
               ka--;
            }
            else {
               M[k] = A[(i - row) * lda + ka];
               jcolM[k] = ka + col;
               k--;
               ka--;
            }
         }
         else if (ka + col > jcolM[km]) {
            // The current element in A comes first.
            if (A[(i - row) * lda + ka] == 0) {
               // The current element of A is zero. Skip it.
               ka--;
            }
            else {
               M[k] = A[(i - row) * lda + ka];
               jcolM[k] = ka + col;
               k--;
               ka--;
            }
         }
         else if (ka + col < jcolM[km]) {
            // The current element in M comes first.
            M[k] = M[km];
            jcolM[k] = jcolM[km];
            k--;
            km--;
         }
         else {
            // The current element in M is overwritten by the element in A.
            M[k] = A[(i - row) * lda + ka];
            jcolM[k] = ka + col;
            k--;
            km--;
            ka--;
         }
      } // End while the current row...are not both empty
   } // End loop over all rows in range.
}

void SparseStorage::shiftRows(int row, int shift, int& info) {
   if (shift == 0) {
      info = 0;
   }
   else if (krowM[m] + shift > len) {
      // Insufficient space
      info = krowM[m] + shift - len;
   }
   else {
      // We perform the copy
      info = 0;
      int lcopy = krowM[m] - krowM[row];
      if (lcopy > 0) {
         // There is anything to copy
         // As a consequence of lcopy > 0, col !== n
         memmove(&jcolM[krowM[row] + shift], &jcolM[krowM[row]], lcopy * sizeof(int));
         memmove(&M[krowM[row] + shift], &M[krowM[row]], lcopy * sizeof(double));
         int i;
         for (i = row; i <= m; i++) {
            krowM[i] += shift;
         }
      }
      else {
         // Still adjust the starts of the rows
         int i;
         int rowStart = krowM[m] + shift;
         for (i = row; i <= m; i++) {
            krowM[i] = rowStart;
         }
      }
   } // end else we perform the copy
}

void SparseStorage::putSparseTriple(const int irow[], int lenA, const int jcol[], const double A[], int& info) {
   if (len < lenA) {
      info = 1;
      return;
   }
   info = 0;
   int i, k;
   krowM[0] = 0;
   i = 0;
   for (k = 0; k < lenA; k++) {
      while (i < irow[k]) {
         i++;
         krowM[i] = k;
      }
      // now i == irow[k] because irow is sorted
      jcolM[k] = jcol[k];
      M[k] = A[k];
   }
   for (i++; i <= m; i++) {
      krowM[i] = lenA;
   }
}

void SparseStorage::fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info) const {
   assert(col >= 0 && col < n);
   assert(row >= 0 && row < m);
   assert(col + colExtent <= n);
   int ka = 0;
   int lastCol = col + colExtent - 1;

   info = 0;

   for (int km = krowM[row]; km < krowM[row + 1]; km++) {
      const int colm = jcolM[km];
      if (colm >= col) {
         if (colm <= lastCol) {
            if (ka < lenA) {
               A[ka] = M[km];
               jcolA[ka] = colm;
               ka++;
            }
            else {
               // Count the number of aditional elements needed in A
               info++;
            }
         }
         else {
            break;
         }
      }
   }
   nnz = ka;
}

void SparseStorage::writeToStream(std::ostream& out) const {
   int i, k;

   for (i = 0; i < m; i++) {
      for (k = krowM[i]; k < krowM[i + 1]; k++) {
         out << i << '\t' << jcolM[k] << '\t' << M[k] << "\n";
      }
   }
}

void SparseStorage::writeNNZpatternToStreamDense(std::ostream& out) const {
   for (int row = 0; row < m; row++) {
      int col = 0;
      for (int k = krowM[row]; k < krowM[row + 1]; k++) {
#ifndef NDEBUG
         // assert that columns are ordered
         assert(k == krowM[row] || jcolM[k - 1] < jcolM[k]);
#endif

         while (jcolM[k] > col) {
            out << "   ";
            col++;
         }
         out << " * ";
         col++;
      }
      while (col < n) {
         out << "   ";
         col++;
      }
      out << "\n";
   }
}

void SparseStorage::writeToStreamDense(std::ostream& out) const {
   //todo: instead of \t, use length of longest value in M

   for (int i = 0; i < m; i++) {
      int j = 0;
      for (int k = krowM[i]; k < krowM[i + 1]; k++) {
#ifndef NDEBUG
         // assert that columns are ordered
         assert(k == krowM[i] || jcolM[k - 1] < jcolM[k]);
#endif
         while (jcolM[k] > j) {
            out << 0 << '\t';
            j++;
         }
         out << M[k] << '\t';
         j++;
      }
      while (j < n) {
         out << 0 << '\t';
         j++;
      }
      out << "\n";
   }
}

void SparseStorage::writeToStreamDenseRow(std::ostream& out, int rowidx) const {
   int j = 0; // Column j
   for (int k = krowM[rowidx]; k < krowM[rowidx + 1]; k++) {
      while (jcolM[k] > j) {
         out << 0 << '\t';
         j++;
      }
      out << M[k] << '\t';
      j++;
   }
   while (j < n) {
      out << 0 << '\t';
      j++;
   }
}

void indexedLexSort(int first[], int n, int swapFirst, int second[], int swapSecond, int index[]) {
   int fi, se, j, k, kinc, inc, ktemp;
   const int nincs = 12;
   const int incs[] = {1, 5, 19, 41, 109, 209, 505, 929, 2161, 3905, 8929, 16001};

   kinc = 0;
   for (k = 0; k < nincs; k++) {
      kinc = k;
      if (incs[kinc] > n / 2) {
         kinc--;
         break;
      }
   }
   // incs[kinc] is the greatest value in the sequence that is also less
   // than n/2.

   //for( k = 0; k < n; k++ ) index[k] = k;

   for (; kinc >= 0; kinc--) {
      // Loop over all increments
      inc = incs[kinc];

      if (!swapFirst && !swapSecond) {
         for (k = inc; k < n; k++) {
            // loop over all subarrays defined by the current increment
            ktemp = index[k];
            fi = first[ktemp];
            se = second[ktemp];
            // Insert element k into the sorted subarray
            for (j = k; j >= inc; j -= inc) {
               // Loop over the elements in the current subarray
               if (fi < first[index[j - inc]] || (fi == first[index[j - inc]] && se < second[index[j - inc]])) {
                  // Swap elements j and j - inc, implicitly use the fact
                  // that ktemp hold element j to avoid having to assign to
                  // element j - inc
                  index[j] = index[j - inc];
               }
               else {
                  // There are no more elements in this sorted subarray which
                  // are less than element j
                  break;
               }
            } // End loop over the elements in the current subarray
            // Move index[j] out of temporary storage
            index[j] = ktemp;
            // The element has been inserted into the subarray.
         } // End loop over all subarrays defined by the current increment
      }
      else if (swapSecond && !swapFirst) {
         for (k = inc; k < n; k++) {
            ktemp = index[k];
            fi = first[ktemp];
            se = second[k];
            for (j = k; j >= inc; j -= inc) {
               if (fi < first[index[j - inc]] || (fi == first[index[j - inc]] && se < second[j - inc])) {
                  index[j] = index[j - inc];
                  second[j] = second[j - inc];
               }
               else {
                  break;
               }
            }
            index[j] = ktemp;
            second[j] = se;
         }
      }
      else if (swapFirst && !swapSecond) {
         for (k = inc; k < n; k++) {
            ktemp = index[k];
            fi = first[k];
            se = second[ktemp];
            for (j = k; j >= inc; j -= inc) {
               if (fi < first[j - inc] || (fi == first[j - inc] && se < second[index[j - inc]])) {
                  index[j] = index[j - inc];
                  first[j] = first[j - inc];
               }
               else {
                  break;
               }
            }
            index[j] = ktemp;
            first[j] = fi;
         }
      }
      else { // Swap both
         for (k = inc; k < n; k++) {
            ktemp = index[k];
            fi = first[k];
            se = second[k];
            for (j = k; j >= inc; j -= inc) {
               if (fi < first[j - inc] || (fi == first[j - inc] && se < second[j - inc])) {
                  index[j] = index[j - inc];
                  first[j] = first[j - inc];
                  second[j] = second[j - inc];
               }
               else {
                  break;
               }
            }
            index[j] = ktemp;
            first[j] = fi;
            second[j] = se;
         }
      }
   } // End loop over all increments
}


void SparseStorage::mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   int i, j, k;
   double temp;
   for (i = 0; i < m; i++) {
      temp = 0;
      for (k = krowM[i]; k < krowM[i + 1]; k++) {
         j = jcolM[k];
         temp += M[k] * x[j * incx];
#ifndef NDEBUG
         assert(j < n);
#endif
      }
      y[i * incy] = beta * y[i * incy] + alpha * temp;
   }
}

void SparseStorage::multSym(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   for (int i = 0; i < m; i++)
      y[i * incy] *= beta;

   for (int i = 0; i < n; i++) {
      for (int k = krowM[i]; k < krowM[i + 1]; k++) {
         const int j = jcolM[k];

         y[i * incy] += alpha * M[k] * x[j * incx];
         // todo fixme won't work with Q or any other "really" symmetric matrix!
         // Necessary because CtDC from sLinsysRootAug is stored as a general matrix
         if (i != j)//&& 0 )
            y[j * incy] += alpha * M[k] * x[i * incx];
      }
   }
}

void SparseStorage::transMult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   if (beta != 1.0)
      for (int j = 0; j < n; j++)
         y[j * incy] *= beta;

   for (int i = 0; i < m; i++) {
      for (int k = krowM[i]; k < krowM[i + 1]; k++) {
         const int j = jcolM[k];

         assert(j < n);

         y[j * incy] += alpha * M[k] * x[i * incx];
      }
   }
}

void SparseStorage::transMultD(double beta, double y[], int incy, double alpha, const double x[], const double d[], int incxd) const {
   if (beta != 1.0)
      for (int j = 0; j < n; j++)
         y[j * incy] *= beta;

   for (int i = 0; i < m; i++) {
      for (int k = krowM[i]; k < krowM[i + 1]; k++) {
         const int j = jcolM[k];

         assert(j < n);
         assert(!PIPSisZero(d[i * incxd]));

         y[j * incy] += alpha * M[k] * x[i * incxd] / d[i * incxd];
      }
   }
}

// adds this * x to the part of row yrow of y that is in the upper half of y
// y is assumed to be symmetric and sorted!
void SparseStorage::multMatSymUpper(double beta, SparseStorage& y, double alpha, const double x[], int yrow, int ycolstart) const {
   assert(yrow >= 0 && yrow < y.m);
   assert(ycolstart >= 0 && ycolstart < y.n);
   assert(y.n == y.m);
   assert(y.n >= m + ycolstart);

   int* const krowM_y = y.krowM;
   int* const jcolM_y = y.jcolM;
   double* const M_y = y.M;

   // assert that yrow is sorted
#ifndef NDEBUG
   for (int ci = krowM_y[yrow] + 1; ci < krowM_y[yrow + 1]; ci++)
      assert(jcolM_y[ci - 1] < jcolM_y[ci]);
#endif

   // scale row yrow
   if (beta != 1.0)
      for (int c_y = krowM_y[yrow]; c_y != krowM_y[yrow + 1]; c_y++)
         M_y[c_y] *= beta;

   // add this * x to yrow (and exploit that y is symmetric)
   int c_y = krowM_y[yrow];
   for (int r = 0; r < m; r++) {
      double yrx;
      const int colplace_y = r + ycolstart; // the column in y where to place yrx

      // not in upper half of y?
      if (colplace_y < yrow)
         continue;

      // compute y_(r,.) * x
      yrx = 0.0;

      for (int c = krowM[r]; c != krowM[r + 1]; c++) {
         const int col = jcolM[c];
         yrx += x[col] * M[c];
      }

      if (PIPSisZero(yrx))
         continue;

      yrx *= alpha;

      for (; c_y != krowM_y[yrow + 1]; c_y++) {
         const int col_y = jcolM_y[c_y];
         if (col_y == colplace_y) {
            M_y[c_y] += yrx;
            break;
         }
      }
      assert(c_y != krowM_y[yrow + 1]);
   }
}

// TODO : ugly
void doubleLexSort(int first[], int n, int second[], double data[]);

void SparseStorage::symmetrize(int& info) {
   int i, k, ku;

   int nnz = krowM[m];
   int* irowM = new int[2 * nnz];

   info = 0;

   ku = nnz;
   for (i = 0; i < m; i++) {
      // For all rows
      for (k = krowM[i]; k < krowM[i + 1]; k++) {
         // Loop over elements of row i
         irowM[k] = i;
         if (i != jcolM[k]) {
            // Not a diagonal element
            if (ku >= len) {
               info++;
            }
            else {
               // Add the element transpose to the scrambled matrix
               irowM[ku] = jcolM[k];
               jcolM[ku] = i;
               M[ku] = M[k];
               ku++;
            }
         } // End not a diagonal element
      } // End loop over elements of column j
   } // End for all columns
   if (info != 0) {
      delete[] irowM;
      return;
   }

   nnz = ku;

   doubleLexSort(irowM, nnz, jcolM, M);
   i = 0;
   krowM[0] = 0;
   for (k = 0; k < nnz; k++) {
      for (; i < irowM[k]; i++) {
         krowM[i + 1] = k;
      }
   }
   for (i++; i <= m; i++) {
      krowM[i] = nnz;
   }
   delete[] irowM;
}

double SparseStorage::inf_norm() const {
   double norm = 0.0;
   const int nnz = numberOfNonZeros();

   for (int i = 0; i < nnz; i++) {
      const double absM = std::fabs(M[i]);
      if (absM > norm)
         norm = absM;
   }
   return norm;
}

double SparseStorage::abminnormNonZero(double tol) const {
   double norm = std::numeric_limits<double>::infinity();
   int nnz = this->numberOfNonZeros();

   for (int i = 0; i < nnz; i++) {
      const double fabsMi = std::fabs(M[i]);
      if (fabsMi < norm && tol < fabsMi)
         norm = fabsMi;
   }
   return norm;
}

void SparseStorage::atPutDiagonal(int idiag, const Vector<double>& vvec) {
   const auto& v = dynamic_cast<const SimpleVector<double>&>(vvec);
   atPutDiagonal(idiag, &v[0], 1, v.length());
}

void SparseStorage::atAddDiagonal(int idiag, const Vector<double>& vvec) {
   const auto& v = dynamic_cast<const SimpleVector<double>&>(vvec);
   atAddDiagonal(idiag, &v[0], 1, v.length());
}

void SparseStorage::atPutDiagonal(int idiag, const double x[], int incx, int extent) {
   int i;
   int info;
   for (i = idiag; i < idiag + extent; i++) {
      // Search for the diagonal elt.
      int lastk;
      lastk = krowM[i + 1];
      for (int k = krowM[i]; k < lastk; k++) {
         if (i >= jcolM[k]) {
            if (i == jcolM[k])
               M[k] = x[incx * (i - idiag)];
            else
               atPutSpRow(i, &x[incx * (i - idiag)], 1, &i, info);
            break;
         }
      }
   }
}

void SparseStorage::atAddDiagonal(int idiag, const double x[], int incx, int extent) {
   assert(idiag + extent <= m);
   assert(idiag + extent <= n);

   for (int row = idiag; row < idiag + extent; row++) {
      int k = krowM[row];
      while (jcolM[k] != row && k < krowM[row + 1])
         ++k;
      if (row == jcolM[k])
         this->M[k] += x[incx * (row - idiag)];
      else
         assert(false);
   }
}

void SparseStorage::diagonal_set_to_constant_from(int from, int length, double value) {
   assert(from + length <= m);
   assert(from + length <= n);

   for (int row = from; row < from + length; ++row) {
      int col_ptr = krowM[row];
      while (jcolM[col_ptr] != row && col_ptr < krowM[row + 1])
         ++col_ptr;
      if (row == jcolM[col_ptr])
         M[col_ptr] = value;
      else
         assert(false);
   }
}

void SparseStorage::diagonal_add_constant_from(int from, int length, double value) {
   assert(from + length <= m);
   assert(from + length <= n);

   for (int row = from; row < from + length; ++row) {
      int col_ptr = krowM[row];
      while (jcolM[col_ptr] != row && col_ptr < krowM[row + 1])
         ++col_ptr;
      if (row == jcolM[col_ptr])
         M[col_ptr] += value;
      else
         assert(false);
   }
}

void SparseStorage::transpose(int* krowMt, int* jcolMt, double* Mt) const {
   int ind, pend;//,ptend;
   /////////////////////////////////////////////////////
   // form the transpose
   /////////////////////////////////////////////////////
   const int nnz = krowM[m];

   //cumulative sum to find the number of elements in each row of At, ie column of A.
   int* w = new int[n];
   for (int i = 0; i < n; i++)
      w[i] = 0;
   for (int i = 0; i < nnz; i++)
      w[jcolM[i]]++;

   krowMt[0] = 0; //double sum=0.0;
   for (int i = 1; i <= n; i++) {
      krowMt[i] = krowMt[i - 1] + w[i - 1];
      w[i - 1] = krowMt[i - 1];
   }

   //populate jcolMt and Mt now
   for (int i = 0; i < m; i++) {
      pend = krowM[i + 1];
      for (int j = krowM[i]; j < pend; j++) {
         ind = w[jcolM[j]];
         jcolMt[ind] = i;
         Mt[ind] = M[j];
         w[jcolM[j]] = ind + 1; //w[jcolM[j]]++
      }
   }
   //we have the transpose
   delete[] w;
}

void SparseStorage::clear() {
   for (int i = 0; i < len; i++)
      M[i] = 0.0;
}

void SparseStorage::matTransDSymbMultMat(const double*, const int* krowMt, const int* jcolMt, const double*, int** krowAtDA, int** jcolAtDA, double** AtDA) const {
   int k, j, nnzAtA = 0, pend, ptend;

   ////////////////////////////////////////////////////
   // count the number of entries in the result AtA  //
   ////////////////////////////////////////////////////
   char* flag = new char[n];
   for (int i = 0; i < n; i++) {
      memset(flag, 0, n);

      ptend = krowMt[i + 1];
      for (int pt = krowMt[i]; pt < ptend; pt++) {
         k = jcolMt[pt]; //At[i,j] is nz

         //add the nonzero pattern of row i of A to AtA
         pend = krowM[k + 1];
         for (int p = krowM[k]; p < pend; p++) {
            j = jcolM[p];

            if (flag[j] == 0) {
               nnzAtA++;
               flag[j] = 1;
            }
         }
      }
   }
   assert(nnzAtA >= 0); //overflow?!?

   delete[] flag;
   ////////////////////////////////////////////////////
   // alocate AtDA
   ///////////////////////////////////////////////////
   *krowAtDA = new int[n + 1];
   *jcolAtDA = new int[nnzAtA];
   *AtDA = new double[nnzAtA];

   std::uninitialized_fill(*krowAtDA, *krowAtDA + n + 1, 0);
   std::uninitialized_fill(*jcolAtDA, *jcolAtDA + nnzAtA, 0);
   std::uninitialized_fill(*AtDA, *AtDA + nnzAtA, 0);
}


void SparseStorage::matTransDMultMat(const double* d, const int* krowMt, const int* jcolMt, const double* Mt, int* krowAtDA, int* jcolAtDA, double* AtDA) const {
   int k, j, pend, ptend;
   double val;
   ////////////////////////////////////////////////////
   // AtDA = A'*D*A
   ////////////////////////////////////////////////////
   auto* W = new double[n];
   auto* flag = new char[n];

   for (int it = 0; it < n; it++)
      W[it] = 0.0;

   int nnzAtA = 0;
   for (int i = 0; i < n; i++) {
      memset(flag, 0, n);

      //start row i of AtDA
      krowAtDA[i] = nnzAtA;

      ptend = krowMt[i + 1];
      for (int pt = krowMt[i]; pt < ptend; pt++) {
         k = jcolMt[pt]; //At[i,k] is non-zero
         val = Mt[pt] * d[k];

         //iterate the row k of A and scatter the values into W
         pend = krowM[k + 1];
         for (int p = krowM[k]; p < pend; p++) {
            j = jcolM[p];
            //we have A[k,j]
            if (flag[j] == 0) {
               jcolAtDA[nnzAtA++] = j;
               flag[j] = 1;
            }

            W[j] += (M[p] * val);
         }
      }
      //gather the values into the i-th row AtDA
      for (int p = krowAtDA[i]; p < nnzAtA; p++) {
         j = jcolAtDA[p];
         AtDA[p] = W[j];
         W[j] = 0.0;
      }
   }
   krowAtDA[n] = nnzAtA;
   delete[] W;
   delete[] flag;
}

void SparseStorage::matTransDinvMultMat(const double* d, const int* krowMt, const int* jcolMt, const double* Mt, int* krowAtDA, int* jcolAtDA, double* AtDA) const {
   int k, j, pend, ptend;
   double val;
   ////////////////////////////////////////////////////
   // AtDA = A'*D*A
   ////////////////////////////////////////////////////
   auto* W = new double[n];
   auto* flag = new char[n];

   for (int it = 0; it < n; it++)
      W[it] = 0.0;

   int nnzAtA = 0;
   for (int i = 0; i < n; i++) {
      memset(flag, 0, n);

      //start row i of AtDA
      krowAtDA[i] = nnzAtA;

      ptend = krowMt[i + 1];
      for (int pt = krowMt[i]; pt < ptend; pt++) {
         k = jcolMt[pt]; //At[i,k] is non-zero
         val = Mt[pt] / d[k];

         //iterate the row k of A and scatter the values into W
         pend = krowM[k + 1];
         for (int p = krowM[k]; p < pend; p++) {
            j = jcolM[p];
            //we have A[k,j]
            if (flag[j] == 0) {
               jcolAtDA[nnzAtA++] = j;
               flag[j] = 1;
            }

            W[j] += (M[p] * val);
         }
      }
      //gather the values into the i-th row AtDA
      for (int p = krowAtDA[i]; p < nnzAtA; p++) {
         j = jcolAtDA[p];
         AtDA[p] = W[j];
         W[j] = 0.0;
      }
   }
   krowAtDA[n] = nnzAtA;
   delete[] W;
   delete[] flag;
}

void SparseStorage::reduceToLower() {
   assert(m == n); //available only for square matrices
   int newNnz = 0; // the new number of nonzeros
   //get the nnz of the reduced matrix
   for (int i = 0; i < n; i++) {
      int jend = krowM[i + 1];
      for (int j = krowM[i]; j < jend; j++)
         if (jcolM[j] <= i)
            newNnz++;
   }

   auto* newjcolM = new int[newNnz];
   auto* newdM = new double[newNnz];
   int it = 0;

   for (int i = 0; i < n; i++) {
      int jend = krowM[i + 1];
      int jstart = krowM[i];
      krowM[i] = it;
      for (int j = jstart; j < jend; j++) {
         if (jcolM[j] <= i) {
            newjcolM[it] = jcolM[j];
            newdM[it] = M[j];
            it++;
         }
      }
   }
   assert(it == newNnz);
   krowM[n] = it;

   delete[] jcolM;
   delete[] M;
   jcolM = newjcolM;
   M = newdM;
   //for(int k=krowC
}

void SparseStorage::dump(const std::string& filename) const {
   std::ofstream fd(filename.c_str());
   fd << std::scientific;
   fd.precision(16);
   fd << m << " " << n << " " << len << "\n";

   int i;
   for (i = 0; i <= m; i++) {
      fd << krowM[i] << " ";
   }
   fd << "\n";
   for (i = 0; i < len; i++) {
      fd << jcolM[i] << " ";
   }
   fd << "\n";
   for (i = 0; i < len; i++) {
      fd << M[i] << " ";
   }
}

void SparseStorage::deleteEmptyRowsCols(const int* nnzRowVec, const int* nnzColVec) {
   assert(nnzRowVec != nullptr && nnzColVec != nullptr);

   int m_new = 0;
   int n_new = 0;

   int* rowsmap = new int[m];

   for (int i = 0; i < m; i++)
      if (nnzRowVec[i] != 0.0)
         rowsmap[i] = m_new++;

   int* colsmap = new int[n];

   for (int i = 0; i < n; i++)
      if (nnzColVec[i] != 0.0)
         colsmap[i] = n_new++;

   assert(m_new >= 0 && m_new <= m);
   assert(n_new >= 0 && n_new <= n);

   assert(krowM[m] == 0 && "method not fully supported");

#if 0
   int nnz = 0;
   int rowcount = 0;
   int len_new = 0;

   for( int r = 0; r < m; r++ )
   {
      const int nnzRow = nnzRowVec[r];

      if( nnzRow == 0.0 )
      {
         nnz += krowM[r + 1] - krowM[r];
         continue;
      }
      krowM[rowcount] = len_new;

      for( int j = krowM[r]; j < krowM[r + 1]; j++ )
      {
         jcolM[len_new] = jcolM[colsmap[nnz]];
         M[len_new++] = M[nnz++];
      }
      rowcount++;
   }


   for( int r = rowcount; r <= m; r++ )
      krowM[r] = len_new;

   assert(nnz == len);
#endif

   m = m_new;
   n = n_new;

   delete[] colsmap;
   delete[] rowsmap;
}


void SparseStorage::getSparseTriplet_c2fortran(int*& irn, int*& jcn, double*& val) const {
   int count = 0;
   assert(len > 0);
   assert(!irn && !jcn && !val);
   assert(!fortranIndexed());

   irn = new int[len];
   jcn = new int[len];
   val = new double[len];

   for (int r = 0; r < m; r++) {
      for (int c = krowM[r]; c < krowM[r + 1]; c++) {
         const int col = jcolM[c];
         const double value = M[c];

         irn[count] = r + 1;
         jcn[count] = col + 1;
         val[count] = value;

         count++;
      }
   }

   assert(count == len);
}


void SparseStorage::getSparseTriplet_fortran2fortran(int*& irn, int*& jcn, double*& val) const {
   int count = 0;
   assert(len > 0);
   assert(!irn && !jcn && !val);
   assert(fortranIndexed());

   irn = new int[len];
   jcn = new int[len];
   val = new double[len];

   for (int row_f = 1; row_f <= m; row_f++) {
      const int row_c = row_f - 1;

      const int row_start_c = krowM[row_c] - 1;
      const int row_end_c = krowM[row_c + 1] - 1;

      for (int col_index = row_start_c; col_index < row_end_c; col_index++) {
         const int col_f = jcolM[col_index];
         const double value = M[col_index];

         irn[count] = row_f;
         jcn[count] = col_f;
         val[count] = value;

         count++;
      }
   }

   assert(count == len);
}

void SparseStorage::deleteEmptyRows(int*& orgIndex) {
   assert(!neverDeleteElts);
   assert(orgIndex == nullptr);

   int m_new = 0;

   // count non-empty rows
   for (int r = 0; r < m; r++)
      if (krowM[r] != krowM[r + 1])
         m_new++;

   int* krowM_new = new int[m_new + 1];
   orgIndex = new int[m_new + 1];

   krowM_new[0] = 0;
   m_new = 0;

   for (int r = 0; r < m; r++)
      if (krowM[r] != krowM[r + 1]) {
         orgIndex[m_new] = r;
         krowM_new[++m_new] = krowM[r + 1];
      }

   assert(krowM_new[m_new] == len);
   m = m_new;

   delete[] krowM;
   krowM = krowM_new;
}


void SparseStorage::c2fortran() {
   assert(krowM[0] == 0 && krowM[m] == len && !isFortranIndexed);

   for (int i = 0; i <= m; i++)
      krowM[i]++;

   for (int i = 0; i < len; i++)
      jcolM[i]++;

   isFortranIndexed = true;
}

void SparseStorage::fortran2c() {
   assert(krowM[0] == 1 && krowM[m] == len + 1 && isFortranIndexed);

   for (int i = 0; i <= m; i++)
      krowM[i]--;

   for (int i = 0; i < len; i++)
      jcolM[i]--;

   isFortranIndexed = false;
}

bool SparseStorage::fortranIndexed() const {
   return isFortranIndexed;
}

void SparseStorage::set2FortranIndexed() {
   assert(krowM[0] == 1 && krowM[m] == len + 1);

   isFortranIndexed = true;
}

void SparseStorage::deleteZeroRowsColsSym(int*& new2orgIdx) {
   assert(m == n);
   assert(!neverDeleteElts);
   assert(new2orgIdx == nullptr);
   assert(this->isValid());

   int* const offset = new int[m];

   for (int r = 0; r < m; r++)
      offset[r] = 0;

   // mark rows (and columns) to be deleted
   for (int r = 0; r < m; r++) {
      const int start = krowM[r];
      const int end = krowM[r + 1];

      if (start == end) {
         offset[r] = -1;
         continue;
      }

      int c = start;

      for (; c < end; c++)
         if (!PIPSisZero(M[c]))
            break;

      // no non-zero found?
      if (c == end) {
         offset[r] = -1;
         continue;
      }

      offset[r] = -2;

      for (c = start; c < end; c++) {
         const int col = jcolM[c];

         assert(offset[col] != 0);

         if (!PIPSisZero(M[c]) && offset[col] == -1)
            offset[col] = -2;
      }
   }

   int rowDeletes = 0;
   int zeroEntryDeletes = 0;

   // count column offsets and entries to be deleted
   for (int r = 0; r < m; r++) {
      const int start = krowM[r];
      const int end = krowM[r + 1];

      // row deleted?
      if (offset[r] == -1) {
         zeroEntryDeletes += end - start;
         rowDeletes++;
         continue;
      }

      for (int c = start; c < end; c++) {
         const int col = jcolM[c];
         assert(col < m);

         if (offset[col] == -1) {
            assert(PIPSisZero(M[c]));
            zeroEntryDeletes++;
         }
      }

      assert(offset[r] == -2);
      offset[r] = rowDeletes;
   }

   const int m_new = m - rowDeletes;
   const int len_new = len - zeroEntryDeletes;
   assert(len_new >= 0 && m_new >= 0);

   new2orgIdx = new int[m_new];
   auto* const krowM_new = new int[m_new + 1];
   auto* const jcolM_new = new int[len_new];
   auto* const M_new = new double[len_new];
   int m_count = 0;
   int len_count = 0;

   // fill the new arrays
   krowM_new[0] = 0;
   for (int r = 0; r < m; r++) {
      if (offset[r] == -1)
         continue;

      for (int c = krowM[r]; c < krowM[r + 1]; c++) {
         const int col = jcolM[c];
         if (offset[col] == -1) {
            assert(PIPSisZero(M[c]));
            continue;
         }

         assert(col - offset[col] >= 0);

         jcolM_new[len_count] = col - offset[col];

         M_new[len_count++] = M[c];
      }

      new2orgIdx[m_count] = r;
      krowM_new[++m_count] = len_count;
      assert(krowM_new[m_count] > krowM_new[m_count - 1]);
   }

   assert(m_count == m_new);
   assert(len_count == len_new);

   delete[] krowM;
   delete[] jcolM;
   delete[] M;
   delete[] offset;

   m = m_new;
   n = m_new;
   len = len_new;
   krowM = krowM_new;
   jcolM = jcolM_new;
   M = M_new;

   assert(this->isValid());
}

void SparseStorage::addNnzPerRow(int* vec, int begin_rows, int end_rows) const {
   assert(0 <= begin_rows && begin_rows <= end_rows && end_rows <= m);

   for (int r = begin_rows; r < end_rows; r++)
      vec[r - begin_rows] += krowM[r + 1] - krowM[r];
}

void SparseStorage::addRowSums(double* vec) const {
   for (int r = 0; r < m; r++) {
      const int end = krowM[r + 1];
      for (int c = krowM[r]; c < end; c++)
         vec[r] += std::fabs(M[c]);
   }
}

void SparseStorage::getRowMinVec(const double* colScaleVec, double* vec) const {
   if (n <= 0)
      return;

   const bool coscale = (colScaleVec != nullptr);

   if (coscale) {
      for (int r = 0; r < m; r++) {
         double minval = vec[r];
         assert(minval >= 0.0);

         for (int i = krowM[r]; i < krowM[r + 1]; i++) {
            const double absval = std::abs(M[i] * colScaleVec[jcolM[i]]);

            if (absval < minval && absval > pips_eps)
               minval = absval;
         }
         vec[r] = minval;
      }
   }
   else {
      for (int r = 0; r < m; r++) {
         double minval = vec[r];
         assert(minval >= 0.0);

         for (int i = krowM[r]; i < krowM[r + 1]; i++) {
            const double absval = std::abs(M[i]);

            if (absval < minval && absval > pips_eps)
               minval = absval;
         }
         vec[r] = minval;
      }
   }
}

void SparseStorage::getRowMaxVec(const double* colScaleVec, double* vec) const {
   if (n <= 0)
      return;

   const bool coscale = (colScaleVec != nullptr);

   if (coscale) {
      for (int r = 0; r < m; r++) {
         double maxval = vec[r];
         assert(maxval >= 0.0);

         for (int i = krowM[r]; i < krowM[r + 1]; i++) {
            const double absval = std::abs(M[i] * colScaleVec[jcolM[i]]);

            if (absval > maxval)
               maxval = absval;
         }
         vec[r] = maxval;
      }
   }
   else {
      for (int r = 0; r < m; r++) {
         double maxval = vec[r];
         assert(maxval >= 0.0);

         for (int i = krowM[r]; i < krowM[r + 1]; i++) {
            const double absval = std::abs(M[i]);

            if (absval > maxval)
               maxval = absval;
         }
         vec[r] = maxval;
      }
   }
}

void SparseStorage::getRowMinMaxVec(bool getMin, const double* colScaleVec, double* vec) const {
   if (getMin)
      getRowMinVec(colScaleVec, vec);
   else
      getRowMaxVec(colScaleVec, vec);
}


void SparseStorage::permuteRows(const std::vector<unsigned int>& permvec) {
   assert(permvec.size() == size_t(m));

   if (len == 0)
      return;

   assert(m > 0 && n > 0);

   auto* jcolM_new = new int[len];
   auto* krowM_new = new int[m + 1];
   auto* M_new = new double[len];

   int len_new = 0;
   krowM_new[0] = 0;

   for (int r = 0; r < m; ++r) {
      const unsigned int r_perm = permvec[r];
      const int rowlength = krowM[r_perm + 1] - krowM[r_perm];

      assert(r_perm < static_cast<unsigned int>(m) && rowlength >= 0);

      if (rowlength > 0) {
         memcpy(jcolM_new + len_new, jcolM + krowM[r_perm], rowlength * sizeof(int));
         memcpy(M_new + len_new, M + krowM[r_perm], rowlength * sizeof(double));

         len_new += rowlength;
      }

      krowM_new[r + 1] = len_new;
   }

   assert(len_new == len);

   delete[] jcolM;
   delete[] krowM;
   delete[] M;

   jcolM = jcolM_new;
   krowM = krowM_new;
   M = M_new;
}


void SparseStorage::permuteCols(const std::vector<unsigned int>& permvec) {
   assert(int(permvec.size()) == n);

   std::vector<int> permvec_rev(n, 0);

   auto* indexvec = new int[n];
   auto* bufferCol = new int[n];
   auto* bufferM = new double[n];

   for (int i = 0; i < n; i++) {
      assert(permvec[i] < unsigned(n));

      permvec_rev[permvec[i]] = i;
   }

   for (int r = 0; r < m; ++r) {
      const int row_start = krowM[r];
      const int row_end = krowM[r + 1];
      const int row_length = row_end - row_start;

      if (row_length == 0)
         continue;

      for (int c = row_start; c < row_end; c++) {
         const int col = jcolM[c];
         assert(col < n);

         jcolM[c] = permvec_rev[col];
      }

      for (int i = 0; i < row_length; i++)
         indexvec[i] = i;

      std::sort(indexvec, indexvec + row_length, index_sort(jcolM + row_start, row_length));

      for (int i = 0; i < row_length; i++) {
         assert(indexvec[i] < row_length);

         bufferCol[i] = jcolM[row_start + indexvec[i]];
         bufferM[i] = M[row_start + indexvec[i]];
      }

      memcpy(jcolM + row_start, bufferCol, row_length * sizeof(int));
      memcpy(M + row_start, bufferM, row_length * sizeof(double));
   }

   delete[] indexvec;
   delete[] bufferM;
   delete[] bufferCol;
}


void SparseStorage::sortCols() {
   auto* indexvec = new int[n];
   auto* bufferCol = new int[n];
   auto* bufferM = new double[n];

   for (int r = 0; r < m; ++r) {
      const int row_start = krowM[r];
      const int row_end = krowM[r + 1];
      const int row_length = row_end - row_start;

      if (row_length == 0)
         continue;

      for (int i = 0; i < row_length; i++)
         indexvec[i] = i;

      std::sort(indexvec, indexvec + row_length, index_sort(jcolM + row_start, row_length));

      for (int i = 0; i < row_length; i++) {
         assert(indexvec[i] < row_length);

         bufferCol[i] = jcolM[row_start + indexvec[i]];
         bufferM[i] = M[row_start + indexvec[i]];
      }

      memcpy(jcolM + row_start, bufferCol, row_length * sizeof(int));
      memcpy(M + row_start, bufferM, row_length * sizeof(double));
   }

   delete[] indexvec;
   delete[] bufferM;
   delete[] bufferCol;
}


/*
 * computes the full sparse matrix representation from a upper triangular symmetric sparse representation
 *
 * Must be square, the storage for the full representation will be allocated within the matrix and must be released later
 */
void SparseStorage::fullMatrixFromUpperTriangular(int*& rowPtrFull, int*& colIdxFull, double*& valuesFull) const {
   assert(n == m);

   /* cout elems per row and assert upper triangular */
   std::vector<int> nelems(n, 0);
   for (int i = 0; i < n; ++i) {
      // diag elem
      for (int j = krowM[i]; j < krowM[i + 1]; ++j) {
         assert(jcolM[j] >= i);

         if (i == jcolM[j])
            nelems[i]++;
         else {
            nelems[jcolM[j]]++;
            nelems[i]++;
         }
      }
   }

   // fill rowptr array
   rowPtrFull = new int[n + 1];

   rowPtrFull[0] = 0;
   for (int i = 0; i < n; ++i)
      rowPtrFull[i + 1] = rowPtrFull[i] + nelems[i];

   colIdxFull = new int[rowPtrFull[n]];
   for (int i = 0; i < rowPtrFull[n]; ++i)
      colIdxFull[i] = -1;

   valuesFull = new double[rowPtrFull[n]];

   // fill in col and value
   for (int i = 0; i < n; ++i) {
      int rowstart = krowM[i];
      int rowend = krowM[i + 1];

      for (int k = rowstart; k < rowend; ++k) {
         double value = M[k];
         int col = jcolM[k];

         int kk = rowPtrFull[i];
         int colfull = colIdxFull[kk];
         while (colfull != -1) {
            kk++;
            colfull = colIdxFull[kk];
         }

         colIdxFull[kk] = col;
         valuesFull[kk] = value;

         if (col != i) {
            kk = rowPtrFull[col];
            colfull = colIdxFull[kk];

            while (colfull != -1) {
               kk++;
               colfull = colIdxFull[kk];
            }

            colIdxFull[kk] = i;
            valuesFull[kk] = value;
         }
      }
   }
}

SparseStorage* SparseStorage::shaveLeft(int n_cols) {
   assert(n_cols <= n);
   assert(0 <= n_cols);

   if (n_cols == 0)
      return new SparseStorage(m, 0, 0);

   // TODO : adjust for when n_cols == n
   const int n_border = n_cols;
   const int m_border = m;

   int count = 0;
   for (int i = 0; i < len; ++i) {
      if (jcolM[i] < n_cols)
         ++count;
   }
   const int len_border = count;
   assert(0 <= len_border && len_border <= len);

   auto* krowM_border = new int[m_border + 1];
   auto* jcolM_border = new int[len_border];
   auto* M_border = new double[len_border];

   count = 0;
   krowM_border[0] = 0;
   for (int row = 0; row < m; ++row) {
      for (int j = krowM[row]; j < krowM[row + 1]; ++j) {
         const int col = jcolM[j];
         if (col < n_cols) {
            jcolM_border[count] = col;
            M_border[count] = M[j];
            ++count;
         }
      }
      krowM_border[row + 1] = count;
   }
   assert(count == len_border);

   const int len_new = len - len_border;
   const int n_new = n - n_cols;
   const int m_new = m;
   auto* krowM_new = new int[m_new + 1];
   auto* jcolM_new = new int[len_new];
   auto* M_new = new double[len_new];
   count = 0;
   krowM_new[0] = 0;
   for (int row = 0; row < m; ++row) {
      for (int j = krowM[row]; j < krowM[row + 1]; ++j) {
         const int col = jcolM[j];
         if (col >= n_cols) {
            jcolM_new[count] = col - n_cols;
            M_new[count] = M[j];
            ++count;
         }
      }
      krowM_new[row + 1] = count;
   }
   assert(count == len_new);

   len = len_new;
   m = m_new;
   n = n_new;
   std::swap(M, M_new);
   std::swap(krowM, krowM_new);
   std::swap(jcolM, jcolM_new);

   delete[] M_new;
   delete[] krowM_new;
   delete[] jcolM_new;

   auto* border = new SparseStorage(m_border, n_border, len_border, krowM_border, jcolM_border, M_border, true);

   assert(this->isValid());
   assert(border->isValid());

   return border;
}

SparseStorage* SparseStorage::shaveBottom(int n_rows) {
   assert(n_rows <= m);
   assert(0 <= n_rows);

   if (n_rows == 0)
      return new SparseStorage(0, n, 0);

   // TODO : adjust for n_rows == m
   const int n_border = n;
   const int m_border = n_rows;
   const int len_border = krowM[m] - krowM[m - n_rows];
   assert(0 <= len_border && len_border <= len);

   auto* krowM_border = new int[m_border + 1];
   auto* jcolM_border = new int[len_border];
   auto* M_border = new double[len_border];

   const int len_new = len - len_border;
   const int m_new = m - m_border;

   std::copy(jcolM + len_new, jcolM + len, jcolM_border);
   std::copy(M + len_new, M + len, M_border);

   for (int i = 0; i <= m_border; ++i) {
      assert(i + m_new <= m);
      assert(krowM[i + m_new] - len_new >= 0);
      krowM_border[i] = krowM[i + m_new] - len_new;
   }

   /* for this mat keep old storages but just set other sizes */
   m = m_new;
   len = len_new;
   assert(krowM[m] == len);

   auto* border = new SparseStorage(m_border, n_border, len_border, krowM_border, jcolM_border, M_border, true);
   assert(this->isValid());
   assert(border->isValid());

   return border;
}

void SparseStorage::dropNEmptyRowsBottom(int n_rows) {
   assert(n_rows <= m);
   assert(0 <= n_rows);

   if (n_rows == 0)
      return;

   // assert rows are empty
   assert(krowM[m] - krowM[m - n_rows] == 0);
   m -= n_rows;
   assert(this->isValid());
}

void SparseStorage::dropNEmptyRowsTop(int n_rows) {
   assert(n_rows <= m);
   assert(0 <= n_rows);
   assert(isValid());

   if (n_rows == 0)
      return;

   assert(krowM[0] - krowM[n_rows] == 0);

   std::move(krowM + n_rows, krowM + m + 1, krowM);
   m -= n_rows;

   assert(isValid());
}

SparseStorage* SparseStorage::shaveSymLeftBottom(int) {
   assert(0 && "TODO : implement");
   // TODO : assert is symmetric - either upper or lower?? then shave of elements from left and top at the same time
   return nullptr;

}

// concatenate matrices
// if "diagonal", make a block diagonal matrix:
// [ A 0 ]
// [ 0 B ]
// if not, stack the matrices:
// [ A ]
// [ B ]
// diagonal will be symmetric if the input is
/*SparseStorage::SparseStorage(const vector<SparseStorage*> &blocks, bool diagonal)
{
  assert(blocks.size() > 0);
  m = n = len = 0;
  neverDeleteElts = 0;
  for (size_t i = 0; i < blocks.size(); i++) {
    m += blocks[i]->m;
    len += blocks[i]->len;
    if (diagonal) {
      n += blocks[i]->n;
    } else {
      assert( blocks[i]->n == blocks[0]->n );
    }
  }
  if (!diagonal) {
    n = blocks[0]->n;
  }


  M     = new double[len];
  jcolM = new int[len];
  krowM = new int[m+1];

  int curnnz = 0, curn = 0, curm = 0;

  for (size_t i = 0; i < blocks.size(); i++) {
    int lnnz = blocks[i]->len;
    int ln = blocks[i]->n;
    int lm = blocks[i]->m;

    memcpy(M+curnnz, blocks[i]->M, lnnz*sizeof(double));
    memcpy(jcolM+curnnz, blocks[i]->jcolM, lnnz*sizeof(int));
    memcpy(krowM+curm, blocks[i]->krowM, lm*sizeof(int)); // skips last element

    if (diagonal) {
      for (int j = 0; j < lnnz; j++) {
        jcolM[curnnz+j] += curn;
      }
    }
    for (int j = 0; j < lm; j++) {
      krowM[curm+j] += curnnz;
    }
    curnnz += lnnz;
    curn += ln;
    curm += lm;
  }
  assert(curm == m);
  krowM[m] = len;

  SparseStorage::instances++;
}
*/