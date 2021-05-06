/*
 * BorderedSymMatrix.C
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */

#include "BorderedSymmetricMatrix.h"
#include "DoubleMatrixTypes.h"
#include "SimpleVector.h"
#include "pipsdef.h"
#include <algorithm>

BorderedSymmetricMatrix::BorderedSymmetricMatrix(DistributedSymmetricMatrix* inner_matrix_, StripMatrix* border_vertical_, SymmetricMatrix* top_left_block_, MPI_Comm mpiComm_)
      : inner_matrix(inner_matrix_), border_vertical(border_vertical_), top_left_block(top_left_block_), mpiComm(mpiComm_),
      iAmDistrib(mpiComm == MPI_COMM_NULL) {
   assert(inner_matrix);
   assert(border_vertical);
   assert(top_left_block);

   assert(inner_matrix->children.size() == border_vertical->children.size());
   assert(border_vertical->is_vertical);

   inner_matrix->getSize(n, n);

   int n_border, m_border;
   border_vertical->getSize(m_border, n_border);

   n += n_border;

#ifndef NDEBUG
   int n_bottom;
   top_left_block->getSize(n_bottom, n_bottom);
   int n_inner;
   inner_matrix->getSize(n_inner, n_inner);

   assert(n_inner == m_border);
   assert(n_bottom == n_border);
#endif
}

BorderedSymmetricMatrix::~BorderedSymmetricMatrix() {
   delete inner_matrix;
   delete border_vertical;
   delete top_left_block;
}

int BorderedSymmetricMatrix::is_a(int type) const {
   return (type == kBorderedSymMatrix || type == kSymMatrix || type == kBorderedMatrix);
}

/** y = beta * y + alpha * this * x */
void BorderedSymmetricMatrix::mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   auto& y = dynamic_cast<DistributedVector<double>&>(y_in);
   const auto& x = dynamic_cast<const DistributedVector<double>&>(x_in);

   assert(x.children.size() == 1 && y.children.size() == 1);
   assert(x.children[0] && y.children[0]);
   assert(x.first);
   assert(y.first);
   assert(!x.last);
   assert(!y.last);

   top_left_block->mult(beta, *y.first, alpha, *x.first);
   border_vertical->transMult(1.0, *y.first, alpha, *x.children[0]);

   border_vertical->mult(beta, *y.children[0], alpha, *x.first);
   inner_matrix->mult(1.0, *y.children[0], alpha, *x.children[0]);
}

/** y = beta * y + alpha * this^T * x */
void BorderedSymmetricMatrix::transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   this->mult(beta, y, alpha, x);
}

void BorderedSymmetricMatrix::fromGetDiagonal(int idiag, Vector<double>& x_in) const {
   assert("The value of the parameter idiag is not supported!" && idiag == 0);

   auto& x = dynamic_cast<DistributedVector<double>&>(x_in);
   assert(x.children.size() == 1);
   assert(x.children[0]);
   assert(!x.last && x.first);

   top_left_block->getDiagonal(*x.first);

   inner_matrix->fromGetDiagonal(idiag, *x.children[0]);
}

double BorderedSymmetricMatrix::inf_norm() const {
   double norm = -std::numeric_limits<double>::infinity();

   norm = std::max(norm, inner_matrix->inf_norm());
   norm = std::max(norm, border_vertical->inf_norm());
   norm = std::max(norm, top_left_block->inf_norm());

   return norm;
}

void BorderedSymmetricMatrix::scalarMult(double num) {
   inner_matrix->scalarMult(num);
   border_vertical->scalarMult(num);
}

void BorderedSymmetricMatrix::getSize(long long& m_, long long& n_) const {
   m_ = n;
   n_ = n;
}

void BorderedSymmetricMatrix::getSize(int& m_, int& n_) const {
   m_ = static_cast<int>(n);
   n_ = static_cast<int>(n);
}

long long BorderedSymmetricMatrix::size() const {
   return n;
}
