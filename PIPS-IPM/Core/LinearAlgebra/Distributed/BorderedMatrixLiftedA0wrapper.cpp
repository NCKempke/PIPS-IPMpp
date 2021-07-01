//
// Created by nils-christian on 02.06.21.
//

#include "BorderedMatrixLiftedA0wrapper.h"
#include "DistributedVector.h"
#include "StripMatrix.h"

BorderedMatrixLiftedA0wrapper::BorderedMatrixLiftedA0wrapper(std::unique_ptr<BorderedMatrix> other) : BorderedMatrix(std::move(other->inner_matrix),
   std::move(other->border_left), std::move(other->border_bottom), std::move(other->bottom_left_block), other->get_mpi_comm()) {};

void BorderedMatrixLiftedA0wrapper::mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   set_column_vector_first_for_child(y_in);
   BorderedMatrix::mult(beta, y_in, alpha, x_in);
   reset_column_vector_first_for_child(y_in);
}

void BorderedMatrixLiftedA0wrapper::transpose_mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   set_column_vector_first_for_child(x_in);
   BorderedMatrix::transpose_mult(beta, y_in, alpha, x_in);
   reset_column_vector_first_for_child(x_in);
}

void BorderedMatrixLiftedA0wrapper::rowScale(const Vector<double>& vec) {
   set_column_vector_first_for_child(vec);
   BorderedMatrix::rowScale(vec);
   reset_column_vector_first_for_child(vec);
}

void BorderedMatrixLiftedA0wrapper::getRowMinMaxVec(bool get_min, bool initialize_vec, const Vector<double>* col_scale_in, Vector<double>& minmax_in) const {
   set_column_vector_first_for_child(minmax_in);
   BorderedMatrix::getRowMinMaxVec(get_min, initialize_vec, col_scale_in, minmax_in);
   reset_column_vector_first_for_child(minmax_in);
}

void BorderedMatrixLiftedA0wrapper::getColMinMaxVec(bool get_min, bool initialize_vec, const Vector<double>* row_scale_in, Vector<double>& minmax_in) const {
   if (row_scale_in)
      set_column_vector_first_for_child(*row_scale_in);
   BorderedMatrix::getColMinMaxVec(get_min, initialize_vec, row_scale_in, minmax_in);
   if (row_scale_in)
      reset_column_vector_first_for_child(*row_scale_in);
}

void BorderedMatrixLiftedA0wrapper::addRowSums(Vector<double>& vec) const {
   set_column_vector_first_for_child(vec);
   BorderedMatrix::addRowSums(vec);
   reset_column_vector_first_for_child(vec);
}

template<typename T>
void BorderedMatrixLiftedA0wrapper::set_column_vector_first_for_child(const Vector<T>& vec_in) const {
   const auto& vec = dynamic_cast<const DistributedVector<double>&>(vec_in);

   assert(vec.first && !vec.children[0]->first);
   vec.children[0]->first = vec.first;
}

template<typename T>
void BorderedMatrixLiftedA0wrapper::reset_column_vector_first_for_child(const Vector<T>& vec_in) const {
   const auto& vec = dynamic_cast<const DistributedVector<double>&>(vec_in);

   assert(vec.children[0]->first);
   vec.children[0]->first = nullptr;
}