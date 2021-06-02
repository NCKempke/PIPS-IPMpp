//
// Created by nils-christian on 02.06.21.
//

#include "BorderedMatrix.h"

#ifndef PIPSIPMPP_BORDEREDMATRIXLIFTEDA0WRAPPER_H
#define PIPSIPMPP_BORDEREDMATRIXLIFTEDA0WRAPPER_H

class BorderedMatrixLiftedA0wrapper : public BorderedMatrix {
public:
   BorderedMatrixLiftedA0wrapper(std::unique_ptr<BorderedMatrix> other);

   void addRowSums(Vector<double>& vec_) const override;
   void getRowMinMaxVec(bool get_min, bool initialize_vec, const Vector<double>* col_scale_in, Vector<double>& minmax_in) const override;
   void getColMinMaxVec(bool get_min, bool initialize_vec, const Vector<double>* row_scale_in, Vector<double>& minmax_in) const override;
   void rowScale(const Vector<double>& vec) override;
   void transMult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const override;
   void mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const override;
private:
   template<typename T>
   void set_column_vector_first_for_child(const Vector<T>& vec) const;

   template<typename T>
   void reset_column_vector_first_for_child(const Vector<T>& vec) const;
};

#endif //PIPSIPMPP_BORDEREDMATRIXLIFTEDA0WRAPPER_H
