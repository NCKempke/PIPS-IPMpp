/*
 * pipschecks.h
 *
 *  Created on: 23.02.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_UTILITIES_PIPSCHECKS_H_
#define PIPS_IPM_CORE_UTILITIES_PIPSCHECKS_H_

#include "../LinearAlgebra/Dense/DenseVector.hpp"
#include "pipsdef.h"

// is the permutation vector valid?
bool permutation_is_valid(const Permutation& perm);

// are the columns of the given sub-matrix ordered?
bool submatrix_is_ordered(const int* rowptr, const int* colidx, int rowstart, int rowend);

// compute residual norms for Ax=rhs with A in CSR with 1-indexing (Fortran)
void compute_fortran_CSR_matrix_residual_norms(const int* rowptr, const int* colidx, const double* vals, const DenseVector<double>& rhs,
      const DenseVector<double>& x, double& res_norm2, double& res_nrmInf, double& sol_inf, double& mat_max);

#endif /* PIPS_IPM_CORE_UTILITIES_PIPSCHECKS_H_ */
