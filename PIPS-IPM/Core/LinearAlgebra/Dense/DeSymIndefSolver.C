/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "DeSymIndefSolver.h"
#include "SimpleVector.h"
#include <cassert>
#include <memory>

#include "DenseSymmetricMatrix.h"
#include "DenseMatrix.h"


#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif


#if defined(TIMING)
#include "mpi.h"
#endif

// declarations for LAPACK functions used to factor/solve:

// TODO : move to ooqpblas ?
// dsytrf_() factors a symmetric indefinite matrix A, see LAPACK
// documentation for more details.
extern "C" void FNAME(dsytrf)(const char* uplo, const int* n, double A[], const int* lda, int ipiv[], double work[], const int* lwork, int* info);

// dsytrs_() solves the system Ax = b using the factor obtained by dsytrf_().
extern "C" void
FNAME(dsytrs)(const char* uplo, const int* n, const int* nrhs, double A[], const int* lda, int ipiv[], double b[], const int* ldb, int* info);

#ifdef TIMING_FLOPS
extern "C" {
    void HPM_Init(void);
    void HPM_Start(char *);
    void HPM_Stop(char *);
    void HPM_Print(void);
    void HPM_Print_Flops(void);
    void HPM_Print_Flops_Agg(void);
    void HPM_Terminate(char*);
}
#endif

DeSymIndefSolver::DeSymIndefSolver(const DenseSymmetricMatrix* dm) : matrix{*dm}, is_mat_sparse{false}, n{static_cast<int>(dm->size())} {
   mStorage = std::make_unique<DenseStorage>(n,n);
   ipiv.resize(this->n);
}

DeSymIndefSolver::DeSymIndefSolver(const SparseSymmetricMatrix* sm) : matrix{*sm}, is_mat_sparse{true}, n{static_cast<int>(sm->size())} {
   mStorage = std::make_unique<DenseStorage>(n,n);
   ipiv.resize(this->n);
}

void DeSymIndefSolver::matrixChanged() {
   int info;

   if (n == 0)
      return;

   if (is_mat_sparse) {
      this->mStorage->fill_from_sparse(dynamic_cast<const SparseSymmetricMatrix&>(matrix).getStorageRef());
   } else {
      this->mStorage->fill_from_dense(dynamic_cast<const DenseSymmetricMatrix&>(matrix).getStorageRef());
   }

#ifdef TIMING
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if( myrank == 0 )
      std::cout << "DENSE_SYM_INDEF: starting factorization" << std::endl;
#endif

   const int query_workspace_size = -1;
   double optimal_workspace_size_output;
   FNAME(dsytrf)(&fortranUplo, &n, &mStorage->M[0][0], &n, ipiv.data(), &optimal_workspace_size_output, &query_workspace_size, &info);

   assert(0 < optimal_workspace_size_output && optimal_workspace_size_output < std::numeric_limits<int>::max());
   const int optimal_workspace_size = static_cast<int>( optimal_workspace_size_output );

   if (optimal_workspace_size > static_cast<int>(work.size()))
      work.resize(optimal_workspace_size);

#ifdef TIMING_FLOPS
      HPM_Start("DSYTRFFact");
#endif
   //factorize
   FNAME(dsytrf)(&fortranUplo, &n, &mStorage->M[0][0], &n, ipiv.data(), work.data(), &optimal_workspace_size, &info);

#ifdef TIMING_FLOPS
   HPM_Stop("DSYTRFFact");
#endif
   if (info != 0)
      printf("DeSymIndefSolver::matrixChanged : error - dsytrf returned info=%d\n", info);
}

void DeSymIndefSolver::solve(Vector<double>& v) {
   int info;
   const int one = 1;

   auto& sv = dynamic_cast<SimpleVector<double>&>(v);

   if (n == 0)
      return;

#ifdef TIMING_FLOPS
      HPM_Start("DSYTRSSolve");
#endif

   FNAME(dsytrs)(&fortranUplo, &n, &one, &mStorage->M[0][0], &n, ipiv.data(), &sv[0], &n, &info);

#ifdef TIMING_FLOPS
   HPM_Stop("DSYTRSSolve");
#endif
   assert(info == 0);
}

void DeSymIndefSolver::solve(GeneralMatrix& rhs_in) {
   auto& rhs = dynamic_cast<DenseMatrix&>(rhs_in);

   int info;
   const int ncols = rhs.n_rows();

   FNAME(dsytrs)(&fortranUplo, &n, &ncols, &mStorage->M[0][0], &n, ipiv.data(), &rhs[0][0], &n, &info);

   assert(info == 0);
}

void DeSymIndefSolver::diagonalChanged(int /* idiag */, int /* extent */) {
   this->matrixChanged();
}

void DeSymIndefSolver::calculate_inertia_from_factorization() const {
   positive_eigenvalues = 0;
   negative_eigenvalues = 0;
   // TODO: what whould be zero evs?
   zero_eigenvalues = 0;

   for(int i = 0; i < this->n; ++i) {
      /* 1x1 diagonal block */
      if(ipiv[i] > 0) {
         if(mStorage->M[i][i] > 0) {
            ++positive_eigenvalues;
         } else if (mStorage->M[i][i] < 0) {
            ++negative_eigenvalues;
         }
         /* 2x2 diagonal block */
      } else {
         assert(i + 1 < mStorage->n);
         assert(ipiv[i + 1] < 0);

         ++i;
         ++positive_eigenvalues;
         ++negative_eigenvalues;
      }
   }
}

std::tuple<unsigned int, unsigned int, unsigned int> DeSymIndefSolver::get_inertia() const {

   calculate_inertia_from_factorization();

   assert(mStorage->n - positive_eigenvalues - negative_eigenvalues >= 0);
   return {positive_eigenvalues, negative_eigenvalues, mStorage->n - positive_eigenvalues - negative_eigenvalues};
}