/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP 
 *
 * dgemm_ definition added by C. Petra 10/09
 */
#ifndef OOQPBLAS_H
#define OOQPBLAS_H

/**
 * Blas routines used in OOQP 
 */

extern "C"


{
void dsymv (char* uplo, int* n, double* alpha, double* A, int* lda, const double* x, int* incx, double* beta, double* y, int* incy);

void dgemv (char* trans, int* m, int* n, double* alpha, double* a, int* lda, const double* x, int* incx, double* beta, double* y, int* incy);

void dgemm (char* transA, char* transB, int* m, int* n, int* k, double* alpha, double* A, int* ldA, double* B, int* ldB, double* beta, double* C,
      int* ldC);

void dpotrf (char* uplo, int* n, double* A, int* lda, int* info);

void dtrsm (char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double* A, int* lda, double* B, int* ldb);

void dsyrk (char* uplo, char* trans, int* n, int* k, double* alpha, double* A, int* lda, double* beta, double* C, int* ldc);

void dsyr (char* uplo, int* n, double* alpha, double* x, int* incx, double* a, int* lda);

void dtrsv (char* uplo, char* trans, char* diag, int* n, double* A, int* lda, double* x, int* incx);

void dscal (const int* n, const double* alpha, double* x, const int* incx);

void daxpy (const int* n, const double* alpha, double x[], const int* incx, const double y[], const int* incy);

double ddot (const int* n, const double dx[], const int* incx, const double dy[], const int* incy);


int idamax (const int* n, const double x[], const int* incx);

/* solving dense linear systems */
void dsytrf(const char* uplo, const int* n, double A[], const int* lda, int ipiv[], double work[], const int* lwork, int* info);

void dsyrfs(const char* uplo, const int* n, const int* n_rhs, const double A[], const int* lda, const double afact[], const int* ldaf, const int ipiv[],
   const double rhs[], const int* lrhs, double x[], const int* ldx, double forward_error[], double backward_error[], double work[], int iwork[], int* info);

void dsytrs(const char* uplo, const int* n, const int* nrhs, double A[], const int* lda, int ipiv[], double b[], const int* ldb, int* info);

}


#endif
