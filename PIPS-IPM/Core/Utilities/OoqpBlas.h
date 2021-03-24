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


#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif

extern "C"


{
   void FNAME(dsymv)( char * uplo, int * n, double * alpha,
            double * A, int * lda,
            const double * x, int * incx,
            double * beta, double * y, int * incy );

   void FNAME(dgemv)( char * trans, int * m, int * n,
            double * alpha, double * a, int * lda,
            const double * x, int * incx,
                 double * beta, double * y, int * incy );

   void FNAME(dgemm) ( char * transA, char * transB,
            int * m, int * n, int * k,
            double * alpha,
            double * A, int * ldA,
            double * B, int * ldB,
            double * beta, double * C, int * ldC);

   void FNAME(dpotrf)(char *uplo, int *n, double *A, int*lda, int*info);

   void FNAME(dtrsm)(char *side,
         char *uplo, char *transa, char *diag,
         int* m, int *n, double *alpha, double *A,
         int *lda, double *B, int *ldb);

   void FNAME(dsyrk)(char *uplo,
         char *trans, int* n, int* k, double *alpha,
         double *A, int *lda, double *beta,
         double *C, int *ldc);

   void FNAME(dsyr)( char * uplo, int * n,
         double * alpha, double * x,
         int * incx, double * a,
         int * lda );

   void FNAME(dtrsv)(char *uplo, char *trans, char* diag, int *n,
         double *A, int *lda, double *x, int *incx);

   void FNAME(dscal)(const int *n, const double *alpha, double *x, const int* incx);

   void FNAME(daxpy)( const int* n, const double* alpha, double x[], const int* incx,
            const double y[], const int* incy );

   double FNAME(ddot)( const int* n, const double dx[], const int* incx, const double dy[], const int* incy );


   int FNAME(idamax)( const int* n, const double x[], const int* incx);
}


#endif
