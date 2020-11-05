/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef MA27LINSYS_H
#define MA27LINSYS_H

#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "SparseStorageHandle.h"
#include "OoqpVectorHandle.h"
#include "pipsport.h"


#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif

extern "C"
{
  void FNAME(ma27id)( int icntl[], double cntl[] );

  void FNAME(ma27ad)( int * n,        int * nz,
      int irn[],      int icn[],
      int iw[],       int * liw,
      int ikeep[],    int iw1[],
      int * nsteps,   int * iflag,
      int icntl[],    double cntl[],
      int info[],     double * ops );

  void FNAME(ma27bd)( int * n,        int * nz,
      int irn[],      int icn[],
      double a[],     int * la,
      int iw[],       int * liw,
      int ikeep[],    int * nsteps,
      int * maxfrt,   int iw1[],
      int icntl[],    double cntl[],
      int  info[] );

  void FNAME(ma27cd)( int * n,
      double a[],     int * la,
      int iw[],       int * liw,
      double w[],     int * maxfrt,
      double rhs[],
      int iw2[],      int * nsteps,
      int icntl[],    int info[] );
}

/** implements the linear solver class using the HSL MA27 solver
 *
 * @ingroup LinearSolvers
 */
class Ma27Solver : public DoubleLinearSolver {

protected:
  int     icntl[30];
  int     info[20];
  double  cntl[5];

  int minimumRealWorkspace() { return info[4]; }
  int minimumIntWorkspace() { return info[5]; }

  const int max_tries = 8;

  int max_n_iter_refinement;

  const int ooqp_print_level_warnings = 10000;

  /** precision we demand from the linear system solver. If it isn't
   * attained on the first solve, we use iterative refinement and
   * possibly refactorization with a higher value of
   * kThresholdPivoting. */
  double precision;

  /** the Threshold Pivoting parameter may need to be increased during
   * the algorithm if poor precision is obtained from the linear
   * solves.  kThresholdPivoting indicates the largest value we are
   * willing to tolerate.  */
  double threshold_pivoting_max;

  /** During factorization entries smaller than small pivot will not
   * be accepted as pivots and the matrix will be treated as singular.
   *
   * This is the max we are willing to go with our pivots.
   */
  const double threshold_pivtol = 1e-14;
  const double threshold_pivtol_factor = 0.1;

  /** the factor in the range (1,inf) by which kThresholdPivoting is
   * increased when it is found to be inadequate.  */
  const double threshold_pivoting_factor = 10;

  /** index array for the factorization */
  int *irowM, *jcolM;

  /** nonzero element of the factors */
  double *fact;

  /** dimension of the whole matrix */
  int n;

  /** number of nonzeros in the matrix */
  int nnz;

  /** length of the array containing factors; may be increased during
   * the numerical factorization if the estimate obtained during the
   * symbolic factorization proves to be inadequate. */
  int la;

  /** pivot sequence and temporary storage information */
  int *ikeep, *iw, liw, *iw1, *iw2, nsteps, maxfrt;

  /** temporary storage for the factorization */
  double *w;

  /** amounts by which to increase allocated factorization space when
   * inadequate space is detected. ipessimism is for array "iw",
   * rpessimism is for the array "fact". */
  double ipessimism = 3.0;
  double rpessimism = 3.0;

  /** called the very first time a matrix is factored. Allocates space
   * for the factorization and performs ordering */
  virtual void firstCall();

public:
  /** the Threshold Pivoting parameter, stored as U in the ma27dd
   *  common block. Takes values in the range [-0.5,0.5]. If negative no
   *  numerical pivoting will be performed - if positive numerical pivoting
   *  will be performed. If no pivoting is applied the subroutine will fail when
   *  a zero pivot is encountered. Larger values enforce greater stability in the
   *  factorization as they insist on
   *  larger pivots. Smaller values preserve sparsity at the cost of
   *  using smaller pivots.
   */
  double thresholdPivoting() const { return cntl[0]; }
  void setThresholdPivoting( double piv  ) { cntl[0] = piv; }

  /** the "small pivot" parameter, stored as pivtol in the common
   * block ma27td. The factorization will not accept a pivot whose
   * absolute value is less than this parameter as a 1x1 pivot or as
   * the off-diagonal in a 2x2 pivot.  */
  double getSmallPivot() const { return cntl[2]; }
  void setSmallPivot( double tol ) { cntl[2] = tol; }

  /** copy elements from matrix into the fact data structure, in
   * preparation for factorization (or refactorization).  */
  virtual void copyMatrixElements ( double fact[], int lfact ) const;

  /** change format for row/column index matrices, in preparation for
   * call to MA27 FORTRAN routines
   *
   * @param irowM array of nnz elements indicating row index (in range
   * 1..n) of the corresponding matrix element
   *
   * @param jcolM array of nnz elements indicating col index (in range
   * 1..n) of the corresponding matrix element */
  virtual void getIndices( int irowM[], int jcolM[] ) const ;

  virtual void diagonalChanged( int idiag, int extent );
  virtual void matrixChanged();

protected:
  const SparseSymMatrix* mat;
  SparseStorageHandle mat_storage;

  void init();
  void freeWorkingArrays();
  bool checkErrorsAndReact();

  void scaleMatrix(); // TODO : implement..
  void orderMatrix(); // TODO : implement..
public:
  /** base class constructor. Allocates values for kTreatAsZero,
   * kThresholdPivoting, kThresholdPivotingMax,
   * kThresholdPivotingFactor, kPrecision, ipessimism,rpessimism */
  Ma27Solver( const SparseSymMatrix * sgm );

  virtual ~Ma27Solver();

  using DoubleLinearSolver::solve;
  void solve( OoqpVector& rhs ) override;
  void solve( int nrhss, double* rhss, int* colSparsity ) override;

};

#endif
