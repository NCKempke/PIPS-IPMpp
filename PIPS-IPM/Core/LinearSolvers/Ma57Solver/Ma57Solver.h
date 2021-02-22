/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef MA57LINSYS_H
#define MA57LINSYS_H

#include "pipsport.h"

#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "SparseStorageHandle.h"
#include "OoqpVectorHandle.h"
#include "SparseStorage.h"


#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif

extern "C"
{
   void FNAME(ma57id)( double cntl[],  int icntl[] );

   void FNAME(ma57ad)( int * n,        int * ne,       int irn[],
      int jcn[],      int * lkeep,    int keep[],
      int iwork[],    int icntl[],    int info[],
      double rinfo[] );
   void FNAME(ma57bd)( int * n,        int * ne,       double a[],
      double fact[],  int * lfact,    int ifact[],
      int * lifact,   int * lkeep,    int keep[],
      int work[],     int * icntl,    double cntl[],
      int info[],     double rinfo[] );
   void FNAME(ma57cd)( int * job,      int * n,        double fact[],
      int * lfact,    int ifact[],    int * lifact,
      int * nrhs,     double rhs[],   int * lrhs,
      double w[],     int * lw,       int iw1[],
      int icntl[],    int info[]);
   void FNAME(ma57dd)( int * job,      int * n,        int * ne,
      double a[],     int irn[],      int jcn[],
      double fact[],  int * lfact,    int ifact[],
      int * lifact,   double rhs[],   double x[],
      double resid[], double w[],     int iw[],
      int icntl[],    double cntl[],  int info[],
      double rinfo[] );
   void FNAME(ma57ed)( int * n,        int * ic,       int keep[],
      double fact[],  int * lfact,    double * newfac,
      int * lnew,     int  ifact[],   int * lifact,
      int newifc[],   int * linew,    int * info );
}

/** implements the linear solver class using the HSL MA57 solver
 *
 * @ingroup LinearSolvers
 */
class Ma57Solver : public DoubleLinearSolver {
protected:
  int icntl[20];
  int info[40];
  double cntl[5];
  double rinfo[20];

  const int n_iterative_refinement = 10;
  const int max_tries = 8;

  const int ooqp_print_level_warnings = 10000;

  /** the Threshold Pivoting parameter, stored as U in the ma27dd
   *  common block. Takes values in the range [0,1]. Larger values
   *  enforce greater stability in the factorization as they insist on
   *  larger pivots. Smaller values preserve sparsity at the cost of
   *  using smaller pivots.  */
  double threshold_pivoting = 0.01;

  /** the Threshold Pivoting parameter may need to be increased during
   * the algorithm if poor precision is obtained from the linear
   * solves. threshold_pivoting indicates the largest value we are
   * willing to tolerate.  */
  double threshold_pivoting_max = 1.0;

  /** the factor in the range (1,inf) by which kThresholdPivoting is
   * increased when it is found to be inadequate.  */
  double threshold_pivoting_factor = 10.0;

  /** the "Treat As Zero" parameter, stored as pivtol in the common
   * block ma27td. The factorization will not accept a pivot whose
   * absolute value is less than this parameter as a 1x1 pivot or as
   * the off-diagonal in a 2x2 pivot.  */
  double treat_pivot_as_zero = 1e-12; // was 1e-10

  /** precision we demand from the linear system solver. If it isn't
   * attained on the first solve, we use iterative refinement and
   * possibly refactorization with a higher value of
   * kThresholdPivoting. */
  double precision = 1e-7;

  /** index array for the factorization */
  int *irowM = nullptr;
  int *jcolM = nullptr;

  /** storage for the original matrix */
  double *M;

  /** dimension of the whole matrix */
  const int n;

  /** number of nonzeros in the matrix */
  int nnz;

  /** temporary storage */
  int lkeep = 0;
  int *keep = nullptr;

  /** temporary storage for the factorization process */
  int lifact = 0, lfact = 0;
  int *ifact = nullptr;

  /* storage for the factors */
  double *fact = nullptr;

  /** amounts by which to increase allocated factorization space when
   * inadequate space is detected. ipessimism is for array "iw",
   * rpessimism is for the array "fact". */
  double ipessimism = 2;
  double rpessimism = 2;

  /** used to indicate when we need a fresh factorization (when
   * iterative refinement has failed to improve the precision of the
   * computed solution satisfactorily */
  bool freshFactor = false;

  /** store as a sparse symmetric matrix */
  SparseStorageHandle mStorage;

  /** called the very first time a matrix is factored. Allocates space
   * for the factorization and performs ordering */
  virtual void firstCall();
public:
  /** sets mStorage to refer to the argument sgm */
  Ma57Solver( SparseSymMatrix * sgm );

  void diagonalChanged( int idiag, int extent ) override;
  void matrixChanged() override;
  using DoubleLinearSolver::solve;
  void solve( OoqpVector& rhs ) override;
  void solve( GenMatrix& rhs) override;
  void solve ( int nrhss, double* rhss, int* ) override;

  //virtual void Lsolve  ( OoqpVector& x );
  //virtual void Dsolve  ( OoqpVector& x );
  //virtual void Ltsolve ( OoqpVector& x );
  //virtual void Refine  ( OoqpVector& x );
 protected:
  void solve(int solveType, OoqpVector& rhs);
  void freeWorkingArrays();
  void init();
  bool checkErrorsAndReact();

  int* iworkn, niworkn;
  int* new_iworkn(int dim);
  double* dworkn;
  int ndworkn;
  double* new_dworkn(int dim);

  /** destructor */
  virtual ~Ma57Solver();

};

#endif
