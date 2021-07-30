/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef MA27LINSYS_H
#define MA27LINSYS_H

#include "../DoubleLinearSolver.h"

#include "../LinearAlgebra/Sparse/SparseSymmetricMatrix.h"
#include "../LinearAlgebra/Dense/DenseVector.hpp"

#include <vector>
#include <string>


#ifndef FNAME
#define FNAME(f) f ## _
#endif

extern "C" {
void FNAME(ma27id)(int icntl[], double cntl[]);

void FNAME(ma27ad)(const int* n, const int* nz, int irn[], int icn[], int iw[], int* liw, int ikeep[], int iw1[], int* nsteps, int* iflag, int icntl[],
   double cntl[], int info[], double* ops);

void FNAME(ma27bd)(const int* n, const int* nz, int irn[], int icn[], double a[], int* la, int iw[], int* liw, int ikeep[], int* nsteps, int* maxfrt,
   int iw1[], int icntl[], double cntl[], int info[]);

void
FNAME(ma27cd)(const int* n, const double a[], const int* la, const int iw[], const int* liw, double w[], const int* maxfrt, double rhs[], int iw1[],
   const int* nsteps, const int icntl[], int info[]);

}

/** implements the linear solver class using the HSL MA27 solver
 *
 * @ingroup LinearSolvers
 */
class Ma27Solver : public DoubleLinearSolver {

protected:
   const SparseStorage* mat_storage;

   /** control structures MA27 */
   std::vector<int> icntl = std::vector<int>(30);
   std::vector<int> info = std::vector<int>(20);
   std::vector<double> cntl = std::vector<double>(5);

   /** precision aimed for with iterative refinement */
   double precision = 1e-7;

   const int max_tries = 8;
   int max_n_iter_refinement = 5;

   const int ooqp_print_level_warnings = 1000;

   /** the Threshold Pivoting parameter may need to be increased during
    * the algorithm if poor precision is obtained from the linear
    * solves. Larger values enforce greater stability in the factorization.
    */
   double threshold_pivoting_max = 0.5;
   const double threshold_pivoting_factor = 10;
   const double default_pivoting = 0.5;

   /** During factorization entries smaller than small pivot will not
    * be accepted as pivots and the matrix will be treated as singular.
    *
    * This is the min we are willing to go with our pivots.
    */
   const double threshold_pivtol = 1e-8;
   const double threshold_pivtol_factor = 0.1;
   const double default_pivtol = 1.0e-7;

   /* detecting dense rows during the factorization to preserve sparsity */
   const double default_fratio = 0.5;

   /** index arrays for the factorization */
   std::vector<int> irowM, jcolM;

   /** nonzero element of the factors */
   std::vector<double> fact;

   /** vectors for iterative refinement */
   std::vector<double> iter;
   std::vector<double> iter_best;
   std::vector<double> resid;

   /** dimension of the whole matrix */
   const int n;

   /** number of nonzeros in the matrix */
   int nnz;

   /** length of the array containing factors; may be increased during
    * the numerical factorization if the estimate obtained during the
    * symbolic factorization proves to be inadequate. */
   int la{-1};

   /** pivot sequence and temporary storage information */
   int* ikeep{}, * iw{}, liw{}, * iw1{}, nsteps{-1}, maxfrt{-1};

   /** temporary storage for the solve stage */
   double* ww{};

   /** amounts by which to increase allocated factorization space when
    * inadequate space is detected. ipessimism is for array "iw",
    * rpessimism is for the array "fact". */
   double ipessimism = 3.0;
   double rpessimism = 3.0;

   /** amount of solvers that might call this routine - used to duplicate data structures
    * where ever necessary to make the solve calls thread safe
    */
   const int n_threads{1};

   /** to distinguish between root and leaf solver in error messages */
   const std::string name;

   /** called the very first time a matrix is factored. Allocates space
    * for the factorization and performs ordering
    */
   virtual void firstCall();

   /** the Threshold Pivoting parameter, stored as U in the ma27dd
    *  common block. Takes values in the range [-0.5,0.5]. If negative no
    *  numerical pivoting will be performed - if positive numerical pivoting
    *  will be performed. If no pivoting is applied the subroutine will fail when
    *  a zero pivot is encountered. Larger values enforce greater stability in the
    *  factorization as they insist on
    *  larger pivots. Smaller values preserve sparsity at the cost of
    *  using smaller pivots.
    */
   [[nodiscard]] double thresholdPivoting() const { return cntl[0]; }
   void setThresholdPivoting(double piv) { cntl[0] = piv; }

   /** the "small pivot" parameter, stored as pivtol in the common
    * block ma27td. The factorization will not accept a pivot whose
    * absolute value is less than this parameter as a 1x1 pivot or as
    * the off-diagonal in a 2x2 pivot.  */
   [[nodiscard]] double getSmallPivot() const { return cntl[2]; }
   void setSmallPivot(double tol) { cntl[2] = tol; }

   void setFratio(double ratio) { cntl[1] = ratio; };

   [[nodiscard]] int minimumRealWorkspace() const { return info[4]; }
   [[nodiscard]] int minimumIntWorkspace() const { return info[5]; }

   void init();
   void freeWorkingArrays();
   bool checkErrorsAndReact();

   /* scaler */
   SymmetricLinearScaler* scaler{};

   /* stuff for MA60 iterative refinement */
   int icntl_ma60[5]{};
   int keep_ma60[10]{};
   double rkeep_ma60[10]{};

   std::vector<double> w_ma60;
   std::vector<int> iw_ma60;

   void copyMatrixElements(std::vector<double>& afact, int lfact) const;
   void getIndices(std::vector<int>& irow, std::vector<int>& jcol) const;

//  void orderMatrix(); // TODO : implement..
public:
   explicit Ma27Solver(const SparseSymmetricMatrix& sgm, std::string name_ = "leaf");

   ~Ma27Solver() override;

   using DoubleLinearSolver::solve;
   /* thread-safe if not called by more OMP_NUM_THREADS than available when calling the ctor */
   void solve(Vector<double>& rhs) override;
   void solve(int nrhss, double* rhss, int* colSparsity) override;

   void diagonalChanged(int idiag, int extent) override;
   void matrixChanged() override;

   [[nodiscard]] bool reports_inertia() const override { return true; };
   [[nodiscard]] std::tuple<unsigned int, unsigned int, unsigned int> get_inertia() const override;

};

#endif
