/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef MA57LINSYS_H
#define MA57LINSYS_H

#include "../DoubleLinearSolver.h"
#include "../../LinearAlgebra/Dense/SimpleVector.h"
#include "../../LinearAlgebra/Sparse/SparseSymmetricMatrix.h"
#include "../../LinearAlgebra/Dense/DenseMatrix.h"

#include <vector>

// TODO : deprecated - we're not computing on BlueGene
#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif

extern "C" {
void FNAME(ma57id)(double cntl[], int icntl[]);

void
FNAME(ma57ad)(const int* n, const int* ne, int irn[], int jcn[], const int* lkeep, int keep[], int iwork[], int icntl[], int info[], double rinfo[]);

void FNAME(ma57bd)(const int* n, const int* ne, const double a[], double fact[], const int* lfact, int ifact[], const int* lifact, const int* lkeep,
      const int keep[], int work[], int* icntl, double cntl[], int info[], double rinfo[]);

void FNAME(ma57cd)(const int* job, const int* n, const double fact[], const int* lfact, const int ifact[], const int* lifact, const int* nrhs,
      double rhs[], const int* lrhs, double w[], const int* lw, int iw1[], int icntl[], int info[]);

void
FNAME(ma57dd)(const int* job, const int* n, const int* ne, const double a[], const int irn[], const int jcn[], const double fact[], const int* lfact,
      const int ifact[], const int* lifact, const double rhs[], double x[], double resid[], double w[], int iw[], int icntl[], double cntl[],
      int info[], double rinfo[]);

void FNAME(ma57ed)(const int* n, int* ic, int keep[], const double fact[], const int* lfact, double* newfac, const int* lnew, const int ifact[],
      const int* lifact, int newifc[], const int* linew, int* info);
}

/** implements the linear solver class using the HSL MA57 solver
 *
 * @ingroup LinearSolvers
 */
class Ma57Solver : public DoubleLinearSolver {

protected:
   SmartPointer<SparseStorage> mat_storage;

   /** control structures MA57 */
   std::vector<int> icntl = std::vector<int>(20);
   std::vector<int> info = std::vector<int>(40);
   std::vector<double> cntl = std::vector<double>(5);
   std::vector<double> rinfo = std::vector<double>(20);

   /** precision aimed for with iterative refinement */
   const double precision = 1e-5;

   const int max_tries = 4;
   const int n_iterative_refinement = 2;

   const int ooqp_print_level_warnings = 1000;

   /** the Threshold Pivoting parameter may need to be increased during
    * the algorithm if poor precision is obtained from the linear
    * solves. Larger values enforce greater stability in the factorization.
    */
   const double threshold_pivoting_max = 1.0;
   const double threshold_pivoting_factor = 10.0;
   const double default_pivoting = 0.01;

   /** During factorization entries smaller than small pivot will not
    * be accepted as pivots and the matrix will be treated as singular.
    *
    * This is the min we are willing to go with our pivots.
    */
   const double threshold_pivtol = 1e-20;
   const double threshold_pivtol_factor = 0.1;
   const double default_pivtol = 1.0e-7;

   /** index arrays for the factorization */
   std::vector<int> irowM, jcolM;

   /** nonzero element of the factors */
   std::vector<double> fact;

   /** vectors for solve */
   std::vector<double> x;
   std::vector<double> resid;

   /** dimension of the whole matrix */
   const int n;

   /** number of nonzeros in the matrix */
   int nnz;

   /** temporary storage */
   int lkeep{-1};
   std::vector<int> keep;

   /** temporary storage for the factorization process */
   int lifact = -1, lfact = -1;
   std::vector<int> ifact;

   std::vector<int> iworkn;
   std::vector<double> dworkn;

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

   bool freshFactor{false};

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
   [[nodiscard]] double getSmallPivot() const { return cntl[1]; }
   void setSmallPivot(double tol) { cntl[1] = tol; }

   [[nodiscard]] int minimumRealWorkspace() const { return info[4]; }
   [[nodiscard]] int minimumIntWorkspace() const { return info[5]; }

   void init();
   bool checkErrorsAndReact();

   void getIndices(std::vector<int>& irowM, std::vector<int>& jcolM) const;

public:
   explicit Ma57Solver(const SparseSymmetricMatrix* sgm, std::string name = "leaf");

   ~Ma57Solver() override = default;

   using DoubleLinearSolver::solve;
   void solve(Vector<double>& rhs) override;
   void solve(int nrhss, double* rhss, int*) override;
   void solve(GeneralMatrix& rhs) override;

   void diagonalChanged(int idiag, int extent) override;
   void matrixChanged() override;

   [[nodiscard]] bool reports_inertia() const override { return true; };
   [[nodiscard]] std::tuple<unsigned int, unsigned int, unsigned int> get_inertia() const override;
protected:
   void solve(int solveType, Vector<double>& rhs);
//   int* new_iworkn(int dim);
//   double* new_dworkn(int dim);
};

#endif
