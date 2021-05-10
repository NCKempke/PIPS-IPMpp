/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
#include <iostream>
#include "PardisoSchurSolver.h"
#include "SparseSymmetricMatrix.h"
#include "DenseMatrix.h"
#include "PIPSIPMppOptions.h"
#include "pipsdef.h"
#include <algorithm>

#ifdef STOCH_TESTING
extern double g_scenNum;
static int rhsCount=0;
#endif

#ifdef STOCH_TESTING
int dumpAugMatrix(int n, int nnz, int nSys, //size, nnz and size of the (1,1) block
        double* eltsA,
        int* rowptr,
        int* colidx,
        const char* fname=nullptr);
int dumpSysMatrix(SparseSymMatrix* Msys,
                  const char* fname=nullptr);
int dumpRhs(SimpleVector<double>& v);
#endif

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

#define SHRINK_SC

PardisoSchurSolver::PardisoSchurSolver(const SparseSymmetricMatrix* sgm) : Msys{sgm} {
   const int myRank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   useSparseRhs = pipsipmpp_options::get_bool_parameter("PARDISO_SPARSE_RHS_LEAF");

   symbFactorInterval = pipsipmpp_options::get_int_parameter("PARDISO_SYMB_INTERVAL");
   if (symbFactorInterval < 0)
      symbFactorInterval = symbFactorIntervalDefault;

   pivotPerturbationExp = pipsipmpp_options::get_int_parameter("PARDISO_PIVOT_PERTURBATION");
   if (pivotPerturbationExp < 0)
      pivotPerturbationExp = pivotPerturbationExpDefault;

   nIterativeRefins = pipsipmpp_options::get_int_parameter("PARDISO_NITERATIVE_REFINS");
   if (nIterativeRefins < 0)
      nIterativeRefins = nIterativeRefinsDefault;

   parallelForwardBackward = parallelForwardBackwardDefault;

   factorizationTwoLevel = factorizationTwoLevelDefault;

   static bool printed = false;

   if (myRank == 0 && !printed) {
      printed = true;
      printf(" using pivot perturbation 10^-%d \n", pivotPerturbationExp);

      printf(" using maximum of %d iterative refinements  \n", nIterativeRefins);

      if (useSparseRhs)
         printf(" using PARDISO_SPARSE_RHS_LEAF \n");
      else
         printf(" NOT using PARDISO_SPARSE_RHS_LEAF \n");
   }
}

void PardisoSchurSolver::firstCall() {
   iparm[0] = 0; /* make init set iparm to default values */

   initPardiso();
   setIparm(iparm);
}

// this function is called only once and creates the augmented system
void
PardisoSchurSolver::firstSolveCall( const SparseMatrix& R, const SparseMatrix& A, const SparseMatrix& C, const SparseMatrix& F, const SparseMatrix& G, int nSC0) {
   int nR, nA, nC, nF, nG, nx;
   nnz = 0;

   F.getSize(nF, nx);
   nnz += F.numberOfNonZeros();
   G.getSize(nG, nx);
   nnz += G.numberOfNonZeros();
   R.getSize(nR, nx);
   nnz += R.numberOfNonZeros();
   A.getSize(nA, nx);
   nnz += A.numberOfNonZeros();
   C.getSize(nC, nx);
   nnz += C.numberOfNonZeros();
   const int Msize = static_cast<int>(Msys->size());

   if (nR == 0)
      assert(R.numberOfNonZeros() == 0);
   if (nA == 0)
      assert(A.numberOfNonZeros() == 0);
   if (nC == 0)
      assert(C.numberOfNonZeros() == 0);
   if (nF == 0)
      assert(F.numberOfNonZeros() == 0);
   if (nG == 0)
      assert(G.numberOfNonZeros() == 0);

   assert(F.getStorageRef().isValid());
   assert(G.getStorageRef().isValid());
   assert(R.getStorageRef().isValid());
   assert(A.getStorageRef().isValid());
   assert(C.getStorageRef().isValid());
   assert(Msys->getStorageRef().isValid());

   // todo not implemented yet
   assert(R.numberOfNonZeros() == 0);

   if (nF > 0 || nG > 0)
      nSC = nSC0;
   else
      nSC = nx;

#ifdef TIMING
   cout << "firstSolveCall: nR=" << nR << " nA=" << nA << " nC=" << nC << " nF=" << nF << " nG=" << nG << " nSC=" << nSC << " sizeKi=" << (nR+nA+nC)<< endl
       << " nnzR=" << R.numberOfNonZeros()
       << " nnzA=" << A.numberOfNonZeros()
       << " nnzC=" << C.numberOfNonZeros()
       << " nnzF=" << F.numberOfNonZeros()
       << " nnzG=" << C.numberOfNonZeros() << endl;
#endif

   n = nR + nA + nC + nSC;

   assert(Msize == nR + nA + nC);
   nnz += Msys->numberOfNonZeros();
   nnz += nSC; //space for the 0 diagonal of 2x2 block

   // the lower triangular part of the augmented system in row-major
   SparseSymmetricMatrix augSys(n, nnz);

   // pointer for augmented system
   int* krowAug = augSys.getStorageRef().krowM;
   int* jcolAug = augSys.getStorageRef().jcolM;
   double* MAug = augSys.getStorageRef().M;

   //
   //put (1,1) block in the augmented system
   //
   memcpy(krowAug, Msys->getStorageRef().krowM, sizeof(int) * Msize);
   memcpy(jcolAug, Msys->getStorageRef().jcolM, sizeof(int) * Msys->numberOfNonZeros());
   memcpy(MAug, Msys->getStorageRef().M, sizeof(double) * Msys->numberOfNonZeros());


   int nnzIt = Msys->numberOfNonZeros();
   //
   //put A and C block in the augmented system as At and Ct in the lower triangular part
   //

   if (nA > 0 || nC > 0 || nF > 0 || nG > 0) {
      const bool putA = A.numberOfNonZeros() > 0;
      const bool putC = C.numberOfNonZeros() > 0;

      // putA = TRUE => nA > 0
      assert(putA == (putA && (nA > 0)));
      assert(putC == (putC && (nC > 0)));

      // initialize variables for At
      SparseMatrix At(putA ? nx : 0, putA ? nA : 0, putA ? A.numberOfNonZeros() : 0);
      int* krowAt = At.getStorageRef().krowM;
      int* jcolAt = At.getStorageRef().jcolM;
      double* MAt = At.getStorageRef().M;

      if (putA)
         A.getStorageRef().transpose(krowAt, jcolAt, MAt);

      const int colShiftA = nR;

      // initialize variables for Ct
      SparseMatrix Ct(putC ? nx : 0, putC ? nC : 0, putC ? C.numberOfNonZeros() : 0);
      int* krowCt = Ct.getStorageRef().krowM;
      int* jcolCt = Ct.getStorageRef().jcolM;
      double* MCt = Ct.getStorageRef().M;

      if (putC)
         C.getStorageRef().transpose(krowCt, jcolCt, MCt);

      const int colShiftC = nR + nA;

      int row = Msize;
      for (; row < Msize + nx; row++) {
         krowAug[row] = nnzIt;

         if (putA) {
            for (int c = krowAt[row - Msize]; c < krowAt[row - Msize + 1]; c++) {
               const int j = jcolAt[c];
               jcolAug[nnzIt] = j + colShiftA;
               MAug[nnzIt] = MAt[c];
               nnzIt++;
            }
         }

         if (putC) {
            for (int c = krowCt[row - Msize]; c < krowCt[row - Msize + 1]; c++) {
               const int j = jcolCt[c];
               jcolAug[nnzIt] = j + colShiftC;
               MAug[nnzIt] = MCt[c];
               nnzIt++;
            }
         }
         //add the zero from the diagonal
         jcolAug[nnzIt] = row;
         MAug[nnzIt] = 0.0;
         nnzIt++;

      }
      krowAug[row] = nnzIt;
   }
   assert(nnzIt == Msys->numberOfNonZeros() + A.numberOfNonZeros() + C.numberOfNonZeros() + nx);

   //
   // add linking constraint matrices F and G
   //
   if (nF > 0 || nG > 0) {
      // put diagonal in zero block
      int row = Msize + nx;
      assert(row <= n - nF - nG);
      for (; row < n - nF - nG; row++) {
         krowAug[row] = nnzIt;
         jcolAug[nnzIt] = row;
         MAug[nnzIt] = 0.0;
         nnzIt++;
      }

      // are there linking equality constraints?
      if (nF > 0) {
         // put F in the lower triangular part below R (and below 0 block)

         int* krowF = F.getStorageRef().krowM;
         int* jcolF = F.getStorageRef().jcolM;
         double* MF = F.getStorageRef().M;

         const bool putF = F.numberOfNonZeros() > 0;

         const int row0 = n - nF - nG;
         int row_F = row0;
         assert(row0 >= Msize + nx);

         for (; row_F < n - nG; row_F++) {
            krowAug[row_F] = nnzIt;

            if (putF) {
               for (int c = krowF[row_F - row0]; c < krowF[row_F - row0 + 1]; c++) {
                  jcolAug[nnzIt] = jcolF[c];
                  MAug[nnzIt] = MF[c];
                  nnzIt++;
               }
            }
            //add the zero from the diagonal
            jcolAug[nnzIt] = row_F;
            MAug[nnzIt] = 0.0;
            nnzIt++;
         }
         krowAug[row_F] = nnzIt;
      }

      // are there linking equality constraints?
      if (nG > 0) {
         // put G in the lower triangular part below F

         int* krowG = G.getStorageRef().krowM;
         int* jcolG = G.getStorageRef().jcolM;
         double* MG = G.getStorageRef().M;

         const bool putG = G.numberOfNonZeros() > 0;

         const int row0 = n - nG;
         int row_G = row0;

         for (; row_G < n; row_G++) {
            krowAug[row_G] = nnzIt;

            if (putG) {
               for (int c = krowG[row_G - row0]; c < krowG[row_G - row0 + 1]; c++) {
                  jcolAug[nnzIt] = jcolG[c];
                  MAug[nnzIt] = MG[c];
                  nnzIt++;
               }
            }
            //add the zero from the diagonal
            jcolAug[nnzIt] = row_G;
            MAug[nnzIt] = 0.0;
            nnzIt++;
         }
         krowAug[row_G] = nnzIt;
      }
   }
   assert(nnzIt == Msys->numberOfNonZeros() + A.numberOfNonZeros() + C.numberOfNonZeros() + F.numberOfNonZeros() + G.numberOfNonZeros() + nSC);

#ifdef SHRINK_SC

   // remove empty or zero rows from augSys and check that in Msys none are removed!
   int* shrinked2orgAug = nullptr;

   augSys.deleteZeroRowsCols(shrinked2orgAug);

   n = static_cast<int>(augSys.size());

#ifndef NDEBUG
   for (int i = 0; i < Msize; i++) {
      if (i != shrinked2orgAug[i]) {
         std::cout << "zero row in (1,1) block of Schur complement!" << std::endl;
         std::cout << "i=" << i << " shrinked2orgAug[i]=" << shrinked2orgAug[i] << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }

      if (Msys->getStorageRef().krowM[i] == Msys->getStorageRef().krowM[i + 1]) {
         std::cout << "(2) zero row in (1,1) block of Schur complement!" << std::endl;
         MPI_Abort(MPI_COMM_WORLD, 1);
      }
   }
#endif

   nSC = static_cast<int>(augSys.size()) - Msize;

   assert(shrinked2orgSC == nullptr);
   shrinked2orgSC = new int[nSC];

   // shrinked2orgSC should only map Schur complement part of augmented saddle-point system
   for (int i = 0; i < nSC; i++) {
      shrinked2orgSC[i] = shrinked2orgAug[i + Msize] - Msize;
      assert(shrinked2orgSC[i] >= i);
   }

   delete[] shrinked2orgAug;
#endif

   nnz = augSys.numberOfNonZeros();

   // we need to transpose to get the augmented system in the row-major upper triangular format of PARDISO
   rowptrAug = new int[n + 1];
   colidxAug = new int[nnz];
   eltsAug = new double[nnz];

   augSys.getStorageRef().transpose(rowptrAug, colidxAug, eltsAug);

   assert(rowptrAug[n] == nnz);

   //save the indices for diagonal entries for a streamlined later update
   int* krowMsys = Msys->getStorageRef().krowM;
   int* jcolMsys = Msys->getStorageRef().jcolM;

   for (int r = 0; r < Msize; r++) {

      // Msys - find the index in jcol for the diagonal (r,r)
      int idxDiagMsys = -1;
      for (int idx = krowMsys[r]; idx < krowMsys[r + 1]; idx++)
         if (jcolMsys[idx] == r) {
            idxDiagMsys = idx;
            break;
         }

      assert(idxDiagMsys >= 0);

      // aug  - find the index in jcol for the diagonal (r,r)
      int idxDiagAug = -1;
      for (int idx = rowptrAug[r]; idx < rowptrAug[r + 1]; idx++)
         if (colidxAug[idx] == r) {
            idxDiagAug = idx;
            break;
         }

      assert(idxDiagAug >= 0);

      diagMap.insert(std::pair<int, int>(idxDiagMsys, idxDiagAug));
   }

   //convert to Fortran indexing
   for (int it = 0; it < n + 1; it++)
      rowptrAug[it]++;
   for (int it = 0; it < nnz; it++)
      colidxAug[it]++;

   //allocate temp vector(s)
   nvec = new double[n];
   nvec2 = new double[n];
   nvec_size = n;
}

bool PardisoSchurSolver::iparmUnchanged() {
   /* put all Parameters that should stay be checked against init into this array */
   static const int check_iparm[] = {1, 30, 10, 12, 23, 24};

   bool unchanged = true;
   bool print = false;

   int iparm_compare[64];
   setIparm(iparm_compare);


   std::vector<int> to_compare(check_iparm, check_iparm + sizeof(check_iparm) / sizeof(check_iparm[0]));

   for (int i = 0; i < 64; ++i) {
      // if entry should be compared
      if (std::find(to_compare.begin(), to_compare.end(), i) != to_compare.end()) {
         if (iparm[i] != iparm_compare[i]) {
            if (print)
               std::cout << "ERROR - PardisoSolver: elements in iparm changed at " << i << ": " << iparm[i] << " != " << iparm_compare[i] << "(new)"
                         << std::endl;
            unchanged = false;
         }
      }
   }
   return unchanged;
}

void PardisoSchurSolver::diagonalChanged(int /* idiag */, int /* extent */) {
   this->matrixChanged();
}

void PardisoSchurSolver::matrixChanged() {
   if (first) {
      firstCall();
      first = false;
   }

   // we don't have the right hand-size, therefore we can't (re)factorize
   // the augmented system at this point.

   //dumpSysMatrix(Msys);
}

void PardisoSchurSolver::schur_solve(const SparseMatrix& R, const SparseMatrix& A, const SparseMatrix& C, const SparseMatrix& F, const SparseMatrix& G,
      DenseSymmetricMatrix& SC0) {
   int* rowptrSC = nullptr;
   int* colidxSC = nullptr;
   double* eltsSC = nullptr;

   computeSC(static_cast<int>(SC0.size()), R, A, C, F, G, rowptrSC, colidxSC, eltsSC);

#ifdef SHRINK_SC
   for (int r = 0; r < nSC; r++) {
      const int r_org = shrinked2orgSC[r];
      assert(r_org >= r && r_org < SC0.size());

      for (int ci = rowptrSC[r]; ci < rowptrSC[r + 1]; ci++) {
         const int c = colidxSC[ci];
         const int c_org = shrinked2orgSC[c];

         assert(c >= r);
         assert(c_org >= c && c_org < SC0.size());

         SC0[c_org][r_org] += eltsSC[ci];

         // NOTE: we only save half of the matrix, so we don't need
         // if( r_org != c_org ) SC0[r_org][c_org] += eltsSC[ci];
      }
   }
#else
   for( int r = 0; r < nSC; r++ )
   {
      for( int ci = rowptrSC[r]; ci < rowptrSC[r + 1]; ci++ )
      {
         const int c = colidxSC[ci];
         assert(c >= r);

         SC0[c][r] += eltsSC[ci];

      }
   }
#endif

   delete[] rowptrSC;
   delete[] colidxSC;
   delete[] eltsSC;
}


void PardisoSchurSolver::schur_solve_sparse(const SparseMatrix& R, const SparseMatrix& A, const SparseMatrix& C, const SparseMatrix& F, const SparseMatrix& G,
      SparseSymmetricMatrix& SC0) {
   int* rowptrSC = nullptr;
   int* colidxSC = nullptr;
   double* eltsSC = nullptr;

   computeSC(static_cast<int>(SC0.size()), R, A, C, F, G, rowptrSC, colidxSC, eltsSC);

   int* rowptrBase = SC0.krowM();
   int* colidxBase = SC0.jcolM();
   double* eltsBase = SC0.M();

#ifdef SHRINK_SC
   assert(SC0.size() >= nSC);

   // add to summed Schur complement
   for (int r = 0; r < nSC; r++) {
      const int r_org = shrinked2orgSC[r];
      assert(r_org >= r && r_org < SC0.size());

      int cbase = rowptrBase[r_org];

      // catch empty diagonal
      if (rowptrSC[r + 1] - rowptrSC[r] == 1 && eltsSC[rowptrSC[r]] == 0.0)
         continue;

      for (int j = rowptrSC[r]; j < rowptrSC[r + 1]; j++) {
         const int c_org = shrinked2orgSC[colidxSC[j]];
         assert(c_org >= colidxSC[j] && c_org < SC0.size());

         while (colidxBase[cbase] != c_org) {
            cbase++;
            assert(cbase < rowptrBase[r_org + 1]);
         }

         eltsBase[cbase] += eltsSC[j];
      }
   }
#else
   assert(SC0.size() == nSC);


   // add to summed Schur complement
   for( int r = 0; r < nSC; r++ )
   {
      int cbase = rowptrBase[r];

      // catch empty diagonal
      if( rowptrSC[r + 1] - rowptrSC[r] == 1 && eltsSC[rowptrSC[r]] == 0.0 )
         continue;

      for( int j = rowptrSC[r]; j < rowptrSC[r + 1]; j++ )
      {
         const int c = colidxSC[j];

         while( colidxBase[cbase] != c )
         {
            cbase++;
            assert(cbase < rowptrBase[r + 1]);
         }

         eltsBase[cbase] += eltsSC[j];
      }
   }
#endif

   delete[] rowptrSC;
   delete[] colidxSC;
   delete[] eltsSC;
}

PardisoSchurSolver::~PardisoSchurSolver() {
   delete[] rowptrAug;
   delete[] colidxAug;
   delete[] eltsAug;
   delete[] shrinked2orgSC;
   delete[] nvec;
   delete[] nvec2;
}

#ifdef STOCH_TESTING
int dumpAugMatrix(int n, int nnz, int nSys,
        double* elts, int* rowptr, int* colidx,
        const char* fname)
{
  char filename[1024];
  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if(fname==nullptr) 
    sprintf(filename, "augMat-r%d-i%g-s%g.dat", myRank, g_iterNumber,  g_scenNum+1);
  else 
    sprintf(filename, "%s-r%d-i%g-s%g.dat", fname, myRank, g_iterNumber,  g_scenNum+1);

  cout << "saving to:" << filename << " ...";

  ofstream fd(filename);
  fd << scientific;
  fd.precision(16);

  fd << n << endl << nSys << endl << nnz << endl;
  for(int it=0; it<n+1; it++) {fd << rowptr[it] << " ";}  fd << endl;
  for(int it=0; it<nnz; it++) {fd << colidx[it] << " ";}  fd << endl;
  for(int it=0; it<nnz; it++) {fd << elts[it]   << " ";}  fd << endl;

  cout << filename << " done!" << endl;
  return 0;
}

int dumpSysMatrix(SparseSymMatrix* Msys, const char* fname)
{
  char filename[1024];
  if(fname==nullptr) 
    sprintf(filename, "Qdump-%g-s%g.dat", g_iterNumber,  g_scenNum+1);
  else 
    sprintf(filename, "%s-%g-s%g.dat", fname, g_iterNumber,  g_scenNum+1);

  cout << "saving to:" << filename << " ...";

  int n  =Msys->size();
  int nnz=Msys->numberOfNonZeros();

  // we need to transpose to get the augmented system in the row-major upper triangular format of  PARDISO 
  int* rowptr  = new int[n+1];
  int* colidx  = new int[nnz];
  double* elts = new double[nnz];
  Msys->getStorageRef().transpose(rowptr,colidx,elts);


  ofstream fd(filename);
  fd << scientific;
  fd.precision(16);

  fd << n << endl << nnz << endl;
  for(int it=0; it<n+1; it++) {fd << rowptr[it]+1 << " ";}  fd << endl;
  for(int it=0; it<nnz; it++) {fd << colidx[it]+1 << " ";}  fd << endl;
  for(int it=0; it<nnz; it++) {fd << elts[it]   << " ";}  fd << endl;

  delete[] rowptr; delete[] colidx; delete[] elts;

  cout << " Done!" << endl;
  

  return 0;
}
int dumpRhs(SimpleVector<double>& v)
{
  rhsCount++;
  char filename[1024];
  sprintf(filename, "rhsDump-%g-%d.dat", g_iterNumber,   rhsCount);
  cout << "saving to: " << filename << " ...";

  ofstream fd(filename);
  fd << scientific;
  fd.precision(16);

  fd << v.length() << endl;
  for(int i=0; i<v.length(); i++) fd << v[i]   << " ";  
  fd << endl;

  cout << "done!" << endl;
  return 0;
}
#endif
