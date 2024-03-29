/*
 * sLinsysLeafMumps.C
 *
 *      Author: bzfrehfe
 */
#include <algorithm>
#include "sLinsysLeafMumps.h"
#include "../LinearSolvers/MumpsSolver/MumpsSolverLeaf.h"


static void appendCsrRow(const int* ia_src, const int* ja_src, const double* a_src, int row_src, int offset, int& nnz, int* ja_dest, double* a_dest) {
   assert(offset >= 0 && nnz >= 0);

   for (int c = ia_src[row_src]; c < ia_src[row_src + 1]; c++) {
      ja_dest[nnz] = ja_src[c] + offset;
      a_dest[nnz++] = a_src[c];
   }
}

sLinsysLeafMumps::~sLinsysLeafMumps() {
   delete[] schurRightNzColId;
   delete[] buffer;
}

void sLinsysLeafMumps::addTermToSparseSchurCompl(SparseSymmetricMatrix& SC) {
   this->addTermToSchurComplMumps(true, SC);
}

void sLinsysLeafMumps::addTermToDenseSchurCompl(DenseSymmetricMatrix& SC) {
   this->addTermToSchurComplMumps(false, SC);
}

// creates new (Fortran indexed, implicit CSC) matrix
//
//  (R^T A^T C^T)
//  (0          )
//  (F          )
//  (G          )
//
// ... but leave out empty rows (columns in the original right Schur factor)
void sLinsysLeafMumps::buildSchurRightMatrix(SymmetricMatrix& SC) {
   assert(prob);
   assert(!schurRightMatrix_csc);
   assert(!schurRightNzColId);
   assert(nSC == -1 && mSchurRight == -1 && nSchurRight == -1);

   const auto& A = data->getLocalA();
   const auto& C = data->getLocalC();
   const auto& F = data->getLocalF();
   const auto& G = data->getLocalG();
   const auto& R = data->getLocalCrossHessian();

   assert(R.numberOfNonZeros() == 0); // currently not supported

   mSchurRight = locnx + locmy + locmz;
   nSchurRight = SC.size();

   const int nnzSchurRight = R.numberOfNonZeros() + A.numberOfNonZeros() + C.numberOfNonZeros() + F.numberOfNonZeros() + G.numberOfNonZeros();

   if (!R.hasTransposed())
      R.updateTransposed();

   schurRightMatrix_csc = std::make_unique<SparseMatrix>(nSchurRight, mSchurRight, nnzSchurRight);

   int* const csc_ia = schurRightMatrix_csc->krowM();
   int* const csc_ja = schurRightMatrix_csc->jcolM();
   double* const csc_a = schurRightMatrix_csc->M();

   int M, nxA, nxC, nxParent;

   A.getSize(M, nxA);

   assert(M == locmy && locmyl >= 0 && locmzl >= 0);

   const int firstRowF = nSchurRight - locmyl - locmzl;
   const int firstRowG = nSchurRight - locmzl;

   C.getSize(M, nxC);

   nxParent = std::max(nxA, nxC);

   if (nxParent == -1)
      nxParent = 0;

   assert(nSchurRight >= nxParent);

   int nnz = 0;
   csc_ia[0] = 0;

   // todo new (variadic?) method in sparseStorage to concatenate matrices

   // add (R^T A^T C^T)
   if (nxParent > 0) {
      const int* const At_ia = A.getTranspose().krowM();
      const int* const At_ja = A.getTranspose().jcolM();
      const double* const At_a = A.getTranspose().M();
      const int* const Ct_ia = C.getTranspose().krowM();
      const int* const Ct_ja = C.getTranspose().jcolM();
      const double* const Ct_a = C.getTranspose().M();

#ifdef DEBUG_WRITE
      std::cout << "At" << std::endl;
      A.getTranspose().write_to_streamDense(cout);
      std::cout << "Ct" << std::endl;
      C.getTranspose().write_to_streamDense(cout);
#endif

      for (int r = 0; r < nxParent; r++) {
         appendCsrRow(At_ia, At_ja, At_a, r, locnx, nnz, csc_ja, csc_a);
         appendCsrRow(Ct_ia, Ct_ja, Ct_a, r, locnx + locmy, nnz, csc_ja, csc_a);
         csc_ia[r + 1] = nnz;
      }
   }

   // add (0 0 0)
   for (int r = nxParent; r < firstRowF; r++)
      csc_ia[r + 1] = csc_ia[nxParent];

   // add (F 0 0)
   if (locmyl > 0) {
      const int* const F_ia = F.krowM();
      const int* const F_ja = F.jcolM();
      const double* const F_a = F.M();

#ifdef DEBUG_WRITE
      std::cout << "F" << std::endl;
      F.write_to_streamDense(cout);
#endif

      for (int r = 0; r < locmyl; r++) {
         appendCsrRow(F_ia, F_ja, F_a, r, 0, nnz, csc_ja, csc_a);
         csc_ia[firstRowF + r + 1] = nnz;
      }
   }

   // add (G 0 0)
   if (locmzl > 0) {
      const int* const G_ia = G.krowM();
      const int* const G_ja = G.jcolM();
      const double* const G_a = G.M();

#ifdef DEBUG_WRITE
      std::cout << "G" << std::endl;
      G.write_to_streamDense(cout);
#endif

      for (int r = 0; r < locmzl; r++) {
         appendCsrRow(G_ia, G_ja, G_a, r, 0, nnz, csc_ja, csc_a);
         csc_ia[firstRowG + r + 1] = nnz;
      }
   }

   assert(nnz == nnzSchurRight);
   assert(nnz == csc_ia[nSchurRight]);

   schurRightMatrix_csc->deleteEmptyRows(schurRightNzColId);

#ifdef DEBUG_WRITE
   std::cout << "\n CSC" << std::endl;
   schurRightMatrix_csc->write_to_streamDense(cout);

   if( schurRightNzColId )
      for( int i = 0; i < schurRightMatrix_csc->getStorageRef().m; i++ )
         std::cout << schurRightNzColId[i] << std::endl;
#endif

   schurRightMatrix_csc->getStorageRef().c2fortran();
}


void sLinsysLeafMumps::addTermToSchurComplMumps(bool sparseSC, SymmetricMatrix& SC) {
   if (!schurRightMatrix_csc)
      buildSchurRightMatrix(SC);

   const int nNzRhs = schurRightMatrix_csc->getStorageRef().m;
   const int solSize = nNzRhs * mSchurRight;

   assert(solSize >= 1);

   if (bufferSize == -1) {
      assert(!buffer);

      if (solSize < bufferMaxSize)
         bufferSize = solSize;
      else
         bufferSize = bufferMaxSize;

      buffer = new double[bufferSize];
   }

   const int bufferNrhs = bufferSize / mSchurRight;
   const int nRuns = nNzRhs / bufferNrhs;
   const int leftoverNrhs = nNzRhs - bufferNrhs * nRuns;

   assert(buffer);
   assert(bufferNrhs >= 1);
   assert(nRuns * bufferNrhs + leftoverNrhs == nNzRhs);

   MumpsSolverLeaf* const solverMumps = dynamic_cast<MumpsSolverLeaf*>(solver.get());
   assert(solverMumps);

   const int n_empty =
      SC.size() - data->getLocalCrossHessian().getStorageRef().n - data->getLocalF().getStorageRef().m - data->getLocalG().getStorageRef().m;

   //  do block-wise computation of
   //
   //     SC +=  B^T K^-1 B
   //
   const BorderBiBlock border_left_transp(data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(),
         data->getLocalC().getTranspose(), n_empty, data->getLocalF(), data->getLocalG());

   for (int i = 0; i < nRuns; i++) {

      solverMumps->solve(*schurRightMatrix_csc, i * bufferNrhs, bufferNrhs, buffer);
      addLeftBorderTimesDenseColsToResTransp(border_left_transp, buffer, schurRightNzColId + i * bufferNrhs, mSchurRight, bufferNrhs, sparseSC, true,
            SC);

//      multLeftSchurComplBlocked(data, buffer, schurRightNzColId + i * bufferNrhs, bufferNrhs, sparseSC, SC);
   }

   if (leftoverNrhs > 0) {
      solverMumps->solve(*schurRightMatrix_csc, nRuns * bufferNrhs, leftoverNrhs, buffer);
      addLeftBorderTimesDenseColsToResTransp(border_left_transp, buffer, schurRightNzColId + nRuns * bufferNrhs, mSchurRight, leftoverNrhs, sparseSC,
            true, SC);
//      multLeftSchurComplBlocked(data, buffer, schurRightNzColId + nRuns * bufferNrhs, leftoverNrhs, sparseSC, SC);
   }


#if DEBUG_WRITE
   // debug stuff
   int myrank;
   static int iteration = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   ofstream myfile;
   char filename[50];
   sprintf(filename, "../A1_%d_%d.txt", myrank, iteration);
   myfile.open(filename);
   iteration++;
   SC.write_to_stream(myfile); // todo write out in each iteration with global counter and MPI rank!
   myfile.close();

   if( iteration >= 100)
      assert(0);
#endif
}
