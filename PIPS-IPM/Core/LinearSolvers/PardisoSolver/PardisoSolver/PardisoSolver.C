/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin, Murat Mut
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "PardisoSolver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"
#include "pipsdef.h"
#include "pipsport.h"
#include <cstdlib>

#include "mpi.h"

PardisoSolver::PardisoSolver( SparseSymMatrix * sgm )
{
  Msys = sgm;
  Mdsys = nullptr;
  n = sgm->size();
  nnz = sgm->numberOfNonZeros();

  krowM = new int[n+1];
  jcolM = new int[nnz];
  M = new double[nnz];

  sol = nullptr;
  sz_sol = 0;

  first = true;
  nvec = new double[n];

  mtype = -2;
  error = 0;
  phase = 0;
  nrhs = -1; // do not specify here
  maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
  mnum = 1; // actual matrix (as in index from 1 to maxfct)
  msglvl = 0; // messaging level (0 = no output, 1 = statistical info to screen)
}


PardisoSolver::PardisoSolver( DenseSymMatrix * m )
{
  Msys = nullptr;
  Mdsys = m;
  n = m->size();
  
  nnz=0; //getNumberOfNonZeros(*Mdsys);

  krowM = nullptr;//new int[n+1];
  jcolM = nullptr;//new int[nnz];
  M = nullptr;//new double[nnz];

  sol = nullptr;
  sz_sol = 0;

  first = true;
  nvec = new double[n];

  mtype = -2;
  error = 0;
  nrhs = -1;
  phase = 0;
  maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
  mnum = 1; // actual matrix (as in index from 1 to maxfct)
  msglvl = 0; // messaging level (0 = no output, 1 = statistical info to screen)
}

bool PardisoSolver::iparmUnchanged() const
{
   bool same = true;
   bool print = true;

   std::vector<int> iparm_compare(64, 0);
   setIparm( iparm_compare.data() );

   static const std::vector<int> to_compare{ 1, 7, 10, 12 };

   for(int i = 0; i < 64; ++i)
   {
      // if entry should be compared
      if(std::find(to_compare.begin(), to_compare.end(), i) != to_compare.end())
      {
         if(iparm[i] != iparm_compare[i])
         {
            if(print)
               std::cout << "ERROR - PardisoSolver: elements in iparm changed at " << i << ": "
                  << iparm[i] << " != " << iparm_compare[i] << "(new)" << std::endl;
            same = false;
         }
      }

   }
   return same;
}

void PardisoSolver::initSystem()
{
   if( Msys )
   {
      //get the matrix in upper triangular
      Msys->getStorageRef().transpose(krowM, jcolM, M);

      //save the indices for diagonal entries for a streamlined later update
      int* krowMsys = Msys->getStorageRef().krowM;
      int* jcolMsys = Msys->getStorageRef().jcolM;
      for( int r = 0; r < n; r++ )
      {
         // Msys - find the index in jcol for the diagonal (r,r)
         int idxDiagMsys = -1;
         for( int idx = krowMsys[r]; idx < krowMsys[r + 1]; idx++ )
            if( jcolMsys[idx] == r )
            {
               idxDiagMsys = idx;
               break;
            }
         assert(idxDiagMsys >= 0);
         //must have all diagonals entries

         // aug - find the index in jcol for the diagonal (r,r)
         int idxDiagAug = -1;
         for( int idx = krowM[r]; idx < krowM[r + 1]; idx++ )
            if( jcolM[idx] == r )
            {
               idxDiagAug = idx;
               break;
            }
         assert(idxDiagAug>=0);

         diagMap.insert(std::pair<int, int>(idxDiagMsys, idxDiagAug));
      }
   }
   else if( Mdsys )
   {
      // the input is a dense matrix
      // the dense matrix is also processed in matrixChanged everytime the method is called

      // the input is a dense matrix
      DenseSymMatrix& Md = (*Mdsys);
      nnz = Md.getNumberOfNonZeros();

      delete[] krowM;
      delete[] jcolM;
      delete[] M;
      krowM = new int[n + 1];
      jcolM = new int[nnz];
      M = new double[nnz];

      int nzIt = 0;
      for( int i = 0; i < n; i++ )
      {

         krowM[i] = nzIt;
         //cout << i << " " << krowM[i] << endl;

         jcolM[nzIt] = i; // the diagonal
         M[nzIt] = Md[i][i];
         assert(nzIt<nnz);
         nzIt++;

         for( int j = i + 1; j < n; j++ )
         {
            if( Md[i][j] != 0.0 )
            {
               jcolM[nzIt] = j;
               M[nzIt] = Md[i][j];
               assert(nzIt<nnz);
               nzIt++;
            }
         }
      }
      //cout << "PardisoSolver::first call nzit=" << nzIt << endl;
      assert(nzIt==nnz);
      krowM[n] = nzIt;

   }
   else
   {
      assert(false);
   }

   // need Fortran indexes
   // extra method?
   for( int i = 0; i < n + 1; i++ )
      krowM[i] += 1;
   for( int i = 0; i < nnz; i++ )
      jcolM[i] += 1;

} 

void PardisoSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

void PardisoSolver::matrixChanged()
{
   if( first )
   {
      firstCall();
      first = false;
   }
   else
   {
      if( Msys )
      {
         //update diagonal entries in the PARDISO aug sys (if the input is sparse)
         double* eltsMsys = Msys->getStorageRef().M;
         std::map<int, int>::iterator it;
         for( it = diagMap.begin(); it != diagMap.end(); it++ )
            M[it->second] = eltsMsys[it->first];
      }
   }

   if( Mdsys )
   {
      // the input is a dense matrix
      DenseSymMatrix& Md = (*Mdsys);
      //double tm=MPI_Wtime();
      int nzIt = 0;
      for( int i = 0; i < n; i++ )
      {

         krowM[i] = nzIt;

         jcolM[nzIt] = i; // the diagonal
         M[nzIt] = Md[i][i];
         nzIt++;

         for( int j = i + 1; j < n; j++ )
         {
            if( Md[i][j] != 0.0 )
            {
               jcolM[nzIt] = j;
               M[nzIt] = Md[i][j];
               nzIt++;
            }
         }
      }
      //cout << "PardisoSolver::matrix changed nnzit=" << nzIt << endl;
      krowM[n] = nzIt;
      assert(nzIt==nnz);
      // need Fortran indexes
      for( int i = 0; i < n + 1; i++ )
         krowM[i] += 1;
      for( int i = 0; i < nnz; i++ )
         jcolM[i] += 1;
      //cout << "Forming the matrix took:" << MPI_Wtime()-tm << endl;
   }


   phase = 12; // Analysis, numerical factorization
   nrhs = 1;

   iparm[30] = 0; // do not specify sparse rhs
                  // if one wants to use the sparse rhs and partial solves according to
                  // the PARDISO user guide the perm vector has to be present during all
                  // phases of pardiso

   assert(iparmUnchanged());

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, nullptr,
         &nrhs, iparm, &msglvl, nullptr, nullptr, &error);

   if( error != 0 )
   {
      printf("PardisoSolver - ERROR during factorization: %d\n", error);
      MPI_Abort(MPI_COMM_WORLD, -1);
   }
}

void PardisoSolver::solve( OoqpVector& rhs_in )
{
   SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
   double * sol_local = nvec;

   //int maxRefinSteps=(gLackOfAccuracy==0?3:6);

   /* same for mkl_pardiso and pardiso */
   phase = 33; // solve and iterative refinement
   nrhs = 1;
   iparm[30] = 0; // do not specify sparse rhs at this point

   assert(iparmUnchanged());

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, nullptr,
         &nrhs, iparm, &msglvl, rhs.elements(), sol_local, &error);

   if( error != 0 )
   {
      printf("PardisoSolver - ERROR during solve: %d", error);
      exit(1);
   }

   rhs.copyFromArray(sol_local);
}

void PardisoSolver::solve( GenMatrix& rhs_in )
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);

  int nrows,ncols; rhs.getSize(ncols,nrows);
  if(sz_sol<nrows*ncols) {
    sz_sol=nrows*ncols;
    if(sol) delete[] sol;
    sol = new double[sz_sol];
  }

  assert(nrows==n);


  phase = 33; // solve and iterative refinement
  nrhs = ncols;
  iparm[30] = 0;
  assert(iparmUnchanged());

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, nullptr,
         &nrhs, iparm, &msglvl, &rhs[0][0], sol, &error);

  if ( error != 0)
    printf ("PardisoSolver - ERROR during solve: %d", error ); 

  memcpy(&rhs[0][0], sol, sz_sol*sizeof(double));
}

void PardisoSolver::solve( GenMatrix& rhs_in, int *colSparsity)
{
  std::cout << "PardisoSolver - using sparse rhs but might lead to numerical troubles .. \n";
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);

  int nrows,ncols; rhs.getSize(ncols,nrows);
  if(sz_sol<nrows*ncols) {
    sz_sol=nrows*ncols;
    if(sol) delete[] sol;
    sol = new double[sz_sol];
  }
  assert(nrows==n);

  /* same for mkl_pardiso and pardiso */
   phase = 33; // solve and iterative refinement
   iparm[30] = 1;

   nrhs = ncols;

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM,
         colSparsity, &nrhs, iparm, &msglvl, &rhs[0][0], sol, &error);

  if ( error != 0)
    printf ("PardisoSolver - ERROR during solve: %d", error );

  memcpy(&rhs[0][0], sol, sz_sol*sizeof(double));
}

void PardisoSolver::solve( int nrhss, double* rhss, int* /*colSparsity*/ )
{
   assert(rhss);
   assert(nrhss >= 1);

   if( sz_sol < nrhss * n )
   {
      sz_sol = nrhss * n;
      delete[] sol;

      sol = new double[sz_sol];
   }

#ifndef NDEBUG
//   if( colSparsity )
//   {
//      for( int nr = 0; nr < nrhss; nr++ )
//      {
//         for( int i = 0; i < n; i++ )
//         {
//            const int rhspos = nr * n + i;
//            if( rhss[rhspos] != 0.0 )
//               assert(colSparsity[i] == 1);
//            else if( nrhss == 1 ) // does not work with zeroes in matrix, e.g. callback example
//               assert(colSparsity[i] == 0);
//         }
//      }
//   }
#endif

   /* same for mkl_pardiso and pardiso */
   phase = 33; // solve and iterative refinement
   int nrhss_local = nrhss;

   assert(iparmUnchanged());

// see notes on [30] earlier - cannot be specified on the go
// seems to detoriate parformance a lot - triggers many BiCGStab steps
//   if( colSparsity )
//   {
//      iparm[30] = 1; //sparse rhs
//   }
//   else
   {
      iparm[30] = 0;
   }

   assert(pt); assert(M); assert(krowM); assert(jcolM); assert(rhss); assert(sol);

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, nullptr,
         &nrhss_local, iparm, &msglvl, rhss, sol, &error);

   if( error != 0 )
   {
      printf("PardisoSolver - ERROR during solve: %d", error);
      exit(1);
   }

   memcpy(rhss, sol, n * nrhss * sizeof(double));
}

PardisoSolver::~PardisoSolver()
{
   delete[] jcolM;
   delete[] krowM;
   delete[] M;
   delete[] nvec;
   if(sol) delete[] sol;
}
