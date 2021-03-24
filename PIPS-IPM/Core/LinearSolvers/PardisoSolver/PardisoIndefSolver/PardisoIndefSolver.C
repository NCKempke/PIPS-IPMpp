/*
 * PardisoIndefSolver.C
 *
 *  Created on: 21.03.2018
 *      Author: Daniel Rehfeldt
 */


#include "PardisoIndefSolver.h"

#include "pipschecks.h"
#include "SimpleVector.h"

#include <cmath>
#include <algorithm>
#include <cassert>

#include "omp.h"
#include "pipsport.h"
#include "pipsdef.h"
#include "StochOptions.h"


PardisoIndefSolver::PardisoIndefSolver( DenseSymMatrix * dm, bool solve_in_parallel, MPI_Comm mpi_comm ) :
      mpi_comm{mpi_comm}, solve_in_parallel{solve_in_parallel}
{
  mStorage = dm->getStorageHandle();
  mStorageSparse = nullptr;

  assert(mStorage);

  n = mStorage->n;

  initPardiso();
}

PardisoIndefSolver::PardisoIndefSolver( SparseSymMatrix * sm, bool solve_in_parallel, MPI_Comm mpi_comm ) :
      mpi_comm{mpi_comm}, solve_in_parallel{solve_in_parallel}
{
  mStorage = nullptr;
  mStorageSparse = sm->getStorageHandle();

  assert(mStorageSparse);

  n = mStorageSparse->n;

  initPardiso();
}

bool PardisoIndefSolver::iparmUnchanged()
{
   /* put all Parameters that should stay be checked against init into this array */
   static const int check_iparm[] =
      { 7, 10, 12, 23, 24 };

   bool unchanged = true;
   bool print = true;

   int iparm_compare[64];
   getIparm(iparm_compare);

   std::vector<int> to_compare(check_iparm,
         check_iparm + sizeof(check_iparm) / sizeof(check_iparm[0]));

   for( int i = 0; i < 64; ++i )
   {
      // if entry should be compared
      if( std::find(to_compare.begin(), to_compare.end(), i)
            != to_compare.end() )
      {
         if( iparm[i] != iparm_compare[i] )
         {
            if( print )
               std::cout
                     << "ERROR - PardisoIndefSolver: elements in iparm changed at "
                     << i << ": " << iparm[i] << " != " << iparm_compare[i]
                     << "(new)" << std::endl;
            unchanged = false;
         }
      }
   }
   return unchanged;
}

void PardisoIndefSolver::initPardiso()
{
   const int myRank = PIPS_MPIgetRank();

   iparm[0] = 0;

   pivotPerturbationExp = pips_options::getIntParameter("PARDISO_PIVOT_PERTURBATION_ROOT");
   if( pivotPerturbationExp < 0 )
	   pivotPerturbationExp = pivotPerturbationExpDefault;

   nIterativeRefins = pips_options::getIntParameter("PARDISO_NITERATIVE_REFINS_ROOT");
   if( nIterativeRefins < 0 )
 	  nIterativeRefins = nIterativeRefinsDefault;

   static bool printed = false;
   if( myRank == 0 && !printed)
   {
      printf("\nPARDISO solver settings\n\n");
      printf("PARDISO root: using pivot perturbation 10^-%d \n", pivotPerturbationExp);

      printf("PARDISO root: using maximum of %d iterative refinements  \n", nIterativeRefins);

      if( parallelForwardBackward )
         printf("PARDISO root: using parallel (forward/backward) solve \n");
      else
         printf("PARDISO root: NOT using parallel (forward/backward) solve \n");

      if( factorizationTwoLevel )
         printf("PARDISO root: using two-level scheduling for numerical factorization \n");
      else
         printf("PARDISO root: NOT using two-level scheduling for numerical factorization \n");

      if( highAccuracy )
         printf("PARDISO root: using high accuracy \n");
      else
         printf("PARDISO root: NOT using high accuracy \n");

      if( useSparseRhs )
         printf("PARDISO root: using sparse rhs \n");
      else
         printf("PARDISO root: NOT using sparse rhs \n");
      if( solve_in_parallel )
         printf("PARDISO root: allreducing SC and solving on every process \n");
      else
         printf("PARDISO root: only rank 0 does the SC solve \n");
      printf("\n");
   }
   printed = true;
}


void PardisoIndefSolver::matrixChanged()
{
   const int my_rank = PIPS_MPIgetRank(mpi_comm);

   if( solve_in_parallel || my_rank == 0 )
   {
      if( !pips_options::getBoolParameter("HIERARCHICAL") && my_rank == 0 )
         printf("\n PardisoIndefSolver: Schur complement factorization is starting ...\n ");

      if( mStorageSparse )
         factorizeFromSparse();
      else
         factorizeFromDense();

      if( !pips_options::getBoolParameter("HIERARCHICAL") && my_rank == 0 )
         printf("\n PardisoIndefSolver: Schur complement factorization completed \n ");
   }
}


void PardisoIndefSolver::matrixRebuild( DoubleMatrix& matrixNew )
{
   const int my_rank = PIPS_MPIgetRank(mpi_comm);
   if( solve_in_parallel || my_rank == 0 )
   {
      SparseSymMatrix& matrixNewSym = dynamic_cast<SparseSymMatrix&>(matrixNew);

      assert(matrixNewSym.getStorageRef().fortranIndexed());

      if( !pips_options::getBoolParameter("HIERARCHICAL") && my_rank == 0 )
         printf("\n Schur complement factorization is starting ...\n ");

      factorizeFromSparse(matrixNewSym);

      if( !pips_options::getBoolParameter("HIERARCHICAL") && my_rank == 0 )
         printf("\n Schur complement factorization completed \n ");
   }
}


void PardisoIndefSolver::factorizeFromSparse(SparseSymMatrix& matrix_fortran)
{
   assert(n == matrix_fortran.size());
   assert(!deleteCSRpointers);
   assert(matrix_fortran.getStorageRef().fortranIndexed());

   ia = matrix_fortran.krowM();
   ja = matrix_fortran.jcolM();
   a = matrix_fortran.M();

   assert(ia[0] == 1);

   // matrix initialized, now do the actual factorization
   factorize();
}


void PardisoIndefSolver::factorizeFromSparse()
{
   assert(mStorageSparse);

   const int nnz = mStorageSparse->len;
   const int* const iaStorage = mStorageSparse->krowM;
   const int* const jaStorage = mStorageSparse->jcolM;
   const double* const aStorage = mStorageSparse->M;
   const bool usePrecondSparse = pips_options::getBoolParameter("PRECONDITION_SPARSE");

   // first call?
   if( ia == nullptr )
   {
      assert(ja == nullptr && a == nullptr);
      deleteCSRpointers = true;

      ia = new int[n + 1];
      ja = new int[nnz];
      a = new double[nnz];
   }

   assert(n >= 0);

   std::vector<double>diag(n);

   // todo the sparse precond. stuff should be moved out and handled by Sparsifier class
   if( usePrecondSparse )
   {
	  const double t = precondDiagDomBound;

	  for( int r = 0; r < n; r++ )
	  {
		  const int j = iaStorage[r];
		  assert(jaStorage[j] == r);

		  diag[r] = fabs(aStorage[j]) * t;
	  }
   }

   ia[0] = 1;
   int nnznew = 0;

   for( int r = 0; r < n; r++ )
   {
      for( int j = iaStorage[r]; j < iaStorage[r + 1]; j++ )
      {
         if( aStorage[j] != 0.0 || jaStorage[j] == r )
         {
            if( usePrecondSparse )
            {
               if( (fabs(aStorage[j]) >= diag[r] || fabs(aStorage[j]) >= diag[jaStorage[j]]) )
               {
                  ja[nnznew] = jaStorage[j] + 1;
                  a[nnznew++] = aStorage[j];
               }
               else
               {
                  assert(jaStorage[j] != r);
               }
            }
            else
            {
               ja[nnznew] = jaStorage[j] + 1;
               a[nnznew++] = aStorage[j];
            }
         }
      }
      ia[r + 1] = nnznew + 1;
   }

   if( !pips_options::getBoolParameter("HIERARCHICAL") && PIPS_MPIgetRank(mpi_comm) == 0 )
      std::cout << "real nnz in KKT: " << nnznew << " (ratio: " << double(nnznew) / double(iaStorage[n]) << ")" << std::endl;

#if 0
   {
      ofstream myfile;
      int mype;  MPI_Comm_rank(mpi_comm, &mype);

      printf("\n\n ...WRITE OUT! \n\n");

      if( mype == 0 )
      {
         myfile.open("../A.txt");

         for( int i = 0; i < n; i++ )
            for( int k = ia[i]; k < ia[i + 1]; k++ )
               myfile << i << '\t' << ja[k - 1] << '\t' << a[k - 1] << endl;

         myfile.close();
      }

      printf("%d...exiting (pardiso) \n", mype);
      exit(1);
  }
#endif

   // matrix initialized, now do the actual factorization
   factorize();
}


void PardisoIndefSolver::factorizeFromDense()
{
   assert(mStorage);

#ifndef NDEBUG
  for( int i = 0; i < n; i++ )
  {
     for( int j = 0; j < n; j++ )
        assert(j <= i || PIPSisZero(mStorage->M[i][j]));
  }
#endif
#ifdef TIMING
  if( PIPS_MPIgetRank(mpi_comm) == 0 )
     std::cout << "from dense, starting factorization" << std::endl;
#endif

   assert(mStorage->n == mStorage->m);
   int nnz = 0;
   for( int i = 0; i < n; i++ )
      for( int j = 0; j <= i; j++ )
         if( mStorage->M[i][j] != 0.0 )
            nnz++;

   if( deleteCSRpointers )
   {
      delete[] ia;
      delete[] ja;
      delete[] a;
   }

   ia = new int[n + 1];
   ja = new int[nnz];
   a = new double[nnz];

   deleteCSRpointers = true;

   nnz = 0;
   for( int j = 0; j < n; j++ )
   {
      ia[j] = nnz;
      for( int i = j; i < n; i++ )
         if( mStorage->M[i][j] != 0.0 )
         {
            ja[nnz] = i;
            a[nnz++] = mStorage->M[i][j];
         }
   }

   ia[n] = nnz;

   for( int i = 0; i < n + 1; i++ )
      ia[i] += 1;
   for( int i = 0; i < nnz; i++ )
      ja[i] += 1;

   factorize();
}


void PardisoIndefSolver::factorize()
{
   int error;
   const int my_rank = PIPS_MPIgetRank(mpi_comm);

   assert(ia && ja && a);
   checkMatrix();

   iparm[17] = -1; /* compute number of nonzeros in factors */

#if 0
   const int nnz = ia[n] - 1;
   double abs_max = 0.0;
   for( int i = 0; i < nnz; i++ )
   {
      const double abs = std::fabs(a[i]);
      if( abs > abs_max)
         abs_max = abs;
   }

   std::cout << "absmax=" << abs_max << " log=" << log10(abs_max) << std::endl;
   if( log10(abs_max) >= 13)
   {
      iparm[9] = min(int(log10(abs_max)), 15);
      std::cout << "new: param " << iparm[9] << std::endl;
   }
else
#endif

   phase = 11;
   assert(iparmUnchanged());

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
         iparm, &msglvl, &ddum, &ddum, &error );

   if( error != 0 )
   {
      printf("\nERROR during symbolic factorization: %d", error);
      exit(1);
   }

   if( !pips_options::getBoolParameter("HIERARCHICAL") && my_rank == 0 )
   {
      printf("\nReordering completed: ");
      printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
   }

   phase = 22;
   assert( iparmUnchanged() );

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
         iparm, &msglvl, &ddum, &ddum, &error);

   if( error != 0 )
   {
      printf("\nERROR during numerical factorization: %d", error);
      exit(2);
   }
}

void PardisoIndefSolver::solveSynchronized( OoqpVector& vec )
{
   const bool solve_in_parallel_old = solve_in_parallel;
   solve_in_parallel = false;

   solve( vec );

   solve_in_parallel = solve_in_parallel_old;
}


void PardisoIndefSolver::solve( OoqpVector& v )
{
   assert( iparmUnchanged() );

   const int size = PIPS_MPIgetSize(mpi_comm);
   const int my_rank = PIPS_MPIgetRank(mpi_comm);

   phase = 33;
   SimpleVector& sv = dynamic_cast<SimpleVector&>(v);

   double* b = sv.elements();

   assert(sv.length() == n);

#ifdef TIMING_FLOPS
   HPM_Start("DSYTRSSolve");
#endif
   if( solve_in_parallel || my_rank == 0 )
   {
      int* rhsSparsity = nullptr;

      int error;

      // first call?
      if( !x )
         x = new double[n];

      if( useSparseRhs )
      {
         iparm[30] = 1; //sparse rhs
         rhsSparsity = new int[n]();

         for( int i = 0; i < n; i++  )
            if( !PIPSisZero(b[i]) )
               rhsSparsity[i] = 1;
      }
      else
      {
         iparm[30] = 0;
      }

      pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, rhsSparsity, &nrhs,
            iparm, &msglvl, b, x, &error);
      if( error != 0 )
      {
         printf("\nERROR during solution: %d", error);
         exit(3);
      }

      iparm[30] = 0;
#if 0
      const double b2norm = sv.twonorm();
      const double binfnorm = sv.infnorm();
      double mat_max = 0.0;
      for( int i = 0; i < n; i++ )
      {
         for( int p = ia[i]; p < ia[i + 1]; p++ )
         {
            const int j = ja[p - 1] - 1;

            sv[i] -= a[p - 1] * x[j]; //r[i] = r[i] - M(i,j)*x(j)

            mat_max = std::max(std::fabs(a[p - 1]), mat_max);

            assert(j >= i);

            if( j != i )
               sv[j] -= a[p - 1] * x[i]; //r[j] = r[j] - M(j,i)*x(i)
         }
      }
      const double res2norm = sv.twonorm();
      const double resinfnorm = sv.infnorm();

      if( my_rank == 0 )
         std::cout << "GLOBAL SCHUR: res.2norm=" << res2norm << " rel.res2norm=" << res2norm / b2norm  <<
            " res.infnorm=" << resinfnorm << " rel.resinfnorm=" << resinfnorm / binfnorm  << " b2norm=" << b2norm <<
            " abs.mat.elem.=" << mat_max << std::endl;
#endif

      for( int i = 0; i < n; i++ )
         b[i] = x[i];

      if( size > 0 && !solve_in_parallel )
         MPI_Bcast(b, n, MPI_DOUBLE, 0, mpi_comm);

      delete[] rhsSparsity;

#ifdef TIMING
      printf("sparse kkt iterative refinement steps: %d \n", iparm[6]);
#endif
   }
   else
   {
      assert( !solve_in_parallel );
      assert(size > 0);
      MPI_Bcast(b, n, MPI_DOUBLE, 0, mpi_comm);
   }


#ifdef TIMING_FLOPS
   HPM_Stop("DSYTRSSolve");
#endif
}

void PardisoIndefSolver::solve( int nrhss, double* rhss, int* /*colSparsity*/ )
{
   assert(rhss);
   assert(nrhss >= 1);

   if( static_cast<int>(sol.size()) < nrhss * n )
      sol.resize( nrhss * n );

   if( static_cast<int>(rhss_nonzero.size()) < nrhss * n )
   {
      rhss_nonzero.resize( nrhss * n );
      map_rhs_nonzero_original.resize( nrhss );
   }

   int nrhss_nnz{0};
   for( int i = 0; i < nrhss; ++i )
   {
      if( std::any_of( rhss + i * n, rhss + (i + 1) * n, [](double d){ return !PIPSisZero(d); }) )
      {
         std::copy( rhss + i * n, rhss + (i + 1) * n, rhss_nonzero.begin() + nrhss_nnz * n);
         map_rhs_nonzero_original[nrhss_nnz] = i;
         ++nrhss_nnz;
      }
   }

#ifndef NDEBUG
//   if( colSparsity )
//   {]
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
   // PARDISO cannot deal well with all zero rhs
   for( int i = 0; i < nrhss_nnz; ++i )
      assert( std::any_of(rhss_nonzero.begin() + i * n, rhss_nonzero.begin() + (i + 1) * n, [](double d){
               return d != 0.0;
         })
      );
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

   assert(pt); assert(a); assert(ia); assert(ja); assert(rhss);
   int error = 0;

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, nullptr,
         &nrhss_local, iparm, &msglvl, rhss_nonzero.data(), sol.data(), &error);

   if( error != 0 )
   {
      printf("PardisoSolver - ERROR during solve: %d", error);
      exit(1);
   }

   for( int i = 0; i < nrhss_nnz; ++i )
      std::copy(sol.begin() + i * n, sol.begin() + (i + 1) * n, rhss + map_rhs_nonzero_original[i] * n );
}

void PardisoIndefSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
   matrixChanged();
}

PardisoIndefSolver::~PardisoIndefSolver()
{
   if( deleteCSRpointers )
   {
      delete[] ia;
      delete[] ja;
      delete[] a;
   }

   delete[] x;
}
