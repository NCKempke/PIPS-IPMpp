/*
 * PardisoProjectSchurSolver.C
 *
 *  Created on: 10.12.2020
 *      Author: bzfkempk
 */

#include "PardisoProjectSchurSolver.h"

#include "pipsdef.h"
#include "pipschecks.h"

#include "SimpleVector.h"

extern "C" void pardisoinit(void*, int*, int*, int*, double*, int*);
extern "C" void pardiso(void*, int*, int*, int*, int*, int*, double*, int*, int*, int*, int*,
   int*, int*, double*, double*, int*, double*);

extern "C" void pardiso_chkmatrix(int*, int*, double*, int*, int*, int*);
extern "C" void pardiso_chkvec(int*, int*, double*, int*);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);
extern "C" void pardiso_get_schur(void*, int*, int*, int*, double*, int*, int*);

extern double g_iterNumber;

PardisoProjectSchurSolver::PardisoProjectSchurSolver( SparseSymMatrix * sgm )
   : PardisoSchurSolver( sgm )
{
   num_threads = PIPSgetnOMPthreads();

   initPardiso();
}

void PardisoProjectSchurSolver::solve( OoqpVector& rhs_in )
{
  SimpleVector& rhs=dynamic_cast<SimpleVector&>(rhs_in);

  int error = 0;
  assert(iparmUnchanged());

  assert(nvec_size == n);
  double* const x_n = nvec;
  double* const rhs_n = nvec2;

  const int dim=rhs.length();
  assert(dim >= 0 && dim <= n);

  memset(x_n, 0, dim * sizeof(double));
  memcpy(rhs_n, rhs.elements(), dim * sizeof(double));

  if( n > dim )
     memset(&rhs_n[dim], 0, (n - dim) * sizeof(double));

#ifdef TIMING_FLOPS
  HPM_Start("PARDISOSolve");
#endif

  // solving phase
  phase = 33; /* solve - iterative refinement */

  int* rhsSparsity = nullptr;
  if( useSparseRhs )
  {
     iparm[30] = 1; //sparse rhs
     rhsSparsity = new int[n]();

     for( int i = 0; i < dim; i++  )
        if( !PIPSisZero(rhs_n[i]) )
           rhsSparsity[i] = 1;
  }
  else
  {
     iparm[30] = 0;
  }

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, eltsAug, rowptrAug,
         colidxAug, rhsSparsity, &nrhs, iparm, &msglvl, rhs_n, x_n, &error, dparm);

  iparm[30] = 0;
  delete[] rhsSparsity;

  assert(error == 0);

#ifdef TIMING_FLOPS
  HPM_Stop("PARDISOSolve");
#endif


#ifdef PRINT_SC_RESIDUAL
  //compute residual (alternative)
  double* tmp_resid=new double[dim];
  memcpy(tmp_resid, rhs.elements(), dim*sizeof(double));
  double mat_max=0;
  for( int i = 0; i < dim; i++ )
  {
     for( int p = rowptrAug[i]; p < rowptrAug[i + 1]; p++ )
     {
        const int j = colidxAug[p - 1] - 1;
        if( j + 1 <= dim )
        {
           //r[i] = r[i] + M(i,j)*x(j)
           tmp_resid[i] -= eltsAug[p - 1] * x_n[j];

           if( abs(eltsAug[p - 1]) > mat_max )
              mat_max = abs(eltsAug[p - 1]);

           if( j != i )
           {
              //r[j] = r[j] + M(j,i)*x(i)
              tmp_resid[j] -= eltsAug[p - 1] * x_n[i];
           }
        }
     }
  }

  double res_norm2=0.0, res_nrmInf=0, sol_inf=0.;
  for( int i = 0; i < dim; i++ )
  {
     res_norm2 += tmp_resid[i] * tmp_resid[i];
     if( res_nrmInf < fabs(tmp_resid[i]) )
        res_nrmInf = tmp_resid[i];
     if( abs(x_n[i]) > sol_inf )
        sol_inf = abs(x_n[i]);
  }
  res_norm2 = sqrt(res_norm2);

  const double rhsNorm=rhs.twonorm();
  //if(min(res_nrmInf/(mat_max*sol_inf),res_norm2/(mat_max*sol_inf))>1e-3) {
  if(min(res_nrmInf/rhsNorm,res_norm2/rhsNorm)>1e-6) {
    cout << "PardisoSchurSolve large residual - norms resid="<< res_norm2 << ":" << res_nrmInf
    << " rhs=" << rhsNorm << " sol="<<sol_inf << " mat="<< mat_max
    << " #refin.=" << iparm[6]
    << " rel.res.nrmInf=" << res_nrmInf/rhsNorm
    << " bicgiter=" << gOuterBiCGIter<< endl;

  }
  delete[] tmp_resid;
#endif

  memcpy(&rhs[0], x_n, dim * sizeof(double));
}


void PardisoProjectSchurSolver::setIparm(int* iparm) const
{
   /* common parameters */
   iparm[1] = 2; // 2 and 3 are for METIS - 2 is METIS 4.1, 3 is METIS 5.1, 0 for min degree ordering

   /* NOTE: if iparm[9] is less than 13 mkl_pardiso will not consistently produce the same schur complement as the other pardiso (on some examples)
    * this might not be an issue should be kept in mind though
    */
   iparm[30] = 0; // do not specify sparse rhs at this point ! MKL_PARDISO can either set iparm[35] or iparm[30]

   iparm[2] = PIPSgetnOMPthreads();

   iparm[7] = nIterativeRefins; // max number of iterative refinement steps
   iparm[9] = pivotPerturbationExp; // pivot perturbation 10^{-x} * |A|_\{\inf}
   iparm[10] = 1; // default, scaling for IPM KKT used with either mtype=11/13 or mtype=-2/-4/6 and iparm[12]=1
   iparm[12] = 2;// 0 disable matching, 1 enable matching, no other settings

   if( factorizationTwoLevel )
      iparm[23] = 1; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
   else
      iparm[23] = 0;

   if( parallelForwardBackward )
      iparm[24] = 1; // parallelization for the forward and backward solve. 0=sequential, 1=parallel solve.
   else
      iparm[24] = 0;

}

void PardisoProjectSchurSolver::initPardiso()
{
   int error = 0;

   #pragma omp critical
   pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);
   if( error != 0 )
   {
      if( error == -10 )
         printf("PARDISO: No license file found \n");
      if( error == -11 )
         printf("PARDISO: License is expired \n");
      if( error == -12 )
         printf("PARDISO: Wrong username or hostname \n");

      PIPS_MPIabortIf(true, "Error in pardisoinit");
   }
}

void PardisoProjectSchurSolver::computeSC(int nSCO,
/*const*/SparseGenMatrix& R,
/*const*/SparseGenMatrix& A,
/*const*/SparseGenMatrix& C,
/*const*/SparseGenMatrix& F,
/*const*/SparseGenMatrix& G, int*& rowptrSC, int*& colidxSC, double*& eltsSC)
{
   assert(!rowptrSC && !colidxSC && !eltsSC);

   bool doSymbFact = false;
   if( firstSolve )
   {
      firstSolveCall(R, A, C, F, G, nSCO);
      firstSolve = false;
      doSymbFact = true;
   }
   else
   {
      //update diagonal entries in the PARDISO aug sys
      const double* eltsMsys = Msys->getStorageRef().M;
      std::map<int, int>::iterator it;

#if 0
      double max = -1e20;
      double min = 1e20;
      double minAbs = 1e20;

      for(it=diagMap.begin(); it!=diagMap.end(); it++)
      {
         const double elem = eltsMsys[it->first];
         if(elem > max)
         max = elem;
         if(elem < min)
         min = elem;
         if(std::fabs(elem) < minAbs && elem > 0.0 )
         minAbs = std::fabs(elem);

      }
      std::cout << "local Schur diag: min/max/minAbs  " << min << " " << max << " " << minAbs << std::endl;
#endif

      for( it = diagMap.begin(); it != diagMap.end(); it++ )
         eltsAug[it->second] = eltsMsys[it->first];
   }

   // call PARDISO
   int error = 0;

#ifndef NDEBUG
   pardiso_chkmatrix(&mtype, &n, eltsAug, rowptrAug, colidxAug, &error);
   if( error != 0 )
   {
      std::cout << "PARDISO matrix error " << error << "\n";
      exit(1);
   }
#endif

   const int nIter = (int) g_iterNumber;
   const int myRank = PIPS_MPIgetRank();

   if( (nIter % symbFactorInterval) == 0 )
   {
      doSymbFact = true;
      if( myRank == 0 )
         printf("PardisoSchur: starting symbolic analysis ... ");
   }

   int phase = 22; // numerical factorization
   if( doSymbFact )
      phase = 12; // numerical factorization & symb analysis

   assert(iparmUnchanged());

   /* compute schur complement */
   iparm[37] = nSC; // compute Schur-complement

   #ifdef TIMING_FLOPS
   HPM_Start("PARDISOFact");
   #endif

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, eltsAug, rowptrAug,
         colidxAug, nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error, dparm);

 #ifdef TIMING_FLOPS
   HPM_Stop("PARDISOFact");
 #endif

   if( doSymbFact && myRank == 0 )
      printf("finished \n");

   const int nnzSC = iparm[38];

#ifdef TIMING
   if(1001*(myRank/1001)==myRank)
   printf("rank %d perturbPiv %d peakmem %d\n", myRank, iparm[13], iparm[14]); // same for MKL and pardiso
   //cout << "NNZ(SCHUR) " << nnzSC << "    SPARSITY " << nnzSC/(1.0*nSC*nSC) << endl;
 #endif

   if( error != 0 )
   {
      printf("PardisoSolver - ERROR during factorization: %d. Phase param=%d\n",
            error, phase);
      exit(1);
   }
   rowptrSC = new int[nSC + 1];
   colidxSC = new int[nnzSC];
   eltsSC = new double[nnzSC];


   pardiso_get_schur(pt, &maxfct, &mnum, &mtype, eltsSC, rowptrSC, colidxSC);

   //convert back to C/C++ indexing
   for( int it = 0; it < nSC + 1; it++ )
      rowptrSC[it]--;
   for( int it = 0; it < nnzSC; it++ )
      colidxSC[it]--;

   assert(subMatrixIsOrdered(rowptrSC, colidxSC, 0, nSC));
}


PardisoProjectSchurSolver::~PardisoProjectSchurSolver()
{
   int phase = -1; /* Release internal memory . */
   msglvl = 0;
   int error = 0;

   pardiso (pt, &maxfct, &mnum, &mtype, &phase,
       &n, nullptr, rowptrAug, colidxAug, nullptr, &nrhs,
       iparm, &msglvl, nullptr, nullptr, &error , dparm );
   if ( error != 0) {
     printf ("PardisoProjectSchurSolver - ERROR in pardiso release: %d", error );
   }
}
