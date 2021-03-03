/*
 * PardisoMKLSchurSolver.C
 *
 *  Created on: 11.12.2020
 *      Author: Nils-Christian Kempke
 */
#include "PardisoMKLSchurSolver.h"

#include "pipschecks.h"
#include "SimpleVector.h"
#include "pipsdef.h"

#include "StochOptions.h"

#include "mkl_pardiso.h"
#include "mkl_types.h"

extern double g_iterNumber;

PardisoMKLSchurSolver::PardisoMKLSchurSolver( const SparseSymMatrix * sm ) :
   PardisoSchurSolver( sm )
{}


void PardisoMKLSchurSolver::setIparm(int* iparm) const
{
   /* common parameters */
   iparm[1] = 2; // 2 and 3 are for METIS - 2 is METIS 4.1, 3 is METIS 5.1, 0 for min degree ordering

   /* NOTE: if iparm[9] is less than 13 mkl_pardiso will not consistently produce the same schur complement as the other pardiso (on some examples)
    * this might not be an issue should be kept in mind though
    */
   iparm[30] = 0; // do not specify sparse rhs at this point ! MKL_PARDISO can either set iparm[35] or iparm[30]

   /* From INTEL (instead of iparm[2] which is not defined there):
    *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
    *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
    *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
    */
   iparm[7] = 0; // MKL_PARDISO runs into troubles otherwise
   iparm[9] = 13; // MKL_PARDISO need this in order to compute same schur decomposition as schenk pardiso
   iparm[10] = 0; // scaling for IPM KKT; used with IPARM(13)=1 or 2
   iparm[12] = 0; // improved accuracy for IPM KKT; used with IPARM(11)=1; use 2 for advanced matchings and higher accuracy.

   /* NOTE: requires iparm[23] = 1 which in return requires iparm[10] = iparm[12] = 0 */
   /* even though the documentation does not tell so setting iparm[23] = 10 is highly unstable and might result in segmentation faults */
   iparm[35] = -2; // compute the schur complement

   /* mkl_pardiso has no chkmatrix method - instead one can set iparm[26] */
#ifndef NDEBUG
   iparm[26] = 1;
#endif
}

void PardisoMKLSchurSolver::initPardiso()
{
   pardisoinit(pt, &mtype, iparm);
}

void PardisoMKLSchurSolver::solve( OoqpVector& rhs_in )
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
  /* pardiso from mkl does not support same functionality as pardiso-project
   *
   * pardiso project:
   * when computing the schur complement S with factorization matrices we will get
   *
   * [A11 A12]   [L11 0] [I 0] [U11 U12]
   * [A21 A22] = [L12 I] [0 S] [0     I]
   *
   * a subsequent solve call will then only solve for A11 x1 = b1 instead of the full
   * system.
   *
   * pardiso mkl:
   * while the schur complement is the same, the factorization computed, stored and
   * used for solve calls is a full factorization. thus pardiso from intel will always
   * solve the full system
   *
   * workaround is to solve
   *
   * (phase 331)
   * [L11   0] [z1] = [b1]
   * [L12   I] [z2] = [b2]
   *
   * (phase 332)
   * [I 0] [y1]   [z1]
   * [0 S] [y2] = [z2]
   *
   * (phase 333)
   * [U11 U12] [x1]   [y1]
   * [0     I] [x2] = [y2]
   *
   */

   // forward substitution
   double* z_n = new double[nvec_size];
   assert(iparm[7] == 0);
   assert(iparm[35] = -2);

   // this is necessary for usage of stage = 331/332/333
   iparm[9] = 0;

   // HACK: keeping iparm[35] = -2 will, for some reason, not compute the correct result
   // iparm[35] will be set to -2 after stage 333
   iparm[35] = 2;

   phase = 331;
   pardiso (pt , &maxfct , &mnum, &mtype, &phase,
         &n, eltsAug, rowptrAug, colidxAug,
         nullptr, &nrhs, iparm, &msglvl, rhs_n, z_n, &error);
   assert(error == 0);

   // diagonal substitution
   phase = 332;
   double* y_n = new double[nvec_size];

   pardiso (pt , &maxfct , &mnum, &mtype, &phase,
         &n, eltsAug, rowptrAug, colidxAug,
         nullptr, &nrhs, iparm, &msglvl, z_n, y_n, &error);
   assert(error == 0);

   // backward substitution
   for(int i = dim; i < nvec_size; ++i)
      y_n[i] = 0.0;
   phase = 333;

   pardiso (pt , &maxfct , &mnum, &mtype, &phase,
         &n, eltsAug, rowptrAug, colidxAug,
         nullptr, &nrhs, iparm, &msglvl, y_n, x_n, &error);
   assert(error == 0);

   iparm[35] = -2;

   delete[] z_n;
   delete[] y_n;


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


void PardisoMKLSchurSolver::computeSC(int nSCO,
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
   assert(iparm[35] == -2);
   assert(iparm[23] == 1);

   /* perm array has to be defined - it is of size of augmented system (dense) and inicates rows/colums we want to have in the schur complement with 1
    * rest set to 0
    * if specified like {1, 0, 0, 1, 0, 0} mkl_pardiso will first reorder row 0 and 3 to the bottom (hopefully stable sort) and then use the last two rows for computing
    * the schur complement
    * setting all entries to 1 will return the input matrix as schur complement
    */
   int perm[n] = {0};
   for(int i = n - nSC; i < n; ++i)
      perm[i] = 1;

   /* reordering and symbolic factorization */
   phase = 11;

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, eltsAug, rowptrAug,
               colidxAug, perm, &nrhs, iparm, &msglvl, nullptr, nullptr, &error);

   /* iparm[35] should now contain the number of non-zero entries in the Schur-complement */
   /* if it does not most likely the general iparm setup is faulty */
   assert(error == 0); assert(iparm[35] >= 0);

   /* preallocation of schur matrix arrays */
   int nnzSC = iparm[35];

   /* get deleted inside of destructor of schur_transposed */
   int* rowptrSCtransp = new int[nSC + 1];
   int* colidxSCtransp = new int[nnzSC];
   double* eltsSCtransp = new double[nnzSC];

   phase = 22;

   /* reset iparm[35] to -2 necessary! */
   iparm[35] = -2;
   int step = 1;

   // mkl pardiso returns arrays with zero-based index via export
   pardiso_export(pt, eltsSCtransp, rowptrSCtransp, colidxSCtransp, &step, iparm, &error);

   #ifdef TIMING_FLOPS
   HPM_Start("PARDISOFact");
   #endif

   /* factorization call */
   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, eltsAug, rowptrAug,
         colidxAug, perm, &nrhs, iparm, &msglvl, nullptr, nullptr, &error);

   #ifdef TIMING_FLOPS
   HPM_Stop("PARDISOFact");
   #endif

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

   /////////////////////////////////////////////////////
   // transpose the matrix since it is given in lower triangular form and we are using upper triangular form
   /////////////////////////////////////////////////////
   SparseStorage schur_transposed( nSC, nSC, nnzSC, rowptrSCtransp, colidxSCtransp, eltsSCtransp, true);

   rowptrSC = new int[nSC + 1];
   colidxSC = new int[nnzSC];
   eltsSC = new double[nnzSC];

   schur_transposed.transpose(rowptrSC, colidxSC, eltsSC);

   assert(subMatrixIsOrdered(rowptrSC, colidxSC, 0, nSC));
}


PardisoMKLSchurSolver::~PardisoMKLSchurSolver()
{
   int phase = -1; /* Release internal memory . */
   msglvl = 0;
   int error = 0;

   pardiso (pt, &maxfct, &mnum, &mtype, &phase,
       &n, nullptr, rowptrAug, colidxAug, nullptr, &nrhs,
       iparm, &msglvl, nullptr, nullptr, &error );
   if ( error != 0) {
     printf ("PardisoMKLSchurSolver - ERROR in pardiso release: %d", error );
   }

}
