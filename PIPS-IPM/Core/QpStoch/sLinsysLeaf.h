/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHLEAFLINSYS
#define STOCHLEAFLINSYS

#include "sLinsys.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"

#include "pipsport.h"

#include "omp.h"

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to QpGenSparseLinsys.
 */
class sLinsysLeaf : public sLinsys
{
 public:
  //sLinsysLeaf(QpGenStoch * factory_, sData * prob_);
  template<class LINSOLVER>
    sLinsysLeaf(sFactory* factory,
		sData* prob_,				    
		OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
		OoqpVector* rhs_, OoqpVector* reg, OoqpVector* primal_reg,
      OoqpVector* dual_y_reg, OoqpVector* dual_z_reg,
      LINSOLVER *linsolver = nullptr);

  ~sLinsysLeaf() override = default;

  void factor2( sData *prob, Variables *vars) override;
  void Lsolve( sData*, OoqpVector& ) override {};
  void Dsolve( sData*, OoqpVector& x ) override;
  void Ltsolve( sData*, OoqpVector& ) override {};

  //virtual void Lsolve2 ( OoqpVector& x );
  //virtual void Dsolve2 ( OoqpVector& x );
  //virtual void solveCompressed( OoqpVector& rhs );

  virtual void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp);

  void putZDiagonal( const OoqpVector& zdiag ) override;
  void putXDiagonal( const OoqpVector& xdiag_ ) override;

  void regularize( const OoqpVector& primal_reg, const OoqpVector& dual_y_reg, const OoqpVector& dual_z_reg ) override;

  //void Ltsolve_internal(  sData *prob, StochVector& x, SimpleVector& xp);
  virtual void deleteChildren();

  using sLinsys::addTermToSchurComplBlocked;
  void addTermToSchurComplBlocked(sData *prob, bool sparseSC, SymMatrix& SC) override;

 protected:

  static void mySymAtPutSubmatrix(SymMatrix& kkt, GenMatrix& B, GenMatrix&, int locnx, int locmy, int);

  template<class LINSOLVER>
  void initBlockedSolvers();

  void addBorderTimesRhsToB0( StochVector& rhs, SimpleVector& b0, BorderLinsys& border ) override;

  void addBorderX0ToRhs( StochVector& rhs, const SimpleVector& x0, BorderLinsys& border ) override;

}; 

template<class LINSOLVER>
sLinsysLeaf::sLinsysLeaf(sFactory *factory_, sData* prob,
			 OoqpVector* dd_,
			 OoqpVector* dq_,
			 OoqpVector* nomegaInv_,
			 OoqpVector* rhs_,
			 OoqpVector* reg,
			 OoqpVector* primal_reg,
			 OoqpVector* dual_y_reg,
			 OoqpVector* dual_z_reg,
			 LINSOLVER* /*thesolver*/)
  : sLinsys(factory_, prob, dd_, dq_, reg, primal_reg, dual_y_reg, dual_z_reg, nomegaInv_, rhs_, false)
{
   static bool printed = false;
   const int n_omp_threads = PIPSgetnOMPthreads();
   if( pips_options::getIntParameter("LINEAR_LEAF_SOLVER") == SolverType::SOLVER_PARDISO )
   {
      n_solvers = std::max( 1, n_omp_threads / 2 );
      n_threads_solvers = ( n_omp_threads > 1 ) ? 2 : 1;
   }
   else
   {
      n_solvers = n_omp_threads;
      n_threads_solvers = 1;
   }

   if( computeBlockwiseSC )
   {
      if( PIPS_MPIgetRank() == 0 && !printed )
      {
         printed = true;
         std::cout << "Using " << n_solvers << " solvers in parallel (with "
            << n_threads_solvers << " threads each) for leaf SC computation - sLinsysLeaf\n";
      }
   }

#ifdef TIMING
  const int myRank = PIPS_MPIgetRank(mpiComm);
  const double t0 = MPI_Wtime();
#endif

  prob->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);
  const int n = locnx + locmy + locmz;

  int nnzQ, nnzB, nnzD;
  prob->getLocalNnz(nnzQ, nnzB, nnzD);

#ifdef TIMING
  if( myRank == 0 )
     std::cout << "Rank 0: building local Schur matrix ..." << std::endl;
#endif

  /* allocate and copy lower triangular:
   *
   * [ Qq BiT DiT ]
   * [ Bi  0   0  ]
   * [ Di  0   0  ]
   *
   * where  Qq = Q + V^-1 Gamma + W^-1 Phi (so we estimate its nnzs as n_1 + nnzQ)
   */

  SparseSymMatrix* kkt_sp = new SparseSymMatrix(n, n + nnzQ + nnzB + nnzD);

  SimpleVectorHandle v( new SimpleVector(n) );
  v->setToZero();
  kkt_sp->setToDiagonal(*v);

  kkt_sp->symAtPutSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx);

  if( locmz > 0 )
  {
     kkt_sp->symAtPutSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx);
     kkt_sp->symAtPutSubmatrix( locnx + locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx);
  }
  else
    mySymAtPutSubmatrix(*kkt_sp, prob->getLocalB(), prob->getLocalD(), locnx, locmy, locmz);

#ifdef TIMING
  if( myRank == 0 ) std::cout << "Rank 0: finished " << std::endl;
#endif

  assert( n_solvers >= 1 );

  solvers_blocked.resize(n_solvers);
  problems_blocked.resize(n_solvers);

  #pragma omp parallel num_threads(n_solvers)
  {
     omp_set_num_threads(n_threads_solvers);
     const int id = omp_get_thread_num();

     if( id == 0 )
        problems_blocked[id].reset( kkt_sp );
     else
        problems_blocked[id].reset( new SparseSymMatrix( *dynamic_cast<SparseSymMatrix*>(kkt_sp) ) );

     solvers_blocked[id].reset( new LINSOLVER( dynamic_cast<SparseSymMatrix*>(problems_blocked[id].get()), reg ) );
  }

  kkt = problems_blocked[0].get();
  solver = solvers_blocked[0].get();

#ifdef TIMING
  const double t1 = MPI_Wtime() - t0;
  if (myRank == 0) printf("Rank 0: new sLinsysLeaf took %f sec\n",t1);
#endif

  mpiComm = (dynamic_cast<StochVector*>(dd_))->mpiComm;
}

#endif
