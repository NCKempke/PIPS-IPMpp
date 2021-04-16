/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysRootAug.h"

#include "DeSymIndefSolver.h"
#include "DeSymIndefSolver2.h"
#include "DeSymPSDSolver.h"

#include "DistributedQP.hpp"
#include "BorderedSymMatrix.h"


#include "pipsport.h"
#include "StochOptions.h"
#include <limits>
#include <unistd.h>
#include "math.h"

//#define DUMPKKT
#ifdef DUMPKKT
#include <iostream>
#include <fstream>
#endif
#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif
extern int gInnerBiCGIter;
extern int gInnerBiCGFails;



static void biCGStabCommunicateStatus(int flag, int it)
{
   gInnerBiCGIter = it;

   if( flag != 0 )
      gInnerBiCGFails++;
}

sLinsysRootAug::sLinsysRootAug(sFactory * factory_, DistributedQP * prob_)
  : sLinsysRoot(factory_, prob_)
{ 
   if( pips_options::getBoolParameter( "HIERARCHICAL" ) )
      assert( false && "should not end up here");

   assert(locmyl >= 0 && locmzl >= 0);

   createSolversAndKKts(data);

   redRhs.reset( new SimpleVector(locnx + locmy + locmz + locmyl + locmzl) );
}

sLinsysRootAug::sLinsysRootAug(sFactory* factory_,
			       DistributedQP* prob_,
			       OoqpVector* dd_, 
			       OoqpVector* dq_,
			       OoqpVector* nomegaInv_,
			       OoqpVector* rhs_, bool create_solvers)
  : sLinsysRoot(factory_, prob_, dd_, dq_, nomegaInv_, rhs_)
{ 
   assert( pips_options::getBoolParameter( "HIERARCHICAL" ) );
   assert(locmyl >= 0 && locmzl >= 0);
   assert( computeBlockwiseSC );


   if( create_solvers )
      createSolversAndKKts(data);

   redRhs.reset( new SimpleVector(locnx + locmy + locmz + locmyl + locmzl) );
}

SymMatrix* sLinsysRootAug::createKKT(DistributedQP* prob) const
{
   const int n = locnx + locmy + locmyl + locmzl;
   if( hasSparseKkt )
   {
      SparseSymMatrix* sparsekkt;

      if( usePrecondDist )
         sparsekkt = prob->createSchurCompSymbSparseUpperDist(childrenProperStart, childrenProperEnd);
      else
         sparsekkt = prob->createSchurCompSymbSparseUpper();

      assert(sparsekkt->size() == n);

      return sparsekkt;
   }
   else
   {
      return new DenseSymMatrix(n);
   }
}

void sLinsysRootAug::createSolversSparse(SolverType solver_type)
{
   SparseSymMatrix* kkt_sp = dynamic_cast<SparseSymMatrix*>(kkt.get());

   if( solver_type == SolverType::SOLVER_MUMPS )
   {
#ifdef WITH_MUMPS
   solver.reset( new MumpsSolverRoot(mpiComm, kkt_sp, allreduce_kkt) );
#endif
   }
   else if( solver_type == SolverType::SOLVER_PARDISO )
   {
#ifdef WITH_PARDISO
      solver.reset( new PardisoProjectIndefSolver(kkt_sp, allreduce_kkt, mpiComm) );
#endif
   }
   else if( solver_type == SolverType::SOLVER_MKL_PARDISO )
   {
#ifdef WITH_MKL_PARDISO
      solver.reset( new PardisoMKLIndefSolver(kkt_sp, allreduce_kkt, mpiComm) );
#endif
   }
   else if( solver_type == SolverType::SOLVER_MA57 )
   {
#ifdef WITH_MA57
      solver.reset( new Ma57SolverRoot(kkt_sp, allreduce_kkt, mpiComm, "sLinsysRootAug") );
#endif
   }
   else
   {
      assert( solver_type == SolverType::SOLVER_MA27 );
#ifdef WITH_MA27
      solver.reset( new Ma27SolverRoot(kkt_sp, allreduce_kkt, mpiComm, "sLinsysRootAug") );
#endif
   }
}

void sLinsysRootAug::createSolversDense()
{
   const SolverTypeDense solver_type = pips_options::getSolverDense();
   DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kkt.get());

   if( solver_type == SolverTypeDense::SOLVER_DENSE_SYM_INDEF )
      solver.reset( new DeSymIndefSolver(kktmat) );
   else if( solver_type == SolverTypeDense::SOLVER_DENSE_SYM_INDEF_SADDLE_POINT )
      solver.reset( new DeSymIndefSolver2(kktmat, locnx) );
   else
   {
      assert( solver_type == SolverTypeDense::SOLVER_DENSE_SYM_PSD );
      solver.reset( new DeSymPSDSolver(kktmat) );
   }
}

void sLinsysRootAug::createSolversAndKKts(DistributedQP* prob)
{
   const SolverType solver_root = pips_options::getSolverRoot();

   static bool printed = false;
   if( !printed && PIPS_MPIgetRank() == 0 )
   {
      if( hasSparseKkt )
         std::cout << "sLinsysRootAug: using " << solver_root << "\n";
      else
         std::cout << "sLinsysRootAug: using " << pips_options::getSolverDense() << "\n";
   }

   kkt.reset(createKKT(prob));

   if( !printed && PIPS_MPIgetRank() == 0 )
   {
      if( hasSparseKkt )
         std::cout << "sLinsysRootAug: getSchurCompMaxNnz " << prob->getSchurCompMaxNnz() << "\n";
      else
      {
         const int n = locnx + locmy + locmyl + locmzl;
         std::cout << "sLinsysRootAug: getSchurCompMaxNnz " << n * n << "\n";
      }
   }
   printed = true;


   if( hasSparseKkt )
      createSolversSparse(solver_root);
   else
      createSolversDense();
}

#ifdef TIMING
static double t_start, troot_total, taux, tchild_total, tcomm_total;
#endif


void sLinsysRootAug::finalizeKKT(DistributedQP* prob, Variables* vars)
{
  stochNode->resMon.recFactTmLocal_start();
  stochNode->resMon.recSchurMultLocal_start();

  if( usePrecondDist )
  {
     // don't do anything, already done previously
  }
  else
  {
     if( hasSparseKkt )
        finalizeKKTsparse(prob, vars);
     else
        finalizeKKTdense(prob, vars);
  }

  stochNode->resMon.recSchurMultLocal_stop();
  stochNode->resMon.recFactTmLocal_stop();
}


void sLinsysRootAug::finalizeKKTdist(DistributedQP* prob)
{
   assert(kkt && hasSparseKkt && prob);

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

   int myRank; MPI_Comm_rank(mpiComm, &myRank);
   int mpiCommSize; MPI_Comm_size(mpiComm, &mpiCommSize);
   const bool iAmLastRank = (myRank == mpiCommSize - 1);
   const int childStart = childrenProperStart;
   const int childEnd = childrenProperEnd;

#ifndef NDEBUG
   int* const jcolKkt = kkts.jcolM();
#endif
   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();

   assert(childStart >= 0 && childStart < childEnd);
   assert(kkts.size() == locnx + locmy + locmyl + locmzl);
   assert(!kkts.isLower);
   assert(locmyl >= 0 && locmzl >= 0);
   assert(prob->getLocalQ().krowM()[locnx] == 0 && "Q currently not supported for dist. sparse kkt");

   /////////////////////////////////////////////////////////////
   // update the KKT with the diagonals
   // xDiag is in fact diag(Q)+X^{-1}S
   /////////////////////////////////////////////////////////////

   if( xDiag && iAmLastRank )
   {
      const SimpleVector& sxDiag = dynamic_cast<const SimpleVector&>(*xDiag);

      for( int i = 0; i < locnx; i++ )
      {
         const int diagIdx = krowKkt[i];
         assert(jcolKkt[diagIdx] == i);

         MKkt[diagIdx] += sxDiag[i];
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with   - C' * diag(zDiag) *C
   /////////////////////////////////////////////////////////////
   if( locmz > 0 && iAmLastRank )
   {
      assert(zDiag);

      SparseGenMatrix& C = prob->getLocalD();
      C.matTransDinvMultMat(*zDiag, &CtDC);
      assert(CtDC->size() == locnx);

      //aliases for internal buffers of CtDC
      SparseSymMatrix* CtDCsp = dynamic_cast<SparseSymMatrix*>(CtDC);
      const int* krowCtDC = CtDCsp->krowM();
      const int* jcolCtDC = CtDCsp->jcolM();
      const double* dCtDC = CtDCsp->M();

      for( int i = 0; i < locnx; i++ )
      {
         const int pend = krowCtDC[i + 1];
         for( int p = krowCtDC[i]; p < pend; p++ )
         {
            const int col = jcolCtDC[p];

            if( col >= i )
            {
               // get start position of dense kkt block
               const int blockStart = krowKkt[i];
               assert(col < locnx && jcolKkt[blockStart + col - i] == col);

               MKkt[blockStart + col - i] -= dCtDC[p];
            }
         }
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with At (symmetric update forced)
   /////////////////////////////////////////////////////////////
   if( locmy > 0 && iAmLastRank )
   {
      SparseGenMatrix& At = prob->getLocalB().getTranspose(); // yes, B
      const double* MAt = At.M();
      const int* krowAt = At.krowM();

      for( int i = 0; i < locnx; ++i )
      {
         const int pstart = krowAt[i];
         const int pend = krowAt[i + 1];

         // get start position of sparse kkt block
         const int blockStart = krowKkt[i] + locnx - i;

         assert(blockStart <= krowKkt[i + 1]);

         for( int p = pstart; p < pend; ++p )
         {
            assert(At.jcolM()[p] < locmy);
            assert(blockStart + (p - pstart) <= krowKkt[i + 1]);
            assert(jcolKkt[blockStart + (p - pstart)] == (locnx + At.jcolM()[p]));

            MKkt[blockStart + (p - pstart)] += MAt[p];
         }
      }
   }

   int local2linksStartEq;
   int local2linksEndEq;
   int local2linksStartIneq;
   int local2linksEndIneq;

   prob->getSCrangeMarkersMy(childStart, childEnd, local2linksStartEq, local2linksEndEq,
         local2linksStartIneq, local2linksEndIneq);

   const int n2linksRowsLocalEq = local2linksEndEq - local2linksStartEq;

   PIPSdebugMessage("rank %d FT local columns: %d-%d \n", myRank, local2linksStartEq, local2linksEndEq);
   PIPSdebugMessage("rank %d GT local columns: %d-%d \n", myRank, local2linksStartIneq, local2linksEndIneq);
   PIPSdebugMessage("rank %d FT local columns: %d-%d \n", myRank, local2linksStartEq, local2linksEndEq);

   /////////////////////////////////////////////////////////////
   // update the KKT with Ft
   /////////////////////////////////////////////////////////////
   if( locmyl > 0 )
   {
      SparseGenMatrix& Ft = prob->getLocalF().getTranspose();

      // add locally owned sparse part of Ft
      addLinkConsBlock0Matrix(prob, Ft, locnx + locmy, 0, local2linksStartEq, local2linksEndEq);

      if( myRank == 0 )
      {
         const int n2linksRowsEq = prob->n2linkRowsEq();
         const int bordersizeEq = locmyl - n2linksRowsEq;
         const int borderstartEq = locnx + locmy + n2linksRowsEq;

         PIPSdebugMessage("rank %d FT border columns: %d-%d \n", myRank, borderstartEq, borderstartEq + bordersizeEq);

         // add (shared) border part of Ft
         addLinkConsBlock0Matrix(prob, Ft, locnx + locmy, n2linksRowsLocalEq, borderstartEq, borderstartEq + bordersizeEq);
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with Gt and add z diagonal
   /////////////////////////////////////////////////////////////
   if( locmzl > 0 )
   {
      SparseGenMatrix& Gt = prob->getLocalG().getTranspose();
      const int n2linksRowsIneq = prob->n2linkRowsIneq();
      const int bordersizeIneq = locmzl - n2linksRowsIneq;
      const int borderstartIneq = locnx + locmy + locmyl + n2linksRowsIneq;

      // add locally owned sparse part of Gt
      addLinkConsBlock0Matrix(prob, Gt, locnx + locmy + locmyl, n2linksRowsLocalEq, local2linksStartIneq, local2linksEndIneq);

      if( myRank == 0 )
      {
         const int n2linksRowsLocalIneq = local2linksEndIneq - local2linksStartIneq;
         PIPSdebugMessage("rank %d GT border columns: %d-%d\n", myRank, borderstartIneq, borderstartIneq + bordersizeIneq);

         // add (shared) border part of Gt
         addLinkConsBlock0Matrix(prob, Gt, locnx + locmy + locmyl, n2linksRowsLocalEq + n2linksRowsLocalIneq,
               borderstartIneq, borderstartIneq + bordersizeIneq);
      }

      assert(zDiagLinkCons);
      const SimpleVector& szDiagLinkCons = dynamic_cast<const SimpleVector&>(*zDiagLinkCons);

      assert(local2linksStartIneq >= locnx + locmy + locmyl);
      assert(local2linksEndIneq <= locnx + locmy + locmyl + locmzl);

      const int szDiagLocalStart = local2linksStartIneq - (locnx + locmy + locmyl);
      assert(szDiagLocalStart >= 0);
      assert(szDiagLocalStart < locmzl || (szDiagLocalStart == locmzl && local2linksStartIneq == local2linksEndIneq));

      // add locally owned part of z diagonal
      for( int i = szDiagLocalStart, iKkt = local2linksStartIneq; iKkt < local2linksEndIneq; ++i, ++iKkt )
      {
         const int idx = krowKkt[iKkt];
         assert(jcolKkt[idx] == iKkt);
         assert(i < locmzl);

         MKkt[idx] += szDiagLinkCons[i];
      }

      if( myRank == 0 )
      {
         const int szDiagBorderStart = borderstartIneq - (locnx + locmy + locmyl);

         assert(szDiagBorderStart >= 0 && szDiagBorderStart <= locmzl);
         assert(szDiagBorderStart + bordersizeIneq == locmzl);

         // add border part of diagonal
         for( int i = szDiagBorderStart, iKkt = borderstartIneq; iKkt < borderstartIneq + bordersizeIneq; ++i, ++iKkt )
         {
            const int idx = krowKkt[iKkt];
            assert(jcolKkt[idx] == iKkt);
            assert(i < locmzl);

            MKkt[idx] += szDiagLinkCons[i];
         }
      }
   }
}

void sLinsysRootAug::assembleLocalKKT( DistributedQP* prob )
{
   const bool is_layer_only_twolinks = prob->isHierarchySparseTopLayerOnlyTwolinks();
   if( !pips_options::getBoolParameter("HIERARCHICAL") )
      assert( !is_layer_only_twolinks );

   for( size_t c = 0; c < children.size(); ++c )
   {
#ifdef STOCH_TESTING
      g_scenNum = c;
#endif
      if( children[c]->mpiComm == MPI_COMM_NULL )
         continue;

      children[c]->stochNode->resMon.recFactTmChildren_start();

      //---------------------------------------------
      addTermToSchurCompl( prob, c, !is_layer_only_twolinks );

      //---------------------------------------------
      children[c]->stochNode->resMon.recFactTmChildren_stop();
   }
}

/* compute Schur rhs b0 - sum Bi^T Ki^-1 bi for all children */
void sLinsysRootAug::Lsolve(DistributedQP *prob, OoqpVector& x)
{
   assert( !is_hierarchy_root );

   StochVector& b = dynamic_cast<StochVector&>(x);
   assert(children.size() == b.children.size() );

   SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.first);
   assert(!b.last);

   if( iAmDistrib && PIPS_MPIgetRank(mpiComm) > 0 )
      b0.setToZero();

   // compute Bi^T Ki^-1 rhs_i and sum it up
   for( size_t it = 0; it < children.size(); it++ )
   {
#ifdef TIMING
      children[it]->stochNode->resMon.eLsolve.clear();
      children[it]->stochNode->resMon.recLsolveTmChildren_start();
#endif
      children[it]->addLniziLinkCons( prob->children[it], b0, *b.children[it], true );

#ifdef TIMING
      children[it]->stochNode->resMon.recLsolveTmChildren_stop();
#endif
   }
#ifdef TIMING
   MPI_Barrier(MPI_COMM_WORLD);
   stochNode->resMon.eReduce.clear();//reset
   stochNode->resMon.recReduceTmLocal_start();
#endif
   if(iAmDistrib)
      PIPS_MPIsumArrayInPlace( b0.elements(), b0.length(), mpiComm );

#ifdef TIMING
   stochNode->resMon.recReduceTmLocal_stop();
#endif
  //dumpRhs(0, "rhs",  b0);
}

/* does Schur Complement solve */
void sLinsysRootAug::Dsolve( DistributedQP *prob, OoqpVector& x )
{
  /* Ki^-1 bi has already been computed in Lsolve */

  /* children have already computed Li^T\Di\Li\bi in Lsolve() */
  StochVector& b = dynamic_cast<StochVector&>(x);
  SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.first);
#ifdef TIMING
  stochNode->resMon.eDsolve.clear();
  stochNode->resMon.recDsolveTmLocal_start();
#endif
  solveReducedLinkCons(prob, b0);
#ifdef TIMING
  stochNode->resMon.recDsolveTmLocal_stop();
#endif
}

void sLinsysRootAug::Ltsolve( DistributedQP *prob, OoqpVector& x )
{
   StochVector& b = dynamic_cast<StochVector&>(x);
   SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.first);

   //dumpRhs(0, "sol",  b0);
   SimpleVector& z0 = b0; //just another name, for clarity

   for(size_t it = 0; it < children.size(); it++)
      children[it]->Ltsolve2(prob->children[it], *b.children[it], z0, true);

#ifdef TIMING
   int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if( 256 * ( myRank / 256 ) == myRank )
   {
      double tTotResChildren=0.0;
      for( size_t it = 0; it < children.size(); it++)
      {
         if( children[it]->mpiComm == MPI_COMM_NULL )
            continue;
         tTotResChildren += children[it]->stochNode->resMon.eLsolve.tmChildren;
         tTotResChildren += children[it]->stochNode->resMon.eLsolve.tmLocal;
      }
      double tComm = stochNode->resMon.eReduce.tmLocal;

      //double tTotChildren = 0.0;
      //for( size_t it = 0; it < children.size(); it++)
      //{
      //   tTotChildren += children[it]->stochNode->resMon.eDsolve.tmChildren;
      //   tTotChildren += children[it]->stochNode->resMon.eDsolve.tmLocal;
      //}
      double tStg1 = stochNode->resMon.eDsolve.tmLocal;

      double tTotStg2Children = 0.0;
      for( size_t it = 0; it < children.size(); it++)
      {
         if( children[it]->mpiComm == MPI_COMM_NULL )
            continue;
         tTotStg2Children += children[it]->stochNode->resMon.eLtsolve.tmChildren;
         tTotStg2Children += children[it]->stochNode->resMon.eLtsolve.tmLocal;
      }
      std::cout << "  rank " << myRank << " " << "Resid comp " << tTotResChildren << " " << "reduce " << tComm << " "
            << "1stStage solve " << tStg1 << " " << "2ndStage solve " << tTotStg2Children << "\n";
   }
#endif
}

/* gets called for computing the dense schur complement*/
void sLinsysRootAug::LsolveHierarchyBorder( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool two_link_border, int begin_cols, int end_cols )
{
   LsolveHierarchyBorder( result, Br, Br_mod_border, true, two_link_border, begin_cols, end_cols );
}

void sLinsysRootAug::LtsolveHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br,
      std::vector<BorderMod>& br_mod_border, bool sym_res, bool sparse_res, int begin_cols, int end_cols )
{
   if( Bl.isEmpty() || (Br.isEmpty() && br_mod_border.empty()) )
      return;

   LtsolveHierarchyBorder( res, X0, Bl, Br, br_mod_border, sym_res, sparse_res, true, begin_cols, end_cols );
}

extern int gLackOfAccuracy;

void sLinsysRootAug::solveReducedLinkCons( DistributedQP*, SimpleVector& b_vec)
{
#ifdef TIMING
   t_start = MPI_Wtime();
   troot_total = tchild_total = tcomm_total = 0.0;
#endif

   assert(locmyl >= 0 && locmzl >= 0);
   assert(locnx + locmy + locmz + locmyl + locmzl == b_vec.length());

   ///////////////////////////////////////////////////////////////////////
   // LOCAL SOLVE WITH SCHUR COMPLEMENT and i-th rhs from buffer
   ///////////////////////////////////////////////////////////////////////
   SparseGenMatrix &C = data->getLocalD();

   double *rhs_reduced = redRhs->elements();
   assert( redRhs->length() >= b_vec.length() );

   ///////////////////////////////////////////////////////////////////////
   // b = [b1; b2; b3; b4; b5] is a locnx + locmy + locmz + locmyl + locmz vector
   // the new rhs should be
   //           r = [b1-C^T*(zDiag)^{-1}*b3; b2; b4; b5]
   ///////////////////////////////////////////////////////////////////////
   double *b = b_vec.elements();

   //copy all elements from b into r except for the the residual values corresponding to z0 = b3
   // copy b1, b2
   std::copy(b, b + locnx + locmy, rhs_reduced);
   // copy b4, b5
   std::copy(b + locnx + locmy + locmz, b + locnx + locmy + locmz + locmyl + locmzl, rhs_reduced + locnx + locmy);
   // copy b3 to the end - used as buffer for reduction computations
   std::copy(b + locnx + locmy, b + locnx + locmy + locmz, rhs_reduced + locnx + locmy + locmyl + locmzl);
   // rhs_reduced now : [ b1; b2; b4; b5; b3]

   // alias to r1 part (no mem allocations)
   SimpleVector rhs1(&rhs_reduced[0], locnx);
   SimpleVector b3(rhs_reduced + locnx + locmy + locmyl + locmzl, locmz);

   ///////////////////////////////////////////////////////////////////////
   // compute r1 = b1 - C^T * (zDiag)^{-1} * b3
   ///////////////////////////////////////////////////////////////////////
   // if we have C part
   if( locmz > 0 )
   {
      assert(b3.length() == zDiag->length());
      b3.componentDiv(*zDiag);
      C.transMult(1.0, rhs1, -1.0, b3);
   }

   ///////////////////////////////////////////////////////////////////////
   // rhs_reduced now contains all components -> solve for it
   ///////////////////////////////////////////////////////////////////////

   // we do not need the last locmz elements of r since they were buffer only
   SimpleVector rhs_short(rhs_reduced, locnx + locmy + locmyl + locmzl);

   if( innerSCSolve == 0 )
   {
      // Option 1. - solve with the factors
      solver->solveSynchronized(rhs_short);
   }
   else if( innerSCSolve == 1 )
   {
      // Option 2 - solve with the factors and perform iter. ref.
      solveWithIterRef(data, rhs_short);
   }
   else
   {
      assert(innerSCSolve == 2);
      // Option 3 - use the factors as preconditioner and apply BiCGStab
      solveWithBiCGStab(data, rhs_short);
   }

   ///////////////////////////////////////////////////////////////////////
   // rhs_small is now the solution to the reduced system
   // the solution to the augmented system can now be computed as
   //      x = [rhs1; rhs2; zDiag^{-1} * (b3 - C * r1); rhs3; rhs4]
   ///////////////////////////////////////////////////////////////////////

   // copy the solution components and calculate r3
   // copy rhs1 and rhs2
   std::copy(rhs_reduced, rhs_reduced + locnx + locmy, b);
   // compute x3
   if( locmz > 0 )
   {
      SimpleVector rhs1(rhs_reduced, locnx);
      SimpleVector b3(b + locnx + locmy, locmz);
      C.mult(1.0, b3, -1.0, rhs1);
      b3.componentDiv(*zDiag);
   }

   // copy rhs3 and rhs4
   std::copy(rhs_reduced + locnx + locmy,
         rhs_reduced + locnx + locmy + locmyl + locmzl,
         b + locnx + locmy + locmz);

#ifdef TIMING
  if( myRank == 0 && innerSCSolve >= 1 )
    std::cout << "Root - Refin times: child=" << tchild_total << " root=" << troot_total
       << " comm=" << tcomm_total << " total=" << MPI_Wtime()-t_start << "\n";
#endif
}

void sLinsysRootAug::solveReducedLinkConsBlocked( DistributedQP* data, DenseGenMatrix& rhs_mat_transp, int rhs_start, int n_rhs )
{
#ifdef TIMING
   t_start = MPI_Wtime();
   troot_total = tchild_total = tcomm_total = 0.0;
#endif

   int m, length_rhs; rhs_mat_transp.getSize( m, length_rhs);

   assert( locmyl >= 0 && locmzl >= 0 );
   assert( locnx + locmy + locmz + locmyl + locmzl == length_rhs );

   const int length_reduced = locnx + locmy + locmyl + locmzl;
   if( reduced_rhss_blocked.size() <= static_cast<unsigned int>(n_rhs * length_reduced) )
      reduced_rhss_blocked.resize(n_rhs * length_reduced);

   ///////////////////////////////////////////////////////////////////////
   // LOCAL SOLVE WITH SCHUR COMPLEMENT and set of buffer rhs
   ///////////////////////////////////////////////////////////////////////
   SparseGenMatrix &C = data->getLocalD();

   #pragma omp parallel for schedule(dynamic, 1)
   for( int rhs_i = rhs_start; rhs_i < rhs_start + n_rhs; ++rhs_i )
   {
      assert( rhs_i < m );

      double* rhs_reduced = reduced_rhss_blocked.data() + (rhs_i - rhs_start) * length_reduced;

      SimpleVector b_vec( rhs_mat_transp[rhs_i], length_rhs );

      double *b = b_vec.elements();

      ///////////////////////////////////////////////////////////////////////
      // b = [b1; b2; b3; b4; b5] is a locnx + locmy + locmz + locmyl + locmz vector
      // the new rhs should be
      //           r = [b1-C^T*(zDiag)^{-1}*b3; b2; b4; b5]
      ///////////////////////////////////////////////////////////////////////

      //copy all elements from b into r except for the the residual values corresponding to z0 = b3
      // copy b1, b2
      std::copy(b, b + locnx + locmy, rhs_reduced);

      // copy b4, b5
      std::copy(b + locnx + locmy + locmz, b + locnx + locmy + locmz + locmyl + locmzl, rhs_reduced + locnx + locmy);

      // rhs_reduced now : [ b1; b2; b4; b5 ]

      // alias to r1 part (no mem allocations)
      SimpleVector rhs1(rhs_reduced, locnx);
      SimpleVector b3(b + locnx + locmy, locmz);

      ///////////////////////////////////////////////////////////////////////
      // compute r1 = b1 - C^T * (zDiag)^{-1} * b3
      ///////////////////////////////////////////////////////////////////////
      // if we have C part
      if( locmz > 0 )
      {
         assert( zDiag );
         assert( b3.length() == zDiag->length( ));
         C.transMultD(1.0, rhs1, -1.0, b3, *zDiag);
      }
   }


   ///////////////////////////////////////////////////////////////////////
   // reduced_rhss_blocked now contains all reduced rhs -> solve for it
   ///////////////////////////////////////////////////////////////////////
   if( innerSCSolve == 0 )
      solver->solve(n_rhs, reduced_rhss_blocked.data(), nullptr);
   else
      assert( false && "bicg and iterref not available for blocked solution");

   ///////////////////////////////////////////////////////////////////////
   // rhs_small is now the solution to the reduced system
   // the solution to the augmented system can now be computed as
   //      x = [rhs1; rhs2; zDiag^{-1} * (b3 - C * r1); rhs3; rhs4]
   ///////////////////////////////////////////////////////////////////////


   // copy the solution components and calculate r3
   // copy rhs1 and rhs2
   #pragma omp parallel for schedule(dynamic, 1)
   for( int rhs_i = rhs_start; rhs_i < rhs_start + n_rhs; ++rhs_i )
   {
      double* rhs_reduced = reduced_rhss_blocked.data() + (rhs_i - rhs_start) * length_reduced;

      SimpleVector b_vec( rhs_mat_transp[rhs_i], length_rhs );
      double *b = b_vec.elements();

      std::copy(rhs_reduced, rhs_reduced + locnx + locmy, b);
      // compute x3
      if( locmz > 0 )
      {
         SimpleVector rhs1(rhs_reduced, locnx);
         SimpleVector b3(b + locnx + locmy, locmz);
         C.mult(1.0, b3, -1.0, rhs1);
         b3.componentDiv(*zDiag);
      }

      // copy rhs3 and rhs4
      std::copy(rhs_reduced + locnx + locmy, rhs_reduced + locnx + locmy + locmyl + locmzl, b + locnx + locmy + locmz);
   }
#ifdef TIMING
   // TODO
  if( myRank == 0 && innerSCSolve >= 1 )
    std::cout << "Root - Refin times: child=" << tchild_total << " root=" << troot_total
       << " comm=" << tcomm_total << " total=" << MPI_Wtime()-t_start << "\n";
#endif
}


/** Ht should be either Ft or Gt */
void sLinsysRootAug::addLinkConsBlock0Matrix( DistributedQP *prob, SparseGenMatrix& Ht, int nHtOffsetCols,
      int nKktOffsetCols, int startCol, int endCol)
{
   assert(startCol >= 0 && startCol <= endCol && nKktOffsetCols >= 0 && nKktOffsetCols <= startCol);

   if( startCol == endCol )
      return;

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

   int* const jcolKkt = kkts.jcolM();
   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();
   const double* const MHt = Ht.M();
   const int* const krowHt = Ht.krowM();
   const int* const jcolHt = Ht.jcolM();
   const int n0Links = prob->getN0LinkVars();

   /* main loop going over all rows of Ht */
   for( int i = 0; i < locnx; ++i )
   {
      const bool sparseRow = (i >= locnx - n0Links);

      // note: upper left block ignores 0-link variables pattern, since CtC pattern is not implemented
      int pKkt = krowKkt[i] + locnx - i;

      if( !sparseRow )
         pKkt += nKktOffsetCols;

      assert(pKkt <= krowKkt[i + 1]);
      assert(sparseRow || pKkt == krowKkt[i + 1] || jcolKkt[pKkt] <= startCol);

      if( jcolKkt[pKkt] >= endCol )
      {
#ifndef NDEBUG
         // make sure that there is no entry of Ht in the given range
         int pHt;

         for( pHt = krowHt[i]; pHt < krowHt[i + 1]; pHt++ )
         {
            const int colHt = jcolHt[pHt] + nHtOffsetCols;
            if( colHt >= startCol && colHt < endCol )
               break;
         }

         assert(pHt == krowHt[i + 1]);
#endif
         return;
      }

      bool hit = false;

      // get first in-range entry of Kkt
      for( ; pKkt < krowKkt[i + 1]; pKkt++ )
      {
         const int colKkt = jcolKkt[pKkt];
         if( colKkt >= startCol && colKkt < endCol )
         {
            hit = true;
            break;
         }

         if( colKkt >= endCol )
            break;
      }

      // no entry of Kkt in range?
      if( !hit )
      {
         assert(startCol == endCol || sparseRow);
         continue;
      }

      assert(pKkt < krowKkt[i + 1]);

      int pHt;
      int colHt = -1;
      hit = false;

      // get first in-range entry of Ht
      for( pHt = krowHt[i]; pHt < krowHt[i + 1]; pHt++ )
      {
         colHt = jcolHt[pHt] + nHtOffsetCols;
         if( colHt >= startCol && colHt < endCol )
         {
            hit = true;
            break;
         }

         if( colHt >= endCol )
            break;
      }

      // no entry of Ht in range?
      if( !hit )
         continue;

      assert(colHt >= startCol && colHt < endCol);

      // add in-range entries of Ht to Kkt
      for( ; pKkt < krowKkt[i + 1]; pKkt++ )
      {
         const int colKkt = jcolKkt[pKkt];

         if( colKkt >= endCol )
            break;

         if( colKkt == colHt )
         {
            assert(pHt < krowHt[i + 1]);

            MKkt[pKkt] += MHt[pHt++];

            // end of Ht row reached?
            if( pHt == krowHt[i + 1] )
               break;

            colHt = jcolHt[pHt] + nHtOffsetCols;
         }
      }

      assert(pHt == krowHt[i + 1] || jcolHt[pHt] + nHtOffsetCols >= endCol); // asserts that no entry of Ht has been missed
   }
}


/** rxy = beta*rxy + alpha * SC * x */
void sLinsysRootAug::SCmult( double beta, SimpleVector& rxy,
              double alpha, SimpleVector& x,
              DistributedQP* prob)
{
  //if (iAmDistrib) {
  //only one process subtracts [ (Q+Dx0+C'*Dz0*C)*xx + A'*xy + F'*xxl + G'*xyl ] from r
  //                           [  A*xx                                         ]
  //                           [  F*xx                                         ]
  //                           [  G*xx                           + Omega * xyl ]

  int myRank; MPI_Comm_rank(mpiComm, &myRank);
  int mpiCommSize; MPI_Comm_size(mpiComm, &mpiCommSize);
  const bool iAmLastRank = (myRank == mpiCommSize - 1);
  assert(mpiCommSize >= 1);

  if( iAmLastRank ) { // needs to be the last rank because only this rank is guaranteed to have CtDC
    assert(rxy.length() == locnx + locmy + locmyl + locmzl);

    //only this proc subtracts from rxy
    rxy.scalarMult(beta);
    SparseSymMatrix& Q = prob->getLocalQ();
    Q.mult(1.0,&rxy[0],1, alpha,&x[0],1);

    if(locmz>0) {
      SparseSymMatrix* CtDC_sp = dynamic_cast<SparseSymMatrix*>(CtDC);
      assert(CtDC_sp);

      CtDC_sp->mult(1.0, &rxy[0], 1, -alpha,&x[0], 1);
    }

    SimpleVector& xDiagv = dynamic_cast<SimpleVector&>(*xDiag);
    assert(xDiagv.length() == locnx);
    for(int i=0; i<xDiagv.length(); i++)
      rxy[i] += alpha*xDiagv[i]*x[i];

    SparseGenMatrix& A=prob->getLocalB();
    A.transMult(1.0,&rxy[0],1, alpha,&x[locnx],1);
    A.mult(1.0,&rxy[locnx],1, alpha,&x[0],1);

    assert(locmyl >= 0 && locmzl >= 0);

    if( locmyl > 0 ) {
       SparseGenMatrix& F = prob->getLocalF();
       F.transMult(1.0,&rxy[0],1, alpha,&x[locnx+locmy],1);
       F.mult(1.0,&rxy[locnx+locmy],1, alpha,&x[0],1);
    }

    if( locmzl > 0 ) {
       SparseGenMatrix& G = prob->getLocalG();
       G.transMult(1.0,&rxy[0],1, alpha,&x[locnx+locmy+locmyl],1);
       G.mult(1.0,&rxy[locnx+locmy+locmyl],1, alpha,&x[0],1);

       SimpleVector& zDiagLinkConsv = dynamic_cast<SimpleVector&>(*zDiagLinkCons);
       assert(zDiagLinkConsv.length() == locmzl);
       const int shift = locnx+locmy+locmyl;
       for(int i=0; i<zDiagLinkConsv.length(); i++)
         rxy[i+shift] += alpha*zDiagLinkConsv[i]*x[i+shift];
    }
  } else {
    //other processes set r to zero since they will get this portion from process 0
    rxy.setToZero();
  }

#ifdef TIMING
    taux=MPI_Wtime();
#endif
    // now children add [0 A^T C^T ]*inv(KKT)*[0;A;C] x
    SimpleVector xx((locmyl || locmzl) ? (locnx + locmy + locmyl + locmzl) : locnx);
    xx.copyFromArray(x.elements());
    xx.scalarMult(-alpha);

    for(size_t it=0; it<children.size(); it++) {
      children[it]->addTermToSchurResidual(prob->children[it],rxy,xx);
    }

#ifdef TIMING
    tchild_total +=  (MPI_Wtime()-taux);
#endif
    //~done computing residual

#ifdef TIMING
    taux=MPI_Wtime();
#endif
    //all-reduce residual
    if(iAmDistrib) {
      SimpleVector buf(rxy.length());
      buf.setToZero(); //we use dx as the recv buffer
      MPI_Allreduce(rxy.elements(), buf.elements(), locnx+locmy+locmyl+locmzl, MPI_DOUBLE, MPI_SUM, mpiComm);
      rxy.copyFrom(buf);
    }
#ifdef TIMING
    tcomm_total += (MPI_Wtime()-taux);
#endif

}


void sLinsysRootAug::solveWithIterRef( DistributedQP *prob, SimpleVector& r)
{
  assert( false && " TODO : not sure if working correctly...");
  SimpleVector r2(&r[locnx],       locmy);
  SimpleVector r1(&r[0],           locnx);

  //SimpleVector realRhs(&r[0], locnx+locmy);
#ifdef TIMING
  taux=MPI_Wtime();
#endif

  double rhsNorm=r.twonorm(); //r== the initial rhs of the reduced system here

  int myRank; MPI_Comm_rank(mpiComm, &myRank);
  SimpleVector rxy(locnx+locmy); rxy.copyFrom(r);
  SimpleVector   x(locnx+locmy); x.setToZero(); //solution
  SimpleVector  dx(locnx+locmy);                //update from iter refinement
  SimpleVector x_prev(locnx+locmy);
  int refinSteps=0;
  std::vector<double> histResid;
  int maxRefinSteps=(gLackOfAccuracy>0?9:8);
  do { //iterative refinement
#ifdef TIMING
    taux=MPI_Wtime();
#endif

    x_prev.copyFrom(x);
    //dx = Ainv * r 
    dx.copyFrom(rxy);
    solver->Dsolve(dx);
    //update x
    x.axpy(1.0,dx);

#ifdef TIMING
    troot_total += (MPI_Wtime()-taux);
#endif  

    if(gLackOfAccuracy<0) break;
    if(refinSteps==maxRefinSteps) break;

    //////////////////////////////////////////////////////////////////////
    //iterative refinement
    //////////////////////////////////////////////////////////////////////
    //compute residual
    
    //if (iAmDistrib) {
    //only one process substracts [ (Q+Dx0+C'*Dz0*C)*xx + A'*xy ] from r
    //                            [  A*xx                       ]
    if(myRank==0) {
      rxy.copyFrom(r);
      if(locmz>0) {
	SparseSymMatrix* CtDC_sp = dynamic_cast<SparseSymMatrix*>(CtDC);
	CtDC_sp->mult(1.0,&rxy[0],1, 1.0,&x[0],1);
      }
      SparseSymMatrix& Q = prob->getLocalQ();
      Q.mult(1.0,&rxy[0],1, -1.0,&x[0],1);
      
      SimpleVector& xDiagv = dynamic_cast<SimpleVector&>(*xDiag);
      assert(xDiagv.length() == locnx);
      for(int i=0; i<xDiagv.length(); i++)
	rxy[i] -= xDiagv[i]*x[i];
      
      SparseGenMatrix& A=prob->getLocalB();
      A.transMult(1.0,&rxy[0],1, -1.0,&x[locnx],1);
      A.mult(1.0,&rxy[locnx],1, -1.0,&x[0],1);
    } else {
      //other processes set r to zero since they will get this portion from process 0
      rxy.setToZero();
    }

#ifdef TIMING
    taux=MPI_Wtime();
#endif  
    // now children add [0 A^T C^T ]*inv(KKT)*[0;A;C] x
    SimpleVector xx(&x[0], locnx);
    for(size_t it=0; it<children.size(); it++) {
      children[it]->addTermToSchurResidual(prob->children[it],rxy,xx);  
    }
#ifdef TIMING
    tchild_total +=  (MPI_Wtime()-taux);
#endif
    //~done computing residual 

#ifdef TIMING
    taux=MPI_Wtime();
#endif
    //all-reduce residual
    if(iAmDistrib) {
      dx.setToZero(); //we use dx as the recv buffer
      MPI_Allreduce(rxy.elements(), dx.elements(), locnx+locmy, MPI_DOUBLE, MPI_SUM, mpiComm);
      rxy.copyFrom(dx);
    }
#ifdef TIMING
    tcomm_total += (MPI_Wtime()-taux);
#endif

    double relResNorm=rxy.twonorm()/rhsNorm;
    
    if(relResNorm<1.0e-10) {
      break;
    } else {
      double prevRelResNorm=1.0e10;
      if(histResid.size()) 
	prevRelResNorm=histResid[histResid.size()-1];

      //check for stop, divergence or slow convergence conditions
      if(relResNorm>prevRelResNorm) {
	// diverging; restore iteration
	if(myRank==0) {
	   std::cout << "1st stg - iter refinement diverges relResNorm=" << relResNorm
	       << "  before was " << prevRelResNorm << "\n";
	   std::cout << "Restoring iterate.\n";
	}
	x.copyFrom(x_prev);
	break;
      }else {
	//check slow convergence for the last xxx iterates.
	// xxx is 1 for now
	//if(relResNorm>0.*prevRelResNorm) {

	//  if(myRank==0) {
	//    cout << "1st stg - iter refinement stuck relResNorm=" << relResNorm 
	//	 << "  before was " << prevRelResNorm << endl;
	//    cout << "exiting refinement." << endl;
	//  }
	//  break;
	//
	//} else {
	//  //really nothing, continue
	//}
      }
      histResid.push_back(relResNorm);
      if(myRank==0)
         std::cout << "1st stg - sol does NOT  have enough accuracy (" << relResNorm << ") after "
         << refinSteps << " refinement steps\n";
    }
    refinSteps++;
  }while(refinSteps<=maxRefinSteps);

#ifdef TIMING
  taux = MPI_Wtime();
#endif

  r1.copyFrom(x);
  r2.copyFromArray(&x[locnx]);

#ifdef TIMING
  troot_total += (MPI_Wtime()-taux);
#endif  
}

void sLinsysRootAug::solveWithBiCGStab( DistributedQP *prob, SimpleVector& b)
{
  assert( false && "TODO: not sure if working correctly");
  int n = b.length();

  const int maxit=75; //500
  const double tol=1e-10, EPS=1e-15; // EPS=2e-16

  int myRank; MPI_Comm_rank(mpiComm, &myRank);

  SimpleVector r(n);           //residual
  SimpleVector s(n);           //residual associated with half iterate
  SimpleVector rt(n);          //shadow residual
  SimpleVector xmin(n);        //minimal residual iterate
  SimpleVector x(n);           //iterate
  SimpleVector xhalf(n);       // half iterate of BiCG
  SimpleVector p(n),paux(n);
  SimpleVector v(n), t(n);
  int flag;
  double n2b;                  //norm of b 
  double normr, normrmin;      //norm of the residual and norm of residual at min-resid iterate
  double normr_act;
  double tolb;                 //relative tolerance
  double rho, omega;
  double alpha = -1;
  int stag, maxmsteps, maxstagsteps, moresteps;
  //double imin;
  //maxit = n/2+1;

  //////////////////////////////////////////////////////////////////
  //  Problem Setup and initialization
  //////////////////////////////////////////////////////////////////

  n2b = b.twonorm();
  tolb = n2b * tol;

  tolb = std::max(tolb, EPS);

#ifdef TIMING
  double relres;
  double iter=0.0;
  if( myRank == 0 )
     std::cout << "initial norm of b " << n2b << "\n";
  taux = MPI_Wtime();
#endif
  //initial guess
  x.copyFrom(b);

  solver->Dsolve(x);
  //initial residual
  r.copyFrom(b);

#ifdef TIMING
  troot_total += (MPI_Wtime()-taux);
  taux = MPI_Wtime();
#endif 

  //applyA(1.0, r, -1.0, x);
  SCmult(1.0,r, -1.0,x, prob);

#ifdef TIMING
    tchild_total +=  (MPI_Wtime()-taux);
#endif

  normr = r.twonorm(); normr_act = normr;

  if( normr <= tolb ) {
    //initial guess is good enough
    b.copyFrom(x); flag=0; return;
  }

  if( myRank == 0 )
      std::cout << "innerBICG starts: " << normr << " > " << tolb << "\n";

  rt.copyFrom(r); //Shadow residual
  double* resvec = new double[2*maxit+1];
  resvec[0] = normr; normrmin=normr;
  rho=1.0; omega=1.0;
  stag=0; maxmsteps=std::min(std::min(n/50, 5), n-maxit);
  maxstagsteps=3; moresteps=0;

  //////////////////////////////////////////////////////////////////
  // loop over maxit iterations
  //////////////////////////////////////////////////////////////////
  int ii=0; while(ii<maxit) {
    //cout << ii << " ";
    flag=-1;
    ///////////////////////////////
    // First half of the iterate
    ///////////////////////////////
    double rho1=rho; double beta;
    rho = rt.dotProductWith(r); 
    //printf("rho=%g\n", rho);
    if(0.0==rho) { flag=4;  break; }

    if(ii==0) p.copyFrom(r);
    else {
      beta = (rho/rho1)*(alpha/omega);
      if(beta==0.0) { flag=4;  break; }

      //-------- p = r + beta*(p - omega*v) --------
      p.axpy(-omega, v); p.scale(beta); p.axpy(1.0, r);
    }

#ifdef TIMING
    taux = MPI_Wtime();
#endif
    //------ v = A*(M2inv*(M1inv*p)) and ph=M2inv*(M1inv*p)
    //first use v as temp storage
    //applyM1(0.0, v,    1.0, p);
    //applyM2(0.0, paux, 1.0, v);
    //applyA (0.0, v,    1.0, paux); 
    paux.copyFrom(p);
    solver->solve(paux);

#ifdef TIMING
  troot_total += (MPI_Wtime()-taux);
#endif 
    
    SCmult(0.0,v, 1.0,paux, prob);
    
    SimpleVector& ph = paux;

    double rtv = rt.dotProductWith(v);
    if(rtv==0.0) { flag=4; break; }

    alpha = rho/rtv;
    if(fabs(alpha)*ph.twonorm()<EPS*x.twonorm()) stag++;
    else                                         stag=0;

    // xhalf = x + alpha*ph and the associated residual
    xhalf.copyFrom(x); xhalf.axpy( alpha, ph);
    s.    copyFrom(r);     s.axpy(-alpha, v);
    normr = s.twonorm(); normr_act = normr;
    resvec[2*ii] = normr;

    //printf("iter %g normr=%g\n", ii+0.5, normr);
    //-------- check for convergence in the middle of the iterate.  -------- 
    if(normr<=tolb || stag>=maxstagsteps || moresteps) {
      s.copyFrom(b);
      //applyA(1.0, s, -1.0, xhalf); // s=b-Ax
      SCmult(1.0,s, -1.0,xhalf, prob);
      normr_act = s.twonorm();
      
      if(normr<=tolb) {
	//converged
	x.copyFrom(xhalf);	
	flag = 0;
#ifdef TIMING
	iter = 0.5+ii;
#endif
	break;
      } else {
	if(stag>=maxstagsteps && moresteps==0) {
	  stag=0;
	}
	moresteps++;
	if(moresteps>=maxmsteps) {
	  //method stagnated
	  flag=3; x.copyFrom(xhalf);
	  break;
	}
      }
    }
    if(stag>=maxstagsteps) { flag=3; break;} //stagnation

    //update quantities related to minimal norm iterate
    if(normr_act<normrmin) {
      xmin.copyFrom(xhalf); normrmin=normr_act;
      //imin=0.5+ii;
    }

#ifdef TIMING
    taux = MPI_Wtime();
#endif
    ///////////////////////////////
    // Second half of the iterate
    //////////////////////////////
    //applyM1(0.0, t,    1.0, s); //applyM1(s,     stemp);
    //applyM2(0.0, paux, 1.0, t); //applyM2(stemp, sh);
    //applyA (0.0, t,    1.0, paux); //applyA (sh, t);
    //kkt->mult(0.0,paux, 1.0,s);
    paux.copyFrom(s);
    solver->solve(paux);
#ifdef TIMING
    troot_total += (MPI_Wtime()-taux);
#endif

    SCmult(0.0,t, 1.0,paux, prob);

    SimpleVector& sh = paux; 
    double tt = t.dotProductWith(t);
    if(tt==0.0) { flag=4; break;}

    omega=t.dotProductWith(s); omega /= tt;

    if(fabs(omega)*sh.twonorm() < EPS*xhalf.twonorm()) stag++;
    else                                               stag=0;

    x.copyFrom(xhalf); x.axpy( omega, sh); // x=xhalf+omega*sh
    r.copyFrom(s);     r.axpy(-omega, t ); // r=s-omega*t

    normr = r.twonorm(); normr_act = normr;
    resvec[2*ii+1] = normr;

    //printf("stag=%d  maxstagsteps=%d moresteps=%d  normr=%g\n",
    //	   stag, maxstagsteps, moresteps, normr);    

    //-------- check for convergence at the end of the iterate.  -------- 
    if(normr<=tolb || stag>=maxstagsteps || moresteps) {
      r.copyFrom(b); 
      //applyA(1.0, r, -1.0, x); //r=b-Ax
      SCmult(1.0,r, -1.0,x, prob);
      normr_act=r.twonorm();

      if(normr<=tolb) {
         flag = 0;
#ifdef TIMING
         iter = 1.0+ii;
#endif
         break;
      }
      else {
	if(stag>=maxstagsteps && moresteps==0) {
	  stag = 0;
	}
	moresteps++;
	if(moresteps>=maxmsteps) {
	  //method stagnated
	  flag=3; break;
	}
      }
    } // end convergence check
    if(stag>=maxstagsteps) { flag=3; break;} //stagnation

    //update quantities related to minimal norm iterate
    if(normr_act<normrmin) {
      xmin.copyFrom(x); normrmin=normr_act;
      //imin=1.5+ii;
    }
    //printf("iter %g normr=%g\n", ii+1.0, normr);
    ///////////////////////////////
    // Next iterate
    ///////////////////////////////
    ii++;
    
  }//end while

  if(ii>=maxit) {
#ifdef TIMING
    iter=ii;
#endif
    flag=10;
  }
  
  if(flag==0 || flag==-1) {
#ifdef TIMING
    relres = normr_act/n2b;
    if(myRank==0) {
      printf("INNER BiCGStab converged: normResid=%g relResid=%g iter=%g\n",
	        normr_act, relres, iter);
    }
#endif
  } 
  else 
  {
    if(ii==maxit) flag=10;//aaa
    //FAILURE -> return minimum resid-norm iterate
    r.copyFrom(b); 
    //applyA(1.0, r, -1.0, xmin);
    SCmult(1.0,r, -1.0,xmin, prob);

    normr=r.twonorm();
    if(normr >= normr_act) {
      x.copyFrom(xmin);
      //iter=imin;
#ifdef TIMING
      relres=normr/n2b;
#endif
    } else {
#ifdef TIMING
      iter=1.0+ii;
      relres = normr/n2b;
#endif
    }

#ifdef TIMING
    if(myRank==0) {
      printf("INNERBiCGStab did not NOT converged after %g[%d] iterations.\n", iter,ii);
      printf("\t - Error code %d\n\t - Act res=%g\n\t - Rel res=%g %g\n\n", 
	     flag, normr, relres, normrmin);
    }
#endif
  }

  if( myRank == 0 )
     std::cout << "innerBICG: " << "ii=" << ii << " flag=" << flag << " normr=" << normr << " normr_act="
        << normr_act << " tolb=" << tolb << "\n";

  biCGStabCommunicateStatus(flag, ii);

  b.copyFrom(x);
  delete[] resvec;
}

void sLinsysRootAug::finalizeKKTsparse(DistributedQP* prob, Variables*)
{
   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

#ifndef NDEBUG
   int* const jcolKkt = kkts.jcolM();
#endif
   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();
   const int n0Links = prob->getN0LinkVars();

   assert(kkts.size() == locnx + locmy + locmyl + locmzl);
   assert(!kkts.isLower);
   assert(locmyl >= 0 && locmzl >= 0);

#if 0
   if( PIPS_MPIgetRank(mpiComm) == 0)
   {
      xDiag->writefToStreamStats(std::cout, "xDiag");
      zDiag->writefToStreamStats(std::cout, "zDiag");

      const SimpleVector& szDiagLinkCons = dynamic_cast<const SimpleVector&>(*zDiagLinkCons);
      assert(szDiagLinkCons.length() == locmzl);
      int zerocount = 0;
      for( int i = 0; i < locmzl; i++ )
      {
         if( szDiagLinkCons[i] == 0)
            zerocount++;
      }

      zDiagLinkCons->writefToStreamStats(std::cout, "zDiagLinkCons");

      std::cout << "zDiagLinkCons zeroes: " << zerocount << "\n";
   }
#endif
   //////////////////////////////////////////////////////
   // compute Q+diag(xdiag) - C' * diag(zDiag) * C
   // and update the KKT
   //////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////
   // update the KKT with Q (DO NOT PUT DIAG)
   /////////////////////////////////////////////////////////////
   assert(prob->getLocalQ().krowM()[locnx] == 0 && "Q currently not supported for sparse kkt");

   /////////////////////////////////////////////////////////////
   // update the KKT with the diagonals
   // xDiag is in fact diag(Q)+X^{-1}S
   /////////////////////////////////////////////////////////////
   if( xDiag )
   {
      const SimpleVector& sxDiag = dynamic_cast<const SimpleVector&>(*xDiag);

      for( int i = 0; i < locnx; i++ )
      {
         const int diagIdx = krowKkt[i];
         assert(jcolKkt[diagIdx] == i);

         MKkt[diagIdx] += sxDiag[i];
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with   - C' * diag(zDiag) *C
   /////////////////////////////////////////////////////////////
   if( locmz > 0 )
   {
      assert(zDiag);

      SparseGenMatrix& C = prob->getLocalD();
      C.matTransDinvMultMat(*zDiag, &CtDC);
      assert(CtDC->size() == locnx);

      //aliases for internal buffers of CtDC
      SparseSymMatrix* CtDCsp = dynamic_cast<SparseSymMatrix*>(CtDC);
      const int* krowCtDC = CtDCsp->krowM();
      const int* jcolCtDC = CtDCsp->jcolM();
      const double* dCtDC = CtDCsp->M();
      for( int i = 0; i < locnx; i++ )
      {
         const int pend = krowCtDC[i + 1];
         for( int p = krowCtDC[i]; p < pend; p++ )
         {
            const int col = jcolCtDC[p];

            if( col >= i )
            {
               // get start position of dense kkt block
               const int blockStart = krowKkt[i];
               assert(col < locnx && jcolKkt[blockStart + col - i] == col);

               MKkt[blockStart + col - i] -= dCtDC[p];
            }
         }
      }
   } //~end if locmz>0

   /////////////////////////////////////////////////////////////
   // update the KKT with At (symmetric update forced)
   /////////////////////////////////////////////////////////////
   if( locmy > 0 )
   {
      SparseGenMatrix& At = prob->getLocalB().getTranspose(); // yes, B
      const double* MAt = At.M();
      const int* krowAt = At.krowM();

      for( int i = 0; i < locnx; ++i )
      {
         const int pstart = krowAt[i];
         const int pend = krowAt[i + 1];

         // get start position of sparse kkt block
         const int blockStart = krowKkt[i] + locnx - i;

         assert(blockStart <= krowKkt[i + 1]);

         for( int p = pstart; p < pend; ++p )
         {
            assert(At.jcolM()[p] < locmy);
            assert(blockStart + (p - pstart) <= krowKkt[i + 1]);
            assert(jcolKkt[blockStart + (p - pstart)] == (locnx + At.jcolM()[p]));

            MKkt[blockStart + (p - pstart)] += MAt[p];
         }
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with Ft
   /////////////////////////////////////////////////////////////
   if( locmyl > 0 )
   {
      SparseGenMatrix& Ft = prob->getLocalF().getTranspose();
      const double* MFt = Ft.M();
      const int* krowFt = Ft.krowM();
      const int* jcolFt = Ft.jcolM();

      int* krowGt = nullptr;

      if( locmzl > 0 )
      {
         SparseGenMatrix& Gt = prob->getLocalG().getTranspose();
         krowGt = Gt.krowM();
      }

      for( int i = 0; i < locnx; ++i )
      {
         const bool sparseRow = (i >= locnx - n0Links);
         const int pend = krowFt[i + 1];

         if( sparseRow )
         {
            int blockStart =  krowKkt[i + 1] - (krowFt[i + 1] - krowFt[i]);

            if( locmzl > 0 )
               blockStart -= (krowGt[i + 1] - krowGt[i]);

            assert(blockStart >= krowKkt[i]);

            for( int p = krowFt[i], shift = 0; p < pend; ++p, ++shift )
            {
               assert(jcolFt[p] < locmyl && jcolKkt[blockStart + shift] == (locnx + locmy + jcolFt[p]));

               MKkt[blockStart + shift] += MFt[p];
            }
         }
         else
         {
            const int blockStart = krowKkt[i + 1] - locmyl - locmzl;
            assert(blockStart >= krowKkt[i]);

            for( int p = krowFt[i]; p < pend; ++p )
            {
               const int col = jcolFt[p];
               assert(col < locmyl && jcolKkt[blockStart + col] == (locnx + locmy + col));

               MKkt[blockStart + col] += MFt[p];
            }
         }
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with Gt and add z diagonal
   /////////////////////////////////////////////////////////////
   if( locmzl > 0 )
   {
      SparseGenMatrix& Gt = prob->getLocalG().getTranspose();
      const double* MGt = Gt.M();
      const int* krowGt = Gt.krowM();
      const int* jcolGt = Gt.jcolM();

      for( int i = 0; i < locnx; ++i )
      {
         const int pend = krowGt[i + 1];
         const bool sparseRow = (i >= locnx - n0Links);

         if( sparseRow )
         {
            const int blockStart = krowKkt[i + 1] - (krowGt[i + 1] - krowGt[i]);

            assert(blockStart >= krowKkt[i]);

            for( int p = krowGt[i], shift = 0; p < pend; ++p, ++shift )
            {
               assert(jcolGt[p] < locmzl && jcolKkt[blockStart + shift] == (locnx + locmy + locmyl + jcolGt[p]));

               MKkt[blockStart + shift] += MGt[p];
            }
         }
         else
         {
            const int blockStart = krowKkt[i + 1] - locmzl;

            assert(blockStart >= krowKkt[i]);

            for( int p = krowGt[i]; p < pend; ++p )
            {
               const int col = jcolGt[p];
               assert(col < locmzl && jcolKkt[blockStart + col] == (locnx + locmy + locmyl + col));

               MKkt[blockStart + col] += MGt[p];
            }
         }
      }

      assert(zDiagLinkCons);

      const SimpleVector& szDiagLinkCons = dynamic_cast<const SimpleVector&>(*zDiagLinkCons);

      for( int i = 0, iKkt = locnx + locmy + locmyl; i < locmzl; ++i, ++iKkt )
      {
         const int idx = krowKkt[iKkt];
         assert(jcolKkt[idx] == iKkt);

         MKkt[idx] += szDiagLinkCons[i];
      }
   }

#ifdef DUMPKKT
   ofstream myfile;
   myfile.open("../sparsekkt");

   int zerocount = 0;
   const int sizeKkt = locnx + locmy + locmyl + locmzl;

   for( int r = 0; r < sizeKkt; r++ )
   {
      for( int i = krowKkt[r]; i < krowKkt[r + 1]; i++ )
      {
         const double val = MKkt[i];
         const double col = jcolKkt[i];
         if( val != 0.0 )
            myfile << r << " " << col << " " << val << "\n";
         else
            zerocount++;
      }
   }

   std::cout << "zero-count " << zerocount << " of " << krowKkt[sizeKkt] << "\n";

   myfile.close();

   assert(0);

#endif
}

void sLinsysRootAug::finalizeKKTdense(DistributedQP* prob, Variables*)
{
   int j, p, pend;

   DenseSymMatrix* const kktd = dynamic_cast<DenseSymMatrix*>(kkt.get());

   //alias for internal buffer of kkt
   double** const dKkt = kktd->Mat();

   //////////////////////////////////////////////////////
   // compute Q+diag(xdiag) - C' * diag(zDiag) * C
   // and update the KKT
   //////////////////////////////////////////////////////


   /////////////////////////////////////////////////////////////
   // update the KKT with Q (DO NOT PUT DIAG)
   /////////////////////////////////////////////////////////////
   SparseSymMatrix& Q = prob->getLocalQ();
   int* krowQ=Q.krowM(); int* jcolQ=Q.jcolM(); double* dQ=Q.M();
   for(int i=0; i<locnx; i++) {
     pend = krowQ[i+1];
     for(p=krowQ[i]; p<pend; p++) {
       j = jcolQ[p];
       if(i==j) continue;
       double val = dQ[p];
       dKkt[i][j] += val;
       dKkt[j][i] += val;

       assert(0 && "non-empty Q currently not supported");
     }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with the diagonals
   // xDiag is in fact diag(Q)+X^{-1}S
   /////////////////////////////////////////////////////////////
   //kktd->atPutDiagonal( 0, *xDiag );
   if( xDiag )
   {
      SimpleVector& sxDiag = dynamic_cast<SimpleVector&>(*xDiag);
      for(int i=0; i<locnx; i++) dKkt[i][i] += sxDiag[i];
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with   - C' * diag(zDiag) *C
   /////////////////////////////////////////////////////////////
   if(locmz>0) {
     assert(zDiag);

     SparseGenMatrix& C = prob->getLocalD();
     C.matTransDinvMultMat(*zDiag, &CtDC);

     assert(CtDC->size() == locnx);

     //aliases for internal buffers of CtDC
     SparseSymMatrix* CtDCsp = reinterpret_cast<SparseSymMatrix*>(CtDC);
     int* krowCtDC=CtDCsp->krowM(); int* jcolCtDC=CtDCsp->jcolM(); double* dCtDC=CtDCsp->M();

     for( int i = 0; i < locnx; i++ )
     {
       pend = krowCtDC[i + 1];
       for( p = krowCtDC[i]; p < pend; p++ )
       {
          j = jcolCtDC[p];

          if( j <= i )
             dKkt[i][j] -= dCtDC[p];
        }
     }
   } //~end if locmz>0
   /////////////////////////////////////////////////////////////
   // update the KKT with A (symmetric update forced)
   /////////////////////////////////////////////////////////////
   if(locmy>0)
   {
     SparseGenMatrix& A = prob->getLocalB(); // yes, B
     const double* dA = A.M();
     const int* krowA = A.krowM();
     const int* jcolA = A.jcolM();

     int iKkt = locnx;
     for( int i = 0; i < locmy; ++i, ++iKkt ) {

       for( p = krowA[i], pend = krowA[i + 1]; p < pend; ++p ) {
         j = jcolA[p];
         assert(j < locnx);

         dKkt[iKkt][j] += dA[p];
       }
     }
   }
   //prob->getLocalB().getStorageRef().dump("stage1eqmat2.dump");


   /////////////////////////////////////////////////////////////
   // update the KKT with F
   /////////////////////////////////////////////////////////////
   if( locmyl > 0 )
   {
     SparseGenMatrix& F = prob->getLocalF();
     const double* dF = F.M();
     const int* krowF = F.krowM();
     const int* jcolF = F.jcolM();

     int iKkt = locnx + locmy;
     for( int i = 0; i < locmyl; ++i, ++iKkt ) {
       for( p = krowF[i], pend = krowF[i+1]; p < pend; ++p ) {
         j = jcolF[p];
         assert(j < locnx);

         const double val = dF[p];
         dKkt[iKkt][j] += val;
       }
     }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with G and put z diagonal
   /////////////////////////////////////////////////////////////
   if( locmzl > 0 )
   {
     SparseGenMatrix& G = prob->getLocalG();
     assert(zDiagLinkCons);
     SimpleVector& szDiagLinkCons = dynamic_cast<SimpleVector&>(*zDiagLinkCons);

     const double* dG = G.M();
     const int* krowG = G.krowM();
     const int* jcolG = G.jcolM();

     int iKkt = locnx + locmy + locmyl;
     for( int i = 0; i < locmzl; ++i, ++iKkt ) {

       dKkt[iKkt][iKkt] += szDiagLinkCons[i];
       for( p = krowG[i], pend = krowG[i+1]; p < pend; ++p ) {
         j = jcolG[p];
         assert(j < locnx);

         const double val = dG[p];
         dKkt[iKkt][j] += val;
       }
     }
   }

#ifdef DUMPKKT
   const int msize = locnx + locmy + locmyl + locmzl;

   ofstream myfile;
   myfile.open("../densekkt");

   for( int col = 0; col < msize; col++ )
      for( int row = col; row < msize; row++ )
         if( dKkt[row][col] != 0.0 )
            myfile << col << " " << row << " " << dKkt[row][col] << "\n";

   myfile.close();

   assert(0);
#endif

   /////////////////////////////////////////////////////////////
   // update the KKT zeros for the lower right block
   /////////////////////////////////////////////////////////////
   //kktd->storage().atPutZeros(locnx, locnx, locmy+locmz, locmy+locmz);
   //myAtPutZeros(kktd, locnx, locnx, locmy, locmy);
}

void sLinsysRootAug::DsolveHierarchyBorder( DenseGenMatrix& rhs_mat_transp, int n_cols )
{
   /* b holds all rhs in transposed form - C part from schur complement is already missing in b */
   const int my_rank = PIPS_MPIgetRank( mpiComm );
#ifdef TIMING
   // TODO
#endif

   int m, n; rhs_mat_transp.getSize(m, n);
#ifndef NDEBUG
   assert(locmyl >= 0 && locmzl >= 0);

   assert( n_cols <= m );
   assert( locnx + locmy + locmz + locmyl + locmzl == n );
#endif

   /*
    * for every right hand side one of the processes now does the SC solve operation and puts it at the corresponding
    * position in b
    * Every process has to do n_rhs / n_procs right hand sides while the first few might have to solve with one additional one
    */
   const int size = PIPS_MPIgetSize( mpiComm );
   const int n_blockrhs = static_cast<int>( n_cols / size );
   const int leftover = n_cols % size;

   const int n_rhs = (my_rank < leftover ) ? n_blockrhs + 1 : n_blockrhs;
   const int rhs_start = my_rank < leftover ? (n_blockrhs + 1) * my_rank :
         (n_blockrhs + 1) * leftover + (my_rank - leftover) * n_blockrhs;

   assert( rhs_start <= n_cols );
   assert( rhs_start + n_rhs <= n_cols );

   // set rhs contributed by other procs to zero
   #pragma omp parallel for schedule(dynamic, 1)
   for( int rhs_i = 0; rhs_i < n_cols; ++rhs_i )
   {
      if( rhs_start <= rhs_i && rhs_i < rhs_start + n_rhs )
         continue;
      SimpleVector b( rhs_mat_transp[rhs_i], n );
      b.setToZero();
   }

   solveReducedLinkConsBlocked( data, rhs_mat_transp, rhs_start, n_rhs );

   if( iAmDistrib )
   {
      // TODO only allreduce relevant part
      // TODO is allreduce even worth it here? Every proc could also compute all its rhs - add if n_rhs big
      int m, n;
      rhs_mat_transp.getSize(m, n);
      submatrixAllReduceFull(&rhs_mat_transp, 0, 0, m, n, mpiComm);
   }
}

void sLinsysRootAug::addBlTKiInvBrToRes( DoubleMatrix& result, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sym_res, bool sparse_res )
{
   assert( !is_hierarchy_root );

   if( Bl.isEmpty() || (Br.isEmpty() && Br_mod_border.empty() ) )
      return;

   int dummy, m_result;
   result.getSize( m_result, dummy );

   assert( m_result > 0 );
   assert( blocksize_hierarchical > 0 );

   // buffer b0 for blockwise computation of Br0 - SUM_i  Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij ), stored in transposed form (for quick access of cols in solve)
   // dense since we have no clue about any structure in the system and Xij are dense
   const int n_buffer = locnx + locmy + locmz + locmyl + locmzl;
   const int m_buffer = allocateAndZeroBlockedComputationsBuffer(m_result, n_buffer);

   const size_t n_chunks = std::ceil( static_cast<double>(m_result) / m_buffer );

   assert( n_chunks > 0 );
   if( !sc_compute_blockwise_hierarchical )
   {
      assert( n_chunks == 1 );
      assert( m_buffer == m_result );
   }

   for( size_t i = 0; i < n_chunks; ++i )
   {
      const int begin_chunk = i * m_buffer;
      const int end_chunk = std::min( static_cast<size_t>(m_result), (i + 1) * m_buffer );

      assert( end_chunk - begin_chunk <= m_buffer );
      addBlTKiInvBrToResBlockwise( result, Bl, Br, Br_mod_border, sym_res, sparse_res, *buffer_blocked_hierarchical, begin_chunk, end_chunk );
   }
}

//TODO: determine start and end of the two-links (exclusively on first and last process);
//TODO: probably separate method for two-links?

/* if Bl^T is a two link border this needs to be done for at most the first and the last child of this communicator (all the other children are guaranteed 0) */
/* compute res += [ Bl^T Ki^-1 (Br - sum_j Bmodj Xj) ]^T = (Br^T - SUM_j Xj^T Bmodj^T) Ki^-1 Bl */
void sLinsysRootAug::addBlTKiInvBrToResBlockwise( DoubleMatrix& result, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
      bool sym_res, bool sparse_res, DenseGenMatrix& buffer_b0, int begin_cols, int end_cols )
{
   /* only called on sLinsysRootBordered and sLinsysRootAugHierInner */
   buffer_b0.putZeros();

   /* buffer_b0 is in transposed for so that we can access its cols (or rows in the storage) quickly in Dsolve */
   assert( !is_hierarchy_root );
   assert( 0 <= begin_cols && begin_cols <= end_cols );
   assert( buffer_b0.getN() == kkt->size() || buffer_b0.getN() == kkt->size() + locmz );
   assert( buffer_b0.getM() >= end_cols - begin_cols );

   const bool two_link_border_left = !(Bl.has_RAC || Bl.use_local_RAC);
   const bool two_link_border_right = !(Br.has_RAC || Br.use_local_RAC);

   if( two_link_border_right );

   /* compute Schur Complement right hand sides SUM_i Bi_{inner} Ki^-1 ( Bri - sum_j Bmodij Xij )
    * (keep in mind that in Bi_{this} and the SC we projected C0 Omega0 out) */
   // buffer_b0 = - [ SUM_i Bi_{inner}^T Ki^{-1} (Bri - SUM_j Bmodij Xij) ]^T = SUM_i (Bri^T - SUM_j Xij^T Bmodij^T) Ki^{-1} Bi_{inner}
   LsolveHierarchyBorder( buffer_b0, Br, Br_mod_border, two_link_border_left, begin_cols, end_cols );

   // buffer_b0 = (Br0 - sum_j Bmod0J X0j ) - buffer_b0 = Br0 - sum_j Bmod0J X0j - SUM_i Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij )}
   finalizeZ0Hierarchical( buffer_b0, Br, Br_mod_border, begin_cols, end_cols );

   // solve with Schur Complement for B0_{outer} - SUM_i Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij ) (stored in transposed form!)
   // buffer_b0 = SC_{inner}^-1 buffer_b0 = X0
   DsolveHierarchyBorder( buffer_b0, end_cols - begin_cols );

   // compute result += -SUM_i Bli^T Ki^{-1} ( ( Bri - sum_j Bmodij Xij )  - Bi_{inner} X0 ) += -SUM_i Bli^T Xi
   LtsolveHierarchyBorder( result, buffer_b0, Bl, Br, Br_mod_border, sym_res, sparse_res, begin_cols, end_cols );

   // compute result += Bl0^T X0
   if( PIPS_MPIgetRank(mpiComm) == 0 )
      finalizeInnerSchurComplementContribution( result, buffer_b0, Bl, sym_res, sparse_res, begin_cols, end_cols );
}
