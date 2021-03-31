/*
 * sLinsysRootAugHierInner.C
 *
 *  Created on: 27.01.2021
 *      Author: bzfkempk
 */

#include "sLinsysRootAugHierInner.h"

sLinsysRootAugHierInner::sLinsysRootAugHierInner(sFactory *factory,
      sData *prob_, OoqpVector *dd_, OoqpVector *dq_, OoqpVector *nomegaInv_, OoqpVector *rhs_) :
      sLinsysRootAug(factory, prob_, dynamic_cast<StochVector*>(dd_)->vec,
            dynamic_cast<StochVector*>(dq_)->vec,
            dynamic_cast<StochVector*>(nomegaInv_)->vec, rhs_, false)
{
   assert( locnx == 0 );
   assert( locmy == 0 );
   assert( locmz == 0 );

   createSolversAndKKts(prob_);
}

void sLinsysRootAugHierInner::createSolversAndKKts(sData* prob)
{
   assert( hasSparseKkt );

   const SolverType solver_sub_root = pips_options::getSolverSubRoot();

   static bool printed = false;
   if( !printed && PIPS_MPIgetRank() == 0 )
      std::cout << "sLinsysRootAugHierInner: using " << solver_sub_root << "\n";

   kkt.reset(createKKT(prob));

   if( !printed && PIPS_MPIgetRank() == 0 )
      std::cout << "sLinsysRootAugHierInner: getSchurCompMaxNnz " << prob->getSchurCompMaxNnz() << "\n";
   printed = true;

   createSolversSparse(solver_sub_root);
}

void sLinsysRootAugHierInner::assembleLocalKKT( sData* prob )
{
   for( size_t c = 0; c < children.size(); ++c )
   {
#ifdef STOCH_TESTING
      g_scenNum = c;
#endif
      if( children[c]->mpiComm == MPI_COMM_NULL )
         continue;

      children[c]->stochNode->resMon.recFactTmChildren_start();
      //---------------------------------------------
      addTermToSchurCompl(prob, c, false);
      //---------------------------------------------
      children[c]->stochNode->resMon.recFactTmChildren_stop();
   }
}

void sLinsysRootAugHierInner::Ltsolve( sData *prob, OoqpVector& x )
{
   StochVector& b = dynamic_cast<StochVector&>(x);
   SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);

   //dumpRhs(0, "sol",  b0);
   SimpleVector& z0 = b0; //just another name, for clarity

   for(size_t it = 0; it < children.size(); it++)
      children[it]->Ltsolve2(prob->children[it], *b.children[it], z0, false);
}

void sLinsysRootAugHierInner::Ltsolve2(sData*, StochVector& x, SimpleVector& x0, bool use_local_RAC )
{
   assert( pips_options::getBoolParameter("HIERARCHICAL") );

   StochVector& b = dynamic_cast<StochVector&>(x);

   computeInnerSystemRightHandSide( b, x0, use_local_RAC );
   solveCompressed( x );
}

void sLinsysRootAugHierInner::LsolveHierarchyBorder( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool two_link_border, int begin_cols, int end_cols  )
{
   LsolveHierarchyBorder( result, Br, Br_mod_border, false, two_link_border, begin_cols, end_cols );
}

void sLinsysRootAugHierInner::LtsolveHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br,
      std::vector<BorderMod>& br_mod_border, bool sym_res, bool sparse_res, bool two_link_border )
{
   if( Bl.isEmpty() || (Br.isEmpty() && br_mod_border.empty()) )
      return;

   LtsolveHierarchyBorder( res, X0, Bl, Br, br_mod_border, sym_res, sparse_res, false, two_link_border );
}

void sLinsysRootAugHierInner::computeInnerSystemRightHandSide( StochVector& rhs_inner, const SimpleVector& b0, bool use_local_RAC )
{
   BorderLinsys Border( dynamic_cast<StringGenMatrix&>(*dynamic_cast<StochGenMatrix&>(*data->A).Blmat),
         dynamic_cast<StringGenMatrix&>(*dynamic_cast<StochGenMatrix&>(*data->C).Blmat), use_local_RAC );

   if( Border.isEmpty() )
      return;

   addBorderX0ToRhs( rhs_inner, b0, Border );
}

/* compute Schur rhs b0 - sum Bi^T Ki^-1 bi for all children */
void sLinsysRootAugHierInner::Lsolve(sData *prob, OoqpVector& x )
{
   assert( !is_hierarchy_root );

   StochVector& b = dynamic_cast<StochVector&>(x);
   assert(children.size() == b.children.size() );

   SimpleVector& b0 = dynamic_cast<SimpleVector&>(*b.vec);
   assert(!b.vecl);

   if( iAmDistrib && PIPS_MPIgetRank(mpiComm) > 0 )
      b0.setToZero();

   // compute Bi^T Ki^-1 rhs_i and sum it up
   for( size_t it = 0; it < children.size(); it++ )
      children[it]->addLniziLinkCons( prob->children[it], b0, *b.children[it], false );

   if(iAmDistrib)
      PIPS_MPIsumArrayInPlace( b0.elements(), b0.length(), mpiComm );
}

void sLinsysRootAugHierInner::addLniziLinkCons( sData*, OoqpVector& z0_, OoqpVector& zi, bool use_local_RAC )
{
   assert( zi.isKindOf(kStochVector) );
   SimpleVector& z0 = dynamic_cast<SimpleVector&>(z0_);

   if( !sol_inner )
      sol_inner.reset(dynamic_cast<StochVector*>(zi.cloneFull()));
   else
      sol_inner->copyFrom(zi);

   /* solve system */
   solveCompressed( *sol_inner );

   BorderLinsys Bl( dynamic_cast<StringGenMatrix&>(*dynamic_cast<StochGenMatrix&>(*data->A).Blmat),
         dynamic_cast<StringGenMatrix&>(*dynamic_cast<StochGenMatrix&>(*data->C).Blmat), use_local_RAC );

   addBorderTimesRhsToB0( *sol_inner, z0, Bl );
}

void sLinsysRootAugHierInner::addBorderTimesRhsToB0( StochVector& rhs, SimpleVector& b0, BorderLinsys& border )
{
   assert( rhs.children.size() == children.size() );
   assert( border.F.children.size() == children.size() );
   for( size_t i = 0; i < children.size(); ++i )
   {
      BorderLinsys child_border = getChild(border, i);
      if( child_border.isEmpty() )
         continue;

      children[i]->addBorderTimesRhsToB0( *rhs.children[i], b0, child_border );
   }

   /* add schur complement part */
   if( (border.has_RAC || border.use_local_RAC) && (PIPS_MPIgetSize(mpiComm) == 0 || PIPS_MPIgetRank(mpiComm) == 0) )
   {
      if( border.has_RAC )
      {
         assert( border.A.mat_link );
         assert( border.C.mat_link );
      }

      SparseGenMatrix& F_border = border.has_RAC ? dynamic_cast<SparseGenMatrix&>(*border.A.mat_link).getTranspose() : data->getLocalF().getTranspose();
      SparseGenMatrix& G_border = border.has_RAC ? dynamic_cast<SparseGenMatrix&>(*border.C.mat_link).getTranspose() : data->getLocalG().getTranspose();

      int mFb, nFb; F_border.getSize(mFb, nFb);
      int mGb, nGb; G_border.getSize(mGb, nGb);

      assert( mFb == mGb );
      assert( rhs.vec );
      assert( rhs.vec->length() == nFb + nGb );

      assert( b0.length() >= mFb );

      SimpleVector& zi = dynamic_cast<SimpleVector&>(*rhs.vec);

      SimpleVector zi1 (&zi[0], nFb );
      SimpleVector zi2 (&zi[nFb], nGb );

      SimpleVector b1( &b0[0], mFb );

      F_border.mult( 1.0, b1, -1.0, zi1 );
      G_border.mult( 1.0, b1, -1.0, zi2 );
   }
}

void sLinsysRootAugHierInner::addBorderX0ToRhs( StochVector& rhs, const SimpleVector& x0, BorderLinsys& border )
{
   assert( rhs.children.size() == children.size() );
   assert( border.F.children.size() == children.size() );

   for( size_t i = 0; i < children.size(); ++i )
   {
      BorderLinsys child_border = getChild(border, i);
      if( child_border.isEmpty() )
         continue;
      children[i]->addBorderX0ToRhs( *rhs.children[i], x0, child_border );
   }

   if( border.has_RAC )
   {
      assert( border.A.mat_link );
      assert( border.C.mat_link );
   }

   /* add schur complement part */
   if( border.has_RAC || border.use_local_RAC )
   {
      SparseGenMatrix& F_border = border.has_RAC ? dynamic_cast<SparseGenMatrix&>(*border.A.mat_link) : data->getLocalF();
      SparseGenMatrix& G_border = border.has_RAC ? dynamic_cast<SparseGenMatrix&>(*border.C.mat_link) : data->getLocalG();
      int mFb, nFb; F_border.getSize(mFb, nFb);
      int mGb, nGb; G_border.getSize(mGb, nGb);

      assert( rhs.vec );
      assert( rhs.vec->length() == mFb + mGb );
      assert( nFb == nGb );
      assert( x0.length() >= nFb );

      SimpleVector& rhs0 = dynamic_cast<SimpleVector&>(*rhs.vec);

      SimpleVector rhs01 (&rhs0[0], mFb );
      SimpleVector rhs02 (&rhs0[mFb], mGb );

      SimpleVector x1( (double*)&x0[0], nFb );

      F_border.mult( 1.0, rhs01, -1.0, x1 );
      G_border.mult( 1.0, rhs02, -1.0, x1 );
   }
}

void sLinsysRootAugHierInner::addInnerBorderKiInvBrToRes( DoubleMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool use_local_RAC, bool sparse_res, bool sym_res, int begin_cols, int end_cols )
{
   assert( dynamic_cast<StochGenMatrix&>(*data->A).Blmat->isKindOf(kStringGenMatrix) );
   assert( dynamic_cast<StochGenMatrix&>(*data->C).Blmat->isKindOf(kStringGenMatrix) );

   BorderLinsys Bl( data->getLocalFBorder(), data->getLocalGBorder(), use_local_RAC );
   if( Bl.isEmpty() || (Br.isEmpty() && Br_mod_border.empty()) )
      return;
   addBTKiInvBToSC( result, Bl, Br, Br_mod_border, sparse_res, sym_res);

}

void sLinsysRootAugHierInner::addTermToSchurComplBlocked(sData* prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC )
{
   assert( data == prob );

   BorderLinsys Bl( prob->getLocalFBorder(), prob->getLocalGBorder(), use_local_RAC );
   BorderLinsys Br( prob->getLocalFBorder(), prob->getLocalGBorder(), use_local_RAC );
   if( Bl.isEmpty() )
      return;

   std::vector<BorderMod> border_mod;

   addBTKiInvBToSC( SC, Bl, Br, border_mod, true, sparseSC );
}

void sLinsysRootAugHierInner::LniTransMultHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0,
      BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sparse_res, bool sym_res, bool use_local_RAC, int begin_cols, int end_cols )
{

   if( Bl.isEmpty() || (Br.isEmpty() && Br_mod_border.empty()) )
      return;

   BorderLinsys B_inner( data->getLocalFBorder(), data->getLocalGBorder(), use_local_RAC );
   std::vector<BorderMod> border_mod( Br_mod_border );
   if( !B_inner.isEmpty() )
   {
      BorderMod B_inner_mod(B_inner, X0);
      border_mod.push_back( B_inner_mod );
   }

   const int n_buffer = locnx + locmy + locmz + locmyl + locmzl;
   const size_t m_buffer = PIPSgetnOMPthreads() * blocksize_hierarchical;

   // buffer b0 for blockwise computation of Br0 - SUM_i  Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij ), stored in transposed form (for quick access of cols in solve)
   // dense since we have no clue about any structure in the system and Xij are dense
   if( !buffer_blocked_hierarchical )
      buffer_blocked_hierarchical.reset( new DenseGenMatrix(m_buffer, n_buffer) );
   buffer_blocked_hierarchical->atPutZeros(0, 0, m_buffer, n_buffer);

   assert( end_cols - begin_cols <= m_buffer );

   addBTKiInvBToSCBlockwise( res, Bl, Br, border_mod, sym_res, sparse_res, *buffer_blocked_hierarchical, begin_cols, end_cols );
}

void sLinsysRootAugHierInner::putXDiagonal( OoqpVector& xdiag_ )
{
  assert( dynamic_cast<StochVector&>(xdiag_).vec->isKindOf(kStochVector) );
  StochVector& xdiag = dynamic_cast<StochVector&>(*dynamic_cast<StochVector&>(xdiag_).vec);

  assert(children.size() == xdiag.children.size());

  xDiag = xdiag.vec;

  for(size_t it = 0; it < children.size(); it++)
    children[it]->putXDiagonal(*xdiag.children[it]);
}


void sLinsysRootAugHierInner::putZDiagonal( OoqpVector& zdiag_ )
{
  assert( dynamic_cast<StochVector&>(zdiag_).vec->isKindOf(kStochVector) );
  StochVector& zdiag = dynamic_cast<StochVector&>(*dynamic_cast<StochVector&>(zdiag_).vec);

  assert(children.size() == zdiag.children.size());

  zDiag = zdiag.vec;
  zDiagLinkCons = zdiag.vecl;

  for(size_t it = 0; it < children.size(); it++)
    children[it]->putZDiagonal(*zdiag.children[it]);
}
