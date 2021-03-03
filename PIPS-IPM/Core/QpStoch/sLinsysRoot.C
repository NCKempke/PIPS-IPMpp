/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#include "sLinsysRoot.h"
#include "sFactory.h"
#include "sData.h"
#include "sDummyLinsys.h"
#include "sLinsysLeaf.h"
#include "StochOptions.h"
#include "math.h"

#include "pipsport.h"

/*********************************************************************/
/************************** ROOT *************************************/
/*********************************************************************/

#ifdef STOCH_TESTING
double g_scenNum;
#endif

sLinsysRoot::sLinsysRoot(sFactory * factory_, sData * prob_, bool is_hierarchy_root)
  : sLinsys(factory_, prob_, is_hierarchy_root)
{
  if( pips_options::getBoolParameter( "HIERARCHICAL" ) )
    assert( is_hierarchy_root );
  init();
}

sLinsysRoot::sLinsysRoot(sFactory* factory_,
			 sData* prob_,
			 OoqpVector* dd_,
			 OoqpVector* dq_,
			 OoqpVector* nomegaInv_,
			 OoqpVector* rhs_)
  : sLinsys(factory_, prob_, dd_, dq_, nomegaInv_, rhs_, true)
{
   init();
}

void sLinsysRoot::init()
{
  createChildren(data);

  precondSC = SCsparsifier(mpiComm);
  usePrecondDist = pips_options::getBoolParameter("PRECONDITION_DISTRIBUTED");

  // use sparse KKT if (enough) 2 links are present
  hasSparseKkt = data->exploitingLinkStructure();
  allreduce_kkt = pips_options::getBoolParameter("ALLREDUCE_SCHUR_COMPLEMENT");
  if( pips_options::getBoolParameter( "HIERARCHICAL" ) )
     assert( allreduce_kkt );

  usePrecondDist = usePrecondDist && hasSparseKkt && iAmDistrib;
  MatrixEntryTriplet_mpi = MPI_DATATYPE_NULL;

  initProperChildrenRange();
}

sLinsysRoot::~sLinsysRoot()
{
  for(size_t c = 0; c < children.size(); c++)
    delete children[c];

  delete kktDist;

  delete[] sparseKktBuffer;
}

//this variable is just reset in this file; children will default to the "safe" linear solver
extern int gLackOfAccuracy;

void sLinsysRoot::assembleKKT(sData* prob, Variables* vars)
{
   if( is_hierarchy_root )
      assert( children.size() == 1 );

   /* set kkt to zero */
   initializeKKT(prob, vars);

   /* important that int separate loops! else block in Allreduce might occur */
   for(size_t c = 0; c < children.size(); c++)
      children[c]->assembleKKT(prob->children[c], vars);
   for(size_t c = 0; c < children.size(); c++)
      children[c]->allreduceAndFactorKKT(prob->children[c], vars);

   /* build KKT from local children */
   assembleLocalKKT( prob );
}

void sLinsysRoot::allreduceAndFactorKKT(sData* prob, Variables* vars)
{
   reduceKKT(prob);

   finalizeKKT(prob, vars);

   factorizeKKT(prob);
}

void sLinsysRoot::factor2(sData *prob, Variables *vars)
{
   if( is_hierarchy_root )
      assert( children.size() == 1 );

   /* set kkt to zero */
   initializeKKT(prob, vars);

   // First tell children to factorize.
   for(size_t c = 0; c < children.size(); c++)
      children[c]->factor2(prob->children[c], vars);

   /* build KKT from local children */
   assembleLocalKKT( prob );

#ifdef TIMING
   MPI_Barrier(MPI_COMM_WORLD);
   stochNode->resMon.recReduceTmLocal_start();
#endif 

   reduceKKT(prob);

 #ifdef TIMING
   stochNode->resMon.recReduceTmLocal_stop();
#endif

   finalizeKKT(prob, vars);

   factorizeKKT(prob);
#ifdef TIMING
   afterFactor();
#endif
}

#ifdef TIMING
void sLinsysRoot::afterFactor()
{
  int mype; MPI_Comm_rank(mpiComm, &mype);

  if( (mype/256)*256==mype) {
    for (size_t c=0; c<children.size(); c++) {
      if (children[c]->mpiComm == MPI_COMM_NULL) continue;
      
      printf("  rank %d NODE %4zu SPFACT %g BACKSOLVE %g SEC ITER %d\n", mype, c,
	     children[c]->stochNode->resMon.eFact.tmLocal,
	     children[c]->stochNode->resMon.eFact.tmChildren, (int)g_iterNumber);
    }
  }
  if( (mype/1024)*1024==mype) {
    for (size_t c=0; c<children.size(); c++) {
      if (children[c]->mpiComm == MPI_COMM_NULL) continue;
      
      double redall = stochNode->resMon.eReduce.tmLocal;
      double redscat = stochNode->resMon.eReduceScatter.tmLocal;
      printf("  rank %d REDUCE %g SEC ITER %d REDSCAT %g DIFF %g\n", mype, redall, 
	     (int)g_iterNumber, redscat, redall-redscat);
    }
  }
}
#endif

/* compute
 *              locnx locmy locmyl locmzl
 * nx_border  [   0    A0T   F0VT   G0VT ]
 * myl_border [  F0C    0     0      0   ]
 * mzl_border [  G0C    0     0      0   ]
 *
 * [  0 F0C^T  G0C^T ]^T
 * [ A0   0     0    ]
 * [ C0   0     0    ]
 * [ F0V  0     0    ]   + buffer
 * [ G0V  0     0    ]
 */
// TODO : refactor! ..
// TODO : move to aug..
void sLinsysRoot::finalizeZ0Hierarchical( DenseGenMatrix& buffer, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border )
{
   // TODO : parallelize over MPI porcs?
   finalizeDenseBorderModBlocked(Br_mod_border, buffer);

   if( !Br.has_RAC && !Br.use_local_RAC )
      return;

   bool has_RAC = Br.has_RAC;

   int mX0, nX0; buffer.getSize( mX0, nX0 );

   SparseGenMatrix* F0cons_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.F.mat) : nullptr;
   SparseGenMatrix* G0cons_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.G.mat) : nullptr;

   SparseGenMatrix* A0_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.A.mat) : nullptr;
   SparseGenMatrix* C0_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.C.mat) : nullptr;

   SparseGenMatrix* F0vec_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.A.mat_link) : &data->getLocalF();
   SparseGenMatrix* G0vec_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.C.mat_link) : &data->getLocalG();

   if( has_RAC )
      assert( F0cons_border && G0cons_border && A0_border && C0_border );
   assert( F0vec_border );
   assert( G0vec_border );

   int mA0{0}; int nA0{0};
   if( A0_border )
      A0_border->getSize(mA0, nA0);

   int mC0{0}; int nC0{0};
   if( C0_border )
      C0_border->getSize(mC0, nC0);

   int mF0C{0}; int nF0C{0};
   if( F0cons_border )
      F0cons_border->getSize( mF0C, nF0C );

   int mG0C{0}; int nG0C{0};
   if( G0cons_border )
      G0cons_border->getSize( mG0C, nG0C );

   int mF0V{0}; int nF0V{0};
   F0vec_border->getSize(mF0V, nF0V);

   int mG0V{0}; int nG0V{0};
   G0vec_border->getSize(mG0V, nG0V);

#ifndef NDEBUG
   assert( nA0 == nC0 );
   assert( nF0V == nG0V );

   if( has_RAC )
      assert( nA0 == nF0V );

   assert( nF0C == nG0C );
   if( mA0 != 0 )
      assert( nF0C + mA0 + mC0 + mF0V + mG0V == nX0 );
   if( has_RAC )
      assert( mX0 == nF0V + mF0C + mG0C );
   else
      assert( mX0 >= nF0V );
#endif

   /* add A0^T, C0^T, F0V^T, G0V^T */
   for( int row = 0; row < nF0V; ++row )
   {
      // TODO refactor and move somewhere else... TODO use addMatAt
      /* A0^T */
      if( mA0 > 0 )
      {
         const SparseGenMatrix& A0_transp = A0_border->getTranspose();

         const double* MAt = A0_transp.M();
         const int* krowAt = A0_transp.krowM();
         const int* jcolAt = A0_transp.jcolM();

         const int rowAt_start = krowAt[row];
         const int rowAt_end = krowAt[row + 1];

         const int col_buffer_start = nF0C;
         for( int k = rowAt_start; k < rowAt_end; ++k )
         {
            const int col_At = jcolAt[k];
            const double val_At = MAt[k];

            const int col_buffer = col_buffer_start + col_At;

            assert( nF0C <= col_buffer && col_buffer < nF0C + mA0 );
            buffer[row][col_buffer] += val_At;
         }
      }

      /* C0^T */
      if( mC0 > 0 )
      {
         const SparseGenMatrix& C0_transp = C0_border->getTranspose();

         const double* MCt = C0_transp.M();
         const int* krowCt = C0_transp.krowM();
         const int* jcolCt = C0_transp.jcolM();

         const int rowCt_start = krowCt[row];
         const int rowCt_end = krowCt[row + 1];

         const int col_buffer_start = nF0C + mA0;
         for( int k = rowCt_start; k < rowCt_end; ++k )
         {
            const int col_At = jcolCt[k];
            const double val_Ct = MCt[k];

            const int col_buffer = col_buffer_start + col_At;

            assert( nF0C + mA0 <= col_buffer && col_buffer < nF0C + mA0 + mC0 );
            buffer[row][col_buffer] += val_Ct;
         }
      }

      /* F0V^T */
      if( mF0V > 0 )
      {
         const SparseGenMatrix& F0V_transp = F0vec_border->getTranspose();

         const double* MF0Vt = F0V_transp.M();
         const int* krowF0Vt = F0V_transp.krowM();
         const int* jcolF0Vt = F0V_transp.jcolM();

         const int rowF0Vt_start = krowF0Vt[row];
         const int rowF0Vt_end = krowF0Vt[row + 1];

         const int col_buffer_start = nF0C + mA0 + mC0;
         for( int k = rowF0Vt_start; k < rowF0Vt_end; ++k )
         {
            const int col_F0Vt = jcolF0Vt[k];
            const double val_F0Vt = MF0Vt[k];

            const int col_buffer = col_buffer_start + col_F0Vt;

            assert( nF0C + mA0 + mC0 <= col_buffer && col_buffer < nF0C + mA0 + mC0 + mF0V );
            buffer[row][col_buffer] += val_F0Vt;
         }
      }

      /* G0V^T */
      if( mG0V > 0 )
      {
         const SparseGenMatrix& G0V_transp = G0vec_border->getTranspose();

         const double* MG0Vt = G0V_transp.M();
         const int* krowG0Vt = G0V_transp.krowM();
         const int* jcolG0Vt = G0V_transp.jcolM();

         const int rowG0Vt_start = krowG0Vt[row];
         const int rowG0Vt_end = krowG0Vt[row + 1];

         const int col_buffer_start = nF0C + mA0 + mC0 + mF0V;
         for( int k = rowG0Vt_start; k < rowG0Vt_end; ++k )
         {
            const int col_G0Vt = jcolG0Vt[k];
            const double val_G0Vt = MG0Vt[k];

            const int col_buffer = col_buffer_start + col_G0Vt;

            assert( nF0C + mA0 + mC0 + mF0V <= col_buffer && col_buffer < nF0C + mA0 + mC0 + mF0V + mG0V );
            buffer[row][col_buffer] += val_G0Vt;
         }
      }
   }

   /* F0C */
   if( mF0C > 0 )
   {
      for( int rowF = 0; rowF < mF0C; ++rowF )
      {
         const double* MF0C = F0cons_border->M();
         const int* krowF0C = F0cons_border->krowM();
         const int* jcolF0C = F0cons_border->jcolM();

         const int rowF0C_start = krowF0C[rowF];
         const int rowF0C_end = krowF0C[rowF + 1];

         for( int k = rowF0C_start; k < rowF0C_end; ++k )
         {
            const int col = jcolF0C[k];
            const double val_F0C = MF0C[k];

            const int row_buffer = rowF + nA0;

            assert( nA0 <= row_buffer && row_buffer < nA0 + mF0C );
            assert( 0 <= col && col < locnx );
            buffer[row_buffer][col] += val_F0C;
         }
      }
   }

   /* G0C */
   if( mG0C > 0 )
   {
      for( int rowG = 0; rowG < mG0C; ++rowG )
      {
         const double* MG0C = G0cons_border->M();
         const int* krowG0C = G0cons_border->krowM();
         const int* jcolG0C = G0cons_border->jcolM();

         const int rowG0C_start = krowG0C[rowG];
         const int rowG0C_end = krowG0C[rowG + 1];

         for( int k = rowG0C_start; k < rowG0C_end; ++k )
         {
            const int col = jcolG0C[k];
            const double val_G0C = MG0C[k];

            const int row_buffer = rowG + nA0 + mF0C;

            assert( nA0 + mF0C <= row_buffer && row_buffer < nA0 + mF0C + mG0C );
            assert( 0 <= col && col < locnx );
            buffer[row_buffer][col] += val_G0C;
         }
      }
   }
}

/* compute SC -= Br0^T X0
 *         [  0  A0T C0T F0VT G0VT ]
 * Br0^T = [ F0C  0   0   0    0   ]
 *         [ G0C  0   0   0    0   ]
 *
 * SC is still stored in transposed form as well as X0
 *
 * SC -= X0 Br0 instead
 *
 * Br0 = [  0  F0CT G0CT ]
 *       [  A   0    0   ]
 *       [  C   0    0   ]
 *       [ F0V  0    0   ]
 *       [ G0V  0    0   ]
 *
 */
void sLinsysRoot::finalizeInnerSchurComplementContribution( DoubleMatrix& SC_, DenseGenMatrix& X0, BorderLinsys& Br, bool is_sym, bool is_sparse )
{
   if( is_sparse )
      assert( is_sym );
   const bool has_RAC = Br.has_RAC;

   if( !has_RAC && !Br.use_local_RAC )
      return;

   int mX0, nX0; X0.getSize(mX0, nX0);
   int mSC, nSC; SC_.getSize(mSC, nSC);

   assert( mSC == mX0 );

   SparseGenMatrix* F0cons_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.F.mat) : nullptr;
   SparseGenMatrix* G0cons_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.G.mat) : nullptr;

   SparseGenMatrix* A0_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.A.mat) : nullptr;
   SparseGenMatrix* C0_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.C.mat) : nullptr;
   SparseGenMatrix* F0vec_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.A.mat_link) : &data->getLocalF();
   SparseGenMatrix* G0vec_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(Br.C.mat_link) : &data->getLocalG();

   assert( F0vec_border );
   assert( G0vec_border );

   int mF0V{0}; int nF0V{0};
   F0vec_border->getSize(mF0V, nF0V);

   int mG0V{0}; int nG0V{0};
   G0vec_border->getSize(mG0V, nG0V);

   if( !has_RAC && nF0V == 0 && nG0V == 0 )
      return;

#ifndef NDEBUG
   int mF0C{0}; int nF0C{0};
   if( F0cons_border )
      F0cons_border->getSize( mF0C, nF0C );

   int mG0C{0}; int nG0C{0};
   if( G0cons_border )
      G0cons_border->getSize( mG0C, nG0C );

   int mA0{0}; int nA0{0};
   if( A0_border )
      A0_border->getSize(mA0, nA0);

   int mC0{0}; int nC0{0};
   if( C0_border )
      C0_border->getSize(mC0, nC0);


   assert( nA0 == nC0 );
   assert( nF0V == nG0V );

   if( has_RAC )
      assert( nA0 == nF0V );

   assert( nF0C + mA0 + mC0 + mF0V + mG0V == nX0 );
   assert( nF0C == nG0C );

   if( has_RAC )
      assert( nSC == nF0V + mF0C + mG0C );
   else
      assert( nSC >= mF0V );

   assert( mX0 == mSC );
#endif

   if( is_sparse )
      finalizeInnerSchurComplementContributionSparse(SC_, X0, A0_border, C0_border, F0vec_border, G0vec_border, F0cons_border, G0cons_border );
   else
      finalizeInnerSchurComplementContributionDense(SC_, X0, A0_border, C0_border, F0vec_border, G0vec_border, F0cons_border, G0cons_border, is_sym );
}

/* SC and X0 stored in transposed form */
void sLinsysRoot::finalizeInnerSchurComplementContributionSparse( DoubleMatrix& SC_, DenseGenMatrix& X0, SparseGenMatrix* A0_border,
      SparseGenMatrix* C0_border, SparseGenMatrix* F0vec_border, SparseGenMatrix* G0vec_border, SparseGenMatrix* F0cons_border, SparseGenMatrix* G0cons_border )
{
   assert( F0vec_border );
   assert( G0vec_border );

   SparseSymMatrix& SC = dynamic_cast<SparseSymMatrix&>(SC_);

   int dummy, mX0; X0.getSize(mX0, dummy);

   int mA0{0}; int nA0{0};
   if( A0_border )
      A0_border->getSize(mA0, nA0);

   int mC0{0};
   if( C0_border )
      C0_border->getSize(mC0, dummy);

   int nF0C{0}; int mF0C{0};
   if( F0cons_border )
      F0cons_border->getSize( mF0C, nF0C );

   int mF0V{0};
   F0vec_border->getSize(mF0V, dummy);

   int mG0V{0};
   G0vec_border->getSize(mG0V, dummy);

   // multiply each column with B_{outer]}^T and add it to res
   // todo: #pragma omp parallel for schedule(dynamic, 10)
   for( int i = 0; i < mX0; i++ )
   {
      const double* const col = X0[i];
      if( A0_border )
         A0_border->transmultMatSymUpper(1.0, SC, -1.0, &col[nF0C], i, 0 );

      if( C0_border )
         C0_border->transmultMatSymUpper(1.0, SC, -1.0, &col[nF0C + mA0], i, 0 );

      if( mF0V > 0 )
         F0vec_border->transmultMatSymUpper(1.0, SC, -1.0, &col[nF0C + mA0 + mC0], i, 0 );

      if( mG0V > 0 )
         G0vec_border->transmultMatSymUpper(1.0, SC, -1.0, &col[nF0C + mA0 + mC0 + mF0V], i, 0 );

      if( F0cons_border )
         F0cons_border->multMatSymUpper(1.0, SC, -1.0, &col[0], i, nA0);

      if( G0cons_border )
         G0cons_border->multMatSymUpper(1.0, SC, -1.0, &col[0], i, nA0 + mF0C );
   }
}

/* SC and X0 stored in transposed form */
void sLinsysRoot::finalizeInnerSchurComplementContributionDense( DoubleMatrix& SC_, DenseGenMatrix& X0, SparseGenMatrix* A0_border,
      SparseGenMatrix* C0_border, SparseGenMatrix* F0vec_border, SparseGenMatrix* G0vec_border, SparseGenMatrix* F0cons_border, SparseGenMatrix* G0cons_border,
      bool is_sym )
{
   assert( F0vec_border );
   assert( G0vec_border );

   double** SC = is_sym ? dynamic_cast<DenseSymMatrix&>(SC_).Mat() : dynamic_cast<DenseGenMatrix&>(SC_).Mat();

   int dummy, mX0; X0.getSize(mX0, dummy);

   int mA0{0}; int nA0{0};
   if( A0_border )
      A0_border->getSize(mA0, nA0);

   int mC0{0};
   if( C0_border )
      C0_border->getSize(mC0, dummy);

   int nF0C{0}; int mF0C{0};
   if( F0cons_border )
      F0cons_border->getSize( mF0C, nF0C );

   int mF0V{0};
   F0vec_border->getSize(mF0V, dummy);

   // multiply each column with B_{outer]}^T and add it to res
   // todo: #pragma omp parallel for schedule(dynamic, 10)
   for( int i = 0; i < mX0; i++ )
   {
      const double* const col = X0[i];
      if( A0_border )
         A0_border->transMult(1.0, &SC[i][0], 1, -1.0, &col[nF0C], 1);

      if( C0_border )
         C0_border->transMult(1.0, &SC[i][0], 1, -1.0, &col[nF0C + mA0], 1);

      F0vec_border->transMult(1.0, &SC[i][0], 1, -1.0, &col[nF0C + mA0 + mC0], 1);

      G0vec_border->transMult(1.0, &SC[i][0], 1, -1.0, &col[nF0C + mA0 + mC0 + mF0V], 1);

      if( F0cons_border )
         F0cons_border->mult(1.0, &SC[i][nA0], 1, -1.0, &col[0], 1);

      if( G0cons_border )
         G0cons_border->mult(1.0, &SC[i][nA0 + mF0C], 1, -1.0, &col[0], 1);
   }
}

/* compute -SUM_i Bi_{inner}^T Ki^{-1} Bri */
void sLinsysRoot::LsolveHierarchyBorder( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool use_local_RAC )
{
   assert( children.size() == Br.F.children.size() );

   /* get contribution to schur_complement from each child */
   for( size_t it = 0; it < children.size(); it++ )
   {
      std::vector<BorderMod> Br_mod_border_child;

      for( auto& br_mod : Br_mod_border )
         Br_mod_border_child.push_back( getChild( br_mod, it ) );

      BorderLinsys border_child = getChild( Br, it );
      children[it]->addInnerBorderKiInvBrToRes(result, border_child, Br_mod_border_child, use_local_RAC);
   }

   /* allreduce the result */
   // TODO : optimize -> do not reduce A_0 part ( all zeros... )
   if( iAmDistrib )
      allreduceMatrix( result, false, false, mpiComm );
}

/* compute SUM_i Bli^T X_i = SUM_i Bli^T Ki^-1 (( Bri - sum_j Bmodij Xij ) - Bi_{inner} X0) */
void sLinsysRoot::LtsolveHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
      bool sym_res, bool sparse_res, bool use_local_RAC )
{
   assert( !is_hierarchy_root );

   /* X0 is still in transposed form */
   assert( children.size() == Br.F.children.size() );
   assert( children.size() == Bl.F.children.size() );
   assert( !Br.F.isKindOf( kStringGenDummyMatrix ) );
   assert( !Bl.F.isKindOf( kStringGenDummyMatrix ) );

   /* for every child - add Bi_{outer}^T Ki^-1 (Bi_{outer} - Bi_{inner} X0) */
   for( size_t it = 0; it < children.size(); it++ )
   {
      BorderLinsys bl_child = getChild( Bl, it );
      BorderLinsys br_child = getChild( Br, it );

      std::vector<BorderMod> border_mod_child;
      for( auto& bm : Br_mod_border )
         border_mod_child.push_back( getChild( bm, it ) );

      children[it]->LniTransMultHierarchyBorder( res, X0, bl_child, br_child, border_mod_child, sparse_res, sym_res, use_local_RAC );
   }
}

void sLinsysRoot::addBorderX0ToRhs( StochVector& rhs, const SimpleVector& x0, BorderLinsys& border )
{
   assert( rhs.children.size() == children.size() );
   assert( border.A.children.size() == children.size() );

   for( size_t i = 0; i < children.size(); ++i )
   {
      BorderLinsys child_border = getChild(border, i);
      children[i]->addBorderX0ToRhs( *rhs.children[i], x0, child_border );
   }

   /* add schur complement part */
   assert( border.A.mat );
   assert( border.C.mat );

   SparseGenMatrix& A0_border = dynamic_cast<SparseGenMatrix&>(*border.A.mat);
   int mA0, nA0; A0_border.getSize(mA0, nA0);

   SparseGenMatrix& C0_border = dynamic_cast<SparseGenMatrix&>(*border.C.mat);
   int mC0, nC0; C0_border.getSize(mC0, nC0);

   assert( border.F.mat );
   assert( border.A.mat_link );
   SparseGenMatrix& F0vec_border = dynamic_cast<SparseGenMatrix&>(*border.A.mat_link);
   int mF0V, nF0V; F0vec_border.getSize(mF0V, nF0V);
   SparseGenMatrix& F0cons_border = dynamic_cast<SparseGenMatrix&>(*border.F.mat);
   int mF0C, nF0C; F0cons_border.getSize(mF0C, nF0C);

   assert( border.C.mat_link );
   assert( border.G.mat );
   SparseGenMatrix& G0vec_border = dynamic_cast<SparseGenMatrix&>(*border.C.mat_link);
   int mG0V, nG0V; G0vec_border.getSize(mG0V, nG0V);
   SparseGenMatrix& G0cons_border = dynamic_cast<SparseGenMatrix&>(*border.G.mat);
   int mG0C, nG0C; G0cons_border.getSize(mG0C, nG0C);

   assert( rhs.vec );
   assert( rhs.vec->length() == nF0C + mA0 + mC0 + mF0V + mG0V );
   assert( x0.length() == nA0 + mF0C + mG0C );

   SimpleVector& rhs0 = dynamic_cast<SimpleVector&>(*rhs.vec);

   double* rhs01 = &rhs0[0];
   double* rhs02 = &rhs0[nF0C];
   double* rhs03 = &rhs0[nF0C + mA0];
   double* rhs04 = &rhs0[nF0C + mA0 + mC0];
   double* rhs05 = &rhs0[nF0C + mA0 + mC0 + mF0V];

   const double* x01 = &x0[0];
   const double* x02 = &x0[nA0];
   const double* x03 = &x0[nA0 + mF0C];

   A0_border.mult(1.0, rhs02, 1, -1.0, x01, 1);
   C0_border.mult(1.0, rhs03, 1, -1.0, x01, 1);
   F0vec_border.mult(1.0, rhs04, 1, -1.0, x01, 1);
   G0vec_border.mult(1.0, rhs05, 1, -1.0, x01, 1);

   F0cons_border.transMult(1.0, rhs01, 1, -1.0, x02, 1 );
   G0cons_border.transMult(1.0, rhs01, 1, -1.0, x03, 1 );
}

void sLinsysRoot::addBorderTimesRhsToB0( StochVector& rhs, SimpleVector& b0, BorderLinsys& border )
{
   assert( rhs.children.size() == children.size() );
   assert( border.A.children.size() == children.size() );

   for( size_t i = 0; i < children.size(); ++i )
   {
      BorderLinsys child_border = getChild( border, i );
      children[i]->addBorderTimesRhsToB0( *rhs.children[i], b0, child_border );
   }

   /* add schur complement part */
   if( PIPS_MPIgetSize(mpiComm) == 0 || PIPS_MPIgetRank(mpiComm) == 0 )
   {
      assert( border.A.mat );
      assert( border.C.mat );

      SparseGenMatrix& A0_border = dynamic_cast<SparseGenMatrix&>(*border.A.mat);
      int mA0, nA0; A0_border.getSize(mA0, nA0);

      SparseGenMatrix& C0_border = dynamic_cast<SparseGenMatrix&>(*border.C.mat);
      int mC0, nC0; C0_border.getSize(mC0, nC0);

      assert( border.F.mat );
      assert( border.A.mat_link );
      SparseGenMatrix& F0vec_border = dynamic_cast<SparseGenMatrix&>(*border.A.mat_link);
      int mF0V, nF0V; F0vec_border.getSize(mF0V, nF0V);
      SparseGenMatrix& F0cons_border = dynamic_cast<SparseGenMatrix&>(*border.F.mat);
      int mF0C, nF0C; F0cons_border.getSize(mF0C, nF0C);

      assert( border.C.mat_link );
      assert( border.G.mat );
      SparseGenMatrix& G0vec_border = dynamic_cast<SparseGenMatrix&>(*border.C.mat_link);
      int mG0V, nG0V; G0vec_border.getSize(mG0V, nG0V);
      SparseGenMatrix& G0cons_border = dynamic_cast<SparseGenMatrix&>(*border.G.mat);
      int mG0C, nG0C; G0cons_border.getSize(mG0C, nG0C);

      assert( rhs.vec );
      assert( rhs.vec->length() == nF0C + mA0 + mC0 + mF0V + mG0V );
      assert( b0.length() == nA0 + mF0C + mG0C );

      SimpleVector& zi = dynamic_cast<SimpleVector&>(*rhs.vec);

      SimpleVector zi1 (&zi[0], nF0C);
      SimpleVector zi2 (&zi[nF0C], mA0 );
      SimpleVector zi3 (&zi[nF0C + mA0], mC0);
      SimpleVector zi4 (&zi[nF0C + mA0 + mC0], mF0V);
      SimpleVector zi5 (&zi[nF0C + mA0 + mC0 + mF0V], mG0V);

      SimpleVector b1( &b0[0], nA0 );
      SimpleVector b2( &b0[nA0], mF0C );
      SimpleVector b3( &b0[nA0 + mF0C], mG0C );

      A0_border.transMult(1.0, b1, -1.0, zi2);
      C0_border.transMult(1.0, b1, -1.0, zi3);
      F0vec_border.transMult(1.0, b1, -1.0, zi4);
      G0vec_border.transMult(1.0, b1, -1.0, zi5);

      F0cons_border.mult(1.0, b2, -1.0, zi1);
      G0cons_border.mult(1.0, b3, -1.0, zi1);
   }
}

void sLinsysRoot::Ltsolve2( sData*, StochVector& x, SimpleVector& x0, bool)
{
   assert( false && "not in use");
   assert( pips_options::getBoolParameter("HIERARCHICAL") );
   assert( children.size() == x.children.size() );

   StochVector& b = dynamic_cast<StochVector&>(x);

   for( size_t i = 0; i < children.size(); ++i )
   {
      children[i]->computeInnerSystemRightHandSide( *b.children[i], x0, true );
      children[i]->solveCompressed( *x.children[i] );
   }

//#ifdef TIMING
//  stochNode->resMon.eLtsolve.clear();
//  stochNode->resMon.recLtsolveTmLocal_start();
//#endif
//   //b_i -= Lni^T x0
//   LniTransMult(prob, bi, -1.0, x0);
//   solver->Ltsolve(bi);
//
//#ifdef TIMING
//  stochNode->resMon.recLtsolveTmLocal_stop();
//#endif
//   SimpleVector& xi = bi;
//   //recursive call in order to get the children to do their part
//   for(size_t it = 0; it < children.size(); it++)
//      children[it]->Ltsolve2(prob->children[it], *b.children[it], xi);
}

void sLinsysRoot::createChildren(sData *prob)
{
   sLinsys* child{};
   assert( dd && dq && nomegaInv && rhs && prob );

   StochVector &ddst = dynamic_cast<StochVector&>(*dd);
   StochVector &dqst = dynamic_cast<StochVector&>(*dq);
   StochVector &nomegaInvst = dynamic_cast<StochVector&>(*nomegaInv);
   StochVector &rhsst = dynamic_cast<StochVector&>(*rhs);

   for( size_t it = 0; it < prob->children.size(); it++ )
   {
      assert( ddst.children[it] != nullptr );
      if( MPI_COMM_NULL == ddst.children[it]->mpiComm )
      {
         child = new sDummyLinsys(dynamic_cast<sFactory*>(factory),
               prob->children[it]);
      }
      else
      {
         assert( prob->children[it] );
         sFactory *stochFactory = dynamic_cast<sFactory*>(factory);
         if( is_hierarchy_root )
         {
            assert( prob->isHierarchyRoot() );
            assert( prob->children.size() == 1 );
            assert( prob->children[0] );
            assert( ddst.children.size() == 1 && dqst.children.size() == 1 && nomegaInvst.children.size() == 1
                  && rhsst.children.size() == 1 );
            assert( MPI_COMM_NULL != ddst.children[0]->mpiComm);
         }

         if( prob->children[it]->children.empty() )
         {
            child = stochFactory->newLinsysLeaf(prob->children[it],
                  ddst.children[it], dqst.children[it],
                  nomegaInvst.children[it], rhsst.children[it]);
         }
         else
         {
            assert( prob->children[it] );
            child = stochFactory->newLinsysRoot(prob->children[it],
                  ddst.children[it], dqst.children[it],
                  nomegaInvst.children[it], rhsst.children[it]);
         }
      }
      assert( child != nullptr );
      AddChild(child);
   }
}

void sLinsysRoot::deleteChildren()
{
  for(size_t it = 0; it < children.size(); it++) {
    children[it]->deleteChildren();
    delete children[it];
  }
  children.clear();
}

void sLinsysRoot::initProperChildrenRange()
{
   assert(children.size() > 0);

   int childStart = -1;
   int childEnd = -1;
   for( size_t it = 0; it < children.size(); it++ )
   {
      if( childEnd != -1 )
         assert(children[it]->isDummy());

      if( children[it]->isDummy() )
      {
         // end of range?
         if( childStart != -1 && childEnd == -1 )
            childEnd = int(it);

         continue;
      }

      // start of range?
      if( childStart == -1 )
         childStart = int(it);
   }

   assert(childStart >= 0);

   if( childEnd == -1 )
   {
      assert(!children[children.size() - 1]->isDummy());
      childEnd = int(children.size());
   }

    assert(childStart < childEnd && childEnd <= int(children.size()));

    childrenProperStart = childStart;
    childrenProperEnd = childEnd;
}

void sLinsysRoot::putXDiagonal( OoqpVector& xdiag_ )
{
  StochVector& xdiag = dynamic_cast<StochVector&>(xdiag_);
  assert(children.size() == xdiag.children.size());

  //kkt->atPutDiagonal( 0, *xdiag.vec );
  xDiag = xdiag.vec;
 
  // propagate it to the subtree
  for(size_t it = 0; it < children.size(); it++)
    children[it]->putXDiagonal(*xdiag.children[it]);
}


void sLinsysRoot::putZDiagonal( OoqpVector& zdiag_ )
{
  StochVector& zdiag = dynamic_cast<StochVector&>(zdiag_);
  assert(children.size() == zdiag.children.size());

  //kkt->atPutDiagonal( locnx+locmy, *zdiag.vec );
  zDiag = zdiag.vec;
  zDiagLinkCons = zdiag.vecl;

  // propagate it to the subtree
  for(size_t it = 0; it < children.size(); it++)
    children[it]->putZDiagonal(*zdiag.children[it]);
}

void sLinsysRoot::AddChild(sLinsys* child)
{
  children.push_back(child);
}

///////////////////////////////////////////////////////////
// ATOMS of FACTOR 2
//////////////////////////////////////////////////////////
/* Atoms methods of FACTOR2 for a non-leaf linear system */
void sLinsysRoot::initializeKKT(sData*, Variables*)
{
   if( hasSparseKkt )
      dynamic_cast<SparseSymMatrix*>(kkt.get())->symPutZeroes();
   else
   {
      DenseSymMatrix* kktd = dynamic_cast<DenseSymMatrix*>(kkt.get());
      myAtPutZeros(kktd);
   }
}

void sLinsysRoot::reduceKKT(sData* prob)
{
   if( usePrecondDist )
      reduceKKTdist(prob);
   else if( hasSparseKkt )
      reduceKKTsparse();
   else
      reduceKKTdense();
}


// collects (reduces) dense global Schur complement
void sLinsysRoot::reduceKKTdense()
{
   DenseSymMatrix* const kktd = dynamic_cast<DenseSymMatrix*>(kkt.get());

   // parallel communication
   if( iAmDistrib )
   {
      if( locnx > 0 )
         submatrixAllReduceDiagLower(kktd, 0, locnx, mpiComm);

      if( locmyl > 0 || locmzl > 0 )
      {
         const int locNxMy = locnx + locmy;
         assert(kktd->size() == locnx + locmy + locmyl + locmzl);

         // reduce lower left part
         if( locnx > 0 )
            submatrixAllReduceFull(kktd, locNxMy, 0, locmyl + locmzl, locnx, mpiComm);

         // reduce lower diagonal linking part
         submatrixAllReduceDiagLower(kktd, locNxMy, locmyl + locmzl, mpiComm);
      }
  }
}


// collects sparse global Schur complement
void sLinsysRoot::reduceKKTsparse()
{
   if( !iAmDistrib )
      return;
   assert(kkt);

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();
   const int sizeKkt = locnx + locmy + locmyl + locmzl;
   const int nnzKkt = krowKkt[sizeKkt];

   assert(kkts.size() == sizeKkt);
   assert(!kkts.isLower);

   if( allreduce_kkt )
      reduceToAllProcs(nnzKkt, MKkt);
   else
      reduceToProc0(nnzKkt, MKkt);
}

#define CHUNK_SIZE (1024*1024*64) //doubles = 128 MBytes (maximum)
void sLinsysRoot::reduceToAllProcs(int size, double* values)
{
   assert(values && values != sparseKktBuffer);
   assert(size > 0);

   if( sparseKktBuffer == nullptr )
      sparseKktBuffer = new double[CHUNK_SIZE];

   const int reps = size / CHUNK_SIZE;
   const int res = size - CHUNK_SIZE * reps;
   assert(res >= 0 && res < CHUNK_SIZE);

   for( int i = 0; i < reps; i++ )
   {
      double* const start = &values[i * CHUNK_SIZE];
      MPI_Allreduce(start, sparseKktBuffer, CHUNK_SIZE, MPI_DOUBLE, MPI_SUM, mpiComm);

      memcpy(start, sparseKktBuffer, size_t(CHUNK_SIZE) * sizeof(double));
   }

   if( res > 0 )
   {
      double* const start = &values[reps * CHUNK_SIZE];
      MPI_Allreduce(start, sparseKktBuffer, res, MPI_DOUBLE, MPI_SUM, mpiComm);

         memcpy(start, sparseKktBuffer, size_t(res) * sizeof(double));
   }
}

#define CHUNK_SIZE (1024*1024*64) //doubles = 128 MBytes (maximum)
void sLinsysRoot::reduceToProc0(int size, double* values)
{
   assert(values && values != sparseKktBuffer);
   assert(size > 0);

   int myRank; MPI_Comm_rank(mpiComm, &myRank);

   if( myRank == 0 && sparseKktBuffer == nullptr )
      sparseKktBuffer = new double[CHUNK_SIZE];

   const int reps = size / CHUNK_SIZE;
   const int res = size - CHUNK_SIZE * reps;
   assert(res >= 0 && res < CHUNK_SIZE);

   for( int i = 0; i < reps; i++ )
   {
      double* const start = &values[i * CHUNK_SIZE];
      MPI_Reduce(start, sparseKktBuffer, CHUNK_SIZE, MPI_DOUBLE, MPI_SUM, 0, mpiComm);

      if( myRank == 0 )
         memcpy(start, sparseKktBuffer, size_t(CHUNK_SIZE) * sizeof(double));
   }

   if( res > 0 )
   {
      double* const start = &values[reps * CHUNK_SIZE];
      MPI_Reduce(start, sparseKktBuffer, res, MPI_DOUBLE, MPI_SUM, 0, mpiComm);

      if( myRank == 0 )
         memcpy(start, sparseKktBuffer, size_t(res) * sizeof(double));
   }
}


void sLinsysRoot::registerMatrixEntryTripletMPI()
{
   assert(MatrixEntryTriplet_mpi == MPI_DATATYPE_NULL);

   const int nitems = 3;
   int blocklengths[3] = { 1, 1, 1 };
   MPI_Datatype Types[3] = { MPI_DOUBLE, MPI_INT, MPI_INT };
   MPI_Aint offsets[3];

   offsets[0] = offsetof(MatrixEntryTriplet, val);
   offsets[1] = offsetof(MatrixEntryTriplet, row);
   offsets[2] = offsetof(MatrixEntryTriplet, col);

   MPI_Type_create_struct(nitems, blocklengths, offsets, Types, &MatrixEntryTriplet_mpi);
   MPI_Type_commit(&MatrixEntryTriplet_mpi);
}

void sLinsysRoot::syncKKTdistLocalEntries(sData* prob)
{
   if( !iAmDistrib )
      return;

   assert(kkt && hasSparseKkt);

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

   int* const krowKkt = kkts.krowM();
   int* const jColKkt = kkts.jcolM();
   double* const MKkt = kkts.M();

   const int childStart = childrenProperStart;
   const int childEnd = childrenProperEnd;
   int myRank; MPI_Comm_rank(mpiComm, &myRank);
   int size; MPI_Comm_size(mpiComm, &size);

   assert(size > 1);

   // MPI matrix entries triplet not registered yet?
   if( MatrixEntryTriplet_mpi == MPI_DATATYPE_NULL )
      registerMatrixEntryTripletMPI();

   // pack the entries that will be send below
   std::vector<MatrixEntryTriplet> prevEntries = this->packKKTdistOutOfRangeEntries(prob, childStart, childEnd);
   std::vector<MatrixEntryTriplet> myEntries(0);

   assert(prevEntries.size() > 0 && prevEntries[0].row == - 1 && prevEntries[0].col == - 1);

   // odd processes send first (one process back)
   if( myRank % 2 != 0 )
   {
      this->sendKKTdistLocalEntries(prevEntries);
   }

   // even processes (except last) receive first
   if( myRank % 2 == 0 && myRank != size - 1 )
   {
      assert(myEntries.size() == 0);
      myEntries = this->receiveKKTdistLocalEntries();
   }

   // even processes (except first) send
   if( myRank % 2 == 0 && myRank > 0 )
   {
      this->sendKKTdistLocalEntries(prevEntries);
   }

   // odd processes (except last) receive
   if( myRank % 2 != 0 && myRank != size - 1 )
   {
      assert(myEntries.size() == 0);
      myEntries = this->receiveKKTdistLocalEntries();
   }

   assert(myEntries.size() > 0 || myRank == size - 1 );

   int lastRow = 0;
   int lastC = -1;

#ifndef NDEBUG
   const std::vector<bool>& rowIsLocal = prob->getSCrowMarkerLocal();
   const std::vector<bool>& rowIsMyLocal = prob->getSCrowMarkerMyLocal();
#endif

   // finally, put received data into Schur complement matrix
   for( size_t i = 1; i < myEntries.size(); i++ )
   {
      const double val = myEntries[i].val;
      const int row = myEntries[i].row;
      const int col = myEntries[i].col;

      assert(myRank != size - 1);
      assert(row >= 0 && row < locnx + locmy + locmyl + locmzl);
      assert(col >= row && col < locnx + locmy + locmyl + locmzl);
      assert(val == val); // catch NaNs

      assert(rowIsMyLocal[row] || (rowIsMyLocal[col] && !rowIsLocal[row]));

      int c;

      // continue from last position?
      if( row == lastRow )
      {
         for( c = lastC + 1; c < krowKkt[row + 1]; c++ )
         {
            const int colKkt = jColKkt[c];

            if( colKkt == col )
            {
               MKkt[c] += val;
               break;
            }
         }

         // found the correct entry in last row?
         if( c != krowKkt[row + 1] )
         {
            assert(c < krowKkt[row + 1]);
            lastRow = row;
            lastC = c;
            continue;
         }
      }

      c = krowKkt[row];
      assert(col >= jColKkt[c]);

      for( ; c < krowKkt[row + 1]; c++ )
      {
         const int colKkt = jColKkt[c];

         if( colKkt == col )
         {
            MKkt[c] += val;
            break;
         }
      }

      assert(c != krowKkt[row + 1]);

      lastRow = row;
      lastC = c;
   }
}


std::vector<sLinsysRoot::MatrixEntryTriplet> sLinsysRoot::receiveKKTdistLocalEntries() const
{
   assert(kkt && hasSparseKkt);
   assert(MatrixEntryTriplet_mpi != MPI_DATATYPE_NULL);

   int myRank; MPI_Comm_rank(mpiComm, &myRank);
   int size; MPI_Comm_size(mpiComm, &size);
   const int nextRank = myRank + 1;
   assert(nextRank < size);

   // receive data from next process

   MPI_Status status;
   MPI_Probe(nextRank, 0, mpiComm, &status);

   int nEntries;
   MPI_Get_count(&status, MatrixEntryTriplet_mpi, &nEntries);

   assert(nEntries >= 1);

   std::vector<MatrixEntryTriplet> entries(nEntries);

   MPI_Recv((void*) &entries[0], nEntries, MatrixEntryTriplet_mpi, nextRank, 0, mpiComm, MPI_STATUS_IGNORE);

   assert(entries[0].row == -1 && entries[0].col == -1); // dummy check

   PIPSdebugMessage("myRank=%d received %d \n", myRank, nEntries);

   return entries;
}


void sLinsysRoot::sendKKTdistLocalEntries(const std::vector<MatrixEntryTriplet>& prevEntries) const
{
   int myRank; MPI_Comm_rank(mpiComm, &myRank);
   const int prevRank = myRank - 1;
   const int nEntries = int(prevEntries.size());

   assert(myRank >= 0);
   assert(nEntries > 0);
   assert(MatrixEntryTriplet_mpi != MPI_DATATYPE_NULL);

   PIPSdebugMessage("myRank=%d sends %d \n", myRank, nEntries);
   MPI_Send(&prevEntries[0], nEntries, MatrixEntryTriplet_mpi, prevRank, 0, mpiComm);
}

std::vector<sLinsysRoot::MatrixEntryTriplet> sLinsysRoot::packKKTdistOutOfRangeEntries(sData* prob, int childStart, int) const
{
   assert(kkt && hasSparseKkt);

   int myRank; MPI_Comm_rank(mpiComm, &myRank);

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);
   const std::vector<bool>& rowIsLocal = prob->getSCrowMarkerLocal();
   const std::vector<bool>& rowIsMyLocal = prob->getSCrowMarkerMyLocal();
   int* const krowKkt = kkts.krowM();
   int* const jColKkt = kkts.jcolM();
   double* const MKkt = kkts.M();
   const int sizeKkt = locnx + locmy + locmyl + locmzl;

   std::vector<MatrixEntryTriplet> packedEntries(0);

   // add dummy value
   const MatrixEntryTriplet entry_zero = {-1.0, -1, -1};
   packedEntries.push_back(entry_zero);

   if( childStart > 0 )
   {
      assert(myRank > 0);

      // pack data
      for( int r = 0; r < sizeKkt; r++ )
      {
         const bool rIsLocal = rowIsLocal[r];

         if( rIsLocal && !rowIsMyLocal[r] )
         {
            for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
            {
               const int col = jColKkt[c];

               if( !rowIsMyLocal[col] )
               {
                  const double val = MKkt[c];

                  if( PIPSisZero(val) )
                     continue;

                  const MatrixEntryTriplet entry = {val, r, col};
                  packedEntries.push_back(entry);
               }
            }
         }

         if( rIsLocal )
            continue;

         for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
         {
            const int col = jColKkt[c];

            if( rowIsLocal[col] && !rowIsMyLocal[col] )
            {
               const double val = MKkt[c];

               if( PIPSisZero(val) )
                  continue;

               const MatrixEntryTriplet entry = {val, r, col};
               packedEntries.push_back(entry);
            }
         }
      }
   }

   return packedEntries;
}


void sLinsysRoot::reduceKKTdist(sData* prob)
{
   assert(prob);
   assert(iAmDistrib);
   assert(kkt);

   const std::vector<bool>& rowIsLocal = prob->getSCrowMarkerLocal();
   const std::vector<bool>& rowIsMyLocal = prob->getSCrowMarkerMyLocal();

   SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

   int* const krowKkt = kkts.krowM();
   int* const jColKkt = kkts.jcolM();
   double* const MKkt = kkts.M();
   const int sizeKkt = locnx + locmy + locmyl + locmzl;
   int nnzDistMyLocal = 0;
   int nnzDistShared = 0;
   int nnzDistLocal;

   std::vector<int> rowSizeMyLocal(sizeKkt, 0);
   std::vector<int> rowSizeShared(sizeKkt, 0);
   std::vector<int> rowSizeLocal(sizeKkt, 0);
   std::vector<int> rowIndexMyLocal(0);
   std::vector<int> colIndexMyLocal(0);

   assert(int(rowIsLocal.size()) == sizeKkt);

   // add up locally owned entries
   this->syncKKTdistLocalEntries(prob);

   // add B_0, F_0, G_0 and diagonals (all scattered)
   this->finalizeKKTdist(prob);

   precondSC.updateDiagDomBound();
   precondSC.unmarkDominatedSCdistLocals(*prob, kkts);

   // compute row lengths
   for( int r = 0; r < sizeKkt; r++ )
   {
      if( rowIsMyLocal[r] )
      {
         for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
         {
            const int col = jColKkt[c];

            if( col < 0 )
               continue;

            nnzDistMyLocal++;
            rowSizeMyLocal[r]++;
            rowIndexMyLocal.push_back(r);
            colIndexMyLocal.push_back(col);
         }

         continue;
      }

      const bool rIsLocal = rowIsLocal[r];

      for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
      {
         const int col = jColKkt[c];

         if( col < 0 )
         {
            assert(-col - 1 >= 0 && -col - 1 < sizeKkt);
            assert((!(!rIsLocal && !rowIsLocal[-col - 1])));
            continue;
         }

         // is (r, col) a shared entry?
         if( !rIsLocal && !rowIsLocal[col] )
         {
            nnzDistShared++;
            rowSizeShared[r]++;
            assert(!rowIsMyLocal[col]);
         }

         if( rowIsMyLocal[col] )
         {
            nnzDistMyLocal++;
            rowSizeMyLocal[r]++;
            rowIndexMyLocal.push_back(r);
            colIndexMyLocal.push_back(col);
         }
      }
   }

   assert(int(rowIndexMyLocal.size()) == nnzDistMyLocal);

   // sum up local sizes
   MPI_Allreduce(&nnzDistMyLocal, &nnzDistLocal, 1, MPI_INT, MPI_SUM, mpiComm);
   MPI_Allreduce(&rowSizeMyLocal[0], &rowSizeLocal[0], sizeKkt, MPI_INT, MPI_SUM, mpiComm);

#ifndef NDEBUG
   {
      int nnzDistSharedMax;
      std::vector<int> rowSizeSharedMax(sizeKkt, 0);

      MPI_Allreduce(&nnzDistShared, &nnzDistSharedMax, 1, MPI_INT, MPI_MAX, mpiComm);
      MPI_Allreduce(&rowSizeShared[0], &rowSizeSharedMax[0], sizeKkt, MPI_INT, MPI_MAX, mpiComm);

      assert(nnzDistSharedMax == nnzDistShared);
      for( int i = 0; i < sizeKkt; i++ )
         assert(rowSizeShared[i] == rowSizeSharedMax[i]);
   }
#endif

   int localGatheredMyStart;
   int localGatheredMyEnd;

   std::vector<int> rowIndexGathered = PIPSallgathervInt(rowIndexMyLocal, mpiComm);
   std::vector<int> colIndexGathered = PIPSallgathervInt(colIndexMyLocal, mpiComm, localGatheredMyStart, localGatheredMyEnd);

#ifndef NDEBUG
   assert(int(rowIndexGathered.size()) == nnzDistLocal);
   assert(int(colIndexGathered.size()) == nnzDistLocal);
   assert(localGatheredMyEnd - localGatheredMyStart == nnzDistMyLocal);

   for( int i = 0; i < int(rowIndexMyLocal.size()); i++ )
   {
      assert(rowIndexMyLocal[i] == rowIndexGathered[i + localGatheredMyStart]);
      assert(colIndexMyLocal[i] == colIndexGathered[i + localGatheredMyStart]);
   }
#endif

   const int nnzDist = nnzDistLocal + nnzDistShared;

   assert(!kktDist || !kktDist->isLower);

   delete kktDist;
   kktDist = new SparseSymMatrix(sizeKkt, nnzDist, false);

   int* const krowDist = kktDist->krowM();
   int* const jColDist  = kktDist->jcolM();
   double* const MDist = kktDist->M();

   assert(krowDist[0] == 0);
   assert(sizeKkt > 0 && krowDist[1] == 0);

   memset(MDist, 0, nnzDist * sizeof(double));

   for( int r = 1; r < sizeKkt; r++ )
      krowDist[r + 1] = krowDist[r] + rowSizeLocal[r - 1] + rowSizeShared[r - 1];

   // fill in global and locally owned positions and values
   for( int r = 0; r < sizeKkt; r++ )
   {
      if( rowIsMyLocal[r] )
      {
         for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
         {
            assert(krowDist[r + 1] < nnzDist);

            const int col = jColKkt[c];

            if( col < 0 )
               continue;

            const double val = MKkt[c];

            MDist[krowDist[r + 1]] = val;
            jColDist[krowDist[r + 1]++] = col;
         }

         continue;
      }

      for( int c = krowKkt[r]; c < krowKkt[r + 1]; c++ )
      {
         const int col = jColKkt[c];

         if( col < 0 )
         {
            assert(-col - 1 >= 0 && -col - 1 < sizeKkt);
            assert(!(!rowIsLocal[r] && !rowIsLocal[-col - 1]));
            continue;
         }

         // is (r, col) a shared entry or locally owned?
         if( (!rowIsLocal[r] && !rowIsLocal[col]) || rowIsMyLocal[col] )
         {
            assert(krowDist[r + 1] < nnzDist);

            const double val = MKkt[c];

            MDist[krowDist[r + 1]] = val;
            jColDist[krowDist[r + 1]++] = col;
         }
      }
   }

   precondSC.resetSCdistEntries(kkts);

   // fill in gathered local pairs not inserted yet
   for( int i = 0; i < nnzDistLocal; i++ )
   {
      const int row = rowIndexGathered[i];
      const int col = colIndexGathered[i];

      assert(row >= 0 && row < sizeKkt);
      assert(col >= row && col < sizeKkt);

      // pair already added?
      if( i >= localGatheredMyStart && i < localGatheredMyEnd )
         continue;

      assert(krowDist[row + 1] < nnzDist);
      assert(MDist[krowDist[row + 1]] == 0.0);
      jColDist[krowDist[row + 1]++] = col;
   }

#ifndef NDEBUG
   assert(krowDist[0] == 0);
   assert(krowDist[sizeKkt] == nnzDist);

   for( int r = 0; r < sizeKkt; r++ )
   {
      assert(krowDist[r + 1] == krowDist[r] + rowSizeLocal[r] + rowSizeShared[r]);
      assert(krowDist[r + 1] >= krowDist[r]);
   }
#endif

   kktDist->getStorageRef().sortCols();

   assert(kktDist->getStorageRef().isValid());

   if( allreduce_kkt )
      reduceToAllProcs(nnzDist, MDist);
   else
      reduceToProc0(nnzDist, MDist);

   assert(kktDist->getStorageRef().isValid());
   assert(kktDist->getStorageRef().isSorted());
}

void sLinsysRoot::factorizeKKT()
{
   factorizeKKT(nullptr);
}

void sLinsysRoot::factorizeKKT(sData* prob)
{
  //stochNode->resMon.recFactTmLocal_start();  
#ifdef TIMING
  MPI_Barrier(mpiComm);
  extern double g_iterNumber;
  double st=MPI_Wtime();
#endif
  if( is_hierarchy_root )
     assert( !usePrecondDist );

  if( usePrecondDist )
  {
     assert( n_solvers == 1 );
     const int myRank = PIPS_MPIgetRank(mpiComm);

     assert(kktDist);
     assert(prob);

     if( allreduce_kkt || myRank == 0)
        precondSC.getSparsifiedSC_fortran(*prob, *kktDist);

     // todo do that properly
     precondSC.updateStats();

#if 0
      {
         ofstream myfile;
         int mype; MPI_Comm_rank(mpiComm, &mype);

         if( mype == 0 )
         {
            printf("\n\n ...WRITE OUT kktDist! \n\n");
            myfile.open("../ADist.txt");
            int* ia = kktDist->krowM(); int* ja = kktDist->jcolM(); double* a = kktDist->M();

            for( int i = 0; i < kktDist->size(); i++ )
               for( int k = ia[i]; k < ia[i + 1]; k++ )
                  myfile << i << '\t' << ja[k - 1] << '\t' << a[k - 1] << endl;

            myfile.close();
         }

         MPI_Barrier(mpiComm);
         printf("...exiting (root) \n");
         exit(1);
      }
#endif

     solver->matrixRebuild(*kktDist);
  }
  else
  {
     if( hasSparseKkt )
     {
       #pragma omp parallel num_threads(n_solvers)
       {
         omp_set_num_threads(n_threads_solvers);

         const int id = omp_get_thread_num();
         solvers_blocked[id]->matrixChanged();
       }
     }
     else
        solver->matrixChanged();
  }

  //stochNode->resMon.recFactTmLocal_stop(); 
#ifdef TIMING
  st = MPI_Wtime()-st;
  MPI_Barrier(mpiComm);
  int mype; MPI_Comm_rank(mpiComm, &mype);
  // note, this will include noop scalapack processors
  if( (mype/512)*512==mype )
    printf("  rank %d 1stSTAGE FACT %g SEC ITER %d\n", mype, st, (int)g_iterNumber);
#endif
}

//faster than DenseSymMatrix::atPutZeros
void sLinsysRoot::myAtPutZeros(DenseSymMatrix* mat, 
			       int row, int col, 
			       int rowExtent, int colExtent)
{
  assert( row >= 0 && row + rowExtent <= mat->size() );
  assert( col >= 0 && col + colExtent <= mat->size() );

  double ** M = mat->getStorageRef().M;

  for(int j=col; j<col+colExtent; j++) {
      M[row][j] = 0.0;
  }

  int nToCopy = colExtent*sizeof(double);

  for(int i=row+1; i<row+rowExtent; i++) {
    memcpy(M[i]+col, M[row]+col, nToCopy);
  }
}

void sLinsysRoot::myAtPutZeros(DenseSymMatrix* mat)
{
  int n = mat->size();
  myAtPutZeros(mat, 0, 0, n, n);
}

void sLinsysRoot::addTermToSchurCompl(sData* prob, size_t childindex, bool use_local_RAC )
{
   assert(childindex < prob->children.size());

   if( computeBlockwiseSC )
      children[childindex]->addTermToSchurComplBlocked(prob->children[childindex], hasSparseKkt, *kkt, use_local_RAC );
   else
   {
	   if( hasSparseKkt )
	   {
	      SparseSymMatrix& kkts = dynamic_cast<SparseSymMatrix&>(*kkt);

	      children[childindex]->addTermToSparseSchurCompl(prob->children[childindex], kkts);
	   }
	   else
	   {
		  DenseSymMatrix& kktd = dynamic_cast<DenseSymMatrix&>(*kkt);
		  children[childindex]->addTermToDenseSchurCompl(prob->children[childindex], kktd);
	   }
   }
}

void sLinsysRoot::submatrixAllReduce(DenseSymMatrix* A,
		             int startRow, int startCol, int nRows, int nCols,
				     MPI_Comm comm)
{
  double ** M = A->mStorage->M;
  int n = A->mStorage->n;

  assert(nRows > 0);
  assert(nCols > 0);
  assert(startRow >= 0);
  assert(startCol >= 0);

  int endRow = startRow + nRows;
  int endCol = startCol + nCols;

  assert(n >= endRow);
  assert(n >= endCol);

  int chunk_size = (CHUNK_SIZE / n) * n;
  chunk_size = std::min(chunk_size, n*nRows);

  double* chunk = new double[chunk_size];

  int rows_in_chunk = chunk_size/n;

  int iRow=startRow;

  // main loop
  do {

    if( iRow + rows_in_chunk > endRow )
      rows_in_chunk = endRow - iRow;

    assert(rows_in_chunk > 0);
#ifndef NDEBUG
    const int iErr=MPI_Allreduce(&M[iRow][0], chunk, rows_in_chunk*n, MPI_DOUBLE, MPI_SUM, comm);
    assert(iErr==MPI_SUCCESS);
#else
    MPI_Allreduce(&M[iRow][0], chunk, rows_in_chunk*n, MPI_DOUBLE, MPI_SUM, comm);
#endif

    int shift = 0;

    // copy into M
    for( int i = iRow; i < iRow + rows_in_chunk; i++ ) {
      for( int j = startCol; j < endCol; j++ )
	    M[i][j] = chunk[shift+j];

      // shift one row forward
      shift += n;
    }
    iRow += rows_in_chunk;

  } while( iRow < endRow );

  delete[] chunk;
}

// TODO: move all this to the respective matrix and storages.......
void sLinsysRoot::allreduceMatrix( DoubleMatrix& mat, bool is_sparse, bool is_sym, MPI_Comm comm )
{
   int n,m; mat.getSize(m,n);

   if( is_sparse )
   {
      if( is_sym )
      {
         SparseSymMatrix& matsp = dynamic_cast<SparseSymMatrix&>(mat);

         int* const krowKkt = matsp.krowM();
         double* const MKkt = matsp.M();
         const int nnzKkt = krowKkt[m];

         assert(!matsp.isLower);

         reduceToAllProcs(nnzKkt, MKkt);
      }
      else
      {
         SparseGenMatrix& matsp = dynamic_cast<SparseGenMatrix&>(mat);

         int* const krowKkt = matsp.krowM();
         double* const MKkt = matsp.M();
         const int nnzKkt = krowKkt[m];

         reduceToAllProcs(nnzKkt, MKkt);
      }
   }
   else
   {
      if( is_sym )
         submatrixAllReduceFull( &dynamic_cast<DenseSymMatrix&>(mat), 0, 0, m, n, comm);
      else
         submatrixAllReduceFull( &dynamic_cast<DenseGenMatrix&>(mat), 0, 0, m, n, comm);
   }
}

void sLinsysRoot::submatrixAllReduceFull(DenseSymMatrix* A, int startRow, int startCol, int nRows, int nCols, MPI_Comm comm)
{
   double** const M = A->mStorage->M;
   assert(A->mStorage->n == A->mStorage->m);
   assert(A->mStorage->n >= startRow + nRows);
   assert(A->mStorage->n >= startCol + nCols);

   submatrixAllReduceFull(M, startRow, startCol, nRows, nCols, comm);
}

void sLinsysRoot::submatrixAllReduceFull(DenseGenMatrix* A,
                   int startRow, int startCol, int nRows, int nCols,
                 MPI_Comm comm)
{
   double** const M = A->mStorage->M;
   assert(A->mStorage->m >= startRow + nRows);
   assert(A->mStorage->n >= startCol + nCols);

   submatrixAllReduceFull(M, startRow, startCol, nRows, nCols, comm);
}

void sLinsysRoot::submatrixAllReduceFull(double** A, int startRow, int startCol, int nRows, int nCols, MPI_Comm comm)
{
   assert(nRows > 0);
   assert(nCols > 0);
   assert(startRow >= 0);
   assert(startCol >= 0);

   const int endRow = startRow + nRows;
   const int buffersize = nRows * nCols;

   double* const bufferSend = new double[buffersize];
   double* const bufferRecv = new double[buffersize];

   // copy into send buffer
   int counter = 0;
   const size_t nColBytes = nCols * sizeof(double);

   for( int r = startRow; r < endRow; r++ )
   {
      memcpy(&bufferSend[counter], &A[r][startCol], nColBytes);
      counter += nCols;
   }

   assert(counter == buffersize);

#ifndef NDEBUG
   const int iErr = MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
   assert(iErr == MPI_SUCCESS);
#else
   MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
#endif

   // copy back
   counter = 0;
   for( int r = startRow; r < endRow; r++ )
   {
      memcpy(&A[r][startCol], &bufferRecv[counter], nColBytes);
      counter += nCols;
   }

   assert(counter == buffersize);

   delete[] bufferRecv;
   delete[] bufferSend;
}


void sLinsysRoot::submatrixAllReduceDiagLower(DenseSymMatrix* A,
                   int substart, int subsize,
                 MPI_Comm comm)
{
   double** const M = A->mStorage->M;

   assert(subsize >= 0);
   assert(substart >= 0);

   if( subsize == 0)
      return;

   const int subend = substart + subsize;
   assert(A->mStorage->n >= subend);

   // number of elements in lower matrix triangle (including diagonal)
   const int buffersize = (subsize * subsize + subsize) / 2;
   assert(buffersize > 0);

   double* const bufferSend = new double[buffersize];
   double* const bufferRecv = new double[buffersize];

   int counter = 0;

   for( int i = substart; i < subend; i++ )
      for( int j = substart; j <= i; j++ )
      {
         assert(counter < buffersize);
         bufferSend[counter++] = M[i][j];
      }

   assert(counter == buffersize);

#ifndef NDEBUG
   const int iErr = MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
   assert(iErr == MPI_SUCCESS);
#else
   MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
#endif

   counter = 0;
   for( int i = substart; i < subend; i++ )
      for( int j = substart; j <= i; j++ )
      {
         assert(counter < buffersize);
         M[i][j] = bufferRecv[counter++];
      }

   delete[] bufferSend;
   delete[] bufferRecv;
}


#ifdef STOCH_TESTING
void sLinsysRoot::dumpMatrix(int scen, int proc, const char* nameToken, DenseSymMatrix& M) 
{
  int n = M.size();
  char szNumber[30];
  string strBuffer="";

  //assert(false);

  int iter = g_iterNumber;

  if(iter!=1 && iter!=5 && iter!=15 && iter!=25 && iter!=35 && iter!=45) return;


  char szFilename[256];
  if(scen==-1)
    sprintf(szFilename, "%s_%d__%d.mat", nameToken, n, iter);
  else 
    sprintf(szFilename, "%s_%03d_%d__%d.mat", nameToken, scen+1, n, iter);
  FILE* file = fopen(szFilename, "w");
  assert(file);
  

  for(int j=0; j<n; j++) {
    for(int i=0; i<n; i++) {
      sprintf(szNumber, "%22.16f ", M[i][j]);
      strBuffer += szNumber;
    }
    strBuffer += "\n";
    
    if(strBuffer.length()>1250000) {
      fwrite(strBuffer.c_str(), 1, strBuffer.length(), file);
      strBuffer = "";
    }
  }

  if(strBuffer.length()>0) {
    fwrite(strBuffer.c_str(), 1, strBuffer.length(), file);
  }
  
  fclose(file);
}

void sLinsysRoot::dumpRhs(int proc, const char* nameToken,  SimpleVector& rhs) 
{
  int n = rhs.length();
  char szNumber[30];
  string strBuffer="";


  int iter = g_iterNumber;
  if(iter!=0 && iter!=2 && iter!=20 && iter!=25 && iter!=55) return;

  char ipmPhase[4];
  if(g_iterNumber-iter>0) strcpy(ipmPhase, "co");
  else strcpy(ipmPhase, "pr");

  char szFilename[256];
  sprintf(szFilename, "%s_%s_%d__%d.mat", nameToken,ipmPhase, n, iter);


  for(int i=0; i<n; i++) {
    sprintf(szNumber, "%22.16f ", rhs[i]);
    strBuffer += szNumber;
  }

  FILE* file = fopen(szFilename, "w");
  assert(file);

  fwrite(strBuffer.c_str(), 1, strBuffer.length(), file);
  fclose(file);
}

#endif
