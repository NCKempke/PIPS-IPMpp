/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#include "StochOptions.h"
#include "BorderedSymMatrix.h"
#include "sLinsys.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "SparseLinearAlgebraPackage.h"
#include "math.h"

#include "pipsport.h"
#include "omp.h"

sLinsys::sLinsys(sFactory* factory_, sData* prob, bool is_hierarchy_root)
  : QpGenLinsys(factory_, prob), data{prob},
    computeBlockwiseSC( pips_options::getBoolParameter("SC_COMPUTE_BLOCKWISE") ),
    blocksizemax( pips_options::getIntParameter("SC_BLOCKWISE_BLOCKSIZE_MAX") ),
    is_hierarchy_root(is_hierarchy_root),
    stochNode{ factory_->tree }
{
  if( pips_options::getBoolParameter( "HIERARCHICAL" ) )
    assert( is_hierarchy_root );

  prob->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);

  //get the communicator from one of the vectors
  StochVector& dds = dynamic_cast<StochVector&>(*dd);
  this->mpiComm = dds.mpiComm;
  this->iAmDistrib = dds.iAmDistrib;
}

sLinsys::sLinsys(sFactory* factory_,
		 sData* prob,				    
		 OoqpVector* dd_, 
		 OoqpVector* dq_,
		 OoqpVector* nomegaInv_,
		 OoqpVector* rhs_,
		 bool create_iter_ref_vecs
		 )
  : QpGenLinsys( factory_, prob, dd_, dq_, nomegaInv_, rhs_, create_iter_ref_vecs ),
    data{prob},
    computeBlockwiseSC( pips_options::getBoolParameter("SC_COMPUTE_BLOCKWISE") ),
    blocksizemax( pips_options::getIntParameter("SC_BLOCKWISE_BLOCKSIZE_MAX") ),
    stochNode{factory_->tree}
{
  prob->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);

  if( dd )
  {
     StochVector& dds = dynamic_cast<StochVector&>(*dd);
     mpiComm = dds.mpiComm;
     iAmDistrib = dds.iAmDistrib;
  }
  else
  {
     mpiComm = MPI_COMM_NULL;
     iAmDistrib = false;
  }

  useRefs = 1;
}


sLinsys::~sLinsys()
{
  if( colSparsity ) delete[] colSparsity;
  if( colId ) delete[] colId;
  if( colsBlockDense ) delete[] colsBlockDense;
}

void sLinsys::joinRHS( OoqpVector& rhs_in, const OoqpVector& rhs1_in,
				const OoqpVector& rhs2_in, const OoqpVector& rhs3_in ) const
{
  StochVector& rhs  = dynamic_cast<StochVector&>(rhs_in);
  const StochVector& rhs1 = dynamic_cast<const StochVector&>(rhs1_in);
  const StochVector& rhs2 = dynamic_cast<const StochVector&>(rhs2_in);
  const StochVector& rhs3 = dynamic_cast<const StochVector&>(rhs3_in);

  rhs.jointCopyFrom(rhs1, rhs2, rhs3);
}

void sLinsys::separateVars( OoqpVector& x_in, OoqpVector& y_in,
				     OoqpVector& z_in, const OoqpVector& vars_in ) const
{
  StochVector& x = dynamic_cast<StochVector&>(x_in);
  StochVector& y = dynamic_cast<StochVector&>(y_in);
  StochVector& z = dynamic_cast<StochVector&>(z_in);
  const StochVector& vars = dynamic_cast<const StochVector&>(vars_in);

  vars.jointCopyTo(x, y, z);
}

void sLinsys::factor(Data *prob_, Variables *vars)
{
#ifdef TIMING
  double tTot = MPI_Wtime();
#endif
  // the call to the the parent's method takes care of all necessary updates
  // to the KKT system (updating diagonals mainly). This is done recursively,
  // we don't have to worry about it anymore. 
  QpGenLinsys::factor(prob_, vars);

  // now DO THE LINEAR ALGEBRA!
  
  sData* prob = dynamic_cast<sData*>(prob_);
  // in order to avoid a call to QpGenLinsys::factor, call factor2 method.
  // factor2(prob, vars);
  assembleKKT(prob, vars);
  allreduceAndFactorKKT(prob, vars);

#ifdef TIMING
  tTot = MPI_Wtime() - tTot;
  MPI_Barrier(MPI_COMM_WORLD);
  const int rank == PIPS_MPIgetRANK(MPI_COMM_WORLD);
  //if( 128 * ( myRank / 128 ) == 0 )
  if( 0 == myRank )
      std::cout << "Outer fact. total time " << tTot << std::endl;
#endif
}
 

/**
 * Computes U = Li\Gi^T.
 *        [ 0 0 0 ]
 * Gi^T = [ A 0 0 ]
 *        [ C 0 0]
 *
 * We have special structure here:
 *             [ 0 ]
 *   U   = Li\ [ A ] ,   U is (nx+my+mz)-by-(np)
 *             [ C ]
 *
 *   V = Di\U
 */
void sLinsys::computeU_V(sData *prob, 
			 DenseGenMatrix* U, DenseGenMatrix* V)
{
  U->scalarMult(0.0);
  V->scalarMult(0.0);
  assert(false); //need code to deal with cross Hessian term
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  int N, nxP;
  A.getSize(N, nxP); assert(N==locmy);

  N = locnx+locmy+locmz;
  SimpleVector uCol(N);

  for(int it=0; it<nxP; it++) {
    
    double* p = &uCol[0];
    for(int it1=0; it1<locnx; it1++) p[it1]=0.0;

    A.fromGetDense(0, it, &uCol[locnx],  1, locmy, 1);    
    C.fromGetDense(0, it, &uCol[locnx+locmy], 1, locmz, 1);

    solver->Lsolve(uCol);
    U->atPutDense(0, it, &uCol[0], 1, N, 1);

    solver->Dsolve(uCol);
    V->atPutDense(0, it, &uCol[0], 1, N, 1);

  }
}

void sLinsys::allocU(DenseGenMatrix ** U, int n0)
{
  int lines,cols;
  if(*U==nullptr) {
    *U = new DenseGenMatrix(locnx + locmy + locmz, n0);
  } else {
    (*U)->getSize(lines,cols);
    
    if(lines!=locnx+locmy+locmz || cols != n0) {

      delete (*U);
      *U = new DenseGenMatrix(locnx+locmy+locmz, n0);
    }
  }
}

void sLinsys::allocV(DenseGenMatrix ** V, int n0)
{
  int lines,cols;
  if(*V==nullptr)
    *V = new DenseGenMatrix(locnx+locmy+locmz, n0);
  else {
    (*V)->getSize(lines,cols);
    
    if(lines!=locnx+locmy+locmz || cols != n0) {
      
      delete (*V);
      *V = new DenseGenMatrix(locnx + locmy + locmz, n0);
    }
  }
}

/**
 *       [ R^i^T Ai^T Ci^T ]          [    ]
 * z0 -= [ 0      0   0    ] * Li\Di\ [ zi ]
 *       [ 0      0   0    ]          [    ]
 *
 * 
 */
void sLinsys::addLnizi(sData *prob, OoqpVector& z0_, OoqpVector& zi_)
{
  SimpleVector& z0 = dynamic_cast<SimpleVector&>(z0_);
  SimpleVector& zi = dynamic_cast<SimpleVector&>(zi_);

  solver->Dsolve(zi);
  solver->Ltsolve(zi);

  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();

  //get n0= nx(parent)= #cols of A or C
  int dummy, n0;
  A.getSize(dummy, n0);

  // zi2 and zi3 are just references to fragments of zi
  SimpleVector zi1 (&zi[0],           locnx);
  SimpleVector zi2 (&zi[locnx],       locmy);
  SimpleVector zi3 (&zi[locnx+locmy], locmz);
  // same for z01 (only the first n0 entries in the output z0 are computed)
  SimpleVector z01 (&z0[0], n0);

  R.transMult(1.0, z01, -1.0, zi1);
  A.transMult(1.0, z01, -1.0, zi2);
  C.transMult(1.0, z01, -1.0, zi3);
}

void sLinsys::finalizeDenseBorderModBlocked( std::vector<BorderMod>& border_mod, DenseGenMatrix& result )
{
   /* compute BiT_buffer += X_j^T Bmodj for all j */
   for( auto& border_mod_block : border_mod )
      finalizeDenseBorderBlocked( border_mod_block.border, border_mod_block.multiplier, result );
}

void sLinsys::multRightDenseBorderModBlocked( std::vector<BorderMod>& border_mod, DenseGenMatrix& result )
{
   /* compute BiT_buffer += X_j^T Bmodj for all j */
   for( auto& border_mod_block : border_mod )
   {
      std::unique_ptr<BorderBiBlock> BiT_mod{};

      BorderLinsys& border = border_mod_block.border;

      if( border.use_local_RAC )
         BiT_mod.reset( new BorderBiBlock( data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(), data->getLocalC().getTranspose(),
               dynamic_cast<SparseGenMatrix&>(*border.F.mat), dynamic_cast<SparseGenMatrix&>(*border.G.mat) ) );
      else if( border.has_RAC )
         BiT_mod.reset( new BorderBiBlock ( dynamic_cast<SparseGenMatrix&>(*border.R.mat).getTranspose(),
               dynamic_cast<SparseGenMatrix&>(*border.A.mat).getTranspose(), dynamic_cast<SparseGenMatrix&>(*border.C.mat).getTranspose(),
               dynamic_cast<SparseGenMatrix&>(*border.F.mat), dynamic_cast<SparseGenMatrix&>(*border.G.mat) ) );
      else
         BiT_mod.reset( new BorderBiBlock ( dynamic_cast<SparseGenMatrix&>(*border.F.mat), dynamic_cast<SparseGenMatrix&>(*border.G.mat), false ) );

      multRightDenseBorderBlocked( *BiT_mod, border_mod_block.multiplier, result );
   }
}

/* compute
 *              locnx locmy     locmyl locmzl
 * nx_border  [   0    A0T  C0T F0VT G0VT ]
 * myl_border [  F0C    0    0   0    0   ]
 * mzl_border [  G0C    0    0   0    0   ]
 *
 *               [  0 F0C^T  G0C^T ]^T
 *               [ A0   0     0    ]
 *               [ C0   0     0    ]
 * buffer -= X * [ F0V  0     0    ]
 *               [ G0V  0     0    ]
 */
void sLinsys::finalizeDenseBorderBlocked( BorderLinsys& B, const DenseGenMatrix& X, DenseGenMatrix& result )
{
   const bool has_RAC = B.has_RAC;

   int mX0, nX0; X.getSize( mX0, nX0 );
   int mRes, nRes; result.getSize( mRes, nRes );
   assert( mX0 == mRes );

   SparseGenMatrix* F0cons_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.F.mat) : nullptr;
   SparseGenMatrix* G0cons_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.G.mat) : nullptr;

   SparseGenMatrix* A0_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.A.mat) : nullptr;
   SparseGenMatrix* C0_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.C.mat) : nullptr;
   SparseGenMatrix* F0vec_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.A.mat_link) : dynamic_cast<SparseGenMatrix*>(B.F.mat);
   SparseGenMatrix* G0vec_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.C.mat_link) : dynamic_cast<SparseGenMatrix*>(B.G.mat);

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

   int mF0V{0}; int nF0V{0};
   F0vec_border->getSize(mF0V, nF0V);

#ifndef NDEBUG
   int mG0C{0}; int nG0C{0};
   if( G0cons_border )
      G0cons_border->getSize( mG0C, nG0C );

   int mG0V{0}; int nG0V{0};
   G0vec_border->getSize(mG0V, nG0V);

   assert( nA0 == nC0 );
   assert( nF0V == nG0V );

   if( has_RAC )
      assert( nA0 == nF0V );

   assert( nF0C + mA0 + mC0 + mF0V + mG0V == nRes );
   assert( nF0C == nG0C );

   if( has_RAC )
      assert( nX0 == nF0V + mF0C + mG0C );
   else
      assert( nX0 >= nF0V );
#endif



   if( has_RAC )
   {
      /*            [  0  ]
       * res -= X * [ F0C ]
       *            [ G0C ]
       */
      X.multMatAt( nF0V, 1.0, 0, result, -1.0, *F0cons_border );
      X.multMatAt( nF0V + mF0C, 1.0, 0, result, -1.0, *G0cons_border );

      /*            [ A0T ]
       * res -= X * [  0  ]
       *            [  0  ]
       */
      X.multMatAt( 0, 1.0, nF0C, result, -1.0, A0_border->getTranspose() );

      /*            [ C0T ]
       * res -= X * [  0  ]
       *            [  0  ]
       */
      X.multMatAt( 0, 1.0, nF0C + mA0, result, -1.0, C0_border->getTranspose() );
   }
   /*            [ F0VT ]
    * res -= X * [  0   ]
    *            [  0   ]
    */
   X.multMatAt( 0, 1.0, nF0C + mA0 + mC0, result, -1.0, F0vec_border->getTranspose() );

   /*            [ G0VT ]
    * res -= X * [  0   ]
    *            [  0   ]
    */
   X.multMatAt( 0, 1.0, nF0C + mA0 + mC0 + mF0V, result, -1.0, G0vec_border->getTranspose() );
}


/* calculate res -= X * BT */
void sLinsys::multRightDenseBorderBlocked( BorderBiBlock& BT, const DenseGenMatrix& X, DenseGenMatrix& result )
{
   /*
    *        [  RiT   AiT   CiT ]
    *        [   0     0     0  ]
    * Bi^T = [   0     0     0  ]
    *        [   Fi    0     0  ]
    *        [   Gi    0     0  ]
    */

   const bool with_RAC = BT.has_RAC;
   int mX, nX; X.getSize(mX, nX);
   int mR, nR; BT.R.getSize( mR, nR);
   int mA, nA; BT.A.getSize( mA, nA );
   int mF, nF; BT.F.getSize( mF, nF );
   int mG, nG; BT.G.getSize( mG, nG );

#ifndef NDEBUG
   int mRes, nRes; result.getSize(mRes, nRes);
   assert( mX == mRes );

   assert( nF == nG );
   if( with_RAC )
   {
      int mC, nC; BT.C.getSize( mC, nC );
      assert( mR == mA );
      assert( mR == mC );
      assert( nR == nF );

      assert( mR + mF + mG <= nX );
   }
   else
      assert( mF + mG == nX );
#endif
   // X from the right with each column of Bi^T todo add OMP to submethods

   /*            [ RiT ]
    *            [  0  ]
    * res -= X * [  0  ]
    *            [  Fi ]
    *            [  Gi ]
    */
   if( with_RAC )
      X.multMatAt( 0, 1.0, 0, result, -1.0, BT.R );

   if( mF > 0 )
      X.multMatAt( nX - mF - mG, 1.0, 0, result, -1.0, BT.F );

   if( mG > 0 )
      X.multMatAt( nX - mG, 1.0, 0, result, -1.0, BT.G );

   if( with_RAC )
   {
      /*            [ AiT ]
       *            [  0  ]
       * res -= X * [  0  ]
       *            [  0  ]
       *            [  0  ]
       */
      X.multMatAt( 0, 1.0, nR, result, -1.0, BT.A );

      /*            [ CiT ]
       *            [  0  ]
       * res -= X * [  0  ]
       *            [  0  ]
       *            [  0  ]
       */
      X.multMatAt( 0, 1.0, nR + nA, result, -1.0, BT.C );
   }
}

void sLinsys::addMatAt( DenseGenMatrix& res, const SparseGenMatrix& mat, int row_0, int col_0 ) const
{
   int mres, nres; res.getSize( mres, nres );
   int mmat, nmat; mat.getSize( mmat, nmat );

   assert( 0 <= row_0 && row_0 + mmat <= mres );
   assert( 0 <= col_0 && col_0 + nmat <= nres );

   for( int row = 0; row < mmat; ++row )
   {
      const int row_start = mat.krowM()[row];
      const int row_end = mat.krowM()[row + 1];

      for( int j = row_start; j < row_end; ++j )
      {
         const int col = mat.jcolM()[j];
         const double val = mat.M()[j];

         assert( col_0 + col < nres );
         res[row_0 + row][col_0 + col] += val;
      }
   }

}

void sLinsys::addBiTBorder( DenseGenMatrix& res, const BorderBiBlock& BiT ) const
{
   /* add (Bri)^T to res
    *
    *                [ RiT AiT CiT ]
    *                [  0   0   0  ]
    * Bri^T        = [  Fi  0   0  ]
    *                [  Gi  0   0  ]
    */

   int mRt, nRt; BiT.R.getSize(mRt, nRt);
   int mAt, nAt; BiT.A.getSize(mAt, nAt);
   int mCt, nCt; BiT.C.getSize(mCt, nCt);
   int mF, nF; BiT.F.getSize(mF, nF);
   int mG, nG; BiT.G.getSize(mG, nG);

   int mres, nres; res.getSize(mres, nres);

#ifndef NDEBUG
   assert( nF == nG );
   if( BiT.has_RAC )
   {
      assert( mRt + mF + mG <= mres );
      assert( nF == nRt );
      assert( nRt + nAt + nCt == nres );
   }
   else
   {
      assert( mres == mF + mG );
      assert( nres == nF );
   }
#endif

   if( BiT.has_RAC )
   {
      addMatAt( res, BiT.R, 0, 0 );
      addMatAt( res, BiT.A, 0, nRt );
      addMatAt( res, BiT.C, 0, nRt + nAt );
   }

   addMatAt( res, BiT.F, mres - mF - mG, 0 );
   addMatAt( res, BiT.G, mres - mG, 0 );
}

void sLinsys::solveCompressed( OoqpVector& rhs_)
{
   StochVector& rhs = dynamic_cast<StochVector&>(rhs_);
#ifdef TIMING
  //double tTot=MPI_Wtime();
#endif
   Lsolve(data, rhs);
   Dsolve(data, rhs);
   Ltsolve(data, rhs);
#ifdef TIMING
  //cout << "SolveCompressed took: " << (MPI_Wtime()-tTot) << endl;
#endif
}


/*
 *  y = alpha*Lni^T x + beta*y
 *
 *                       ( [ R 0 0 ]     )
 *  y = beta*y + Di\Li\ (  [ A 0 0 ] * x )
 *                      (  [ C 0 0 ]    )
 */
void sLinsys::LniTransMult(sData *prob, 
			   SimpleVector& y, 
			   double alpha, SimpleVector& x)
{
  SparseGenMatrix& A = prob->getLocalA();
  int N{0}, nx0{0};

  //get nx(parent) from the number of cols of A (or C). Use N as dummy
  if( data->hasRAC() )
     A.getSize(N,nx0);
  // a mild assert
  assert(nx0 <= x.length());

  N = locnx+locmy+locmz;
  assert(y.length() == N);
  
  //!memopt
  SimpleVector LniTx(N);

  // shortcuts
  SimpleVector x1(&x[0], nx0);
  SimpleVector LniTx1(&LniTx[0], locnx);

  LniTx1.setToZero();
  if( data->hasRAC() )
  {
     SimpleVector LniTx2(&LniTx[locnx], locmy);
     SimpleVector LniTx3(&LniTx[locnx+locmy], locmz);

     SparseGenMatrix& C = prob->getLocalC();
     SparseGenMatrix& R = prob->getLocalCrossHessian();
     R.mult(0.0, LniTx1, 1.0, x1);
     A.mult(0.0, LniTx2, 1.0, x1);
     C.mult(0.0, LniTx3, 1.0, x1);

  }

  if( locmyl > 0 )
  {
	 int nxMyMzP = x.length() - locmyl - locmzl;

	 SparseGenMatrix& F = prob->getLocalF();
    SimpleVector xlink(&x[nxMyMzP], locmyl);

    F.transMult(1.0, LniTx1, 1.0, xlink);
  }

  if( locmzl > 0 )
  {
    int nxMyMzMylP = x.length() - locmzl;

    SparseGenMatrix& G = prob->getLocalG();
    SimpleVector xlink(&x[nxMyMzMylP], locmzl);

    G.transMult(1.0, LniTx1, 1.0, xlink);
  }

//  solver->Lsolve(LniTx); -> empty
  solver->Dsolve(LniTx);
  y.axpy(alpha,LniTx);
}


/*
 * Computes res += [R^T A^T C^T ] * inv(KKT) * [R 0 F^T G^T ] x
 *                 [0         ]              [A             ]
 *                 [F         ]              [C             ]
 *                 [G         ]
 */

void sLinsys::addTermToSchurResidual(sData* prob, 
				     SimpleVector& res, 
				     SimpleVector& x)
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& F = prob->getLocalF();
  SparseGenMatrix& G = prob->getLocalG();
  SparseGenMatrix& R = prob->getLocalCrossHessian();

  int nxP, aux;
  A.getSize(aux,nxP); assert(aux==locmy);
  C.getSize(aux,nxP); assert(aux==locmz);
  F.getSize(aux,nxP); assert(aux==locmyl);
  G.getSize(aux,nxP); assert(aux==locmzl);
  R.getSize(aux,nxP); assert(aux==locnx);

  // res contains mz buffer part
  assert(res.length() >= x.length());
  assert(x.length() >= nxP);

  int N=locnx+locmy+locmz;
  SimpleVector y(N);

  R.mult( 0.0,&y[0],1,           1.0,&x[0],1);
  A.mult( 0.0,&y[locnx],1,       1.0,&x[0],1);
  C.mult( 0.0,&y[locnx+locmy],1, 1.0,&x[0],1);

  if( locmyl > 0 )
  {
     assert(res.length() == x.length());
     F.transMult( 1.0,&y[0],1,       1.0,&x[x.length() - locmyl - locmzl],1);
  }

  if( locmzl > 0 )
  {
     assert(res.length() == x.length());
     G.transMult( 1.0,&y[0],1,       1.0,&x[x.length() - locmzl],1);
  }

  //cout << "4 - y norm:" << y.twonorm() << endl;
  //printf("%g  %g  %g  %g\n", y[locnx+locmy+0], y[locnx+locmy+1], y[locnx+locmy+2], y[locnx+locmy+3]);
  solver->solve(y);

  R.transMult(1.0,&res[0],1, 1.0,&y[0],1);
  A.transMult(1.0,&res[0],1, 1.0,&y[locnx],1);
  C.transMult(1.0,&res[0],1, 1.0,&y[locnx+locmy],1);

  if( locmyl > 0 )
     F.mult(1.0,&res[res.length() - locmyl - locmzl],1, 1.0,&y[0],1);

  if( locmzl > 0 )
     G.mult(1.0,&res[res.length() - locmzl],1, 1.0,&y[0],1);
}

void sLinsys::addTermToDenseSchurCompl(sData *prob,
				       DenseSymMatrix& SC)
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& F = prob->getLocalF();
  SparseGenMatrix& G = prob->getLocalG();
  SparseGenMatrix& R = prob->getLocalCrossHessian();

  int N, nxP;

  R.getSize(N, nxP);
  const bool withR = (nxP != -1);

  A.getSize(N, nxP);
  const bool withA = (nxP != -1);

  assert(N==locmy);
  assert(locmyl >= 0);
  assert(locmzl >= 0);

  const int NP = SC.size();
  assert(NP>=nxP);

  const int nxMyP = NP - locmyl - locmzl;
  const int nxMyMzP = NP - locmzl;

  if(nxP==-1)
    C.getSize(N,nxP);

  int N2, nxP2;
  C.getSize(N2,nxP2);
  const bool withC = (nxP2 != -1);

  if(nxP==-1)
    nxP = NP;

  N = locnx+locmy+locmz;

  SimpleVector col(N);
  SimpleVectorBase<int> nnzPerColRAC(nxP);

  if( withR )
     R.addNnzPerCol(nnzPerColRAC);

  if( withA )
     A.addNnzPerCol(nnzPerColRAC);

  if( withC )
     C.addNnzPerCol(nnzPerColRAC);

  const int withMyl = (locmyl > 0);
  const int withMzl = (locmzl > 0);

  for(int it=0; it<nxP; it++) {
    if( nnzPerColRAC[it] == 0 )
      continue;

    double* const pcol = &col[0];

    for(int it1=0; it1<locnx; it1++) pcol[it1]=0.0;

    R.fromGetDense(0, it, &col[0],           1, locnx, 1);
    A.fromGetDense(0, it, &col[locnx],       1, locmy, 1);
    C.fromGetDense(0, it, &col[locnx+locmy], 1, locmz, 1);

    solver->solve(col);

    //here we have colGi = inv(H_i)* it-th col of Gi^t
    //now do colSC = Gi * inv(H_i)* it-th col of Gi^t

    // SC+=R*x
    R.transMult( 1.0, &SC[it][0],     1,  -1.0, &col[0],      1);

    // SC+=At*y
    A.transMult( 1.0, &SC[it][0],     1,  -1.0, &col[locnx],  1);

    // SC+=Ct*z
    C.transMult( 1.0, &SC[it][0],     1,  -1.0, &col[locnx+locmy], 1);

    // do we have linking equality constraints? If so, set SC+=F*x
    if( withMyl )
       F.mult( 1.0, &SC[it][nxMyP],     1, -1.0, &col[0],      1);

    // do we have linking inequality constraints? If so, set SC+=G*x
    if( withMzl )
       G.mult( 1.0, &SC[it][nxMyMzP],     1,  -1.0, &col[0],      1);
  }

  // do we have linking equality constraints?
  if( withMyl )
  {
    SimpleVectorBase<int> nnzPerColFt(locmyl);
    F.addNnzPerRow(nnzPerColFt);

    // do column-wise multiplication for columns containing Ft (F transposed)
    for(int it=0; it<locmyl; it++) {

      if( nnzPerColFt[it] == 0 )
         continue;

      double* pcol = &col[0];

      // get it'th column from Ft (i.e., it'th row from F)
      F.fromGetDense(it, 0, &col[0],           1, 1, locnx);

      for(int it1=locnx; it1 < locnx+locmy+locmz; it1++) pcol[it1]=0.0;

      solver->solve(col);

      R.transMult( 1.0, &SC[it + nxMyP][0],   1,  -1.0, &col[0],      1);
      A.transMult( 1.0, &SC[it + nxMyP][0],   1,  -1.0, &col[locnx],  1);
      C.transMult( 1.0, &SC[it + nxMyP][0],   1,  -1.0, &col[locnx+locmy], 1);

      // here we have colGi = inv(H_i)* (it + locnx + locmy)-th col of Gi^t
      // now do colSC = Gi * inv(H_i)* (it + locnx + locmy)-th col of Gi^t

      F.mult( 1.0, &SC[it + nxMyP][nxMyP],   1, -1.0, &col[0],  1);

      if( withMzl )
         G.mult( 1.0, &SC[it + nxMyP][nxMyMzP],   1, -1.0, &col[0],  1);
    }
  }

  // do we have linking inequality constraints?
  if( withMzl )
  {
    SimpleVectorBase<int> nnzPerColGt(locmzl);
    G.addNnzPerRow(nnzPerColGt);

    // do column-wise multiplication for columns containing Gt (G transposed)
    for(int it=0; it<locmzl; it++) {
      if( nnzPerColGt[it] == 0 )
         continue;

      double* pcol = &col[0];

      // get it'th column from Gt (i.e., it'th row from G)
      G.fromGetDense(it, 0, &col[0],           1, 1, locnx);

      for(int it1=locnx; it1 < locnx+locmy+locmz; it1++) pcol[it1]=0.0;

      solver->solve(col);

      R.transMult( 1.0, &SC[it + nxMyMzP][0],   1,  -1.0, &col[0],      1);
      A.transMult( 1.0, &SC[it + nxMyMzP][0],   1,  -1.0, &col[locnx],  1);
      C.transMult( 1.0, &SC[it + nxMyMzP][0],   1,  -1.0, &col[locnx+locmy], 1);

      // here we have colGi = inv(H_i)* (it + locnx + locmy + locmyl)-th col of Gi^t
      // now do colSC = Gi * inv(H_i)* (it + locnx + locmy + locmyl)-th col of Gi^t

      if( withMyl )
        F.mult( 1.0, &SC[it + nxMyMzP][nxMyP],   1,  -1.0, &col[0],  1);

      G.mult( 1.0, &SC[it + nxMyMzP][nxMyMzP],   1,  -1.0, &col[0],  1);
    }
  }
}

//#define TIME_SCHUR

/* Solving with
 *  [ Ri 0 FiT GiT]
 *  [ Ai 0  0   0 ]
 *  [ Ci 0  0   0 ]
 *
 * where the size of the zero part depends on the size of the factor K^-1:
 * l = n_K^-1 - nR - nF - nG
 *
 * computes (B_left_transp^T K^{-1} B_right)^T and adds it to the SC
 * matrices in border_left_transp and border_right are already transposed/not depending on what is wanted
 */
//void sLinsys::addBiTLeftKiBiRightToResBlocked( bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
//      /* const */ BorderBiBlock& border_right, DoubleMatrix& result)
//{
//   int m_res, n_res; result.getSize(m_res, n_res);
//   assert( m_res >= 0 && n_res >= 0 );
//   if( sym_res )
//      assert( m_res == n_res );
//
//   int mR_right, nR_right; border_right.R.getSize(mR_right, nR_right);
//
//   int mA_right, nA_right; border_right.A.getSize(mA_right, nA_right);
//   assert( nR_right == nA_right);
//
//   int mC_right, nC_right; border_right.C.getSize(mC_right, nC_right);
//   assert( nR_right == nC_right);
//
//   int mF_right, nF_right; border_right.F.getSize(mF_right, nF_right);
//   assert( mR_right == mF_right );
//
//   int mG_right, nG_right; border_right.G.getSize(mG_right, nG_right);
//   assert( mR_right == mG_right );
//
//   /* we add to the transposed of res */
//   assert( nR_right + nF_right + nG_right <= m_res);
//
//   SimpleVectorBase<int> nnzPerColRAC(nR_right);
//
//   border_right.R.addNnzPerCol(nnzPerColRAC);
//   border_right.A.addNnzPerCol(nnzPerColRAC);
//   border_right.C.addNnzPerCol(nnzPerColRAC);
//
//   const int withF = ( nF_right > 0 );
//   const int withG = ( nG_right > 0 );
//   const int length_col = mR_right + mA_right + mC_right;
//
//   if( colsBlockDense == nullptr )
//      colsBlockDense = new double[blocksizemax * length_col];
////   else
////   {
////      delete[] colsBlockDense;
////      colsBlockDense = new double[blocksizemax * length_col];
////   }
//
//   if( colId == nullptr )
//      colId = new int[blocksizemax];
//
//   // indicating whether a right hand side is zero - deactivated since problems in Pardiso
//#if 0
//   if( colSparsity == nullptr )
//      colSparsity = new int[length_col * blocksizemax];
//#else
//   colSparsity = nullptr;
//#endif
//
//#ifdef TIME_SCHUR
//   const double t_start = omp_get_wtime();
//#endif
//
//   int colpos = 0;
//
//   //                       (R)
//   //     SC +=  B^T  K^-1  (A)
//   //                       (C)
//   while( colpos < nR_right )
//   {
//      int blocksize = 0;
//
//      for( ; colpos < nR_right && blocksize < blocksizemax; colpos++ )
//         if( nnzPerColRAC[colpos] != 0 )
//            colId[blocksize++] = colpos;
//
//      if( blocksize == 0 )
//         break;
//
//      memset(colsBlockDense, 0, blocksize * length_col * sizeof(double));
//
//      if( colSparsity )
//         memset(colSparsity, 0, length_col * sizeof(int));
//
//      border_right.R.fromGetColsBlock(colId, blocksize, length_col, 0, colsBlockDense, colSparsity);
//      border_right.A.fromGetColsBlock(colId, blocksize, length_col, mR_right, colsBlockDense, colSparsity);
//      border_right.C.fromGetColsBlock(colId, blocksize, length_col, (mR_right + mA_right), colsBlockDense, colSparsity);
//
//      solvers_blocked[0]->solve(blocksize, colsBlockDense, colSparsity);
//
//      addLeftBorderTimesDenseColsToResTransp( border_left_transp, colsBlockDense, colId, length_col, blocksize, sparse_res, sym_res, result);
//   }
//
//#ifdef TIME_SCHUR
//   const double t_end = omp_get_wtime();
//   std::cout << "t_end - t_start:" << (t_end - t_start) << std::endl;
//   assert(0);
//#endif
//
//   // do we have linking equality constraints?
//   if( withF )
//   {
//      //                       (F^T)
//      //     SC +=  B^T  K^-1  (0  )
//      //                       (0  )
//
//      SimpleVectorBase<int> nnzPerColFt(nF_right);
//      border_right.F.addNnzPerCol(nnzPerColFt);
//
//      colpos = 0;
//      // do block-wise multiplication for columns of F^T part
//      while( colpos < nF_right )
//      {
//         int blocksize = 0;
//
//         for( ; colpos < nF_right && blocksize < blocksizemax; colpos++ )
//            if( nnzPerColFt[colpos] != 0 )
//               colId[blocksize++] = colpos;
//
//         if( blocksize == 0 )
//            break;
//
//         if( colSparsity )
//            memset(colSparsity, 0, length_col * sizeof(int));
//
//         memset(colsBlockDense, 0, blocksize * length_col * sizeof(double));
//
//         // get column block from Ft (i.e., row block from F)
//         border_right.F.fromGetColsBlock(colId, blocksize, length_col, 0, colsBlockDense, colSparsity);
//
//         solvers_blocked[0]->solve(blocksize, colsBlockDense, colSparsity);
//
//         for( int i = 0; i < blocksize; i++ )
//            colId[i] += (m_res - nF_right - nG_right);
//
//         addLeftBorderTimesDenseColsToResTransp( border_left_transp, colsBlockDense, colId, length_col, blocksize, sparse_res, sym_res, result);
//      }
//   }
//
//   // do we have linking inequality constraints?
//   if( withG )
//   {
//      //                       (G^T)
//      //     SC +=  B^T  K^-1  (0  )
//      //                       (0  )
//
//      SimpleVectorBase<int> nnzPerColGt(nG_right);
//      border_right.G.addNnzPerCol(nnzPerColGt);
//
//      colpos = 0;
//
//      // do block-wise multiplication for columns of G^T part
//      while( colpos < nG_right )
//      {
//         int blocksize = 0;
//
//         for( ; colpos < nG_right && blocksize < blocksizemax; colpos++ )
//            if( nnzPerColGt[colpos] != 0 )
//               colId[blocksize++] = colpos;
//
//         if( blocksize == 0 )
//            break;
//
//         if( colSparsity )
//            memset(colSparsity, 0, length_col * sizeof(int));
//
//         memset(colsBlockDense, 0, blocksize * length_col * sizeof(double));
//
//         border_right.G.fromGetColsBlock(colId, blocksize, length_col, 0, colsBlockDense, colSparsity);
//
//         solvers_blocked[0]->solve(blocksize, colsBlockDense, colSparsity);
//
//         for( int i = 0; i < blocksize; i++ )
//             colId[i] += (m_res - nG_right);
//
//         addLeftBorderTimesDenseColsToResTransp( border_left_transp, colsBlockDense, colId, length_col, blocksize, sparse_res, sym_res, result);
//      }
//   }
//
//#if 0
//   // debug stuff
//   int myrank;
//   static int iteration = 0;
//   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//   ofstream myfile;
//   char filename[50];
//   sprintf(filename, "../blocked_%d_%d.txt", myrank, iteration);
//   myfile.open(filename);
//   iteration++;
//   SC.writeToStream(myfile); // todo write out in each iteration with global counter and MPI rank!
//   myfile.close();
//
//   assert(0);
//#endif
//}

void sLinsys::addBiTLeftKiDenseToResBlockedParallelSolvers( bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
      /* const */ DenseGenMatrix& BT, DoubleMatrix& result)
{
   int m_res, n_res; result.getSize(m_res, n_res);
   assert( m_res >= 0 && n_res >= 0 );
   if( sym_res )
      assert( m_res == n_res );

   int mB, nB; BT.getSize(mB, nB);

   assert( n_solvers >= 1 && n_threads_solvers >= 1 );

   std::vector<int> col_id_cont(blocksizemax * n_solvers);

#ifdef TIME_SCHUR
   const double t_start = omp_get_wtime();
#endif

   //
   //     res +=  B^T  K^-1 B
   //
   const int chunks = std::ceil( static_cast<double>(mB) / blocksizemax );

   #pragma omp parallel for schedule(dynamic, 1) num_threads(n_solvers)
   for( int i = 0; i < chunks; i++ )
   {
      omp_set_num_threads(n_threads_solvers);

      const int actual_blocksize = std::min( (i + 1) * blocksizemax, mB) - i * blocksizemax;

      const int id = omp_get_thread_num();

      int* colId_loc = col_id_cont.data() + id * blocksizemax;

      for( int j = 0; j < actual_blocksize; ++j )
            colId_loc[j] = j + i * blocksizemax;

      double* colsBlockDense_loc = BT[ i * blocksizemax ];

      solvers_blocked[id]->solve(actual_blocksize, colsBlockDense_loc, nullptr);

      addLeftBorderTimesDenseColsToResTransp(border_left_transp, colsBlockDense_loc, colId_loc, nB, actual_blocksize, sparse_res, sym_res, result);
   }

#ifdef TIME_SCHUR
   const double t_end = omp_get_wtime();
   std::cout << "t_end - t_start:" << (t_end - t_start) << std::endl;
   assert(0);
#endif
}

/* res is in transposed form here */
void sLinsys::addBiTLeftKiBiRightToResBlockedParallelSolvers( bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
      /* const */ BorderBiBlock& border_right, DoubleMatrix& result)
{
   if( sparse_res )
      assert( sym_res );

   int m_res_tp, n_res_tp;
   result.getSize(n_res_tp, m_res_tp);
   if( sym_res )
      assert( m_res_tp == n_res_tp );

   int mF_right, nF_right; border_right.F.getSize(mF_right, nF_right);
   int mG_right, nG_right; border_right.G.getSize(mG_right, nG_right);

   const bool with_RAC = border_right.has_RAC;
   const bool withF = ( nF_right > 0 );
   const bool withG = ( nG_right > 0 );

   const int length_col = dynamic_cast<SparseSymMatrix&>(*kkt).size();

#ifndef NDEBUG
   assert( n_solvers >= 1 && n_threads_solvers >= 1 );
   int mF_left, nF_left; border_left_transp.F.getSize(mF_left, nF_left);
   int mG_left, nG_left; border_left_transp.G.getSize(mG_left, nG_left);

   if( border_left_transp.has_RAC )
   {
      int mR_left, nR_left; border_left_transp.R.getSize(mR_left, nR_left);
      int mA_left, nA_left; border_left_transp.A.getSize(mA_left, nA_left);
      int mC_left, nC_left; border_left_transp.C.getSize(mC_left, nC_left);
      assert( mR_left == mA_left );
      assert( mR_left == mC_left );
      assert( nR_left == nF_left );
      assert( nR_left == nG_left );

      assert( mR_left + mF_left + mG_left <= m_res_tp);
      assert( nR_left + nA_left + nC_left == length_col );
   }
   else
   {
      assert( nF_left == nG_left );
      assert( mF_left <= length_col );
      assert( mF_left + mG_left == m_res_tp );
   }

   if( with_RAC )
   {
      int mR_right, nR_right; border_right.R.getSize(mR_right, nR_right);
      int mA_right, nA_right; border_right.A.getSize(mA_right, nA_right);
      int mC_right, nC_right; border_right.C.getSize(mC_right, nC_right);
      assert( nR_right == nA_right);
      assert( nR_right == nC_right);
      assert( mR_right == mF_right );
      assert( mR_right == mG_right );

      assert( nR_right + nF_right + nG_right <= n_res_tp );
      assert( mR_right + mA_right + mC_right == length_col);
   }
   else
   {
      assert( mF_right == mG_right );
      assert( mF_right <= length_col );
      assert( nF_right + nG_right == n_res_tp );
   }
#endif


   if( !colsBlockDense )
      colsBlockDense = new double[blocksizemax * length_col * n_solvers];

   if( !colId )
      colId = new int[blocksizemax * n_solvers];

   // indicating whether a right hand side is zero - deactivated since problems in Pardiso
#if 0
   if( !colSparsity )
      colSparsity = new int[length_col * blocksizemax * n_solvers];
#else
   colSparsity = nullptr;
#endif

#ifdef TIME_SCHUR
   const double t_start = omp_get_wtime();
#endif

   if( with_RAC )
   {
      //                       (R)
      //     SC +=  B^T  K^-1  (A)
      //                       (C)

      int nR_r, mR_r;
      border_right.R.getSize(mR_r, nR_r);
      int nA_r, mA_r;
      border_right.A.getSize(mA_r, nA_r);

      SimpleVectorBase<int> nnzPerColRAC(nR_r);

      border_right.R.addNnzPerCol(nnzPerColRAC);
      border_right.A.addNnzPerCol(nnzPerColRAC);
      border_right.C.addNnzPerCol(nnzPerColRAC);

      const int chunks_RAC = std::ceil( static_cast<double>(nR_r) / blocksizemax );

      #pragma omp parallel for schedule(dynamic, 1) num_threads(n_solvers)
      for( int i = 0; i < chunks_RAC; i++ )
      {
         omp_set_num_threads(n_threads_solvers);

         const int actual_blocksize = std::min( (i + 1) * blocksizemax, nR_r) - i * blocksizemax;

         int nrhs = 0;
         const int id = omp_get_thread_num();

         int * colId_loc = colId + id * blocksizemax;

         for( int j = 0; j < actual_blocksize; ++j )
            if( nnzPerColRAC[j + i * blocksizemax] != 0 )
               colId_loc[nrhs++] = j + i * blocksizemax;

         if( nrhs == 0 )
            continue;

         double* colsBlockDense_loc = colsBlockDense + id * length_col * blocksizemax;
         memset(colsBlockDense_loc, 0, blocksizemax * length_col * sizeof(double));

         int* colSparsity_loc = nullptr;
         if( colSparsity )
            colSparsity_loc = colSparsity + id * length_col * blocksizemax;

         border_right.R.fromGetColsBlock(colId_loc, nrhs, length_col, 0, colsBlockDense_loc, colSparsity_loc);
         border_right.A.fromGetColsBlock(colId_loc, nrhs, length_col, mR_r, colsBlockDense_loc, colSparsity_loc);
         border_right.C.fromGetColsBlock(colId_loc, nrhs, length_col, (mR_r + mA_r), colsBlockDense_loc, colSparsity_loc);

         solvers_blocked[id]->solve(nrhs, colsBlockDense_loc, colSparsity_loc);

         addLeftBorderTimesDenseColsToResTransp(border_left_transp, colsBlockDense_loc, colId_loc, length_col, nrhs, sparse_res, sym_res, result);
      }
   }

#ifdef TIME_SCHUR
   const double t_end = omp_get_wtime();
   std::cout << "t_end - t_start:" << (t_end - t_start) << std::endl;
   assert(0);
#endif

   // do we have linking equality constraints?
   if( withF )
   {
      //                       (F^T)
      //     SC +=  B^T  K^-1  (0  )
      //                       (0  )

      SimpleVectorBase<int> nnzPerColFt(nF_right);
      border_right.F.addNnzPerCol(nnzPerColFt);

      const int chunks_F = std::ceil( static_cast<double>(nF_right) / blocksizemax );

      // do block-wise multiplication for columns of F^T part
      #pragma omp parallel for schedule(dynamic, 1) num_threads(n_solvers)
      for(int i = 0; i < chunks_F; ++i )
      {
         omp_set_num_threads(n_threads_solvers);

         const int actual_blocksize = std::min( (i + 1) * blocksizemax, nF_right) - i * blocksizemax;

         int nrhs = 0;
         const int id = omp_get_thread_num();

         int * colId_loc = colId + id * blocksizemax;

         for( int j = 0; j < actual_blocksize; ++j )
            if( nnzPerColFt[j + i * blocksizemax] != 0 )
               colId_loc[nrhs++] = j + i * blocksizemax;

         if( nrhs == 0 )
            continue;

         double* colsBlockDense_loc = colsBlockDense + id * length_col * blocksizemax;
         memset(colsBlockDense_loc, 0, blocksizemax * length_col * sizeof(double));

         int* colSparsity_loc = nullptr;
         if( colSparsity )
            colSparsity_loc = colSparsity + id * length_col * blocksizemax;

         // get column block from Ft (i.e., row block from F)
         border_right.F.fromGetColsBlock(colId_loc, nrhs, length_col, 0, colsBlockDense_loc, colSparsity_loc);

         solvers_blocked[id]->solve(nrhs, colsBlockDense_loc, colSparsity_loc);

         for( int j = 0; j < nrhs; ++j )
            colId_loc[j] += n_res_tp - nF_right - nG_right;

         addLeftBorderTimesDenseColsToResTransp(border_left_transp, colsBlockDense_loc, colId_loc, length_col, nrhs, sparse_res, sym_res, result);
      }
   }

   // do we have linking inequality constraints?
   if( withG )
   {
      //                       (G^T)
      //     SC +=  B^T  K^-1  (0  )
      //                       (0  )

      SimpleVectorBase<int> nnzPerColGt(nG_right);
      border_right.G.addNnzPerCol( nnzPerColGt );

      const int chunks_G = std::ceil( static_cast<double>(nG_right) / blocksizemax );

      // do block-wise multiplication for columns of G^T part
      #pragma omp parallel for schedule(dynamic, 1) num_threads(n_solvers)
      for( int i = 0; i < chunks_G; ++i )
      {
         omp_set_num_threads(n_threads_solvers);

         const int actual_blocksize = std::min( (i + 1) * blocksizemax, nG_right) - i * blocksizemax;

         int nrhs = 0;
         const int id = omp_get_thread_num();

         int * colId_loc = colId + id * blocksizemax;

         for( int j = 0; j < actual_blocksize; ++j )
            if( nnzPerColGt[j + i * blocksizemax] != 0 )
               colId_loc[nrhs++] = j + i * blocksizemax;

         if( nrhs == 0 )
            continue;

         double* colsBlockDense_loc = colsBlockDense + id * length_col * blocksizemax;
         memset(colsBlockDense_loc, 0, blocksizemax * length_col * sizeof(double));

         int* colSparsity_loc = nullptr;
         if( colSparsity )
            colSparsity_loc = colSparsity + id * length_col * blocksizemax;

         border_right.G.fromGetColsBlock(colId_loc, nrhs, length_col, 0, colsBlockDense_loc, colSparsity_loc);

         solvers_blocked[id]->solve(nrhs, colsBlockDense_loc, colSparsity_loc);

         for( int j = 0; j < nrhs; ++j )
            colId_loc[j] += n_res_tp - nG_right;

         addLeftBorderTimesDenseColsToResTransp(border_left_transp, colsBlockDense_loc, colId_loc, length_col, nrhs, sparse_res, sym_res, result);
      }
   }

#if 0
   // debug stuff
   const int myrank = PIPS_MPIgetRank();

   if( myrank == 0 )
   {
      static int iteration = 0;
      ofstream myfile;
      char filename[50];
      sprintf(filename, "../blocked_%d_%d.txt", myrank, iteration);
      myfile.open(filename);
      iteration++;
      SC.writeToStream(myfile); // todo write out in each iteration with global counter and MPI rank!
      myfile.close();

      assert(0);
   }
#endif
}

void sLinsys::addLeftBorderTimesDenseColsToResTranspSparse( const BorderBiBlock& Bl, const double* cols,
      const int* cols_id, int length_col, int n_cols, SparseSymMatrix& res ) const
{
   /*                  [ R A C ]
    * compute res^T += [ 0 0 0 ] * cols = border_left * cols
    *                  [ F 0 0 ]
    *                  [ G 0 0 ]
    *  the size of the zero rows in border_left is determined by res and can be zero
    */
   int mF, nF; Bl.F.getSize(mF, nF);
   int mG, nG; Bl.G.getSize(mG, nG);
   int mRes, nRes; res.getSize(mRes, nRes);

   const bool with_RAC = Bl.has_RAC;
   const bool with_F = mF > 0;
   const bool with_G = mG > 0;

#ifndef NDEBUG
   assert( mRes == nRes );
   if( with_RAC )
   {
      int mR, nR; Bl.R.getSize(mR, nR);
      int mA, nA; Bl.A.getSize(mA, nA);
      int mC, nC; Bl.C.getSize(mC, nC);
      assert( nF == nG && nF == nR );
      assert( length_col == nR + nA + nC );
      assert( nRes >= mR + mF + mG );
      assert( mR == mA && mA == mC );
   }
   else
      assert( nRes == mF + mG );
#endif

   // multiply each column with left_border and add if to res
   // todo: #pragma omp parallel for schedule(dynamic, 10)
   for( int it_col = 0; it_col < n_cols; it_col++ )
   {
      const double* const col = &cols[it_col * length_col];
      const int row_res = cols_id[it_col];

      assert( row_res < mRes );
      if( with_RAC )
      {
         int mR, nR; Bl.R.getSize(mR, nR);
         int mA, nA; Bl.A.getSize(mA, nA);

         Bl.R.multMatSymUpper(1.0, res, -1.0, &col[0], row_res, 0);
         Bl.A.multMatSymUpper(1.0, res, -1.0, &col[nR], row_res, 0);
         Bl.C.multMatSymUpper(1.0, res, -1.0, &col[nR + nA], row_res, 0);
      }

      if( with_F )
         Bl.F.multMatSymUpper(1.0, res, -1.0, &col[0], row_res, nRes - mF - mG);

      if( with_G )
         Bl.G.multMatSymUpper(1.0, res, -1.0, &col[0], row_res, nRes - mG);
   }
}

void sLinsys::addLeftBorderTimesDenseColsToResTranspDense( const BorderBiBlock& Bl, const double* cols,
      const int* cols_id, int length_col, int n_cols, int n_cols_res, double** res) const
{
   /*                  [ R A C ]
    * compute res^T += [ 0 0 0 ] * colsBlockDense = border_left * colsBlockDense
    *                  [ F 0 0 ]
    *                  [ G 0 0 ]
    *
    *  cols lie as rows in storage
    *  the size of the zero rows in border_left is determined by res and can be zero
    */
   int mF, nF; Bl.F.getSize(mF, nF);
   int mG, nG; Bl.G.getSize(mG, nG);

   const bool with_RAC = Bl.has_RAC;
   const bool with_F = mF > 0;
   const bool with_G = mG > 0;

#ifndef NDEBUG
   if( with_RAC )
   {
      int mR, nR; Bl.R.getSize(mR, nR);
      int mA, nA; Bl.A.getSize(mA, nA);
      int mC, nC; Bl.C.getSize(mC, nC);
      assert( mR == mA && mA == mC );
      assert( nF == nG && nF == nR );
      assert( length_col == nR + nA + nC );
      assert( n_cols_res >= mR + mF + mG );
   }
   else
      assert( n_cols_res >= mF + mG );

   assert( n_cols_res >= 1 );
#endif

   // multiply each column with left factor of SC todo add OMP
   for( int it_col = 0; it_col < n_cols; it_col++ )
   {
      const double* const col = &cols[it_col * length_col];
      const int row_res = cols_id[it_col];

      if( with_RAC )
      {
         int mR, nR; Bl.R.getSize(mR, nR);
         int mA, nA; Bl.A.getSize(mA, nA);
         Bl.R.mult(1.0, &res[row_res][0], 1, -1.0, &col[0], 1);
         Bl.A.mult(1.0, &res[row_res][0], 1, -1.0, &col[nR], 1);
         Bl.C.mult(1.0, &res[row_res][0], 1, -1.0, &col[nR + nA], 1);
      }

      if( with_F )
         Bl.F.mult(1.0, &res[row_res][n_cols_res - mF - mG], 1, -1.0, &col[0], 1);
      if( with_G )
         Bl.G.mult(1.0, &res[row_res][n_cols_res - mG], 1, -1.0, &col[0], 1);
   }
}

void sLinsys::addLeftBorderTimesDenseColsToResTransp( const BorderBiBlock& border_left, const double* cols,
      const int* cols_id, int length_col, int blocksize, bool sparse_res, bool sym_res, DoubleMatrix& res ) const
{
   assert(cols_id && cols);

   if( sparse_res )
   {
      assert( sym_res );
      SparseSymMatrix& res_sparse = dynamic_cast<SparseSymMatrix&>(res);
      addLeftBorderTimesDenseColsToResTranspSparse(border_left, cols, cols_id, length_col, blocksize, res_sparse);
   }
   else
   {
      double** res_array;
      int res_ncols;

#ifndef NDEBUG
      int res_mrows;
      if( sym_res )
         res_mrows = dynamic_cast<DenseSymMatrix&>(res).size();
      else
         res_mrows = dynamic_cast<DenseGenMatrix&>(res).mStorage->m;
      for( int i = 0; i < blocksize; ++i )
         assert( cols_id[i] < res_mrows );
#endif

      if( sym_res )
      {
         DenseSymMatrix& res_dense = dynamic_cast<DenseSymMatrix&>(res);
         res_array = res_dense.mStorage->M;
         res_ncols = res_dense.size();
      }
      else
      {
         DenseGenMatrix& res_dense = dynamic_cast<DenseGenMatrix&>(res);
         res_array = res_dense.mStorage->M;
         res_ncols = res_dense.mStorage->n;
      }

      addLeftBorderTimesDenseColsToResTranspDense(border_left, cols, cols_id, length_col, blocksize, res_ncols, res_array);
   }

}
