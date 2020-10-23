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

extern int gOuterIterRefin;

sLinsys::sLinsys(sFactory* factory_, sData* prob, bool is_hierarchy_root)
  : QpGenLinsys(), kkt(nullptr), solver(nullptr), nThreads(PIPSgetnOMPthreads()),
        blocksizemax( pips_options::getIntParameter("SC_BLOCKWISE_BLOCKSIZE_MAX") ),
        colsBlockDense(nullptr), colId(nullptr), colSparsity(nullptr), is_hierarchy_root(is_hierarchy_root)
{
  prob->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);
  factory = factory_;

#ifdef PARDISO_BLOCKSC
  computeBlockwiseSC = true;
#else
  computeBlockwiseSC = false;
#endif
  // compute schur complement blcokwise when neiher MUMPS nor PARDISO are available
#if !defined(WTIH_MUMPS_ROOT) && !defined(WITH_PARDISO)
   computeBlockwiseSC = true;
#endif

   if( computeBlockwiseSC )
      if( PIPS_MPIgetRank() == 0 )
         std::cout << "Using " << nThreads << " solvers parallely for blockwise SC computation" << std::endl;

  nx = prob->nx; my = prob->my; mz = prob->mz;
  ixlow = prob->ixlow;
  ixupp = prob->ixupp;
  iclow = prob->iclow;
  icupp = prob->icupp;

  nxlow = prob->nxlow;
  nxupp = prob->nxupp;
  mclow = prob->mclow;
  mcupp = prob->mcupp;

  dd = factory_->tree->newPrimalVector();
  assert( dd != nullptr );
  
  dq = factory_->tree->newPrimalVector();
  assert(dq != nullptr);
  prob->getDiagonalOfQ( *dq );

  nomegaInv = factory_->tree->newDualZVector();
  rhs = factory_->tree->newRhs();

  //get the communicator from one of the vectors
  StochVector& dds = dynamic_cast<StochVector&>(*dd);
  this->mpiComm = dds.mpiComm;
  this->iAmDistrib = dds.iAmDistrib;

  useRefs=0;
  data = prob;
  stochNode = factory_->tree;
}

sLinsys::sLinsys(sFactory* factory_,
		 sData* prob,				    
		 OoqpVector* dd_, 
		 OoqpVector* dq_,
		 OoqpVector* nomegaInv_,
		 OoqpVector* rhs_)
  : QpGenLinsys(), kkt(nullptr), solver(nullptr), nThreads(PIPSgetnOMPthreads()),
    blocksizemax( pips_options::getIntParameter("SC_BLOCKWISE_BLOCKSIZE_MAX") ),
    colsBlockDense(nullptr), colId(nullptr), colSparsity(nullptr), is_hierarchy_root(false)
{
  assert( prob );

  prob->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);
  factory = factory_;
#ifdef PARDISO_BLOCKSC
  computeBlockwiseSC = true;
#else
  computeBlockwiseSC = false;
#endif
  // compute schur complement blockkwise when neiher MUMPS nor PARDISO are available
#if !defined(WTIH_MUMPS_ROOT) && !defined(WITH_PARDISO)
   computeBlockwiseSC = true;
#endif

   if( computeBlockwiseSC )
      if( PIPS_MPIgetRank() == 0 )
         std::cout << "Using " << nThreads << " solvers parallely for blockwise SC computation" << std::endl;

  nx = prob->nx; my = prob->my; mz = prob->mz;
  ixlow = prob->ixlow;
  ixupp = prob->ixupp;
  iclow = prob->iclow;
  icupp = prob->icupp;

  nxlow = prob->nxlow;
  nxupp = prob->nxupp;
  mclow = prob->mclow;
  mcupp = prob->mcupp;

  dd = dd_;
  dq = dq_;
  nomegaInv = nomegaInv_;
  rhs = rhs_;

  if( dd )
  {
     StochVector& dds = dynamic_cast<StochVector&>(*dd);
     this->mpiComm = dds.mpiComm;
     this->iAmDistrib = dds.iAmDistrib;
  }
  else
  {
     this->mpiComm = MPI_COMM_NULL;
     this->iAmDistrib = false;
  }

  useRefs = 1;
  data = prob;
  stochNode = factory_->tree;
}


sLinsys::~sLinsys()
{
  if( colSparsity ) delete[] colSparsity;
  if( colId ) delete[] colId;
  if( colsBlockDense ) delete[] colsBlockDense;
  if( !computeBlockwiseSC )
  {
     if (solver) delete solver;
  }
  if (kkt)    delete kkt;
}

void sLinsys::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
				OoqpVector& rhs2_in, OoqpVector& rhs3_in )
{
  StochVector& rhs  = dynamic_cast<StochVector&>(rhs_in);
  StochVector& rhs1 = dynamic_cast<StochVector&>(rhs1_in);
  StochVector& rhs2 = dynamic_cast<StochVector&>(rhs2_in);
  StochVector& rhs3 = dynamic_cast<StochVector&>(rhs3_in);

  rhs.jointCopyFromLinkCons(rhs1, rhs2, rhs3);
}

void sLinsys::separateVars( OoqpVector& x_in, OoqpVector& y_in,
				     OoqpVector& z_in, OoqpVector& vars_in )
{
  StochVector& x    = dynamic_cast<StochVector&>(x_in);
  StochVector& y    = dynamic_cast<StochVector&>(y_in);
  StochVector& z    = dynamic_cast<StochVector&>(z_in);
  StochVector& vars = dynamic_cast<StochVector&>(vars_in);

  vars.jointCopyToLinkCons(x, y, z);
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
  factor2(prob, vars);

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
    *U = new DenseGenMatrix(locnx+locmy+locmz, n0);
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
      *V = new DenseGenMatrix(locnx+locmy+locmz, n0);
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

  solver->Dsolve (zi);  
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

/** sum up right hand side for (current) scenario i and add it to right hand side of scenario 0 */
void sLinsys::addLniziLinkCons(sData *prob, OoqpVector& z0_, OoqpVector& zi_, int parentmy, int parentmz)
{
  SimpleVector& z0 = dynamic_cast<SimpleVector&>(z0_);
  SimpleVector& zi = dynamic_cast<SimpleVector&>(zi_);

  solver->Dsolve(zi);
//  solver->Ltsolve(zi); -> empty

  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();

  int dummy, nx0;
  A.getSize(dummy, nx0);

  // zi2, zi3 are just references to fragments of zi
  SimpleVector zi1 (&zi[0],           locnx);
  SimpleVector zi2 (&zi[locnx],       locmy);
  SimpleVector zi3 (&zi[locnx+locmy], locmz);
  // same for z01 (only the first n0 entries in the output z0 are computed)
  SimpleVector z01 (&z0[0], nx0);

  R.transMult(1.0, z01, -1.0, zi1);
  A.transMult(1.0, z01, -1.0, zi2);
  C.transMult(1.0, z01, -1.0, zi3);

  if( locmyl > 0 )
  {
    assert(locmyl >= 0);
    const int nxMyMz = z0.length() - locmyl - locmzl;

    assert(nxMyMz == nx0 + parentmy + parentmz);

    SimpleVector z0myl (&z0[nxMyMz], locmyl);
    SparseGenMatrix& F = prob->getLocalF();
    F.mult(1.0, z0myl, -1.0, zi1);
  }

  if( locmzl > 0 )
  {
    assert(locmyl >= 0);
    const int nxMyMzMyl = z0.length() - locmzl;

    assert(nxMyMzMyl == nx0 + parentmy + parentmz + locmyl);

    SimpleVector z0mzl (&z0[nxMyMzMyl], locmzl);
    SparseGenMatrix& G = prob->getLocalG();
    G.mult(1.0, z0mzl, -1.0, zi1);
  }
}

void sLinsys::addLniZiHierarchyBorder( DenseGenMatrix& result, BorderLinsys& border )
{
   assert( border.R.children.size() == 0 );
   assert( border.A.mat );
   assert( border.C.mat );
   assert( border.F.mat );
   assert( border.G.mat );
   assert( border.R.mat );

   const bool result_sparse = false;
   const bool result_sym = false;

   BorderBiBlock border_right( *border.R.mat, *border.A.mat, *border.C.mat, border.F.mat->getTranspose(), border.G.mat->getTranspose() );
   BorderBiBlock border_left_transp( data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(),
         data->getLocalC().getTranspose(), data->getLocalF(), data->getLocalG() );

   addBiTLeftKiBiRightToResBlocked( result_sparse, result_sym, border_left_transp, border_right, result);
}

/* calculate res += X_i * B_i^T */
void sLinsys::multRightDenseSchurComplBlocked( /*const*/sData* prob, const DenseGenMatrix& X, DenseGenMatrix& result, int parent_nx, int parent_my, int parent_mz )
{
   /*          locnx locmy locmz
    *        [  RiT   AiT   CiT ] pnx
    *        [   0     0     0  ] pmy
    * Bi^T = [   0     0     0  ] pmz
    *        [   Fi    0     0  ] locmyl
    *        [   Gi    0     0  ] locmzl
    */

#ifndef NDEBUG
   int mX, nX; X.getSize(mX, nX);
   int mRes, nRes; result.getSize(mRes, nRes);
   assert( locnx + locmy + locmz == nRes );
   assert( mX == mRes );
   assert( parent_nx + parent_my + parent_mz + locmyl + locmzl == nX );
#endif

   SparseGenMatrix& A = prob->getLocalA();
   SparseGenMatrix& C = prob->getLocalC();
   SparseGenMatrix& F = prob->getLocalF();
   SparseGenMatrix& G = prob->getLocalG();
   SparseGenMatrix& R = prob->getLocalCrossHessian();

#ifndef NDEBUG
   int mR, nR; R.getSize(mR, nR);
   assert( nR == parent_nx );
   assert( mR == locnx );

   int mA, nA; A.getSize( mA, nA );
   assert( nA == parent_nx );
   assert( mA == locmy );

   int mC, nC; C.getSize( mC, nC );
   assert( nC == parent_nx );
   assert( mC == locmz );

   int mF, dummy; F.getSize( mF, dummy );
   assert( locmyl == mF );

   int mG; G.getSize( mG, dummy );
   assert( locmzl == mG );
#endif

   // X from the right with each column of Bi^T todo add OMP to submethods
   /*            [ RiT ]
    *            [  0  ]
    * res += X * [  0  ]
    *            [  Fi ]
    *            [  Gi ]
    */
   X.multMatAt( 0, 1.0, 0, result, 1.0, R.getTranspose() );

   X.multMatAt( parent_nx + parent_my + parent_mz, 1.0, 0, result, 1.0, F );

   X.multMatAt( parent_nx + parent_my + parent_mz + locmyl, 1.0, 0, result, 1.0, G );

   /*            [ AiT ]
    *            [  0  ]
    * res += X * [  0  ]
    *            [  0  ]
    *            [  0  ]
    */
   X.multMatAt( 0, 1.0, locnx, result, 1.0, A.getTranspose() );

   /*            [ CiT ]
    *            [  0  ]
    * res += X * [  0  ]
    *            [  0  ]
    *            [  0  ]
    */
   X.multMatAt( 0, 1.0, locnx + locmy, result, 1.0, C.getTranspose() );
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
   int mRt, nRt; BiT.R.getSize(mRt, nRt);
   int mAt, nAt; BiT.A.getSize(mAt, nAt);
   int mCt, nCt; BiT.C.getSize(mCt, nCt);
   int mF, nF; BiT.F.getSize(mF, nF);
   int mG, nG; BiT.G.getSize(mG, nG);

#ifndef NDEBUG
   int mres, nres; res.getSize(mres, nres);

   assert( mres == mRt + mF + mG );
   assert( nF == nG );
   assert( nF == nRt );
   assert( nRt + nAt + nCt == nres );
#endif

   addMatAt( res, BiT.R, 0, 0 );
   addMatAt( res, BiT.A, 0, nRt );
   addMatAt( res, BiT.C, 0, nRt + nAt );

   addMatAt( res, BiT.F, mRt, 0 );
   addMatAt( res, BiT.G, mRt + mF, 0 );
}

/* compute Bi_{outer}^T X_i = Bi_{outer}^T Ki^-1 (Bi_{outer} - Bi_{inner} X0) and add it to SC */
void sLinsys::LniTransMultHierarchyBorder( DenseSymMatrix& SC, const DenseGenMatrix& X0, BorderLinsys& border, int parent_nx, int parent_my, int parent_mz )
{
   int nx_border, myl_border, mzl_border, dummy;

   border.R.getSize(dummy, nx_border);
   border.F.getSize(myl_border, dummy);
   border.G.getSize(mzl_border, dummy);

   /* buffer for (Bi_{outer} - Bi_{inner} X0)^T = Bi_{outer}^T - X0^T Bi_{inner}^T */
   // TODO : reuse an make member
   DenseGenMatrix* BiT_buffer = new DenseGenMatrix( nx_border + myl_border + mzl_border, locnx + locmy + locmz );

   /* Bi buffer and X0 are in transposed form for memory alignment reasons when solving with K_i */
   int m, n; BiT_buffer->getSize(m, n);
   BiT_buffer->atPutZeros(0, 0, m, n );

   /* put (Bi_{outer})^T into buffer
    *
    *                  nxb myb mzb
    *                [ RiT AiT CiT ]
    * Bi_{outer}^T = [  Fi  0   0  ]
    *                [  Gi  0   0  ]
    */
   assert( border.R.mat );
   assert( border.A.mat );
   assert( border.C.mat );
   assert( border.F.mat );
   assert( border.G.mat );

   const BorderBiBlock BiT_outer( border.R.mat->getTranspose(), border.A.mat->getTranspose(), border.C.mat->getTranspose(), *border.F.mat, *border.G.mat);
   addBiTBorder( *BiT_buffer, BiT_outer);

   /* compute (Bi_{outer} - Bi_{inner} * X0)^T = Bi_{outer}^T - X0^T * Bi_{inner}^T
    *
    *                     [ Ri 0 0 FiT GiT ]^T
    * Bi_{inner} = X0^T * [ Ai 0 0  0   0  ]
    *                     [ Ci 0 0  0   0  ]
    */
   multRightDenseSchurComplBlocked( data, X0, *BiT_buffer, parent_nx, parent_my, parent_mz );

   /* solve blockwise Ki X = Bi_buffer and multiply from left with Bi_{outer} and add to SC */
   addBiTLeftKiDenseToResBlockedParallelSolvers( false, true, BiT_outer, *BiT_buffer, SC );
}

void sLinsys::solveCompressed( OoqpVector& rhs_ )
{
  StochVector& rhs = dynamic_cast<StochVector&>(rhs_);
#ifdef TIMING
  //double tTot=MPI_Wtime();
#endif
  Lsolve (data,rhs); 
  Dsolve (data,rhs);
  Ltsolve(data,rhs);
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
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();
  int N, nx0;

  //get nx(parent) from the number of cols of A (or C). Use N as dummy
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
  SimpleVector LniTx2(&LniTx[locnx], locmy);
  SimpleVector LniTx3(&LniTx[locnx+locmy], locmz);
  
  LniTx1.setToZero();
  R.mult(0.0, LniTx1, 1.0, x1);
  A.mult(0.0, LniTx2, 1.0, x1);
  C.mult(0.0, LniTx3, 1.0, x1);

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

#include "PardisoSolver.h"
/**
 * Computes U = Gi * inv(H_i) * Gi^T.
 *        [ R 0 0 ]
 * Gi^T = [ A 0 0 ]
 *        [ C 0 0]
 *
 * A and C are the recourse eq. and ineq. matrices, R is the cross
 * Hessian term.
 */

/*void sLinsys::addTermToDenseSchurCompl(sData *prob, 
				       DenseSymMatrix& SC) 
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();


  int N, nxP, NP;
  A.getSize(N, nxP); assert(N==locmy);
  NP = SC.size(); assert(NP>=nxP);

  if(nxP==-1) C.getSize(N,nxP);
  if(nxP==-1) nxP = NP;
  N = locnx+locmy+locmz;

  int blocksize = 64;
  DenseGenMatrix cols(blocksize,N);

  bool ispardiso=false;
  PardisoSolver* pardisoSlv = dynamic_cast<PardisoSolver*>(solver);
  int* colSparsity=nullptr;
  if(pardisoSlv) {
    ispardiso=true;
    colSparsity=new int[N];
    //blocksize=32;
  }

  for (int it=0; it < nxP; it += blocksize) {
    int start=it;
    int end = MIN(it+blocksize,nxP);
    int numcols = end-start;
    cols.getStorageRef().m = numcols; // avoid extra solves


    bool allzero = true;
    memset(&cols[0][0],0,N*blocksize*sizeof(double));

    if(ispardiso) {
      for(int i=0; i<N; i++) colSparsity[i]=0;
      R.getStorageRef().fromGetColBlock(start, &cols[0][0], N, numcols, colSparsity, allzero);
      A.getStorageRef().fromGetColBlock(start, &cols[0][locnx], N, numcols, colSparsity, allzero);
      C.getStorageRef().fromGetColBlock(start, &cols[0][locnx+locmy], N, numcols, colSparsity, allzero);

    } else {
      R.getStorageRef().fromGetColBlock(start, &cols[0][0], N, numcols, allzero);
      A.getStorageRef().fromGetColBlock(start, &cols[0][locnx], N, numcols, allzero);
      C.getStorageRef().fromGetColBlock(start, &cols[0][locnx+locmy], N, numcols, allzero);
    }

    if(!allzero) {
      
      if(ispardiso)
	pardisoSlv->solve(cols,colSparsity);
      else 
	solver->solve(cols);

      R.getStorageRef().transMultMat( 1.0, SC[start], numcols, NP,  
       				      -1.0, &cols[0][0], N);
      A.getStorageRef().transMultMat( 1.0, SC[start], numcols, NP,  
       				      -1.0, &cols[0][locnx], N);
      C.getStorageRef().transMultMat( 1.0, SC[start], numcols, NP,
       				      -1.0, &cols[0][locnx+locmy], N);

      // this code seems to have problems
       //R.getStorageRef().transMultMatLower( 1.0, SC[start], numcols, NP,
       //					   -1.0, &cols[0][0], N, start);
       //A.getStorageRef().transMultMatLower( 1.0, SC[start], numcols, NP,
       //					   -1.0, &cols[0][locnx], N, start);
       //C.getStorageRef().transMultMatLower( 1.0, SC[start], numcols, NP,
       //				   -1.0, &cols[0][locnx+locmy], N, start); 

    } //end !allzero
  }

  if(ispardiso) delete[] colSparsity;
}

*/
/* this is the original code that was doing one column at a time. */

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
void sLinsys::addBiTLeftKiBiRightToResBlocked( bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
      /* const */ BorderBiBlock& border_right, DoubleMatrix& result)
{
   int m_res, n_res; result.getSize(m_res, n_res);
   assert( m_res >= 0 && n_res >= 0 );
   if( sym_res )
      assert( m_res == n_res );

   int mR_right, nR_right; border_right.R.getSize(mR_right, nR_right);

   int mA_right, nA_right; border_right.A.getSize(mA_right, nA_right);
   assert( nR_right == nA_right);

   int mC_right, nC_right; border_right.C.getSize(mC_right, nC_right);
   assert( nR_right == nC_right);

   int mF_right, nF_right; border_right.F.getSize(mF_right, nF_right);
   assert( mR_right == mF_right );

   int mG_right, nG_right; border_right.G.getSize(mG_right, nG_right);
   assert( mR_right == mG_right );

   /* we add to the transposed of res */
   assert( nR_right + nF_right + nG_right <= m_res);

   SimpleVectorBase<int> nnzPerColRAC(nR_right);

   border_right.R.addNnzPerCol(nnzPerColRAC);
   border_right.A.addNnzPerCol(nnzPerColRAC);
   border_right.C.addNnzPerCol(nnzPerColRAC);

   const int withF = ( nF_right > 0 );
   const int withG = ( nG_right > 0 );
   const int length_col = mR_right + mA_right + mC_right;

   assert(nThreads >= 1);

   if( colsBlockDense == nullptr )
      colsBlockDense = new double[blocksizemax * length_col];
//   else
//   {
//      delete[] colsBlockDense;
//      colsBlockDense = new double[blocksizemax * length_col];
//   }

   if( colId == nullptr )
      colId = new int[blocksizemax];

   // indicating whether a right hand side is zero - deactivated since problems in Pardiso
#if 0
   if( colSparsity == nullptr )
      colSparsity = new int[length_col * blocksizemax];
#else
   colSparsity = nullptr;
#endif

#ifdef TIME_SCHUR
   const double t_start = omp_get_wtime();
#endif

   int colpos = 0;

   //                       (R)
   //     SC +=  B^T  K^-1  (A)
   //                       (C)
   while( colpos < nR_right )
   {
      int blocksize = 0;

      for( ; colpos < nR_right && blocksize < blocksizemax; colpos++ )
         if( nnzPerColRAC[colpos] != 0 )
            colId[blocksize++] = colpos;

      if( blocksize == 0 )
         break;

      memset(colsBlockDense, 0, blocksize * length_col * sizeof(double));

      if( colSparsity )
         memset(colSparsity, 0, length_col * sizeof(int));

      border_right.R.fromGetColsBlock(colId, blocksize, length_col, 0, colsBlockDense, colSparsity);
      border_right.A.fromGetColsBlock(colId, blocksize, length_col, mR_right, colsBlockDense, colSparsity);
      border_right.C.fromGetColsBlock(colId, blocksize, length_col, (mR_right + mA_right), colsBlockDense, colSparsity);

      solver->solve(blocksize, colsBlockDense, colSparsity);

      addLeftBorderTimesDenseColsToResTransp( border_left_transp, colsBlockDense, colId, length_col, blocksize, sparse_res, sym_res, result);
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

      colpos = 0;
      // do block-wise multiplication for columns of F^T part
      while( colpos < nF_right )
      {
         int blocksize = 0;

         for( ; colpos < nF_right && blocksize < blocksizemax; colpos++ )
            if( nnzPerColFt[colpos] != 0 )
               colId[blocksize++] = colpos;

         if( blocksize == 0 )
            break;

         if( colSparsity )
            memset(colSparsity, 0, length_col * sizeof(int));

         memset(colsBlockDense, 0, blocksize * length_col * sizeof(double));

         // get column block from Ft (i.e., row block from F)
         border_right.F.fromGetColsBlock(colId, blocksize, length_col, 0, colsBlockDense, colSparsity);

         solver->solve(blocksize, colsBlockDense, colSparsity);

         for( int i = 0; i < blocksize; i++ )
            colId[i] += (m_res - nF_right - nG_right);

         addLeftBorderTimesDenseColsToResTransp( border_left_transp, colsBlockDense, colId, length_col, blocksize, sparse_res, sym_res, result);
      }
   }

   // do we have linking inequality constraints?
   if( withG )
   {
      //                       (G^T)
      //     SC +=  B^T  K^-1  (0  )
      //                       (0  )

      SimpleVectorBase<int> nnzPerColGt(nG_right);
      border_right.G.addNnzPerCol(nnzPerColGt);

      colpos = 0;

      // do block-wise multiplication for columns of G^T part
      while( colpos < nG_right )
      {
         int blocksize = 0;

         for( ; colpos < nG_right && blocksize < blocksizemax; colpos++ )
            if( nnzPerColGt[colpos] != 0 )
               colId[blocksize++] = colpos;

         if( blocksize == 0 )
            break;

         if( colSparsity )
            memset(colSparsity, 0, length_col * sizeof(int));

         memset(colsBlockDense, 0, blocksize * length_col * sizeof(double));

         border_right.G.fromGetColsBlock(colId, blocksize, length_col, 0, colsBlockDense, colSparsity);

         solver->solve(blocksize, colsBlockDense, colSparsity);

         for( int i = 0; i < blocksize; i++ )
             colId[i] += (m_res - nG_right);

         addLeftBorderTimesDenseColsToResTransp( border_left_transp, colsBlockDense, colId, length_col, blocksize, sparse_res, sym_res, result);
      }
   }

#if 0
   // debug stuff
   int myrank;
   static int iteration = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   ofstream myfile;
   char filename[50];
   sprintf(filename, "../blocked_%d_%d.txt", myrank, iteration);
   myfile.open(filename);
   iteration++;
   SC.writeToStream(myfile); // todo write out in each iteration with global counter and MPI rank!
   myfile.close();

   assert(0);
#endif
}

void sLinsys::addBiTLeftKiDenseToResBlockedParallelSolvers( bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
      /* const */ DenseGenMatrix& BT, DoubleMatrix& result)
{
   int m_res, n_res; result.getSize(m_res, n_res);
   assert( m_res >= 0 && n_res >= 0 );
   if( sym_res )
      assert( m_res == n_res );

   int mB, nB; BT.getSize(mB, nB);

   assert(nThreads >= 1);

   int * col_id_cont = new int[blocksizemax * n_solvers];

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
      const int actual_blocksize = std::min( (i + 1) * blocksizemax, mB) - i * blocksizemax;

      const int id = omp_get_thread_num();

      int * colId_loc = col_id_cont + id * blocksizemax;

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


void sLinsys::addBiTLeftKiBiRightToResBlockedParallelSolvers( bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
      /* const */ BorderBiBlock& border_right, DoubleMatrix& result)
{
   int m_res, n_res; result.getSize(m_res, n_res);
   assert( m_res >= 0 && n_res >= 0 );
   if( sym_res )
      assert( m_res == n_res );

   int mR_right, nR_right; border_right.R.getSize(mR_right, nR_right);

   int mA_right, nA_right; border_right.A.getSize(mA_right, nA_right);
   assert( nR_right == nA_right);

   int mC_right, nC_right; border_right.C.getSize(mC_right, nC_right);
   assert( nR_right == nC_right);

   int mF_right, nF_right; border_right.F.getSize(mF_right, nF_right);
   assert( mR_right == mF_right );

   int mG_right, nG_right; border_right.G.getSize(mG_right, nG_right);
   assert( mR_right == mG_right );

   assert( nR_right + nF_right + nG_right <= m_res);

   SimpleVectorBase<int> nnzPerColRAC(nR_right);

   border_right.R.addNnzPerCol(nnzPerColRAC);
   border_right.A.addNnzPerCol(nnzPerColRAC);
   border_right.C.addNnzPerCol(nnzPerColRAC);

   const int withF = ( nF_right > 0 );
   const int withG = ( nG_right > 0 );
   const int length_col = mR_right + mA_right + mC_right;

   assert(nThreads >= 1);

   if( colsBlockDense == nullptr )
      colsBlockDense = new double[blocksizemax * length_col * n_solvers];

   if( colId == nullptr )
      colId = new int[blocksizemax * n_solvers];

   // indicating whether a right hand side is zero - deactivated since problems in Pardiso
#if 0
   if( colSparsity == nullptr )
      colSparsity = new int[length_col * blocksizemax * n_solvers];
#else
   colSparsity = nullptr;
#endif

#ifdef TIME_SCHUR
   const double t_start = omp_get_wtime();
#endif


   //                       (R)
   //     SC +=  B^T  K^-1  (A)
   //                       (C)
   const int chunks_RAC = std::ceil( static_cast<double>(nR_right) / blocksizemax );

   #pragma omp parallel for schedule(dynamic, 1) num_threads(n_solvers)
   for( int i = 0; i < chunks_RAC; i++ )
   {
      const int actual_blocksize = std::min( (i + 1) * blocksizemax, nR_right) - i * blocksizemax;

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
      border_right.A.fromGetColsBlock(colId_loc, nrhs, length_col, mR_right, colsBlockDense_loc, colSparsity_loc);
      border_right.C.fromGetColsBlock(colId_loc, nrhs, length_col, (mR_right + mA_right), colsBlockDense_loc, colSparsity_loc);

      solvers_blocked[id]->solve(nrhs, colsBlockDense_loc, colSparsity_loc);

      addLeftBorderTimesDenseColsToResTransp(border_left_transp, colsBlockDense_loc, colId_loc, length_col, nrhs, sparse_res, sym_res, result);
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
            colId_loc[j] += m_res - nF_right - nG_right;

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
            colId_loc[j] += m_res - nG_right;

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

void sLinsys::addLeftBorderTimesDenseColsToResTranspSparse( const BorderBiBlock& border_left, const double* cols,
      const int* cols_id, int length_col, int n_cols, SparseSymMatrix& res) const
{
   /*                  [ R A C ]
    * compute res^T += [ 0 0 0 ] * cols = border_left * cols
    *                  [ F 0 0 ]
    *                  [ G 0 0 ]
    *  the size of the zero rows in border_left is determined by res and can be zero
    */

   int mR, nR; border_left.R.getSize(mR, nR);
   int mA, nA; border_left.A.getSize(mA, nA);
   int mC, nC; border_left.C.getSize(mC, nC);
   int mF, nF; border_left.F.getSize(mF, nF);
   int mG, nG; border_left.G.getSize(mG, nG);
   assert( mR == mA && mA == mC );
   assert( nF == nG && nF == nR );
   assert( length_col == nR + nA + nC );

   int mRes, nRes; res.getSize(mRes, nRes);
   assert( mRes == nRes );
   assert( nRes >= mR + mF + mG );

   // multiply each column with left_border and add if to res
   // todo: #pragma omp parallel for schedule(dynamic, 10)
   for( int it_col = 0; it_col < n_cols; it_col++ )
   {
      const double* const col = &cols[it_col * length_col];
      const int row_res = cols_id[it_col];

      assert( row_res < mRes );
      border_left.R.multMatSymUpper(1.0, res, -1.0, &col[0], row_res, 0);

      border_left.A.multMatSymUpper(1.0, res, -1.0, &col[nR], row_res, 0);

      border_left.C.multMatSymUpper(1.0, res, -1.0, &col[nR + nA], row_res, 0);

      if( mF > 0 )
         border_left.F.multMatSymUpper(1.0, res, -1.0, &col[0], row_res, nRes - mF - mG);

      if( mG > 0 )
         border_left.G.multMatSymUpper(1.0, res, -1.0, &col[0], row_res, nRes - mG);
   }
}

void sLinsys::addLeftBorderTimesDenseColsToResTranspDense( const BorderBiBlock& border_left, const double* cols,
      const int* cols_id, int length_col, int n_cols, int m_rows_res, int n_cols_res, double** res) const
{
   /*                  [ R A C ]
    * compute res^T += [ 0 0 0 ] * colsBlockDense = border_left * colsBlockDense
    *                  [ F 0 0 ]
    *                  [ G 0 0 ]
    *
    *  cols lie as rows in storage
    *  the size of the zero rows in border_left is determined by res and can be zero
    */

   int mR, nR; border_left.R.getSize(mR, nR);
   int mA, nA; border_left.A.getSize(mA, nA);
   int mC, nC; border_left.C.getSize(mC, nC);
   int mF, nF; border_left.F.getSize(mF, nF);
   int mG, nG; border_left.G.getSize(mG, nG);
   assert( mR == mA && mA == mC );
   assert( nF == nG && nF == nR );
   assert( length_col == nR + nA + nC );

   assert( n_cols_res >= 1 );
   assert( n_cols_res >= mR + mF + mG );

   // multiply each column with left factor of SC todo add OMP
   for( int it_col = 0; it_col < n_cols; it_col++ )
   {
      const double* const col = &cols[it_col * length_col];
      const int row_res = cols_id[it_col];
      assert( row_res < m_rows_res );

      border_left.R.mult(1.0, &res[row_res][0], 1, -1.0, &col[0], 1);

      border_left.A.mult(1.0, &res[row_res][0], 1, -1.0, &col[nR], 1);

      border_left.C.mult(1.0, &res[row_res][0], 1, -1.0, &col[nR + nA], 1);

      border_left.F.mult(1.0, &res[row_res][n_cols_res - mF - mG], 1, -1.0, &col[0], 1);

      border_left.G.mult(1.0, &res[row_res][n_cols_res - mG], 1, -1.0, &col[0], 1);
   }
}

void sLinsys::addLeftBorderTimesDenseColsToResTransp( const BorderBiBlock& border_left, const double* cols,
      const int* cols_id, int length_col, int blocksize, bool sparse_res, bool sym_res, DoubleMatrix& res) const
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
      int res_mrows;

      if( sym_res )
      {
         DenseSymMatrix& res_dense = dynamic_cast<DenseSymMatrix&>(res);
         res_array = res_dense.mStorage->M;
         res_ncols = res_dense.size();
         res_mrows = res_dense.size();
      }
      else
      {
         DenseGenMatrix& res_dense = dynamic_cast<DenseGenMatrix&>(res);
         res_array = res_dense.mStorage->M;
         res_ncols = res_dense.mStorage->n;
         res_mrows = res_dense.mStorage->m;
      }

      addLeftBorderTimesDenseColsToResTranspDense(border_left, cols, cols_id, length_col, blocksize, res_mrows, res_ncols, res_array);
   }

}
