/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysLeaf.h"


sLinsysLeaf::sLinsysLeaf(sFactory *factory_, sData* prob,
          OoqpVector* dd_,
          OoqpVector* dq_,
          OoqpVector* nomegaInv_,
          OoqpVector* rhs_
       )
  : sLinsys(factory_, prob, dd_, dq_, nomegaInv_, rhs_, false)
{
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

  if( computeBlockwiseSC )
  {
     const int n_omp_threads = PIPSgetnOMPthreads();
     const int n_threads_for_pardiso = 2;
     if( pips_options::getIntParameter("LINEAR_LEAF_SOLVER") == SolverType::SOLVER_PARDISO )
     {
        n_solvers = std::max( 1, n_omp_threads / n_threads_for_pardiso );
        n_threads_solvers = ( n_omp_threads > 1 ) ? 2 : 1;
     }
     else
     {
        n_solvers = n_omp_threads;
        n_threads_solvers = 1;
     }

     static bool printed = false;
     if( PIPS_MPIgetRank() == 0 && !printed )
     {
        printed = true;
        std::cout << "Using " << n_solvers << " solvers in parallel (with "
           << n_threads_solvers << " threads each) for leaf SC computation - sLinsysLeaf\n";
     }
  }

  kkt.reset(kkt_sp);
  solver.reset( factory_->newLeafSolver( kkt_sp ) );

#ifdef TIMING
  const double t1 = MPI_Wtime() - t0;
  if (myRank == 0) printf("Rank 0: new sLinsysLeaf took %f sec\n",t1);
#endif

  mpiComm = (dynamic_cast<StochVector*>(dd_))->mpiComm;
}

void sLinsysLeaf::factor2(sData*, Variables*)
{
   // Diagonals were already updated, so
   // just trigger a local refactorization (if needed, depends on the type of lin solver).
   stochNode->resMon.recFactTmLocal_start();

   solver->matrixChanged();

   stochNode->resMon.recFactTmLocal_stop();
}

void sLinsysLeaf::putXDiagonal( OoqpVector& xdiag_ )
{
  StochVector& xdiag = dynamic_cast<StochVector&>(xdiag_);
  kkt->atPutDiagonal( 0, *xdiag.vec );
}

void sLinsysLeaf::putZDiagonal( OoqpVector& zdiag_)
{
  StochVector& zdiag = dynamic_cast<StochVector&>(zdiag_);
  kkt->atPutDiagonal( locnx + locmy, *zdiag.vec );
}

void sLinsysLeaf::Dsolve( sData*, OoqpVector& x_in )
{
   StochVector& x = dynamic_cast<StochVector&>(x_in);
   assert(x.children.size()==0);
   stochNode->resMon.recDsolveTmChildren_start();
   solver->Dsolve(*x.vec);
   stochNode->resMon.recDsolveTmChildren_stop();
}

void sLinsysLeaf::Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp, bool)
{
   StochVector& b = dynamic_cast<StochVector&>(x);
   SimpleVector& bi = dynamic_cast<SimpleVector&>(*b.vec);
   assert( 0 == b.children.size() );

#ifdef TIMING
   stochNode->resMon.eLtsolve.clear();
   stochNode->resMon.recLtsolveTmLocal_start();
#endif

   //b_i -= Lni^T x0
   LniTransMult(prob, bi, -1.0, xp);
   //  solver->Ltsolve(bi); -> empty
#ifdef TIMING
   stochNode->resMon.recLtsolveTmChildren_stop();
#endif
}

void sLinsysLeaf::deleteChildren()
{ }

/** sum up right hand side for (current) scenario i and add it to right hand side of scenario 0 */
void sLinsysLeaf::addLniziLinkCons(sData *prob, OoqpVector& z0_, OoqpVector& zi_, bool /*use_local_RAC*/ )
{
  SimpleVector& z0 = dynamic_cast<SimpleVector&>(z0_);
  SimpleVector& zi = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(zi_).vec);

  solver->solve(zi);

  int dummy{0};
  int nx0{0};

  SparseGenMatrix& A = prob->getLocalA();
  if( data->hasRAC() )
     A.getSize(dummy, nx0);

  SimpleVector z01 (&z0[0], nx0);
  SimpleVector zi1 (&zi[0], locnx);

  if( data->hasRAC() )
  {
     SparseGenMatrix& C = prob->getLocalC();
     SparseGenMatrix& R = prob->getLocalCrossHessian();

     SimpleVector zi2 (&zi[locnx], locmy);
     SimpleVector zi3 (&zi[locnx+locmy], locmz);

     R.transMult(1.0, z01, -1.0, zi1);
     A.transMult(1.0, z01, -1.0, zi2);
     C.transMult(1.0, z01, -1.0, zi3);
  }

  if( locmyl > 0 )
  {
    assert(locmyl >= 0);
    const int nxMyMz = z0.length() - locmyl - locmzl;

    SimpleVector z0myl (&z0[nxMyMz], locmyl);
    SparseGenMatrix& F = prob->getLocalF();
    F.mult(1.0, z0myl, -1.0, zi1);
  }

  if( locmzl > 0 )
  {
    assert(locmyl >= 0);
    const int nxMyMzMyl = z0.length() - locmzl;

    SimpleVector z0mzl (&z0[nxMyMzMyl], locmzl);
    SparseGenMatrix& G = prob->getLocalG();
    G.mult(1.0, z0mzl, -1.0, zi1);
  }
}

/* compute Bli^T X_i = Bli^T Ki^-1 (Bri - Br_mod_border - Bi_{inner} X0) and add it to res */
void sLinsysLeaf::LniTransMultHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br,
      std::vector<BorderMod>& Br_mod_border, bool sparse_res, bool sym_res, bool )
{
   int mres, nres; res.getSize(mres, nres);

#ifndef NDEBUG
   int nx_border, myl_border, mzl_border, dummy;
   if( Br.has_RAC )
      Br.R.getSize(dummy, nx_border);
   else if( Br.use_local_RAC )
      data->getLocalCrossHessian().getSize(dummy, nx_border);
   else
      nx_border = 0;
   Br.F.getSize(myl_border, dummy);
   Br.G.getSize(mzl_border, dummy);
   assert( nx_border + myl_border + mzl_border <= mres );
#endif
   /* buffer for (Bri - (sum_j Brmodj * Xmodj)_i - Bi_{inner} X0)^T = Bri^T - X0^T Bi_{inner}^T - (sum_j Xmodj^T Brmodj^T)_i */
   // TODO : reuse and make member ? possible?

   std::unique_ptr<DenseGenMatrix> BiT_buffer( new DenseGenMatrix( mres, dynamic_cast<SparseSymMatrix&>(*kkt).size() ) );

   /* Bi buffer and X0 are in transposed form for memory alignment reasons when solving with K_i */
   int m, n; BiT_buffer->getSize(m, n);
   BiT_buffer->atPutZeros(0, 0, m, n );

   /* put (Bri)^T into buffer
    *
    *                [ RiT AiT CiT ]
    *                [  0   0   0  ]
    * Bri^T        = [  Fi  0   0  ]
    *                [  Gi  0   0  ]
    */
   std::unique_ptr<BorderBiBlock> BriT{};

   if( Br.has_RAC )
      BriT.reset( new BorderBiBlock( dynamic_cast<SparseGenMatrix&>(*Br.R.mat).getTranspose(),
            dynamic_cast<SparseGenMatrix&>(*Br.A.mat).getTranspose(), dynamic_cast<SparseGenMatrix&>(*Br.C.mat).getTranspose(),
            dynamic_cast<SparseGenMatrix&>(*Br.F.mat), dynamic_cast<SparseGenMatrix&>(*Br.G.mat) ) );
   else if( Br.use_local_RAC )
      BriT.reset( new BorderBiBlock( data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(), data->getLocalC().getTranspose(),
            dynamic_cast<SparseGenMatrix&>(*Br.F.mat), dynamic_cast<SparseGenMatrix&>(*Br.G.mat) ) );
   else
      BriT.reset( new BorderBiBlock( dynamic_cast<SparseGenMatrix&>(*Br.F.mat), dynamic_cast<SparseGenMatrix&>(*Br.G.mat), false ) );
   addBiTBorder( *BiT_buffer, *BriT);

   /* compute (Bri - Bi_{inner} * X0)^T = Bri^T - X0^T * Bi_{inner}^T
    *
    *                     [ Ri 0 0 FiT GiT ]^T
    * Bi_{inner} = X0^T * [ Ai 0 0  0   0  ]
    *                     [ Ci 0 0  0   0  ]
    */
   BorderBiBlock BiT_inner = data->hasRAC() ? BorderBiBlock( data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(), data->getLocalC().getTranspose(),
         data->getLocalF(), data->getLocalG() ) : BorderBiBlock( data->getLocalF(), data->getLocalG(), false );
   multRightDenseBorderBlocked( BiT_inner, X0, *BiT_buffer );

   /* now similarly compute BiT_buffer += X_j^T Bmodj for all j */
   multRightDenseBorderModBlocked( Br_mod_border, *BiT_buffer );

   /* compute Bli^T Ki^-1 Bi_buffe = Bli^T Ki^-1 (Bri^T - X0^T * Bi_{inner}^T - sumj Xj^T Bmodj^T) */
   std::unique_ptr<BorderBiBlock> BliT{};
   if( Bl.has_RAC )
      BliT.reset( new BorderBiBlock ( dynamic_cast<SparseGenMatrix&>(*Bl.R.mat).getTranspose(),
         dynamic_cast<SparseGenMatrix&>(*Bl.A.mat).getTranspose(), dynamic_cast<SparseGenMatrix&>(*Bl.C.mat).getTranspose(),
         dynamic_cast<SparseGenMatrix&>(*Bl.F.mat), dynamic_cast<SparseGenMatrix&>(*Bl.G.mat) ) );
   else if( Bl.use_local_RAC )
      BliT.reset( new BorderBiBlock( data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(), data->getLocalC().getTranspose(),
         dynamic_cast<SparseGenMatrix&>(*Bl.F.mat), dynamic_cast<SparseGenMatrix&>(*Bl.G.mat) ) );
   else
      BliT.reset( new BorderBiBlock( dynamic_cast<SparseGenMatrix&>(*Bl.F.mat), dynamic_cast<SparseGenMatrix&>(*Bl.G.mat), false ) );

   addBiTLeftKiDenseToResBlockedParallelSolvers( sparse_res, sym_res, *BliT, *BiT_buffer, res );
}

void sLinsysLeaf::addTermToSchurComplBlocked( sData *prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC )
{
   assert( prob == data );

   const bool sc_is_sym = true;

   if( use_local_RAC )
   {
      assert( prob->hasRAC() );
      BorderBiBlock border_right( prob->getLocalCrossHessian(), prob->getLocalA(), prob->getLocalC(),
            prob->getLocalF().getTranspose(), prob->getLocalG().getTranspose() );
      BorderBiBlock border_left_transp( prob->getLocalCrossHessian().getTranspose(), prob->getLocalA().getTranspose(), prob->getLocalC().getTranspose(),
            prob->getLocalF(), prob->getLocalG() );

      addBiTLeftKiBiRightToResBlockedParallelSolvers( sparseSC, sc_is_sym, border_left_transp, border_right, SC );
   }
   else
   {
      assert( !prob->hasRAC() );

      BorderBiBlock border_right( prob->getLocalF().getTranspose(), prob->getLocalG().getTranspose(), use_local_RAC );
      BorderBiBlock border_left_transp( prob->getLocalF(), prob->getLocalG(), use_local_RAC );

      addBiTLeftKiBiRightToResBlockedParallelSolvers( sparseSC, sc_is_sym, border_left_transp, border_right, SC );
   }
}

void sLinsysLeaf::mySymAtPutSubmatrix(SymMatrix& kkt_, 
					     GenMatrix& B_, GenMatrix&,
					     int locnx, int locmy, int)
{
  SparseSymMatrix& kkt = dynamic_cast<SparseSymMatrix&>(kkt_);
  SparseGenMatrix& B   = dynamic_cast<SparseGenMatrix&>(B_);
  //SparseGenMatrix& D   = reinterpret_cast<SparseGenMatrix&>(D_);

  int* jcolK = kkt.jcolM(); int* jcolB = B.jcolM(); //int* jcolD = D.jcolM(); 
  int* krowK = kkt.krowM(); int* krowB = B.krowM(); //int* krowD =  D.krowM();
  double* MK = kkt.M();     double* MB = B.M();

  for(int i=0; i<locmy; i++) {
    int itK = krowK[i+locnx];
    int j = krowB[i];

    for(; j<krowB[i+1]; j++) { 

      if(jcolB[j]<i+locnx) {
	jcolK[itK]=jcolB[j]; 
	MK[itK]=MB[j]; 
	itK++;
      }
    }
    jcolK[itK]=i+locnx; MK[itK] = 0.0; itK++;

    assert(j==krowB[i+1]);

    krowK[i+locnx+1]=itK;
  }
}

/* compute result += B_inner^T K^-1 Br */
void sLinsysLeaf::addInnerBorderKiInvBrToRes( DenseGenMatrix& result, BorderLinsys& Br )
{
   assert( Br.A.children.size() == 0 );

   const bool result_sparse = false;
   const bool result_sym = false;

   std::unique_ptr<BorderBiBlock> border_inner_tp{};
   std::unique_ptr<BorderBiBlock> border_right{};

   if( data->hasRAC() )
      border_inner_tp.reset( new BorderBiBlock( data->getLocalCrossHessian().getTranspose(),
         data->getLocalA().getTranspose(), data->getLocalC().getTranspose(), data->getLocalF(), data->getLocalG() ) );
   else
      border_inner_tp.reset( new BorderBiBlock( data->getLocalF(), data->getLocalG(), false ) );

   if( Br.has_RAC )
      border_right.reset(
            new BorderBiBlock(dynamic_cast<SparseGenMatrix&>(*Br.R.mat),
                  dynamic_cast<SparseGenMatrix&>(*Br.A.mat),
                  dynamic_cast<SparseGenMatrix&>(*Br.C.mat),
                  dynamic_cast<SparseGenMatrix&>(*Br.F.mat).getTranspose(),
                  dynamic_cast<SparseGenMatrix&>(*Br.G.mat).getTranspose()));
   else if( Br.use_local_RAC )
      border_right.reset(
            new BorderBiBlock(data->getLocalCrossHessian(), data->getLocalA(),
                  data->getLocalC(),
                  dynamic_cast<SparseGenMatrix&>(*Br.F.mat).getTranspose(),
                  dynamic_cast<SparseGenMatrix&>(*Br.G.mat).getTranspose()));
   else
      border_right.reset(
            new BorderBiBlock(
                  dynamic_cast<SparseGenMatrix&>(*Br.F.mat).getTranspose(),
                  dynamic_cast<SparseGenMatrix&>(*Br.G.mat).getTranspose(),
                  false));

   addBiTLeftKiBiRightToResBlockedParallelSolvers( result_sparse, result_sym, *border_inner_tp, *border_right, result);
}

/* compute result += B_inner^T K^-1 ( Br - Br_mod_border ) */
// TODO merge with other super similar method further up
void sLinsysLeaf::addInnerBorderKiInvBrToResDense( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border )
{
   int mRes, nRes;
   result.getSize(mRes, nRes);
#ifndef NDEBUG
   int nx_border, myl_border, mzl_border, dummy;
   if( Br.has_RAC )
      Br.R.getSize(dummy, nx_border);
   else if( Br.use_local_RAC )
      data->getLocalCrossHessian().getSize(dummy, nx_border);
   else
      nx_border = 0;
   Br.F.mat->getSize(myl_border, dummy);
   Br.G.mat->getSize(mzl_border, dummy);
   assert( nx_border + myl_border + mzl_border <= mRes );
#endif
   /* buffer for (Bri - (sum_j Brmodj * Xmodj)_i - Bi_{inner} X0)^T = Bri^T - X0^T Bi_{inner}^T - (sum_j Xmodj^T Brmodj^T)_i */
   // TODO : reuse and make member ? possible?

   std::unique_ptr<DenseGenMatrix> BiT_buffer( new DenseGenMatrix( mRes, dynamic_cast<SparseSymMatrix&>(*kkt).size() ) );

   /* Bi buffer and X0 are in transposed form for memory alignment reasons when solving with K_i */
   int m, n; BiT_buffer->getSize(m, n);
   BiT_buffer->atPutZeros(0, 0, m, n );

   /* put (Bri)^T into buffer
    *
    *                  nxb myb mzb
    *                [ RiT AiT CiT ]
    * Bri^T        = [  Fi  0   0  ]
    *                [  Gi  0   0  ]
    */
   std::unique_ptr<BorderBiBlock> BriT{};

   if( Br.has_RAC )
      BriT.reset( new BorderBiBlock( dynamic_cast<SparseGenMatrix&>(*Br.R.mat).getTranspose(),
            dynamic_cast<SparseGenMatrix&>(*Br.A.mat).getTranspose(), dynamic_cast<SparseGenMatrix&>(*Br.C.mat).getTranspose(),
            dynamic_cast<SparseGenMatrix&>(*Br.F.mat), dynamic_cast<SparseGenMatrix&>(*Br.G.mat) ) );
   else if( Br.use_local_RAC )
      BriT.reset( new BorderBiBlock( data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(), data->getLocalC().getTranspose(),
            dynamic_cast<SparseGenMatrix&>(*Br.F.mat), dynamic_cast<SparseGenMatrix&>(*Br.G.mat) ) );
   else
      BriT.reset( new BorderBiBlock( dynamic_cast<SparseGenMatrix&>(*Br.F.mat), dynamic_cast<SparseGenMatrix&>(*Br.G.mat), false ) );

   addBiTBorder( *BiT_buffer, *BriT);

   /* BiT_buffer += X_j^T Bmodj for all j */
   multRightDenseBorderModBlocked( Br_mod_border, *BiT_buffer );

   /* compute B_{inner}^T Ki^-1 Bi_buffer = B_{inner}^T Ki^-1 (Bri^T - sumj Xj^T Bmodj^T) */
   BorderBiBlock BiT_inner = data->hasRAC() ? BorderBiBlock( data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(), data->getLocalC().getTranspose(),
         data->getLocalF(), data->getLocalG() ) : BorderBiBlock( data->getLocalF(), data->getLocalG(), false );

   addBiTLeftKiDenseToResBlockedParallelSolvers( false, false, BiT_inner, *BiT_buffer, result );
}

/* compute result += B_inner^T K^-1 ( Br - Br_mod_border ) */
void sLinsysLeaf::addInnerBorderKiInvBrToRes( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool )
{
   if( Br_mod_border.empty() )
      addInnerBorderKiInvBrToRes( result, Br );
   else
      addInnerBorderKiInvBrToResDense( result, Br, Br_mod_border );
}

void sLinsysLeaf::addBorderTimesRhsToB0( StochVector& rhs, SimpleVector& b0, BorderLinsys& border )
{
   assert( border.F.children.size() == 0 );
   assert( rhs.children.size() == 0 );

   assert( rhs.vec );
   if( border.has_RAC )
   {
      assert( border.R.mat );
      assert( border.A.mat );
      assert( border.C.mat );
   }
   assert( border.F.mat );
   assert( border.G.mat );

   std::unique_ptr<BorderBiBlock> border_block{};

   if( border.has_RAC )
      border_block.reset(
            new BorderBiBlock(dynamic_cast<SparseGenMatrix&>(*border.R.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.A.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.C.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.F.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.G.mat) ) );
   else if( border.use_local_RAC )
      border_block.reset(
            new BorderBiBlock(data->getLocalCrossHessian(), data->getLocalA(),
                  data->getLocalC(),
                  dynamic_cast<SparseGenMatrix&>(*border.F.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.G.mat)));
   else
      border_block.reset(
            new BorderBiBlock(
                  dynamic_cast<SparseGenMatrix&>(*border.F.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.G.mat),
                  false));

   addBorderTimesRhsToB0( dynamic_cast<SimpleVector&>(*rhs.vec), b0, *border_block);
}

void sLinsysLeaf::addBorderTimesRhsToB0( SimpleVector& rhs, SimpleVector& b0, BorderBiBlock& border )
{
   int mFi, nFi; border.F.getSize(mFi, nFi);
   int mGi, nGi; border.G.getSize(mGi, nGi);

   int mRi{0}; int nRi{0};
   int mAi{0}; int nAi{0};
   int mCi{0}; int nCi{0};
   if( border.has_RAC )
   {
      border.R.getSize(mRi, nRi);
      border.A.getSize(mAi, nAi);
      border.C.getSize(mCi, nCi);

      assert( nFi == mRi );
      assert( rhs.length() == mRi + mAi + mCi );
   }
   else
      assert( rhs.length() >= nFi );

   assert( b0.length() >= nRi + mFi + mGi );
   const int nb0 = b0.length();

   SimpleVector& zi = dynamic_cast<SimpleVector&>(rhs);
   SimpleVector zi1 (&zi[0], nFi);

   if( border.has_RAC )
   {
      SimpleVector zi2 (&zi[mRi], mAi );
      SimpleVector zi3 (&zi[mRi + mAi], mCi);

      SimpleVector b1( &b0[0], nRi );

      border.R.transMult(1.0, b1, -1.0, zi1);
      border.A.transMult(1.0, b1, -1.0, zi2);
      border.C.transMult(1.0, b1, -1.0, zi3);
   }

   SimpleVector b2( &b0[nb0 - mFi - mGi], mFi );
   SimpleVector b3( &b0[nb0 - mGi], mGi );

   border.F.mult(1.0, b2, -1.0, zi1);
   border.G.mult(1.0, b3, -1.0, zi1);
}

void sLinsysLeaf::addBorderX0ToRhs( StochVector& rhs, const SimpleVector& x0, BorderLinsys& border )
{
   assert( border.F.children.size() == 0 );
   assert( rhs.children.size() == 0 );

   assert( rhs.vec );
   if( border.has_RAC )
   {
      assert( border.R.mat );
      assert( border.A.mat );
      assert( border.C.mat );
   }
   assert( border.F.mat );
   assert( border.G.mat );

   std::unique_ptr<BorderBiBlock> border_block{};

   if( border.has_RAC )
      border_block.reset(
            new BorderBiBlock(dynamic_cast<SparseGenMatrix&>(*border.R.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.A.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.C.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.F.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.G.mat) ) );
   else if( border.use_local_RAC )
      border_block.reset(
            new BorderBiBlock(data->getLocalCrossHessian(), data->getLocalA(),
                  data->getLocalC(),
                  dynamic_cast<SparseGenMatrix&>(*border.F.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.G.mat)));
   else
      border_block.reset(
            new BorderBiBlock(
                  dynamic_cast<SparseGenMatrix&>(*border.F.mat),
                  dynamic_cast<SparseGenMatrix&>(*border.G.mat),
                  false));

   addBorderX0ToRhs( dynamic_cast<SimpleVector&>(*rhs.vec), x0, *border_block);
}

void sLinsysLeaf::addBorderX0ToRhs( SimpleVector& rhs, const SimpleVector& x0, BorderBiBlock& border )
{
   int mFi, nFi; border.F.getSize(mFi, nFi);
   int mGi, nGi; border.G.getSize(mGi, nGi);
   int mRi{0}; int nRi{0};
   int mAi{0}; int nAi{0};
   int mCi{0}; int nCi{0};

   if( border.has_RAC )
   {
      border.R.getSize(mRi, nRi);
      border.A.getSize(mAi, nAi);
      border.C.getSize(mCi, nCi);

      assert( nFi == mRi );
      assert( rhs.length() == mRi + mAi + mCi );
   }
   else
      assert( rhs.length() >= nFi );

   assert( x0.length() >= nRi + mFi + mGi );
   const int nb0 = x0.length();

   double* rhsi1 = &rhs[0];

   if( border.has_RAC )
   {
      const double* x1 = &x0[0];

      double* rhsi2 = &rhs[mRi];
      double* rhsi3 = &rhs[mRi + mAi];
      border.R.mult( 1.0, rhsi1, 1, -1.0, x1, 1 );
      border.A.mult( 1.0, rhsi2, 1, -1.0, x1, 1 );
      border.C.mult( 1.0, rhsi3, 1, -1.0, x1, 1 );
   }

   const double* x2 = &x0[nb0 - mFi - mGi];
   const double* x3 = &x0[nb0 - mGi];

   border.F.transMult( 1.0, rhsi1, 1, -1.0, x2, 1 );
   border.G.transMult( 1.0, rhsi1, 1, -1.0, x3, 1 );
}
