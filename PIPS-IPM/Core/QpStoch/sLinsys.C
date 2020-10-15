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

void sLinsys::addLniZiHierarchyBorder( DenseGenMatrix& result, StringGenMatrix& R_border, StringGenMatrix& A_border,
      StringGenMatrix& C_border, StringGenMatrix& F_border, StringGenMatrix& G_border)
{
   assert( R_border.children.size() == 0 );
   assert( A_border.mat );
   assert( C_border.mat );
   assert( F_border.mat );
   assert( G_border.mat );
   assert( R_border.mat );

   SparseGenMatrix& A = *A_border.mat;
   SparseGenMatrix& C = *C_border.mat;
   SparseGenMatrix& F = *F_border.mat;
   SparseGenMatrix& G = *G_border.mat;
   SparseGenMatrix& R = *R_border.mat;

   addTermToSchurComplBlocked(data, false, false, R, A, C, F, G, result);
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
 *        nx_  l_  myl_ mzl_
 * nxl_ [ Ri   0   FiT  GiT ]
 * my_  [ Ai   0    0    0  ]
 * mz_  [ Ci   0    0    0  ]
 *
 * where the size of the zero part depends on the size of the factor K^-1:
 * l = n_K^-1 - nx - myl - mzl
 *
 * computes (B^T K^{-1} B_right)^T and adds it to the SC
 */
void sLinsys::addTermToSchurComplBlocked(/*const*/sData* prob, bool sparseSC, bool symSC,
      SparseGenMatrix& R_right, SparseGenMatrix& A_right, SparseGenMatrix& C_right,
      SparseGenMatrix& F_right, SparseGenMatrix& G_right, DoubleMatrix& result)
{
   int tmp;

   int m_res, n_res;
   result.getSize(m_res, n_res);
   assert( m_res >= 0 && n_res >= 0 );

   int nxl_right, nx_right;
   R_right.getSize(nxl_right, nx_right);
   assert( nx_right >= 0 && nxl_right >= 0 );

   int my_right;
   A_right.getSize(my_right, tmp);
   assert( nx_right == tmp );
   assert( my_right >= 0 );

   int mz_right;
   C_right.getSize(mz_right, tmp);
   assert( nx_right == tmp );
   assert( mz_right >= 0 );

   int myl_right, mzl_right;

   F_right.getSize(myl_right, tmp);
   assert( tmp == nxl_right );
   assert( myl_right >= 0 );

   G_right.getSize(mzl_right, tmp);
   assert( tmp == nxl_right );
   assert( mzl_right >= 0 );

   assert( nx_right + myl_right + mzl_right <= n_res);
   assert( nxl_right == locnx && my_right == locmy && mz_right == locmz );

#ifndef HIERARCHICAL
   assert( myl_right == locmyl && mzl_right == locmzl );
#endif

   SimpleVectorBase<int> nnzPerColRAC(nx_right);

   R_right.addNnzPerCol(nnzPerColRAC);
   A_right.addNnzPerCol(nnzPerColRAC);
   C_right.addNnzPerCol(nnzPerColRAC);

   const int withF = (locmyl > 0);
   const int withG = (locmzl > 0);
   const int N = nxl_right + my_right + mz_right;

   assert(nThreads >= 1);

   // storage for sets of dense cols
   if( colsBlockDense == nullptr )
      colsBlockDense = new double[blocksizemax * N];
   else
   {
      delete[] colsBlockDense;
      colsBlockDense = new double[blocksizemax * N];
   }
   // to save original column index of each column in colsBlockTrans
   if( colId == nullptr )
      colId = new int[blocksizemax];
   // indicating whether a right hand side is zero
#if 0
   if( colSparsity == nullptr )
      colSparsity = new int[N * blocksizemax];
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
   while( colpos < nx_right )
   {
      int blocksize = 0;

      for( ; colpos < nx_right && blocksize < blocksizemax; colpos++ )
         if( nnzPerColRAC[colpos] != 0 )
            colId[blocksize++] = colpos;

      if( blocksize == 0 )
         break;

      memset(colsBlockDense, 0, blocksize * N * sizeof(double));

      if( colSparsity )
         memset(colSparsity, 0, N * sizeof(int));

      R_right.fromGetColsBlock(colId, blocksize, N, 0, colsBlockDense, colSparsity);
      A_right.fromGetColsBlock(colId, blocksize, N, nxl_right, colsBlockDense, colSparsity);
      C_right.fromGetColsBlock(colId, blocksize, N, (nxl_right + my_right), colsBlockDense, colSparsity);

      solver->solve(blocksize, colsBlockDense, colSparsity);

      multLeftSchurComplBlocked( prob, colsBlockDense, colId, blocksize, sparseSC, symSC, result);
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

      SimpleVectorBase<int> nnzPerColFt(myl_right);
      F_right.addNnzPerRow(nnzPerColFt);

      const int nxMySC = m_res - myl_right - mzl_right;

      colpos = 0;
      // do block-wise multiplication for columns of F^T part
      while( colpos < myl_right )
      {
         int blocksize = 0;

         for( ; colpos < myl_right && blocksize < blocksizemax; colpos++ )
            if( nnzPerColFt[colpos] != 0 )
               colId[blocksize++] = colpos;

         if( blocksize == 0 )
            break;

         if( colSparsity )
            memset(colSparsity, 0, N * sizeof(int));

         memset(colsBlockDense, 0, blocksize * N * sizeof(double));

         // get column block from Ft (i.e., row block from F)
         F_right.fromGetRowsBlock(colId, blocksize, N, 0, colsBlockDense, colSparsity);

         solver->solve(blocksize, colsBlockDense, colSparsity);

         for( int i = 0; i < blocksize; i++ )
            colId[i] += nxMySC;

         multLeftSchurComplBlocked( prob, colsBlockDense, colId, blocksize, sparseSC, symSC, result);
      }
   }

   // do we have linking inequality constraints?
   if( withG )
   {
      //                       (G^T)
      //     SC +=  B^T  K^-1  (0  )
      //                       (0  )

      SimpleVectorBase<int> nnzPerColGt(mzl_right);
      G_right.addNnzPerRow(nnzPerColGt);
      const int nxMyMzSC = m_res - mzl_right;

      colpos = 0;

      // do block-wise multiplication for columns of G^T part
      while( colpos < mzl_right )
      {
         int blocksize = 0;

         for( ; colpos < mzl_right && blocksize < blocksizemax; colpos++ )
            if( nnzPerColGt[colpos] != 0 )
               colId[blocksize++] = colpos;

         if( blocksize == 0 )
            break;

         if( colSparsity )
            memset(colSparsity, 0, N * sizeof(int));

         memset(colsBlockDense, 0, blocksize * N * sizeof(double));

         G_right.fromGetRowsBlock(colId, blocksize, N, 0, colsBlockDense, colSparsity);

         solver->solve(blocksize, colsBlockDense, colSparsity);

         for( int i = 0; i < blocksize; i++ )
             colId[i] += nxMyMzSC;

          multLeftSchurComplBlocked( prob, colsBlockDense, colId, blocksize, sparseSC, symSC, result);
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

void sLinsys::addTermToSchurComplBlockedParallelSolvers(sData *prob, bool sparseSC, bool symSC, SymMatrix& SC)
{
   // TODO make member?
   const int blocksize = 20;

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

   assert(N == locmy);
   assert(locmyl >= 0);
   assert(locmzl >= 0);

   const int NP = SC.size();
   assert(NP >= nxP);

   const int nxMyP = NP - locmyl - locmzl;
   const int nxMyMzP = NP - locmzl;

   assert( nxMyP >= 0 );
   assert( nxMyMzP >= 0 );

   if( nxP == -1 )
      C.getSize(N, nxP);

   int N2, nxP2;
   C.getSize(N2, nxP2);
   const bool withC = (nxP2 != -1);

   // TODO : is this correct?
   if( nxP == -1 )
      nxP = NP;

   N = locnx + locmy + locmz;

   SimpleVectorBase<int> nnzPerColRAC(nxP);

   if( withR )
      R.addNnzPerCol(nnzPerColRAC);

   if( withA )
      A.addNnzPerCol(nnzPerColRAC);

   if( withC )
      C.addNnzPerCol(nnzPerColRAC);

   const int withF = (locmyl > 0);
   const int withG = (locmzl > 0);

   assert(nThreads >= 1);

   // indicating whether a right hand side is zero
   if( colSparsity == nullptr )
   {
      colSparsity = new int[N * blocksize * n_solvers];
      memset(colSparsity, 0, N * blocksize * n_solvers * sizeof(int));
   }

   if( colsBlockDense == nullptr )
      colsBlockDense = new double[N * blocksize * n_solvers];

   // to save original column index of each column in colsBlockTrans
   if( colId == nullptr )
      colId = new int[N * blocksizemax];

#ifdef TIME_SCHUR
   const double t_start = omp_get_wtime();
#endif


   R.updateTransposed();
   A.updateTransposed();
   C.updateTransposed();

   const int chunks_RAC = std::ceil( static_cast<double>(nxP) / blocksize );
   //                       (R)
   //     SC +=  B^T  K^-1  (A)
   //                       (C)
   #pragma omp parallel for schedule(dynamic, 1) num_threads(n_solvers)
   for( int i = 0; i < chunks_RAC; i++ )
   {
      const int actual_blocksize = std::min( (i + 1) * blocksize, nxP) - i * blocksize;

      int nrhs = 0;
      const int id = omp_get_thread_num();

      int * colId_loc = colId + id * blocksize;

      for( int j = 0; j < actual_blocksize; ++j )
         if( nnzPerColRAC[j + i * blocksize] != 0 )
            colId_loc[nrhs++] = j + i * blocksize;

      if( nrhs == 0 )
         continue;

      double* colsBlockDense_loc = colsBlockDense + id * N * blocksize;
      memset(colsBlockDense_loc, 0, blocksize * N * sizeof(double));

      int* colSparsity_loc = colSparsity + id * N * blocksize;

      R.fromGetColsBlock(colId_loc, nrhs, N, 0, colsBlockDense_loc, colSparsity_loc);
      A.fromGetColsBlock(colId_loc, nrhs, N, locnx, colsBlockDense_loc, colSparsity_loc);
      C.fromGetColsBlock(colId_loc, nrhs, N, (locnx + locmy), colsBlockDense_loc, colSparsity_loc);

      solvers_blocked[id]->solve(nrhs, colsBlockDense_loc, colSparsity_loc);

      // TODO: thread-safe?
//      #pragma omp critical
      {
         multLeftSchurComplBlocked(prob, colsBlockDense_loc, colId_loc, nrhs, sparseSC, symSC, SC);
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

      SimpleVectorBase<int> nnzPerColFt(locmyl);
      F.addNnzPerRow(nnzPerColFt);

      const int chunks_F = std::ceil( static_cast<double>(locmyl) / blocksize );

      // do block-wise multiplication for columns of F^T part
      #pragma omp parallel for schedule(dynamic, 1) num_threads(n_solvers)
      for(int i = 0; i < chunks_F; ++i )
      {
         const int actual_blocksize = std::min( (i + 1) * blocksize, locmyl) - i * blocksize;

         int nrhs = 0;
         const int id = omp_get_thread_num();

         int * colId_loc = colId + id * blocksize;

         for( int j = 0; j < actual_blocksize; ++j )
            if( nnzPerColFt[j + i * blocksize] != 0 )
               colId_loc[nrhs++] = j + i * blocksize;

         if( nrhs == 0 )
            continue;

         double* colsBlockDense_loc = colsBlockDense + id * N * blocksize;
         memset(colsBlockDense_loc, 0, blocksize * N * sizeof(double));

         int* colSparsity_loc = colSparsity + id * N * blocksize;

         // get column block from Ft (i.e., row block from F)
         F.fromGetRowsBlock(colId_loc, nrhs, N, 0, colsBlockDense_loc, colSparsity_loc);

         solvers_blocked[id]->solve(nrhs, colsBlockDense_loc, colSparsity_loc);

         for( int j = 0; j < nrhs; ++j )
            colId_loc[j] += nxMyP;

//         #pragma omp critical
         {
            // TODO : thread-safe?
            multLeftSchurComplBlocked(prob, colsBlockDense_loc, colId_loc, nrhs, sparseSC, symSC, SC);
         }
      }
   }

   // do we have linking inequality constraints?
   if( withG )
   {
      //                       (G^T)
      //     SC +=  B^T  K^-1  (0  )
      //                       (0  )

      SimpleVectorBase<int> nnzPerColGt(locmzl);
      G.addNnzPerRow(nnzPerColGt);

      const int chunks_G = std::ceil( static_cast<double>(locmzl) / blocksize );

      // do block-wise multiplication for columns of G^T part
      #pragma omp parallel for schedule(dynamic, 1) num_threads(n_solvers)
      for( int i = 0; i < chunks_G; ++i )
      {
         const int actual_blocksize = std::min( (i + 1) * blocksize, locmzl) - i * blocksize;

         int nrhs = 0;
         const int id = omp_get_thread_num();

         int * colId_loc = colId + id * blocksize;

         for( int j = 0; j < actual_blocksize; ++j )
            if( nnzPerColGt[j + i * blocksize] != 0 )
               colId_loc[nrhs++] = j + i * blocksize;

         if( nrhs == 0 )
            continue;

         double* colsBlockDense_loc = colsBlockDense + id * N * blocksize;
         memset(colsBlockDense_loc, 0, blocksize * N * sizeof(double));

         int* colSparsity_loc = colSparsity + id * N * blocksize;

         G.fromGetRowsBlock(colId_loc, nrhs, N, 0, colsBlockDense_loc, colSparsity_loc);

         solvers_blocked[id]->solve(nrhs, colsBlockDense_loc, colSparsity_loc);

         for( int j = 0; j < nrhs; ++j )
            colId_loc[j] += nxMyMzP;

//         #pragma omp critical
         {
            // TODO : thread-safe?
            multLeftSchurComplBlocked(prob, colsBlockDense_loc, colId_loc, nrhs, sparseSC, symSC, SC);
         }
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


#include <set>
#include <algorithm>

// we load the calculated columns into rows of out
// to match the column-major scalapack format

void sLinsys::addColsToDenseSchurCompl(sData *prob, 
				       DenseGenMatrix& out, 
				       int startcol, int endcol) 
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();

  int ncols = endcol-startcol;
  int N, nxP, ncols_t, N_out;
  A.getSize(N, nxP); assert(N==locmy);
  out.getSize(ncols_t, N_out); 
  assert(N_out == nxP);
  assert(endcol <= nxP &&  ncols_t >= ncols);

  if(nxP==-1) C.getSize(N,nxP);
  //if(nxP==-1) nxP = NP;

  N = locnx+locmy+locmz;
  DenseGenMatrix cols(ncols,N);
  bool allzero = true;
  memset(cols[0],0,N*ncols*sizeof(double));

  R.getStorageRef().fromGetColBlock(startcol, &cols[0][0],
				    N, endcol-startcol, allzero);
  A.getStorageRef().fromGetColBlock(startcol, &cols[0][locnx], 
				    N, endcol-startcol, allzero);
  C.getStorageRef().fromGetColBlock(startcol, &cols[0][locnx+locmy], 
				    N, endcol-startcol, allzero);

  //int mype; MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  //printf("solving with multiple RHS %d \n", mype);	
  solver->solve(cols);
  //printf("done solving %d \n", mype);
  
  
  const int blocksize = 20;
  
  for (int it=0; it < ncols; it += blocksize) {
    int end = std::min(it+blocksize,ncols);
    int numcols = end-it;
    assert(false); //add Rt*x -- and test the code
    // SC-=At*y
    A.getStorageRef().transMultMat( 1.0, out[it], numcols, N_out,  
				  -1.0, &cols[it][locnx], N);
    // SC-=Ct*z
    C.getStorageRef().transMultMat( 1.0, out[it], numcols, N_out,
				  -1.0, &cols[it][locnx+locmy], N);
  }
}

/* solve own linear system with border data
 *
 * rhs-block^T looks like
 *  [ R1 F1T G1T ]^T       [ RN FNT GNT ]^T [  0   A0 (C0) F0V G0V ]
 *  [ A1  0   0  ]    ...  [ AN  0   0  ]   [ F0C  0  (0 )  0   0  ]
 *  [ C1  0   0  ]         [ CN  0   0  ]   [ G0C  0  (0 )  0   0  ]
 *
 */
void sLinsys::addInnerToHierarchicalSchurComplement( DenseSymMatrix& schur_comp, sData* data_border )
{
   /* only called in sLinsysRootAug */
   assert( !is_hierarchy_root );

   /* get right hand side parts */
   StringGenMatrix& R_border = *dynamic_cast<BorderedSymMatrix&>(*data_border->Q).border_vertical;
   StringGenMatrix& A_border = *dynamic_cast<BorderedGenMatrix&>(*data_border->A).border_left;
   StringGenMatrix& C_border = *dynamic_cast<BorderedGenMatrix&>(*data_border->C).border_left;

   StringGenMatrix& F_border = *dynamic_cast<BorderedGenMatrix&>(*data_border->A).border_bottom;
   StringGenMatrix& G_border = *dynamic_cast<BorderedGenMatrix&>(*data_border->C).border_bottom;

   /* compute Schur Complement right hand sides SUM_i Bi_{inner} K^-1 Bi_{border}
    * (keep in mind that in Bi_{inner} and the SC we projected C0 Omega0 out)
    */
   int nx_border, myl_border, mzl_border, dummy;
   A_border.getSize(dummy, nx_border);
   F_border.getSize(myl_border, dummy);
   G_border.getSize(mzl_border, dummy);

   const int m_buffer = nx_border + myl_border + mzl_border;
   const int n_buffer = locnx + locmy + locmyl + locmzl;

   // buffer for B0_{outer} - SUM_i Bi_{inner}^T Ki^{-1} Bi_{outer}, stored in transposed form (for quick access of cols in solve)
   DenseGenMatrix* buffer = new DenseGenMatrix(m_buffer, n_buffer);
   LsolveHierarchyBorder(*buffer, R_border, A_border, C_border, F_border, G_border);

   SparseGenMatrix& A0_border = *dynamic_cast<BorderedGenMatrix&>(*data_border->A).border_left->mat;
   SparseGenMatrix& F0vec_border = *dynamic_cast<BorderedGenMatrix&>(*data_border->A).border_left->mat_link;
   SparseGenMatrix& F0cons_border = *dynamic_cast<BorderedGenMatrix&>(*data_border->A).border_bottom->mat;

   SparseGenMatrix& G0vec_border = *dynamic_cast<BorderedGenMatrix&>(*data_border->C).border_left->mat_link;
   SparseGenMatrix& G0cons_border = *dynamic_cast<BorderedGenMatrix&>(*data_border->C).border_bottom->mat;

   finalizeZ0Hierarchical(*buffer, A0_border, F0vec_border, F0cons_border, G0vec_border, G0cons_border);

   // TODO : solve with Schur Complement



   // TODO : check that in hierarchical mode Schur Complement gets allreduced

   // TODO : each process only solves a part of the hat{B}0 rhs with the Schur Complement -> then allreduce

   // TODO : for each child solve Xi = Ki^-1 (Bi_{outer} - Bi_{inner} X0)

   // TODO : for each child multiply and add Bi_{outer}^T X_i to Schur Complement

   assert( false && "TODO : implement" );
}

// adds only lower triangular elements to out

void sLinsys::symAddColsToDenseSchurCompl(sData *prob,
				       double *out,
				       int startcol, int endcol)
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  int N, nxP;
  A.getSize(N, nxP); assert(N==locmy);
  //out.getSize(ncols, N); assert(N == nxP);
  assert(endcol <= nxP);

  if( nxP == -1 ) C.getSize(N,nxP);

  N = locnx+locmy+locmz;

  const int BLOCKSIZE = 40;
  
  DenseGenMatrix cols(BLOCKSIZE,N);
  int outi = 0;
  

  for (int col = startcol; col < endcol; col += BLOCKSIZE) {
    int ecol = std::min(col+BLOCKSIZE,endcol);
    int nbcols = ecol-col;
    
    
    memset(cols[0],0,BLOCKSIZE*N*sizeof(double));
    
 
    bool allzero = true;
    assert(false); //! code needs to consider the cross Hessian; fixme
    A.getStorageRef().fromGetColBlock(col, &cols[0][locnx], N, nbcols, allzero);
    C.getStorageRef().fromGetColBlock(col, &cols[0][locnx+locmy], N, nbcols, allzero);
    
    if (!allzero) {
      solver->solve(cols);
      
      A.getStorageRef().transMultMatLower(out+outi, nbcols, col,
					-1.0, &cols[0][locnx], N);
      C.getStorageRef().transMultMatLower(out+outi, nbcols, col,
					  -1.0, &cols[0][locnx+locmy], N); 
    }
    for (int c = col; c < ecol; c++) {
      outi += nxP-c;
    }
  }
}

void sLinsys:: multLeftSparseSchurComplBlocked( /*const*/ sData* prob, /*const*/double* colsBlockDense,
      const int* colId, int blocksize, SparseSymMatrix& SC)
{
   const int withF = (locmyl > 0);
   const int withG = (locmzl > 0);
   const int N = locnx + locmy + locmz;
   const int NP = SC.size();
   const int nxMyP = NP - locmyl - locmzl;
   const int nxMyMzP = NP - locmzl;

   SparseGenMatrix& A = prob->getLocalA();
   SparseGenMatrix& C = prob->getLocalC();
   SparseGenMatrix& F = prob->getLocalF();
   SparseGenMatrix& G = prob->getLocalG();
   SparseGenMatrix& R = prob->getLocalCrossHessian();

   // multiply each column with left factor of SC
   // todo: #pragma omp parallel for schedule(dynamic, 10)
   for( int it_col = 0; it_col < blocksize; it_col++ )
   {
      double* const col = &colsBlockDense[it_col * N];
      const int row_sc = colId[it_col];

      // SC+=R^T*x
      R.transmultMatSymUpper(1.0, SC, -1.0, &col[0], row_sc, 0);

      // SC+=A^T*y
      A.transmultMatSymUpper(1.0, SC, -1.0, &col[locnx], row_sc, 0);

      // SC+=C^T*z
      C.transmultMatSymUpper(1.0, SC, -1.0, &col[locnx + locmy], row_sc, 0);

      // do we have linking equality constraints? If so, set SC+=F*x
      if( withF )
         F.multMatSymUpper(1.0, SC, -1.0, &col[0], row_sc, nxMyP);

      // do we have linking inequality constraints? If so, set SC+=G*x
      if( withG )
         G.multMatSymUpper(1.0, SC, -1.0, &col[0], row_sc, nxMyMzP);
   }
}

void sLinsys:: multLeftDenseSchurComplBlocked( /*const*/sData* prob, /*const*/double* colsBlockDense,
      const int* colId, int blocksize, int ncolsSC, double** SC)
{
   const int withF = (locmyl > 0);
   const int withG = (locmzl > 0);
   const int nxMySC = ncolsSC - locmyl - locmzl;
   const int nxMyMzSC = ncolsSC - locmzl;
   const int blockcolLength = locnx + locmy + locmz;

   SparseGenMatrix& A = prob->getLocalA();
   SparseGenMatrix& C = prob->getLocalC();
   SparseGenMatrix& F = prob->getLocalF();
   SparseGenMatrix& G = prob->getLocalG();
   SparseGenMatrix& R = prob->getLocalCrossHessian();

   assert( ncolsSC >= 1 );
   assert(nxMySC >= 0 && nxMyMzSC >= 0);

   // multiply each column with left factor of SC todo add OMP
   for( int it_col = 0; it_col < blocksize; it_col++ )
   {
      double* const col = &colsBlockDense[it_col * blockcolLength];
      const int row_sc = colId[it_col];

      // SC+=R^T*x
      R.transMult(1.0, &SC[row_sc][0], 1, -1.0, &col[0], 1);

      // SC+=A^T*y
      A.transMult(1.0, &SC[row_sc][0], 1, -1.0, &col[locnx], 1);

      // SC+=C^T*z
      C.transMult(1.0, &SC[row_sc][0], 1, -1.0, &col[locnx + locmy], 1);

      // do we have linking equality constraints? If so, set SC+=F*x
      if( withF )
         F.mult(1.0, &SC[row_sc][nxMySC], 1, -1.0, &col[0], 1);

      // do we have linking inequality constraints? If so, set SC+=G*x
      if( withG )
         G.mult(1.0, &SC[row_sc][nxMyMzSC], 1, -1.0, &col[0], 1);
   }
}

void sLinsys::multLeftSchurComplBlocked( /*const*/sData* prob, /*const*/double* colsBlockDense,
      const int* colId, int blocksize, bool sparseSC, bool symSC, DoubleMatrix& SC)
{
   assert(colId && colsBlockDense);

   if( sparseSC )
   {
      assert( symSC );
      SparseSymMatrix& SC_sparse = dynamic_cast<SparseSymMatrix&>(SC);
      multLeftSparseSchurComplBlocked( prob, colsBlockDense, colId, blocksize, SC_sparse);
   }
   else
   {
      double** SCarray;
      int ncolsSC;
      if( symSC )
      {
         DenseSymMatrix& SC_dense = dynamic_cast<DenseSymMatrix&>(SC);
         SCarray = SC_dense.mStorage->M;
         ncolsSC = SC_dense.size();
      }
      else
      {
         DenseGenMatrix& SC_dense = dynamic_cast<DenseGenMatrix&>(SC);
         SCarray = SC_dense.mStorage->M;
         ncolsSC = SC_dense.mStorage->n;
      }
      assert(ncolsSC >= locmyl + locmzl);

      multLeftDenseSchurComplBlocked(prob, colsBlockDense, colId, blocksize, ncolsSC, SCarray);
   }

}

/*

void sLinsys::symAddColsToDenseSchurCompl(sData *prob, 
				       double *out, 
				       int startcol, int endcol) 
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();

  int N, nxP, NP, ncols;
  A.getSize(N, nxP); assert(N==locmy);
  N = nxP;
  //out.getSize(ncols, N); assert(N == nxP);
  assert(endcol <= nxP);

  if(nxP==-1) C.getSize(N,nxP);
  if(nxP==-1) nxP = NP;

  N = locnx+locmy+locmz;
  SimpleVector col(N);
  int outi = 0;

  for(int it=startcol; it<endcol; it++) {
    
    double* pcol = &col[0];
    for(int it1=0; it1<locnx; it1++) pcol[it1]=0.0;

    A.fromGetDense(0, it, &col[locnx],  1, locmy, 1);    
    C.fromGetDense(0, it, &col[locnx+locmy], 1, locmz, 1);

    solver->solve(col);

    //here we have colGi = inv(H_i)* it-th col of Gi^t
    
    //now do colSC = Gi * inv(H_i)* it-th col of Gi^t
    //for(int it1=0; it1<locnx; it1++) pcol[it1]=0.0;
    
    // SC-=At*y
    A.getStorageRef().transMultLower( 1.0, out+outi,
                        -1.0, &col[locnx], it);
    // SC-=Ct*z
    C.getStorageRef().transMultLower( 1.0, out+outi, 
                        -1.0, &col[locnx+locmy], it);

    int nelts = nxP-it;
    outi += nelts;
  }
}*/
