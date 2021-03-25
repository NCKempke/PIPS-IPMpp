/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sFactory.h"

#include "sData.h"
#include "sTreeCallbacks.h"
#include "StochInputTree.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"

#include "sVars.h"
#include "sResiduals.h"

#include "sLinsysRoot.h"
#include "sLinsysLeaf.h"

#include "DeSymIndefSolver.h"

#include "pipsport.h"
#include "StochOptions.h"
#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>

class PardisoProjectSolver;
class PardisoMKLSolver;
class MumpsSolverLeaf;
class PardisoSchurSolver;
class Ma27Solver;
class Ma57Solver;

#include "sLinsysLeafSchurSlv.h"
#include "sLinsysLeafMumps.h"

#ifdef WITH_MA57
#include "Ma57Solver.h"
#endif

#ifdef WITH_MA27
#include "Ma27Solver.h"
#endif

#ifdef WITH_PARDISO
#include "PardisoProjectSolver.h"
#include "../LinearSolvers/PardisoSolver/PardisoSchurSolver/PardisoProjectSchurSolver.h"
#endif

#ifdef WITH_MKL_PARDISO
#include "PardisoMKLSolver.h"
#include "../LinearSolvers/PardisoSolver/PardisoSchurSolver/PardisoMKLSchurSolver.h"
#endif

#ifdef WITH_MUMPS
#include "../LinearSolvers/MumpsSolver/MumpsSolverLeaf.h"
#endif


sFactory::sFactory( StochInputTree* inputTree, MPI_Comm comm)
  : tree( new sTreeCallbacks(inputTree) )
{
  tree->assignProcesses(comm);

  tree->computeGlobalSizes();
  //now the sizes of the problem are available, set them for the parent class
  tree->getGlobalSizes(nx, my, mz);
}

sFactory::~sFactory()
{
   if(tree)
      delete tree;
}

DoubleLinearSolver* sFactory::newLeafSolver( const DoubleMatrix* kkt_ )
{
   const SparseSymMatrix* kkt = dynamic_cast<const SparseSymMatrix*>(kkt_);
   assert( kkt );

   const SolverType leaf_solver = pips_options::getSolverLeaf();

   if( !pips_options::getBoolParameter( "SC_COMPUTE_BLOCKWISE" ) )
   {
      if( leaf_solver == SolverType::SOLVER_MUMPS )
      {
#ifdef WITH_MUMPS
         return new MumpsSolverLeaf(kkt);
#endif
      }
      else if( leaf_solver == SolverType::SOLVER_PARDISO )
      {
#ifdef WITH_PARDISO
         return new PardisoProjectSchurSolver(kkt);
#endif
      }
      else if( leaf_solver == SolverType::SOLVER_MKL_PARDISO )
      {
#ifdef WITH_MKL_PARDISO
         return new PardisoMKLSchurSolver(kkt);
#endif
      }

      PIPS_MPIabortIf(true, "No leaf solver for Schur Complement computation could be found - should not happen..");
   }
   else
   {
      if( leaf_solver == SolverType::SOLVER_PARDISO )
      {
#ifdef WITH_PARDISO
         return new PardisoProjectSolver(kkt);
#endif
      }
      else if( leaf_solver == SolverType::SOLVER_MKL_PARDISO )
      {
#ifdef WITH_MKL_PARDISO
         return new PardisoMKLSolver(kkt);
#endif
      }
      else if( leaf_solver == SolverType::SOLVER_MA57 )
      {
#ifdef WITH_MA57
         return new Ma57Solver(kkt);
#endif
      }
      else if( leaf_solver == SolverType::SOLVER_MA27 )
      {
#ifdef WITH_MA27
         return new Ma27Solver(kkt);
#endif
      }
      else if( leaf_solver == SolverType::SOLVER_MUMPS )
      {
#ifdef WITH_MUMPS
         return new MumpsSolverLeaf(kkt);
#endif
      }

      PIPS_MPIabortIf(true, "No leaf solver for Blockwise Schur Complement computation could be found - should not happen..");
   }
   return nullptr;
}

sLinsysLeaf* sFactory::newLinsysLeaf(sData* prob,
			OoqpVector* dd, OoqpVector* dq,
			OoqpVector* nomegaInv, OoqpVector* rhs)
{
   assert( prob );
   static bool printed = false;
   const SolverType leaf_solver = pips_options::getSolverLeaf();

   if( !pips_options::getBoolParameter( "SC_COMPUTE_BLOCKWISE" ) )
   {
      if( PIPS_MPIgetRank() == 0 && !printed )
          std::cout << "Using " << leaf_solver << " for leaf schur complement computation - sFactory\n";
      printed = true;

      if( leaf_solver == SolverType::SOLVER_MUMPS )
      {
#ifdef WITH_MUMPS
         return new sLinsysLeafMumps(this, prob, dd, dq, nomegaInv, rhs);
#endif
      }
      else if( leaf_solver == SolverType::SOLVER_PARDISO || leaf_solver == SolverType::SOLVER_MKL_PARDISO )
         return new sLinsysLeafSchurSlv(this, prob, dd, dq, nomegaInv, rhs);
      else
      {
         if( PIPS_MPIgetRank() == 0 )
         {
            std::cout << "WARNING: did not specify SC_COMPUTE_BLOCKWISE but " << leaf_solver << " which can only compute the Schur Complement blockwise\n";
            std::cout << "WARNING: checking for PARDISO or MUMPS instead...\n";
         }

         SolverType solver = SolverType::SOLVER_NONE;
         if( pips_options::isSolverAvailable( SolverType::SOLVER_PARDISO ) )
            solver = SolverType::SOLVER_PARDISO;
         else if( pips_options::isSolverAvailable( SolverType::SOLVER_MKL_PARDISO) )
            solver = SolverType::SOLVER_MKL_PARDISO;
         else if( pips_options::isSolverAvailable( SolverType::SOLVER_MUMPS) )
            solver = SolverType::SOLVER_MUMPS;

         if( solver != SolverType::SOLVER_NONE )
         {
            if( PIPS_MPIgetRank() == 0 )
               std::cout << " Found solver " << solver << " - using that for leaf computations\n";
            pips_options::setIntParameter( "LINEAR_LEAF_SOLVER", solver );
            return new sLinsysLeafSchurSlv(this, prob, dd, dq, nomegaInv, rhs);
         }

         PIPS_MPIabortIf(true, "Error: Could not find suitable solver - please specify SC_COMPUTE_BLOCKWISE");
         return nullptr;
      }
   }
   else
   {
      if( PIPS_MPIgetRank() == 0 && !printed )
         std::cout << "Using " << leaf_solver << " for blockwise Schur Complement computation - deactivating distributed preconditioner - sFactory\n";
      pips_options::setBoolParameter("PRECONDITION_DISTRIBUTED", false);
      printed = true;

      if( leaf_solver == SolverType::SOLVER_MUMPS )
      {
#ifdef WITH_MUMPS
         return new sLinsysLeafMumps(this, prob, dd, dq, nomegaInv, rhs);
#endif
      }
      else
         return new sLinsysLeaf(this, prob, dd, dq, nomegaInv, rhs);
   }
   return nullptr;
}


void dumpaug(int nx, SparseGenMatrix &A, SparseGenMatrix &C) {

    long long my, mz, nx_1, nx_2;
    A.getSize(my,nx_1);
    C.getSize(mz,nx_2);
    assert(nx_1 == nx_2);

    int nnzA = A.numberOfNonZeros();
    int nnzC = C.numberOfNonZeros();
    std::cout << "augdump  nx=" << nx << "\n";
    std::cout << "A: " << my << "x" << nx_1 << "   nnz=" << nnzA << "\n"
              << "C: " << mz << "x" << nx_1 << "   nnz=" << nnzC << "\n";

	std::vector<double> eltsA(nnzA), eltsC(nnzC), elts(nnzA+nnzC);
	std::vector<int> colptrA(nx_1+1),colptrC(nx_1+1), colptr(nx_1+1), rowidxA(nnzA), rowidxC(nnzC), rowidx(nnzA+nnzC);
	A.getStorageRef().transpose(&colptrA[0],&rowidxA[0],&eltsA[0]);
	C.getStorageRef().transpose(&colptrC[0],&rowidxC[0],&eltsC[0]);

	int nnz = 0;
	for (int col = 0; col < nx_1; col++) {
		colptr[col] = nnz;
		for (int r = colptrA[col]; r < colptrA[col+1]; r++) {
			int row = rowidxA[r]+nx+1; // +1 for fortran
			rowidx[nnz] = row;
			elts[nnz++] = eltsA[r];
		}
		for (int r = colptrC[col]; r < colptrC[col+1]; r++) {
			int row = rowidxC[r]+nx+my+1;
			rowidx[nnz] = row;
			elts[nnz++] = eltsC[r];
		}
	}
	colptr[nx_1] = nnz;
	assert(nnz == nnzA + nnzC);

	std::ofstream fd("augdump.dat");
	fd << std::scientific;
	fd.precision(16);
	fd << (nx + my + mz) << "\n";
	fd << nx_1 << "\n";
	fd << nnzA+nnzC << "\n";

   for( int i = 0; i <= nx_1; i++ )
      fd << colptr[i] << " ";
   fd << "\n";
   for( int i = 0; i < nnz; i++ )
      fd << rowidx[i] << " ";
   fd << "\n";
   for( int i = 0; i < nnz; i++ )
      fd << elts[i] << " ";
   fd << "\n";

   std::cout << "finished dumping aug\n";
}

Data* sFactory::makeData()
{
#ifdef TIMING
   double t2 = MPI_Wtime();
#endif

   StochGenMatrixHandle A(tree->createA());
   StochVectorHandle b(tree->createb());

   StochGenMatrixHandle C(tree->createC());
   StochVectorHandle clow(tree->createclow());
   StochVectorHandle iclow(tree->createiclow());
   StochVectorHandle cupp(tree->createcupp());
   StochVectorHandle icupp(tree->createicupp());

   StochSymMatrixHandle Q(tree->createQ());
   StochVectorHandle c(tree->createc());

   StochVectorHandle xlow(tree->createxlow());
   StochVectorHandle ixlow(tree->createixlow());
   StochVectorHandle xupp(tree->createxupp());
   StochVectorHandle ixupp(tree->createixupp());

#ifdef TIMING
   MPI_Barrier( tree->getCommWrkrs() );
   t2 = MPI_Wtime() - t2;
   if ( PIPS_MPIgetRank(tree->getCommWorkers() == 0) )
         std::cout << "IO second part took " << t2 << " sec\n";
#endif

   data = new sData(tree, c, Q, xlow, ixlow, xupp, ixupp,
         A, b, C, clow, iclow, cupp, icupp );
   return data;
}

// TODO adjust this for hierarchical approach
Variables* sFactory::makeVariables( Data * prob_in )
{
  sData* prob = dynamic_cast<sData*>(prob_in);

  OoqpVectorHandle x      = OoqpVectorHandle( makePrimalVector() );
  OoqpVectorHandle s      = OoqpVectorHandle( makeDualZVector() );
  OoqpVectorHandle y      = OoqpVectorHandle( makeDualYVector() );
  OoqpVectorHandle z      = OoqpVectorHandle( makeDualZVector() );
  OoqpVectorHandle v      = OoqpVectorHandle( makePrimalVector() );
  OoqpVectorHandle gamma  = OoqpVectorHandle( makePrimalVector() );
  OoqpVectorHandle w      = OoqpVectorHandle( makePrimalVector() );
  OoqpVectorHandle phi    = OoqpVectorHandle( makePrimalVector() );
  OoqpVectorHandle t      = OoqpVectorHandle( makeDualZVector() );
  OoqpVectorHandle lambda = OoqpVectorHandle( makeDualZVector() );
  OoqpVectorHandle u      = OoqpVectorHandle( makeDualZVector() );
  OoqpVectorHandle pi     = OoqpVectorHandle( makeDualZVector() );

  sVars* vars = new sVars( tree, x, s, y, z,
			   v, gamma, w, phi,
			   t, lambda, u, pi,
			   prob->ixlow, prob->ixlow->numberOfNonzeros(),
			   prob->ixupp, prob->ixupp->numberOfNonzeros(),
			   prob->iclow, prob->iclow->numberOfNonzeros(),
			   prob->icupp, prob->icupp->numberOfNonzeros());
  registeredVars.push_back(vars);
  return vars;
}

Residuals* sFactory::makeResiduals( Data * prob_in )
{
  sData* prob = dynamic_cast<sData*>(prob_in);
  resid =  new sResiduals( tree, prob->ixlow, prob->ixupp, prob->iclow, prob->icupp);
  return resid;
}

LinearSystem* sFactory::makeLinsys( Data* )
{
   if( pips_options::getBoolParameter( "HIERARCHICAL" ) )
      linsys = newLinsysRootHierarchical();
   else
      linsys = newLinsysRoot();

   return linsys;
}

OoqpVector* sFactory::makePrimalVector() const
{
   assert( !la );
   return tree->newPrimalVector();
}

OoqpVector* sFactory::makeDualYVector() const
{
   assert( !la );
   return tree->newDualYVector();
}

OoqpVector* sFactory::makeDualZVector() const
{
   assert( !la );
   return tree->newDualZVector();
}

OoqpVector* sFactory::makeRhs() const
{
   assert( !la );
   return tree->newRhs();
}

void sFactory::iterateStarted()
{
  iterTmMonitor.recIterateTm_start();
  tree->startMonitors();
}

void sFactory::iterateEnded()
{
  tree->stopMonitors();

  if(tree->balanceLoad()) {
    printf("Should not get here! OMG OMG OMG\n");
  }
  //logging and monitoring
  iterTmMonitor.recIterateTm_stop();
  m_tmTotal += iterTmMonitor.tmIterate;

  if( PIPS_MPIgetRank() == 0 )
  {
#ifdef TIMING
    extern double g_iterNumber;
    printf("TIME %g SOFAR %g ITER %d\n", iterTmMonitor.tmIterate, m_tmTotal, (int)g_iterNumber);
    //#elseif STOCH_TESTING
    //printf("ITERATION WALLTIME: iter=%g  Total=%g\n", iterTmMonitor.tmIterate, m_tmTotal);
#endif
  }
}
