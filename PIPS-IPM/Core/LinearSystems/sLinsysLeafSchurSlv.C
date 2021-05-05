/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysLeafSchurSlv.h"
#include "DistributedQP.hpp"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"
#include "PardisoSolver.h"
#include "PardisoSchurSolver.h"

extern int gLackOfAccuracy;

/**
 * Computes U = Gi * inv(H_i) * Gi^T.
 *        [ R 0 0 ]
 * Gi^T = [ A 0 0 ]
 *        [ C 0 0]
 *
 * A and C are the recourse eq. and ineq. matrices, R is the cross
 * Hessian term.
 */
void sLinsysLeafSchurSlv::addTermToDenseSchurCompl(DistributedQP* prob, DenseSymMatrix& SC) {
   SparseGenMatrix& A = prob->getLocalA();
   SparseGenMatrix& C = prob->getLocalC();
   SparseGenMatrix& F = prob->getLocalF();
   SparseGenMatrix& G = prob->getLocalG();
   SparseGenMatrix& R = prob->getLocalCrossHessian();

   //if(!gLackOfAccuracy && !switchedToSafeSlv) {
   PardisoSchurSolver* scSolver = dynamic_cast<PardisoSchurSolver*>(solver.get());
   scSolver->schur_solve(R, A, C, F, G, SC);
   //} else {
   ////cout << "\tdefaulting to sLinsysLeaf::addTermToDenseSchurCompl ...";
   //sLinsysLeaf::addTermToDenseSchurCompl(prob, SC);
   ////cout << "done" << endl;
   //}
}

void sLinsysLeafSchurSlv::addTermToSparseSchurCompl(DistributedQP* prob, SparseSymMatrix& SC) {
   SparseGenMatrix& A = prob->getLocalA();
   SparseGenMatrix& C = prob->getLocalC();
   SparseGenMatrix& F = prob->getLocalF();
   SparseGenMatrix& G = prob->getLocalG();
   SparseGenMatrix& R = prob->getLocalCrossHessian();

   PardisoSchurSolver* scSolver = dynamic_cast<PardisoSchurSolver*>(solver.get());
   scSolver->schur_solve_sparse(R, A, C, F, G, SC);
}

void sLinsysLeafSchurSlv::factor2(DistributedQP* prob, Variables* vars) {
   // if(gLackOfAccuracy) {
   //   cout << "sLinsysLeafSchurSlv -> accuracy lost, switching to vanilla PARDISO" << endl;
   //   delete solver;
   //   //cout << "\tsolver deleted\n";
   //   SparseSymMatrix* kktsp = dynamic_cast<SparseSymMatrix*>(kkt);
   //   solver = new PardisoSolver(kktsp);
   //   //solver = new Ma57Solver(kktsp);
   //   //cout << "\tnew solver created." << endl;
   //   switchedToSafeSlv=true;
   // }

   DistributedLeafLinearSystem::factor2(prob, vars);
}
