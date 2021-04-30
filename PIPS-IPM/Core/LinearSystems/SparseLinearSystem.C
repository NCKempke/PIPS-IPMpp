/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "SparseLinearSystem.h"
#include "DoubleLinearSolver.h"
#include "SparseSymMatrix.h"

SparseLinearSystem::SparseLinearSystem(ProblemFactory* factory_in, Problem* problem, SparseSymMatrix* Mat_in, DoubleLinearSolver* solver_in) : LinearSystem(
      factory_in, problem), solver(solver_in) {
   SpReferTo(Mat, Mat_in);
}


void SparseLinearSystem::put_primal_diagonal() {
   assert(primal_diagonal);
   Mat->atPutDiagonal(0, *primal_diagonal);
}


void SparseLinearSystem::put_dual_inequalites_diagonal() {
   assert(nomegaInv);
   Mat->atPutDiagonal(nx + my, *nomegaInv);
}


void SparseLinearSystem::solveCompressed(Vector<double>& arhs) {
   //printf("-----\n");arhs.writeToStream(cout);
   solver->solve(arhs);
   //printf("~~~~~\n");arhs.writeToStream(cout);
}


void SparseLinearSystem::factorize(Problem* prob, Variables* vars) {
   this->LinearSystem::factorize(prob, vars);
   solver->matrixChanged();
}

SparseLinearSystem::~SparseLinearSystem() {
   delete solver;
}
