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


void SparseLinearSystem::putXDiagonal(OoqpVector& xdiag) {
   Mat->atPutDiagonal(0, xdiag);
}


void SparseLinearSystem::putZDiagonal(OoqpVector& zdiag) {
   Mat->atPutDiagonal(nx + my, zdiag);
   //zdiag.writeToStream(cout);
   //!assert(false);

   /*//!log
   printf("QpGenSparseLinsys::putZDiagonal\n");
   SparseSymMatrix& M = *Mat;

   int* krowM = M.krowM();
   int* jcolM = M.jcolM();
   double* dM = M.M();

   int nn; M.getSize(nn, nn);
   for(int i=0; i<nn; i++) {
     printf("row %d \n\t", i);

     for(int j=krowM[i];j<krowM[i+1]; j++)
       printf("%9d ", jcolM[j]);
     printf("\n\t");

     for(int j=krowM[i];j<krowM[i+1]; j++)
       printf("%9.6f ", dM[j]);
     printf("\n");
   }
   */
}


void SparseLinearSystem::solveCompressed(OoqpVector& arhs) {
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
