/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGenSparseSeq.h"
#include "QuadraticProblem.h"
#include "SimpleVector.h"
#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"
#include "SparseLinearAlgebraPackage.h"
#include "QpGenVars.h"

//Problem * QpGenSparseSeq::create_problem()
//{
//  return new QuadraticProblem( la, nx, my, mz, nnzQ, nnzA, nnzC );
//}

void QpGenSparseSeq::joinRHS(OoqpVector& rhs_in, const OoqpVector& rhs1_in, const OoqpVector& rhs2_in, const OoqpVector& rhs3_in) const {
   SimpleVector& rhs = dynamic_cast<SimpleVector&>(rhs_in);
   const SimpleVector& rhs1 = dynamic_cast<const SimpleVector&>(rhs1_in);
   const SimpleVector& rhs2 = dynamic_cast<const SimpleVector&>(rhs2_in);
   const SimpleVector& rhs3 = dynamic_cast<const SimpleVector&>(rhs3_in);

   memcpy(&rhs[0], &rhs1[0], nx * sizeof(double));
   if (my > 0)
      memcpy(&rhs[nx], &rhs2[0], my * sizeof(double));
   if (mz > 0)
      memcpy(&rhs[nx + my], &rhs3[0], mz * sizeof(double));
}

void QpGenSparseSeq::separateVars(OoqpVector& x_in, OoqpVector& y_in, OoqpVector& z_in, const OoqpVector& vars_in) const {
   const SimpleVector& vars = dynamic_cast<const SimpleVector&>(vars_in);
   SimpleVector& x = dynamic_cast<SimpleVector&>(x_in);
   SimpleVector& y = dynamic_cast<SimpleVector&>(y_in);
   SimpleVector& z = dynamic_cast<SimpleVector&>(z_in);

   memcpy(&x[0], &vars[0], nx * sizeof(double));
   if (my > 0)
      memcpy(&y[0], &vars[nx], my * sizeof(double));
   if (mz > 0)
      memcpy(&z[0], &vars[nx + my], mz * sizeof(double));
}


Problem* QpGenSparseSeq::create_problem(double c_[], int krowQ[], int jcolQ[], double dQ[], double xlow_[], char ixlow_[], double xupp_[], char ixupp_[],
      int krowA[], int jcolA[], double dA[], double b_[], int krowC[], int jcolC[], double dC[], double clow_[], char iclow_[], double cupp_[],
      char icupp_[]) {
   // Objective funcition
   SimpleVectorHandle c(new SimpleVector(c_, nx));

   nnzQ = krowQ[nx];
   SparseSymMatrixHandle Q(new SparseSymMatrix(nx, nnzQ, krowQ, jcolQ, dQ));

   // Bounds on variables
   SimpleVectorHandle xlow(new SimpleVector(xlow_, nx));
   SimpleVectorHandle ixlow(new SimpleVector(nx));
   ixlow->copyFromArray(ixlow_);

   SimpleVectorHandle xupp(new SimpleVector(xupp_, nx));
   SimpleVectorHandle ixupp(new SimpleVector(nx));
   ixupp->copyFromArray(ixupp_);

   // Equality constraints
   nnzA = 0;
   if (my > 0)
      nnzA = krowA[my];
   SparseGenMatrixHandle A(new SparseGenMatrix(my, nx, nnzA, krowA, jcolA, dA));

   SimpleVectorHandle b(new SimpleVector(b_, my));

   // Inequality constraints
   nnzC = 0;
   if (mz > 0)
      nnzC = krowC[mz];
   SparseGenMatrixHandle C(new SparseGenMatrix(mz, nx, nnzC, krowC, jcolC, dC));

   SimpleVectorHandle clow(new SimpleVector(clow_, mz));
   SimpleVectorHandle iclow(new SimpleVector(mz));
   iclow->copyFromArray(iclow_);

   SimpleVectorHandle cupp(new SimpleVector(cupp_, mz));
   SimpleVectorHandle icupp(new SimpleVector(mz));
   icupp->copyFromArray(icupp_);

   QuadraticProblem* data = new QuadraticProblem(SparseLinearAlgebraPackage::soleInstance(), c, Q, xlow, ixlow, xupp, ixupp, A, b, C, clow, iclow,
         cupp, icupp);

   return data;
}