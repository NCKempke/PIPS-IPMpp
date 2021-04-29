/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGenSparseSeq.h"
#include "QP.hpp"
#include "SimpleVector.h"
#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"
#include "SparseLinearAlgebraPackage.h"
#include "Variables.h"

//Problem * QpGenSparseSeq::create_problem()
//{
//  return new QP( la, nx, my, mz, nnzQ, nnzA, nnzC );
//}

void QpGenSparseSeq::join_right_hand_side(OoqpVector& rhs_in, const OoqpVector& rhs1_in, const OoqpVector& rhs2_in, const OoqpVector& rhs3_in) const {
   SimpleVector<double>& rhs = dynamic_cast<SimpleVector<double>&>(rhs_in);
   const SimpleVector<double>& rhs1 = dynamic_cast<const SimpleVector<double>&>(rhs1_in);
   const SimpleVector<double>& rhs2 = dynamic_cast<const SimpleVector<double>&>(rhs2_in);
   const SimpleVector<double>& rhs3 = dynamic_cast<const SimpleVector<double>&>(rhs3_in);

   memcpy(&rhs[0], &rhs1[0], nx * sizeof(double));
   if (my > 0)
      memcpy(&rhs[nx], &rhs2[0], my * sizeof(double));
   if (mz > 0)
      memcpy(&rhs[nx + my], &rhs3[0], mz * sizeof(double));
}

void QpGenSparseSeq::separate_variables(OoqpVector& x_in, OoqpVector& y_in, OoqpVector& z_in, const OoqpVector& vars_in) const {
   const SimpleVector<double>& vars = dynamic_cast<const SimpleVector<double>&>(vars_in);
   SimpleVector<double>& x = dynamic_cast<SimpleVector<double>&>(x_in);
   SimpleVector<double>& y = dynamic_cast<SimpleVector<double>&>(y_in);
   SimpleVector<double>& z = dynamic_cast<SimpleVector<double>&>(z_in);

   memcpy(&x[0], &vars[0], nx * sizeof(double));
   if (my > 0)
      memcpy(&y[0], &vars[nx], my * sizeof(double));
   if (mz > 0)
      memcpy(&z[0], &vars[nx + my], mz * sizeof(double));
}


Problem*
QpGenSparseSeq::create_problem(double c_[], int krowQ[], int jcolQ[], double dQ[], double xlow_[], char ixlow_[], double xupp_[], char ixupp_[],
      int krowA[], int jcolA[], double dA[], double b_[], int krowC[], int jcolC[], double dC[], double clow_[], char iclow_[], double cupp_[],
      char icupp_[]) {
   // objective function
   SimpleVector<double> c(c_, nx);
   nnzQ = krowQ[nx];
   SparseSymMatrixHandle Q(new SparseSymMatrix(nx, nnzQ, krowQ, jcolQ, dQ));

   // bound constraints
   SimpleVector<double> xlow(xlow_, nx);
   SimpleVector<double> ixlow(nx);
   ixlow.copyFromArray(ixlow_);
   SimpleVector<double> xupp(xupp_, nx);
   SimpleVector<double> ixupp(nx);
   ixupp.copyFromArray(ixupp_);

   // equality constraints
   nnzA = (my > 0) ? krowA[my] : 0;
   SparseGenMatrixHandle A(new SparseGenMatrix(my, nx, nnzA, krowA, jcolA, dA));
   SimpleVector<double> b(b_, my);

   // inequality constraints
   nnzC = (mz > 0) ? krowC[mz] : 0;
   SparseGenMatrixHandle C(new SparseGenMatrix(mz, nx, nnzC, krowC, jcolC, dC));
   SimpleVector<double> clow(clow_, mz);
   SimpleVector<double> iclow(mz);
   iclow.copyFromArray(iclow_);
   SimpleVector<double> cupp(cupp_, mz);
   SimpleVector<double> icupp(mz);
   icupp.copyFromArray(icupp_);

   QP* data = new QP(SparseLinearAlgebraPackage::soleInstance(), &c, Q, &xlow, &ixlow, &xupp, &ixupp, A, &b, C, &clow, &iclow, &cupp, &icupp);
   return data;
}