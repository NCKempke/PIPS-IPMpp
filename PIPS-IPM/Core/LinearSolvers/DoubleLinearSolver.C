#include "DoubleLinearSolver.h"
#include "SimpleVector.hpp"

DoubleIterativeLinearSolver::DoubleIterativeLinearSolver(MatTimesVec* Ain, MatTimesVec* M1in, MatTimesVec* M2in) : A(Ain), ML(M1in), MR(M2in) {

}

void DoubleIterativeLinearSolver::applyA(double beta, Vector<double>& res, double alpha, Vector<double>& x) {
   A->doIt(beta, res, alpha, x);
}


void DoubleIterativeLinearSolver::applyM1(double beta, Vector<double>& res, double alpha, Vector<double>& x) {
   ML->doIt(beta, res, alpha, x);
}

void DoubleIterativeLinearSolver::applyM2(double beta, Vector<double>& res, double alpha, Vector<double>& x) {
   if (nullptr == MR) {
      //just a identity precond
      if (beta == 0.0) {
         if (alpha == 1.0) {
            res.copyFrom(x);
            return;
         }
         //alpha not zero
         res.setToZero();

      }
      else if (beta != 1.0)
         res.scale(beta);

      //beta not 0.0 and alpha not 1.0
      if (alpha != 0)
         res.axpy(alpha, x);

   }
   else
      MR->doIt(beta, res, alpha, x);
}
