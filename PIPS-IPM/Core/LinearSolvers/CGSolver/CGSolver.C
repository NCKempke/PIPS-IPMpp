#include "CGSolver.h"
#include "DenseVector.hpp"

extern int print_level;

#define EPS 2.220e-16

CGSolver::CGSolver(MatTimesVec* A, MatTimesVec* M1, MatTimesVec* M2) : DoubleIterativeLinearSolver(A, M1, M2) {
   tol = 5.0e-13;
   maxit = 10;
   iter = -1;
   flag = -1;
   //firstSolve = 1;
};

void CGSolver::solve(Vector<double>& rhs_) {
   auto& b = dynamic_cast<DenseVector<double>&>(rhs_);
   int n = b.length();

   int flag, imin;
   int stag, maxmsteps, maxstagsteps, moresteps;
   double normrmin;
   double alpha, beta, rho, rho1, pq;
   int iter;

   double n2b = b.two_norm();
   double tolb = n2b * tol;

   if (tmpVec1.size() < static_cast<size_t>(n)) {
      tmpVec1.resize(n);
      tmpVec2.resize(n);
      tmpVec3.resize(n);
      tmpVec4.resize(n);
      tmpVec5.resize(n);
      tmpVec6.resize(n);
   }

   DenseVector<double> x(tmpVec1.data(), n);      //iterate
   DenseVector<double> r(tmpVec2.data(), n);      //residual
   DenseVector<double> xmin(tmpVec3.data(), n);   //minimal residual iterate
   DenseVector<double> y(tmpVec4.data(), n);      //work vectors
   DenseVector<double> z(tmpVec5.data(), n);      //work vectors
   DenseVector<double> p(tmpVec6.data(), n);
   //if(firstSolve)
   //  //initial guess is 0, the previous found solution otherwise
   x.setToZero();

   xmin.copyFrom(x);
   y.setToZero();
   flag = 1;
   imin = 0;

   r.copyFrom(b);
   applyA(1.0, r, -1.0, x);
   double normr = r.two_norm();
   double normr_act = normr;

   maxit = n / 2 + 20;
   if (normr < tolb) {
      printf("lucky!!!!!!!!!!!!\n");
      //initial guess is good enough
      b.copyFrom(x);
      return;
   }
   normrmin = normr;
   rho = 1.0;
   stag = 0;
   maxmsteps = std::min(std::min(n / 50, 5), n - maxit);
   maxstagsteps = 2;
   moresteps = 0;
   iter = 0;

   //////////////////////////////////////////////////////////////////
   // loop over maxit iterations
   //////////////////////////////////////////////////////////////////
   int ii = 0;
   while (ii < maxit) {
      applyM1(0.0, y, 1.0, r);
      applyM2(0.0, z, 1.0, y);

      rho1 = rho;
      rho = r.dotProductWith(z);
      if (rho == 0.0) {
         flag = 4;
         break;
      }

      if (ii == 0)
         p.copyFrom(z);
      else {
         beta = rho / rho1;
         if (beta == 0.0) {
            flag = 4;
            break;
         }
         p.scale(beta);
         p.add(1.0, z); // p=z + beta*p
      }

      DenseVector<double>& q = y;
      applyA(0.0, q, 1.0, p); //q=A*p
      pq = p.dotProductWith(q);
      if (pq <= 0) {
         flag = 4;
         break;
      }
      alpha = rho / pq;

      //check for stagnation
      if (p.two_norm() * fabs(alpha) < EPS * x.two_norm())
         stag++;
      else
         stag = 0;

      //---- updates ----
      x.add(alpha, p);
      r.add(-alpha, q);
      normr = r.two_norm();
      normr_act = normr;

      //printf("stag=%d  maxstagsteps=%d moresteps=%d  normr=%g\n",
      //   stag, maxstagsteps, moresteps, normr);
      // check for convergence
      if (normr <= tolb || stag >= maxstagsteps || moresteps) {
         r.copyFrom(b);
         applyA(1.0, r, -1.0, x);
         normr_act = r.two_norm();

         if (normr_act <= tolb) {
            flag = 0;
            iter = 1 + ii;
            break;
         }
         else {
            if (stag >= maxstagsteps && moresteps == 0)
               stag = 0;
            moresteps++;
            if (moresteps >= maxmsteps) {
               flag = 3;
               iter = 1 + ii;
               break;
            }
         }
      }

      if (normr_act < normrmin) {
         normrmin = normr_act;
         xmin.copyFrom(x);
         imin = ii;
      }

      if (stag >= maxstagsteps) {
         flag = 3;
         break;
      }

      ii++;
   }

   //////////////////////////////////////////////////////////
   // status/error output
   /////////////////////////////////////////////////////////
   if (flag == 0) {
      double relres = normr_act / n2b;
      printf("CG converged: actual normResid=%g relResid=%g iter=%d\n", normr_act, relres, iter);
   }
   else {
      double relres = normr_act / n2b;
      r.copyFrom(b);
      applyA(1.0, r, -1.0, xmin);
      normr = r.two_norm();
      if (normr < normr_act) {
         x.copyFrom(xmin);
         iter = imin;
         relres = normr / n2b;
      }
      else {
         iter = ii;
         relres = normr_act / n2b;
      }

      if (print_level >= 1) {
         printf("CG did not NOT converged after %d  max of %d iters were made.\n", iter, ii);
         printf("\t - Error code %d\n\t - Act res=%g\n\t - Rel res=%g %g\n\n", flag, normr, relres, normrmin);
      }

   }
   b.copyFrom(x);
}