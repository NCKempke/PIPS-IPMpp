#include "PCGSolver.h"
#include "DenseVector.hpp"

extern int print_level;

#define EPS 2.220e-16

PCGSolver::PCGSolver(MatTimesVec* H_in, MatTimesVec* P_in, MatTimesVec* At_in, int n_in, int m_in) : DoubleIterativeLinearSolver(H_in, P_in, nullptr),
      n(n_in), m(m_in) {
   tol = 4.0e-12;
   maxit = 10;
   iter = -1;
   flag = -1;
   At = At_in;
};

PCGSolver::~PCGSolver() {

   delete[] tmpVec1;

   delete[] tmpVec2;
   delete[] tmpVec3;
   delete[] tmpVec4;
   delete[] tmpVec5;
   delete[] tmpVec6;
}

void PCGSolver::solve(Vector<double>& rhs_) {
   solvefull(rhs_);
}

void PCGSolver::solvefull(Vector<double>& rhs_) {
   auto& b = dynamic_cast<DenseVector<double>&>(rhs_);
   assert(n + m == b.length());

   int flag, imin;
   int stagsteps, maxstagsteps;
   double normr, normr_act, normrmin;
   double alpha, beta, rg, pHp;

   double n2b = b.two_norm();
   double tolb = n2b * tol;

   if (tmpVec1 == nullptr)
      tmpVec1 = new double[n + m];
   if (tmpVec2 == nullptr)
      tmpVec2 = new double[n + m];
   if (tmpVec3 == nullptr)
      tmpVec3 = new double[n];
   if (tmpVec4 == nullptr)
      tmpVec4 = new double[n];
   if (tmpVec5 == nullptr)
      tmpVec5 = new double[n];
   if (tmpVec6 == nullptr)
      tmpVec6 = new double[n + m];

   DenseVector<double> xy(tmpVec1, n + m);      //iterate
   DenseVector<double> auxnm(tmpVec2, n + m);      //auxiliary
   DenseVector<double> xmin(tmpVec3, n);   //minimal residual iterate
   DenseVector<double> g(tmpVec4, n);      //work vectors
   DenseVector<double> p(tmpVec5, n);
   DenseVector<double> res(tmpVec6, n + m);

   DenseVector<double> x(&xy[0], n); //y-part of the iterate
   DenseVector<double> y(&xy[n], m); //x-part of the iterate
   DenseVector<double> rx(&res[0], n); //residual
   DenseVector<double> ry(&res[n], m); //residual corresponding to last m eqn
   DenseVector<double> auxn(&auxnm[0], n);
   DenseVector<double> auxm(&auxnm[n], m);
   //////////////////////////////////////////////////////////////////
   // Starting procedure
   /////////////////////////////////////////////////////////////////

   //find starting point x satisfying Ax=b_2
   applyM1(0.0, xy, 1.0, b);

   //compute the x-residual for the starting point rx=Hx-b_1
   rx.copyFromArray(&b[0]);
   applyA(-1.0, rx, 1.0, x);

   //find y such that it minimizes ||r-A'y||_Ginv
   //this is done by a preconditioner solve with rhs=[rx;0]
   ry.setToZero();
   applyM1(0.0, auxnm, 1.0, res);
   y.copyFromArray(&auxnm[n]);

   //remove A'y from residual
   At->doIt(1.0, rx, -1.0, y);

   //initialize projected residual g=Pr and update p=-g
   ry.setToZero();
   applyM1(0.0, auxnm, 1.0, res);
   g.copyFromArray(&auxnm[0]);
   p.copyFrom(g);
   p.negate();

   normr = rx.two_norm();
   rg = rx.dotProductWith(g);

   xmin.copyFrom(x);
   flag = 1;
   imin = 0;

   maxit = n / 2 + 10;
   if (normr < tolb) {
      //initial guess is good enough
      for (int i = 0; i < n; i++)
         b[i] = x[i];
      for (int i = n; i < n + m; i++)
         b[i] = y[i - n];
      return;
   }
   stagsteps = 0;
   maxstagsteps = 5;
   normrmin = normr;

   //////////////////////////////////////////////////////////////////
   // loop over maxit iterations
   //////////////////////////////////////////////////////////////////
   int ii = 0;
   while (ii < maxit) {
      ii++;
      // compute Hp and p'Hp
      DenseVector<double> Hp(&auxn[0], n);
      applyA(0.0, Hp, 1.0, p);
      pHp = p.dotProductWith(Hp);

      //check for negative curvature
      if (pHp < 0.0) {
         flag = 2;
         break;
      }

      alpha = rg / pHp;

      //update x=x+alpha*p and r=r+alpha*H*p
      x.add(alpha, p);
      rx.add(alpha, Hp);

      normr = rx.two_norm();
      ///////////////////////////////////////
      //convergence tests
      ///////////////////////////////////////
      if (normr <= tolb) {
         //compute actual residual
         DenseVector<double> rx_act(&auxnm[0], n);
         rx_act.copyFromArray(&b[0]);
         applyA(-1.0, rx_act, 1.0, x);
         normr_act = rx_act.two_norm();

         //if(normr_act/n2b<tolb) { flag=0; break; }
         {
            flag = 0;
            break;
         }
      }

      if (normr < normrmin) {
         imin = ii;
         xmin.copyFrom(x);
         normrmin = normr;
         stagsteps = 0;
      }
      else
         stagsteps++;
      //check for stagnation!!!
      if (stagsteps > maxstagsteps) {
         flag = 4;
         break;
      }
      //------- end convergence tests -------

      //compute y that minimizes ||r-A'y||_Ginv
      ry.setToZero();
      applyM1(0.0, auxnm, 1.0, res);
      y.copyFromArray(&auxnm[n]);
      //substract A'y from r
      At->doIt(1.0, rx, -1.0, y);

      //projected residual g=Pr
      ry.setToZero();
      applyM1(0.0, auxnm, 1.0, res);
      g.copyFromArray(&auxnm[0]);

      double rpgp = rx.dotProductWith(g);
      beta = rpgp / rg;

      //p = -g+beta*p
      p.scale(beta);
      p.add(-1.0, g);

      rg = rpgp;

      //rounding error
      if (rg < 0.0) {
         flag = 3;
         break;
      }
   }


   //////////////////////////////////////////////////////////
   // status/error output
   /////////////////////////////////////////////////////////
   if (flag == 0) {
      double relres = normr_act / n2b;
      printf("CG converged: actual normResid=%g relResid=%g iter=%d\n", normr_act, relres, ii);

      b.setToZero();
      for (int i = 0; i < n; i++)
         b[i] = x[i];
      for (int i = n; i < m + n; i++)
         b[i] = y[i - n];
   }
   else {

      if (flag == 4)
         x.copyFrom(xmin);
      //compute actual residual
      DenseVector<double> rx(&auxnm[0], n);
      rx.copyFromArray(&b[0]);
      applyA(1.0, rx, -1.0, x);
      normr_act = rx.two_norm();

      for (int i = 0; i < n; i++)
         b[i] = x[i];
      for (int i = n; i < m + n; i++)
         b[i] = y[i - n];

      if (print_level >= 1) {
         printf("Projected CG did not NOT converged after %d iters. The solution from iter %d was returned.\n", ii, imin);
         printf("\t - Error code %d\n\t - Act res=%g\n\t - Rel res=%g\n\n", flag, normr_act, normrmin);
      }

   }
   //b.copyFrom(x);
}

