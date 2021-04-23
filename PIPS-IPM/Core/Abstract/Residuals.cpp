#include "Residuals.h"
#include "Variables.h"
#include "Problem.h"
#include "LinearAlgebraPackage.h"
#include "pipsdef.h"

#include <iostream>

Residuals::Residuals(const Residuals& residuals) : mResidualNorm{residuals.mResidualNorm}, mDualityGap{residuals.mDualityGap},
      primal_objective{residuals.primal_objective}, dual_objective{residuals.dual_objective}, nx{residuals.nx}, my{residuals.my}, mz{residuals.mz},
      nxupp{residuals.nxupp}, nxlow{residuals.nxlow}, mcupp{residuals.mcupp}, mclow{residuals.mclow} {
   ixlow = OoqpVectorHandle(residuals.ixlow->cloneFull());
   ixupp = OoqpVectorHandle(residuals.ixupp->cloneFull());
   iclow = OoqpVectorHandle(residuals.iclow->cloneFull());
   icupp = OoqpVectorHandle(residuals.icupp->cloneFull());

   rQ = OoqpVectorHandle(residuals.rQ->cloneFull());
   rA = OoqpVectorHandle(residuals.rA->cloneFull());
   rC = OoqpVectorHandle(residuals.rC->cloneFull());

   rz = OoqpVectorHandle(residuals.rz->cloneFull());

   rt = OoqpVectorHandle(residuals.rt->cloneFull());
   rlambda = OoqpVectorHandle(residuals.rlambda->cloneFull());

   ru = OoqpVectorHandle(residuals.ru->cloneFull());
   rpi = OoqpVectorHandle(residuals.rpi->cloneFull());

   rv = OoqpVectorHandle(residuals.rv->cloneFull());
   rgamma = OoqpVectorHandle(residuals.rgamma->cloneFull());

   rw = OoqpVectorHandle(residuals.rw->cloneFull());
   rphi = OoqpVectorHandle(residuals.rphi->cloneFull());
}

double updateNormAndPrint(double norm, const OoqpVector& vec, bool print, std::string&& name) {
   const double infnorm = vec.infnorm();

   if (print) {
      const double twonorm = vec.twonorm();

      if (0 == PIPS_MPIgetRank())
         std::cout << name << " infnorm = " << infnorm << " | twonorm = " << twonorm << "\n";
   }

   return std::max(norm, infnorm);
}

void Residuals::evaluate(Problem& problem, Variables* iterate, bool print_resids) {
#ifdef TIMING
   print_resids = true;
#endif
   const int myRank = PIPS_MPIgetRank();

   double norm = 0.0, gap = 0.0;
   primal_objective = 0.0;
   dual_objective = 0.0;

   primal_objective = problem.objective_value(iterate);

   /*** rQ = Qx + g - A^T y - C^T z - gamma + phi ***/
   problem.objective_gradient(iterate, *rQ); // Qx + g

   // contribution calculate x^T (g + Qx) to duality gap */
   gap = rQ->dotProductWith(*iterate->x);

   problem.ATransmult(1.0, *rQ, -1.0, *iterate->y);
   problem.CTransmult(1.0, *rQ, -1.0, *iterate->z);

   iterate->gamma->selectNonZeros(*ixlow);
   iterate->phi->selectNonZeros(*ixupp);
   if (nxlow > 0)
      rQ->axpy(-1.0, *iterate->gamma);
   if (nxupp > 0)
      rQ->axpy(1.0, *iterate->phi);

   norm = updateNormAndPrint(norm, *rQ, print_resids, "rQ");

   /*** rA = Ax - b ***/
   problem.getbA(*rA);
   problem.Amult(-1.0, *rA, 1.0, *iterate->x);

   norm = updateNormAndPrint(norm, *rA, print_resids, "rA");

   // contribution -d^T y to duality gap
   const double ba_y = problem.bA->dotProductWith(*iterate->y);
   gap -= ba_y;
   dual_objective += ba_y;

   /*** rC = Cx - s ***/
   rC->copyFrom(*iterate->s);
   problem.Cmult(-1.0, *rC, 1.0, *iterate->x);

   norm = updateNormAndPrint(norm, *rC, print_resids, "rC");

   /*** rz = z - lambda + pi ***/
   rz->copyFrom(*iterate->z);

   if (mclow > 0)
      rz->axpy(-1.0, *iterate->lambda);
   if (mcupp > 0)
      rz->axpy(1.0, *iterate->pi);

   norm = updateNormAndPrint(norm, *rz, print_resids, "rz");

   if (mclow > 0) {
      /*** rt = s - d - t ***/
      rt->copyFrom(*iterate->s);
      rt->axpy(-1.0, problem.slowerBound());
      rt->selectNonZeros(*iclow);
      rt->axpy(-1.0, *iterate->t);

      // contribution - d^T lambda to duality gap
      const double bl_lambda = problem.bl->dotProductWith(*iterate->lambda);
      gap -= bl_lambda;
      dual_objective += bl_lambda;

      norm = updateNormAndPrint(norm, *rt, print_resids, "rt");
   }

   if (mcupp > 0) {
      /*** ru = s - f + u ***/
      ru->copyFrom(*iterate->s);
      ru->axpy(-1.0, problem.supperBound());
      ru->selectNonZeros(*icupp);
      ru->axpy(1.0, *iterate->u);

      // contribution - f^T pi to duality gap
      const double bu_pi = problem.bu->dotProductWith(*iterate->pi);
      gap += bu_pi;
      dual_objective -= bu_pi;

      norm = updateNormAndPrint(norm, *ru, print_resids, "ru");
   }

   if (nxlow > 0) {
      /*** rv = x - lx - v ***/
      rv->copyFrom(*iterate->x);
      rv->axpy(-1.0, problem.xlowerBound());
      rv->selectNonZeros(*ixlow);
      rv->axpy(-1.0, *iterate->v);

      norm = updateNormAndPrint(norm, *rv, print_resids, "rv");
      // contribution - lx^T gamma to duality gap
      const double blx_gamma = problem.blx->dotProductWith(*iterate->gamma);
      gap -= blx_gamma;
      dual_objective += blx_gamma;

   }

   if (nxupp > 0) {
      /*** rw = x - ux + w ***/
      rw->copyFrom(*iterate->x);
      rw->axpy(-1.0, problem.xupperBound());
      rw->selectNonZeros(*ixupp);
      rw->axpy(1.0, *iterate->w);

      norm = updateNormAndPrint(norm, *rw, print_resids, "rw");

      // contribution + bu^T phi to duality gap
      const double bux_phi = problem.bux->dotProductWith(*iterate->phi);
      gap += bux_phi;
      dual_objective -= bux_phi;
   }

   mDualityGap = gap;
   mResidualNorm = norm;

   if (print_resids && myRank == 0) {
      std::cout << "Norm residuals: " << mResidualNorm << "\tduality gap: " << mDualityGap << "\n";
   }
}

double Residuals::recomputeResidualNorm() {
   mResidualNorm = 0.0;

   double componentNorm = 0.0;
   componentNorm = rQ->infnorm();

   if (componentNorm > mResidualNorm)
      mResidualNorm = componentNorm;

   componentNorm = rA->infnorm();
   if (componentNorm > mResidualNorm)
      mResidualNorm = componentNorm;

   componentNorm = rC->infnorm();
   if (componentNorm > mResidualNorm)
      mResidualNorm = componentNorm;

   if (mclow > 0) {
      componentNorm = rt->infnorm();
      if (componentNorm > mResidualNorm)
         mResidualNorm = componentNorm;
   }

   if (mcupp > 0) {
      componentNorm = ru->infnorm();
      if (componentNorm > mResidualNorm)
         mResidualNorm = componentNorm;
   }

   componentNorm = rz->infnorm();
   if (componentNorm > mResidualNorm)
      mResidualNorm = componentNorm;

   if (nxlow > 0) {
      componentNorm = rv->infnorm();
      if (componentNorm > mResidualNorm)
         mResidualNorm = componentNorm;
   }

   if (nxupp > 0) {
      componentNorm = rw->infnorm();
      if (componentNorm > mResidualNorm)
         mResidualNorm = componentNorm;
   }
   return mResidualNorm;
}

void Residuals::add_r3_xz_alpha(const Variables* variables, double alpha) {
   if (mclow > 0)
      rlambda->axzpy(1.0, *variables->t, *variables->lambda);
   if (mcupp > 0)
      rpi->axzpy(1.0, *variables->u, *variables->pi);
   if (nxlow > 0)
      rgamma->axzpy(1.0, *variables->v, *variables->gamma);
   if (nxupp > 0)
      rphi->axzpy(1.0, *variables->w, *variables->phi);

   if (alpha != 0.0) {
      if (mclow > 0)
         rlambda->addSomeConstants(alpha, *iclow);
      if (mcupp > 0)
         rpi->addSomeConstants(alpha, *icupp);
      if (nxlow > 0)
         rgamma->addSomeConstants(alpha, *ixlow);
      if (nxupp > 0)
         rphi->addSomeConstants(alpha, *ixupp);
   }
}

void Residuals::set_r3_xz_alpha(const Variables* vars, double alpha) {
   this->clear_r3();
   this->add_r3_xz_alpha(vars, alpha);
}

void Residuals::clear_r3() {
   if (mclow > 0)
      rlambda->setToZero();
   if (mcupp > 0)
      rpi->setToZero();
   if (nxlow > 0)
      rgamma->setToZero();
   if (nxupp > 0)
      rphi->setToZero();
}

void Residuals::clear_r1r2() {
   rQ->setToZero();
   rA->setToZero();
   rC->setToZero();
   rz->setToZero();
   if (nxlow > 0)
      rv->setToZero();
   if (nxupp > 0)
      rw->setToZero();
   if (mclow > 0)
      rt->setToZero();
   if (mcupp > 0)
      ru->setToZero();
}

void Residuals::project_r3(double rmin, double rmax) {
   if (mclow > 0) {
      rlambda->gondzioProjection(rmin, rmax);
      rlambda->selectNonZeros(*iclow);
   }
   if (mcupp > 0) {
      rpi->gondzioProjection(rmin, rmax);
      rpi->selectNonZeros(*icupp);
   }
   if (nxlow > 0) {
      rgamma->gondzioProjection(rmin, rmax);
      rgamma->selectNonZeros(*ixlow);
   }
   if (nxupp > 0) {
      rphi->gondzioProjection(rmin, rmax);
      rphi->selectNonZeros(*ixupp);
   }

}


int Residuals::validNonZeroPattern() {
   if (nxlow > 0 && (!rv->matchesNonZeroPattern(*ixlow) || !rgamma->matchesNonZeroPattern(*ixlow))) {
      return 0;
   }

   if (nxupp > 0 && (!rw->matchesNonZeroPattern(*ixupp) || !rphi->matchesNonZeroPattern(*ixupp))) {
      return 0;
   }
   if (mclow > 0 && (!rt->matchesNonZeroPattern(*iclow) || !rlambda->matchesNonZeroPattern(*iclow))) {
      return 0;
   }

   if (mcupp > 0 && (!ru->matchesNonZeroPattern(*icupp) || !rpi->matchesNonZeroPattern(*icupp))) {
      return 0;
   }

   return 1;
}

void Residuals::copyFrom(const Residuals& residuals) {
   mResidualNorm = residuals.mResidualNorm;
   mDualityGap = residuals.mDualityGap;
   m = residuals.m;
   n = residuals.n;

   nx = residuals.nx;
   my = residuals.my;
   mz = residuals.mz;

   nxupp = residuals.nxupp;
   ixupp->copyFrom(*residuals.ixupp);

   nxlow = residuals.nxlow;
   ixlow->copyFrom(*residuals.ixlow);

   mcupp = residuals.mcupp;
   icupp->copyFrom(*residuals.icupp);

   mclow = residuals.mclow;
   iclow->copyFrom(*residuals.iclow);

   rQ->copyFrom(*residuals.rQ);
   rA->copyFrom(*residuals.rA);
   rC->copyFrom(*residuals.rC);

   rz->copyFrom(*residuals.rz);
   rv->copyFrom(*residuals.rv);
   rw->copyFrom(*residuals.rw);
   rt->copyFrom(*residuals.rt);
   ru->copyFrom(*residuals.ru);
   rgamma->copyFrom(*residuals.rgamma);
   rphi->copyFrom(*residuals.rphi);
   rlambda->copyFrom(*residuals.rlambda);
   rpi->copyFrom(*residuals.rpi);
}

void Residuals::writeToStream(std::ostream& out) {
   /*
   printf("--------------rQ\n");
   rQ->writeToStream(out);printf("---------------------------\n");

   printf("rA\n");
   rA->writeToStream(out);printf("---------------------------\n");
   */
   printf("rC\n");
   rC->writeToStream(out);
   printf("---------------------------\n");


   printf("rz\n");
   rz->writeToStream(out);
   printf("---------------------------\n");
}
