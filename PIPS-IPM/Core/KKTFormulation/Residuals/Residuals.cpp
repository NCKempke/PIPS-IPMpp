#include "Residuals.h"
#include "Variables.h"
#include "Problem.h"
#include "pipsdef.h"

#include <iostream>

Residuals::Residuals(const Residuals& residuals) : mResidualNorm{residuals.mResidualNorm}, mDualityGap{residuals.mDualityGap},
      primal_objective{residuals.primal_objective}, dual_objective{residuals.dual_objective}, nx{residuals.nx}, my{residuals.my}, mz{residuals.mz},
      nxupp{residuals.nxupp}, nxlow{residuals.nxlow}, mcupp{residuals.mcupp}, mclow{residuals.mclow} {
   ixlow = SmartPointer<Vector<double> >(residuals.ixlow->cloneFull());
   ixupp = SmartPointer<Vector<double> >(residuals.ixupp->cloneFull());
   iclow = SmartPointer<Vector<double> >(residuals.iclow->cloneFull());
   icupp = SmartPointer<Vector<double> >(residuals.icupp->cloneFull());

   lagrangian_gradient = SmartPointer<Vector<double> >(residuals.lagrangian_gradient->cloneFull());
   rA = SmartPointer<Vector<double> >(residuals.rA->cloneFull());
   rC = SmartPointer<Vector<double> >(residuals.rC->cloneFull());

   rz = SmartPointer<Vector<double> >(residuals.rz->cloneFull());

   rt = SmartPointer<Vector<double> >(residuals.rt->cloneFull());
   rlambda = SmartPointer<Vector<double> >(residuals.rlambda->cloneFull());

   ru = SmartPointer<Vector<double> >(residuals.ru->cloneFull());
   rpi = SmartPointer<Vector<double> >(residuals.rpi->cloneFull());

   rv = SmartPointer<Vector<double> >(residuals.rv->cloneFull());
   rgamma = SmartPointer<Vector<double> >(residuals.rgamma->cloneFull());

   rw = SmartPointer<Vector<double> >(residuals.rw->cloneFull());
   rphi = SmartPointer<Vector<double> >(residuals.rphi->cloneFull());
}

double updateNormAndPrint(double norm, const Vector<double>& vec, bool print, std::string&& name) {
   const double infnorm = vec.inf_norm();

   if (print) {
      const double twonorm = vec.two_norm();

      if (0 == PIPS_MPIgetRank())
         std::cout << name << " inf_norm = " << infnorm << " | two_norm = " << twonorm << "\n";
   }

   return std::max(norm, infnorm);
}

void Residuals::evaluate(Problem& problem, Variables& iterate, bool print_resids) {
#ifdef TIMING
   print_resids = true;
#endif
   const int myRank = PIPS_MPIgetRank();

   double norm = 0.0, gap = 0.0;
   primal_objective = 0.0;
   dual_objective = 0.0;

   primal_objective = problem.objective_value(iterate);

   /*** rQ = Qx + g - A^T y - C^T z - gamma + phi ***/
   problem.objective_gradient(iterate, *lagrangian_gradient); // Qx + g

   // contribution calculate x^T (g + Qx) to duality gap */
   gap = lagrangian_gradient->dotProductWith(*iterate.x);

   problem.ATransmult(1.0, *lagrangian_gradient, -1.0, *iterate.y);
   problem.CTransmult(1.0, *lagrangian_gradient, -1.0, *iterate.z);

   iterate.gamma->selectNonZeros(*ixlow);
   iterate.phi->selectNonZeros(*ixupp);
   if (nxlow > 0)
      lagrangian_gradient->axpy(-1.0, *iterate.gamma);
   if (nxupp > 0)
      lagrangian_gradient->axpy(1.0, *iterate.phi);

   norm = updateNormAndPrint(norm, *lagrangian_gradient, print_resids, "rQ");

   /*** rA = Ax - b ***/
   problem.getbA(*rA);
   problem.Amult(-1.0, *rA, 1.0, *iterate.x);

   norm = updateNormAndPrint(norm, *rA, print_resids, "rA");

   // contribution -d^T y to duality gap
   const double ba_y = problem.bA->dotProductWith(*iterate.y);
   gap -= ba_y;
   dual_objective += ba_y;

   /*** rC = Cx - s ***/
   rC->copyFrom(*iterate.s);
   problem.Cmult(-1.0, *rC, 1.0, *iterate.x);

   norm = updateNormAndPrint(norm, *rC, print_resids, "rC");

   /*** rz = z - lambda + pi ***/
   rz->copyFrom(*iterate.z);

   if (mclow > 0)
      rz->axpy(-1.0, *iterate.lambda);
   if (mcupp > 0)
      rz->axpy(1.0, *iterate.pi);

   norm = updateNormAndPrint(norm, *rz, print_resids, "rz");

   if (mclow > 0) {
      /*** rt = s - d - t ***/
      rt->copyFrom(*iterate.s);
      rt->axpy(-1.0, problem.slowerBound());
      rt->selectNonZeros(*iclow);
      rt->axpy(-1.0, *iterate.t);

      // contribution - d^T lambda to duality gap
      const double bl_lambda = problem.bl->dotProductWith(*iterate.lambda);
      gap -= bl_lambda;
      dual_objective += bl_lambda;

      norm = updateNormAndPrint(norm, *rt, print_resids, "rt");
   }

   if (mcupp > 0) {
      /*** ru = s - f + u ***/
      ru->copyFrom(*iterate.s);
      ru->axpy(-1.0, problem.supperBound());
      ru->selectNonZeros(*icupp);
      ru->axpy(1.0, *iterate.u);

      // contribution - f^T pi to duality gap
      const double bu_pi = problem.bu->dotProductWith(*iterate.pi);
      gap += bu_pi;
      dual_objective -= bu_pi;

      norm = updateNormAndPrint(norm, *ru, print_resids, "ru");
   }

   if (nxlow > 0) {
      /*** rv = x - lx - v ***/
      rv->copyFrom(*iterate.x);
      rv->axpy(-1.0, problem.xlowerBound());
      rv->selectNonZeros(*ixlow);
      rv->axpy(-1.0, *iterate.v);

      norm = updateNormAndPrint(norm, *rv, print_resids, "rv");
      // contribution - lx^T gamma to duality gap
      const double blx_gamma = problem.blx->dotProductWith(*iterate.gamma);
      gap -= blx_gamma;
      dual_objective += blx_gamma;

   }

   if (nxupp > 0) {
      /*** rw = x - ux + w ***/
      rw->copyFrom(*iterate.x);
      rw->axpy(-1.0, problem.xupperBound());
      rw->selectNonZeros(*ixupp);
      rw->axpy(1.0, *iterate.w);

      norm = updateNormAndPrint(norm, *rw, print_resids, "rw");

      // contribution + bu^T phi to duality gap
      const double bux_phi = problem.bux->dotProductWith(*iterate.phi);
      gap += bux_phi;
      dual_objective -= bux_phi;
   }

   mDualityGap = gap;
   mResidualNorm = norm;

   if (print_resids && myRank == 0) {
      std::cout << "Norm residuals: " << mResidualNorm << "\tduality gap: " << mDualityGap << "\n";
   }
}

double Residuals::recompute_residual_norm() {
   mResidualNorm = 0.0;

   double componentNorm = 0.0;
   componentNorm = lagrangian_gradient->inf_norm();

   if (componentNorm > mResidualNorm)
      mResidualNorm = componentNorm;

   componentNorm = rA->inf_norm();
   if (componentNorm > mResidualNorm)
      mResidualNorm = componentNorm;

   componentNorm = rC->inf_norm();
   if (componentNorm > mResidualNorm)
      mResidualNorm = componentNorm;

   if (mclow > 0) {
      componentNorm = rt->inf_norm();
      if (componentNorm > mResidualNorm)
         mResidualNorm = componentNorm;
   }

   if (mcupp > 0) {
      componentNorm = ru->inf_norm();
      if (componentNorm > mResidualNorm)
         mResidualNorm = componentNorm;
   }

   componentNorm = rz->inf_norm();
   if (componentNorm > mResidualNorm)
      mResidualNorm = componentNorm;

   if (nxlow > 0) {
      componentNorm = rv->inf_norm();
      if (componentNorm > mResidualNorm)
         mResidualNorm = componentNorm;
   }

   if (nxupp > 0) {
      componentNorm = rw->inf_norm();
      if (componentNorm > mResidualNorm)
         mResidualNorm = componentNorm;
   }
   return mResidualNorm;
}

void Residuals::add_to_complementarity_residual(const Variables& variables, double alpha) {
   if (mclow > 0)
      rlambda->axzpy(1.0, *variables.t, *variables.lambda);
   if (mcupp > 0)
      rpi->axzpy(1.0, *variables.u, *variables.pi);
   if (nxlow > 0)
      rgamma->axzpy(1.0, *variables.v, *variables.gamma);
   if (nxupp > 0)
      rphi->axzpy(1.0, *variables.w, *variables.phi);

   if (alpha != 0.0) {
      if (mclow > 0)
         rlambda->add_constant(alpha, *iclow);
      if (mcupp > 0)
         rpi->add_constant(alpha, *icupp);
      if (nxlow > 0)
         rgamma->add_constant(alpha, *ixlow);
      if (nxupp > 0)
         rphi->add_constant(alpha, *ixupp);
   }
}

void Residuals::set_complementarity_residual(const Variables& variables, double alpha) {
   this->clear_complementarity_residual();
   this->add_to_complementarity_residual(variables, alpha);
}

void Residuals::clear_complementarity_residual() {
   if (mclow > 0)
      rlambda->setToZero();
   if (mcupp > 0)
      rpi->setToZero();
   if (nxlow > 0)
      rgamma->setToZero();
   if (nxupp > 0)
      rphi->setToZero();
}

void Residuals::clear_linear_residuals() {
   lagrangian_gradient->setToZero();
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


int Residuals::valid_non_zero_pattern() {
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

void Residuals::copy(const Residuals& residuals) {
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

   lagrangian_gradient->copyFrom(*residuals.lagrangian_gradient);
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

double Residuals::constraint_violation() {
   return rA->one_norm() + rC->one_norm();
}

double Residuals::optimality_measure(double mu) {
   return this->lagrangian_gradient->inf_norm();// + mu;
}
double Residuals::feasibility_measure(double mu) {
   //double complementarity = 0.;
   // compute componentwise products of complementary variables
//   if (mclow > 0) {
//      SmartPointer<Vector<double> > rt_copy = SmartPointer<Vector<double> >(rt->cloneFull());
//      rt_copy->componentMult(*rlambda);
//      rt_copy->addConstant(-mu);
//      complementarity += rt_copy->onenorm();
//   }
//   if (mcupp > 0) {
//      SmartPointer<Vector<double> > ru_copy = SmartPointer<Vector<double> >(ru->cloneFull());
//      ru_copy->componentMult(*rphi);
//      ru_copy->addConstant(-mu);
//      complementarity += ru_copy->onenorm();
//   }
//   if (nxlow > 0) {
//      SmartPointer<Vector<double> > rv_copy = SmartPointer<Vector<double> >(rv->cloneFull());
//      rv_copy->componentMult(*rgamma);
//      rv_copy->addConstant(-mu);
//      complementarity += rv_copy->onenorm();
//   }
//   if (nxupp > 0) {
//      SmartPointer<Vector<double> > rw_copy = SmartPointer<Vector<double> >(rw->cloneFull());
//      rw_copy->componentMult(*rphi);
//      rw_copy->addConstant(-mu);
//      complementarity += rw_copy->onenorm();
//   }
   return this->constraint_violation();
}