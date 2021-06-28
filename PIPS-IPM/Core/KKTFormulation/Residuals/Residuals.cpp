#include "Residuals.h"
#include "Variables.h"
#include "Problem.hpp"
#include "pipsdef.h"

#include <iostream>
#include <utility>

Residuals::Residuals(std::unique_ptr<Vector<double>> rQ_, std::unique_ptr<Vector<double>> rA_,
   std::unique_ptr<Vector<double>> rC_, std::unique_ptr<Vector<double>> rz_, std::unique_ptr<Vector<double>> rt_,
   std::unique_ptr<Vector<double>> rlambda_, std::unique_ptr<Vector<double>> ru_, std::unique_ptr<Vector<double>> rpi_,
   std::unique_ptr<Vector<double>> rv_, std::unique_ptr<Vector<double>> rgamma_, std::unique_ptr<Vector<double>> rw_,
   std::unique_ptr<Vector<double>> rphi_, std::shared_ptr<Vector<double>> ixlow_,
   std::shared_ptr<Vector<double>> ixupp_, std::shared_ptr<Vector<double>> iclow_,
   std::shared_ptr<Vector<double>> icupp_) :
   ixupp{std::move(ixupp_)},
   nxupp{ixupp->number_nonzeros()},
   ixlow{std::move(ixlow_)}, nxlow{ixlow->number_nonzeros()},
   icupp{std::move(icupp_)}, mcupp{icupp->number_nonzeros()}, iclow{std::move(iclow_)}, mclow{iclow->number_nonzeros()},
   lagrangian_gradient{std::move(rQ_)},
   equality_residuals{std::move(rA_)}, inequality_residuals{std::move(rC_)}, inequality_dual_residuals{std::move(rz_)}, rv{std::move(rv_)},
   rw{std::move(rw_)},
   rt{std::move(rt_)}, ru{std::move(ru_)}, rgamma{std::move(rgamma_)},
   rphi{std::move(rphi_)},
   rlambda{std::move(rlambda_)}, rpi{std::move(rpi_)} {
}

Residuals::Residuals(const Residuals& residuals) : residual_norm{residuals.residual_norm},
   duality_gap{residuals.duality_gap},
   primal_objective{residuals.primal_objective}, dual_objective{residuals.dual_objective}, nx{residuals.nx},
   my{residuals.my}, mz{residuals.mz},
   ixupp{residuals.ixupp}, nxupp{residuals.nxupp}, ixlow{residuals.ixlow}, nxlow{residuals.nxlow},
   icupp{residuals.icupp}, mcupp{residuals.mcupp}, iclow{residuals.iclow}, mclow{residuals.mclow},
   lagrangian_gradient{residuals.lagrangian_gradient->cloneFull()}, equality_residuals{residuals.equality_residuals->cloneFull()},
   inequality_residuals{residuals.inequality_residuals->cloneFull()}, inequality_dual_residuals{residuals.inequality_dual_residuals->cloneFull()}, rv{residuals.rv->cloneFull()},
   rw{residuals.rw->cloneFull()},
   rt{residuals.rt->cloneFull()}, ru{residuals.ru->cloneFull()}, rgamma{residuals.rgamma->cloneFull()},
   rphi{residuals.rphi->cloneFull()},
   rlambda{residuals.rlambda->cloneFull()}, rpi{residuals.rpi->cloneFull()} {}

std::unique_ptr<Residuals> Residuals::cloneFull() const {
   return std::make_unique<Residuals>(*this);
}

double compute_inf_norm(const Vector<double>& vec, bool print, std::string&& name) {
   const double infnorm = vec.inf_norm();

   if (print) {
      const double twonorm = vec.two_norm();

      if (0 == PIPS_MPIgetRank())
         std::cout << name << " inf_norm = " << infnorm << " | two_norm = " << twonorm << "\n";
   }

   return infnorm;
}

void Residuals::evaluate(Problem& problem, Variables& iterate, bool print_residuals) {
#ifdef TIMING
   print_resids = true;
#endif
   const int myRank = PIPS_MPIgetRank();

   this->residual_norm = 0.0;
   this->duality_gap = 0.0;
   this->primal_objective = problem.evaluate_objective(iterate);
   this->dual_objective = 0.0;

   /*** rQ = Qx + g - A^T y - C^T z - gamma + phi ***/
   problem.evaluate_objective_gradient(iterate, *this->lagrangian_gradient); // Qx + g

   // contribution calculate x^T (g + Qx) to duality gap */
   this->duality_gap = this->lagrangian_gradient->dotProductWith(*iterate.primals);

   problem.ATransmult(1.0, *this->lagrangian_gradient, -1.0, *iterate.equality_duals);
   problem.CTransmult(1.0, *this->lagrangian_gradient, -1.0, *iterate.inequality_duals);

   iterate.primal_lower_bound_gap_dual->selectNonZeros(*ixlow);
   iterate.primal_upper_bound_gap_dual->selectNonZeros(*ixupp);
   if (nxlow > 0)
      this->lagrangian_gradient->axpy(-1.0, *iterate.primal_lower_bound_gap_dual);
   if (nxupp > 0)
      this->lagrangian_gradient->axpy(1.0, *iterate.primal_upper_bound_gap_dual);

   this->residual_norm = std::max(this->residual_norm, compute_inf_norm(*this->lagrangian_gradient, print_residuals, "rQ"));

   /*** rA = Ax - b ***/
   problem.getbA(*this->equality_residuals);
   problem.Amult(-1.0, *this->equality_residuals, 1.0, *iterate.primals);
   this->residual_norm = std::max(this->residual_norm, compute_inf_norm(*this->equality_residuals, print_residuals, "rA"));

   // contribution -d^T y to duality gap
   const double ba_y = problem.equality_rhs->dotProductWith(*iterate.equality_duals);
   this->duality_gap -= ba_y;
   this->dual_objective += ba_y;

   /*** rC = Cx - s ***/
   this->inequality_residuals->copyFrom(*iterate.slacks);
   problem.Cmult(-1.0, *this->inequality_residuals, 1.0, *iterate.primals);

   this->residual_norm = std::max(this->residual_norm, compute_inf_norm(*this->inequality_residuals, print_residuals, "rC"));

   /*** rz = z - lambda + pi ***/
   this->inequality_dual_residuals->copyFrom(*iterate.inequality_duals);

   if (mclow > 0)
      this->inequality_dual_residuals->axpy(-1.0, *iterate.slack_lower_bound_gap_dual);
   if (mcupp > 0)
      this->inequality_dual_residuals->axpy(1.0, *iterate.slack_upper_bound_gap_dual);

   this->residual_norm = std::max(this->residual_norm, compute_inf_norm(*this->inequality_dual_residuals, print_residuals, "rz"));

   if (mclow > 0) {
      /*** rt = s - d - t ***/
      this->rt->copyFrom(*iterate.slacks);
      this->rt->axpy(-1.0, problem.s_lower_bound());
      this->rt->selectNonZeros(*iclow);
      this->rt->axpy(-1.0, *iterate.slack_lower_bound_gap);

      // contribution - d^T lambda to duality gap
      const double bl_lambda = problem.inequality_lower_bounds->dotProductWith(*iterate.slack_lower_bound_gap_dual);
      this->duality_gap -= bl_lambda;
      this->dual_objective += bl_lambda;
      this->residual_norm = std::max(this->residual_norm, compute_inf_norm(*this->rt, print_residuals, "rt"));
   }

   if (mcupp > 0) {
      /*** ru = s - f + u ***/
      this->ru->copyFrom(*iterate.slacks);
      this->ru->axpy(-1.0, problem.s_upper_bound());
      this->ru->selectNonZeros(*icupp);
      this->ru->axpy(1.0, *iterate.slack_upper_bound_gap);

      // contribution - f^T pi to duality gap
      const double bu_pi = problem.inequality_upper_bounds->dotProductWith(*iterate.slack_upper_bound_gap_dual);
      this->duality_gap += bu_pi;
      this->dual_objective -= bu_pi;
      this->residual_norm = std::max(this->residual_norm, compute_inf_norm(*this->ru, print_residuals, "ru"));
   }

   if (nxlow > 0) {
      /*** rv = x - lx - v ***/
      this->rv->copyFrom(*iterate.primals);
      this->rv->axpy(-1.0, problem.x_lower_bound());
      this->rv->selectNonZeros(*ixlow);
      this->rv->axpy(-1.0, *iterate.primal_lower_bound_gap);

      this->residual_norm = std::max(this->residual_norm, compute_inf_norm(*this->rv, print_residuals, "rv"));
      // contribution - lx^T gamma to duality gap
      const double blx_gamma = problem.primal_lower_bounds->dotProductWith(*iterate.primal_lower_bound_gap_dual);
      this->duality_gap -= blx_gamma;
      this->dual_objective += blx_gamma;
   }

   if (nxupp > 0) {
      /*** rw = x - ux + w ***/
      this->rw->copyFrom(*iterate.primals);
      this->rw->axpy(-1.0, problem.x_upper_bound());
      this->rw->selectNonZeros(*ixupp);
      this->rw->axpy(1.0, *iterate.primal_upper_bound_gap);

      this->residual_norm = std::max(this->residual_norm, compute_inf_norm(*this->rw, print_residuals, "rw"));

      // contribution + bu^T phi to duality gap
      const double bux_phi = problem.primal_upper_bounds->dotProductWith(*iterate.primal_upper_bound_gap_dual);
      this->duality_gap += bux_phi;
      this->dual_objective -= bux_phi;
   }

   if (print_residuals && myRank == 0) {
      std::cout << "Norm residuals: " << this->residual_norm << "\tduality gap: " << this->duality_gap << "\n";
   }
}

double Residuals::compute_residual_norm() {
   residual_norm = 0.0;

   double componentNorm = 0.0;
   componentNorm = lagrangian_gradient->inf_norm();

   if (componentNorm > residual_norm)
      residual_norm = componentNorm;

   componentNorm = equality_residuals->inf_norm();
   if (componentNorm > residual_norm)
      residual_norm = componentNorm;

   componentNorm = inequality_residuals->inf_norm();
   if (componentNorm > residual_norm)
      residual_norm = componentNorm;

   if (mclow > 0) {
      componentNorm = rt->inf_norm();
      if (componentNorm > residual_norm)
         residual_norm = componentNorm;
   }

   if (mcupp > 0) {
      componentNorm = ru->inf_norm();
      if (componentNorm > residual_norm)
         residual_norm = componentNorm;
   }

   componentNorm = inequality_dual_residuals->inf_norm();
   if (componentNorm > residual_norm)
      residual_norm = componentNorm;

   if (nxlow > 0) {
      componentNorm = rv->inf_norm();
      if (componentNorm > residual_norm)
         residual_norm = componentNorm;
   }

   if (nxupp > 0) {
      componentNorm = rw->inf_norm();
      if (componentNorm > residual_norm)
         residual_norm = componentNorm;
   }
   return residual_norm;
}

void Residuals::add_to_complementarity_residual(const Variables& variables, double alpha) {
   if (mclow > 0)
      rlambda->axzpy(1.0, *variables.slack_lower_bound_gap, *variables.slack_lower_bound_gap_dual);
   if (mcupp > 0)
      rpi->axzpy(1.0, *variables.slack_upper_bound_gap, *variables.slack_upper_bound_gap_dual);
   if (nxlow > 0)
      rgamma->axzpy(1.0, *variables.primal_lower_bound_gap, *variables.primal_lower_bound_gap_dual);
   if (nxupp > 0)
      rphi->axzpy(1.0, *variables.primal_upper_bound_gap, *variables.primal_upper_bound_gap_dual);

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
   equality_residuals->setToZero();
   inequality_residuals->setToZero();
   inequality_dual_residuals->setToZero();
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


int Residuals::valid_non_zero_pattern() const {
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
   residual_norm = residuals.residual_norm;
   duality_gap = residuals.duality_gap;

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
   equality_residuals->copyFrom(*residuals.equality_residuals);
   inequality_residuals->copyFrom(*residuals.inequality_residuals);

   inequality_dual_residuals->copyFrom(*residuals.inequality_dual_residuals);
   rv->copyFrom(*residuals.rv);
   rw->copyFrom(*residuals.rw);
   rt->copyFrom(*residuals.rt);
   ru->copyFrom(*residuals.ru);
   rgamma->copyFrom(*residuals.rgamma);
   rphi->copyFrom(*residuals.rphi);
   rlambda->copyFrom(*residuals.rlambda);
   rpi->copyFrom(*residuals.rpi);
}

void Residuals::write_to_stream(std::ostream& out) const {
   /*
   printf("--------------rQ\n");
   rQ->write_to_stream(out);printf("---------------------------\n");

   printf("rA\n");
   rA->write_to_stream(out);printf("---------------------------\n");
   */
   printf("rC\n");
   inequality_residuals->write_to_stream(out);
   printf("---------------------------\n");


   printf("rz\n");
   inequality_dual_residuals->write_to_stream(out);
   printf("---------------------------\n");
}

double Residuals::constraint_violation() const {
   return equality_residuals->one_norm() + inequality_residuals->one_norm();
}

double Residuals::optimality_measure(/*Problem& problem, Variables& iterate, Scaler& scaler*/) const {
   //return this->lagrangian_gradient->inf_norm();// + mu;
   //return problem.objective_value(iterate);
   return this->primal_objective;
}

double Residuals::feasibility_measure() const {
   return this->constraint_violation();
}

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