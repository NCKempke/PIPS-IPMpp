/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "LinearSystem.h"
#include "Residuals.h"
#include "Problem.hpp"
#include "Variables.h"
#include "Vector.hpp"
#include "DoubleLinearSolver.h"
#include "PIPSIPMppOptions.h"
#include "ProblemFactory.h"
#include <utility>
#include <vector>
#include <functional>
#include <type_traits>
#include <memory>

extern int gOuterBiCGIter;
extern double gOuterBiCGIterAvg;
extern double g_iterNumber;

static std::vector<int> bicgIters;

LinearSystem::LinearSystem(const ProblemFactory& factory_, const Problem& problem, bool create_iter_ref_vecs) : factory(factory_),
      apply_regularization(options::get_bool_parameter("REGULARIZATION")), outerSolve(options::get_int_parameter("OUTER_SOLVE")),
      innerSCSolve(options::get_int_parameter("INNER_SC_SOLVE")),
      outer_solve_refine_original_system{options::get_bool_parameter("OUTER_SOLVE_REFINE_ORIGINAL_SYSTEM")},
      outer_bicg_print_statistics(options::get_bool_parameter("OUTER_BICG_PRINT_STATISTICS")),
      print_linear_system_diagonal_statistics{options::get_bool_parameter("IPM_PRINT_LINEAR_SYSTEM_DIAGONAL_STATISTICS")},
      outer_bicg_eps(options::get_double_parameter("OUTER_BICG_EPSILON")), outer_bicg_max_iter(options::get_int_parameter("OUTER_BICG_MAX_ITER")),
      outer_bicg_max_normr_divergences(options::get_int_parameter("OUTER_BICG_MAX_NORMR_DIVERGENCES")),
      outer_bicg_max_stagnations(options::get_int_parameter("OUTER_BICG_MAX_STAGNATIONS")),
      xyzs_solve_print_residuals(options::get_bool_parameter("XYZS_SOLVE_PRINT_RESISDUAL")), problem(problem) {

   ixlow = problem.primal_lower_bound_indicators;
   ixupp = problem.primal_upper_bound_indicators;
   iclow = problem.inequality_lower_bound_indicators;
   icupp = problem.inequality_upper_bound_indicators;

   nxlow = problem.number_primal_lower_bounds;
   nxupp = problem.number_primal_upper_bounds;
   mclow = problem.number_inequality_lower_bounds;
   mcupp = problem.number_inequality_upper_bounds;

   if (create_iter_ref_vecs) {
      if (outerSolve || xyzs_solve_print_residuals) {
         //for iterative refinement or BICGStab
         sol = factory.make_right_hand_side();
         sol2 = factory.make_right_hand_side();
         res = factory.make_right_hand_side();
         resx = factory.make_primal_vector();
         resy = factory.make_equalities_dual_vector();
         resz = factory.make_inequalities_dual_vector();

         if (outerSolve == 2) {
            //BiCGStab; additional vectors needed
            sol3 = factory.make_right_hand_side();
            res2 = factory.make_right_hand_side();
            res3 = factory.make_right_hand_side();
            res4 = factory.make_right_hand_side();
            res5 = factory.make_right_hand_side();
         }
      }
   }
}

LinearSystem::LinearSystem(const ProblemFactory& factory_, const Problem& problem, std::shared_ptr<Vector<double>> primal_diagonal_,
      std::shared_ptr<Vector<double>> dq_,
      std::shared_ptr<Vector<double>> nomegaInv_, std::shared_ptr<Vector<double>> primal_regularization_,
      std::shared_ptr<Vector<double>> dual_equality_regularization_,
      std::shared_ptr<Vector<double>> dual_inequality_regularization_, std::shared_ptr<Vector<double>> rhs_, bool create_iter_ref_vecs) :
      LinearSystem(factory_, problem,
            create_iter_ref_vecs) {
   primal_diagonal = std::move(primal_diagonal_);
   dq = std::move(dq_);
   nomegaInv = std::move(nomegaInv_);
   primal_regularization_diagonal = std::move(primal_regularization_);
   dual_equality_regularization_diagonal = std::move(dual_equality_regularization_);
   dual_inequality_regularization_diagonal = std::move(dual_inequality_regularization_);
   rhs = std::move(rhs_);
}

LinearSystem::LinearSystem(const ProblemFactory& factory_, const Problem& problem) : LinearSystem(factory_, problem, true) {
   if (nxupp + nxlow > 0) {
      primal_diagonal = factory.make_primal_vector();
      dq = factory.make_primal_vector();
      problem.hessian_diagonal(*dq);
   }

   nomegaInv = factory.make_inequalities_dual_vector();
   rhs = factory.make_right_hand_side();


   primal_regularization_diagonal = factory.make_primal_vector();
   dual_equality_regularization_diagonal = factory.make_equalities_dual_vector();
   dual_inequality_regularization_diagonal = factory.make_inequalities_dual_vector();
}

int LinearSystem::getIntValue(const std::string& s) const {
   if (s == "BICG_NITERATIONS")
      return bicg_niterations;
   else if (s == "BICG_CONV_FLAG")
      return static_cast<std::underlying_type<IterativeSolverSolutionStatus>::type>(bicg_conv_flag);
   else {
      std::cout << "Unknown observer int request in LinearSystem.C: " << s << "\n";
      return -1;
   }
}

bool LinearSystem::getBoolValue(const std::string& s) const {
   if (s == "BICG_CONVERGED")
      return bicg_conv_flag == IterativeSolverSolutionStatus::CONVERGED;
   else if (s == "BICG_SKIPPED")
      return bicg_conv_flag == IterativeSolverSolutionStatus::SKIPPED;
   else if (s == "BICG_DIVERGED")
      return bicg_conv_flag == IterativeSolverSolutionStatus::DIVERGED;
   else if (s == "BICG_BREAKDOWN")
      return bicg_conv_flag == IterativeSolverSolutionStatus::BREAKDOWN;
   else if (s == "BICG_STAGNATION")
      return bicg_conv_flag == IterativeSolverSolutionStatus::STAGNATION;
   else if (s == "BICG_EXCEED_MAX_ITER")
      return bicg_conv_flag == IterativeSolverSolutionStatus::NOT_CONVERGED_MAX_ITERATIONS;
   else {
      std::cout << "Unknown observer bool request in LinearSystem.C: " << s << "\n";
      return false;
   }
}

double LinearSystem::getDoubleValue(const std::string& s) const {
   if (s == "BICG_RESNORM")
      return bicg_resnorm;
   else if (s == "BICG_RELRESNORM")
      return bicg_relresnorm;
   else {
      std::cout << "Unknown observer double request in LinearSystem.C: " << s << "\n";
      return 0.0;
   }
}


static void biCGStabCommunicateStatus(int flag, int it) {
   double iterAvg = 0.0;

   /* IP algorithm started? */
   if (g_iterNumber >= 0.5) {
      bicgIters.push_back(it);

      for (int bicgIter : bicgIters)
         iterAvg += double(bicgIter);

      iterAvg /= bicgIters.size();
   }
   else {
      iterAvg = it;
   }

   gOuterBiCGIterAvg = iterAvg;
   gOuterBiCGIter = it;
}

static bool isZero(double val, LinearSystem::IterativeSolverSolutionStatus& status) {
   if (PIPSisZero(val)) {
      status = LinearSystem::IterativeSolverSolutionStatus::BREAKDOWN;
      return true;
   }

   return false;
}

void LinearSystem::factorize(const Variables& iterate) {

   assert(iterate.valid_non_zero_pattern());
   assert(iterate.valid_non_zero_pattern());

   put_barrier_parameter(iterate.mu());

   computeDiagonals(*iterate.slack_lower_bound_gap, *iterate.slack_lower_bound_gap_dual, *iterate.slack_upper_bound_gap,
         *iterate.slack_upper_bound_gap_dual, *iterate.primal_lower_bound_gap, *iterate.primal_lower_bound_gap_dual,
         *iterate.primal_upper_bound_gap, *iterate.primal_upper_bound_gap_dual);

   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL_TESTING")) {
      std::cout << "Setting diags to 1.0 for Hierarchical debugging\n";
      primal_diagonal->setToConstant(1.0);
      nomegaInv->setToConstant(-1.0);
   }

   reset_regularization();

   if (nxlow + nxupp > 0) {
      put_primal_diagonal();
   }

   if (mclow + mcupp > 0) {
      put_dual_inequalites_diagonal();
   }

   if (print_linear_system_diagonal_statistics)
      print_diagonal_statistics(*iterate.slack_lower_bound_gap, *iterate.slack_lower_bound_gap_dual, *iterate.slack_upper_bound_gap,
         *iterate.slack_upper_bound_gap_dual, *iterate.primal_lower_bound_gap, *iterate.primal_lower_bound_gap_dual,
         *iterate.primal_upper_bound_gap, *iterate.primal_upper_bound_gap_dual);
}

void LinearSystem::reset_regularization() {

   primal_regularization_diagonal->setToZero();
   dual_equality_regularization_diagonal->setToZero();
   dual_inequality_regularization_diagonal->setToZero();

   clear_dual_equality_diagonal();
}

void LinearSystem::print_diagonal_statistics(Vector<double>& t, Vector<double>& lambda, Vector<double>& u, Vector<double>& pi, Vector<double>& v,
   Vector<double>& gamma, Vector<double>& w, Vector<double>& phi) const {
   assert(primal_diagonal);
   assert(nomegaInv);

   const double infnorm_primal_diagonal = primal_diagonal->inf_norm();
   const double twonorm_primal_diagonal = primal_diagonal->two_norm();

   /* find biggest Gamma/V + Phi/W */
   const bool find_max = false;
   const auto [dual_lower_bound_primal, slack_lower_bound_primal, dual_upper_bound_primal, slack_upper_bound_primal] =
      primal_diagonal->find_abs_nonzero_max_min_pair_a_by_b_plus_c_by_d(gamma, v, *ixlow, nxlow > 0, phi, w, *ixupp, nxupp > 0, find_max);

   const double infnorm_inequalities_diagonal = nomegaInv->inf_norm();
   const double twonorm_inequalities_diagonal = nomegaInv->two_norm();

   /*** find smallest Lambda/T + Pi/U (nOmegaInv is the negative inverse of that value) ***/
   const bool find_min = true;
   const auto [dual_lower_bound_inequalities, slack_lower_bound_inequalitites, dual_upper_bound_inequalities, slack_upper_bound_inequalities] =
      nomegaInv->find_abs_nonzero_max_min_pair_a_by_b_plus_c_by_d(lambda, t, *iclow, mclow > 0, pi, u, *icupp, mcupp > 0, find_min);

   if (PIPS_MPIgetRank() == 0) {
      std::cout << "Primal diagonal: ||.||_inf=" << infnorm_primal_diagonal << ", ||.||_2=" << twonorm_primal_diagonal << "\n";
      std::cout << "\tBiggest entry in diagonal Gamma/V + Phi/W :\n";
      std::cout << "\tGamma/V = " << dual_lower_bound_primal << " / " << slack_lower_bound_primal << " = " << dual_lower_bound_primal/slack_lower_bound_primal << "\n";
      std::cout << "\tPhi/W = " << dual_upper_bound_primal << " / " << slack_upper_bound_primal << " = " << dual_upper_bound_primal/slack_upper_bound_primal << "\n\n";
      std::cout << "Inequalities diagonal: ||.||_inf=" << infnorm_inequalities_diagonal << ", ||.||_2=" << twonorm_inequalities_diagonal << "\n";
      std::cout << "\tSmallest entry in diagonal Lambda/T + Pi/U (gets inverted) :\n";
      std::cout << "\tLambda/T = " << dual_lower_bound_inequalities << " / " << slack_lower_bound_inequalitites << " = " << dual_lower_bound_inequalities/slack_lower_bound_inequalitites << "\n";
      std::cout << "\tPi/U = " << dual_upper_bound_inequalities << " / " << slack_upper_bound_inequalities << " = " << dual_upper_bound_inequalities/slack_upper_bound_inequalities << "\n";
      std::cout << "\t(Lambda/T + Pi/U)^-1 = " << 1.0 / (dual_lower_bound_inequalities/slack_lower_bound_inequalitites + dual_upper_bound_inequalities/slack_upper_bound_inequalities) << std::endl;
   }
}

void LinearSystem::print_regularization_statistics() const {
   assert(primal_regularization_diagonal);
   assert(dual_equality_regularization_diagonal);
   assert(dual_inequality_regularization_diagonal);

   const double infnorm_primal_regularization = primal_regularization_diagonal->inf_norm();
   const double infnorm_dual_equality_regularization = dual_equality_regularization_diagonal->inf_norm();
   const double infnorm_dual_inequality_regularization = dual_inequality_regularization_diagonal->inf_norm();

   if (PIPS_MPIgetRank() == 0) {
      std::cout << "Regularized system with (||dP||_inf, ||dDy||_inf, ||dDz||_inf) = (" << infnorm_primal_regularization << ", "
                << infnorm_dual_equality_regularization << ", " << infnorm_dual_inequality_regularization << ")\n";
   }
}

void LinearSystem::computeDiagonals(Vector<double>& t, Vector<double>& lambda, Vector<double>& u, Vector<double>& pi, Vector<double>& v,
      Vector<double>& gamma, Vector<double>& w, Vector<double>& phi) {

   /*** dd = dQ + Gamma/V + Phi/W ***/
   if (nxlow + nxupp > 0) {
      primal_diagonal->copyFrom(*dq);
   }
   if (nxupp + nxlow > 0) {
      if (nxlow > 0) {
         primal_diagonal->add_quotient(1.0, gamma, v, *ixlow);
      }
      if (nxupp > 0)
         primal_diagonal->add_quotient(1.0, phi, w, *ixupp);
   }
   assert(primal_diagonal->all_of([](const double& d) {
      return d >= 0;
   }));

   nomegaInv->setToZero();
   /*** omega = Lambda/T + Pi/U ***/
   if (mclow > 0)
      nomegaInv->add_quotient(1.0, lambda, t, *iclow);
   if (mcupp > 0)
      nomegaInv->add_quotient(1.0, pi, u, *icupp);

   assert(nomegaInv->all_of([](const double& d) {
      return d >= 0;
   }));

   /*** omega = -omega^-1 ***/
   nomegaInv->safe_invert();
   nomegaInv->negate();
}

void LinearSystem::factorize_with_correct_inertia() {
   regularization_strategy->notify_new_step();

   auto[last_primal_regularization, last_dual_equality_regularization, last_dual_inequality_regularization] =
   this->regularization_strategy->get_default_regularization();

   this->add_regularization_local_kkt(last_primal_regularization,
         last_dual_equality_regularization, last_dual_inequality_regularization);

   /* factor once without applying regularization */
   solver->matrixChanged();
   if (!solver->reports_inertia()) {
      return;
   }

   // TODO : add max tries..
   while (!regularization_strategy->is_inertia_correct(solver->get_inertia())) {
      auto[primal_regularization_value, dual_equality_regularization_value, dual_inequality_regularization_value] =
      this->regularization_strategy->get_regularization_parameters(solver->get_inertia(), barrier_parameter_current_iterate);

      assert(primal_regularization_value >= last_primal_regularization);
      assert(dual_equality_regularization_value >= last_dual_equality_regularization);
      assert(dual_inequality_regularization_value >= last_dual_inequality_regularization);

      this->add_regularization_local_kkt(primal_regularization_value - last_primal_regularization,
            dual_equality_regularization_value - last_dual_equality_regularization,
            dual_inequality_regularization_value - last_dual_inequality_regularization);
      solver->matrixChanged();
   }
}

void LinearSystem::solve(const Variables& variables, const Residuals& residuals, Variables& step) {
   assert(variables.valid_non_zero_pattern());
   assert(residuals.valid_non_zero_pattern());

   /*** compute rX ***/
   /* rx = rQ */
   step.primals->copyFrom(*residuals.lagrangian_gradient);
   if (nxlow > 0) {
      Vector<double>& gamma_by_v = *step.primal_lower_bound_gap;
      gamma_by_v.copyFrom(*variables.primal_lower_bound_gap_dual);
      gamma_by_v.divideSome(*variables.primal_lower_bound_gap, *ixlow);

      /* rx = rQ + Gamma/V rv */
      step.primals->add_product(1.0, gamma_by_v, *residuals.rv);
      /* rx = rQ + Gamma/V rv + rGamma/V */
      step.primals->add_quotient(1.0, *residuals.rgamma, *variables.primal_lower_bound_gap, *ixlow);
   }

   if (nxupp > 0) {
      Vector<double>& phi_by_w = *step.primal_upper_bound_gap;
      phi_by_w.copyFrom(*variables.primal_upper_bound_gap_dual);
      phi_by_w.divideSome(*variables.primal_upper_bound_gap, *ixupp);

      /* rx = rQ + Gamma/V * rv + rGamma/V + Phi/W * rw */
      step.primals->add_product(1.0, phi_by_w, *residuals.rw);
      /* rx = rQ + Gamma/V * rv + rGamma/V + Phi/W * rw - rphi/W */
      step.primals->add_quotient(-1.0, *residuals.rphi, *variables.primal_upper_bound_gap, *ixupp);
   }

   // start by partially computing step.s
   /*** compute rs ***/
   /* step->s = rz */
   step.slacks->copyFrom(*residuals.inequality_dual_residuals);
   if (mclow > 0) {
      Vector<double>& lambda_by_t = *step.slack_lower_bound_gap;
      lambda_by_t.copyFrom(*variables.slack_lower_bound_gap_dual);
      lambda_by_t.divideSome(*variables.slack_lower_bound_gap, *iclow);

      /* step.s = rz + Lambda/T * rt */
      step.slacks->add_product(1.0, lambda_by_t, *residuals.rt);
      /* step.s = rz + Lambda/T * rt + rlambda/T */
      step.slacks->add_quotient(1.0, *residuals.rlambda, *variables.slack_lower_bound_gap, *iclow);
   }

   if (mcupp > 0) {
      Vector<double>& pi_by_u = *step.slack_upper_bound_gap;
      pi_by_u.copyFrom(*variables.slack_upper_bound_gap_dual);
      pi_by_u.divideSome(*variables.slack_upper_bound_gap, *icupp);

      /* step.s = rz + Lambda/T * rt + rlambda/T + Pi/U *ru */
      step.slacks->add_product(1.0, pi_by_u, *residuals.ru);
      /* step.s = rz + Lambda/T * rt + rlambda/T + Pi/U *ru - rpi/U */
      step.slacks->add_quotient(-1.0, *residuals.rpi, *variables.slack_upper_bound_gap, *icupp);
   }

   /*** ry = rA ***/
   step.equality_duals->copyFrom(*residuals.equality_residuals);
   /*** rz = rC ***/
   step.inequality_duals->copyFrom(*residuals.inequality_residuals);

   {
      solveXYZS(*step.primals, *step.equality_duals, *step.inequality_duals, *step.slacks);
   }

   if (mclow > 0) {
      /* Dt = Ds - rt */
      step.slack_lower_bound_gap->copyFrom(*step.slacks);
      step.slack_lower_bound_gap->add(-1.0, *residuals.rt);
      step.slack_lower_bound_gap->selectNonZeros(*iclow);

      /* Dlambda = T^-1 (rlambda - Lambda * Dt ) */
      step.slack_lower_bound_gap_dual->copyFrom(*residuals.rlambda);
      step.slack_lower_bound_gap_dual->add_product(-1.0, *variables.slack_lower_bound_gap_dual, *step.slack_lower_bound_gap);
      step.slack_lower_bound_gap_dual->divideSome(*variables.slack_lower_bound_gap, *iclow);
      //!
      step.slack_lower_bound_gap_dual->selectNonZeros(*iclow);
   }

   if (mcupp > 0) {
      /* Du = ru - Ds */
      step.slack_upper_bound_gap->copyFrom(*residuals.ru);
      step.slack_upper_bound_gap->add(-1.0, *step.slacks);
      step.slack_upper_bound_gap->selectNonZeros(*icupp);

      /* Dpi = U^-1 ( rpi - Pi * Du ) */
      step.slack_upper_bound_gap_dual->copyFrom(*residuals.rpi);
      step.slack_upper_bound_gap_dual->add_product(-1.0, *variables.slack_upper_bound_gap_dual, *step.slack_upper_bound_gap);
      step.slack_upper_bound_gap_dual->divideSome(*variables.slack_upper_bound_gap, *icupp);
      //!
      step.slack_upper_bound_gap_dual->selectNonZeros(*icupp);
   }

   if (nxlow > 0) {
      /* Dv = Dx - rv */
      step.primal_lower_bound_gap->copyFrom(*step.primals);
      step.primal_lower_bound_gap->add(-1.0, *residuals.rv);
      step.primal_lower_bound_gap->selectNonZeros(*ixlow);

      /* Dgamma = V^-1 ( rgamma - Gamma * Dv ) */
      step.primal_lower_bound_gap_dual->copyFrom(*residuals.rgamma);
      step.primal_lower_bound_gap_dual->add_product(-1.0, *variables.primal_lower_bound_gap_dual, *step.primal_lower_bound_gap);
      step.primal_lower_bound_gap_dual->divideSome(*variables.primal_lower_bound_gap, *ixlow);
      //!
      step.primal_lower_bound_gap_dual->selectNonZeros(*ixlow);
   }

   if (nxupp > 0) {
      /* Dw = rw - Dx */
      step.primal_upper_bound_gap->copyFrom(*residuals.rw);
      step.primal_upper_bound_gap->add(-1.0, *step.primals);
      step.primal_upper_bound_gap->selectNonZeros(*ixupp);

      /* Dphi = W^-1 ( rphi - Phi * Dw ) */
      step.primal_upper_bound_gap_dual->copyFrom(*residuals.rphi);
      step.primal_upper_bound_gap_dual->add_product(-1.0, *variables.primal_upper_bound_gap_dual, *step.primal_upper_bound_gap);
      step.primal_upper_bound_gap_dual->divideSome(*variables.primal_upper_bound_gap, *ixupp);
      //!
      step.primal_upper_bound_gap_dual->selectNonZeros(*ixupp);
   }
   assert(step.valid_non_zero_pattern());
}

void LinearSystem::solveXYZS(Vector<double>& stepx, Vector<double>& stepy, Vector<double>& stepz, Vector<double>& steps) {
   /* step->z = rC */
   /* step->s = rz + Lambda/T * rt + rlambda/T + Pi/U *ru - rpi/U */

   /* rx = rQ + Gamma/V * rv + rGamma/V + Phi/W * rw - rphi/W */
   /* ry = rA */
   /* rz = rC + Omega^-1 ( rz + Lambda/T * rt + rlambda/T + Pi/U *ru - rpi/U ) */
   stepz.add_product(-1.0, *nomegaInv, steps);

   std::unique_ptr<Vector<double>> residual;
   if (xyzs_solve_print_residuals) {
      residual.reset(rhs->clone_full());
      joinRHS(*residual, stepx, stepy, stepz);

      const double xinf = stepx.inf_norm();
      const double yinf = stepy.inf_norm();
      const double zinf = stepz.inf_norm();

      if (PIPS_MPIgetRank() == 0)
         std::cout << "rhsx norm : " << xinf << ",\trhsy norm : " << yinf << ",\trhsz norm : " << zinf << "\n";
   }

   assert(rhs);
   LinearSystem::joinRHS(*rhs, stepx, stepy, stepz);

   if (outerSolve == 1) {
      ///////////////////////////////////////////////////////////////
      // Iterative refinement
      ///////////////////////////////////////////////////////////////
      if (outer_solve_refine_original_system) {
         auto computeResiduals = [this, &stepx, &stepy, &stepz, &capture0 = problem](auto&& PH1, auto&& PH2) {
            compute_system_residuals(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2), stepx, stepy, stepz, capture0);
         };

         solveCompressedIterRefin(computeResiduals);
      }
      else {
         auto computeResiduals = [this, &stepx, &stepy, &stepz, &capture0 = problem](auto&& PH1, auto&& PH2) {
            compute_regularized_system_residuals(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2), stepx, stepy, stepz, capture0);
         };

         solveCompressedIterRefin(computeResiduals);
      }
   }
   else if (outerSolve == 0) {
      ///////////////////////////////////////////////////////////////
      // Default solve - Schur complement based decomposition
      ///////////////////////////////////////////////////////////////
      solveCompressed(*rhs);
   }
   else {
      assert(outerSolve == 2);
      ///////////////////////////////////////////////////////////////
      // BiCGStab
      ///////////////////////////////////////////////////////////////

      const bool use_regularized_system = !outer_solve_refine_original_system;
      auto matMult = [this, &capture0 = problem, &stepx, &stepy, &stepz, use_regularized_system](auto&& PH1, auto&& PH2, auto&& PH3, auto&& PH4) {
         system_mult(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2), std::forward<decltype(PH3)>(PH3),
               std::forward<decltype(PH4)>(PH4), capture0, stepx, stepy, stepz, use_regularized_system);
      };

      auto matInfnorm = [this, &capture0 = problem, &stepx, &stepy, &stepz, use_regularized_system] {
         return matXYZinfnorm(capture0, stepx, stepy, stepz, use_regularized_system);
      };

      solveCompressedBiCGStab(matMult, matInfnorm);
      /* notify observers about result of BiCGStab */
      notifyObservers();
   }

   LinearSystem::separateVars(stepx, stepy, stepz, *sol);

   if (xyzs_solve_print_residuals) {
      assert(sol);
      const double bnorm = residual->inf_norm();
      LinearSystem::joinRHS(*sol, stepx, stepy, stepz);
      this->system_mult(1.0, *residual, -1.0, *sol, problem, stepx, stepy, stepz, true);

      LinearSystem::separateVars(*resx, *resy, *resz, *residual);
      const double resxnorm = resx->inf_norm();
      const double resynorm = resy->inf_norm();
      const double resznorm = resz->inf_norm();

      if (PIPS_MPIgetRank() == 0) {
         std::cout << "bnorm " << bnorm << "\n";
         std::cout << "resx norm: " << resxnorm << "\tnorm/bnorm " << resxnorm / bnorm << "\n";
         std::cout << "resy norm: " << resynorm << "\tnorm/bnorm " << resynorm / bnorm << "\n";
         std::cout << "resz norm: " << resznorm << "\tnorm/bnorm " << resznorm / bnorm << "\n";
      }
   }

   stepy.negate();
   stepz.negate();

   /* Ds = Omega^-1 (rz + Lambda/T * rt + rlambda/T + Pi/U *ru - rpi/U - Dz ) */
   steps.add(-1.0, stepz);
   steps.componentMult(*nomegaInv);
   steps.negate();
}

void check_stagnation(int& n_stagnations, const double& step_length, const double& step_norm, const double& convergence_tolerance, const double& iterate_norm) {
   if (std::fabs(step_length) * step_norm <= convergence_tolerance * iterate_norm)
     ++n_stagnations;
   else
      n_stagnations = 0;
}

void LinearSystem::solveCompressedBiCGStab(const std::function<void(double, Vector<double>&, double, Vector<double>&)>& matMult,
      const std::function<double()>& matInfnorm) {
   //aliases
   Vector<double>& r0 = *res2, & dx = *sol2, & best_x = *sol3, & v = *res3, & t = *res4, & p = *res5;
   Vector<double>& x = *sol, & r = *res, & b = *rhs;

   const double tol = options::get_double_parameter("OUTER_BICG_TOL");
   const double right_hand_side_two_nord = b.two_norm();
   const double convergence_tolerance_scaled_by_rhs_norm = std::max(right_hand_side_two_nord * tol, outer_bicg_eps);

   gOuterBiCGIter = 0;
   bicg_niterations = 0;

   assert(right_hand_side_two_nord >= 0);

   const int myRank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   //starting guess/point
   x.copyFrom(b);

   //solution to the approx. system
   solveCompressed(x);

   //initial residual: res = b - Ax
   r.copyFrom(b);
   matMult(1.0, r, -1.0, x);

   double residual_two_norm = r.two_norm(), min_residual_two_norm = residual_two_norm;
   best_x.copyFrom(x);

   bicg_resnorm = residual_two_norm;
   bicg_relresnorm = bicg_resnorm / right_hand_side_two_nord;

   //quick return if solve is accurate enough
   if (residual_two_norm <= convergence_tolerance_scaled_by_rhs_norm) {

      bicg_conv_flag = IterativeSolverSolutionStatus::SKIPPED;
      biCGStabCommunicateStatus(static_cast<std::underlying_type<IterativeSolverSolutionStatus>::type>(bicg_conv_flag), bicg_niterations);

      if (myRank == 0)
         std::cout << "BiCGStab (it=" << bicg_niterations << ", rel.res.norm=" << bicg_relresnorm << ", rel.r.norm=" << residual_two_norm / right_hand_side_two_nord << ", avg.iter="
                   << gOuterBiCGIterAvg << ") " << bicg_conv_flag << "\n";
      return;
   }

   if (outer_bicg_print_statistics) {
      const double infb = b.inf_norm();
      const double glbinfnorm = matInfnorm();
      const double xonenorm = x.one_norm();

      if (myRank == 0) {
         std::cout << "global system inf_norm=" << glbinfnorm << " x 1norm=" << xonenorm << " convergence_tolerance_scaled_by_rhs_norm/tolnew: " << convergence_tolerance_scaled_by_rhs_norm << " "
                   << (tol * xonenorm * glbinfnorm) << "\n";
         std::cout << "outer BiCGStab starts: " << residual_two_norm << " > " << convergence_tolerance_scaled_by_rhs_norm << " normb2=" << right_hand_side_two_nord << " normbinf=" << infb << " (tolerance=" << tol << ")"
                   << "\n";
      }
   }

   r0.copyFrom(r);

   //normalize
   r0.scale(1 / residual_two_norm);

   bicg_conv_flag = IterativeSolverSolutionStatus::NOT_CONVERGED_MAX_ITERATIONS;
   int normrNDiv = 0;
   int n_half_iterate_stagnations = 0;
   double rho = 1., omega = 1., alpha = 1.;

   //main loop
   for (bicg_niterations = 0; bicg_niterations < outer_bicg_max_iter; bicg_niterations++) {

      //first half of the iterate
      {
         const double rho_last_iter = rho;
         rho = r0.dotProductWith(r);

         if (isZero(rho, bicg_conv_flag))
            break;

         if (bicg_niterations == 0)
            p.copyFrom(r);
         else {
            assert(rho_last_iter != 0.0 && omega != 0.0);
            const double beta = (rho / rho_last_iter) * (alpha / omega);

            if (isZero(beta, bicg_conv_flag))
               break;

            //-------- p = r_i-1 + beta * (p_i-1 - omega_i-1 * v_i-1) --------
            p.add(-omega, v);
            p.scale(beta);
            p.add(1.0, r);
         }

         //precond: dx = \tilde{K}^{-1} p
         dx.copyFrom(p);
         solveCompressed(dx);

         //mat-first: v_i = K * dx
         matMult(0.0, v, 1.0, dx);

         const double rtv = r0.dotProductWith(v);

         if (isZero(rtv, bicg_conv_flag))
            break;

         /* alpha is half step for x */
         alpha = rho / rtv;

         check_stagnation(n_half_iterate_stagnations, alpha, dx.two_norm(), outer_bicg_eps, x.two_norm());

         // x_i/2 = x + alpha * dx ( x = x + alpha * ph)
         x.add(alpha, dx); // half way iterate
         // r = r - alpha * v ( s = r - alpha * v = r - alpha * A * y)
         r.add(-alpha, v); // half way (expected) residual

         // check for convergence
         residual_two_norm = r.two_norm();

         if (residual_two_norm <= convergence_tolerance_scaled_by_rhs_norm) {

            //compute the actual residual using dx as buffer
            dx.copyFrom(b);
            matMult(1.0, dx, -1.0, x);

            bicg_resnorm = dx.two_norm();
            if (bicg_resnorm <= convergence_tolerance_scaled_by_rhs_norm) {
               bicg_conv_flag = IterativeSolverSolutionStatus::CONVERGED;
               break;
            }
         }

         /* for rollback to half-iterate */
         if (residual_two_norm < min_residual_two_norm) {
            // update best for rollback
            min_residual_two_norm = residual_two_norm;
            best_x.copyFrom(x);
         }
      }

      //second half of the iterate now
      {
         // precondition half_step residual
         dx.copyFrom(r);
         solveCompressed(dx);

         //mat-first
         matMult(0.0, t, 1.0, dx);

         const double tt = t.dotProductSelf(1.0);

         if (isZero(tt, bicg_conv_flag))
            break;

         // omega_i = <K_1^-1 t | K1^-1 r> / <K_1^-1 t, K_1^-1 t> (for us K_1 = 1 -> else this is 4 solves with the PIPS system..)
         omega = t.dotProductWith(r) / tt;

         check_stagnation(n_half_iterate_stagnations, omega, dx.two_norm(), outer_bicg_eps, x.two_norm());

         // x=x+omega*dx  ( x =x+omega*sh)
         x.add(omega, dx);
         // r = r-omega*t (r=s-omega*sh)
         r.add(-omega, t);
         //check for convergence
         residual_two_norm = r.two_norm();

         if (residual_two_norm <= convergence_tolerance_scaled_by_rhs_norm || n_half_iterate_stagnations >= outer_bicg_max_stagnations) {
            //compute the actual residual
            dx.copyFrom(b); // dx as buffer
            matMult(1.0, dx, -1.0, x);

            bicg_resnorm = dx.two_norm();

            if (bicg_resnorm <= convergence_tolerance_scaled_by_rhs_norm) {
               //converged
               bicg_conv_flag = IterativeSolverSolutionStatus::CONVERGED;
               break;
            }
         } else {
            if (residual_two_norm >= min_residual_two_norm)
               normrNDiv++;
            else
               normrNDiv = 0;

            if (normrNDiv > outer_bicg_max_normr_divergences) {
               // rollback to best iterate
               x.copyFrom(best_x);
               residual_two_norm = min_residual_two_norm;

               //compute the actual residual
               dx.copyFrom(b); // dx as buffer
               matMult(1.0, dx, -1.0, x);

               bicg_resnorm = dx.two_norm();
               bicg_conv_flag = IterativeSolverSolutionStatus::DIVERGED;
               break;
            }
         } //~end of convergence test

         if (residual_two_norm < min_residual_two_norm) {
            // update best for rollback
            min_residual_two_norm = residual_two_norm;
            best_x.copyFrom(x);
         }
      } //~end of scoping

      if (n_half_iterate_stagnations >= outer_bicg_max_stagnations) {
         // rollback to best iterate
         if (min_residual_two_norm < residual_two_norm) {
            residual_two_norm = min_residual_two_norm;
            x.copyFrom(best_x);
         }

         //compute the actual residual
         dx.copyFrom(b); // dx as buffer
         matMult(1.0, dx, -1.0, x);

         bicg_resnorm = dx.two_norm();

         bicg_conv_flag = IterativeSolverSolutionStatus::STAGNATION;
         break;
      }

      if (isZero(omega, bicg_conv_flag))
         break;

   } //~ end of BiCGStab loop

   bicg_relresnorm = bicg_resnorm / right_hand_side_two_nord;

   biCGStabCommunicateStatus(static_cast<std::underlying_type<IterativeSolverSolutionStatus>::type>(bicg_conv_flag), std::max(bicg_niterations, 1));
   if (myRank == 0) {
      std::cout << "BiCGStab (it=" << std::max(1, bicg_niterations) << ", relative residual norm=" << bicg_relresnorm
         << ", predicted relative residual norm=" << residual_two_norm / right_hand_side_two_nord
         << ", avg.iter=" << gOuterBiCGIterAvg << ") " << bicg_conv_flag << "\n";
   }
}

/**
 * res = beta * res + alpha * mat * sol
 *       [ Q + dq + gamma/ v + phi/w + regP                   AT                                 CT               ]
 * mat = [            A                      dual_equality_regularization_diagonal               0                ]
 *       [            C                                       0                     -(lambda/V + pi/u)^-1 + regDz ]
 * stepx, stepy, stepz are used as temporary buffers
 * if use_regularized_system == false primal_regularization_diagonal, dual_equality_regularization_diagonal and dual_inequality_regularization_diagonal are not used
 */
void
LinearSystem::system_mult(double beta, Vector<double>& res, double alpha, const Vector<double>& sol, const Problem& problem, Vector<double>& solx,
      Vector<double>& soly, Vector<double>& solz, bool use_regularized_system) {
   assert(resx);
   assert(resy);
   assert(resz);
   assert(nomegaInv);
   assert(primal_diagonal);
   assert((primal_regularization_diagonal && dual_equality_regularization_diagonal && dual_inequality_regularization_diagonal) ||
          (!primal_regularization_diagonal && !dual_equality_regularization_diagonal && !dual_inequality_regularization_diagonal));

   separateVars(solx, soly, solz, sol);
   separateVars(*resx, *resy, *resz, res);

   /* resx = beta resx + alpha Q solx + alpha dd solx + alpha primal_regularization_diagonal solx */
   problem.hessian_multiplication(beta, *resx, alpha, solx);
   resx->add_product(alpha, *primal_diagonal, solx);
   if (use_regularized_system && primal_regularization_diagonal)
      resx->add_product(alpha, *primal_regularization_diagonal, solx);

   /* resx = beta resx + alpha Q solx + alpha dd solx + alpha AT soly + alpha CT solz */
   problem.ATransmult(1.0, *resx, alpha, soly);
   problem.CTransmult(1.0, *resx, alpha, solz);

   /* resy = beta resy + alpha A solx */
   problem.Amult(beta, *resy, alpha, solx);
   if (use_regularized_system && dual_equality_regularization_diagonal)
      resy->add_product(alpha, *dual_equality_regularization_diagonal, soly);

   /* resz = beta resz + alpha C solx + alpha nomegaInv solz */
   problem.Cmult(beta, *resz, alpha, solx);
   resz->add_product(alpha, *nomegaInv, solz);
   if (use_regularized_system && dual_inequality_regularization_diagonal)
      resz->add_product(alpha, *dual_inequality_regularization_diagonal, solz);

   LinearSystem::joinRHS(res, *resx, *resy, *resz);
}

/* computes infinity norm of entire system; solx, soly, solz are used as temporary buffers */
double
LinearSystem::matXYZinfnorm(const Problem& problem, Vector<double>& solx, Vector<double>& soly, Vector<double>& solz, bool use_regularized_system) {
   solx.copyFromAbs(*primal_diagonal);
   if (use_regularized_system && primal_regularization_diagonal)
      solx.add(1.0, *primal_regularization_diagonal);

   problem.equality_jacobian->addColSums(solx);
   problem.inequality_jacobian->addColSums(solx);
   double infnorm = solx.inf_norm();

   soly.setToZero();
   if (use_regularized_system && dual_equality_regularization_diagonal)
      soly.add(1.0, *dual_equality_regularization_diagonal);
   soly.negate();

   problem.equality_jacobian->addRowSums(soly);
   infnorm = std::max(infnorm, soly.inf_norm());

   solz.copyFromAbs(*nomegaInv);
   if (use_regularized_system && dual_inequality_regularization_diagonal) {
      solz.add(1.0, *dual_inequality_regularization_diagonal);
   }
   solz.negate();

   problem.inequality_jacobian->addRowSums(solz);
   infnorm = std::max(infnorm, solz.inf_norm());

   return infnorm;
}

void LinearSystem::solveCompressedIterRefin(const std::function<void(Vector<double>& solution, Vector<double>& residual)>& computeResidual) {
   assert(this->rhs);
   assert(this->res);
   assert(this->sol);
   assert(this->sol2);

   /* aliases */
   Vector<double>& residual = *this->res;
   Vector<double>& solution = *this->sol;
   Vector<double>& best_solution = *this->sol2;
   const Vector<double>& right_hand_side = *this->rhs;

   residual.copyFrom(right_hand_side);
   solution.setToZero();

   const double rhs_norm = std::max(1.0, residual.two_norm());
   const double tolerance_iterative_refinement = 1e-8;
   const int max_iterative_refinement_steps = 50;
   const int max_steps_stalling = 5;
   double residual_norm_best_solution = std::numeric_limits<double>::max();
   double residual_norm_current_solution{};
   double residual_last_iterate = std::numeric_limits<double>::max();

   const double min_improvement_factor = 1e-1;
   const double max_divergence_factor = 1e4;
   int n_steps_no_improvement = 0;

   IterativeSolverSolutionStatus convergence_status = IterativeSolverSolutionStatus::DID_NOT_RUN;

   int current_iterative_refinement_step = -1;

   do {
      current_iterative_refinement_step++;

      /* solve A'x = b */
      this->solveCompressed(residual);

      /* x = x + dx */
      solution.add(1.0, residual);

      residual.copyFrom(right_hand_side);

      computeResidual(solution, residual);
      residual_norm_current_solution = residual.two_norm();

      if (residual_norm_current_solution > min_improvement_factor * residual_last_iterate) {
         ++n_steps_no_improvement;
      }
      else {
         n_steps_no_improvement = 0;
      }

      if (residual_norm_current_solution < tolerance_iterative_refinement * rhs_norm) {
         convergence_status = IterativeSolverSolutionStatus::CONVERGED;
         break;
      }
      else if (current_iterative_refinement_step == max_iterative_refinement_steps) {
         convergence_status = IterativeSolverSolutionStatus::NOT_CONVERGED_MAX_ITERATIONS;
         break;
      }
      else if (max_divergence_factor * residual_norm_best_solution < residual_norm_current_solution) {
         convergence_status = IterativeSolverSolutionStatus::DIVERGED;
         break;
      }
      else if (n_steps_no_improvement > max_steps_stalling) {
         convergence_status = IterativeSolverSolutionStatus::STAGNATION;
         break;
      }

      /* store best solution if we compute one more in the next round */
      if (residual_norm_current_solution < residual_norm_best_solution) {
         residual_norm_best_solution = residual_norm_current_solution;
         best_solution.copyFrom(solution);
      }
      residual_last_iterate = residual_norm_current_solution;

   } while (true);

   assert(convergence_status != IterativeSolverSolutionStatus::DID_NOT_RUN);

   /* backtrack to best solution if it was not the last one */
   if (residual_norm_best_solution < residual_norm_current_solution)
      solution.copyFrom(best_solution);

   if (PIPS_MPIgetRank() == 0) {
      const double norm_final = std::min(residual_norm_best_solution, residual_norm_current_solution);
      std::cout << "IterRef (it=" << current_iterative_refinement_step << ", residual norm=" << norm_final << ", relative residual norm="
                << norm_final / rhs_norm << ") " << convergence_status << "\n";
   }
}

/**
 * res = res - mat * sol
 *
 * [ resx ]   [ resx ]   [ Q + dd + regularization           AT                    CT               ] [ solx ]
 * [ resy ] = [ resy ] - [       A                     regularization              0                ] [ soly ]
 * [ resz ]   [ resz ]   [       C                            0          nOmegaInv + regularization ] [ solz ]
 *
 * stepx, stepy, stepz are used as temporary buffers
 */
void LinearSystem::compute_regularized_system_residuals(const Vector<double>& sol, Vector<double>& res, Vector<double>& solx, Vector<double>& soly,
      Vector<double>& solz, const Problem& problem) {
   const double alpha = -1.0;
   const double beta = 1.0;

   system_mult(beta, res, alpha, sol, problem, solx, soly, solz, true);
}

void LinearSystem::compute_system_residuals(const Vector<double>& sol, Vector<double>& res, Vector<double>& solx, Vector<double>& soly,
      Vector<double>& solz, const Problem& problem) {
   const double alpha = -1.0;
   const double beta = 1.0;

   system_mult(beta, res, alpha, sol, problem, solx, soly, solz, false);
}

void
LinearSystem::joinRHS(Vector<double>& rhs_in, const Vector<double>& rhs1_in, const Vector<double>& rhs2_in, const Vector<double>& rhs3_in) {
   rhs_in.jointCopyFrom(rhs1_in, rhs2_in, rhs3_in);
}

void LinearSystem::separateVars(Vector<double>& x_in, Vector<double>& y_in, Vector<double>& z_in, const Vector<double>& vars_in) {
   vars_in.jointCopyTo(x_in, y_in, z_in);
}
