#include <iostream>
#include <utility>
#include <DenseVector.hpp>
#include "Variables.h"
#include "Vector.hpp"
#include "Problem.hpp"
#include "MpsReader.h"

Variables::Variables(std::unique_ptr<Vector<double>> x_in, std::unique_ptr<Vector<double>> s_in, std::unique_ptr<Vector<double>> y_in, std::unique_ptr<Vector<double>> z_in, std::unique_ptr<Vector<double>> v_in,
   std::unique_ptr<Vector<double>> gamma_in, std::unique_ptr<Vector<double>> w_in, std::unique_ptr<Vector<double>> phi_in, std::unique_ptr<Vector<double>> t_in, std::unique_ptr<Vector<double>> lambda_in, std::unique_ptr<Vector<double>> u_in,
   std::unique_ptr<Vector<double>> pi_in, std::shared_ptr<Vector<double>> ixlow_in, std::shared_ptr<Vector<double>> ixupp_in, std::shared_ptr<Vector<double>> iclow_in,
   std::shared_ptr<Vector<double>> icupp_in) :
      primal_lower_bound_indicators{std::move(ixlow_in)}, primal_upper_bound_indicators{std::move(ixupp_in)}, inequality_lower_bound_indicators{std::move(iclow_in)}, inequality_upper_bound_indicators{std::move(icupp_in)}, primals{std::move(x_in)}, slacks{std::move(s_in)},
      equality_duals{std::move(y_in)}, inequality_duals{std::move(z_in)}, primal_lower_bound_gap{std::move(v_in)}, primal_lower_bound_gap_dual{std::move(gamma_in)},
      primal_upper_bound_gap{std::move(w_in)}, primal_upper_bound_gap_dual{std::move(phi_in)}, slack_lower_bound_gap{std::move(t_in)},
      slack_lower_bound_gap_dual{std::move(lambda_in)}, slack_upper_bound_gap{std::move(u_in)}, slack_upper_bound_gap_dual{std::move(pi_in)}{

   nx = primals->length();
   my = equality_duals->length();
   mz = inequality_duals->length();

   assert(nx == primal_lower_bound_indicators->length() || 0 == primal_lower_bound_indicators->length());
   assert(nx == primal_lower_bound_indicators->length() || 0 == primal_lower_bound_indicators->length());
   assert(mz == inequality_lower_bound_indicators->length() || 0 == inequality_lower_bound_indicators->length());
   assert(mz == inequality_upper_bound_indicators->length() || 0 == inequality_upper_bound_indicators->length());

   nxlow = primal_lower_bound_indicators->number_nonzeros();
   nxupp = primal_upper_bound_indicators->number_nonzeros();
   mclow = inequality_lower_bound_indicators->number_nonzeros();
   mcupp = inequality_upper_bound_indicators->number_nonzeros();
   number_complementarity_pairs = mclow + mcupp + nxlow + nxupp;

   assert(mz == slacks->length());
   assert(nx == primal_lower_bound_gap->length() || (0 == primal_lower_bound_gap->length() && nxlow == 0));
   assert(nx == primal_lower_bound_gap_dual->length() || (0 == primal_lower_bound_gap_dual->length() && nxlow == 0));

   assert(nx == primal_upper_bound_gap->length() || (0 == primal_upper_bound_gap->length() && nxupp == 0));
   assert(nx == primal_upper_bound_gap_dual->length() || (0 == primal_upper_bound_gap_dual->length() && nxupp == 0));

   assert(mz == slack_lower_bound_gap->length() || (0 == slack_lower_bound_gap->length() && mclow == 0));
   assert(mz == slack_lower_bound_gap_dual->length() || (0 == slack_lower_bound_gap_dual->length() && mclow == 0));

   assert(mz == slack_upper_bound_gap->length() || (0 == slack_upper_bound_gap->length() && mcupp == 0));
   assert(mz == slack_upper_bound_gap_dual->length() || (0 == slack_upper_bound_gap_dual->length() && mcupp == 0));
}

Variables::Variables(const Variables& other)  :
   number_complementarity_pairs{other.number_complementarity_pairs}, nx{other.nx}, nxupp{other.nxupp}, nxlow{other.nxlow}, my{other.my}, mz{other.mz}, mcupp{other.mcupp}, mclow{other.mclow},
   primal_lower_bound_indicators{other.primal_lower_bound_indicators}, primal_upper_bound_indicators{other.primal_upper_bound_indicators}, inequality_lower_bound_indicators{other.inequality_lower_bound_indicators}, inequality_upper_bound_indicators{other.inequality_upper_bound_indicators},
   primals{other.primals->clone_full()}, slacks{other.slacks->clone_full()}, equality_duals{other.equality_duals->clone_full()}, inequality_duals{other.inequality_duals->clone_full()},
   primal_lower_bound_gap{other.primal_lower_bound_gap->clone_full()}, primal_lower_bound_gap_dual{other.primal_lower_bound_gap_dual->clone_full()},
   primal_upper_bound_gap{other.primal_upper_bound_gap->clone_full()}, primal_upper_bound_gap_dual{other.primal_upper_bound_gap_dual->clone_full()},
   slack_lower_bound_gap{other.slack_lower_bound_gap->clone_full()}, slack_lower_bound_gap_dual{other.slack_lower_bound_gap_dual->clone_full()},
   slack_upper_bound_gap{other.slack_upper_bound_gap->clone_full()}, slack_upper_bound_gap_dual{other.slack_upper_bound_gap_dual->clone_full()}
   {}

std::unique_ptr<Variables> Variables::clone_full() const {
   return std::make_unique<Variables>(*this);
}

double Variables::get_average_distance_to_bound_for_converged_vars(const Problem&, double tol) const {
   assert(0 < tol);

   double sum_small_distance = 0.0;
   int n_close = 0;
   primal_lower_bound_gap->getSumCountIfSmall(tol, sum_small_distance, n_close, &*primal_lower_bound_indicators);
   primal_upper_bound_gap->getSumCountIfSmall(tol, sum_small_distance, n_close, &*primal_upper_bound_indicators);
   slack_upper_bound_gap->getSumCountIfSmall(tol, sum_small_distance, n_close, &*inequality_upper_bound_indicators);
   slack_lower_bound_gap->getSumCountIfSmall(tol, sum_small_distance, n_close, &*inequality_lower_bound_indicators);

   if (n_close == 0)
      return std::numeric_limits<double>::infinity();
   else
      return sum_small_distance / (double) n_close;
}

void Variables::push_slacks_from_bound(double tol, double amount) {
   if (nxlow > 0)
      primal_lower_bound_gap->pushAwayFromZero(tol, amount, &*primal_lower_bound_indicators);
   if (nxupp > 0)
      primal_upper_bound_gap->pushAwayFromZero(tol, amount, &*primal_upper_bound_indicators);
   if (mclow > 0)
      slack_lower_bound_gap->pushAwayFromZero(tol, amount, &*inequality_lower_bound_indicators);
   if (mcupp > 0)
      slack_upper_bound_gap->pushAwayFromZero(tol, amount, &*inequality_upper_bound_indicators);
}

double Variables::mu() const {
   double mu = 0.;
   if (number_complementarity_pairs == 0) {
      return 0.;
   }
   else {

      if (mclow > 0)
         mu += slack_lower_bound_gap->dotProductWith(*slack_lower_bound_gap_dual);
      if (mcupp > 0)
         mu += slack_upper_bound_gap->dotProductWith(*slack_upper_bound_gap_dual);
      if (nxlow > 0)
         mu += primal_lower_bound_gap->dotProductWith(*primal_lower_bound_gap_dual);
      if (nxupp > 0)
         mu += primal_upper_bound_gap->dotProductWith(*primal_upper_bound_gap_dual);

      mu /= (double) number_complementarity_pairs;
      return mu;
   }
}

double Variables::mustep_pd(const Variables& iterate, double alpha_primal, double alpha_dual) const {
   double mu = 0.;
   if (number_complementarity_pairs == 0) {
      return 0.;
   }
   else {
      if (mclow > 0) {
         mu += slack_lower_bound_gap->shiftedDotProductWith(alpha_primal, *iterate.slack_lower_bound_gap, *slack_lower_bound_gap_dual, alpha_dual, *iterate.slack_lower_bound_gap_dual);
      }
      if (mcupp > 0) {
         mu += slack_upper_bound_gap->shiftedDotProductWith(alpha_primal, *iterate.slack_upper_bound_gap, *slack_upper_bound_gap_dual, alpha_dual, *iterate.slack_upper_bound_gap_dual);
      }
      if (nxlow > 0) {
         mu += primal_lower_bound_gap->shiftedDotProductWith(alpha_primal, *iterate.primal_lower_bound_gap, *primal_lower_bound_gap_dual, alpha_dual, *iterate.primal_lower_bound_gap_dual);
      }
      if (nxupp > 0) {
         mu += primal_upper_bound_gap->shiftedDotProductWith(alpha_primal, *iterate.primal_upper_bound_gap, *primal_upper_bound_gap_dual, alpha_dual, *iterate.primal_upper_bound_gap_dual);
      }
      mu /= (double) number_complementarity_pairs;
      return mu;
   }
}

void Variables::add(const Variables& step, double alpha_primal, double alpha_dual) {
   primals->add(alpha_primal, *step.primals);
   equality_duals->add(alpha_dual, *step.equality_duals);
   inequality_duals->add(alpha_dual, *step.inequality_duals);
   slacks->add(alpha_primal, *step.slacks);
   if (mclow > 0) {
      assert(step.slack_lower_bound_gap->matchesNonZeroPattern(*inequality_lower_bound_indicators) && step.slack_lower_bound_gap_dual->matchesNonZeroPattern(*inequality_lower_bound_indicators));

      slack_lower_bound_gap->add(alpha_primal, *step.slack_lower_bound_gap);
      slack_lower_bound_gap_dual->add(alpha_dual, *step.slack_lower_bound_gap_dual);
   }
   if (mcupp > 0) {
      assert(step.slack_upper_bound_gap->matchesNonZeroPattern(*inequality_upper_bound_indicators) && step.slack_upper_bound_gap_dual->matchesNonZeroPattern(*inequality_upper_bound_indicators));

      slack_upper_bound_gap->add(alpha_primal, *step.slack_upper_bound_gap);
      slack_upper_bound_gap_dual->add(alpha_dual, *step.slack_upper_bound_gap_dual);
   }
   if (nxlow > 0) {
      assert(step.primal_lower_bound_gap->matchesNonZeroPattern(*primal_lower_bound_indicators) && step.primal_lower_bound_gap_dual->matchesNonZeroPattern(*primal_lower_bound_indicators));

      primal_lower_bound_gap->add(alpha_primal, *step.primal_lower_bound_gap);
      primal_lower_bound_gap_dual->add(alpha_dual, *step.primal_lower_bound_gap_dual);
   }
   if (nxupp > 0) {
      assert(step.primal_upper_bound_gap->matchesNonZeroPattern(*primal_upper_bound_indicators) && step.primal_upper_bound_gap_dual->matchesNonZeroPattern(*primal_upper_bound_indicators));

      primal_upper_bound_gap->add(alpha_primal, *step.primal_upper_bound_gap);
      primal_upper_bound_gap_dual->add(alpha_dual, *step.primal_upper_bound_gap_dual);
   }
}

void Variables::add(const Variables& step, double alpha) {
   // equal primal and dual coefficients
   this->add(step, alpha, alpha);
}

void Variables::negate() {
   slacks->negate();
   primals->negate();
   equality_duals->negate();
   inequality_duals->negate();
   if (mclow > 0) {
      slack_lower_bound_gap->negate();
      slack_lower_bound_gap_dual->negate();
   }
   if (mcupp > 0) {
      slack_upper_bound_gap->negate();
      slack_upper_bound_gap_dual->negate();
   }
   if (nxlow > 0) {
      primal_lower_bound_gap->negate();
      primal_lower_bound_gap_dual->negate();
   }
   if (nxupp > 0) {
      primal_upper_bound_gap->negate();
      primal_upper_bound_gap_dual->negate();
   }
}

double Variables::fraction_to_boundary(const Variables& iterate, double fraction) const {
   double length = 1.;
   if (mclow > 0) {
      assert(slack_lower_bound_gap->are_positive(*inequality_lower_bound_indicators));
      assert(slack_lower_bound_gap_dual->are_positive(*inequality_lower_bound_indicators));

      length = std::min(length, slack_lower_bound_gap->fraction_to_boundary(*iterate.slack_lower_bound_gap, fraction));
      length = std::min(length, slack_lower_bound_gap_dual->fraction_to_boundary(*iterate.slack_lower_bound_gap_dual, fraction));
   }
   if (mcupp > 0) {
      assert(slack_upper_bound_gap->are_positive(*inequality_upper_bound_indicators));
      assert(slack_upper_bound_gap_dual->are_positive(*inequality_upper_bound_indicators));

      length = std::min(length, slack_upper_bound_gap->fraction_to_boundary(*iterate.slack_upper_bound_gap, fraction));
      length = std::min(length, slack_upper_bound_gap_dual->fraction_to_boundary(*iterate.slack_upper_bound_gap_dual, fraction));
   }
   if (nxlow > 0) {
      assert(primal_lower_bound_gap->are_positive(*primal_lower_bound_indicators));
      assert(primal_lower_bound_gap_dual->are_positive(*primal_lower_bound_indicators));

      length = std::min(length, primal_lower_bound_gap->fraction_to_boundary(*iterate.primal_lower_bound_gap, fraction));
      length = std::min(length, primal_lower_bound_gap_dual->fraction_to_boundary(*iterate.primal_lower_bound_gap_dual, fraction));
   }
   if (nxupp > 0) {
      assert(primal_upper_bound_gap->are_positive(*primal_upper_bound_indicators));
      assert(primal_upper_bound_gap_dual->are_positive(*primal_upper_bound_indicators));

      length = std::min(length, primal_upper_bound_gap->fraction_to_boundary(*iterate.primal_upper_bound_gap, fraction));
      length = std::min(length, primal_upper_bound_gap_dual->fraction_to_boundary(*iterate.primal_upper_bound_gap_dual, fraction));
   }
   assert(length <= 1.);
   return length;
}

std::pair<double, double> Variables::stepbound_pd(const Variables& iterate) const {
   double primal_length = 1.;
   double dual_length = 1.;

   if (mclow > 0) {
      assert(slack_lower_bound_gap->are_positive(*inequality_lower_bound_indicators));
      assert(slack_lower_bound_gap_dual->are_positive(*inequality_lower_bound_indicators));

      primal_length = std::min(primal_length, slack_lower_bound_gap->fraction_to_boundary(*iterate.slack_lower_bound_gap));
      dual_length = std::min(dual_length, slack_lower_bound_gap_dual->fraction_to_boundary(*iterate.slack_lower_bound_gap_dual));
   }
   if (mcupp > 0) {
      assert(slack_upper_bound_gap->are_positive(*inequality_upper_bound_indicators));
      assert(slack_upper_bound_gap_dual->are_positive(*inequality_upper_bound_indicators));

      primal_length = std::min(primal_length, slack_upper_bound_gap->fraction_to_boundary(*iterate.slack_upper_bound_gap));
      dual_length = std::min(dual_length, slack_upper_bound_gap_dual->fraction_to_boundary(*iterate.slack_upper_bound_gap_dual));
   }
   if (nxlow > 0) {
      assert(primal_lower_bound_gap->are_positive(*primal_lower_bound_indicators));
      assert(primal_lower_bound_gap_dual->are_positive(*primal_lower_bound_indicators));

      primal_length = std::min(primal_length, primal_lower_bound_gap->fraction_to_boundary(*iterate.primal_lower_bound_gap));
      dual_length = std::min(dual_length, primal_lower_bound_gap_dual->fraction_to_boundary(*iterate.primal_lower_bound_gap_dual));
   }
   if (nxupp > 0) {
      assert(primal_upper_bound_gap->are_positive(*primal_upper_bound_indicators));
      assert(primal_upper_bound_gap_dual->are_positive(*primal_upper_bound_indicators));

      primal_length = std::min(primal_length, primal_upper_bound_gap->fraction_to_boundary(*iterate.primal_upper_bound_gap));
      dual_length = std::min(dual_length, primal_upper_bound_gap_dual->fraction_to_boundary(*iterate.primal_upper_bound_gap_dual));
   }
   assert(primal_length <= 1.);
   assert(dual_length <= 1.);
   return std::make_pair(primal_length, dual_length);
}

double
Variables::find_blocking(const Variables& step_in, double& primalValue, double& primalStep, double& dualValue, double& dualStep, int& firstOrSecond) const {
   double alpha = 1.;
   firstOrSecond = 0;

   if (mclow > 0) {
      alpha = slack_lower_bound_gap->find_blocking(*step_in.slack_lower_bound_gap, *slack_lower_bound_gap_dual, *step_in.slack_lower_bound_gap_dual, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   if (mcupp > 0) {
      alpha = slack_upper_bound_gap->find_blocking(*step_in.slack_upper_bound_gap, *slack_upper_bound_gap_dual, *step_in.slack_upper_bound_gap_dual, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   if (nxlow > 0) {
      alpha = primal_lower_bound_gap->find_blocking(*step_in.primal_lower_bound_gap, *primal_lower_bound_gap_dual, *step_in.primal_lower_bound_gap_dual, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   if (nxupp > 0) {
      alpha = primal_upper_bound_gap->find_blocking(*step_in.primal_upper_bound_gap, *primal_upper_bound_gap_dual, *step_in.primal_upper_bound_gap_dual, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   return alpha;
}

void
Variables::find_blocking(const Variables& step_in, double& primalValue, double& primalStep, double& dualValue, double& dualStep, double& primalValue_d,
      double& primalStep_d, double& dualValue_d, double& dualStep_d, double& alphaPrimal, double& alphaDual, bool& primalBlocking,
      bool& dualBlocking) const {
   alphaPrimal = 1., alphaDual = 1.;
   primalBlocking = false, dualBlocking = false;

   if (mclow > 0) {
      slack_lower_bound_gap->find_blocking_pd(*step_in.slack_lower_bound_gap, *slack_lower_bound_gap_dual, *step_in.slack_lower_bound_gap_dual, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d,
            primalStep_d, dualValue_d, dualStep_d, primalBlocking, dualBlocking);
   }

   if (mcupp > 0) {
      slack_upper_bound_gap->find_blocking_pd(*step_in.slack_upper_bound_gap, *slack_upper_bound_gap_dual, *step_in.slack_upper_bound_gap_dual, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d, primalStep_d,
            dualValue_d, dualStep_d, primalBlocking, dualBlocking);
   }

   if (nxlow > 0) {
      primal_lower_bound_gap->find_blocking_pd(*step_in.primal_lower_bound_gap, *primal_lower_bound_gap_dual, *step_in.primal_lower_bound_gap_dual, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d,
            primalStep_d, dualValue_d, dualStep_d, primalBlocking, dualBlocking);
   }

   if (nxupp > 0) {
      primal_upper_bound_gap->find_blocking_pd(*step_in.primal_upper_bound_gap, *primal_upper_bound_gap_dual, *step_in.primal_upper_bound_gap_dual, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d,
            primalStep_d, dualValue_d, dualStep_d, primalBlocking, dualBlocking);
   }
}

void Variables::push_to_interior(double alpha, double beta) {
   slacks->setToZero();
   primals->setToZero();
   equality_duals->setToZero();
   inequality_duals->setToZero();

   if (nxlow > 0) {
      primal_lower_bound_gap->setToConstant(alpha);
      primal_lower_bound_gap->selectNonZeros(*primal_lower_bound_indicators);
      primal_lower_bound_gap_dual->setToConstant(beta);
      primal_lower_bound_gap_dual->selectNonZeros(*primal_lower_bound_indicators);
   }
   if (nxupp > 0) {
      primal_upper_bound_gap->setToConstant(alpha);
      primal_upper_bound_gap->selectNonZeros(*primal_upper_bound_indicators);
      primal_upper_bound_gap_dual->setToConstant(beta);
      primal_upper_bound_gap_dual->selectNonZeros(*primal_upper_bound_indicators);
   }

   if (mclow > 0) {
      slack_lower_bound_gap->setToConstant(alpha);
      slack_lower_bound_gap->selectNonZeros(*inequality_lower_bound_indicators);
      slack_lower_bound_gap_dual->setToConstant(beta);
      slack_lower_bound_gap_dual->selectNonZeros(*inequality_lower_bound_indicators);
   }
   if (mcupp > 0) {
      slack_upper_bound_gap->setToConstant(alpha);
      slack_upper_bound_gap->selectNonZeros(*inequality_upper_bound_indicators);
      slack_upper_bound_gap_dual->setToConstant(beta);
      slack_upper_bound_gap_dual->selectNonZeros(*inequality_upper_bound_indicators);
   }
}

double Variables::violation() const {
   double viol = 0., cmin = 0.;
   int iblock;

   if (nxlow > 0) {
      primal_lower_bound_gap->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;

      primal_lower_bound_gap_dual->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;
   }
   if (nxupp > 0) {
      primal_upper_bound_gap->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;

      primal_upper_bound_gap_dual->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;
   }
   if (mclow > 0) {
      slack_lower_bound_gap->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;

      slack_lower_bound_gap_dual->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;
   }
   if (mcupp > 0) {
      slack_upper_bound_gap->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;

      slack_upper_bound_gap_dual->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;
   }
   return -viol;
}

void Variables::shift_bound_variables(double alpha, double beta) {
   if (nxlow > 0) {
      primal_lower_bound_gap->add_constant(alpha, *primal_lower_bound_indicators);
      primal_lower_bound_gap_dual->add_constant(beta, *primal_lower_bound_indicators);
   }
   if (nxupp > 0) {
      primal_upper_bound_gap->add_constant(alpha, *primal_upper_bound_indicators);
      primal_upper_bound_gap_dual->add_constant(beta, *primal_upper_bound_indicators);
   }
   if (mclow > 0) {
      slack_lower_bound_gap->add_constant(alpha, *inequality_lower_bound_indicators);
      slack_lower_bound_gap_dual->add_constant(beta, *inequality_lower_bound_indicators);
   }
   if (mcupp > 0) {
      slack_upper_bound_gap->add_constant(alpha, *inequality_upper_bound_indicators);
      slack_upper_bound_gap_dual->add_constant(beta, *inequality_upper_bound_indicators);
   }
}

void Variables::copy(const Variables& b) {

   slacks->copyFrom(*b.slacks);
   if (nxlow > 0) {
      primal_lower_bound_gap->copyFrom(*b.primal_lower_bound_gap);
      primal_lower_bound_gap_dual->copyFrom(*b.primal_lower_bound_gap_dual);
   }
   if (nxupp > 0) {
      primal_upper_bound_gap->copyFrom(*b.primal_upper_bound_gap);
      primal_upper_bound_gap_dual->copyFrom(*b.primal_upper_bound_gap_dual);
   }
   if (mclow > 0) {
      slack_lower_bound_gap->copyFrom(*b.slack_lower_bound_gap);
      slack_lower_bound_gap_dual->copyFrom(*b.slack_lower_bound_gap_dual);
   }
   if (mcupp > 0) {
      slack_upper_bound_gap->copyFrom(*b.slack_upper_bound_gap);
      slack_upper_bound_gap_dual->copyFrom(*b.slack_upper_bound_gap_dual);
   }
   primals->copyFrom(*b.primals);
   equality_duals->copyFrom(*b.equality_duals);
   inequality_duals->copyFrom(*b.inequality_duals);

}

double Variables::one_norm() const {
   double norm = primals->one_norm();
   norm += slacks->one_norm();
   norm += equality_duals->one_norm();
   norm += inequality_duals->one_norm();

   norm += primal_lower_bound_gap->one_norm();
   norm += primal_upper_bound_gap_dual->one_norm();
   norm += primal_upper_bound_gap->one_norm();
   norm += primal_lower_bound_gap_dual->one_norm();
   norm += slack_lower_bound_gap->one_norm();
   norm += slack_lower_bound_gap_dual->one_norm();
   norm += slack_upper_bound_gap->one_norm();
   norm += slack_upper_bound_gap_dual->one_norm();

   return norm;
}

double Variables::inf_norm() const {
   double norm = 0.0;
   double temp = primals->inf_norm();
   if (temp > norm)
      norm = temp;
   temp = slacks->inf_norm();
   if (temp > norm)
      norm = temp;
   temp = equality_duals->inf_norm();
   if (temp > norm)
      norm = temp;
   temp = inequality_duals->inf_norm();
   if (temp > norm)
      norm = temp;

   temp = primal_lower_bound_gap->inf_norm();
   if (temp > norm)
      norm = temp;
   temp = primal_upper_bound_gap_dual->inf_norm();
   if (temp > norm)
      norm = temp;

   temp = primal_upper_bound_gap->inf_norm();
   if (temp > norm)
      norm = temp;
   temp = primal_lower_bound_gap_dual->inf_norm();
   if (temp > norm)
      norm = temp;

   temp = slack_lower_bound_gap->inf_norm();
   if (temp > norm)
      norm = temp;
   temp = slack_lower_bound_gap_dual->inf_norm();
   if (temp > norm)
      norm = temp;

   temp = slack_upper_bound_gap->inf_norm();
   if (temp > norm)
      norm = temp;
   temp = slack_upper_bound_gap_dual->inf_norm();
   if (temp > norm)
      norm = temp;

   return norm;
}

void Variables::set_to_zero() {
   primals->setToZero();
   slacks->setToZero();
   equality_duals->setToZero();
   inequality_duals->setToZero();

   primal_lower_bound_gap->setToZero();
   primal_lower_bound_gap_dual->setToZero();

   primal_upper_bound_gap->setToZero();
   primal_upper_bound_gap_dual->setToZero();

   slack_lower_bound_gap->setToZero();
   slack_lower_bound_gap_dual->setToZero();

   slack_upper_bound_gap->setToZero();
   slack_upper_bound_gap_dual->setToZero();
}

int Variables::valid_non_zero_pattern() const {
   if (nxlow > 0 && (!primal_lower_bound_gap->matchesNonZeroPattern(*primal_lower_bound_indicators) || !primal_lower_bound_gap_dual->matchesNonZeroPattern(*primal_lower_bound_indicators))) {

      if (!primal_lower_bound_gap->matchesNonZeroPattern(*primal_lower_bound_indicators))
         printf("invalidNonZeroPattern v\n");
      if (!primal_lower_bound_gap_dual->matchesNonZeroPattern(*primal_lower_bound_indicators))
         printf("invalidNonZeroPattern gamma\n");
      return 0;
   }

   if (nxupp > 0 && (!primal_upper_bound_gap->matchesNonZeroPattern(*primal_upper_bound_indicators) || !primal_upper_bound_gap_dual->matchesNonZeroPattern(*primal_upper_bound_indicators))) {
      if (!primal_upper_bound_gap->matchesNonZeroPattern(*primal_upper_bound_indicators))
         printf("invalidNonZeroPattern w\n");
      if (!primal_upper_bound_gap_dual->matchesNonZeroPattern(*primal_upper_bound_indicators))
         printf("invalidNonZeroPattern phi\n");
      return 0;
   }
   if (mclow > 0 && (!slack_lower_bound_gap->matchesNonZeroPattern(*inequality_lower_bound_indicators) || !slack_lower_bound_gap_dual->matchesNonZeroPattern(*inequality_lower_bound_indicators))) {
      if (!slack_lower_bound_gap->matchesNonZeroPattern(*inequality_lower_bound_indicators))
         printf("invalidNonZeroPattern t\n");
      if (!slack_lower_bound_gap_dual->matchesNonZeroPattern(*inequality_lower_bound_indicators))
         printf("invalidNonZeroPattern lambda\n");
      return 0;
   }

   if (mcupp > 0 && (!slack_upper_bound_gap->matchesNonZeroPattern(*inequality_upper_bound_indicators) || !slack_upper_bound_gap_dual->matchesNonZeroPattern(*inequality_upper_bound_indicators))) {
      if (!slack_upper_bound_gap->matchesNonZeroPattern(*inequality_upper_bound_indicators))
         printf("invalidNonZeroPattern u\n");
      if (!primal_upper_bound_gap_dual->matchesNonZeroPattern(*inequality_upper_bound_indicators))
         printf("invalidNonZeroPattern phi\n");
      return 0;
   }

   return 1;
}

void Variables::unscale_solution(Problem* problem) {
// Modifying sx is equivalent to modifying x
   auto& sx = (DenseVector<double>&) *this->primals;

// x = D * x'
   sx.componentMult(problem->scale());
}

void Variables::unscale_bounds(Problem* problem) {
   auto& sxlow = (DenseVector<double>&) problem->x_lower_bound();
   auto& sxupp = (DenseVector<double>&) problem->x_upper_bound();

// l = D * l'
   sxlow.componentMult(problem->scale());

// u = D * u'
   sxupp.componentMult(problem->scale());
}

void Variables::print_solution(MpsReader* reader, Problem* problem, int& iErr) {
   assert(primals->isKindOf(kDenseVector)); // Otherwise this routine

   DenseVector<double> g(nx);
   problem->get_objective_gradient(g);
   problem->hessian_multiplication(1.0, g, 0.5, *primals);
   double objective = g.dotProductWith(*primals);

   auto& sx = dynamic_cast<DenseVector<double>&>(*this->primals);
   auto& sxlow = dynamic_cast<DenseVector<double>&>(problem->x_lower_bound());
   auto& sixlow = dynamic_cast<DenseVector<double>&>(problem->has_x_lower_bound());
   auto& sxupp = dynamic_cast<DenseVector<double>&>(problem->x_upper_bound());
   auto& sixupp = dynamic_cast<DenseVector<double>&>(problem->has_x_upper_bound());
   auto& sgamma = dynamic_cast<DenseVector<double>&>(*this->primal_lower_bound_gap_dual);
   auto& sphi = dynamic_cast<DenseVector<double>&>(*this->primal_upper_bound_gap_dual);
   auto& sy = dynamic_cast<DenseVector<double>&>(*this->equality_duals);
   auto& ss = dynamic_cast<DenseVector<double>&>(*this->slacks);
   auto& slambda = dynamic_cast<DenseVector<double>&>(*this->slack_lower_bound_gap_dual);
   auto& spi = dynamic_cast<DenseVector<double>&>(*this->slack_upper_bound_gap_dual);
   auto& sz = dynamic_cast<DenseVector<double>&>(*this->inequality_duals);
   auto& sclow = dynamic_cast<DenseVector<double>&>(problem->s_lower_bound());
   auto& siclow = dynamic_cast<DenseVector<double>&>(problem->has_s_lower_bound());
   auto& scupp = dynamic_cast<DenseVector<double>&>(problem->s_upper_bound());
   auto& sicupp = dynamic_cast<DenseVector<double>&>(problem->has_s_upper_bound());

   char* cxupp = new char[nx];
   char* cxlow = new char[nx];
   for (int j = 0; j < nx; j++) {
      if (nxupp > 0 && sixupp[j] != 0) {
         cxupp[j] = 1;
      }
      else {
         cxupp[j] = 0;
      }
      if (nxlow > 0 && sixlow[j] != 0) {
         cxlow[j] = 1;
      }
      else {
         cxlow[j] = 0;
      }
   }
   char* cclow, * ccupp;
   if (mz <= 0) {
      cclow = nullptr;
      ccupp = nullptr;
   }
   else {
      cclow = new char[mz];
      ccupp = new char[mz];
      for (int i = 0; i < mz; i++) {
         if (mclow > 0 && siclow[i] != 0.0) {
            cclow[i] = 1;
         }
         else {
            cclow[i] = 0;
         }
         if (mcupp > 0 && sicupp[i] != 0.0) {
            ccupp[i] = 1;
         }
         else {
            ccupp[i] = 0;
         }
      }
   }

   if (reader->scalingOption == 1) {
      // Unscale the solution and bounds before printing
      this->unscale_solution(problem);
      this->unscale_bounds(problem);
   }

   reader->printSolution(sx.elements(), nx, sxlow.elements(), cxlow, sxupp.elements(), cxupp, sgamma.elements(), sphi.elements(), sy.elements(), my,
         ss.elements(), mz, sclow.elements(), cclow, scupp.elements(), ccupp, slambda.elements(), spi.elements(), sz.elements(), objective, iErr);
   delete[] cclow;
   delete[] ccupp;
   delete[] cxlow;
   delete[] cxupp;
}

// default implementation for Variables::print() prints abusive
// message. Since we don't have any knowledge of how the variables are
// stored at this top level, we can't do much except print their
// dimensions.

void Variables::print() const {
   std::cout << " Complementary Variables = " << number_complementarity_pairs << std::endl;
   std::cout << "(Cannot tell you more at this level)" << std::endl;
}

void Variables::print_norms(bool print_bound_gaps_and_duals) const {

   const double primals_infnorm = this->primals->inf_norm();
   const double slacks_infnorm = this->slacks->inf_norm();
   const double equality_duals_infnorm = this->equality_duals->inf_norm();
   const double inequality_duals_infnorm = this->inequality_duals->inf_norm();

   if (PIPS_MPIgetRank() == 0) {
      std::cout << "Inf norm of x: " << primals_infnorm << "\n";
      std::cout << "Inf norm of s: " << slacks_infnorm << "\n";
      std::cout << "Inf norm of eq. duals: " << equality_duals_infnorm << "\n";
      std::cout << "Inf norm of ineq. duals: " << inequality_duals_infnorm << "\n";
   }

   if (print_bound_gaps_and_duals) {
      const double primal_lower_bound_gap_infnorm = primal_lower_bound_gap->inf_norm();
      const double primal_lower_bound_gap_dual_infnorm = primal_lower_bound_gap_dual->inf_norm();

      const double primal_upper_bound_gap_infnorm = primal_upper_bound_gap->inf_norm();
      const double primal_upper_bound_gap_dual_infnorm = primal_upper_bound_gap_dual->inf_norm();

      const double slack_lower_bound_gap_infnorm = slack_lower_bound_gap->inf_norm();
      const double slack_lower_bound_gap_dual_infnorm = slack_lower_bound_gap_dual->inf_norm();

      const double slack_upper_bound_gap_infnorm = slack_upper_bound_gap->inf_norm();
      const double slack_upper_bound_gap_dual_infnorm = slack_upper_bound_gap_dual->inf_norm();

      if (PIPS_MPIgetRank() == 0) {
         std::cout << "Inf norm of x lower bound: " << primal_lower_bound_gap_infnorm << "\n";
         std::cout << "Inf norm of x lower bound dual: " << primal_lower_bound_gap_dual_infnorm << "\n";

         std::cout << "Inf norm of x upper bound: " << primal_upper_bound_gap_infnorm << "\n";
         std::cout << "Inf norm of x upper bound dual: " << primal_upper_bound_gap_dual_infnorm << "\n";

         std::cout << "Inf norm of s lower bound gap: " << slack_lower_bound_gap_infnorm << "\n";
         std::cout << "Inf norm of s lower bound gap dual: " << slack_lower_bound_gap_dual_infnorm << "\n";

         std::cout << "Inf norm of s upper bound gap: " << slack_upper_bound_gap_infnorm << "\n";
         std::cout << "Inf norm of s upper bound gap dual: " << slack_upper_bound_gap_dual_infnorm << "\n";
      }
   }
}