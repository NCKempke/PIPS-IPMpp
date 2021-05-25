#include <iostream>
#include <SimpleVector.h>
#include "Variables.h"
#include "Vector.hpp"
#include "SmartPointer.h"

#include "Problem.h"
#include "MpsReader.h"

Variables::Variables(Vector<double>* x_in, Vector<double>* s_in, Vector<double>* y_in, Vector<double>* z_in, Vector<double>* v_in,
      Vector<double>* gamma_in, Vector<double>* w_in, Vector<double>* phi_in, Vector<double>* t_in, Vector<double>* lambda_in, Vector<double>* u_in,
      Vector<double>* pi_in, Vector<double>* ixlow_in, Vector<double>* ixupp_in, Vector<double>* iclow_in, Vector<double>* icupp_in) {
   SpReferTo(primals, x_in);
   SpReferTo(slacks, s_in);
   SpReferTo(equality_duals, y_in);
   SpReferTo(inequality_duals, z_in);
   SpReferTo(primal_lower_bound_gap, v_in);
   SpReferTo(primal_upper_bound_gap_dual, phi_in);
   SpReferTo(primal_upper_bound_gap, w_in);
   SpReferTo(primal_lower_bound_gap_dual, gamma_in);
   SpReferTo(slack_lower_bound_gap, t_in);
   SpReferTo(slack_lower_bound_gap_dual, lambda_in);
   SpReferTo(slack_upper_bound_gap, u_in);
   SpReferTo(slack_upper_bound_gap_dual, pi_in);
   SpReferTo(ixlow, ixlow_in);
   SpReferTo(ixupp, ixupp_in);
   SpReferTo(iclow, iclow_in);
   SpReferTo(icupp, icupp_in);

   nx = primals->length();
   my = equality_duals->length();
   mz = inequality_duals->length();

   assert(nx == ixlow->length() || 0 == ixlow->length());
   assert(nx == ixlow->length() || 0 == ixlow->length());
   assert(mz == iclow->length() || 0 == iclow->length());
   assert(mz == icupp->length() || 0 == icupp->length());

   nxlow = ixlow->number_nonzeros();
   nxupp = ixupp->number_nonzeros();
   mclow = iclow->number_nonzeros();
   mcupp = icupp->number_nonzeros();
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

Variables::Variables(const Variables& vars) {
   ixlow = SmartPointer<Vector<double> >(vars.ixlow->cloneFull());
   ixupp = SmartPointer<Vector<double> >(vars.ixupp->cloneFull());
   iclow = SmartPointer<Vector<double> >(vars.iclow->cloneFull());
   icupp = SmartPointer<Vector<double> >(vars.icupp->cloneFull());

   nx = vars.nx;
   my = vars.my;
   mz = vars.mz;

   nxlow = ixlow->number_nonzeros();
   nxupp = ixupp->number_nonzeros();
   mclow = iclow->number_nonzeros();
   mcupp = icupp->number_nonzeros();

   slacks = SmartPointer<Vector<double> >(vars.slacks->cloneFull());

   slack_lower_bound_gap = SmartPointer<Vector<double> >(vars.slack_lower_bound_gap->cloneFull());
   slack_lower_bound_gap_dual = SmartPointer<Vector<double> >(vars.slack_lower_bound_gap_dual->cloneFull());

   slack_upper_bound_gap = SmartPointer<Vector<double> >(vars.slack_upper_bound_gap->cloneFull());
   slack_upper_bound_gap_dual = SmartPointer<Vector<double> >(vars.slack_upper_bound_gap_dual->cloneFull());

   primal_lower_bound_gap = SmartPointer<Vector<double> >(vars.primal_lower_bound_gap->cloneFull());
   primal_lower_bound_gap_dual = SmartPointer<Vector<double> >(vars.primal_lower_bound_gap_dual->cloneFull());

   primal_upper_bound_gap = SmartPointer<Vector<double> >(vars.primal_upper_bound_gap->cloneFull());
   primal_upper_bound_gap_dual = SmartPointer<Vector<double> >(vars.primal_upper_bound_gap_dual->cloneFull());

   primals = SmartPointer<Vector<double> >(vars.primals->cloneFull());
   equality_duals = SmartPointer<Vector<double> >(vars.equality_duals->cloneFull());
   inequality_duals = SmartPointer<Vector<double> >(vars.inequality_duals->cloneFull());
   number_complementarity_pairs = mclow + mcupp + nxlow + nxupp;
}

double Variables::get_average_distance_to_bound_for_converged_vars(const Problem&, double tol) const {
   assert(0 < tol);

   double sum_small_distance = 0.0;
   int n_close = 0;
   primal_lower_bound_gap->getSumCountIfSmall(tol, sum_small_distance, n_close, &*ixlow);
   primal_upper_bound_gap->getSumCountIfSmall(tol, sum_small_distance, n_close, &*ixupp);
   slack_upper_bound_gap->getSumCountIfSmall(tol, sum_small_distance, n_close, &*icupp);
   slack_lower_bound_gap->getSumCountIfSmall(tol, sum_small_distance, n_close, &*iclow);

   if (n_close == 0)
      return std::numeric_limits<double>::infinity();
   else
      return sum_small_distance / (double) n_close;
}


void Variables::push_slacks_from_bound(double tol, double amount) {
   if (nxlow > 0)
      primal_lower_bound_gap->pushAwayFromZero(tol, amount, &*ixlow);
   if (nxupp > 0)
      primal_upper_bound_gap->pushAwayFromZero(tol, amount, &*ixupp);
   if (mclow > 0)
      slack_lower_bound_gap->pushAwayFromZero(tol, amount, &*iclow);
   if (mcupp > 0)
      slack_upper_bound_gap->pushAwayFromZero(tol, amount, &*icupp);
}


double Variables::mu() {
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

      mu /= number_complementarity_pairs;
      return mu;
   }
}

double Variables::mustep_pd(const Variables* iterate, double alpha_primal, double alpha_dual) {
   double mu = 0.;
   if (number_complementarity_pairs == 0) {
      return 0.;
   }
   else {
      if (mclow > 0) {
         mu += slack_lower_bound_gap->shiftedDotProductWith(alpha_primal, *iterate->slack_lower_bound_gap, *slack_lower_bound_gap_dual, alpha_dual, *iterate->slack_lower_bound_gap_dual);
      }
      if (mcupp > 0) {
         mu += slack_upper_bound_gap->shiftedDotProductWith(alpha_primal, *iterate->slack_upper_bound_gap, *slack_upper_bound_gap_dual, alpha_dual, *iterate->slack_upper_bound_gap_dual);
      }
      if (nxlow > 0) {
         mu += primal_lower_bound_gap->shiftedDotProductWith(alpha_primal, *iterate->primal_lower_bound_gap, *primal_lower_bound_gap_dual, alpha_dual, *iterate->primal_lower_bound_gap_dual);
      }
      if (nxupp > 0) {
         mu += primal_upper_bound_gap->shiftedDotProductWith(alpha_primal, *iterate->primal_upper_bound_gap, *primal_upper_bound_gap_dual, alpha_dual, *iterate->primal_upper_bound_gap_dual);
      }
      mu /= number_complementarity_pairs;
      return mu;
   }
}

void Variables::saxpy(const Variables* iterate, double alpha) {
   primals->axpy(alpha, *iterate->primals);
   equality_duals->axpy(alpha, *iterate->equality_duals);
   inequality_duals->axpy(alpha, *iterate->inequality_duals);
   slacks->axpy(alpha, *iterate->slacks);
   if (mclow > 0) {
      assert(iterate->slack_lower_bound_gap->matchesNonZeroPattern(*iclow) && iterate->slack_lower_bound_gap_dual->matchesNonZeroPattern(*iclow));

      slack_lower_bound_gap->axpy(alpha, *iterate->slack_lower_bound_gap);
      slack_lower_bound_gap_dual->axpy(alpha, *iterate->slack_lower_bound_gap_dual);
   }
   if (mcupp > 0) {
      assert(iterate->slack_upper_bound_gap->matchesNonZeroPattern(*icupp) && iterate->slack_upper_bound_gap_dual->matchesNonZeroPattern(*icupp));

      slack_upper_bound_gap->axpy(alpha, *iterate->slack_upper_bound_gap);
      slack_upper_bound_gap_dual->axpy(alpha, *iterate->slack_upper_bound_gap_dual);
   }
   if (nxlow > 0) {
      assert(iterate->primal_lower_bound_gap->matchesNonZeroPattern(*ixlow) && iterate->primal_lower_bound_gap_dual->matchesNonZeroPattern(*ixlow));

      primal_lower_bound_gap->axpy(alpha, *iterate->primal_lower_bound_gap);
      primal_lower_bound_gap_dual->axpy(alpha, *iterate->primal_lower_bound_gap_dual);
   }
   if (nxupp > 0) {
      assert(iterate->primal_upper_bound_gap->matchesNonZeroPattern(*ixupp) && iterate->primal_upper_bound_gap_dual->matchesNonZeroPattern(*ixupp));

      primal_upper_bound_gap->axpy(alpha, *iterate->primal_upper_bound_gap);
      primal_upper_bound_gap_dual->axpy(alpha, *iterate->primal_upper_bound_gap_dual);
   }
}

void Variables::saxpy_pd(const Variables* iterate, double alpha_primal, double alpha_dual) {
   primals->axpy(alpha_primal, *iterate->primals);
   equality_duals->axpy(alpha_dual, *iterate->equality_duals);
   inequality_duals->axpy(alpha_dual, *iterate->inequality_duals);
   slacks->axpy(alpha_primal, *iterate->slacks);
   if (mclow > 0) {
      assert(iterate->slack_lower_bound_gap->matchesNonZeroPattern(*iclow) && iterate->slack_lower_bound_gap_dual->matchesNonZeroPattern(*iclow));

      slack_lower_bound_gap->axpy(alpha_primal, *iterate->slack_lower_bound_gap);
      slack_lower_bound_gap_dual->axpy(alpha_dual, *iterate->slack_lower_bound_gap_dual);
   }
   if (mcupp > 0) {
      assert(iterate->slack_upper_bound_gap->matchesNonZeroPattern(*icupp) && iterate->slack_upper_bound_gap_dual->matchesNonZeroPattern(*icupp));

      slack_upper_bound_gap->axpy(alpha_primal, *iterate->slack_upper_bound_gap);
      slack_upper_bound_gap_dual->axpy(alpha_dual, *iterate->slack_upper_bound_gap_dual);
   }
   if (nxlow > 0) {
      assert(iterate->primal_lower_bound_gap->matchesNonZeroPattern(*ixlow) && iterate->primal_lower_bound_gap_dual->matchesNonZeroPattern(*ixlow));

      primal_lower_bound_gap->axpy(alpha_primal, *iterate->primal_lower_bound_gap);
      primal_lower_bound_gap_dual->axpy(alpha_dual, *iterate->primal_lower_bound_gap_dual);
   }
   if (nxupp > 0) {
      assert(iterate->primal_upper_bound_gap->matchesNonZeroPattern(*ixupp) && iterate->primal_upper_bound_gap_dual->matchesNonZeroPattern(*ixupp));

      primal_upper_bound_gap->axpy(alpha_primal, *iterate->primal_upper_bound_gap);
      primal_upper_bound_gap_dual->axpy(alpha_dual, *iterate->primal_upper_bound_gap_dual);
   }
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

double Variables::stepbound(const Variables* iterate) {
   double max_step = 1.0;
   if (mclow > 0) {
      assert(slack_lower_bound_gap->somePositive(*iclow));
      assert(slack_lower_bound_gap_dual->somePositive(*iclow));

      max_step = slack_lower_bound_gap->stepbound(*iterate->slack_lower_bound_gap, max_step);
      max_step = slack_lower_bound_gap_dual->stepbound(*iterate->slack_lower_bound_gap_dual, max_step);
   }

   if (mcupp > 0) {
      assert(slack_upper_bound_gap->somePositive(*icupp));
      assert(slack_upper_bound_gap_dual->somePositive(*icupp));

      max_step = slack_upper_bound_gap->stepbound(*iterate->slack_upper_bound_gap, max_step);
      max_step = slack_upper_bound_gap_dual->stepbound(*iterate->slack_upper_bound_gap_dual, max_step);
   }

   if (nxlow > 0) {
      assert(primal_lower_bound_gap->somePositive(*ixlow));
      assert(primal_lower_bound_gap_dual->somePositive(*ixlow));

      max_step = primal_lower_bound_gap->stepbound(*iterate->primal_lower_bound_gap, max_step);
      max_step = primal_lower_bound_gap_dual->stepbound(*iterate->primal_lower_bound_gap_dual, max_step);
   }

   if (nxupp > 0) {
      assert(primal_upper_bound_gap->somePositive(*ixupp));
      assert(primal_upper_bound_gap_dual->somePositive(*ixupp));

      max_step = primal_upper_bound_gap->stepbound(*iterate->primal_upper_bound_gap, max_step);
      max_step = primal_upper_bound_gap_dual->stepbound(*iterate->primal_upper_bound_gap_dual, max_step);
   }

   assert(max_step <= 1.0);
   return max_step;
}

std::pair<double, double> Variables::stepbound_pd(const Variables* iterate) {
   double maxStep_primal = 1.0;
   double maxStep_dual = 1.0;

   if (mclow > 0) {
      assert(slack_lower_bound_gap->somePositive(*iclow));
      assert(slack_lower_bound_gap_dual->somePositive(*iclow));

      maxStep_primal = slack_lower_bound_gap->stepbound(*iterate->slack_lower_bound_gap, maxStep_primal);
      maxStep_dual = slack_lower_bound_gap_dual->stepbound(*iterate->slack_lower_bound_gap_dual, maxStep_dual);
   }

   if (mcupp > 0) {
      assert(slack_upper_bound_gap->somePositive(*icupp));
      assert(slack_upper_bound_gap_dual->somePositive(*icupp));

      maxStep_primal = slack_upper_bound_gap->stepbound(*iterate->slack_upper_bound_gap, maxStep_primal);
      maxStep_dual = slack_upper_bound_gap_dual->stepbound(*iterate->slack_upper_bound_gap_dual, maxStep_dual);
   }

   if (nxlow > 0) {
      assert(primal_lower_bound_gap->somePositive(*ixlow));
      assert(primal_lower_bound_gap_dual->somePositive(*ixlow));

      maxStep_primal = primal_lower_bound_gap->stepbound(*iterate->primal_lower_bound_gap, maxStep_primal);
      maxStep_dual = primal_lower_bound_gap_dual->stepbound(*iterate->primal_lower_bound_gap_dual, maxStep_dual);
   }

   if (nxupp > 0) {
      assert(primal_upper_bound_gap->somePositive(*ixupp));
      assert(primal_upper_bound_gap_dual->somePositive(*ixupp));

      maxStep_primal = primal_upper_bound_gap->stepbound(*iterate->primal_upper_bound_gap, maxStep_primal);
      maxStep_dual = primal_upper_bound_gap_dual->stepbound(*iterate->primal_upper_bound_gap_dual, maxStep_dual);
   }

   assert(maxStep_primal <= 1.0);
   assert(maxStep_dual <= 1.0);
   return std::make_pair(maxStep_primal, maxStep_dual);
}

double
Variables::find_blocking(const Variables* step_in, double& primalValue, double& primalStep, double& dualValue, double& dualStep, int& firstOrSecond) {
   double alpha = 1.;
   firstOrSecond = 0;

   if (mclow > 0) {
      alpha = slack_lower_bound_gap->find_blocking(*step_in->slack_lower_bound_gap, *slack_lower_bound_gap_dual, *step_in->slack_lower_bound_gap_dual, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   if (mcupp > 0) {
      alpha = slack_upper_bound_gap->find_blocking(*step_in->slack_upper_bound_gap, *slack_upper_bound_gap_dual, *step_in->slack_upper_bound_gap_dual, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   if (nxlow > 0) {
      alpha = primal_lower_bound_gap->find_blocking(*step_in->primal_lower_bound_gap, *primal_lower_bound_gap_dual, *step_in->primal_lower_bound_gap_dual, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   if (nxupp > 0) {
      alpha = primal_upper_bound_gap->find_blocking(*step_in->primal_upper_bound_gap, *primal_upper_bound_gap_dual, *step_in->primal_upper_bound_gap_dual, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   return alpha;
}

void
Variables::find_blocking(const Variables* step_in, double& primalValue, double& primalStep, double& dualValue, double& dualStep, double& primalValue_d,
      double& primalStep_d, double& dualValue_d, double& dualStep_d, double& alphaPrimal, double& alphaDual, bool& primalBlocking,
      bool& dualBlocking) {
   alphaPrimal = 1.0, alphaDual = 1.0;
   primalBlocking = false, dualBlocking = false;

   if (mclow > 0) {
      slack_lower_bound_gap->find_blocking_pd(*step_in->slack_lower_bound_gap, *slack_lower_bound_gap_dual, *step_in->slack_lower_bound_gap_dual, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d,
            primalStep_d, dualValue_d, dualStep_d, primalBlocking, dualBlocking);
   }

   if (mcupp > 0) {
      slack_upper_bound_gap->find_blocking_pd(*step_in->slack_upper_bound_gap, *slack_upper_bound_gap_dual, *step_in->slack_upper_bound_gap_dual, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d, primalStep_d,
            dualValue_d, dualStep_d, primalBlocking, dualBlocking);
   }

   if (nxlow > 0) {
      primal_lower_bound_gap->find_blocking_pd(*step_in->primal_lower_bound_gap, *primal_lower_bound_gap_dual, *step_in->primal_lower_bound_gap_dual, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d,
            primalStep_d, dualValue_d, dualStep_d, primalBlocking, dualBlocking);
   }

   if (nxupp > 0) {
      primal_upper_bound_gap->find_blocking_pd(*step_in->primal_upper_bound_gap, *primal_upper_bound_gap_dual, *step_in->primal_upper_bound_gap_dual, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d,
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
      primal_lower_bound_gap->selectNonZeros(*ixlow);
      primal_lower_bound_gap_dual->setToConstant(beta);
      primal_lower_bound_gap_dual->selectNonZeros(*ixlow);
   }
   if (nxupp > 0) {
      primal_upper_bound_gap->setToConstant(alpha);
      primal_upper_bound_gap->selectNonZeros(*ixupp);
      primal_upper_bound_gap_dual->setToConstant(beta);
      primal_upper_bound_gap_dual->selectNonZeros(*ixupp);
   }

   if (mclow > 0) {
      slack_lower_bound_gap->setToConstant(alpha);
      slack_lower_bound_gap->selectNonZeros(*iclow);
      slack_lower_bound_gap_dual->setToConstant(beta);
      slack_lower_bound_gap_dual->selectNonZeros(*iclow);
   }
   if (mcupp > 0) {
      slack_upper_bound_gap->setToConstant(alpha);
      slack_upper_bound_gap->selectNonZeros(*icupp);
      slack_upper_bound_gap_dual->setToConstant(beta);
      slack_upper_bound_gap_dual->selectNonZeros(*icupp);
   }
}

double Variables::violation() {
   double viol = 0.0, cmin = 0.0;
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
      primal_lower_bound_gap->add_constant(alpha, *ixlow);
      primal_lower_bound_gap_dual->add_constant(beta, *ixlow);
   }
   if (nxupp > 0) {
      primal_upper_bound_gap->add_constant(alpha, *ixupp);
      primal_upper_bound_gap_dual->add_constant(beta, *ixupp);
   }
   if (mclow > 0) {
      slack_lower_bound_gap->add_constant(alpha, *iclow);
      slack_lower_bound_gap_dual->add_constant(beta, *iclow);
   }
   if (mcupp > 0) {
      slack_upper_bound_gap->add_constant(alpha, *icupp);
      slack_upper_bound_gap_dual->add_constant(beta, *icupp);
   }
}

void Variables::copy(const Variables* b_in) {
   const Variables* b = (const Variables*) b_in;

   slacks->copyFrom(*b->slacks);
   if (nxlow > 0) {
      primal_lower_bound_gap->copyFrom(*b->primal_lower_bound_gap);
      primal_lower_bound_gap_dual->copyFrom(*b->primal_lower_bound_gap_dual);
   }
   if (nxupp > 0) {
      primal_upper_bound_gap->copyFrom(*b->primal_upper_bound_gap);
      primal_upper_bound_gap_dual->copyFrom(*b->primal_upper_bound_gap_dual);
   }
   if (mclow > 0) {
      slack_lower_bound_gap->copyFrom(*b->slack_lower_bound_gap);
      slack_lower_bound_gap_dual->copyFrom(*b->slack_lower_bound_gap_dual);
   }
   if (mcupp > 0) {
      slack_upper_bound_gap->copyFrom(*b->slack_upper_bound_gap);
      slack_upper_bound_gap_dual->copyFrom(*b->slack_upper_bound_gap_dual);
   }
   primals->copyFrom(*b->primals);
   equality_duals->copyFrom(*b->equality_duals);
   inequality_duals->copyFrom(*b->inequality_duals);

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

int Variables::valid_non_zero_pattern() {
   if (nxlow > 0 && (!primal_lower_bound_gap->matchesNonZeroPattern(*ixlow) || !primal_lower_bound_gap_dual->matchesNonZeroPattern(*ixlow))) {

      if (!primal_lower_bound_gap->matchesNonZeroPattern(*ixlow))
         printf("invalidNonZeroPattern v\n");
      if (!primal_lower_bound_gap_dual->matchesNonZeroPattern(*ixlow))
         printf("invalidNonZeroPattern gamma\n");
      return 0;
   }

   if (nxupp > 0 && (!primal_upper_bound_gap->matchesNonZeroPattern(*ixupp) || !primal_upper_bound_gap_dual->matchesNonZeroPattern(*ixupp))) {
      if (!primal_upper_bound_gap->matchesNonZeroPattern(*ixupp))
         printf("invalidNonZeroPattern w\n");
      if (!primal_upper_bound_gap_dual->matchesNonZeroPattern(*ixupp))
         printf("invalidNonZeroPattern phi\n");
      return 0;
   }
   if (mclow > 0 && (!slack_lower_bound_gap->matchesNonZeroPattern(*iclow) || !slack_lower_bound_gap_dual->matchesNonZeroPattern(*iclow))) {
      if (!slack_lower_bound_gap->matchesNonZeroPattern(*iclow))
         printf("invalidNonZeroPattern t\n");
      if (!slack_lower_bound_gap_dual->matchesNonZeroPattern(*iclow))
         printf("invalidNonZeroPattern lambda\n");
      return 0;
   }

   if (mcupp > 0 && (!slack_upper_bound_gap->matchesNonZeroPattern(*icupp) || !slack_upper_bound_gap_dual->matchesNonZeroPattern(*icupp))) {
      if (!slack_upper_bound_gap->matchesNonZeroPattern(*icupp))
         printf("invalidNonZeroPattern u\n");
      if (!primal_upper_bound_gap_dual->matchesNonZeroPattern(*icupp))
         printf("invalidNonZeroPattern phi\n");
      return 0;
   }

   return 1;
}

void Variables::unscale_solution(Problem* problem) {
// Modifying sx is equivalent to modifying x
   SimpleVector<double>& sx = (SimpleVector<double>&) *this->primals;

// x = D * x'
   sx.componentMult(problem->scale());
}

void Variables::unscale_bounds(Problem* problem) {
   SimpleVector<double>& sxlow = (SimpleVector<double>&) problem->xlowerBound();
   SimpleVector<double>& sxupp = (SimpleVector<double>&) problem->xupperBound();

// l = D * l'
   sxlow.componentMult(problem->scale());

// u = D * u'
   sxupp.componentMult(problem->scale());
}

void Variables::print_solution(MpsReader* reader, Problem* problem, int& iErr) {
   assert(primals->isKindOf(kSimpleVector)); // Otherwise this routine

   SimpleVector<double> g(nx);
   problem->getg(g);
   problem->hessian_multiplication(1.0, g, 0.5, *primals);
   double objective = g.dotProductWith(*primals);

   SimpleVector<double>& sx = (SimpleVector<double>&) *this->primals;
   SimpleVector<double>& sxlow = (SimpleVector<double>&) problem->xlowerBound();
   SimpleVector<double>& sixlow = (SimpleVector<double>&) problem->ixlowerBound();
   SimpleVector<double>& sxupp = (SimpleVector<double>&) problem->xupperBound();
   SimpleVector<double>& sixupp = (SimpleVector<double>&) problem->ixupperBound();
   SimpleVector<double>& sgamma = (SimpleVector<double>&) *this->primal_lower_bound_gap_dual;
   SimpleVector<double>& sphi = (SimpleVector<double>&) *this->primal_upper_bound_gap_dual;
   SimpleVector<double>& sy = (SimpleVector<double>&) *this->equality_duals;
   SimpleVector<double>& ss = (SimpleVector<double>&) *this->slacks;
   SimpleVector<double>& slambda = (SimpleVector<double>&) *this->slack_lower_bound_gap_dual;
   SimpleVector<double>& spi = (SimpleVector<double>&) *this->slack_upper_bound_gap_dual;
   SimpleVector<double>& sz = (SimpleVector<double>&) *this->inequality_duals;
   SimpleVector<double>& sclow = (SimpleVector<double>&) problem->slowerBound();
   SimpleVector<double>& siclow = (SimpleVector<double>&) problem->islowerBound();
   SimpleVector<double>& scupp = (SimpleVector<double>&) problem->supperBound();
   SimpleVector<double>& sicupp = (SimpleVector<double>&) problem->isupperBound();

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
      cclow = 0;
      ccupp = 0;
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

void Variables::print() {
   std::cout << " Complementary Variables = " << number_complementarity_pairs << std::endl;
   std::cout << "(Cannot tell you more at this level)" << std::endl;
}
