/*
 * QpScaler.C
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

//#define PIPS_DEBUG
#include <algorithm>
#include "Scaler.hpp"
#include "PIPSIPMppOptions.h"
#include "DistributedVector.h"
#include "QP.hpp"
#include "Variables.h"
#include "Residuals.h"
#include "pipsdef.h"

Scaler::Scaler(const Problem& problem, bool bitshifting, bool usesides) : do_bitshifting(bitshifting), with_sides(usesides),
      dnorm_orig(problem.datanorm()),
      scaling_output{pipsipmpp_options::get_bool_parameter("SCALER_OUTPUT")} {
   A = problem.A;
   C = problem.C;
   obj = problem.g;
   bA = problem.bA;
   bux = problem.bux; // upper bound of x
   blx = problem.blx; // lower bound of x
   rhsC = problem.bu; // RHS of C
   lhsC = problem.bl; // LHS of C
   factor_objscale = 1.;
}

double Scaler::get_unscaled_objective(double objval) const {
   assert(vec_colscale != nullptr);
   assert(factor_objscale > 0.0);

   if (scaling_applied)
      return (objval / factor_objscale);
   else
      return objval;
}

Variables* Scaler::get_unscaled_variables(const Variables& variables) const {
   auto* unscaled_variables = new Variables(variables);
   unscale_variables(*unscaled_variables);
   return unscaled_variables;
};

Residuals* Scaler::get_unscaled_residuals(const Residuals& residuals) const {
   Residuals* unscaled_residuals = new Residuals(residuals);
   this->unscale_residuals(*unscaled_residuals);
   return unscaled_residuals;
};

void Scaler::unscale_variables(Variables& variables) const {
   if (!scaling_applied)
      return;

   // todo : Q
   assert(vec_colscale);
   assert(vec_rowscaleA);
   assert(vec_rowscaleC);

   variables.primals->componentMult(*vec_colscale);
   variables.slacks->componentDiv(*vec_rowscaleC);
   variables.equality_duals->componentMult(*vec_rowscaleA);
   variables.inequality_duals->componentMult(*vec_rowscaleC);
   variables.primal_lower_bound_gap->componentMult(*vec_colscale);
   variables.primal_lower_bound_gap_dual->componentDiv(*vec_colscale);
   variables.primal_upper_bound_gap->componentMult(*vec_colscale);
   variables.primal_upper_bound_gap_dual->componentDiv(*vec_colscale);
   variables.slack_lower_bound_gap->componentDiv(*vec_rowscaleC);
   variables.slack_lower_bound_gap_dual->componentMult(*vec_rowscaleC);
   variables.slack_upper_bound_gap->componentDiv(*vec_rowscaleC);
   variables.slack_upper_bound_gap_dual->componentMult(*vec_rowscaleC);
}

void Scaler::unscale_residuals(Residuals& residuals) const {
   if (!scaling_applied)
      return;

   assert(vec_colscale);
   assert(vec_rowscaleA);
   assert(vec_rowscaleC);

   residuals.lagrangian_gradient->componentDiv(*vec_colscale);
   residuals.rA->componentDiv(*vec_rowscaleA);
   residuals.rC->componentDiv(*vec_rowscaleC);
   residuals.rz->componentMult(*vec_rowscaleC);

   if (residuals.getNxlow() > 0)
      residuals.rv->componentMult(*vec_colscale);

   if (residuals.getNxupp() > 0)
      residuals.rw->componentMult(*vec_colscale);

   if (residuals.getMclow() > 0)
      residuals.rt->componentDiv(*vec_rowscaleC);

   if (residuals.getMcupp() > 0)
      residuals.ru->componentDiv(*vec_rowscaleC);
   // nothing to to for rgamma, rphi, rlambda, rpi;

   // gap is scaling resistant

   residuals.recompute_residual_norm();
}

Vector<double>* Scaler::get_primal_unscaled(const Vector<double>& primal_solution) const {
   assert(vec_colscale);
   Vector<double>* unscaledprimal = primal_solution.cloneFull();

   // unscale primal
   if (scaling_applied)
      unscaledprimal->componentMult(*vec_colscale);

   return unscaledprimal;
}

Vector<double>* Scaler::get_dual_eq_unscaled(const Vector<double>& dual_solution) const {
   assert(vec_rowscaleA);
   Vector<double>* unscaleddual = dual_solution.cloneFull();

   // unscale dual
   if (scaling_applied)
      unscaleddual->componentMult(*vec_rowscaleA);

   return unscaleddual;
}

Vector<double>* Scaler::get_dual_ineq_unscaled(const Vector<double>& dual_solution) const {
   assert(vec_rowscaleC);
   Vector<double>* unscaleddual = dual_solution.cloneFull();

   // unscale dual
   if (scaling_applied)
      unscaleddual->componentMult(*vec_rowscaleC);

   return unscaleddual;
}

Vector<double>* Scaler::get_dual_var_bounds_upp_unscaled(const Vector<double>& dual_solution) const {
   assert(vec_colscale);
   Vector<double>* unscaleddual = dual_solution.cloneFull();

   // unscale primal
   if (scaling_applied)
      unscaleddual->componentDiv(*vec_colscale);

   return unscaleddual;
}

Vector<double>* Scaler::get_dual_var_bounds_low_unscaled(const Vector<double>& dual_solution) const {
   assert(vec_colscale);
   Vector<double>* unscaleddual = dual_solution.cloneFull();

   // unscale primal
   if (scaling_applied)
      unscaleddual->componentDiv(*vec_colscale);

   return unscaleddual;
}


void Scaler::applyScaling() {
   PIPSdebugMessage("before scaling: \n "
                    "objnorm: %f \n Anorm:  %f \n Cnorm  %f \n bAnorm %f \n rhsCnorm %f \n lhsCnorm %f \n buxnorm %f \n blxnorm %f \n  ",
            obj->inf_norm(), A->inf_norm(), C->inf_norm(), bA->inf_norm(), rhsC->inf_norm(), lhsC->inf_norm(), bux->inf_norm(), blx->inf_norm());

   // todo scale Q
   doObjScaling();

   // scale A and rhs
   A->columnScale(*vec_colscale);
   A->rowScale(*vec_rowscaleA);
   bA->componentMult(*vec_rowscaleA);

   // scale C and lhs, rhs
   C->columnScale(*vec_colscale);
   C->rowScale(*vec_rowscaleC);
   rhsC->componentMult(*vec_rowscaleC);
   lhsC->componentMult(*vec_rowscaleC);

   // scale ub and lb of x
   bux->componentDiv(*vec_colscale);
   blx->componentDiv(*vec_colscale);

   scaling_applied = true;

   PIPSdebugMessage("after scaling: \n "
                    "objnorm: %f \n Anorm:  %f \n Cnorm  %f \n bAnorm %f \n rhsCnorm %f \n lhsCnorm %f \n buxnorm %f \n blxnorm %f \n  ",
            obj->inf_norm(), A->inf_norm(), C->inf_norm(), bA->inf_norm(), rhsC->inf_norm(), lhsC->inf_norm(), bux->inf_norm(), blx->inf_norm());
}

double Scaler::maxRowRatio(Vector<double>& maxvecA, Vector<double>& maxvecC, Vector<double>& minvecA, Vector<double>& minvecC,
      const Vector<double>* colScalevec) {
   A->getRowMinMaxVec(true, true, colScalevec, minvecA);
   A->getRowMinMaxVec(false, true, colScalevec, maxvecA);
   C->getRowMinMaxVec(true, true, colScalevec, minvecC);
   C->getRowMinMaxVec(false, true, colScalevec, maxvecC);

#ifndef NDEBUG

   if (!colScalevec) {
      int j = -1;
      double max = -std::numeric_limits<double>::max();

      maxvecA.max(max, j);
      assert(max < 0 || max == A->inf_norm());

      j = -1;
      max = -std::numeric_limits<double>::max();
      maxvecC.max(max, j);

      assert(max < 0 || max == C->inf_norm());
   }
#endif

   if (with_sides) {
      bA->absminVecUpdate(minvecA);
      rhsC->absminVecUpdate(minvecC);
      lhsC->absminVecUpdate(minvecC);

      bA->absmaxVecUpdate(maxvecA);
      rhsC->absmaxVecUpdate(maxvecC);
      lhsC->absmaxVecUpdate(maxvecC);
   }

   Vector<double>* const ratiovecA = maxvecA.clone();
   Vector<double>* const ratiovecC = maxvecC.clone();

   ratiovecA->copyFrom(maxvecA);
   ratiovecC->copyFrom(maxvecC);

   ratiovecA->divideSome(minvecA, minvecA);
   ratiovecC->divideSome(minvecC, minvecC);

   int i = -1;
   double maxratio = -std::numeric_limits<double>::infinity();
   ratiovecA->max(maxratio, i);

   PIPSdebugMessage("max row ratio A: %f \n", maxratio);

   double maxvalC = -std::numeric_limits<double>::infinity();
   ratiovecC->max(maxvalC, i);

   PIPSdebugMessage("max row ratio C: %f \n", maxvalC);

   if (maxvalC > maxratio)
      maxratio = maxvalC;

   delete ratiovecA;
   delete ratiovecC;

   return maxratio;
}

double Scaler::maxColRatio(Vector<double>& maxvec, Vector<double>& minvec, const Vector<double>* rowScaleVecA, const Vector<double>* rowScaleVecC) {
   A->getColMinMaxVec(true, true, rowScaleVecA, minvec);
   C->getColMinMaxVec(true, false, rowScaleVecC, minvec);

   A->getColMinMaxVec(false, true, rowScaleVecA, maxvec);
   C->getColMinMaxVec(false, false, rowScaleVecC, maxvec);

#ifndef NDEBUG
   if (!rowScaleVecA || !rowScaleVecC) {
      int j;
      double max;

      maxvec.max(max, j);

      assert(max < 0 || max == std::max(A->inf_norm(), C->inf_norm()));
   }
#endif

   Vector<double>* const ratiovec = maxvec.clone();

   ratiovec->copyFrom(maxvec);

   ratiovec->divideSome(minvec, minvec);

   int i;
   double maxratio;
   ratiovec->max(maxratio, i);
   assert(maxratio >= 0.0);

   PIPSdebugMessage("max column ratio: %f \n", maxratio);

   delete ratiovec;

   return maxratio;
}

void Scaler::scaleObjVector(double scaling_factor) {
   if (scaling_factor > 0.0)
      factor_objscale = 1.0 / scaling_factor;
   else
      factor_objscale = 1.0;

   if (do_bitshifting) {
      int exp;
      const double mantissa = std::frexp(factor_objscale, &exp);

      if (mantissa >= 0.75)
         factor_objscale = std::ldexp(0.5, exp + 1);
      else
         factor_objscale = std::ldexp(0.5, exp);
   }

   assert(obj);

   if (factor_objscale != 1.0)
      obj->scalarMult(factor_objscale);

   scaling_applied = true;
}

void Scaler::printRowColRatio() {
   if (scaling_output) {
      std::unique_ptr<DistributedVector<double>> xrowmaxA(dynamic_cast<DistributedVector<double>*>(bA->clone()));
      std::unique_ptr<DistributedVector<double>> xrowminA(dynamic_cast<DistributedVector<double>*>(bA->clone()));
      std::unique_ptr<DistributedVector<double>> xrowmaxC(dynamic_cast<DistributedVector<double>*>(rhsC->clone()));
      std::unique_ptr<DistributedVector<double>> xrowminC(dynamic_cast<DistributedVector<double>*>(rhsC->clone()));
      std::unique_ptr<DistributedVector<double>> xcolmax(dynamic_cast<DistributedVector<double>*>(bux->clone()));
      std::unique_ptr<DistributedVector<double>> xcolmin(dynamic_cast<DistributedVector<double>*>(bux->clone()));

      const double rowratio = maxRowRatio(*xrowmaxA, *xrowmaxC, *xrowminA, *xrowminC, nullptr);
      const double colratio = maxColRatio(*xcolmax, *xcolmin, nullptr, nullptr);

      if (PIPS_MPIgetRank() == 0) {
         printf("rowratio after scaling %f \n", rowratio);
         printf("colratio after scaling %f \n", colratio);
      }
   }
}

void Scaler::setScalingVecsToOne() {
   assert(vec_rowscaleA && vec_rowscaleC && vec_colscale);
   vec_rowscaleA->setToConstant(1.0);
   vec_rowscaleC->setToConstant(1.0);
   vec_colscale->setToConstant(1.0);
}
