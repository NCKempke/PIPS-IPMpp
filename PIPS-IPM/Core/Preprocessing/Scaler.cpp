/*
 * QpScaler.C
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

//#define PIPS_DEBUG
#include "Scaler.hpp"

#include <algorithm>

#include "PIPSIPMppOptions.h"
#include "ProblemFactory.h"
#include "Variables.h"
#include "Problem.hpp"
#include "Residuals.h"
#include "pipsdef.h"

Scaler::Scaler(const ProblemFactory& problem_factory_, const Problem& problem, bool bitshifting, bool usesides) : problem_factory{problem_factory_},
   A{problem.equality_jacobian}, C{problem.inequality_jacobian}, obj{problem.objective_gradient}, bA{problem.equality_rhs}, bux{problem.primal_upper_bounds}, blx{problem.primal_lower_bounds},
   rhsC{problem.inequality_upper_bounds}, lhsC{problem.inequality_lower_bounds}, factor_objscale{1.},
   dnorm_orig(problem.datanorm()), do_bitshifting(bitshifting), with_sides(usesides),
   scaling_output{pipsipmpp_options::get_bool_parameter("SCALER_OUTPUT")} {}

double Scaler::get_unscaled_objective(double objval) const {
   assert(vec_colscale);
   assert(factor_objscale > 0.0);

   if (scaling_applied)
      return (objval / factor_objscale);
   else
      return objval;
}

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
   residuals.equality_residuals->componentDiv(*vec_rowscaleA);
   residuals.inequality_residuals->componentDiv(*vec_rowscaleC);
   residuals.inequality_dual_residuals->componentMult(*vec_rowscaleC);

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

void Scaler::invertAndRound(bool round, Vector<double>& vector) {
   vector.safe_invert(1.0);
   if (round)
      vector.roundToPow2();
}

void Scaler::create_scaling_vectors() {
   assert(!vec_rowscaleA && !vec_rowscaleC && !vec_colscale);

   vec_rowscaleA = problem_factory.make_equalities_dual_vector();
   vec_rowscaleC = problem_factory.make_inequalities_dual_vector();
   vec_colscale = problem_factory.make_primal_vector();
}

void Scaler::applyScaling() const {
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
      const Vector<double>* colScalevec) const {
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

   std::unique_ptr<Vector<double>> ratiovecA{maxvecA.clone()};
   std::unique_ptr<Vector<double>> ratiovecC{maxvecC.clone()};

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

   return maxratio;
}

double Scaler::maxColRatio(Vector<double>& maxvec, Vector<double>& minvec, const Vector<double>* rowScaleVecA, const Vector<double>* rowScaleVecC) const {
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

   std::unique_ptr<Vector<double>> ratiovec{maxvec.clone()};

   ratiovec->copyFrom(maxvec);
   ratiovec->divideSome(minvec, minvec);

   int i;
   double maxratio;
   ratiovec->max(maxratio, i);
   assert(maxratio >= 0.0);

   PIPSdebugMessage("max column ratio: %f \n", maxratio);

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

void Scaler::printRowColRatio() const {
   if (scaling_output) {
      std::unique_ptr<Vector<double>> xrowmaxA(problem_factory.make_equalities_dual_vector());
      std::unique_ptr<Vector<double>> xrowminA(problem_factory.make_equalities_dual_vector());
      std::unique_ptr<Vector<double>> xrowmaxC(problem_factory.make_inequalities_dual_vector());
      std::unique_ptr<Vector<double>> xrowminC(problem_factory.make_inequalities_dual_vector());
      std::unique_ptr<Vector<double>> xcolmax(problem_factory.make_primal_vector());
      std::unique_ptr<Vector<double>> xcolmin(problem_factory.make_primal_vector());

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
