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
   assert(scaling_factors_columns);
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
   assert(scaling_factors_columns);
   assert(scaling_factors_equalities);
   assert(scaling_factors_inequalities);

   variables.primals->componentMult(*scaling_factors_columns);
   variables.slacks->componentDiv(*scaling_factors_inequalities);
   variables.equality_duals->componentMult(*scaling_factors_equalities);
   variables.inequality_duals->componentMult(*scaling_factors_inequalities);
   variables.primal_lower_bound_gap->componentMult(*scaling_factors_columns);
   variables.primal_lower_bound_gap_dual->componentDiv(*scaling_factors_columns);
   variables.primal_upper_bound_gap->componentMult(*scaling_factors_columns);
   variables.primal_upper_bound_gap_dual->componentDiv(*scaling_factors_columns);
   variables.slack_lower_bound_gap->componentDiv(*scaling_factors_inequalities);
   variables.slack_lower_bound_gap_dual->componentMult(*scaling_factors_inequalities);
   variables.slack_upper_bound_gap->componentDiv(*scaling_factors_inequalities);
   variables.slack_upper_bound_gap_dual->componentMult(*scaling_factors_inequalities);
}

void Scaler::unscale_residuals(Residuals& residuals) const {
   if (!scaling_applied)
      return;

   assert(scaling_factors_columns);
   assert(scaling_factors_equalities);
   assert(scaling_factors_inequalities);

   residuals.lagrangian_gradient->componentDiv(*scaling_factors_columns);
   residuals.equality_residuals->componentDiv(*scaling_factors_equalities);
   residuals.inequality_residuals->componentDiv(*scaling_factors_inequalities);
   residuals.inequality_dual_residuals->componentMult(*scaling_factors_inequalities);

   if (residuals.getNxlow() > 0)
      residuals.rv->componentMult(*scaling_factors_columns);

   if (residuals.getNxupp() > 0)
      residuals.rw->componentMult(*scaling_factors_columns);

   if (residuals.getMclow() > 0)
      residuals.rt->componentDiv(*scaling_factors_inequalities);

   if (residuals.getMcupp() > 0)
      residuals.ru->componentDiv(*scaling_factors_inequalities);
   // nothing to to for rgamma, rphi, rlambda, rpi;

   // gap is scaling resistant

   residuals.recompute_residual_norm();
}

Vector<double>* Scaler::get_primal_unscaled(const Vector<double>& primal_solution) const {
   assert(scaling_factors_columns);
   Vector<double>* unscaledprimal = primal_solution.cloneFull();

   // unscale primal
   if (scaling_applied)
      unscaledprimal->componentMult(*scaling_factors_columns);

   return unscaledprimal;
}

Vector<double>* Scaler::get_dual_eq_unscaled(const Vector<double>& dual_solution) const {
   assert(scaling_factors_equalities);
   Vector<double>* unscaleddual = dual_solution.cloneFull();

   // unscale dual
   if (scaling_applied)
      unscaleddual->componentMult(*scaling_factors_equalities);

   return unscaleddual;
}

Vector<double>* Scaler::get_dual_ineq_unscaled(const Vector<double>& dual_solution) const {
   assert(scaling_factors_inequalities);
   Vector<double>* unscaleddual = dual_solution.cloneFull();

   // unscale dual
   if (scaling_applied)
      unscaleddual->componentMult(*scaling_factors_inequalities);

   return unscaleddual;
}

Vector<double>* Scaler::get_dual_var_bounds_upp_unscaled(const Vector<double>& dual_solution) const {
   assert(scaling_factors_columns);
   Vector<double>* unscaleddual = dual_solution.cloneFull();

   // unscale primal
   if (scaling_applied)
      unscaleddual->componentDiv(*scaling_factors_columns);

   return unscaleddual;
}

Vector<double>* Scaler::get_dual_var_bounds_low_unscaled(const Vector<double>& dual_solution) const {
   assert(scaling_factors_columns);
   Vector<double>* unscaleddual = dual_solution.cloneFull();

   // unscale primal
   if (scaling_applied)
      unscaleddual->componentDiv(*scaling_factors_columns);

   return unscaleddual;
}

void Scaler::invertAndRound(bool round, Vector<double>& vector) {
   vector.safe_invert(1.0);
   if (round)
      vector.roundToPow2();
}

void Scaler::create_scaling_vectors() {
   assert(!scaling_factors_equalities && !scaling_factors_inequalities && !scaling_factors_columns);
   std::tie(scaling_factors_equalities, scaling_factors_inequalities, scaling_factors_columns) = create_primal_dual_vector_triplet();
}

PrimalDualTriplet Scaler::create_primal_dual_vector_triplet() const {
   return {problem_factory.make_equalities_dual_vector(), problem_factory.make_inequalities_dual_vector(), problem_factory.make_primal_vector()};
}

void Scaler::applyScaling() const {
   PIPSdebugMessage("before scaling: \n "
                    "objnorm: %f \n Anorm:  %f \n Cnorm  %f \n bAnorm %f \n rhsCnorm %f \n lhsCnorm %f \n buxnorm %f \n blxnorm %f \n  ",
            obj->inf_norm(), A->inf_norm(), C->inf_norm(), bA->inf_norm(), rhsC->inf_norm(), lhsC->inf_norm(), bux->inf_norm(), blx->inf_norm());

   // todo scale Q
   scale_objective();

   // scale A and rhs
   A->columnScale(*scaling_factors_columns);
   A->rowScale(*scaling_factors_equalities);
   bA->componentMult(*scaling_factors_equalities);

   // scale C and lhs, rhs
   C->columnScale(*scaling_factors_columns);
   C->rowScale(*scaling_factors_inequalities);
   rhsC->componentMult(*scaling_factors_inequalities);
   lhsC->componentMult(*scaling_factors_inequalities);

   // scale ub and lb of x
   bux->componentDiv(*scaling_factors_columns);
   blx->componentDiv(*scaling_factors_columns);

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
   assert(scaling_factors_equalities && scaling_factors_inequalities && scaling_factors_columns);
   scaling_factors_equalities->setToConstant(1.0);
   scaling_factors_inequalities->setToConstant(1.0);
   scaling_factors_columns->setToConstant(1.0);
}
