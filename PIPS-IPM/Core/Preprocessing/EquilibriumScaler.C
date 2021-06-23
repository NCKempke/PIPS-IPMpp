/*
 * EquilibriumScaler.C
 *
 *  Created on: 20.12.2017
 *      Author: bzfrehfe
 */

//#define PIPS_DEBUG
#include "EquilibriumScaler.h"
#include "Vector.hpp"
#include "ProblemFactory.h"
#include "AbstractMatrix.h"

#include "pipsdef.h"
#include <cmath>
#include <memory>

EquilibriumScaler::EquilibriumScaler(const ProblemFactory& problem_factory, const Problem& problem, bool bitshifting) :
   Scaler(problem_factory, problem, bitshifting) {
   if (PIPS_MPIgetRank() == 0 && scaling_output)
      std::cout << "Creating EquilibriumScaler...\n";
}

void EquilibriumScaler::scale_objective() const {
   assert(scaling_factors_columns);

   obj->componentMult(*scaling_factors_columns);

#if 0 // note: seems to deteriorate performance and stability
   const double absmax = obj->infnorm();
   assert(absmax >= 0);

   scaleObjVector(absmax);
#endif
}

// todo scale Q
void EquilibriumScaler::scale() {
   create_scaling_vectors();

   /* We want to do the direction with lower maximal ratio first,
    * since the absolute smallest value in the scaled matrix is bounded from below by
    * the inverse of the maximum ratio of the direction that is done first
    */
   std::unique_ptr<Vector<double>> rowminA{problem_factory.make_equalities_dual_vector()};
   std::unique_ptr<Vector<double>> rowminC{problem_factory.make_inequalities_dual_vector()};
   std::unique_ptr<Vector<double>> colmin{problem_factory.make_primal_vector()};

   const double rowratio = maxRowRatio(*scaling_factors_equalities, *scaling_factors_inequalities, *rowminA, *rowminC, nullptr);
   const double colratio = maxColRatio(*scaling_factors_columns, *colmin, nullptr, nullptr);

   const int myRank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   if (myRank == 0 && scaling_output) {
      printf("rowratio before scaling %f \n", rowratio);
      printf("colratio before scaling %f \n", colratio);
   }

   // column scaling first?
   if (colratio < rowratio && !with_sides) {
      invertAndRound(do_bitshifting, *scaling_factors_columns);

      A->getRowMinMaxVec(false, true, scaling_factors_columns.get(), *scaling_factors_equalities);
      C->getRowMinMaxVec(false, true, scaling_factors_columns.get(), *scaling_factors_inequalities);

      invertAndRound(do_bitshifting, *scaling_factors_equalities);
      invertAndRound(do_bitshifting, *scaling_factors_inequalities);
   }
   else // row first
   {
      invertAndRound(do_bitshifting, *scaling_factors_equalities);
      invertAndRound(do_bitshifting, *scaling_factors_inequalities);

      A->getColMinMaxVec(false, true, scaling_factors_equalities.get(), *scaling_factors_columns);
      C->getColMinMaxVec(false, false, scaling_factors_inequalities.get(), *scaling_factors_columns);

      invertAndRound(do_bitshifting, *scaling_factors_columns);
   }

   applyScaling();

   printRowColRatio();

   if (!scaling_applied) {
      setScalingVecsToOne();
#ifndef NDEBUG
      scaling_factors_equalities->setToConstant(NAN);
      scaling_factors_inequalities->setToConstant(NAN);
      scaling_factors_columns->setToConstant(NAN);
#endif
   }

   assert(A->inf_norm() <= 2.0 && C->inf_norm() <= 2.0);
}
