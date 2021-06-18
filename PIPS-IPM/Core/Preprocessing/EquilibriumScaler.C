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

void EquilibriumScaler::doObjScaling() const {
   assert(vec_colscale != nullptr);

   obj->componentMult(*vec_colscale);

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

   const double rowratio = maxRowRatio(*vec_rowscaleA, *vec_rowscaleC, *rowminA, *rowminC, nullptr);
   const double colratio = maxColRatio(*vec_colscale, *colmin, nullptr, nullptr);

   const int myRank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   if (myRank == 0 && scaling_output) {
      printf("rowratio before scaling %f \n", rowratio);
      printf("colratio before scaling %f \n", colratio);
   }

   // column scaling first?
   if (colratio < rowratio && !with_sides) {
      invertAndRound(do_bitshifting, *vec_colscale);

      A->getRowMinMaxVec(false, true, vec_colscale.get(), *vec_rowscaleA);
      C->getRowMinMaxVec(false, true, vec_colscale.get(), *vec_rowscaleC);

      invertAndRound(do_bitshifting, *vec_rowscaleA);
      invertAndRound(do_bitshifting, *vec_rowscaleC);
   }
   else // row first
   {
      invertAndRound(do_bitshifting, *vec_rowscaleA);
      invertAndRound(do_bitshifting, *vec_rowscaleC);

      A->getColMinMaxVec(false, true, vec_rowscaleA.get(), *vec_colscale);
      C->getColMinMaxVec(false, false, vec_rowscaleC.get(), *vec_colscale);

      invertAndRound(do_bitshifting, *vec_colscale);
   }

   applyScaling();

   printRowColRatio();

   if (!scaling_applied) {
      setScalingVecsToOne();
#ifndef NDEBUG
      vec_rowscaleA->setToConstant(NAN);
      vec_rowscaleC->setToConstant(NAN);
      vec_colscale->setToConstant(NAN);
#endif
   }

   assert(A->inf_norm() <= 2.0 && C->inf_norm() <= 2.0);
}
