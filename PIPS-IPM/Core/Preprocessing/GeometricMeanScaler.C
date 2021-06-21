/*
 * GeometricMeanScaler.C
 *
 *  Created on: 16.07.2018
 *      Author: Svenja Uslu
 */

//#define PIPS_DEBUG
#include "GeometricMeanScaler.h"
#include "Vector.hpp"
#include "AbstractMatrix.h"
#include "ProblemFactory.h"
#include "pipsdef.h"
#include <memory>
#include <cmath>

static const double maxobjscale = 100.0;


GeometricMeanScaler::GeometricMeanScaler(const ProblemFactory& problem_factory, const Problem& problem, bool equiScaling, bool bitshifting) :
   Scaler(problem_factory, problem, bitshifting), minImpr{0.85}, goodEnough{500.}, maxIters{10}, equilibrate{equiScaling} {
   if (PIPS_MPIgetRank() == 0 && scaling_output)
      std::cout << "Creating GeometricMeanScaler... bitshifting=" << bitshifting << " equiscaling=" << equiScaling << "\n";
}

void GeometricMeanScaler::doObjScaling() const {
   assert(scaling_factors_columns != nullptr);

   obj->componentMult(*scaling_factors_columns);

   assert(factor_objscale == 1.0);

#if 0   // note: seems to deteriorate performance and stability
   const double absmax = obj->infnorm();
   double absmin = 0.0;

   obj->absminNonZero( absmin, pips_eps );

   // all elements of scaled obj smaller than pips_eps?
   if( PIPSisEQ(absmin, -1.0) )
   {
      int myRank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

      if( myRank == 0 )
         std::cout << "Almost zero objective after geometric scaling!" << "\n";
   }
   else
   {
      const double scaleFactor = std::sqrt(absmax * absmin);
      PIPSdebugMessage("Objective Scaling: absmin=%f, absmax=%f, scaleFactor=%f \n", absmin, absmax, scaleFactor);

      assert( scaleFactor >= 0.0 );
      scaleObjVector(scaleFactor);

      if( equilibrate )
      {
         const double absmax2 = obj->infnorm();
         assert(absmax2 >= 0);
         scaleObjVector(absmax2);
      }
   }
#endif
}

void GeometricMeanScaler::scale() {
   create_scaling_vectors();

   /* We want to do the direction with lower maximal ratio first,
    * since the absolute smallest value in the scaled matrix is bounded from below by
    * the inverse of the maximum ratio of the direction that is done first */
   std::unique_ptr<Vector<double>> rowminA{problem_factory.make_equalities_dual_vector()};
   std::unique_ptr<Vector<double>> rowminC{problem_factory.make_inequalities_dual_vector()};
   std::unique_ptr<Vector<double>> colmin{problem_factory.make_primal_vector()};

   const double rowratio = maxRowRatio(*scaling_factors_equalities, *scaling_factors_inequalitites, *rowminA, *rowminC, nullptr);
   const double colratio = maxColRatio(*scaling_factors_columns, *colmin, nullptr, nullptr);

   const int myRank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   if (myRank == 0 && scaling_output) {
      printf("rowratio before scaling %f \n", rowratio);
      printf("colratio before scaling %f \n", colratio);
   }

   double p0start, p1start;
   if (colratio < rowratio && !with_sides) {
      p0start = colratio;
      p1start = rowratio;
   }
   else {
      p0start = rowratio;
      p1start = colratio;
   }

   setScalingVecsToOne();

   bool geoscale = p1start > goodEnough;

   if (!geoscale) {
      if (myRank == 0 && scaling_output)
         std::cout << "No geometric scaling done, ratio already good enough.\n";
      if (!equilibrate)
         return;
      if (myRank == 0 && scaling_output)
         std::cout << "But will still perform equilibrium scaling.\n";
   }

   double p0 = 0.0;
   double p1 = 0.0;

   if (geoscale) {
      double p0prev = p0start;
      double p1prev = p1start;

      for (int i = 0; i < maxIters; i++) {
         // column scaling first?
         if (colratio < rowratio && !with_sides) {
            p0 = maxColRatio(*scaling_factors_columns, *colmin, scaling_factors_equalities.get(), scaling_factors_inequalitites.get());
            applyGeoMean(*scaling_factors_columns, *colmin);
            invertAndRound(do_bitshifting, *scaling_factors_columns);

            p1 = maxRowRatio(*scaling_factors_equalities, *scaling_factors_inequalitites, *rowminA, *rowminC, scaling_factors_columns.get());
            applyGeoMean(*scaling_factors_equalities, *rowminA);
            applyGeoMean(*scaling_factors_inequalitites, *rowminC);

            invertAndRound(do_bitshifting, *scaling_factors_equalities);
            invertAndRound(do_bitshifting, *scaling_factors_inequalitites);

            PIPSdebugMessage("Geometric Scaling round %d. colratio=%f, rowratio=%f \n", i, p0, p1);
         }
         else // row first
         {

            p0 = maxRowRatio(*scaling_factors_equalities, *scaling_factors_inequalitites, *rowminA, *rowminC, scaling_factors_columns.get());
            applyGeoMean(*scaling_factors_equalities, *rowminA);
            applyGeoMean(*scaling_factors_inequalitites, *rowminC);
            invertAndRound(do_bitshifting, *scaling_factors_equalities);
            invertAndRound(do_bitshifting, *scaling_factors_inequalitites);

            p1 = maxColRatio(*scaling_factors_columns, *colmin, scaling_factors_equalities.get(), scaling_factors_inequalitites.get());
            applyGeoMean(*scaling_factors_columns, *colmin);

            invertAndRound(do_bitshifting, *scaling_factors_columns);

            PIPSdebugMessage("Geometric Scaling round %d. colratio=%f, rowratio=%f \n", i, p1, p0);
         }
         // if ratio improvement is not good enough, then break:
         PIPSdebugMessage("p0=%f, p0prev=%f, p1=%f, p1prev=%f \n", p0, p0prev, p1, p1prev);
         if (p0 > minImpr * p0prev && p1 > minImpr * p1prev)
            break;

         p0prev = p0;
         p1prev = p1;
      }
      // perform geometric scaling only if there is enough (default 15%) improvement:
      geoscale = (p0 <= minImpr * p0start || p1 <= minImpr * p1start);
   }

   if (!geoscale && !equilibrate) {
      if (myRank == 0)
         std::cout << "No geometric scaling done, improvement was not good enough..." << "\n";
   }
   else {
      if (equilibrate) {
         if (!geoscale)
            setScalingVecsToOne();

         // equiScaling using the scaling vectors from GeoScaling:
         postEquiScale();
      }

      applyScaling();

#if 0
      double absmaxAll = bA->infnorm();
      absmaxAll = std::max(absmaxAll, bux->infnorm());
      absmaxAll = std::max(absmaxAll, blx->infnorm());
      absmaxAll = std::max(absmaxAll, rhsC->infnorm());
      absmaxAll = std::max(absmaxAll, lhsC->infnorm());

      std::cout << "absmax: " << absmaxAll <<  "\n\n\n\n\n" <<"\n";

      bA->scalarMult(1.0 / absmaxAll);
      bux->scalarMult(1.0 / absmaxAll);
      blx->scalarMult(1.0 / absmaxAll);
      rhsC->scalarMult(1.0 / absmaxAll);
      lhsC->scalarMult(1.0 / absmaxAll);
#endif

      printRowColRatio();

      if (equilibrate)
         assert(A->inf_norm() <= 2.0 && C->inf_norm() <= 2.0);
   }

   if (!scaling_applied) {
      setScalingVecsToOne();
#ifndef NDEBUG
      scaling_factors_equalities->setToConstant(NAN);
      scaling_factors_inequalitites->setToConstant(NAN);
      scaling_factors_columns->setToConstant(NAN);
#endif
   }
}

/** apply an approximation to the geometric mean to Vector maxvec:
 * Multiply maxvec and minvec componentwise and take the square root of the result.
 * Return result in maxvec.
 * */
void GeometricMeanScaler::applyGeoMean(Vector<double>& maxvec, const Vector<double>& minvec) {
   assert(maxvec.length() == minvec.length());

   maxvec.componentMult(minvec);
   maxvec.sqrt();
}

/** apply Equilibrium Scaling after having done Geometric Scaling.
 * The scaling vectors scaling_factors_equalities, scaling_factors_inequalitites and scaling_factors_columns should contain
 * the previously determined scaling factors.
 */
void GeometricMeanScaler::postEquiScale() {
   assert(scaling_factors_equalities && scaling_factors_inequalitites && scaling_factors_columns);

   std::unique_ptr<Vector<double>> rowmaxA{problem_factory.make_equalities_dual_vector()};
   std::unique_ptr<Vector<double>> rowminA{problem_factory.make_equalities_dual_vector()};
   std::unique_ptr<Vector<double>> rowmaxC{problem_factory.make_inequalities_dual_vector()};
   std::unique_ptr<Vector<double>> rowminC{problem_factory.make_inequalities_dual_vector()};
   std::unique_ptr<Vector<double>> colmax{problem_factory.make_primal_vector()};
   std::unique_ptr<Vector<double>> colmin{problem_factory.make_primal_vector()};

   const double rowratio = maxRowRatio(*rowmaxA, *rowmaxC, *rowminA, *rowminC, scaling_factors_columns.get());
   const double colratio = maxColRatio(*colmax, *colmin, scaling_factors_equalities.get(), scaling_factors_inequalitites.get());

   PIPSdebugMessage("rowratio before Post-EquiScale %f \n", rowratio);
   PIPSdebugMessage("colratio before Post-EquiScale %f \n", colratio);

   std::swap(scaling_factors_columns, colmax);
   std::swap(scaling_factors_equalities, rowmaxA);
   std::swap(scaling_factors_inequalitites, rowmaxC);

   // column scaling first?
   if (colratio < rowratio && !with_sides) {
      invertAndRound(do_bitshifting, *scaling_factors_columns);

      A->getRowMinMaxVec(false, true, scaling_factors_columns.get(), *scaling_factors_equalities);
      C->getRowMinMaxVec(false, true, scaling_factors_columns.get(), *scaling_factors_inequalitites);

      invertAndRound(do_bitshifting, *scaling_factors_equalities);
      invertAndRound(do_bitshifting, *scaling_factors_inequalitites);
   }
   else // row first
   {
      invertAndRound(do_bitshifting, *scaling_factors_equalities);
      invertAndRound(do_bitshifting, *scaling_factors_inequalitites);

      A->getColMinMaxVec(false, true, scaling_factors_equalities.get(), *scaling_factors_columns);
      C->getColMinMaxVec(false, false, scaling_factors_inequalitites.get(), *scaling_factors_columns);

      invertAndRound(do_bitshifting, *scaling_factors_columns);
   }
}


