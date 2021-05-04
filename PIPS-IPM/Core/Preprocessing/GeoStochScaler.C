/*
 * GeoStochScaler.C
 *
 *  Created on: 16.07.2018
 *      Author: Svenja Uslu
 */

//#define PIPS_DEBUG
#include "GeoStochScaler.h"
#include "pipsdef.h"
#include <memory>
#include <cmath>

static const double maxobjscale = 100.0;


GeoStochScaler::GeoStochScaler(Problem* prob, bool equiScaling, bool bitshifting) : StochScaler(prob, bitshifting) {
   if (PIPS_MPIgetRank() == 0 && scaling_output)
      std::cout << "Creating GeoStochScaler... bitshifting=" << bitshifting << " equiscaling=" << equiScaling << "\n";
   equilibrate = equiScaling;

   // todo: adjust parameters
   maxIters = 10;
   minImpr = 0.85;
   goodEnough = 500;
}

void GeoStochScaler::doObjScaling() {
   assert(vec_colscale != nullptr);

   obj->componentMult(*vec_colscale);

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

void GeoStochScaler::scale() {
   assert(!vec_rowscaleA && !vec_rowscaleC && !vec_colscale);

   /* We want to do the direction with lower maximal ratio first,
    * since the absolute smallest value in the scaled matrix is bounded from below by
    * the inverse of the maximum ratio of the direction that is done first */
   vec_rowscaleA.reset(dynamic_cast<DistributedVector<double>*>(bA->clone()));
   std::unique_ptr<DistributedVector<double>> rowminA{dynamic_cast<DistributedVector<double>*>(bA->clone())};
   vec_rowscaleC.reset(dynamic_cast<DistributedVector<double>*>(rhsC->clone()));
   std::unique_ptr<DistributedVector<double>> rowminC{dynamic_cast<DistributedVector<double>*>(rhsC->clone())};
   vec_colscale.reset(dynamic_cast<DistributedVector<double>*>(bux->clone()));
   std::unique_ptr<DistributedVector<double>> colmin{dynamic_cast<DistributedVector<double>*>(bux->clone())};

   const double rowratio = maxRowRatio(*vec_rowscaleA, *vec_rowscaleC, *rowminA, *rowminC, nullptr);
   const double colratio = maxColRatio(*vec_colscale, *colmin, nullptr, nullptr);

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
            p0 = maxColRatio(*vec_colscale, *colmin, vec_rowscaleA.get(), vec_rowscaleC.get());
            applyGeoMean(*vec_colscale, *colmin);
            invertAndRound(do_bitshifting, *vec_colscale);

            p1 = maxRowRatio(*vec_rowscaleA, *vec_rowscaleC, *rowminA, *rowminC, vec_colscale.get());
            applyGeoMean(*vec_rowscaleA, *rowminA);
            applyGeoMean(*vec_rowscaleC, *rowminC);

            invertAndRound(do_bitshifting, *vec_rowscaleA);
            invertAndRound(do_bitshifting, *vec_rowscaleC);

            PIPSdebugMessage("Geometric Scaling round %d. colratio=%f, rowratio=%f \n", i, p0, p1);
         }
         else // row first
         {

            p0 = maxRowRatio(*vec_rowscaleA, *vec_rowscaleC, *rowminA, *rowminC, vec_colscale.get());
            applyGeoMean(*vec_rowscaleA, *rowminA);
            applyGeoMean(*vec_rowscaleC, *rowminC);
            invertAndRound(do_bitshifting, *vec_rowscaleA);
            invertAndRound(do_bitshifting, *vec_rowscaleC);

            p1 = maxColRatio(*vec_colscale, *colmin, vec_rowscaleA.get(), vec_rowscaleC.get());
            applyGeoMean(*vec_colscale, *colmin);

            invertAndRound(do_bitshifting, *vec_colscale);

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
         assert(A->abmaxnorm() <= 2.0 && C->abmaxnorm() <= 2.0);
   }

   if (!scaling_applied) {
      setScalingVecsToOne();
#ifndef NDEBUG
      vec_rowscaleA->setToConstant(NAN);
      vec_rowscaleC->setToConstant(NAN);
      vec_colscale->setToConstant(NAN);
#endif
   }
}

/** apply an approximation to the geometric mean to Vector maxvec:
 * Multiply maxvec and minvec componentwise and take the square root of the result.
 * Return result in maxvec.
 * */
void GeoStochScaler::applyGeoMean(Vector<double>& maxvec, const Vector<double>& minvec) {
   assert(maxvec.length() == minvec.length());

   maxvec.componentMult(minvec);
   maxvec.applySqrt();
}

/** apply Equilibrium Scaling after having done Geometric Scaling.
 * The scaling vectors vec_rowscaleA, vec_rowscaleC and vec_colscale should contain
 * the previously determined scaling factors.
 */
void GeoStochScaler::postEquiScale() {
   assert(vec_rowscaleA != nullptr && vec_rowscaleC != nullptr && vec_colscale != nullptr);

   DistributedVector<double>* rowmaxA = dynamic_cast<DistributedVector<double>*>(bA->clone());
   std::unique_ptr<DistributedVector<double>> rowminA{dynamic_cast<DistributedVector<double>*>(bA->clone())};
   DistributedVector<double>* rowmaxC = dynamic_cast<DistributedVector<double>*>(rhsC->clone());
   std::unique_ptr<DistributedVector<double>> rowminC{dynamic_cast<DistributedVector<double>*>(rhsC->clone())};
   DistributedVector<double>* colmax = dynamic_cast<DistributedVector<double>*>(bux->clone());
   std::unique_ptr<DistributedVector<double>> colmin{dynamic_cast<DistributedVector<double>*>(bux->clone())};

   const double rowratio = maxRowRatio(*rowmaxA, *rowmaxC, *rowminA, *rowminC, vec_colscale.get());
   const double colratio = maxColRatio(*colmax, *colmin, vec_rowscaleA.get(), vec_rowscaleC.get());

   PIPSdebugMessage("rowratio before Post-EquiScale %f \n", rowratio);
   PIPSdebugMessage("colratio before Post-EquiScale %f \n", colratio);

   vec_colscale.reset(colmax);
   vec_rowscaleA.reset(rowmaxA);
   vec_rowscaleC.reset(rowmaxC);

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

}


