/*
 * EquiStochScaler.C
 *
 *  Created on: 20.12.2017
 *      Author: bzfrehfe
 */

//#define PIPS_DEBUG
#include "EquiStochScaler.h"
#include "StochVector.h"
#include "pipsdef.h"
#include "pipsport.h"

#include <cmath>
#include <memory>

EquiStochScaler::EquiStochScaler(Data* prob, bool bitshifting)
  : StochScaler(prob, bitshifting)
{
   if( PIPS_MPIgetRank() == 0 && scaling_output )
      std::cout << "Creating EquiStochScaler...\n";
}

void EquiStochScaler::doObjScaling()
{
   assert(vec_colscale != nullptr);

   obj->componentMult(*vec_colscale);

#if 0 // note: seems to deteriorate performance and stability
   const double absmax = obj->infnorm();
   assert(absmax >= 0);

   scaleObjVector(absmax);
#endif
}

// todo scale Q
void EquiStochScaler::scale()
{
   assert( !vec_rowscaleA && !vec_rowscaleC && !vec_colscale );

   /* We want to do the direction with lower maximal ratio first,
    * since the absolute smallest value in the scaled matrix is bounded from below by
    * the inverse of the maximum ratio of the direction that is done first */
   vec_rowscaleA.reset( dynamic_cast<StochVector*>(bA->clone()) );
   std::unique_ptr<StochVector> rowminA{ dynamic_cast<StochVector*>(bA->clone()) };
   vec_rowscaleC.reset( dynamic_cast<StochVector*>(rhsC->clone()) );
   std::unique_ptr<StochVector> rowminC{ dynamic_cast<StochVector*>(rhsC->clone()) };
   vec_colscale.reset( dynamic_cast<StochVector*>(bux->clone()) );
   std::unique_ptr<StochVector> colmin{ dynamic_cast<StochVector*>(bux->clone()) };

   const double rowratio = maxRowRatio(*vec_rowscaleA, *vec_rowscaleC, *rowminA, *rowminC, nullptr);
   const double colratio = maxColRatio(*vec_colscale, *colmin, nullptr, nullptr);

   const int myRank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   if( myRank == 0 && scaling_output )
   {
      printf("rowratio before scaling %f \n", rowratio);
      printf("colratio before scaling %f \n", colratio);
   }

   // column scaling first?
   if( colratio < rowratio && !with_sides )
   {
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

   if( !scaling_applied )
   {
      setScalingVecsToOne();
#ifndef NDEBUG
      vec_rowscaleA->setToConstant(NAN);
      vec_rowscaleC->setToConstant(NAN);
      vec_colscale->setToConstant(NAN);
#endif
   }

   assert(A->abmaxnorm() <= 2.0 && C->abmaxnorm() <= 2.0);
}
