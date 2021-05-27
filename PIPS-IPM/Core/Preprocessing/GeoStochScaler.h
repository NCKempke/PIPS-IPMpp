/*
 * GeoStochScaler.h
 *
 *  Created on: 16.07.2018
 *      Author: Svenja Uslu
 */

#ifndef GEOSTOCHSCALER_H
#define GEOSTOCHSCALER_H

#include "StochScaler.h"
#include "DistributedVector.h"

class Problem;

/**  * @defgroup QpPreprocess
 *
 * Geometric scaler
 * @{
 */
class GeoStochScaler : public StochScaler {
   bool equilibrate;
   int maxIters;
   double minImpr;
   double goodEnough;

protected:
   void doObjScaling() override;

   void applyGeoMean(Vector<double>& maxvec, const Vector<double>& minvec);
   void postEquiScale();
public:

   GeoStochScaler(const Problem& problem, bool equiScaling, bool bitshifting = false);
   virtual ~GeoStochScaler() = default;

   /** scale */
   void scale() override;
};

#endif /* GEOSTOCHSCALER_H */
