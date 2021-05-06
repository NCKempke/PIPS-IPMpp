/*
 * GeoStochScaler.h
 *
 *  Created on: 16.07.2018
 *      Author: Svenja Uslu
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_GEOSTOCHSCALER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_GEOSTOCHSCALER_H_

#include "StochScaler.h"
#include "DistributedVector.h"

class Problem;


/**  * @defgroup QpPreprocess
 *
 * Geometric scaler
 * @{
 */

/**
 * Derived class for Geometric scaler.
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

   GeoStochScaler(Problem* prob, bool equiScaling, bool bitshifting = false);
   virtual ~GeoStochScaler() = default;

   /** scale */
   void scale() override;
};

//@}

#endif /* PIPS_IPM_CORE_QPPREPROCESS_GEOSTOCHSCALER_H_ */