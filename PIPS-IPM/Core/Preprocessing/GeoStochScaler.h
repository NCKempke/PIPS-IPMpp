/*
 * GeoStochScaler.h
 *
 *  Created on: 16.07.2018
 *      Author: Svenja Uslu
 */

#ifndef GEOSTOCHSCALER_H
#define GEOSTOCHSCALER_H

#include "Scaler.hpp"

class Problem;

/**  * @defgroup QpPreprocess
 *
 * Geometric scaler
 * @{
 */
class GeoStochScaler : public Scaler {

private:
   const double minImpr;
   const double goodEnough;

   const int maxIters;
   const bool equilibrate;

protected:
   void doObjScaling() const override;

   static void applyGeoMean(Vector<double>& maxvec, const Vector<double>& minvec);
   void postEquiScale();

public:
   GeoStochScaler(const ProblemFactory& problem_factory, const Problem& problem, bool equiScaling, bool bitshifting = false);
   ~GeoStochScaler() override = default;

   /** scale */
   void scale() override;
};

#endif /* GEOSTOCHSCALER_H */
