/*
 * EquiStochScaler.h
 *
 *  Created on: 20.12.2017
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_EQUISTOCHSCALER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_EQUISTOCHSCALER_H_

#include "StochScaler.h"

class Problem;

/**  * @defgroup QpPreprocess
 *
 * scaler
 * @{
 */

class EquiStochScaler : public StochScaler {
protected:
   void doObjScaling() override;

public:
   EquiStochScaler(const Problem& problem, bool bitshifting = true);
   ~EquiStochScaler() override = default;

   /** scale */
   void scale() override;
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_EQUISTOCHSCALER_H_ */
