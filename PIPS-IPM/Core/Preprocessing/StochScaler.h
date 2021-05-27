/*
 * StochScaler.h
 *
 *  Created on: 18.07.2019
 *      Author: Nils-Christian Kempke
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHSCALER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHSCALER_H_

#include "Scaler.hpp"

class Problem;

class StochScaler : public Scaler {
public:
   StochScaler(const Problem& problem, bool bitshifting);
   ~StochScaler() override = default;

   Variables* get_unscaled_variables(const Variables& variables) const override;
   Residuals* get_unscaled_residuals(const Residuals& residuals) const override;
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHSCALER_H_ */
