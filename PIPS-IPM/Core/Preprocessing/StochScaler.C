/*
 * StochScaler.C
 *
 *  Created on: 18.07.2019
 *      Author: Nils-Christian
 */

#include "StochScaler.h"
#include "pipsdef.h"
#include "DistributedVariables.h"
#include "DistributedResiduals.hpp"


StochScaler::StochScaler(const Problem& problem, bool bitshifting) : QpScaler(problem, bitshifting) {
}

Variables* StochScaler::get_unscaled_variables(const Variables& variables) const {
   Variables* s_vars = new DistributedVariables(dynamic_cast<const DistributedVariables&>(variables));
   assert(s_vars);
   assert(dynamic_cast<DistributedVariables*>(s_vars)->primals);

   unscale_variables(*s_vars);

   return s_vars;
};

Residuals* StochScaler::get_unscaled_residuals(const Residuals& residuals) const {
   DistributedResiduals* s_resids = new DistributedResiduals(dynamic_cast<const DistributedResiduals&>(residuals));
   assert(s_resids);

   unscale_residuals(*s_resids);

   return s_resids;
};
