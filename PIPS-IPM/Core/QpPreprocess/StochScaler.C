/*
 * StochScaler.C
 *
 *  Created on: 18.07.2019
 *      Author: Nils-Christian
 */

#include "StochScaler.h"
#include "pipsdef.h"
#include "sVars.h"
#include "sResiduals.h"


StochScaler::StochScaler(Problem* prob, bool bitshifting)
	: QpScaler(prob, bitshifting)
{
}

Variables* StochScaler::getVariablesUnscaled(const Variables& vars) const
{
   Variables* s_vars = new sVars(dynamic_cast<const sVars&>(vars)); 
   assert(s_vars);
   assert( dynamic_cast<sVars*>(s_vars)->x);

   unscaleVariables(*s_vars);

   return s_vars;
};

Residuals* StochScaler::getResidualsUnscaled(const Residuals& resids) const
{
   sResiduals* s_resids = new sResiduals(dynamic_cast<const sResiduals&>(resids));
   assert(s_resids);

   unscaleResiduals(*s_resids);

   return s_resids;
};
