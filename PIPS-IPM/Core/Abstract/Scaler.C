/*
 * Scaler.C
 *
 *  Created on: 19.12.2017
 *      Authors: Daniel Rehfeldt, Svenja Uslu
 */

#include "Scaler.h"
#include "Problem.h"

Scaler::Scaler(Problem* prob, bool bitshifting, bool usesides)
: problem(prob), do_bitshifting(bitshifting), with_sides(usesides),
  dnorm_orig(prob->datanorm())
{
}
