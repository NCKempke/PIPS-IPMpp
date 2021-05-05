/*
 * Scaler.C
 *
 *  Created on: 19.12.2017
 *      Authors: Daniel Rehfeldt, Svenja Uslu
 */

#include "Scaler.h"
#include "Problem.h"

Scaler::Scaler(Problem* problem, bool bitshifting, bool usesides) : problem(problem), do_bitshifting(bitshifting), with_sides(usesides),
      dnorm_orig(problem->datanorm()) {
}
