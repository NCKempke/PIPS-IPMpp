//
// Created by nils-christian on 15.06.21.
//

#ifndef PIPSIPMPP_BORDERMOD_BLOCK_H
#define PIPSIPMPP_BORDERMOD_BLOCK_H

#include "RACFG_BLOCK.h"
#include "DenseMatrix.h"

/* stores a border and the corresponding multiplier - needed for hierarchical Schur Complement computation and Ltsolves */
template<typename T>
struct BorderMod_Block {
public:
   BorderLinsys border;
   const T& multiplier;

   BorderMod_Block(BorderLinsys& border_, const T& multiplier) : border{border_}, multiplier{multiplier} {};
};

using BorderMod = BorderMod_Block<DenseMatrix>;

BorderMod getChild(BorderMod& bordermod, unsigned int i);


#endif //PIPSIPMPP_BORDERMOD_BLOCK_H
