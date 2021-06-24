//
// Created by nils-christian on 15.06.21.
//

#include "BorderMod_Block.h"

BorderMod getChild(BorderMod& bordermod, unsigned int i) {
   BorderLinsys child = getChild(bordermod.border, i);
   return BorderMod(child, bordermod.multiplier);
}
