//
// Created by charlie on 26.04.21.
//


#include "ProblemFormulation.h"
#include "QP.hpp"
#include "Variables.h"
#include "Residuals.h"
#include "LinearAlgebraPackage.h"


ProblemFormulation::ProblemFormulation(int nx_, int my_, int mz_) : nx(nx_), my(my_), mz(mz_) {
}

OoqpVector* ProblemFormulation::makePrimalVector() const {
   assert(la);
   return la->newVector(nx);
}

OoqpVector* ProblemFormulation::makeDualYVector() const {
   assert(la);
   return la->newVector(my);

}

OoqpVector* ProblemFormulation::makeDualZVector() const {
   assert(la);
   return la->newVector(mz);
}

OoqpVector* ProblemFormulation::makeRhs() const {
   assert(la);
   return la->newVector(nx + my + mz);
}

