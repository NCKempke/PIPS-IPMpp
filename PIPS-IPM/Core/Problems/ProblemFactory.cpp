//
// Created by charlie on 26.04.21.
//


#include "ProblemFactory.h"
#include "QP.hpp"
#include "Variables.h"
#include "Residuals.h"
#include "LinearAlgebraPackage.h"


ProblemFactory::ProblemFactory(int nx_, int my_, int mz_) : nx(nx_), my(my_), mz(mz_) {
}

Vector<double>* ProblemFactory::make_primal_vector() const {
   assert(la);
   return la->newVector(nx);
}

Vector<double>* ProblemFactory::make_equalities_dual_vector() const {
   assert(la);
   return la->newVector(my);

}

Vector<double>* ProblemFactory::make_inequalities_dual_vector() const {
   assert(la);
   return la->newVector(mz);
}

Vector<double>* ProblemFactory::make_right_hand_side() const {
   assert(la);
   return la->newVector(nx + my + mz);
}

