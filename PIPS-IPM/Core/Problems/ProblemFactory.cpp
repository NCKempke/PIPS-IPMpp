//
// Created by charlie on 26.04.21.
//


#include "ProblemFactory.h"
#include "SimpleVector.h"

ProblemFactory::ProblemFactory(int nx_, int my_, int mz_) : nx(nx_), my(my_), mz(mz_) {
}

Vector<double>* ProblemFactory::make_primal_vector() const {
   return new SimpleVector<double>(nx);
}

Vector<double>* ProblemFactory::make_equalities_dual_vector() const {
   return new SimpleVector<double>(my);
}

Vector<double>* ProblemFactory::make_inequalities_dual_vector() const {
   return new SimpleVector<double>(mz);
}

Vector<double>* ProblemFactory::make_right_hand_side() const {
   return new SimpleVector<double>(nx + my + mz);
}

