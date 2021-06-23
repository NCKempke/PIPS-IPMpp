/*
 * PreprocessFactory.h
 *
 *  Created on: Dec 30, 2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_
#define PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_

#include "EquilibriumScaler.h"
#include "GeometricMeanScaler.h"
#include "CurtisReidScaler.h"
#include "StochPresolver.h"
#include "PreprocessType.h"
#include "StochPostsolver.h"
#include "DistributedProblem.hpp"

class ProblemFactory;

class PreprocessFactory {
public:

   static Scaler* make_scaler(const ProblemFactory& problem_factory, const Problem& problem, ScalerType type) {
      switch (type) {
         case ScalerType::EQUILIBRIUM:
            return new EquilibriumScaler(problem_factory, problem, false);
         case ScalerType::GEOMETRIC_MEAN:
            return new GeometricMeanScaler(problem_factory, problem, false, false);
         case ScalerType::GEOMETRIC_MEAN_EQUILIBRIUM:
            return new GeometricMeanScaler(problem_factory, problem, true, false);
         case ScalerType::CURTIS_REID:
            return new CurtisReidScaler(problem_factory, problem, false);
         default:
            return nullptr;
      }
   };

   static Presolver* make_presolver(DistributedTree& tree, const Problem* data, PresolverType type, Postsolver* postsolver = nullptr) {
      assert(data);
      switch (type) {
         case PresolverType::PRESOLVE:
            return new StochPresolver(tree, *data, postsolver);
         default:
            return nullptr;
      }
   };

   // todo : not sure about the factory design here
   static Postsolver* make_postsolver(const Problem* original_problem) {
      return new StochPostsolver(dynamic_cast<const DistributedProblem&>(*original_problem));
   }
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_ */
