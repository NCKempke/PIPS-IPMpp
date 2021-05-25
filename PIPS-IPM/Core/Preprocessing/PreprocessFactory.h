/*
 * PreprocessFactory.h
 *
 *  Created on: Dec 30, 2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_
#define PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_

#include "EquiStochScaler.h"
#include "GeoStochScaler.h"
#include "StochPresolver.h"
#include "IotrRefCount.h"
#include "PreprocessType.h"
#include "StochPostsolver.h"
#include "DistributedQP.hpp"

class PreprocessFactory : public IotrRefCount {
public:

   static Scaler* makeScaler(Problem* data, ScalerType type) {
      switch (type) {
         case ScalerType::SCALER_EQUI_STOCH:
            return new EquiStochScaler(data, false);
         case ScalerType::SCALER_GEO_STOCH:
            return new GeoStochScaler(data, false, false);
         case ScalerType::SCALER_GEO_EQUI_STOCH:
            return new GeoStochScaler(data, true, false);
         default:
            return nullptr;
      }
   };

   static Presolver* makePresolver(DistributedTree* tree, const Problem* data, PresolverType type, Postsolver* postsolver = nullptr) {
      assert(data);
      switch (type) {
         case PresolverType::PRESOLVER_STOCH:
            return new StochPresolver(tree, *data, postsolver);
         default:
            return nullptr;
      }
   };

   // todo : not sure about the factory design here
   static Postsolver* makePostsolver(const Problem* original_problem) {
      return new StochPostsolver(dynamic_cast<const DistributedQP&>(*original_problem));
   }
};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_PREPROCESSFACTORY_H_ */
