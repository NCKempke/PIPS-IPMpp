/*
 * sFactoryHierarchical.h
 *
 *  Created on: 10.09.2020
 *      Author: bzfkempk
 */


#ifndef PIPS_IPM_CORE_QPSTOCH_SFACTORYHIERARCHICAL_H_
#define PIPS_IPM_CORE_QPSTOCH_SFACTORYHIERARCHICAL_H_

#include "sFactoryAug.h"

class sFactoryHierarchical : public sFactory {
 public:

      sFactoryHierarchical( StochInputTree* inputTree, MPI_Comm comm=MPI_COMM_WORLD )
      : sFactory(inputTree, comm) {};

      sFactoryHierarchical( stochasticInput& in, MPI_Comm comm=MPI_COMM_WORLD )
      : sFactory(in,comm) {};

   virtual ~sFactoryHierarchical() {};

   sLinsysRoot* newLinsysRoot() override;
   sLinsysRoot* newLinsysRoot(sData* prob, OoqpVector* dd,OoqpVector* dq,
         OoqpVector* nomegaInv, OoqpVector* rhs) override;

   Data* switchToHierarchicalData(Data* prob_in) override;
};



#endif /* PIPS_IPM_CORE_QPSTOCH_SFACTORYHIERARCHICAL_H_ */
