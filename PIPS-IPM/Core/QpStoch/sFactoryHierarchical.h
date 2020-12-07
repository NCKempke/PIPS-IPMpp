/*
 * sFactoryHierarchical.h
 *
 *  Created on: 10.09.2020
 *      Author: bzfkempk
 */


#ifndef PIPS_IPM_CORE_QPSTOCH_SFACTORYHIERARCHICAL_H_
#define PIPS_IPM_CORE_QPSTOCH_SFACTORYHIERARCHICAL_H_

#include "sFactoryAug.h"

class sFactoryHierarchical : public sFactoryAug {
   public:

      sFactoryHierarchical( StochInputTree* inputTree, MPI_Comm comm=MPI_COMM_WORLD )
      : sFactoryAug(inputTree, comm) {};

//      sFactoryHierarchical( stochasticInput& in, MPI_Comm comm=MPI_COMM_WORLD )
//      : sFactoryAug(in, comm) {};

      virtual ~sFactoryHierarchical() {};

      sLinsysRoot* newLinsysRootHierarchical() override;

      Data* switchToHierarchicalData(Data* prob_in) override;
};

#endif /* PIPS_IPM_CORE_QPSTOCH_SFACTORYHIERARCHICAL_H_ */
