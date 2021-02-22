/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHACTORYAUG
#define STOCHACTORYAUG

#include "sFactory.h"

class sFactoryAug : public sFactory {

public:
  sFactoryAug( StochInputTree*, MPI_Comm comm = MPI_COMM_WORLD );
 ~sFactoryAug() override = default;

  sLinsysRoot* newLinsysRoot() override;
  sLinsysRoot* newLinsysRootHierarchical() override;

  sLinsysRoot* newLinsysRoot(sData* prob, OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs) override;

  Data* switchToHierarchicalData(Data* prob_in) override;
  void switchToOriginalTree() override;

};
#endif
