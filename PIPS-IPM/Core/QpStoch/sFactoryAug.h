/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHACTORYAUG
#define STOCHACTORYAUG

#include "sFactory.h"

class sFactoryAug : public sFactory {
 public:
  sFactoryAug( StochInputTree*, MPI_Comm comm=MPI_COMM_WORLD );
//  sFactoryAug( stochasticInput&, MPI_Comm comm=MPI_COMM_WORLD );
 private:
  sFactoryAug();
 public:
  virtual ~sFactoryAug();

  sLinsysRoot* newLinsysRoot() override;
  sLinsysRoot* newLinsysRoot(sData* prob, sTree* tree_,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs) override;
};
#endif
