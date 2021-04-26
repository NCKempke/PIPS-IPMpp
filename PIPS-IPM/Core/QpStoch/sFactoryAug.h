/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHACTORYAUG
#define STOCHACTORYAUG

#include "DistributedFactory.h"

class sFactoryAug : public DistributedFactory {
public:
   sFactoryAug(StochInputTree* tree, MPI_Comm comm = MPI_COMM_WORLD);

   ~sFactoryAug() override = default;

   DoubleLinearSolver* newRootSolver() override;

   sLinsysRoot* newLinsysRoot() override;

   sLinsysRoot* newLinsysRootHierarchical() override;

   sLinsysRoot* newLinsysRoot(DistributedQP* prob, OoqpVector* dd, OoqpVector* dq, OoqpVector* nomegaInv, OoqpVector* rhs) override;

   Problem* switchToHierarchicalData(Problem* prob_in) override;

   void switchToOriginalTree() override;

};

#endif
