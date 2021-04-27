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

   DoubleLinearSolver* make_root_solver() override;

   sLinsysRoot* make_linear_system_root() override;

   sLinsysRoot* newLinsysRootHierarchical() override;

   sLinsysRoot* make_linear_system_root(DistributedQP* problem, OoqpVector* dd, OoqpVector* dq, OoqpVector* nomegaInv, OoqpVector* rhs) override;

   Problem* switchToHierarchicalData(Problem* problem);

   void switchToOriginalTree();

};

#endif
