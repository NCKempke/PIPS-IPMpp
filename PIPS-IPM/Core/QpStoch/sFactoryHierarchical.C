/*
 * sFactoryHierarchical.C
 *
 *  Created on: 10.09.2020
 *      Author: bzfkempk
 */

#include "sFactoryHierarchical.h"
#include "sData.h"
#include "sTree.h"
#include <cassert>

sLinsysRoot* sFactoryHierarchical::newLinsysRoot()
{
   // TODO return new root linear system (linear system with inner matrix and border only)
   //return new sLinsysRootBordered(this, data);

   assert("not implemented" && 0 );
   return nullptr;
}

sLinsysRoot* sFactoryHierarchical::newLinsysRoot(sData *prob, OoqpVector *dd, OoqpVector *dq,
      OoqpVector *nomegaInv, OoqpVector *rhs)
{
   // TODO return new root linear system (linear system with inner matrix and border only)
   //return new sLinsysRootBordered(this, prob, dd, dq, nomegaInv, rhs);

   assert("not implemented" && 0 );
   return nullptr;
}

Data* sFactoryHierarchical::switchToHierarchicalData( Data* prob_in )
{
   data = dynamic_cast<sData*>(prob_in);

   assert( data->exploitingLinkStructure() );
   const int nx_to_shave = data->getNGlobalVars();
   const int myl_to_shave = data->getNGlobalEQConss();
   const int mzl_to_shave = data->getNGlobalINEQConss();

   // adjust tree
   // TODO : decide on how to split the lower tree levels
   tree = tree->switchToHierarchicalTree(nx_to_shave, myl_to_shave, mzl_to_shave);
   assert( tree->children.size() == 1 );
   assert( tree->isHierarchicalRoot() );

   // adjust data
   data = data->switchToHierarchicalData( tree );

   return data;
}
