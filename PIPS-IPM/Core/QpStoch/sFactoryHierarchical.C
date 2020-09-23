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
#include "sLinsysRootBordered.h"

sLinsysRoot* sFactoryHierarchical::newLinsysRoot()
{
   return new sLinsysRootBordered(this, data);
}

sLinsysRoot* sFactoryHierarchical::newLinsysRoot(sData *prob, OoqpVector *dd, OoqpVector *dq,
      OoqpVector *nomegaInv, OoqpVector *rhs)
{
   return new sLinsysRootBordered(this, prob, dd, dq, nomegaInv, rhs);
}

Data* sFactoryHierarchical::switchToHierarchicalData( Data* prob_in )
{
   data = dynamic_cast<sData*>(prob_in);

   assert( data->exploitingLinkStructure() );
   const int nx_to_shave = data->getNGlobalVars();
   const int myl_to_shave = data->getNGlobalEQConss();
   const int mzl_to_shave = data->getNGlobalINEQConss();

   // TODO : lower hierarchy missing
   // adjust tree
   tree = tree->switchToHierarchicalTree(nx_to_shave, myl_to_shave, mzl_to_shave);

   assert( tree->children.size() == 1 );
   assert( tree->isHierarchicalRoot() );

   // adjust data
   data = data->switchToHierarchicalData( tree );

   return data;
}
