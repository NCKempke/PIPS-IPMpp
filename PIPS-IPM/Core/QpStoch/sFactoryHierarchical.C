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

sLinsysRoot* sFactoryHierarchical::newLinsysRootHierarchical()
{
   return new sLinsysRootBordered(this, data);
}

Data* sFactoryHierarchical::switchToHierarchicalData( Data* prob_in )
{
   data = dynamic_cast<sData*>(prob_in);

   assert( data->exploitingLinkStructure() );
   const int nx_to_shave = data->getNGlobalVars();
   const int myl_to_shave = data->getNGlobalEQConss();
   const int mzl_to_shave = data->getNGlobalINEQConss();

   // adjust tree
   const std::vector<int>& twoLinksStartBlockA = data->getTwoLinksStartBlockA();
   const std::vector<int>& twoLinksStartBlockC = data->getTwoLinksStartBlockC();

   tree = tree->switchToHierarchicalTree( nx_to_shave, myl_to_shave, mzl_to_shave, twoLinksStartBlockA, twoLinksStartBlockC );

   assert( tree->children.size() == 1 );
   assert( tree->isHierarchicalRoot() );

   // adjust data
   data = data->switchToHierarchicalData( tree );

   assert( data->isHierarchieRoot() );

   return data;
}

void sFactoryHierarchical::collapseHierarchicalTree()
{
   sTree* new_top = tree->collapseHierarchicalTree();
   assert( tree->children.size() == 0 );
   delete( tree );
   tree = new_top;

}
