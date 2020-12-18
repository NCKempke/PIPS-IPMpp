/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sFactoryAug.h"

#include "sData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "sLinsysRootAug.h"
#include "sLinsysRootBordered.h"

sFactoryAug::sFactoryAug( StochInputTree* inputTree, MPI_Comm comm)
  : sFactory(inputTree, comm)
{ };

sLinsysRoot* sFactoryAug::newLinsysRoot()
{
   assert( data );
   return new sLinsysRootAug(this, data);
}

sLinsysRoot* sFactoryAug::newLinsysRoot(sData* prob, sTree* tree_,
			   OoqpVector* dd,OoqpVector* dq,
			   OoqpVector* nomegaInv, OoqpVector* rhs)
{

   return new sLinsysRootAug(this, tree_, prob, dd, dq, nomegaInv, rhs);
}

sLinsysRoot* sFactoryAug::newLinsysRootHierarchical()
{
   return new sLinsysRootBordered(this, data);
}

Data* sFactoryAug::switchToHierarchicalData( Data* prob_in )
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

   assert( tree->nChildren() == 1 );
   assert( tree->isHierarchicalRoot() );

   // adjust data
   data = data->switchToHierarchicalData( tree );

   assert( data->isHierarchieRoot() );

   return data;
}

void sFactoryAug::collapseHierarchicalTree()
{
   sTree* new_top = tree->collapseHierarchicalTree();
   assert( tree->nChildren() == 0 );
   delete( tree );
   tree = new_top;

}
