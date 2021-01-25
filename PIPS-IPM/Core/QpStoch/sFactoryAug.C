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

   // TODO : DELETEME
//   OoqpVector* x_bef = tree->newPrimalVector();
//   OoqpVector* y_bef = tree->newDualYVector();
//   OoqpVector* z_bef = tree-   >newDualZVector();
//   x_bef->setToConstant(2.0);
//   y_bef->setToConstant(2.0);
//   z_bef->setToConstant(2.0);
//
//   data->A->mult(2.0, *y_bef, 3.0, *x_bef);
//   data->C->mult(2.0, *z_bef, 3.0, *x_bef);
//
//   const double A2norm_bef = y_bef->twonorm();
//   const double A1norm_bef = y_bef->onenorm();
//
//   const double C2norm_bef = z_bef->twonorm();
//   const double C1norm_bef = z_bef->onenorm();

   assert( data->exploitingLinkStructure() );
   const int nx_to_shave = data->getNGlobalVars();
   const int myl_to_shave = data->getNGlobalEQConss();
   const int mzl_to_shave = data->getNGlobalINEQConss();

   // adjust tree
   const std::vector<int>& twoLinksStartBlockA = data->getTwoLinksStartBlockA();
   const std::vector<int>& twoLinksStartBlockC = data->getTwoLinksStartBlockC();

   tree = tree->switchToHierarchicalTree( nx_to_shave, myl_to_shave, mzl_to_shave, twoLinksStartBlockA, twoLinksStartBlockC );

   assert( tree->getChildren().size() == 1 );
   assert( tree->isHierarchicalRoot() );

   // adjust data
   data = data->switchToHierarchicalData( tree );

   assert( data->isHierarchyRoot() );

// TODO : DELETEME
//   OoqpVector* x_after = tree->newPrimalVector();
//   OoqpVector* y_after = tree->newDualYVector();
//   OoqpVector* z_after = tree->newDualZVector();
//   x_after->setToConstant(2.0);
//   y_after->setToConstant(2.0);
//   z_after->setToConstant(2.0);
//   data->A->mult(2.0, *y_after, 3.0, *x_after);
//   data->C->mult(2.0, *z_after, 3.0, *x_after);
//   const double A2norm_after = y_after->twonorm();
//   const double A1norm_after = y_after->onenorm();
//
//   const double C2norm_after = z_after->twonorm();
//   const double C1norm_after = z_after->onenorm();
//
//   std::cout << "A2norm before : " << A2norm_bef << " vs A2norm after : " << A2norm_after << " difference " << A2norm_bef - A2norm_after << "\n";
//   std::cout << "C2norm before : " << C2norm_bef << " vs C2norm after : " << C2norm_after << " difference " << C2norm_bef - C2norm_after << "\n";
//   std::cout << "\n";
//   std::cout << "A1norm before : " << A1norm_bef << " vs A1norm after : " << A1norm_after << " difference " << A1norm_bef - A1norm_after << "\n";
//   std::cout << "C1norm before : " << C1norm_bef << " vs C1norm after : " << C1norm_after << " difference " << C1norm_bef - C1norm_after << "\n";

   return data;
}

void sFactoryAug::collapseHierarchicalTree()
{
   sTree* new_top = tree->collapseHierarchicalTree();
   assert( tree->getChildren().size() == 0 );
   delete( tree );
   tree = new_top;
}
