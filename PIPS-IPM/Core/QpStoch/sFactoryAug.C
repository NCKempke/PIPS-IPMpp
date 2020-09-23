/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sFactoryAug.h"

#include "sData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "sLinsysRootAug.h"

sFactoryAug::sFactoryAug( StochInputTree* inputTree, MPI_Comm comm)
  : sFactory(inputTree, comm)
{ };

sFactoryAug::sFactoryAug( stochasticInput& in, MPI_Comm comm)
  : sFactory(in,comm)
{ };

sFactoryAug::sFactoryAug()
{ };

sFactoryAug::~sFactoryAug()
{ };


sLinsysRoot* sFactoryAug::newLinsysRoot()
{
  return new sLinsysRootAug(this, data);
}

sLinsysRoot* 
sFactoryAug::newLinsysRoot(sData* prob,
			   OoqpVector* dd,OoqpVector* dq,
			   OoqpVector* nomegaInv, OoqpVector* rhs)
{
  return new sLinsysRootAug(this, prob, dd, dq, nomegaInv, rhs);
}
