/*
 * sFactoryHierarchical.C
 *
 *  Created on: 10.09.2020
 *      Author: bzfkempk
 */

#include "sFactoryHierarchical.h"
#include "sData.h"
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

   // TODO: get sizes n0Links, global F and G constraints -> then shave Matrices and Vectors and put everything in higher level of hierarchy

   // TODO: adjust tree so that Residuals and Variables (so Stochvectors) can be mades with it (one more layer in tree with sizes accordingly

   // TODO: adjust Data -> border level will have the border matrices and on bottom right matrix - need Amul Cmul and other functions here
   data->switchToHierarchicalData();

   // TODO: second implementation step for hierarchical approach : add recursive linear systems inbetween border layer and the bottom layers..
   assert(0 && "todo : implement ;(");
   return data;
}

