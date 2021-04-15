/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGen.h"
#include "QuadraticProblem.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "LinearAlgebraPackage.h"


QpGen::QpGen( int nx_, int my_, int mz_ ) :
   nx( nx_ ), my( my_ ), mz( mz_ )
{
}

OoqpVector* QpGen::makePrimalVector() const
{
   assert(la);
   return la->newVector( nx );
}

OoqpVector* QpGen::makeDualYVector() const
{
   assert(la);
   return la->newVector( my );

}

OoqpVector* QpGen::makeDualZVector() const
{
   assert(la);
   return la->newVector( mz );
}

OoqpVector* QpGen::makeRhs() const
{
   assert(la);
   return la->newVector( nx + my + mz );
}

