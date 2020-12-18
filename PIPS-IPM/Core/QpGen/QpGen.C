/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGen.h"
#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"

QpGen::QpGen( int nx_, int my_, int mz_ ) :
   nx( nx_ ), my( my_ ), mz( mz_ )
{
}
