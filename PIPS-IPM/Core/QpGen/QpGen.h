/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENFACTORY
#define QPGENFACTORY

#include "ProblemFormulation.h"
#include "OoqpVector_fwd.h"

class Data;
class Residuals;
class LinearSystem;
class Variables;
class LinearAlgebraPackage;

/**
 * @defgroup QpGen
 * 
 * OOQP's default general problem formulation:
 *
 *  <pre>
 *  minimize    c' x + ( 1/2 ) x' * Q x         ; 
 *  subject to                      A x  = b    ;
 *                         clow <=  C x <= cupp ;
 *                         xlow <=    x <= xupp ;
 *  </pre> 
 *
 *  The general linear equality constraints must have either an upper
 *  or lower bound, but need not have both bounds. The variables may have 
 *  no bounds; an upper bound; a lower bound or both an upper and lower
 *  bound.
*/
class QpGen : public ProblemFormulation {
protected:
  LinearAlgebraPackage * la;
  /** number of elements in x */
  long long nx;

  /** number of rows in A and b */
  long long my;

  /** number of rows in C */
  long long mz;

  QpGen( int nx_, int my_, int mz_ );
public:
   Residuals     * makeResiduals( Data * prob_in ) override;
   Variables     * makeVariables( Data * prob_in ) override;

   virtual void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			OoqpVector& rhs2_in, OoqpVector& rhs3_in ) = 0;

   virtual void separateVars( OoqpVector& x_in, OoqpVector& y_in,
			     OoqpVector& z_in, OoqpVector& vars_in ) = 0;

   virtual void writeProblemToStream(std::ostream& out) const;

   virtual ~QpGen() {};
};

#endif



