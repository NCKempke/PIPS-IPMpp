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
  LinearAlgebraPackage * la{};
  /** number of elements in x */
  long long nx{0};

  /** number of rows in A and b including linking rows (sFactory..) */
  long long my{0};

  /** number of rows in C including linking rows */
  long long mz{0};

  QpGen() = default;
  QpGen( int nx_, int my_, int mz_ );
public:

  /** create x shaped vector using LinearAlgebraPackage */
  virtual OoqpVector* makePrimalVector() const;
  /** create dual A shaped vector using LinearAlgebraPackage */
  virtual OoqpVector* makeDualYVector() const;
  /** create dual C shaped vector using LinearAlgebraPackage */
  virtual OoqpVector* makeDualZVector() const;
  /** create a rhs vector for the augmented system */
  virtual OoqpVector* makeRhs() const;

  virtual void joinRHS( OoqpVector& rhs_in, const OoqpVector& rhs1_in,
        const OoqpVector& rhs2_in, const OoqpVector& rhs3_in ) const = 0;

  virtual void separateVars( OoqpVector& x_in, OoqpVector& y_in,
        OoqpVector& z_in, const OoqpVector& vars_in ) const = 0;

  ~QpGen() override = default;
};

#endif



