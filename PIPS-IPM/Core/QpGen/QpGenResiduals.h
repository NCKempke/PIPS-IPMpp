/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENRESIDUALS
#define QPGENRESIDUALS

#include <iostream>
#include <fstream>

#include "Residuals.h"
#include "OoqpVectorHandle.h"

class QpGen;
class QpGenData;
class Variables;
class LinearAlgebraPackage;

/** 
 * Residuals for the general QP formulation 
 *
 * @ingroup QpGen
 */

class QpGenResiduals : public Residuals {
protected:
  long long nx{0};
  long long my{0};
  long long mz{0};

  long long nxupp{0};
  OoqpVectorHandle ixupp;

  long long nxlow{0};
  OoqpVectorHandle ixlow;

  long long mcupp{0};
  OoqpVectorHandle icupp;

  long long mclow{0};
  OoqpVectorHandle iclow;

  QpGenResiduals() = default;

public:
  OoqpVectorHandle rQ;
  OoqpVectorHandle rA;
  OoqpVectorHandle rC;
  OoqpVectorHandle rz;
  OoqpVectorHandle rv;
  OoqpVectorHandle rw;
  OoqpVectorHandle rt;
  OoqpVectorHandle ru;
  OoqpVectorHandle rgamma;
  OoqpVectorHandle rphi;
  OoqpVectorHandle rlambda;
  OoqpVectorHandle rpi;

  QpGenResiduals( LinearAlgebraPackage * la,
		  long long nx, long long my, long long mz,
		  OoqpVector * ixlow, OoqpVector * ixupp,
		  OoqpVector * iclow, OoqpVector * icupp );

  QpGenResiduals( const QpGenResiduals& res);

  const long long& getNxupp() { return nxupp; };
  const long long& getNxlow() { return nxlow; };
  const long long& getMcupp() { return mcupp; };
  const long long& getMclow() { return mclow; };

  virtual ~QpGenResiduals() = default;
  
  void calcresids(Problem* problem, Variables *vars, bool print_resids = false) override;

  void add_r3_xz_alpha(const Variables *vars, double alpha) override;

  double recomputeResidualNorm() override;

  void set_r3_xz_alpha(const Variables *vars, double alpha) override;
  
  void clear_r3() override;
  
  void clear_r1r2() override;

  void project_r3(double rmin, double rmax) override;

  virtual int  validNonZeroPattern();
  
  virtual void writeToStream(std::ostream& out);

  void copyFrom( const Residuals& other ) override;
};

#endif





