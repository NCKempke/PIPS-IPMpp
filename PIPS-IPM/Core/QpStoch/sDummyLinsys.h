#ifndef SDUMMYLINSYS
#define SDUMMYLINSYS

#include "sLinsys.h"
#include "pipsport.h"

/** 
 * DUMMY Linear system class
 */
class sDummyLinsys : public sLinsys
{
 public:
  sDummyLinsys(sFactory* factory, sData* prob)
    : sLinsys(factory, prob, nullptr, nullptr, nullptr, nullptr) 
    {
      mpiComm = MPI_COMM_NULL;
    };


  virtual ~sDummyLinsys(){};

  void factor2( sData *prob, Variables *vars) override {};
  void Lsolve ( sData *prob, OoqpVector& x ) override {};
  void Dsolve ( sData *prob, OoqpVector& x ) override {};
  void Ltsolve( sData *prob, OoqpVector& x ) override {};

  void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp) override {};

  void putZDiagonal( OoqpVector& zdiag ) override {};
  void solveCompressed( OoqpVector& rhs ) override {};
  void putXDiagonal( OoqpVector& xdiag_ ) override {};

  void joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
		OoqpVector& rhs2_in, OoqpVector& rhs3_in ) override {};

  void separateVars( OoqpVector& x_in, OoqpVector& y_in,
		     OoqpVector& z_in, OoqpVector& vars_in ) override {};
  

  void addLnizi(sData *prob, OoqpVector& z0, OoqpVector& zi) override {};
  void addLniziLinkCons(sData *prob, OoqpVector& z0, OoqpVector& zi, int parentmy, int parentmz) override {};

  /** y += alpha * Lni^T * x */
  void LniTransMult(sData *prob, 
		    SimpleVector& y, 
		    double alpha, SimpleVector& x) override {};

  void addTermToSchurResidual(sData* prob, 
			      SimpleVector& res, 
			      SimpleVector& x) override {};


  void allocU(DenseGenMatrix ** Ut, int np) override {};
  void allocV (DenseGenMatrix ** V, int np) override {};
  void computeU_V(sData *prob, DenseGenMatrix* U, DenseGenMatrix* V) override {};
  void sync() override {};
  void deleteChildren() override {};

  bool isDummy() const override { return true; };
}; 

#endif
