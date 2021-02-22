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
    : sLinsys(factory, prob, nullptr, nullptr, nullptr, nullptr, false)
 {
     mpiComm = MPI_COMM_NULL;
 };

  ~sDummyLinsys() override = default;

  void factor2( sData*, Variables*) override {};
  void Lsolve( sData*, OoqpVector& ) override {};
  void Dsolve( sData*, OoqpVector& ) override {};
  void Ltsolve( sData*, OoqpVector& ) override {};

  void Ltsolve2( sData*, StochVector&, SimpleVector& ) override {};

  void putZDiagonal( OoqpVector& ) override {};
  void solveCompressed( OoqpVector& ) override {};
  void putXDiagonal( OoqpVector& ) override {};

  void joinRHS( OoqpVector&, const OoqpVector&, const OoqpVector&, const OoqpVector& ) const override {};

  void separateVars( OoqpVector&, OoqpVector&, OoqpVector&, const OoqpVector& ) const override {};

  void addLnizi(sData*, OoqpVector&, OoqpVector& ) override {};
  void addLniziLinkCons(sData*, OoqpVector&, OoqpVector&, bool ) override {};

  /** y += alpha * Lni^T * x */
  //  void LniTransMult(sData *prob, SimpleVector& y, double alpha, SimpleVector& x) override {};

  void addTermToSchurResidual( sData*, SimpleVector&, SimpleVector& ) override {};

  void LsolveHierarchyBorder( DenseGenMatrix&, BorderLinsys&, std::vector<BorderMod>& ) override {};
  void addInnerBorderKiInvBrToRes( DenseGenMatrix&, BorderLinsys&, std::vector<BorderMod>& ) override {};
  void LniTransMultHierarchyBorder( DoubleMatrix&, const DenseGenMatrix&, BorderLinsys&, BorderLinsys&, std::vector<BorderMod>&, bool, bool ) override {};

  void allocU( DenseGenMatrix**, int ) override {};
  void allocV( DenseGenMatrix**, int ) override {};
  void computeU_V( sData*, DenseGenMatrix*, DenseGenMatrix* ) override {};
  void deleteChildren() override {};

  bool isDummy() const override { return true; };

  void addBorderTimesRhsToB0( StochVector&, SimpleVector&, BorderLinsys& ) override {};
  void addBorderX0ToRhs( StochVector&, const SimpleVector&, BorderLinsys& ) override {};
  void computeInnerSystemRightHandSide( StochVector&, const SimpleVector& ) override {};
};

#endif
