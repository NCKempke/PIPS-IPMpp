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
  sDummyLinsys(sFactory* factory, DistributedQP* prob)
    : sLinsys(factory, prob, nullptr, nullptr, nullptr, nullptr, false)
 {
     mpiComm = MPI_COMM_NULL;
 };

  ~sDummyLinsys() override = default;

  void factor2( DistributedQP*, Variables*) override {};
  void allreduceAndFactorKKT(DistributedQP*, Variables*) override {};
  void assembleKKT(DistributedQP*, Variables*) override {}

  void Lsolve( DistributedQP*, OoqpVector& ) override {};
  void Dsolve( DistributedQP*, OoqpVector& ) override {};
  void Ltsolve( DistributedQP*, OoqpVector& ) override {};
  void Ltsolve2( DistributedQP*, StochVector&, SimpleVector&, bool ) override {};

  void putZDiagonal( OoqpVector& ) override {};
  void solveCompressed( OoqpVector& ) override {};
  void putXDiagonal( OoqpVector& ) override {};

  void joinRHS( OoqpVector&, const OoqpVector&, const OoqpVector&, const OoqpVector& ) const override {};

  void separateVars( OoqpVector&, OoqpVector&, OoqpVector&, const OoqpVector& ) const override {};

  void addLnizi(DistributedQP*, OoqpVector&, OoqpVector& ) override {};
  void addLniziLinkCons(DistributedQP*, OoqpVector&, OoqpVector&, bool ) override {};

  /** y += alpha * Lni^T * x */
  //  void LniTransMult(DistributedQP *prob, SimpleVector& y, double alpha, SimpleVector& x) override {};

  void addTermToSchurResidual( DistributedQP*, SimpleVector&, SimpleVector& ) override {};

  void LsolveHierarchyBorder( DenseGenMatrix&, BorderLinsys&, std::vector<BorderMod>&, bool, int, int ) override {};
  void LsolveHierarchyBorder( DenseGenMatrix&, BorderLinsys&, std::vector<BorderMod>&, bool, bool, int, int ) override {};
  void addInnerBorderKiInvBrToRes( DoubleMatrix&, BorderLinsys&, std::vector<BorderMod>&, bool, bool, bool, int, int, int ) override {};
  void LniTransMultHierarchyBorder( DoubleMatrix&, const DenseGenMatrix&, BorderLinsys&, BorderLinsys&, std::vector<BorderMod>&, bool, bool, bool, int, int, int ) override {};

  void allocU( DenseGenMatrix**, int ) override {};
  void allocV( DenseGenMatrix**, int ) override {};
  void computeU_V( DistributedQP*, DenseGenMatrix*, DenseGenMatrix* ) override {};
  void deleteChildren() override {};

  bool isDummy() const override { return true; };

  void addBorderTimesRhsToB0( StochVector&, SimpleVector&, BorderLinsys& ) override {};
  void addBorderX0ToRhs( StochVector&, const SimpleVector&, BorderLinsys& ) override {};
  void computeInnerSystemRightHandSide( StochVector&, const SimpleVector&, bool ) override {};
};

#endif
