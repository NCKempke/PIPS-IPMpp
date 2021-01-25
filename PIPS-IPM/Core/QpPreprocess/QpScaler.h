/*
 * QpScaler.h
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_QPSCALER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_QPSCALER_H_


#include "Scaler.h"
#include "OoqpVector.h"
#include "DoubleMatrix.h"

#include <memory>

/**  * @defgroup QpPreprocess
 *
 * QP scaler
 * @{
 */

/**
 * Abstract base class for QP scalers.
 */
class QpScaler : public Scaler
{
protected:
  static void invertAndRound(bool round, OoqpVector& vector)
  {
     vector.invertSave(1.0);
     if( round )
        vector.roundToPow2();
  }

  const bool scaling_output{false};

  // has scaling been applied
  bool scaling_applied{false};

  // scaling vector
  std::unique_ptr<OoqpVector> vec_rowscaleQ{};
  std::unique_ptr<OoqpVector> vec_rowscaleA{};
  std::unique_ptr<OoqpVector> vec_rowscaleC{};
  std::unique_ptr<OoqpVector> vec_colscale{};

  // problem data
  SymMatrixHandle Q;
  GenMatrixHandle A;
  GenMatrixHandle C;
  OoqpVectorHandle obj;
  OoqpVectorHandle bA;
  OoqpVectorHandle bux;
  OoqpVectorHandle blx;
  OoqpVectorHandle rhsC;
  OoqpVectorHandle lhsC;

  // scaling factor for objective
  double factor_objscale;

  void applyScaling();

  virtual void doObjScaling() = 0;

  /** get maximum absolute row ratio and write maximum row entries into vectors */
  double maxRowRatio(OoqpVector& maxvecA, OoqpVector& maxvecC, OoqpVector& minvecA, OoqpVector& minvecC, const OoqpVector* colScalevec);

  /** get maximum absolute column ratio and write maximum column entries into vectors */
  double maxColRatio(OoqpVector& maxvec, OoqpVector& minvec, const OoqpVector* rowScaleVecA,  const OoqpVector* rowScaleVecC);

  void scaleObjVector(double scaling_factor);

  /** print row col ratios */
  void printRowColRatio();

  void setScalingVecsToOne();
public:

  QpScaler(Data* prob, bool bitshifting = false);
  ~QpScaler() override = default;

  /** scale */
  void scale() override = 0;

  double getObjUnscaled(double objval) const override;

  Variables* getVariablesUnscaled(const Variables& vars) const override;
  Residuals* getResidualsUnscaled(const Residuals& resids) const override;

  void unscaleVariables( Variables& vars ) const override;
  void unscaleResiduals( Residuals& resids ) const override;

  OoqpVector* getPrimalUnscaled(const OoqpVector& solprimal) const override;
  OoqpVector* getDualEqUnscaled(const OoqpVector& soldual) const override;
  OoqpVector* getDualIneqUnscaled(const OoqpVector& soldual) const override;
  OoqpVector* getDualVarBoundsUppUnscaled(const OoqpVector& soldual) const override;
  OoqpVector* getDualVarBoundsLowUnscaled(const OoqpVector& soldual) const override;
};

//@}


#endif /* PIPS_IPM_CORE_QPPREPROCESS_QPSCALER_H_ */
