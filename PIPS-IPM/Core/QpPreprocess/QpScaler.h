/*
 * QpScaler.h
 *
 *  Created on: 19.12.2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_QPSCALER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_QPSCALER_H_


#include "Scaler.h"
#include "Vector.hpp"
#include "SmartPointer.h"
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
class QpScaler : public Scaler {
protected:
   static void invertAndRound(bool round, Vector<double>& vector) {
      vector.invertSave(1.0);
      if (round)
         vector.roundToPow2();
   }

   const bool scaling_output{false};

   // has scaling been applied
   bool scaling_applied{false};

   // scaling vector
   std::unique_ptr<Vector<double>> vec_rowscaleQ{};
   std::unique_ptr<Vector<double>> vec_rowscaleA{};
   std::unique_ptr<Vector<double>> vec_rowscaleC{};
   std::unique_ptr<Vector<double>> vec_colscale{};

   // problem data
   SymMatrixHandle Q;
   GenMatrixHandle A;
   GenMatrixHandle C;
   SmartPointer<Vector<double> > obj;
   SmartPointer<Vector<double> > bA;
   SmartPointer<Vector<double> > bux;
   SmartPointer<Vector<double> > blx;
   SmartPointer<Vector<double> > rhsC;
   SmartPointer<Vector<double> > lhsC;

   // scaling factor for objective
   double factor_objscale;

   void applyScaling();

   virtual void doObjScaling() = 0;

   /** get maximum absolute row ratio and write maximum row entries into vectors */
   double
   maxRowRatio(Vector<double>& maxvecA, Vector<double>& maxvecC, Vector<double>& minvecA, Vector<double>& minvecC, const Vector<double>* colScalevec);

   /** get maximum absolute column ratio and write maximum column entries into vectors */
   double maxColRatio(Vector<double>& maxvec, Vector<double>& minvec, const Vector<double>* rowScaleVecA, const Vector<double>* rowScaleVecC);

   void scaleObjVector(double scaling_factor);

   /** print row col ratios */
   void printRowColRatio();

   void setScalingVecsToOne();
public:

   QpScaler(Problem* problem, bool bitshifting = false);
   ~QpScaler() override = default;

   /** scale */
   void scale() override = 0;

   double getObjUnscaled(double objval) const override;

   Variables* getVariablesUnscaled(const Variables& variables) const override;
   Residuals* getResidualsUnscaled(const Residuals& residuals) const override;

   void unscaleVariables(Variables& vars) const override;
   void unscaleResiduals(Residuals& resids) const override;

   Vector<double>* getPrimalUnscaled(const Vector<double>& solprimal) const override;
   Vector<double>* getDualEqUnscaled(const Vector<double>& soldual) const override;
   Vector<double>* getDualIneqUnscaled(const Vector<double>& soldual) const override;
   Vector<double>* getDualVarBoundsUppUnscaled(const Vector<double>& soldual) const override;
   Vector<double>* getDualVarBoundsLowUnscaled(const Vector<double>& soldual) const override;
};

//@}


#endif /* PIPS_IPM_CORE_QPPREPROCESS_QPSCALER_H_ */
