/*
 * StochPresolver.h
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_

#include "QpPresolver.h"
#include "StochVector.h"
#include "StochGenMatrix.h"
#include "SmartPointer.h"
#include "sData.h"

class Data;


/**  * @defgroup QpPreprocess
 *
 * QP presolver
 * @{
 */

/**
 * Derived class for QP presolvers.
 */
class StochPresolver : public QpPresolver
{
private:
  static const double feastol = 1.0e-6;
  static const double infinity = 10e30;

  // number of non-zero elements of each row
  StochVectorHandle nRowElemsA;
  StochVectorHandle nRowElemsC;

  // number of non-zero elements of each column
  StochVectorHandle nColElems;

  // initialize row and column nnz counter
  void initNnzCounter();

  // remove small matrix entries and return number of eliminations
  int cleanUp();

  int cleanUpRowsC(StochGenMatrix* matrixC, StochVector* clow, StochVector* cupp,
        StochVector* xlow, StochVector* xupp, StochVector* nnzRowC);

  sData* presProb;

public:

  StochPresolver(const Data* prob);
  virtual ~StochPresolver();

  virtual Data* presolve();
};

//@}




#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_ */
