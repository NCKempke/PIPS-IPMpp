/*
 * QpPresolver.h
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_QPPRESOLVER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_QPPRESOLVER_H_

#include "Presolver.h"
#include "pipsport.h"

/**  * @defgroup QpPreprocess
 *
 * QP scaler
 * @{
 */

/**
 * Abstract base class for QP Presolvers.
 */

class QpPresolver : public Presolver {
public:
   QpPresolver(const Problem& prob, Postsolver* postsolver = nullptr);
   ~QpPresolver() override = default;
};


//@}


#endif /* PIPS_IPM_CORE_QPPREPROCESS_QPPRESOLVER_H_ */
