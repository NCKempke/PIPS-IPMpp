/*
 * QpPostsolver.h
 *
 *  Created on: 03.05.2019
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_QPPOSTSOLVER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_QPPOSTSOLVER_H_

#include "Postsolver.h"

class Problem;

/**
 * Abstract base class for QP Postsolvers.
 */

class QpPostsolver : public Postsolver {
protected:
   const Problem& original_problem;

public:
   QpPostsolver(const Problem& prob);
   ~QpPostsolver() override = default;
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_QPPOSTSOLVER_H_ */
