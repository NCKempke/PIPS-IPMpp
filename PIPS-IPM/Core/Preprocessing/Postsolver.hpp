/*
 * QpPostsolver.h
 *
 *  Created on: 03.05.2019
 *      Author: bzfkempk
 */
#ifndef QPPOSTSOLVER_H
#define QPPOSTSOLVER_H

#include "TerminationStatus.hpp"

class Problem;

class Variables;

enum PostsolveStatus {
   PRESOLVE_OK, PRESOLVE_FAIL // TODO
};

/**
 * Abstract base class for QP Postsolvers.
 */

class Postsolver {
public:
   explicit Postsolver(const Problem& problem);
   virtual ~Postsolver() = default;

   /** postsolve reduced solution and set original solution accordingly */
   virtual PostsolveStatus postsolve(const Variables& reduced_solution, Variables& original_solution, TerminationStatus result_code) = 0;

protected:
   const Problem& original_problem;
};

#endif /* QPPOSTSOLVER_H */
