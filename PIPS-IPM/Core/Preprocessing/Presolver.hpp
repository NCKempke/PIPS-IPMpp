/*
 * QpPresolver.h
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#ifndef PRESOLVER_H
#define PRESOLVER_H

#include <cstddef>

class Problem;

class Postsolver;

/**
 * Abstract base class for Presolvers.
 */

class Presolver {
public:
   explicit Presolver(const Problem& problem, Postsolver* postsolver = nullptr);
   virtual ~Presolver() = default;

   /** presolve and return pointer to presolved data */
   virtual Problem* presolve() = 0;

protected:
   const Problem& original_problem;
   Postsolver* const postsolver{};
};

#endif /* PRESOLVER_H */
