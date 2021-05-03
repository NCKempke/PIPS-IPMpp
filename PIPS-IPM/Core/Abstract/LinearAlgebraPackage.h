/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef LINEARALGEBRA
#define LINEARALGEBRA
/**
 * @defgroup AbstractLinearAlgebra
 *
 * Abstract base classes for linear algebra object (vectors/matrices/solvers).
 * @{
 */

#include "IotrRefCount.h"
#include "Vector.hpp"

class DoubleLinearSolver;
class GenMatrix;
class SymMatrix;

/**
 * A class whose instances creates matrices and vectors of
 * an appropriate type. */
class LinearAlgebraPackage : public IotrRefCount {
protected:
   LinearAlgebraPackage() = default;
   ~LinearAlgebraPackage() override = default;
public:
   /** Create a new symmetric matrix (of appropriate type). */
   virtual SymMatrix* newSymMatrix(int size, int nnz) const = 0;
   /** Create a new non-symmetric matrix (of appropriate type). */
   virtual GenMatrix* newGenMatrix(int m, int n, int nnz) const = 0;
   /** Create a new vector (of appropriate type.) */
   virtual Vector<double>* newVector(int n) const = 0;

   /** Get a string indicating the type of this object
    *  (for debugging purposes.)
    */
   virtual void whatami(char type[32]) const = 0;
};

/**
 * @}
 */
#endif
