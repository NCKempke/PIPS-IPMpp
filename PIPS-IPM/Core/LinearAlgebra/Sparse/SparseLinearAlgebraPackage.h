/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SPARSELINEARALGEBRA
#define SPARSELINEARALGEBRA

#include "LinearAlgebraPackage.h"
#include "Vector.hpp"
#include "SmartPointer.h"
/**
 * @defgroup SparseLinearAlgebra
 *
 * A module for sparse linear operators and equations
 * @{
 */

/**
 * Creates sparse matrices and SimpleVectors.
 *
 * Singleton class. Only accessible through 
 * SparseLinearAlgebraPackage::soleInstance
 */
class SparseLinearAlgebraPackage : public LinearAlgebraPackage {
protected:
   SparseLinearAlgebraPackage() = default;
   ~SparseLinearAlgebraPackage() override = default;
public:
   /** Return the sole instance of this class. This instance must not
    * ever be deleted. */
   static SparseLinearAlgebraPackage* soleInstance();
   SymMatrix* newSymMatrix(int size, int nnz) const override;
   GenMatrix* newGenMatrix(int m, int n, int nnz) const override;
   Vector<double>* newVector(int n) const override;
   void whatami(char type[32]) const override;

};

/**
 * @}
 */

#endif
