/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENDATA
#define QPGENDATA

#include "Problem.h"
#include "OoqpVectorHandle.h"
#include "OoqpVector.h"
#include "DoubleMatrixHandle.h"

class MpsReader;

class LinearAlgebraPackage;

class Variables;

#ifdef TESTING
class QpGenDataTester;
#endif

/**
 * Data for the general QP formulation.
 *
 * @ingroup QpGen
 */

class QP : public Problem {
#ifdef TESTING
   friend QpGenDataTester;
#endif

protected:

   QP() = default;

public:
   SymMatrixHandle Q;

   /** constructor that sets up pointers to the data objects that are
       passed as arguments */
   QP(LinearAlgebraPackage* la, OoqpVector* c, SymMatrix* Q, OoqpVector* xlow, OoqpVector* ixlow, OoqpVector* xupp, OoqpVector* ixupp,
         GenMatrix* A, OoqpVector* bA, GenMatrix* C, OoqpVector* clow, OoqpVector* iclow, OoqpVector* cupp, OoqpVector* ciupp);

   /** insert the Hessian Q into the matrix M for the fundamental linear system, where M is stored as a SymMatrix */
   virtual void putQIntoAt(SymMatrix& M, int row, int col);

   /** insert the Hessian Q into the matrix M for the fundamental linear system, where M is stored as a GenMatrix */
   virtual void putQIntoAt(GenMatrix& M, int row, int col);

   /** y = beta * y + alpha * Q * x */
   virtual void hessian_multiplication(double beta, OoqpVector& y, double alpha, const OoqpVector& x) const override;

   /** extract the diagonal of the Hessian and put it in the OoqpVector hessian_diagonal */
   void hessian_diagonal(OoqpVector& hessian_diagonal) override;

   void createScaleFromQ();

   void scaleQ();

   void flipQ();

   double datanorm() const override;

   virtual void datainput() {};

   virtual void datainput(MpsReader* reader, int& iErr);

   void print() override;

   virtual void objective_gradient(const Variables* vars, OoqpVector& gradient) const override;

   virtual double objective_value(const Variables* vars) const override;

   ~QP() override = default;
};

#endif
