/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef PROBLEM_H
#define PROBLEM_H

#include "Variables.h"
#include "OoqpVector.h"
#include "DoubleMatrixHandle.h"
#include "DoubleMatrix.h"

class LinearAlgebraPackage;

class Problem {
protected:
   Problem() = default;

   LinearAlgebraPackage* la{};

public:
   GenMatrixHandle A;
   GenMatrixHandle C;
   OoqpVectorHandle g; // objective
   OoqpVectorHandle bA; // rhs equality
   OoqpVectorHandle bux; // upper bounds x
   OoqpVectorHandle ixupp; // index for upper bounds
   OoqpVectorHandle blx; // lower bounds x
   OoqpVectorHandle ixlow; // index for lower bounds
   OoqpVectorHandle bu; // upper bounds C
   OoqpVectorHandle icupp; // index upper bounds
   OoqpVectorHandle bl; // lower bounds C
   OoqpVectorHandle iclow; // index lower bounds
   OoqpVectorHandle sc; // scale (and diag of Q) -> not maintained currently

   long long nx{0};
   long long my{0};
   long long mz{0};

   long long nxlow{0};
   long long nxupp{0};
   long long mclow{0};
   long long mcupp{0};

   Problem(LinearAlgebraPackage* la, OoqpVector* c, OoqpVector* xlow, OoqpVector* ixlow, OoqpVector* xupp, OoqpVector* ixupp, GenMatrix* A,
         OoqpVector* bA, GenMatrix* C, OoqpVector* clow, OoqpVector* iclow, OoqpVector* cupp, OoqpVector* ciupp);

   virtual ~Problem() = default;

   virtual double objective_value(const Variables& x) const = 0;

   virtual void objective_gradient(const Variables& vars, OoqpVector& gradient) const = 0;

   /** compute the norm of the problem data */
   virtual double datanorm() const;

   /** print the problem data */
   virtual void print();

   OoqpVector& xupperBound() { return *bux; };

   const OoqpVector& xupperBound() const { return *bux; };

   OoqpVector& ixupperBound() { return *ixupp; };

   OoqpVector& xlowerBound() { return *blx; };

   const OoqpVector& xlowerBound() const { return *blx; };

   OoqpVector& ixlowerBound() { return *ixlow; };

   OoqpVector& supperBound() { return *bu; };

   OoqpVector& isupperBound() { return *icupp; };

   OoqpVector& slowerBound() { return *bl; };

   OoqpVector& islowerBound() { return *iclow; };

   OoqpVector& scale() { return *sc; };

   virtual void hessian_multiplication(double beta, OoqpVector& y, double alpha, const OoqpVector& x) const = 0;

   virtual void hessian_diagonal(OoqpVector& hessian_diagonal) = 0;

   /** insert the constraint matrix A into the matrix M for the
    fundamental linear system, where M is stored as a GenMatrix */
   virtual void putAIntoAt(GenMatrix& M, int row, int col);

   /** insert the constraint matrix C into the matrix M for the
       fundamental linear system, where M is stored as a GenMatrix */
   virtual void putCIntoAt(GenMatrix& M, int row, int col);

   /** insert the constraint matrix A into the matrix M for the
       fundamental linear system, where M is stored as a SymMatrix */
   virtual void putAIntoAt(SymMatrix& M, int row, int col);

   /** insert the constraint matrix C into the matrix M for the
       fundamental linear system, where M is stored as a SymMatrix */
   virtual void putCIntoAt(SymMatrix& M, int row, int col);

   /** y = beta * y + alpha * A * x */
   virtual void Amult(double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;

   /** y = beta * y + alpha * C * x   */
   virtual void Cmult(double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;

   /** y = beta * y + alpha * A\T * x */
   virtual void ATransmult(double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;

   /** y = beta * y + alpha * C\T * x */
   virtual void CTransmult(double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;

   void getg(OoqpVector& cout) const;

   void getbA(OoqpVector& bout) const;

   void scaleA();

   void scaleC();

   void scaleg();

   void scalexupp();

   void scalexlow();

   void flipg();
};

#endif
