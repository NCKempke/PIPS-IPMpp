/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef PROBLEM_H
#define PROBLEM_H

#include "Variables.h"
#include "Vector.hpp"
#include "SmartPointer.h"
#include "DoubleMatrix.h"

class LinearAlgebraPackage;

class Problem {
protected:
   Problem() = default;

   LinearAlgebraPackage* la{};

public:
   SmartPointer<GenMatrix> A;
   SmartPointer<GenMatrix> C;
   SmartPointer<Vector<double> > g; // objective
   SmartPointer<Vector<double> > bA; // rhs equality
   SmartPointer<Vector<double> > bux; // upper bounds x
   SmartPointer<Vector<double> > ixupp; // index for upper bounds
   SmartPointer<Vector<double> > blx; // lower bounds x
   SmartPointer<Vector<double> > ixlow; // index for lower bounds
   SmartPointer<Vector<double> > bu; // upper bounds C
   SmartPointer<Vector<double> > icupp; // index upper bounds
   SmartPointer<Vector<double> > bl; // lower bounds C
   SmartPointer<Vector<double> > iclow; // index lower bounds
   SmartPointer<Vector<double> > sc; // scale (and diag of Q) -> not maintained currently

   long long nx{0};
   long long my{0};
   long long mz{0};

   long long nxlow{0};
   long long nxupp{0};
   long long mclow{0};
   long long mcupp{0};

   Problem(LinearAlgebraPackage* la, Vector<double>* c, Vector<double>* xlow, Vector<double>* ixlow, Vector<double>* xupp, Vector<double>* ixupp,
         GenMatrix* A, Vector<double>* bA, GenMatrix* C, Vector<double>* clow, Vector<double>* iclow, Vector<double>* cupp, Vector<double>* ciupp);

   virtual ~Problem() = default;

   virtual double objective_value(const Variables& x) const = 0;

   virtual void objective_gradient(const Variables& vars, Vector<double>& gradient) const = 0;

   /** compute the norm of the problem data */
   virtual double datanorm() const;

   /** print the problem data */
   virtual void print();

   Vector<double>& xupperBound() { return *bux; };

   const Vector<double>& xupperBound() const { return *bux; };

   Vector<double>& ixupperBound() { return *ixupp; };

   Vector<double>& xlowerBound() { return *blx; };

   const Vector<double>& xlowerBound() const { return *blx; };

   Vector<double>& ixlowerBound() { return *ixlow; };

   Vector<double>& supperBound() { return *bu; };

   Vector<double>& isupperBound() { return *icupp; };

   Vector<double>& slowerBound() { return *bl; };

   Vector<double>& islowerBound() { return *iclow; };

   Vector<double>& scale() { return *sc; };

   virtual void hessian_multiplication(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const = 0;

   virtual void hessian_diagonal(Vector<double>& hessian_diagonal) = 0;

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
   virtual void Amult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const;

   /** y = beta * y + alpha * C * x   */
   virtual void Cmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const;

   /** y = beta * y + alpha * A\T * x */
   virtual void ATransmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const;

   /** y = beta * y + alpha * C\T * x */
   virtual void CTransmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const;

   void getg(Vector<double>& cout) const;

   void getbA(Vector<double>& bout) const;

   void scaleA();

   void scaleC();

   void scaleg();

   void scalexupp();

   void scalexlow();

   void flipg();
};

#endif
