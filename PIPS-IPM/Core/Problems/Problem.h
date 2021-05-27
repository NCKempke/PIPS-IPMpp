/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef PROBLEM_H
#define PROBLEM_H

#include "Variables.h"
#include "Vector.hpp"
#include "SimpleVector.h"
#include "AbstractMatrix.h"

class Problem {
protected:
   Problem() = default;

public:
   std::shared_ptr<GeneralMatrix> A;
   std::shared_ptr<GeneralMatrix> C;
   std::shared_ptr<Vector<double>> g; // objective
   std::shared_ptr<Vector<double>> bA; // rhs equality
   std::shared_ptr<Vector<double>> bux; // upper bounds x
   std::shared_ptr<Vector<double>> ixupp; // index for upper bounds
   std::shared_ptr<Vector<double>> blx; // lower bounds x
   std::shared_ptr<Vector<double>> ixlow; // index for lower bounds
   std::shared_ptr<Vector<double>> bu; // upper bounds C
   std::shared_ptr<Vector<double>> icupp; // index upper bounds
   std::shared_ptr<Vector<double>> bl; // lower bounds C
   std::shared_ptr<Vector<double>> iclow; // index lower bounds
   std::shared_ptr<Vector<double>> sc; // scale (and diag of Q) -> not maintained currently

   long long nx{0};
   long long my{0};
   long long mz{0};

   long long nxlow{0};
   long long nxupp{0};
   long long mclow{0};
   long long mcupp{0};

   Problem(std::shared_ptr<Vector<double>> g_in, std::shared_ptr<Vector<double>> xlow_in,
      std::shared_ptr<Vector<double>> ixlow_in, std::shared_ptr<Vector<double>> xupp_in, std::shared_ptr<Vector<double>> ixupp_in,
      std::shared_ptr<GeneralMatrix> A_in, std::shared_ptr<Vector<double>> bA_in, std::shared_ptr<GeneralMatrix> C_in,
      std::shared_ptr<Vector<double>> clow_in, std::shared_ptr<Vector<double>> iclow_in,
      std::shared_ptr<Vector<double>> cupp_in, std::shared_ptr<Vector<double>> icupp_in);

   virtual ~Problem() = default;

   [[nodiscard]] virtual double objective_value(const Variables& x) const = 0;

   virtual void objective_gradient(const Variables& vars, Vector<double>& gradient) const = 0;

   /** compute the norm of the problem data */
   [[nodiscard]] virtual double datanorm() const;

   /** print the problem data */
   virtual void print();

   Vector<double>& xupperBound() { return *bux; };

   [[nodiscard]] const Vector<double>& xupperBound() const { return *bux; };

   Vector<double>& ixupperBound() { return *ixupp; };

   Vector<double>& xlowerBound() { return *blx; };

   [[nodiscard]] const Vector<double>& xlowerBound() const { return *blx; };

   Vector<double>& ixlowerBound() { return *ixlow; };

   Vector<double>& supperBound() { return *bu; };

   Vector<double>& isupperBound() { return *icupp; };

   Vector<double>& slowerBound() { return *bl; };

   Vector<double>& islowerBound() { return *iclow; };

   Vector<double>& scale() { return *sc; };

   virtual void hessian_multiplication(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const = 0;

   virtual void hessian_diagonal(Vector<double>& hessian_diagonal) const = 0;

   /** insert the constraint matrix A into the matrix M for the
    fundamental linear system, where M is stored as a GenMatrix */
   virtual void putAIntoAt(GeneralMatrix& M, int row, int col);

   /** insert the constraint matrix C into the matrix M for the
       fundamental linear system, where M is stored as a GenMatrix */
   virtual void putCIntoAt(GeneralMatrix& M, int row, int col);

   /** insert the constraint matrix A into the matrix M for the
       fundamental linear system, where M is stored as a SymMatrix */
   virtual void putAIntoAt(SymmetricMatrix& M, int row, int col);

   /** insert the constraint matrix C into the matrix M for the
       fundamental linear system, where M is stored as a SymMatrix */
   virtual void putCIntoAt(SymmetricMatrix& M, int row, int col);

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
