/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef PROBLEM_H
#define PROBLEM_H

#include "Variables.h"
#include "Vector.hpp"
#include "AbstractMatrix.h"

class MpsReader;

class Variables;

class Problem {
protected:
   Problem() = default;

public:
   std::shared_ptr<SymmetricMatrix> hessian;
   std::shared_ptr<GeneralMatrix> equality_jacobian;
   std::shared_ptr<GeneralMatrix> inequality_jacobian;
   std::shared_ptr<Vector<double>> objective_gradient; // objective
   std::shared_ptr<Vector<double>> equality_rhs; // rhs equality
   std::shared_ptr<Vector<double>> primal_upper_bounds; // upper bounds x
   std::shared_ptr<Vector<double>> primal_upper_bound_indicators; // index for upper bounds
   std::shared_ptr<Vector<double>> primal_lower_bounds; // lower bounds x
   std::shared_ptr<Vector<double>> primal_lower_bound_indicators; // index for lower bounds
   std::shared_ptr<Vector<double>> inequality_upper_bounds; // upper bounds C
   std::shared_ptr<Vector<double>> inequality_upper_bound_indicators; // index upper bounds
   std::shared_ptr<Vector<double>> inequality_lower_bounds; // lower bounds C
   std::shared_ptr<Vector<double>> inequality_lower_bound_indicators; // index lower bounds
   std::shared_ptr<Vector<double>> sc; // scale (and diag of Q) -> not maintained currently

   long long nx{0};
   long long my{0};
   long long mz{0};

   long long number_primal_lower_bounds{0};
   long long number_primal_upper_bounds{0};
   long long number_inequality_lower_bounds{0};
   long long number_inequality_upper_bounds{0};

   Problem(std::shared_ptr<Vector<double>> g_in, std::shared_ptr<SymmetricMatrix> Q_in,
         std::shared_ptr<Vector<double>> xlow_in,
      std::shared_ptr<Vector<double>> ixlow_in, std::shared_ptr<Vector<double>> xupp_in, std::shared_ptr<Vector<double>> ixupp_in,
      std::shared_ptr<GeneralMatrix> A_in, std::shared_ptr<Vector<double>> bA_in, std::shared_ptr<GeneralMatrix> C_in,
      std::shared_ptr<Vector<double>> clow_in, std::shared_ptr<Vector<double>> iclow_in,
      std::shared_ptr<Vector<double>> cupp_in, std::shared_ptr<Vector<double>> icupp_in);

   virtual ~Problem() = default;

   [[nodiscard]] virtual double evaluate_objective(const Variables& x) const;

   virtual void evaluate_objective_gradient(const Variables& vars, Vector<double>& gradient) const;

   /** compute the norm of the problem data */
   [[nodiscard]] virtual double datanorm() const;

   /** print the problem data */
   virtual void print();

   [[nodiscard]] Vector<double>& x_lower_bound() const { return *primal_lower_bounds; };
   [[nodiscard]] Vector<double>& x_upper_bound() const { return *primal_upper_bounds; };
   [[nodiscard]] Vector<double>& s_lower_bound() const { return *inequality_lower_bounds; };
   [[nodiscard]] Vector<double>& s_upper_bound() const { return *inequality_upper_bounds; };

   [[nodiscard]] Vector<double>& has_x_lower_bound() const { return *primal_lower_bound_indicators; };
   [[nodiscard]] Vector<double>& has_x_upper_bound() const { return *primal_upper_bound_indicators; };
   [[nodiscard]] Vector<double>& has_s_lower_bound() const { return *inequality_lower_bound_indicators; };
   [[nodiscard]] Vector<double>& has_s_upper_bound() const { return *inequality_upper_bound_indicators; };

   Vector<double>& scale() { return *sc; };

   [[nodiscard]] virtual std::unique_ptr<Problem> clone_full() const = 0;
   virtual void write_to_streamDense(std::ostream& out) const = 0;

   /** y = beta * y + alpha * Q * x */
   void hessian_multiplication(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const;

   /** extract the diagonal of the Hessian and put it in the Vector<double> hessian_diagonal */
   virtual void hessian_diagonal(Vector<double>& hessian_diagonal) const;

   /** y = beta * y + alpha * A * x */
   virtual void Amult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const;

   /** y = beta * y + alpha * C * x   */
   virtual void Cmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const;

   /** y = beta * y + alpha * A\T * x */
   virtual void ATransmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const;

   /** y = beta * y + alpha * C\T * x */
   virtual void CTransmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const;

   void get_objective_gradient(Vector<double>& myG) const;

   void getbA(Vector<double>& bout) const;

   void flip_objective_gradient();

   void flip_hessian();

   virtual void datainput(MpsReader* reader, int& iErr);

   void print_ranges() const;
};

#endif
