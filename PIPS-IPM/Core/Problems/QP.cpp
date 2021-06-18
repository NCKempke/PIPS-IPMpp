/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QP.hpp"
#include "Variables.h"
#include "AbstractMatrix.h"
#include <cmath>
#include "SimpleVector.hpp"
#include "MpsReader.h"

QP::QP(std::shared_ptr<Vector<double>> c_in, std::shared_ptr<SymmetricMatrix> Q_in,
   std::shared_ptr<Vector<double>> xlow_in,
   std::shared_ptr<Vector<double>> ixlow_in, std::shared_ptr<Vector<double>> xupp_in,
   std::shared_ptr<Vector<double>> ixupp_in, std::shared_ptr<GeneralMatrix> A_in, std::shared_ptr<Vector<double>> bA_in,
   std::shared_ptr<GeneralMatrix> C_in, std::shared_ptr<Vector<double>> clow_in,
   std::shared_ptr<Vector<double>> iclow_in,
   std::shared_ptr<Vector<double>> cupp_in,
   std::shared_ptr<Vector<double>> icupp_in) :
   Problem(std::move(c_in), std::move(xlow_in), std::move(ixlow_in), std::move(xupp_in), std::move(ixupp_in),
      std::move(A_in), std::move(bA_in), std::move(C_in), std::move(clow_in), std::move(iclow_in), std::move(cupp_in),
      std::move(icupp_in)) {
   Q = std::move(Q_in);
}

void QP::hessian_multiplication(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   // y = beta * y + alpha * Q * x
   Q->mult(beta, y, alpha, x);
}

double QP::datanorm() const {
   return std::max(Problem::datanorm(), Q->inf_norm());
}

void QP::datainput(MpsReader* reader, int& iErr) {
   reader->readQpGen(*objective_gradient, *Q, *primal_lower_bounds, *primal_lower_bound_indicators, *primal_upper_bounds, *primal_upper_bound_indicators, *equality_jacobian, *equality_rhs, *inequality_jacobian, *inequality_lower_bounds, *inequality_lower_bound_indicators, *inequality_upper_bounds, *inequality_upper_bound_indicators, iErr);

   if (reader->scalingOption == 1) {
      // Create the scaling vector
      this->create_scale_from_hessian();

      //Scale the variables
      this->scale_hessian();
      this->scaleA();
      this->scaleC();
      this->scaleg();
      this->scalexlow();
      this->scalexupp();
   }

   /* If objective sense is "MAX", flip the C and Q matrices */
   if (!strncmp(reader->objectiveSense, "MAX", 3)) {
      this->flipg();
      this->flip_hessian();
   }
}

void QP::print() {
   std::cout << "begin Q\n";
   Q->write_to_stream(std::cout);
   std::cout << "end Q\n";
   Problem::print();
}

void QP::putQIntoAt(SymmetricMatrix& M, int row, int col) {
   M.symAtPutSubmatrix(row, col, *Q, 0, 0, nx, nx);
}

void QP::putQIntoAt(GeneralMatrix& M, int row, int col) {
   M.atPutSubmatrix(row, col, *Q, 0, 0, nx, nx);
}

void QP::hessian_diagonal(Vector<double>& hessian_diagonal) const {
   Q->fromGetDiagonal(0, hessian_diagonal);
}

void QP::evaluate_objective_gradient(const Variables& variables, Vector<double>& gradient) const {
   this->get_objective_gradient(gradient);
   this->hessian_multiplication(1., gradient, 1., *variables.primals);
}

double QP::evaluate_objective(const Variables& variables) const {
   SimpleVector<double> gradient(nx);
   this->get_objective_gradient(gradient);
   this->hessian_multiplication(1., gradient, 0.5, *variables.primals);
   return gradient.dotProductWith(*variables.primals);
}

void QP::create_scale_from_hessian() {
   // Stuff the diagonal elements of Q into the vector "sc"
   this->hessian_diagonal(*sc);

   // Modifying scVector is equivalent to modifying sc
   auto& scVector = dynamic_cast<SimpleVector<double>&>(*sc);
   for (int i = 0; i < scVector.length(); i++) {
      if (scVector[i] > 1)
         scVector[i] = 1.0 / sqrt(scVector[i]);
      else
         scVector[i] = 1.0;
   }
}

void QP::scale_hessian() {
   Q->symmetricScale(*sc);
}

void QP::flip_hessian() {
   // Multiply Q matrix by -1
   Q->scalarMult(-1.0);
}