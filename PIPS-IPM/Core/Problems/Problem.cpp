#include "Problem.hpp"
#include "SimpleVector.hpp"
#include "Variables.h"
#include "AbstractMatrix.h"
#include <cmath>
#include "MpsReader.h"

Problem::Problem(std::shared_ptr<Vector<double>> g_in, std::shared_ptr<SymmetricMatrix> Q_in,
      std::shared_ptr<Vector<double>> xlow_in,
   std::shared_ptr<Vector<double>> ixlow_in, std::shared_ptr<Vector<double>> xupp_in,
   std::shared_ptr<Vector<double>> ixupp_in,
   std::shared_ptr<GeneralMatrix> A_in, std::shared_ptr<Vector<double>> bA_in, std::shared_ptr<GeneralMatrix> C_in,
   std::shared_ptr<Vector<double>> clow_in, std::shared_ptr<Vector<double>> iclow_in,
   std::shared_ptr<Vector<double>> cupp_in, std::shared_ptr<Vector<double>> icupp_in) :
   hessian(std::move(Q_in)),
   equality_jacobian{std::move(A_in)},
   inequality_jacobian{std::move(C_in)},
   objective_gradient{std::move(g_in)}, equality_rhs{std::move(bA_in)}, primal_upper_bounds{std::move(xupp_in)}, primal_upper_bound_indicators{std::move(ixupp_in)},
   primal_lower_bounds{std::move(xlow_in)},
   primal_lower_bound_indicators{std::move(ixlow_in)}, inequality_upper_bounds{std::move(cupp_in)}, inequality_upper_bound_indicators{std::move(icupp_in)}, inequality_lower_bounds{std::move(clow_in)},
   inequality_lower_bound_indicators{std::move(iclow_in)},
   nx{objective_gradient->length()}, my{equality_jacobian->n_rows()}, mz{inequality_jacobian->n_rows()}, number_primal_lower_bounds{primal_lower_bound_indicators->number_nonzeros()},
   number_primal_upper_bounds{primal_upper_bound_indicators->number_nonzeros()},
   number_inequality_lower_bounds{inequality_lower_bound_indicators->number_nonzeros()}, number_inequality_upper_bounds{inequality_upper_bound_indicators->number_nonzeros()} {
   assert(primal_lower_bound_indicators && primal_upper_bound_indicators && inequality_lower_bound_indicators && inequality_upper_bound_indicators);
}

void Problem::Amult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   equality_jacobian->mult(beta, y, alpha, x);
}

void Problem::Cmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   inequality_jacobian->mult(beta, y, alpha, x);
}

void Problem::ATransmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   equality_jacobian->transMult(beta, y, alpha, x);
}

void Problem::CTransmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   inequality_jacobian->transMult(beta, y, alpha, x);
}

void Problem::get_objective_gradient(Vector<double>& myG) const {
   myG.copyFrom(*objective_gradient);
}

void Problem::getbA(Vector<double>& bout) const {
   bout.copyFrom(*equality_rhs);
}

void Problem::putAIntoAt(GeneralMatrix& M, int row, int col) {
   M.atPutSubmatrix(row, col, *equality_jacobian, 0, 0, my, nx);
}

void Problem::putAIntoAt(SymmetricMatrix& M, int row, int col) {
   M.symAtPutSubmatrix(row, col, *equality_jacobian, 0, 0, my, nx);
}

void Problem::putCIntoAt(GeneralMatrix& M, int row, int col) {
   M.atPutSubmatrix(row, col, *inequality_jacobian, 0, 0, mz, nx);
}

void Problem::putCIntoAt(SymmetricMatrix& M, int row, int col) {
   M.symAtPutSubmatrix(row, col, *inequality_jacobian, 0, 0, mz, nx);
}

void Problem::scaleA() {
   equality_jacobian->columnScale(*sc);
}

void Problem::scaleC() {
   inequality_jacobian->columnScale(*sc);
}

void Problem::scaleg() {
   auto& scVector = dynamic_cast<SimpleVector<double>&>(*sc);
   assert (scVector.length() == objective_gradient->length());

   // D * g
   objective_gradient->componentMult(scVector);
}

void Problem::scalexupp() {
   auto& scVector = dynamic_cast<SimpleVector<double>&>(*sc);

   assert (scVector.length() == primal_upper_bounds->length());

   // inverse(D) * bux
   primal_upper_bounds->componentDiv(scVector);
}


void Problem::scalexlow() {
   auto& scVector = dynamic_cast<SimpleVector<double>&>(*sc);

   assert (scVector.length() == primal_lower_bounds->length());

   // inverse(D) * blx
   primal_lower_bounds->componentDiv(scVector);
}

void Problem::flip_objective_gradient() {
   // Multiply g by -1
   objective_gradient->scalarMult(-1.0);
}

double Problem::datanorm() const {
   double norm = objective_gradient->inf_norm();

   norm = std::max(norm, equality_rhs->inf_norm());

   norm = std::max(norm, equality_jacobian->inf_norm());

   norm = std::max(norm, inequality_jacobian->inf_norm());

   assert(primal_lower_bounds->matchesNonZeroPattern(*primal_lower_bound_indicators));
   norm = std::max(norm, primal_lower_bounds->inf_norm());

   assert(primal_upper_bounds->matchesNonZeroPattern(*primal_upper_bound_indicators));
   norm = std::max(norm, primal_upper_bounds->inf_norm());

   assert(inequality_lower_bounds->matchesNonZeroPattern(*inequality_lower_bound_indicators));
   norm = std::max(norm, inequality_lower_bounds->inf_norm());

   assert(inequality_upper_bounds->matchesNonZeroPattern(*inequality_upper_bound_indicators));
   norm = std::max(norm, inequality_upper_bounds->inf_norm());

   norm = std::max(norm, hessian->inf_norm());

   return norm;
}

void Problem::print() {
   std::cout << "begin Q\n";
   hessian->write_to_stream(std::cout);
   std::cout << "end Q\n";

   std::cout << "begin c\n";
   objective_gradient->write_to_stream(std::cout);
   std::cout << "end c\n";

   std::cout << "begin xlow\n";
   primal_lower_bounds->write_to_stream(std::cout);
   std::cout << "end xlow\n";
   std::cout << "begin ixlow\n";
   primal_lower_bound_indicators->write_to_stream(std::cout);
   std::cout << "end ixlow\n";

   std::cout << "begin xupp\n";
   primal_upper_bounds->write_to_stream(std::cout);
   std::cout << "end xupp\n";
   std::cout << "begin ixupp\n";
   primal_upper_bound_indicators->write_to_stream(std::cout);
   std::cout << "end ixupp\n";
   std::cout << "begin A\n";

   equality_jacobian->write_to_stream(std::cout);
   std::cout << "end A\n";
   std::cout << "begin b\n";
   equality_rhs->write_to_stream(std::cout);
   std::cout << "end b\n";
   std::cout << "begin C\n";
   inequality_jacobian->write_to_stream(std::cout);
   std::cout << "end C\n";

   std::cout << "begin clow\n";
   inequality_lower_bounds->write_to_stream(std::cout);
   std::cout << "end clow\n";
   std::cout << "begin iclow\n";
   inequality_lower_bound_indicators->write_to_stream(std::cout);
   std::cout << "end iclow\n";

   std::cout << "begin cupp\n";
   inequality_upper_bounds->write_to_stream(std::cout);
   std::cout << "end cupp\n";
   std::cout << "begin icupp\n";
   inequality_upper_bound_indicators->write_to_stream(std::cout);
   std::cout << "end icupp\n";
}

void Problem::hessian_multiplication(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   // y = beta * y + alpha * Q * x
   hessian->mult(beta, y, alpha, x);
}

void Problem::datainput(MpsReader* reader, int& iErr) {
   reader->readQpGen(*objective_gradient, *hessian, *primal_lower_bounds, *primal_lower_bound_indicators, *primal_upper_bounds,
         *primal_upper_bound_indicators, *equality_jacobian, *equality_rhs, *inequality_jacobian, *inequality_lower_bounds,
         *inequality_lower_bound_indicators, *inequality_upper_bounds, *inequality_upper_bound_indicators, iErr);

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
      this->flip_objective_gradient();
      this->flip_hessian();
   }
}

void Problem::put_hessian_into_At(SymmetricMatrix& M, int row, int col) {
   M.symAtPutSubmatrix(row, col, *hessian, 0, 0, nx, nx);
}

void Problem::put_hessian_into_At(GeneralMatrix& M, int row, int col) {
   M.atPutSubmatrix(row, col, *hessian, 0, 0, nx, nx);
}

void Problem::hessian_diagonal(Vector<double>& hessian_diagonal) const {
   hessian->fromGetDiagonal(0, hessian_diagonal);
}

void Problem::evaluate_objective_gradient(const Variables& variables, Vector<double>& gradient) const {
   this->get_objective_gradient(gradient);
   this->hessian_multiplication(1., gradient, 1., *variables.primals);
}

double Problem::evaluate_objective(const Variables& variables) const {
   SimpleVector<double> gradient(nx);
   this->get_objective_gradient(gradient);
   this->hessian_multiplication(1., gradient, 0.5, *variables.primals);
   return gradient.dotProductWith(*variables.primals);
}

void Problem::create_scale_from_hessian() {
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

void Problem::scale_hessian() {
   hessian->symmetricScale(*sc);
}

void Problem::flip_hessian() {
   // Multiply Q matrix by -1
   hessian->scalarMult(-1.0);
}