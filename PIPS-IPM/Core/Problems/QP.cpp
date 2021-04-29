/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QP.hpp"
#include "Variables.h"
#include "DoubleMatrix.h"
#include <cmath>
#include "SimpleVector.h"
#include "LinearAlgebraPackage.h"
#include "MpsReader.h"

QP::QP(LinearAlgebraPackage* la_in, OoqpVector* c_in, SymMatrix* Q_in, OoqpVector* xlow_in, OoqpVector* ixlow_in,
      OoqpVector* xupp_in, OoqpVector* ixupp_in, GenMatrix* A_in, OoqpVector* bA_in, GenMatrix* C_in, OoqpVector* clow_in, OoqpVector* iclow_in,
      OoqpVector* cupp_in, OoqpVector* icupp_in) :
// superclass constructor
      Problem(la_in, c_in, xlow_in, ixlow_in, xupp_in, ixupp_in, A_in, bA_in, C_in, clow_in, iclow_in, cupp_in, icupp_in) {
   SpReferTo(Q, Q_in);
}

void QP::hessian_multiplication(double beta, OoqpVector& y, double alpha, const OoqpVector& x) const {
   // y = beta * y + alpha * Q * x
   Q->mult(beta, y, alpha, x);
}

double QP::datanorm() const {
   double norm = Problem::datanorm();
   double componentNorm;

   componentNorm = Q->abmaxnorm();
   if (componentNorm > norm)
      norm = componentNorm;
   return norm;
}

void QP::datainput(MpsReader* reader, int& iErr) {
   reader->readQpGen(*g, *Q, *blx, *ixlow, *bux, *ixupp, *A, *bA, *C, *bl, *iclow, *bu, *icupp, iErr);

   if (reader->scalingOption == 1) {
      // Create the scaling vector
      this->createScaleFromQ();

      //Scale the variables
      this->scaleQ();
      this->scaleA();
      this->scaleC();
      this->scaleg();
      this->scalexlow();
      this->scalexupp();
   }

   /* If objective sense is "MAX", flip the C and Q matrices */
   if (!strncmp(reader->objectiveSense, "MAX", 3)) {
      this->flipg();
      this->flipQ();
   }
}

void QP::print() {
   std::cout << "begin Q\n";
   Q->writeToStream(std::cout);
   std::cout << "end Q\n";
   Problem::print();
}

void QP::putQIntoAt(SymMatrix& M, int row, int col) {
   M.symAtPutSubmatrix(row, col, *Q, 0, 0, nx, nx);
}

void QP::putQIntoAt(GenMatrix& M, int row, int col) {
   M.atPutSubmatrix(row, col, *Q, 0, 0, nx, nx);
}

void QP::hessian_diagonal(OoqpVector& hessian_diagonal) {
   Q->fromGetDiagonal(0, hessian_diagonal);
}

void QP::objective_gradient(const Variables& variables, OoqpVector& gradient) const {
   this->getg(gradient);
   this->hessian_multiplication(1., gradient, 1., *variables.x);
   return;
}

double QP::objective_value(const Variables& variables) const {
   OoqpVectorHandle gradient(la->newVector(nx));
   this->getg(*gradient);
   this->hessian_multiplication(1., *gradient, 0.5, *variables.x);

   return gradient->dotProductWith(*variables.x);
}

void QP::createScaleFromQ() {
   // Stuff the diagonal elements of Q into the vector "sc"
   this->hessian_diagonal(*sc);

   // Modifying scVector is equivalent to modifying sc
   SimpleVector<double>& scVector = dynamic_cast<SimpleVector<double>&>(*sc);

   int scLength = scVector.length();

   for (int i = 0; i < scLength; i++) {
      if (scVector[i] > 1)
         scVector[i] = 1.0 / sqrt(scVector[i]);
      else
         scVector[i] = 1.0;
   }
}

void QP::scaleQ() {
   Q->symmetricScale(*sc);
}

void QP::flipQ() {
   // Multiply Q matrix by -1
   Q->scalarMult(-1.0);
}