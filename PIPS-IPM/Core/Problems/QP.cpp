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

QP::QP(LinearAlgebraPackage* la_in, Vector<double>* c_in, SymMatrix* Q_in, Vector<double>* xlow_in, Vector<double>* ixlow_in,
      Vector<double>* xupp_in, Vector<double>* ixupp_in, GenMatrix* A_in, Vector<double>* bA_in, GenMatrix* C_in, Vector<double>* clow_in, Vector<double>* iclow_in,
      Vector<double>* cupp_in, Vector<double>* icupp_in) :
// superclass constructor
      Problem(la_in, c_in, xlow_in, ixlow_in, xupp_in, ixupp_in, A_in, bA_in, C_in, clow_in, iclow_in, cupp_in, icupp_in) {
   SpReferTo(Q, Q_in);
}

void QP::hessian_multiplication(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
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

void QP::hessian_diagonal(Vector<double>& hessian_diagonal) {
   Q->fromGetDiagonal(0, hessian_diagonal);
}

void QP::objective_gradient(const Variables& variables, Vector<double>& gradient) const {
   this->getg(gradient);
   this->hessian_multiplication(1., gradient, 1., *variables.x);
   return;
}

double QP::objective_value(const Variables& variables) const {
   SmartPointer<Vector<double> > gradient(la->newVector(nx));
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