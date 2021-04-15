/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QuadraticProblem.h"
#include "QpGenVars.h"
#include "DoubleMatrix.h"
#include <cmath>
#include "SimpleVector.h"
#include "LinearAlgebraPackage.h"
#include "MpsReader.h"

QuadraticProblem::QuadraticProblem(LinearAlgebraPackage *la_in, OoqpVector *c_in, SymMatrix *Q_in, OoqpVector *xlow_in,
                                   OoqpVector *ixlow_in,
                                   OoqpVector *xupp_in, OoqpVector *ixupp_in, GenMatrix *A_in, OoqpVector *bA_in,
                                   GenMatrix *C_in, OoqpVector *clow_in, OoqpVector *iclow_in, OoqpVector *cupp_in,
                                   OoqpVector *icupp_in) :
   // superclass constructor
   Problem(la_in, c_in, xlow_in, ixlow_in, xupp_in, ixupp_in, A_in, bA_in, C_in, clow_in, iclow_in, cupp_in, icupp_in) {
   SpReferTo(Q, Q_in);
}

void QuadraticProblem::Qmult(double beta, OoqpVector &y, double alpha, const OoqpVector &x) const {
   Q->mult(beta, y, alpha, x);
}


double QuadraticProblem::datanorm() const {
   double norm = Problem::datanorm();
   double componentNorm;

   componentNorm = Q->abmaxnorm();
   if (componentNorm > norm) norm = componentNorm;
   return norm;
}

void QuadraticProblem::datainput(MpsReader *reader, int &iErr) {
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

void QuadraticProblem::print() {
   std::cout << "begin Q\n";
   Q->writeToStream(std::cout);
   std::cout << "end Q\n";
   Problem::print();
}

void QuadraticProblem::putQIntoAt(SymMatrix &M, int row, int col) {
   M.symAtPutSubmatrix(row, col, *Q, 0, 0, nx, nx);
}

void QuadraticProblem::putQIntoAt(GenMatrix &M, int row, int col) {
   M.atPutSubmatrix(row, col, *Q, 0, 0, nx, nx);
}

void QuadraticProblem::getDiagonalOfQ(OoqpVector &q_diagonal) {
   Q->fromGetDiagonal(0, q_diagonal);
}

void QuadraticProblem::objective_gradient(const QpGenVars *vars, OoqpVector &gradient) {
   this->getg(gradient);
   this->Qmult(1., gradient, 1., *vars->x);
   return;
}

double QuadraticProblem::objective_value(const QpGenVars *vars) const {
   OoqpVectorHandle gradient(la->newVector(nx));
   this->getg(*gradient);
   this->Qmult(1., *gradient, 0.5, *vars->x);

   return gradient->dotProductWith(*vars->x);
}

void QuadraticProblem::createScaleFromQ() {
   // Stuff the diagonal elements of Q into the vector "sc"
   this->getDiagonalOfQ(*sc);

   // Modifying scVector is equivalent to modifying sc
   SimpleVector &scVector = dynamic_cast<SimpleVector &>(*sc);

   int scLength = scVector.length();

   for (int i = 0; i < scLength; i++) {
      if (scVector[i] > 1)
         scVector[i] = 1.0 / sqrt(scVector[i]);
      else
         scVector[i] = 1.0;
   }
}

void QuadraticProblem::scaleQ() {
   Q->symmetricScale(*sc);
}

void QuadraticProblem::flipQ() {
   // Multiply Q matrix by -1
   Q->scalarMult(-1.0);
}