#include <SimpleVector.h>
#include "Problem.h"

Problem::Problem(LinearAlgebraPackage *la_in, OoqpVector *c_in, OoqpVector *xlow_in, OoqpVector *ixlow_in,
   OoqpVector *xupp_in, OoqpVector *ixupp_in, GenMatrix *A_in, OoqpVector *bA_in, GenMatrix *C_in, OoqpVector *clow_in,
   OoqpVector *iclow_in, OoqpVector *cupp_in, OoqpVector *icupp_in) :
      la{la_in},
      nxlow{ixlow_in->numberOfNonzeros()},
      nxupp{ixupp_in->numberOfNonzeros()},
      mclow{iclow_in->numberOfNonzeros()},
      mcupp{icupp_in->numberOfNonzeros()} {
   SpReferTo(g, c_in);
   SpReferTo(bA, bA_in);
   SpReferTo(blx, xlow_in);
   SpReferTo(ixlow, ixlow_in);
   SpReferTo(bux, xupp_in);
   SpReferTo(ixupp, ixupp_in);
   SpReferTo(bl, clow_in);
   SpReferTo(iclow, iclow_in);
   SpReferTo(bu, cupp_in);
   SpReferTo(icupp, icupp_in);

   long long dummy;

   nx = g->length();

   SpReferTo(A, A_in);
   A->getSize(my, dummy);

   SpReferTo(C, C_in);
   C->getSize(mz, dummy);
}

void Problem::Amult(double beta, OoqpVector &y, double alpha, const OoqpVector &x) const {
   A->mult(beta, y, alpha, x);
}

void Problem::Cmult(double beta, OoqpVector &y, double alpha, const OoqpVector &x) const {
   C->mult(beta, y, alpha, x);
}

void Problem::ATransmult(double beta, OoqpVector &y, double alpha, const OoqpVector &x) const {
   A->transMult(beta, y, alpha, x);
}

void Problem::CTransmult(double beta, OoqpVector &y, double alpha, const OoqpVector &x) const {
   C->transMult(beta, y, alpha, x);
}

void Problem::getg(OoqpVector &myG) const {
   myG.copyFrom(*g);
}

void Problem::getbA(OoqpVector &bout) const {
   bout.copyFrom(*bA);
}

// precondition: hessian_diagonal is a vector of zeros
void Problem::hessian_diagonal(OoqpVector& hessian_diagonal) {
   // if there is no Hessian, do nothing
   return;
}

void Problem::putAIntoAt(GenMatrix &M, int row, int col) {
   M.atPutSubmatrix(row, col, *A, 0, 0, my, nx);
}

void Problem::putAIntoAt(SymMatrix &M, int row, int col) {
   M.symAtPutSubmatrix(row, col, *A, 0, 0, my, nx);
}

void Problem::putCIntoAt(GenMatrix &M, int row, int col) {
   M.atPutSubmatrix(row, col, *C, 0, 0, mz, nx);
}

void Problem::putCIntoAt(SymMatrix &M, int row, int col) {
   M.symAtPutSubmatrix(row, col, *C, 0, 0, mz, nx);
}

void Problem::scaleA() {
   A->columnScale(*sc);
}

void Problem::scaleC() {
   C->columnScale(*sc);
}

void Problem::scaleg() {
   SimpleVector &scVector = dynamic_cast<SimpleVector &>(*sc);
   assert (scVector.length() == g->length());

   // D * g
   g->componentMult(scVector);
}

void Problem::scalexupp() {
   SimpleVector &scVector = dynamic_cast<SimpleVector &>(*sc);

   assert (scVector.length() == bux->length());

   // inverse(D) * bux
   bux->componentDiv(scVector);

}


void Problem::scalexlow() {
   SimpleVector &scVector = dynamic_cast<SimpleVector &>(*sc);

   assert (scVector.length() == blx->length());

   // inverse(D) * blx
   blx->componentDiv(scVector);

}

void Problem::flipg() {
   // Multiply C matrix by -1
   g->scalarMult(-1.0);
}

double Problem::datanorm() const {
   double norm = 0.0;
   double componentNorm;

   componentNorm = g->infnorm();
   if (componentNorm > norm) norm = componentNorm;

   componentNorm = bA->infnorm();
   if (componentNorm > norm) norm = componentNorm;

   componentNorm = A->abmaxnorm();
   if (componentNorm > norm) norm = componentNorm;

   componentNorm = C->abmaxnorm();
   if (componentNorm > norm) norm = componentNorm;

   assert(blx->matchesNonZeroPattern(*ixlow));
   componentNorm = blx->infnorm();
   if (componentNorm > norm) norm = componentNorm;

   assert(bux->matchesNonZeroPattern(*ixupp));
   componentNorm = bux->infnorm();
   if (componentNorm > norm) norm = componentNorm;

   assert(bl->matchesNonZeroPattern(*iclow));
   componentNorm = bl->infnorm();
   if (componentNorm > norm) norm = componentNorm;

   assert(bu->matchesNonZeroPattern(*icupp));
   componentNorm = bu->infnorm();
   if (componentNorm > norm) norm = componentNorm;

   return norm;
}

void Problem::print() {
   std::cout << "begin c\n";
   g->writeToStream(std::cout);
   std::cout << "end c\n";

   std::cout << "begin xlow\n";
   blx->writeToStream(std::cout);
   std::cout << "end xlow\n";
   std::cout << "begin ixlow\n";
   ixlow->writeToStream(std::cout);
   std::cout << "end ixlow\n";

   std::cout << "begin xupp\n";
   bux->writeToStream(std::cout);
   std::cout << "end xupp\n";
   std::cout << "begin ixupp\n";
   ixupp->writeToStream(std::cout);
   std::cout << "end ixupp\n";
   std::cout << "begin A\n";

   A->writeToStream(std::cout);
   std::cout << "end A\n";
   std::cout << "begin b\n";
   bA->writeToStream(std::cout);
   std::cout << "end b\n";
   std::cout << "begin C\n";
   C->writeToStream(std::cout);
   std::cout << "end C\n";

   std::cout << "begin clow\n";
   bl->writeToStream(std::cout);
   std::cout << "end clow\n";
   std::cout << "begin iclow\n";
   iclow->writeToStream(std::cout);
   std::cout << "end iclow\n";

   std::cout << "begin cupp\n";
   bu->writeToStream(std::cout);
   std::cout << "end cupp\n";
   std::cout << "begin icupp\n";
   icupp->writeToStream(std::cout);
   std::cout << "end icupp\n";
}