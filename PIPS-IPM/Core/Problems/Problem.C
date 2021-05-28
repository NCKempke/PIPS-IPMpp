#include <SimpleVector.h>
#include "Problem.h"

Problem::Problem(std::shared_ptr<Vector<double>> g_in, std::shared_ptr<Vector<double>> xlow_in,
   std::shared_ptr<Vector<double>> ixlow_in, std::shared_ptr<Vector<double>> xupp_in,
   std::shared_ptr<Vector<double>> ixupp_in,
   std::shared_ptr<GeneralMatrix> A_in, std::shared_ptr<Vector<double>> bA_in, std::shared_ptr<GeneralMatrix> C_in,
   std::shared_ptr<Vector<double>> clow_in, std::shared_ptr<Vector<double>> iclow_in,
   std::shared_ptr<Vector<double>> cupp_in, std::shared_ptr<Vector<double>> icupp_in) : A{std::move(A_in)},
   C{std::move(C_in)},
   g{std::move(g_in)}, bA{std::move(bA_in)}, bux{std::move(xupp_in)}, ixupp{std::move(ixupp_in)},
   blx{std::move(xlow_in)},
   ixlow{std::move(ixlow_in)}, bu{std::move(cupp_in)}, icupp{std::move(icupp_in)}, bl{std::move(clow_in)},
   iclow{std::move(iclow_in)},
   nx{g->length()}, my{A->n_rows()}, mz{C->n_rows()}, nxlow{ixlow->number_nonzeros()},
   nxupp{ixupp->number_nonzeros()},
   mclow{iclow->number_nonzeros()}, mcupp{icupp->number_nonzeros()} {
   assert(ixlow && ixupp && iclow && icupp);
}

void Problem::Amult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   A->mult(beta, y, alpha, x);
}

void Problem::Cmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   C->mult(beta, y, alpha, x);
}

void Problem::ATransmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   A->transMult(beta, y, alpha, x);
}

void Problem::CTransmult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   C->transMult(beta, y, alpha, x);
}

void Problem::getg(Vector<double>& myG) const {
   myG.copyFrom(*g);
}

void Problem::getbA(Vector<double>& bout) const {
   bout.copyFrom(*bA);
}

void Problem::putAIntoAt(GeneralMatrix& M, int row, int col) {
   M.atPutSubmatrix(row, col, *A, 0, 0, my, nx);
}

void Problem::putAIntoAt(SymmetricMatrix& M, int row, int col) {
   M.symAtPutSubmatrix(row, col, *A, 0, 0, my, nx);
}

void Problem::putCIntoAt(GeneralMatrix& M, int row, int col) {
   M.atPutSubmatrix(row, col, *C, 0, 0, mz, nx);
}

void Problem::putCIntoAt(SymmetricMatrix& M, int row, int col) {
   M.symAtPutSubmatrix(row, col, *C, 0, 0, mz, nx);
}

void Problem::scaleA() {
   A->columnScale(*sc);
}

void Problem::scaleC() {
   C->columnScale(*sc);
}

void Problem::scaleg() {
   auto& scVector = dynamic_cast<SimpleVector<double>&>(*sc);
   assert (scVector.length() == g->length());

   // D * g
   g->componentMult(scVector);
}

void Problem::scalexupp() {
   auto& scVector = dynamic_cast<SimpleVector<double>&>(*sc);

   assert (scVector.length() == bux->length());

   // inverse(D) * bux
   bux->componentDiv(scVector);

}


void Problem::scalexlow() {
   auto& scVector = dynamic_cast<SimpleVector<double>&>(*sc);

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

   componentNorm = g->inf_norm();
   if (componentNorm > norm)
      norm = componentNorm;

   componentNorm = bA->inf_norm();
   if (componentNorm > norm)
      norm = componentNorm;

   componentNorm = A->inf_norm();
   if (componentNorm > norm)
      norm = componentNorm;

   componentNorm = C->inf_norm();
   if (componentNorm > norm)
      norm = componentNorm;

   assert(blx->matchesNonZeroPattern(*ixlow));
   componentNorm = blx->inf_norm();
   if (componentNorm > norm)
      norm = componentNorm;

   assert(bux->matchesNonZeroPattern(*ixupp));
   componentNorm = bux->inf_norm();
   if (componentNorm > norm)
      norm = componentNorm;

   assert(bl->matchesNonZeroPattern(*iclow));
   componentNorm = bl->inf_norm();
   if (componentNorm > norm)
      norm = componentNorm;

   assert(bu->matchesNonZeroPattern(*icupp));
   componentNorm = bu->inf_norm();
   if (componentNorm > norm)
      norm = componentNorm;

   return norm;
}

void Problem::print() {
   std::cout << "begin c\n";
   g->write_to_stream(std::cout);
   std::cout << "end c\n";

   std::cout << "begin xlow\n";
   blx->write_to_stream(std::cout);
   std::cout << "end xlow\n";
   std::cout << "begin ixlow\n";
   ixlow->write_to_stream(std::cout);
   std::cout << "end ixlow\n";

   std::cout << "begin xupp\n";
   bux->write_to_stream(std::cout);
   std::cout << "end xupp\n";
   std::cout << "begin ixupp\n";
   ixupp->write_to_stream(std::cout);
   std::cout << "end ixupp\n";
   std::cout << "begin A\n";

   A->write_to_stream(std::cout);
   std::cout << "end A\n";
   std::cout << "begin b\n";
   bA->write_to_stream(std::cout);
   std::cout << "end b\n";
   std::cout << "begin C\n";
   C->write_to_stream(std::cout);
   std::cout << "end C\n";

   std::cout << "begin clow\n";
   bl->write_to_stream(std::cout);
   std::cout << "end clow\n";
   std::cout << "begin iclow\n";
   iclow->write_to_stream(std::cout);
   std::cout << "end iclow\n";

   std::cout << "begin cupp\n";
   bu->write_to_stream(std::cout);
   std::cout << "end cupp\n";
   std::cout << "begin icupp\n";
   icupp->write_to_stream(std::cout);
   std::cout << "end icupp\n";
}