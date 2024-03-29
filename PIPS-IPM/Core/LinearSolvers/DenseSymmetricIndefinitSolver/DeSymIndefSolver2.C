/* PIPS                                                               *
 * Authors: Miles Lubin                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "DeSymIndefSolver2.h"
#include "DenseVector.hpp"
#include <cassert>

#include "OoqpBlas.h"
#include "DenseSymmetricMatrix.h"

DeSymIndefSolver2::DeSymIndefSolver2(const DenseSymmetricMatrix& dm, int nx) : mStorage(std::make_unique<DenseStorage>(dm.getStorage())), nx(nx) {
   n = mStorage->n;
   ny = n - nx;
}
//#include "mpi.h"
void DeSymIndefSolver2::matrixChanged() {

   char fortranUplo = 'U';
   char fortranL = 'L';
   char fortranN = 'N';
   char fortranT = 'T';
   int info;
   double* mat = &mStorage->M[0][0];
   double one = 1, zero = 0;

   // treat matrix as column-major and assume upper
   // triangle is filled
   // this is different from the elemental version
   // which uses the lower triangle

   /*
   start with:
   [Q A^T
    *  0]

   cholesky on Q to get M^T
   trsm to get M^-1A^T
   syrk to form (M^-1A^T)^T(M^-1A^T)
   cholesky on bottom right
   */


   dpotrf_(&fortranUplo, &nx, mat, &n, &info);
   if (info != 0) {
      std::cerr << "error factoring Q block: info = " << info << "\n";
   }
   assert(info == 0);
   if (ny == 0)
      return;

   printf("dtrsm_\n");

   dtrsm_(&fortranL, &fortranUplo, &fortranT, &fortranN, &nx, &ny, &one, mat, &n, mat + nx * n, &n);

   printf("dsyrk_\n");

   dsyrk_(&fortranUplo, &fortranT, &ny, &nx, &one, mat + nx * n, &n, &zero, mat + nx * n + nx, &n);

   printf("dpotrf2\n");

   dpotrf_(&fortranUplo, &ny, mat + nx * n + nx, &n, &info);
   if (info != 0) {
      std::cerr << "error factoring AQ^-1A^T block: info = " << info << "\n";
   }
   assert(info == 0);

   printf("finished factorization\n");


}

void DeSymIndefSolver2::solve(Vector<double>& v) {
   char fortranUplo = 'U';
   char fortranT = 'T';
   char fortranN = 'N';

   double* mat = &mStorage->M[0][0];
   int one = 1;
   double minus1 = -1;

   DenseVector<double>& sv = dynamic_cast<DenseVector<double>&>(v);
   double* rhs = &sv[0];

   dtrsv_(&fortranUplo, &fortranT, &fortranN, &n, mat, &n, rhs, &one);

   if (ny > 0)
      dscal_(&ny, &minus1, rhs + nx, &one);

   dtrsv_(&fortranUplo, &fortranN, &fortranN, &n, mat, &n, rhs, &one);

   printf("finished rhs solve\n");

}

void DeSymIndefSolver2::diagonalChanged(int /* idiag */, int /* extent */) {
   this->matrixChanged();
}

DeSymIndefSolver2::~DeSymIndefSolver2() {
}
