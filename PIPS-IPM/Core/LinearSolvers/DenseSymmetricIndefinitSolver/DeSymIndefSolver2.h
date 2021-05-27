/* PIPS                                                               *
 * Authors: Miles Lubin                                               *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DESYMINDEFSOLVER2_H
#define DESYMINDEFSOLVER2_H

#include "DoubleLinearSolver.h"
#include "DenseSymmetricMatrix.h"
#include <memory>

/** Specialized LDL^T solver for saddle point systems
 * @ingroup DenseLinearAlgebra
 * @ingroup LinearSolvers
 */
class DeSymIndefSolver2 : public DoubleLinearSolver {
protected:
   std::unique_ptr<DenseStorage> mStorage;
   int nx, ny, n;
public:
   DeSymIndefSolver2(const DenseSymmetricMatrix& storage, int nx);
   void diagonalChanged(int idiag, int extent) override;
   void matrixChanged() override;
   using DoubleLinearSolver::solve;
   void solve(Vector<double>& vec) override;
   ~DeSymIndefSolver2() override;


   [[nodiscard]] bool reports_inertia() const override { return true; };
   [[nodiscard]] std::tuple<unsigned int, unsigned int, unsigned int> get_inertia() const override {
      assert(false && "TODO : implement");
      return {0, 0, 0};
   };

};

#endif
