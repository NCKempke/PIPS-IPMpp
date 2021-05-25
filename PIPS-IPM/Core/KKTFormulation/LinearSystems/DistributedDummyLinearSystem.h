#ifndef DISTRIBUTEDDUMMYLINEARSYSTEM_H
#define DISTRIBUTEDDUMMYLINEARSYSTEM_H

#include "DistributedLinearSystem.h"

/** 
 * DUMMY Linear system class
 */
class DistributedDummyLinearSystem : public DistributedLinearSystem {
public:
   DistributedDummyLinearSystem(DistributedFactory* factory, DistributedQP* prob) : DistributedLinearSystem(factory, prob, nullptr, nullptr, nullptr, nullptr,
         nullptr, nullptr, nullptr, false) {
      mpiComm = MPI_COMM_NULL;
   };

   ~DistributedDummyLinearSystem() override = default;

   void factor2(DistributedQP*, Variables*) override {};
   void allreduceAndFactorKKT(DistributedQP*, Variables*) override {};
   void assembleKKT(DistributedQP*, Variables*) override {}

   void Lsolve(DistributedQP*, Vector<double>&) override {};
   void Dsolve(DistributedQP*, Vector<double>&) override {};
   void Ltsolve(DistributedQP*, Vector<double>&) override {};
   void Ltsolve2(DistributedQP*, DistributedVector<double>&, SimpleVector<double>&, bool) override {};

   void solveCompressed(Vector<double>&) override {};

   void put_primal_diagonal() override {};
   void put_dual_inequalites_diagonal() override {};
   void put_barrier_parameter(double) override {};
   void clear_dual_equality_diagonal() override {};

   void add_regularization_local_kkt(double, double, double) override {};

   void joinRHS(Vector<double>&, const Vector<double>&, const Vector<double>&, const Vector<double>&) const override {};

   void separateVars(Vector<double>&, Vector<double>&, Vector<double>&, const Vector<double>&) const override {};

   void addLnizi(DistributedQP*, Vector<double>&, Vector<double>&) override {};
   void addLniziLinkCons(DistributedQP*, Vector<double>&, Vector<double>&, bool) override {};

   /** y += alpha * Lni^T * x */
   //  void LniTransMult(DistributedQP *prob, SimpleVector<double>& y, double alpha, SimpleVector<double>& x) override {};

   void addTermToSchurResidual(DistributedQP*, SimpleVector<double>&, SimpleVector<double>&) override {};

   void LsolveHierarchyBorder(DenseMatrix&, BorderLinsys&, std::vector<BorderMod>&, bool, int, int) override {};
   void LsolveHierarchyBorder(DenseMatrix&, BorderLinsys&, std::vector<BorderMod>&, bool, bool, int, int) override {};
   void addInnerBorderKiInvBrToRes(AbstractMatrix&, BorderLinsys&, std::vector<BorderMod>&, bool, bool, bool, int, int, int) override {};
   void
   LniTransMultHierarchyBorder(AbstractMatrix&, const DenseMatrix&, BorderLinsys&, BorderLinsys&, std::vector<BorderMod>&, bool, bool, bool, int,
         int, int) override {};

   void deleteChildren() override {};

   [[nodiscard]] bool isDummy() const override { return true; };

   void addBorderTimesRhsToB0(DistributedVector<double>&, SimpleVector<double>&, BorderLinsys&) override {};
   void addBorderX0ToRhs(DistributedVector<double>&, const SimpleVector<double>&, BorderLinsys&) override {};
   void computeInnerSystemRightHandSide(DistributedVector<double>&, const SimpleVector<double>&, bool) override {};
};

#endif
