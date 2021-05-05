/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#ifndef DISTRIBUTEDROOTLINEARSYSTEM_H
#define DISTRIBUTEDROOTLINEARSYSTEM_H

#include "DistributedLinearSystem.h"
#include "DistributedMatrix.h"
#include "SCsparsifier.h"

class SCsparsifier;

class DistributedFactory;

class DistributedQP;

/** 
 * ROOT (= NON-leaf) linear system
 */
class DistributedRootLinearSystem : public DistributedLinearSystem {
   struct MatrixEntryTriplet {
      double val;
      int row;
      int col;
   };

protected:
   void createChildren(DistributedQP* prob);
   void deleteChildren() override;

private:
   void init();

public:
   std::vector<DistributedLinearSystem*> children;

   DistributedRootLinearSystem(DistributedFactory* factory_, DistributedQP* prob_, bool is_hierarchy_root = false);
   DistributedRootLinearSystem(DistributedFactory* factory, DistributedQP* prob_, Vector<double>* dd_, Vector<double>* dq_, Vector<double>* nomegaInv_,
         Vector<double>* primal_reg_, Vector<double>* dual_y_reg_, Vector<double>* dual_z_reg_, Vector<double>* rhs_);

   void factor2(DistributedQP* prob, Variables* vars) override;
   void assembleKKT(DistributedQP* prob, Variables* vars) override;
   void allreduceAndFactorKKT(DistributedQP* prob, Variables* vars) override;

   /* Atoms methods of FACTOR2 for a non-leaf linear system */
   virtual void initializeKKT(DistributedQP* prob, Variables* vars);
   virtual void assembleLocalKKT(DistributedQP* prob) = 0;
   void addTermToSchurCompl(DistributedQP* prob, size_t childindex, bool use_local_RAC);
   virtual void reduceKKT(DistributedQP* prob);
   virtual void factorizeKKT();
   virtual void factorizeKKT(DistributedQP* prob);
   virtual void finalizeKKT(DistributedQP* prob, Variables* vars) = 0;
   virtual void finalizeKKTdist(DistributedQP* /*prob*/ ) { assert("not implemented here \n" && 0); };

   void Ltsolve2(DistributedQP* prob, DistributedVector<double>& x, SimpleVector<double>& xp, bool) override;

   /* compute (Br0 - sum_j Br_mod_border) - buffer */
   virtual void finalizeZ0Hierarchical(DenseMatrix& buffer, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, int begin_rows, int end_rows);
   /* compute SC += B0_{outer}^T X0 */
   virtual void
   finalizeInnerSchurComplementContribution(AbstractMatrix& SC, DenseMatrix& X0, BorderLinsys& Br, bool is_sym, bool is_sparse, int begin_rows,
         int end_rows);

   /* compute -SUM_i Bi_{inner} Ki^-1 Bi_{outer} */
   using DistributedLinearSystem::LsolveHierarchyBorder;
   void
   LsolveHierarchyBorder(DenseMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool use_local_RAC, bool two_link_border,
         int begin_cols, int end_cols) override;

   /* compute SUM_i Bli^T X_i = Bli^T Ki^-1 ( ( Bri - sum_j Bmodij Xij ) - Bi_{inner} X0) */
   using DistributedLinearSystem::LtsolveHierarchyBorder;
   void LtsolveHierarchyBorder(AbstractMatrix& res, const DenseMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& br_mod_border,
         bool sym_res, bool sparse_res, bool use_local_RAC, int begin_cols, int end_cols) override;

   void addBorderTimesRhsToB0(DistributedVector<double>& rhs, SimpleVector<double>& b0, BorderLinsys& border) override;

   void addBorderX0ToRhs(DistributedVector<double>& rhs, const SimpleVector<double>& x0, BorderLinsys& border) override;

   void put_primal_diagonal() override;
   void put_dual_inequalites_diagonal() override;
   void clear_dual_equality_diagonal() override;

   void put_barrier_parameter(double barrier) override;

   virtual void AddChild(DistributedLinearSystem* child);

   virtual bool usingSparseKkt() { return hasSparseKkt; };

   ~DistributedRootLinearSystem() override;

   //utilities
   static void myAtPutZeros(DenseSymMatrix* mat);
   static void myAtPutZeros(DenseSymMatrix* mat, int row, int col, int rowExtent, int colExtent);

   // all_reduces specified submatrix (in chunks)
   static void submatrixAllReduce(DenseSymMatrix* A, int startRow, int startCol, int nRows, int nCols, MPI_Comm comm);

   void allreduceMatrix(AbstractMatrix& mat, bool is_sparse, bool is_sym, MPI_Comm comm);

   static void submatrixAllReduceFull(DenseSymMatrix* A, int startRow, int startCol, int nRows, int nCols, MPI_Comm comm);
   static void submatrixAllReduceFull(DenseMatrix* A, int startRow, int startCol, int nRows, int nCols, MPI_Comm comm);

   // all_reduces specified submatrix as a while
   static void submatrixAllReduceFull(double** A, int startRow, int startCol, int nRows, int nCols, MPI_Comm comm);

   // all_reducees lower half (including diagonal) of specified submatrix
   static void submatrixAllReduceDiagLower(DenseSymMatrix* A, int substart, int subsize, MPI_Comm comm);

   SCsparsifier precondSC;

protected: //buffers

   SparseSymmetricMatrix* kktDist{};

   Vector<double>* xDiag{};

   Vector<double>* zDiag{};
   Vector<double>* zDiagLinkCons{};

   double* sparseKktBuffer{};

   std::unique_ptr<DistributedVector<double>> sol_inner{};

   int childrenProperStart; // first non-dummy child
   int childrenProperEnd;   // end of non-dummy children range (not included)
   bool hasSparseKkt;
   bool usePrecondDist;
   bool allreduce_kkt;

private:
   void initProperChildrenRange();
   void registerMatrixEntryTripletMPI();
   void reduceKKTdist(DistributedQP* prob);
   void reduceKKTdense();
   void reduceKKTsparse();
   void reduceToProc0(int size, double* values);
   void reduceToAllProcs(int size, double* values);
   void syncKKTdistLocalEntries(DistributedQP* prob);
   void sendKKTdistLocalEntries(const std::vector<MatrixEntryTriplet>& prevEntries) const;
   std::vector<MatrixEntryTriplet> receiveKKTdistLocalEntries() const;
   std::vector<MatrixEntryTriplet> packKKTdistOutOfRangeEntries(DistributedQP* prob, int childStart, int childEnd) const;

   static void finalizeInnerSchurComplementContributionDense(AbstractMatrix& SC_, DenseMatrix& X0, SparseMatrix* A0_border, SparseMatrix* C0_border,
         SparseMatrix* F0vec_border, SparseMatrix* G0vec_border, SparseMatrix* F0cons_border, SparseMatrix* G0cons_border, bool is_sym,
         int begin_rows, int end_rows);

   static void finalizeInnerSchurComplementContributionSparse(AbstractMatrix& SC_, DenseMatrix& X0, SparseMatrix* A0_border, SparseMatrix* C0_border,
         SparseMatrix* F0vec_border, SparseMatrix* G0vec_border, SparseMatrix* F0cons_border, SparseMatrix* G0cons_border, int begin_rows,
         int end_rows);

   MPI_Datatype MatrixEntryTriplet_mpi;

#ifdef STOCH_TESTING
   protected:
    static void dumpRhs(int proc, const char* nameToken,  SimpleVector<double>& rhs);
    static void dumpMatrix(int scen, int proc, const char* nameToken, DenseSymMatrix& M);
#endif
#ifdef TIMING
   protected:
    void afterFactor();
#endif
};

#endif

