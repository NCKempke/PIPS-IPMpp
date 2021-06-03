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
   void createChildren();

   void print_solver_regularization_and_sc_info(const std::string&& name) const;

private:
   void init();

public:
   std::vector<DistributedLinearSystem*> children;

   DistributedRootLinearSystem(DistributedFactory* factory_, DistributedQP* prob_, bool is_hierarchy_root = false);

   DistributedRootLinearSystem(DistributedFactory* factory, DistributedQP* prob_, std::shared_ptr<Vector<double>> dd_,
      std::shared_ptr<Vector<double>> dq_, std::shared_ptr<Vector<double>> nomegaInv_,
      std::shared_ptr<Vector<double>> primal_reg_, std::shared_ptr<Vector<double>> dual_y_reg_,
      std::shared_ptr<Vector<double>> dual_z_reg_, std::shared_ptr<Vector<double>> rhs_, bool create_sub_root_solver);

   void factor2() override;

   void assembleKKT() override;

   void allreduceAndFactorKKT() override;

   /* Atoms methods of FACTOR2 for a non-leaf linear system */
   virtual void initializeKKT();

   virtual void assembleLocalKKT() = 0;

   void addTermToSchurCompl(size_t childindex, bool use_local_RAC);

   virtual void reduceKKT();

   virtual void factorizeKKT();

   virtual void finalizeKKT() = 0;

   virtual void finalizeKKTdist() { assert("not implemented here \n" && 0); };

   void Ltsolve2(DistributedVector<double>& x, SimpleVector<double>& xp, bool) override;

   /* compute (Br0 - sum_j Br_mod_border) - buffer */
   virtual void
   finalizeZ0Hierarchical(DenseMatrix& buffer, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, int begin_rows,
      int end_rows);

   /* compute SC += B0_{outer}^T X0 */
   virtual void
   finalizeInnerSchurComplementContribution(AbstractMatrix& SC, DenseMatrix& X0, BorderLinsys& Br, bool is_sym,
      bool is_sparse, int begin_rows,
      int end_rows);

   /* compute -SUM_i Bi_{inner} Ki^-1 Bi_{outer} */
   using DistributedLinearSystem::LsolveHierarchyBorder;

   void
   LsolveHierarchyBorder(DenseMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
      bool use_local_RAC, bool two_link_border,
      int begin_cols, int end_cols) override;

   /* compute SUM_i Bli^T X_i = Bli^T Ki^-1 ( ( Bri - sum_j Bmodij Xij ) - Bi_{inner} X0) */
   using DistributedLinearSystem::LtsolveHierarchyBorder;

   void LtsolveHierarchyBorder(AbstractMatrix& res, const DenseMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br,
      std::vector<BorderMod>& br_mod_border,
      bool sym_res, bool sparse_res, bool use_local_RAC, int begin_cols, int end_cols) override;

   void put_primal_diagonal() override;

   void put_dual_inequalites_diagonal() override;

   void clear_dual_equality_diagonal() override;

   void put_barrier_parameter(double barrier) override;

   virtual void AddChild(DistributedLinearSystem* child);

   virtual bool usingSparseKkt() { return hasSparseKkt; };

   ~DistributedRootLinearSystem() override;

   //utilities
   static void myAtPutZeros(DenseSymmetricMatrix* mat);

   static void myAtPutZeros(DenseSymmetricMatrix* mat, int row, int col, int rowExtent, int colExtent);

   // all_reduces specified submatrix (in chunks)
   static void
   submatrixAllReduce(DenseSymmetricMatrix& A, int startRow, int startCol, int nRows, int nCols, MPI_Comm comm);

   void allreduceMatrix(AbstractMatrix& mat, bool is_sparse, bool is_sym, MPI_Comm comm);

   static void
   submatrixAllReduceFull(DenseSymmetricMatrix& A, int startRow, int startCol, int nRows, int nCols, MPI_Comm comm);

   static void submatrixAllReduceFull(DenseMatrix& A, int startRow, int startCol, int nRows, int nCols, MPI_Comm comm);

   // all_reduces specified submatrix as a while
   static void submatrixAllReduceFull(double** A, int startRow, int startCol, int nRows, int nCols, MPI_Comm comm);

   // all_reducees lower half (including diagonal) of specified submatrix
   static void submatrixAllReduceDiagLower(DenseSymmetricMatrix& A, int substart, int subsize, MPI_Comm comm);

   SCsparsifier precondSC;

protected: //buffers

   void createSolverAndSchurComplement(bool sub_root);
   void createSparseSolver(bool sub_root);
   void createDenseSolver();
   [[nodiscard]] virtual std::unique_ptr<SymmetricMatrix> createKKT() const;

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

   void reduceKKTdist();

   void reduceKKTdense();

   void reduceKKTsparse();

   void reduceToProc0(int size, double* values);

   void reduceToAllProcs(int size, double* values);

   void syncKKTdistLocalEntries();

   void sendKKTdistLocalEntries(const std::vector<MatrixEntryTriplet>& prevEntries) const;

   [[nodiscard]] std::vector<MatrixEntryTriplet> receiveKKTdistLocalEntries() const;

   [[nodiscard]] std::vector<MatrixEntryTriplet> packKKTdistOutOfRangeEntries(int childStart, int childEnd) const;

   static void finalizeInnerSchurComplementContributionDense(AbstractMatrix& SC_, const DenseMatrix& X0,
      const SparseMatrix* A0_border, const SparseMatrix* C0_border, const SparseMatrix* A00_border,
      const SparseMatrix* F0vec_border, const SparseMatrix* G0vec_border, const SparseMatrix* F0cons_border,
      const SparseMatrix* G0cons_border, bool is_sym,
      int begin_rows, int end_rows);

   static void finalizeInnerSchurComplementContributionSparse(AbstractMatrix& SC_, const DenseMatrix& X0,
      const SparseMatrix* A0_border, const SparseMatrix* C0_border, const SparseMatrix* A00_border,
      const SparseMatrix* F0vec_border, const SparseMatrix* G0vec_border, const SparseMatrix* F0cons_border,
      const SparseMatrix* G0cons_border, int begin_rows,
      int end_rows);

   MPI_Datatype MatrixEntryTriplet_mpi;

#ifdef STOCH_TESTING
   protected:
    static void dumpRhs(int proc, const char* nameToken,  SimpleVector<double>& rhs);
    static void dumpMatrix(int scen, int proc, const char* nameToken, DenseSymmetricMatrix& M);
#endif
#ifdef TIMING
   protected:
    void afterFactor();
#endif
};

#endif

