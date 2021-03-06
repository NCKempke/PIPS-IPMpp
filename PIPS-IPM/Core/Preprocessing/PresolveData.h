/*
 * PresolveData.h
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_
#define PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_

#include "DistributedProblem.hpp"
#include "StochPostsolver.h"
#include "SparseStorageDynamic.h"
#include "SystemType.h"

#include <algorithm>
#include <list>
#include <limits>
#include <queue>
#include <memory>

class PresolveData {
private:
   DistributedProblem* presProb;

   StochPostsolver* const postsolver;

   const double limit_max_bound_accepted;

   const int length_array_outdated_indicators{6};
   bool* array_outdated_indicators;
   bool& outdated_lhsrhs;
   bool& outdated_nnzs;
   bool& outdated_linking_var_bounds;
   bool& outdated_activities;
   bool& outdated_obj_vector;
   bool& postsolve_linking_row_propagation_needed;

   /* counter to indicate how many linking row bounds got changed locally and thus need activity recomputation */
   int linking_rows_need_act_computation{0};

   /* number of non-zero elements of each row / column */
   std::unique_ptr<DistributedVector<int>> nnzs_row_A;
   std::unique_ptr<DistributedVector<int>> nnzs_row_C;
   std::unique_ptr<DistributedVector<int>> nnzs_col;

   /* size of non-zero changes array = #linking rows A + #linking rows C + # linking variables */
   std::vector<int> array_nnz_chgs;
   std::unique_ptr<DenseVector<int>> nnzs_row_A_chgs{};
   std::unique_ptr<DenseVector<int>> nnzs_row_C_chgs{};
   std::unique_ptr<DenseVector<int>> nnzs_col_chgs{};

   /* In the constructor all unbounded entries will be counted.
    * Unbounded entries mean variables with non-zero multiplier that are unbounded in either upper or lower direction.
    * Activities will be computed once the amount of unbounded variables in upper or lower direction falls below 2 so
    * that bound strengthening becomes possible.
    */
   /* StochVecs for upper and lower activities and unbounded entries */
   std::unique_ptr<DistributedVector<double>> actmax_eq_part{};
   std::unique_ptr<DistributedVector<double>> actmin_eq_part{};

   std::unique_ptr<DistributedVector<int>> actmax_eq_ubndd{};
   std::unique_ptr<DistributedVector<int>> actmin_eq_ubndd{};

   std::unique_ptr<DistributedVector<double>> actmax_ineq_part{};
   std::unique_ptr<DistributedVector<double>> actmin_ineq_part{};

   std::unique_ptr<DistributedVector<int>> actmax_ineq_ubndd{};
   std::unique_ptr<DistributedVector<int>> actmin_ineq_ubndd{};

   /// changes in boundedness and activities of linking rows get stored and synchronized
   std::vector<double> array_act_chgs;
   std::unique_ptr<DenseVector<double>> actmax_eq_chgs{};
   std::unique_ptr<DenseVector<double>> actmin_eq_chgs{};
   std::unique_ptr<DenseVector<double>> actmax_ineq_chgs{};
   std::unique_ptr<DenseVector<double>> actmin_ineq_chgs{};

   std::vector<int> array_act_unbounded_chgs;
   std::unique_ptr<DenseVector<int>> actmax_eq_ubndd_chgs{};
   std::unique_ptr<DenseVector<int>> actmin_eq_ubndd_chgs{};
   std::unique_ptr<DenseVector<int>> actmax_ineq_ubndd_chgs{};
   std::unique_ptr<DenseVector<int>> actmin_ineq_ubndd_chgs{};

   /* handling changes in bounds */
   std::vector<double> array_bound_chgs;
   std::unique_ptr<DenseVector<double>> bound_chgs_A{};
   std::unique_ptr<DenseVector<double>> bound_chgs_C{};

   /* storing so far found singleton rows and columns */
   std::queue<INDEX> singleton_rows;
   std::queue<INDEX> singleton_cols;

   const int my_rank{PIPS_MPIgetRank()};
   const bool distributed{PIPS_MPIgetDistributed()};

   const double INF_NEG;
   const double INF_POS;

   // number of children
   const int nChildren;

   /* should we track a row/column through the presolving process - set in StochOptions */
   const bool track_row;
   const bool track_col;

   const INDEX tracked_row;
   const INDEX tracked_col;

   // objective offset created by presolving
   double objOffset{0.0};
   double obj_offset_chgs{0.0};
   std::unique_ptr<DenseVector<double>> objective_vec_chgs{};

   // store free variables which bounds are only implied by bound tightening to remove bounds later again
   std::unique_ptr<DistributedVector<int>> lower_bound_implied_by_system{};
   std::unique_ptr<DistributedVector<int>> lower_bound_implied_by_row{};
   std::unique_ptr<DistributedVector<int>> lower_bound_implied_by_node{};

   // TODO a vector of INDEX would be nicer
   std::unique_ptr<DistributedVector<int>> upper_bound_implied_by_system{};
   std::unique_ptr<DistributedVector<int>> upper_bound_implied_by_row{};
   std::unique_ptr<DistributedVector<int>> upper_bound_implied_by_node{};

   /* storing biggest and smallest absolute nonzero-coefficient in system matrix (including objective vector) */
   std::unique_ptr<DistributedVector<double>> absmin_col{};
   std::unique_ptr<DistributedVector<double>> absmax_col{};

   bool in_bound_tightening{false};
   std::vector<int> store_linking_row_boundTightening_A;
   std::vector<int> store_linking_row_boundTightening_C;

public :

   PresolveData(const DistributedProblem& sorigprob, StochPostsolver* postsolver);
   ~PresolveData();

   [[nodiscard]] const DistributedProblem& getPresProb() const { return *presProb; };

   [[nodiscard]] double getObjOffset() const { return objOffset; };
   [[nodiscard]] int getNChildren() const { return nChildren; };

   void getRowActivities(const INDEX& row, double& max_act, double& min_act, int& max_ubndd, int& min_ubndd) const;
   [[nodiscard]] std::pair<double,double> getRowBounds(const INDEX& row) const;
   [[nodiscard]] std::pair<double,double> getColBounds(const INDEX& col) const;


   [[nodiscard]] double getRowCoeff(const INDEX& row, const INDEX& col) const;

   [[nodiscard]] const DistributedVector<int>& getNnzsRow(SystemType system_type) const { return (system_type == EQUALITY_SYSTEM) ? *nnzs_row_A : *nnzs_row_C; }
   [[nodiscard]] const DistributedVector<int>& getNnzsRowA() const { return *nnzs_row_A; }; // todo maybe this is a problem - these counters might not be up to date
   [[nodiscard]] const DistributedVector<int>& getNnzsRowC() const { return *nnzs_row_C; };
   [[nodiscard]] const DistributedVector<int>& getNnzsCol() const { return *nnzs_col; };

   [[nodiscard]] int getNnzsRow(const INDEX& row) const;
   [[nodiscard]] int getNnzsCol(const INDEX& col) const;

   std::queue<INDEX>& getSingletonRows() { return singleton_rows; };
   std::queue<INDEX>& getSingletonCols() { return singleton_cols; };

   void delete_transposed();
   DistributedProblem* finalize();

   /* reset originally free variables' bounds to +- inf iff their current bounds are still implied by the problem */
   void resetOriginallyFreeVarsBounds(const DistributedProblem& orig_prob);

   /* whether or not there is currently changes buffered that need synchronization among all procs */
   bool reductionsEmpty();

   /* checks activities, non-zeros and root node */
   [[nodiscard]] bool presolve_dataInSync() const;

   /// synchronizing the problem over all mpi processes if necessary
   // TODO : add a allreduceEverything method that simply calls all the others
   void allreduceLinkingVarBounds();
   void allreduceAndApplyLinkingRowActivities();
   void allreduceAndApplyNnzChanges();
   void allreduceAndApplyBoundChanges();
   void allreduceAndApplyObjVecChanges();
   void allreduceObjOffset();

   [[nodiscard]] bool wasColumnRemoved(const INDEX& col) const;
   [[nodiscard]] bool wasRowRemoved(const INDEX& row) const;

   /// interface methods called from the presolvers when they detect a possible modification
   void startColumnFixation();
   void fixColumn(const INDEX& col, double value);
   void fixEmptyColumn(const INDEX& col, double val);

   void removeSingletonRow(const INDEX& row, const INDEX& col, double xlow_new, double xupp_new, double coeff);
   void removeSingletonRowSynced(const INDEX& row, const INDEX& col, double xlow_new, double xupp_new, double coeff);

   void syncPostsolveOfBoundsPropagatedByLinkingRows();

   void startBoundTightening();
   bool rowPropagatedBounds(const INDEX& row, const INDEX& col, double ubx, double lbx);
   void endBoundTightening();

   void startParallelRowPresolve();
   void
   substituteVariableNearlyParallelRows(const INDEX& row1, const INDEX& row2, const INDEX& col1, const INDEX& col2, double scalar, double translation,
         double parallelity);
   void tightenBoundsNearlyParallelRows(const INDEX& row1, const INDEX& row2, const INDEX& col1, const INDEX& col2, double xlow_new, double xupp_new,
         double scalar, double translation, double parallel_factor);

   void removeRedundantParallelRow(const INDEX& rm_row, const INDEX& par_row);
   void removeRedundantRow(const INDEX& row);
   void removeRedundantSide(const INDEX& row, bool is_upper_side);

   void startSingletonColumnPresolve();
   void fixColumnInequalitySingleton(const INDEX& col, const INDEX& row, double value, double coeff);
   void removeImpliedFreeColumnSingletonEqualityRow(const INDEX& row, const INDEX& col);
   void removeImpliedFreeColumnSingletonEqualityRowSynced(const INDEX& row, const INDEX& col);

   void removeFreeColumnSingletonInequalityRow(const INDEX& row, const INDEX& col, double coeff);
   void removeFreeColumnSingletonInequalityRowSynced(const INDEX& row, const INDEX& col, double coeff);

   void tightenRowBoundsParallelRow(const INDEX& row_tightened, const INDEX& row_tightening, double clow_new, double cupp_new, double factor);

   void transfrom_ineqalities_to_equalities();

   /* call whenever a single entry has been deleted from the matrix */
   void deleteEntryAtIndex(const INDEX& row, const INDEX& col, int col_index);

   /* methods for verifying state of presolve_data or querying the problem */
   [[nodiscard]] bool verifyNnzcounters() const;
   [[nodiscard]] bool verifyActivities() const;

   [[nodiscard]] bool nodeIsDummy(int node) const;
   [[nodiscard]] bool hasLinking(SystemType system_type) const;

   /* compute and update activities */
   void recomputeActivities() { recomputeActivities(false); }

   bool varBoundImpliedFreeBy(bool upper, const INDEX& col, const INDEX& row);
private:
   [[nodiscard]] bool iTrackColumn() const;
   [[nodiscard]] bool iTrackRow() const;

   void setRowBounds(const INDEX& row, double clow, double cupp);
   bool updateColBounds(const INDEX& col, double xlow, double xupp);

   void setRowUpperBound(const INDEX& row, double rhs) {
      row.inEqSys() ? setRowBounds(row, rhs, rhs) : setRowBounds(row, INF_NEG, rhs);
   }

   void setRowLowerBound(const INDEX& row, double lhs) {
      row.inEqSys() ? setRowBounds(row, lhs, lhs) : setRowBounds(row, lhs, INF_POS);
   }

   bool updateColLowerBound(const INDEX& col, double xlow) {
      return updateColBounds(col, xlow, INF_POS);
   }

   bool updateColUpperBound(const INDEX& col, double xupp) {
      return updateColBounds(col, INF_NEG, xupp);
   }

   void adaptObjectiveSubstitutedRow(const INDEX& row, const INDEX& col, double obj_coeff, double col_coeff);
   void addCoeffColToRow(double coeff, const INDEX& col, const INDEX& row);

   INDEX getRowMarkedAsImplyingColumnBound(const INDEX& col, bool upper_bound);
   void markRowAsImplyingColumnBound(const INDEX& col, const INDEX& row, bool upper_bound);

   void markColumnRemoved(const INDEX& col);

   void varboundImpliedFreeFullCheck(bool& upper_implied, bool& lower_implied, const INDEX& col, const INDEX& row) const;

   /// methods for printing debug information
   // initialize row and column nnz counter
   void initNnzCounter(DistributedVector<int>& nnzs_row_A, DistributedVector<int>& nnzs_row_C, DistributedVector<int>& nnzs_col) const;
   void initSingletons();

   void initAbsminAbsmaxInCols(DistributedVector<double>& absmin, DistributedVector<double>& absmax) const;

   void setUndefinedVarboundsTo(double value);
   void setUndefinedRowboundsTo(double value);

   static void addActivityOfBlock(const SparseStorageDynamic& matrix, DenseVector<double>& min_partact, DenseVector<int>& unbounded_min,
         DenseVector<double>& max_partact, DenseVector<int>& unbounded_max, const DenseVector<double>& xlow, const DenseVector<double>& ixlow,
         const DenseVector<double>& xupp, const DenseVector<double>& ixupp) ;

   long resetOriginallyFreeVarsBounds(const DenseVector<double>& ixlow_orig, const DenseVector<double>& ixupp_orig, int node);

   void adjustMatrixRhsLhsBy(const INDEX& row, double value, bool at_root);
   /// methods for modifying the problem
   void adjustRowActivityFromDeletion(const INDEX& row, const INDEX& col, double coeff);
   /// set bounds if new bound is better than old bound
   void updateRowActivities(const INDEX& col, double xlow_new, double xupp_new, double xlow_old, double xupp_old);

   void updateRowActivitiesBlock(const INDEX& row, const INDEX& col, double xlow_new, double xupp_new, double xlow_old, double xupp_old);

   void updateRowActivitiesBlock(const INDEX& row, const INDEX& col, double bound, double old_bound, bool upper);

   /* computes all row activities and number of unbounded variables per row
    * If there is more than one unbounded variable in the min/max activity of a row
    * +/-infinity() is stored. Else the actual partial activity is computed and stored.
    * For rows with one unbounded variable we store the partial activity without that
    * one variable, for rows with zero unbounded vars the stored activity is the actual
    * activity of that row.
    */
   void recomputeActivities(bool linking_only);

   void recomputeActivities(bool linkinig_only, DistributedVector<double>& actmax_eq_part, DistributedVector<double>& actmin_eq_part,
         DistributedVector<int>& actmax_eq_ubndd, DistributedVector<int>& actmin_eq_ubndd, DistributedVector<double>& actmax_ineq_part,
         DistributedVector<double>& actmin_ineq_part, DistributedVector<int>& actmax_ineq_ubndd, DistributedVector<int>& actmin_ineq_ubndd) const;

   [[nodiscard]] double computeLocalLinkingRowMinOrMaxActivity(const INDEX& row, bool upper) const;
   void computeRowMinOrMaxActivity(const INDEX& row, bool upper);

   void removeColumn(const INDEX& col, double fixation);
   void removeColumnFromMatrix(const INDEX& row, const INDEX& col, double fixation);
   void removeRow(const INDEX& row);
   void removeRowFromMatrix(const INDEX& row, const INDEX& col);

   void reduceNnzCounterRowBy(const INDEX& row, int amount, bool at_root);
   void increaseNnzCounterRowBy(const INDEX& row, int amount, bool at_root);

   void changeNnzCounterRow(const INDEX& row, int amount, bool at_root);

   void reduceNnzCounterColumnBy(const INDEX& col, int amount, bool at_root);
   void increaseNnzCounterColumnBy(const INDEX& col, int amount, bool at_root);

   void changeNnzCounterColumn(const INDEX& col, int amount, bool at_root);

   /// methods for querying the problem in order to get certain structures etc.
   [[nodiscard]] DistributedMatrix& getSystemMatrix(SystemType system_type) const;
   [[nodiscard]] SparseMatrix* getSparseGenMatrix(const INDEX& row, const INDEX& col) const;

   void checkBoundsInfeasible(const INDEX& col, double xlow_new, double xupp_new) const;

   void transform_inequalities_into_equalities(int node);
   void transform_inequalities_into_equalities(int node, bool linking);

   void append_bounds_inequalities_to_equalities_transformation(int node, bool linking, int n_slack_variables);
   void append_new_slacks_to_objective_vector(int node, int n_new_slack_variables);
   void extend_q_matrix_by_new_variables(int node, int n_variables);
   void adjust_nonzeros_after_inequalities_to_equalities_transformation(int node, bool linking, const DenseVector<int>& nonzero_pattern_slacks);
   void transform_matrices_inequalities_into_equalities(int node, int n_new_slack_variables, const std::vector<int>& diagonal_for_identity);
   void transform_linking_matrices_inequalities_into_equalities(int n_new_slack_variables, const std::vector<int>& diagonal_for_identity);
   void extend_linking_variable_child_matrices_by(int n_new_slack_variables);
   std::vector<int> get_slack_diagonal_for_inequality_equality_tranformation(int node, bool linking);

public:
   void writeRowLocalToStreamDense(std::ostream& out, const INDEX& row) const;
   void printRowColStats() const;
   [[nodiscard]] int countEmptyRowsBDmat() const;

private:
   void writeMatrixRowToStreamDense(std::ostream& out, const SparseMatrix& mat, int node, int row, const DenseVector<double>& ixupp,
         const DenseVector<double>& xupp, const DenseVector<double>& ixlow, const DenseVector<double>& xlow) const;
   void printVarBoundStatistics(std::ostream& out) const;
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_ */
