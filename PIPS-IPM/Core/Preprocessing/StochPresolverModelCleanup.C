/*
 * StochPresolverModelCleanup.C
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#include "StochPresolverModelCleanup.h"

#include "PIPSIPMppOptions.h"
#include <cmath>
#include <vector>

StochPresolverModelCleanup::StochPresolverModelCleanup(PresolveData& presolve_data, const DistributedProblem& origProb) : StochPresolverBase(presolve_data,
      origProb), limit_min_mat_entry(pipsipmpp_options::get_double_parameter("PRESOLVE_MODEL_CLEANUP_MIN_MATRIX_ENTRY")),
      limit_max_matrix_entry_impact(pipsipmpp_options::get_double_parameter("PRESOLVE_MODEL_CLEANUP_MAX_MATRIX_ENTRY_IMPACT")),
      limit_matrix_entry_impact_feasdist(pipsipmpp_options::get_double_parameter("PRESOLVE_MODEL_CLEANUP_MATRIX_ENTRY_IMPACT_FEASDIST")),
      removed_entries_total(0), fixed_empty_cols_total(0), removed_rows_total(0) {
}

bool StochPresolverModelCleanup::applyPresolving() {
   assert(presolve_data.reductionsEmpty());
   assert(presolve_data.presolve_dataInSync());

#ifndef NDEBUG
   if (my_rank == 0 && verbosity > 1) {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
      std::cout << "--- Before model cleanup:\n";
   }
   countRowsCols();
#endif

   std::vector<int> counts(3, 0);
   int& n_removed_entries = counts[0];
   int& n_fixed_empty_columns = counts[1];
   int& n_removed_rows = counts[2];

   // removal of redundant constraints
   int n_removed_rows_eq = removeRedundantRows(EQUALITY_SYSTEM);
   int n_removed_rows_ineq = removeRedundantRows(INEQUALITY_SYSTEM);
   n_removed_rows = n_removed_rows_eq + n_removed_rows_ineq;

   /* remove entries from A and C matrices and updates transposed systems */
   int n_removed_entries_eq = removeTinyEntriesFromSystem(EQUALITY_SYSTEM);
   int n_removed_entries_ineq = removeTinyEntriesFromSystem(INEQUALITY_SYSTEM);
   n_removed_entries = n_removed_entries_eq + n_removed_entries_ineq;

   int local_count_empty_brows = presolve_data.countEmptyRowsBDmat();
   PIPS_MPIgetSumInPlace(local_count_empty_brows);
   if (my_rank == 0)
      std::cout << "empty brows " << local_count_empty_brows << std::endl;

   presolve_data.allreduceAndApplyNnzChanges();

   n_fixed_empty_columns = fixEmptyColumns();

   // update all nnzCounters - set reductionStochvecs to zero afterwards
   presolve_data.allreduceAndApplyBoundChanges();
   presolve_data.allreduceObjOffset();
   presolve_data.allreduceAndApplyLinkingRowActivities();

   if (distributed)
      PIPS_MPIsumArrayInPlace(counts, MPI_COMM_WORLD);

   fixed_empty_cols_total += n_fixed_empty_columns;
   removed_entries_total += n_removed_entries;
   removed_rows_total += n_removed_rows;

#ifndef NDEBUG
   if (my_rank == 0 && verbosity > 1) {
      std::cout << "\tRemoved redundant rows in model cleanup: " << removed_rows_total << "\n";
      std::cout << "\tRemoved tiny entries in model cleanup: " << removed_entries_total << "\n";
      std::cout << "\tFixed empty columns in model cleanup: " << fixed_empty_cols_total << "\n";
      std::cout << "--- After model cleanup:" << "\n";
   }
   else if (my_rank == 0 && verbosity == 1)
      std::cout << "Clean:\t removed " << removed_rows_total << " rows, " << fixed_empty_cols_total << " cols, " << removed_entries_total << " entries\n";

   countRowsCols();
   if (my_rank == 0 && verbosity > 1)
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
#endif

   assert(presolve_data.reductionsEmpty());
   assert(presolve_data.presolve_dataInSync());

   if (n_removed_entries != 0 || n_removed_rows != 0 || n_fixed_empty_columns != 0)
      return true;
   else
      return false;
}

/** Remove redundant rows in the constraint system. Compares the minimal and maximal row activity
 * with the row bounds (lhs and rhs). If a row is found to be redundant, it is removed.
 * If infeasiblity is detected, then Abort.
 */
int StochPresolverModelCleanup::removeRedundantRows(SystemType system_type) {
   int nRemovedRows = 0;

   // root:
   int nRemovedRowsRoot = removeRedundantRows(system_type, -1);

   if (my_rank == 0)
      nRemovedRows += nRemovedRowsRoot;

   // children:
   for (int node = 0; node < nChildren; node++)
      if (!presolve_data.nodeIsDummy(node))
         nRemovedRows += removeRedundantRows(system_type, node);

   return nRemovedRows;
}

int StochPresolverModelCleanup::removeRedundantRows(SystemType system_type, int node) {
   assert(!presolve_data.nodeIsDummy(node));
   int n_removed_rows = 0;
   int n_removed_rows_link = 0;
   if (node == -1)
      n_removed_rows_link = removeRedundantRows(system_type, node, true);

   n_removed_rows += removeRedundantRows(system_type, node, false);

   return n_removed_rows + n_removed_rows_link;
}

int StochPresolverModelCleanup::removeRedundantRows(SystemType system_type, int node, bool linking) {
   assert(-1 <= node && node < nChildren);
   assert((linking && node == -1) || !linking);
   assert(!presolve_data.nodeIsDummy(node));

   if (linking && !presolve_data.hasLinking(system_type))
      return 0;
   updatePointersForCurrentNode(node, system_type);

   int n_removed_rows = 0;

   const DenseVector<int>& nnzs = !linking ? *currNnzRow : *currNnzRowLink;

   for (int row_index = 0; row_index < nnzs.length(); ++row_index) {

      const INDEX row(ROW, node, row_index, linking, system_type);

      if (presolve_data.wasRowRemoved(row)) {
         assert(nnzs[row_index] == 0);
         continue;
      }

      double actmin_part, actmax_part;
      int actmin_ubndd, actmax_ubndd;
      presolve_data.getRowActivities(row, actmax_part, actmin_part, actmax_ubndd, actmin_ubndd);
      const auto [clow, cupp] = presolve_data.getRowBounds(row);

      const bool has_clow = clow != INF_NEG;
      const bool has_cupp = cupp != INF_POS;
      const bool has_max_activitiy = actmax_ubndd == 0;
      const bool has_min_activity = actmin_ubndd == 0;

      if (row.inEqSys()) {
         assert(clow == cupp);
         if (actmin_ubndd != 0 || actmax_ubndd != 0)
            continue;

         if ((PIPSisLTFeas(clow, actmin_part) && has_min_activity) ||
             (PIPSisLTFeas(actmax_part, cupp) && has_max_activitiy)) {
            PIPS_MPIabortInfeasible("Found row that cannot meet it's rhs with it's computed activities", "StochPresolverModelCleanup.C",
                  "removeRedundantRows");
         }
         else if (PIPSisLEFeas(clow, actmin_part) && PIPSisLEFeas(actmax_part, cupp)) {
            presolve_data.removeRedundantRow(row);
            n_removed_rows++;
         }
      }
      else {
         assert(row.inInEqSys());

         if ((has_clow && has_max_activitiy && PIPSisLTFeas(actmax_part, clow)) ||
             (has_cupp && has_min_activity && PIPSisLTFeas(cupp, actmin_part)))
            PIPS_MPIabortInfeasible("Found row that cannot meet it's lhs or rhs with it's computed activities", "StochPresolverModelCleanup.C",
                  "removeRedundantRows");
         else if ((!has_clow || PIPSisLE(clow, -infinity)) &&
                  (!has_cupp || PIPSisLE(infinity, cupp))) {
            presolve_data.removeRedundantRow(row);
            n_removed_rows++;
         }
         else if ((!has_clow || (has_min_activity && PIPSisLEFeas(clow, actmin_part))) &&
                  (!has_cupp || (has_max_activitiy && PIPSisLEFeas(actmax_part, cupp)))) {
            presolve_data.removeRedundantRow(row);
            n_removed_rows++;
         }
         else if (has_cupp && has_max_activitiy && PIPSisLEFeas(actmax_part, cupp)) {
            presolve_data.removeRedundantSide(row, true);
         }
         else if (has_clow && has_min_activity && PIPSisLEFeas(clow, actmin_part)) {
            presolve_data.removeRedundantSide(row, false);
         }
      }
   }
   return n_removed_rows;
}

/* removes all small entries from the specified system
 *
 * While removing the reduction vectors in presolve_data get set - nnzVectors are not updated an that must be
 * done using updateNnzFromReductions if needed.
 * Transposed matrices get updated in a subroutine - so after calling this method, that matrix should
 * be in a consistent state.
 */
int StochPresolverModelCleanup::removeTinyEntriesFromSystem(SystemType system_type) {
   assert(dynamic_cast<const DistributedMatrix&>(*(presolve_data.getPresProb().equality_jacobian)).children.size() == (size_t) nChildren);
   assert(dynamic_cast<const DistributedMatrix&>(*(presolve_data.getPresProb().inequality_jacobian)).children.size() == (size_t) nChildren);

   int n_elims = 0;

   assert(!presolve_data.nodeIsDummy(-1));

   /* reductions in root node */
   /* process B0 and Bl0 */
   n_elims += removeTinyInnerLoop(system_type, -1, B_MAT);
   if (presolve_data.hasLinking(system_type))
      n_elims += removeTinyInnerLoop(system_type, -1, BL_MAT);

   /* count eliminations in B0 and Bl0 only once */
   if (distributed && my_rank != 0)
      n_elims = 0;

   // go through the children
   for (int node = 0; node < nChildren; node++) {
      if (!presolve_data.nodeIsDummy(node)) {
         /* Amat */
         n_elims += removeTinyInnerLoop(system_type, node, A_MAT);

         /* Bmat */
         n_elims += removeTinyInnerLoop(system_type, node, B_MAT);

         /* this has to be synchronized */
         /* Blmat */
         if (presolve_data.hasLinking(system_type))
            n_elims += removeTinyInnerLoop(system_type, node, BL_MAT);
      }
   }

   return n_elims;
}

/** Removes tiny entries in storage and adapts the lhs/rhs accordingly.
 * system type indicates matrix A or C, block_type indicates the block
 */
int StochPresolverModelCleanup::removeTinyInnerLoop(SystemType system_type, int node, BlockType block_type) {
   if (presolve_data.nodeIsDummy(node))
      return 0;

   updatePointersForCurrentNode(node, system_type);

   const SparseStorageDynamic* mat = nullptr;

   const DenseVector<double>* x_lower = nullptr;
   const DenseVector<double>* x_lower_idx = nullptr;
   const DenseVector<double>* x_upper = nullptr;
   const DenseVector<double>* x_upper_idx = nullptr;
   const DenseVector<int>* nnzRow = nullptr;

   /* set matrix */
   if (block_type == B_MAT) {
      mat = currBmat;
      assert(currBmatTrans);
   }
   else if (block_type == A_MAT) {
      assert(node != -1);
      mat = currAmat;
      assert(currAmatTrans);
   }
   else if (block_type == BL_MAT) {
      mat = currBlmat;
      assert(currBlmatTrans);
   }

   /* set variables */
   x_lower = (block_type == A_MAT || node == -1) ? currxlowParent : currxlowChild;
   x_lower_idx = (block_type == A_MAT || node == -1) ? currIxlowParent : currIxlowChild;
   x_upper = (block_type == A_MAT || node == -1) ? currxuppParent : currxuppChild;
   x_upper_idx = (block_type == A_MAT || node == -1) ? currIxuppParent : currIxuppChild;

   /* set reduction vectors */

   /* set non-zero row vectors */
   nnzRow = (block_type == BL_MAT) ? currNnzRowLink : currNnzRow;

   int n_elims = 0;

   const bool linking_row = (block_type == BL_MAT);
   const int node_row = linking_row ? -1 : node;
   const int node_col = (block_type == A_MAT) ? -1 : node;
   const SparseStorageDynamic* storage = mat;

   /* for every row in row in matrix */
   for (int r = 0; r < storage->n_rows(); r++) {
      double total_sum_modifications_row = 0.0;

      int start = storage->getRowPtr(r).start;
      int end = storage->getRowPtr(r).end;

      /* for every nonzero column in that row */
      for (int col_index = start; col_index < end; ++col_index) {
         const int col = storage->getJcolM(col_index);
         const double mat_abs = std::fabs(storage->getMat(col_index));

         /* remove all small entries */
         if (mat_abs < limit_min_mat_entry) {
            const INDEX row_INDEX(ROW, node_row, r, linking_row, system_type);
            const INDEX col_INDEX(COL, node_col, col);
            presolve_data.deleteEntryAtIndex(row_INDEX, col_INDEX, col_index);

            /* since the current entry got deleted we have to step back one entry */
            --col_index;
            --end;
            if (my_rank == 0 || !(node_row == -1 && node_col == -1))
               ++n_elims;
         }
         /* remove entries where their corresponding variables have valid lower and upper bounds, that overall do not have a real influence though */
         else if (!PIPSisZero((*x_upper_idx)[col]) && !PIPSisZero((*x_lower_idx)[col])) {
            const double bux = (*x_upper)[col];
            const double blx = (*x_lower)[col];

            /* don't remove entries in rows that need to be fixed */
            if (PIPSisEQ(bux, blx))
               continue;

            const int nnz = (*nnzRow)[r];
            assert(nnz != 0);

            if (mat_abs < limit_max_matrix_entry_impact && mat_abs * (bux - blx) * nnz < limit_matrix_entry_impact_feasdist * feastol) {
               const INDEX row_INDEX(ROW, node_row, r, linking_row, system_type);
               const INDEX col_INDEX(COL, node_col, col);
               presolve_data.deleteEntryAtIndex(row_INDEX, col_INDEX, col_index);

               /* since the current entry got deleted we have to step back one entry */
               --col_index;
               --end;
               if (my_rank == 0 || !(node_row == -1 && node_col == -1))
                  ++n_elims;
            }
               /* for linking constraints this is a slight modification of criterion three that does not require communication but is only a slight relaxation to
                * criterion two
                */
            else if ((block_type == BL_MAT && mat_abs * (bux - blx) * nnz < 1.0e-1 * feastol) ||
                     (block_type != BL_MAT && total_sum_modifications_row + mat_abs * (bux - blx) < 1.0e-1 * feastol / 2.0)) {
               total_sum_modifications_row += mat_abs * (bux - blx);
               const INDEX row_INDEX(ROW, node_row, r, linking_row, system_type);
               const INDEX col_INDEX(COL, node_col, col);
               presolve_data.deleteEntryAtIndex(row_INDEX, col_INDEX, col_index);

               /* since the current entry got deleted we have to step back one entry */
               --col_index;
               --end;
               if (my_rank == 0 || !(node_row == -1 && node_col == -1))
                  ++n_elims;
            }
         }
         /* not removed */
      }
   }

   return n_elims;
}


/* Go through columns and fix all empty ones to the current variables lower/upper bound (depending on objective)
 * Might detect unboundedness of problem.
 */
int StochPresolverModelCleanup::fixEmptyColumns() {
   int fixations = 0;

   for (int node = -1; node < nChildren; ++node) {
      if (presolve_data.nodeIsDummy(node))
         continue;

      updatePointersForCurrentNode(node, EQUALITY_SYSTEM);

      const DenseVector<double>& g = (node == -1) ? *currgParent : *currgChild;
      const DenseVector<double>& ixupp = (node == -1) ? *currIxuppParent : *currIxuppChild;
      const DenseVector<double>& ixlow = (node == -1) ? *currIxlowParent : *currIxlowChild;
      const DenseVector<double>& xupp = (node == -1) ? *currxuppParent : *currxuppChild;
      const DenseVector<double>& xlow = (node == -1) ? *currxlowParent : *currxlowChild;
      const DenseVector<int>& nnzs_col = (node == -1) ? *currNnzColParent : *currNnzColChild;

      for (int col_index = 0; col_index < nnzs_col.length(); ++col_index) {
         const INDEX col(COL, node, col_index);

         if (presolve_data.wasColumnRemoved(col)) {
            assert(nnzs_col[col_index] == 0);
            assert(PIPSisZero(ixlow[col_index]));
            assert(PIPSisZero(ixupp[col_index]));
            assert(PIPSisZero(xlow[col_index]));
            assert(PIPSisZero(xupp[col_index]));
            assert(PIPSisZero(g[col_index]));

            continue;
         }
         /* column fixation candidate */
         if (nnzs_col[col_index] == 0) {
            if (PIPSisLT(g[col_index], 0.0)) {
               if (!PIPSisZero(ixupp[col_index])) {
                  presolve_data.fixEmptyColumn(col, xupp[col_index]);
               }
               else {
                  PIPS_MPIabortInfeasible("Found empty column with non-zero objective vector and no bounds in objective direction! Unbounded!",
                        "StochPresolverModelCleanup.C", "fixEmptyColumns");
               }
            }
            else if (PIPSisLT(0.0, g[col_index])) {
               if (!PIPSisZero(ixlow[col_index])) {
                  presolve_data.fixEmptyColumn(col, xlow[col_index]);
               }
               else {
                  PIPS_MPIabortInfeasible("Found empty column with non-zero objective vector and no bounds in objective direction! Unbounded!",
                        "StochPresolverModelCleanup.C", "fixEmptyColumns");
               }
            }
            else {
               assert(PIPSisEQ(g[col_index], 0.0));
               if (!PIPSisZero(ixlow[col_index]))
                  presolve_data.fixEmptyColumn(col, xlow[col_index]);
               else if (!PIPSisZero(ixlow[col_index]))
                  presolve_data.fixEmptyColumn(col, xlow[col_index]);
               else
                  presolve_data.fixEmptyColumn(col, 0.0);
            }

            if (my_rank == 0 || !col.isLinkingCol())
               ++fixations;
         }
      }
   }

   return fixations;
}


