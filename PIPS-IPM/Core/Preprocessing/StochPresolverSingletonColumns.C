/*
 * StochPresolverSingletonColumns.C
 *
 *  Created on: 08.05.2018
 *      Author: bzfuslus
 */

// TODO : distinguish between fixing and removals when printing the stats

#include "StochPresolverSingletonColumns.h"
#include "DistributedVectorUtilities.h"

StochPresolverSingletonColumns::StochPresolverSingletonColumns(PresolveData& presolve_data, const DistributedProblem& origProb) : StochPresolverBase(
      presolve_data, origProb), removed_cols(0), local_singletons(false),
      n_linking_rows_eq(dynamic_cast<const DistributedVector<double>&>(*origProb.equality_rhs).last->length()),
      n_linking_rows_ineq(dynamic_cast<const DistributedVector<double>&>(*origProb.inequality_upper_bounds).last->length()),
      local_linking_column_for_row_in_proc(n_linking_rows_eq + n_linking_rows_ineq), cols(n_linking_rows_eq + n_linking_rows_ineq),
      coeffs(n_linking_rows_ineq) {
}

bool StochPresolverSingletonColumns::applyPresolving() {
   assert(presolve_data.reductionsEmpty());
   assert(presolve_data.presolve_dataInSync());

   resetArrays();
   presolve_data.startSingletonColumnPresolve();

#ifndef NDEBUG
   if (my_rank == 0 && verbosity > 1) {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
      std::cout << "--- Before singleton columns presolving:\n";
   }
   countRowsCols();
#endif

#ifndef NDEBUG
   if (presolve_data.getSingletonCols().empty()) {
      if (my_rank == 0 && verbosity > 1)
         std::cout << "No more singletons left - exiting\n";
   }
#endif

   int removed_cols_run = 0;

   // main loop:
   while (!presolve_data.getSingletonCols().empty()) {
      const INDEX& col = presolve_data.getSingletonCols().front();
      const bool removed = removeSingletonColumn(col);

      if (removed)
         ++removed_cols_run;

      presolve_data.getSingletonCols().pop();
   }

   /* check for local singletons and communicate */
   PIPS_MPIgetLogicOrInPlace(local_singletons, MPI_COMM_WORLD);
   if (local_singletons) {
      /* allreduce local singleton columns and coeffs */

      // TODO: change this - it should be the process that provides most numerical stability..
      /* allreduce the procs that found a local singleton column in a linking row - the lowest ranking one will get to remove the column */
      PIPS_MPImaxArrayInPlace(local_linking_column_for_row_in_proc, MPI_COMM_WORLD);

      /* remove local singleton columns */
      for (unsigned int i = 0; i < local_linking_column_for_row_in_proc.size(); ++i) {
         const int proc_that_removes = local_linking_column_for_row_in_proc[i];

         if (proc_that_removes == -1)
            continue;

         const SystemType system_type = (i < n_linking_rows_eq) ? EQUALITY_SYSTEM : INEQUALITY_SYSTEM;
         const int row_index = (i < n_linking_rows_eq) ? i : i - n_linking_rows_eq;
         const INDEX row(ROW, -1, row_index, true, system_type);
         if (my_rank == proc_that_removes) {
            assert(!cols[i].isLinkingCol());
            if (row.inEqSys())
               presolve_data.removeImpliedFreeColumnSingletonEqualityRowSynced(row, cols[i]);
            else {
               const double iclow = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_lower_bound_indicators, row);
               const double icupp = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_upper_bound_indicators, row);
               const double clow = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_lower_bounds, row);
               const double cupp = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_upper_bounds, row);

               if (!PIPSisZero(iclow) && !PIPSisZero(icupp) && PIPSisEQ(clow, cupp))
                  presolve_data.removeImpliedFreeColumnSingletonEqualityRowSynced(row, cols[i]);
               else
                  presolve_data.removeFreeColumnSingletonInequalityRowSynced(row, cols[i], coeffs.at(row_index));
            }
            ++removed_cols_run;
         }
         else {
            if (row.inEqSys())
               presolve_data.removeImpliedFreeColumnSingletonEqualityRowSynced(row, INDEX());
            else {
               const double iclow = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_lower_bound_indicators, row);
               const double icupp = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_upper_bound_indicators, row);
               const double clow = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_lower_bounds, row);
               const double cupp = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_upper_bounds, row);

               if (!PIPSisZero(iclow) && !PIPSisZero(icupp) && PIPSisEQ(clow, cupp))
                  presolve_data.removeImpliedFreeColumnSingletonEqualityRowSynced(row, INDEX());
               else
                  presolve_data.removeFreeColumnSingletonInequalityRowSynced(row, INDEX(), 0.0);
            }
         }
      }
   }

   presolve_data.allreduceAndApplyBoundChanges();
   presolve_data.allreduceAndApplyNnzChanges();
   presolve_data.allreduceAndApplyLinkingRowActivities();
   presolve_data.allreduceLinkingVarBounds();
   presolve_data.allreduceAndApplyObjVecChanges();
   presolve_data.allreduceObjOffset();

   PIPS_MPIgetSumInPlace(removed_cols_run, MPI_COMM_WORLD);
   removed_cols += removed_cols_run;

#ifndef NDEBUG
   if (my_rank == 0 && verbosity > 1) {
      std::cout << "--- After singleton columns presolving:\n";
      std::cout << "\tRemoved columns during singleton column elimination: " << removed_cols << "\n";
   }
   else if (my_rank == 0 && verbosity == 1)
      std::cout << "SinCol:\t removed " << removed_cols << " cols\n";

   countRowsCols();
   if (my_rank == 0 && verbosity > 1)
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
#endif

   assert(presolve_data.reductionsEmpty());
   assert(presolve_data.presolve_dataInSync());

   if (removed_cols_run != 0)
      return true;
   else
      return false;
}

bool StochPresolverSingletonColumns::removeSingletonColumn(const INDEX& col) {
   assert(col.isCol());
   assert(col.hasValidNode(nChildren));
   assert(!presolve_data.nodeIsDummy(col.getNode()));

   updatePointersForCurrentNode(col.getNode(), EQUALITY_SYSTEM);
   const DenseVector<int>& nnzs_col = col.isLinkingCol() ? *currNnzColParent : *currNnzColChild;

   if (nnzs_col[col.getIndex()] == 0 || presolve_data.wasColumnRemoved(col))
      return false;

   assert(nnzs_col[col.getIndex()] == 1);

   bool found = false;
   /* find the associated row via checking transposed matrices */
   INDEX row = findRowForColumnSingleton(col, found);

   /* the singleton is a linking variable located on another process */
   if (!found) {
      assert(row.isEmpty());
      assert(col.isLinkingCol());
      return false;
   }

   assert(row.isRow());
   assert(row.hasValidNode(nChildren));

   if (row.isLinkingRow())
      updatePointersForCurrentNode(col.getNode(), row.getSystemType());
   else
      updatePointersForCurrentNode(row.getNode(), row.getSystemType());

   /* check whether col is free/implied free */
   bool lb_implied_free = false;
   bool ub_implied_free = false;

   // TODO : we should probably not do anything if the coefficient of the column we are looking at is too small
   checkColImpliedFree(col, row, lb_implied_free, ub_implied_free);
   bool implied_free = lb_implied_free && ub_implied_free;

   /* equality singleton variables */
   if (row.inEqSys() && implied_free) {
      /* store local singleton cols for later */
      if (row.isLinkingRow() && !col.isLinkingCol()) {
         local_singletons = true;

         const int index_row = row.getIndex();
         local_linking_column_for_row_in_proc.at(index_row) = my_rank;
         cols.at(index_row) = col;

         return false;
      }
      else
         presolve_data.removeImpliedFreeColumnSingletonEqualityRow(row, col);
   }
   else if (row.inInEqSys()) {
      /* inequality singleton variables */
      const double iclow = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_lower_bound_indicators, row);
      const double icupp = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_upper_bound_indicators, row);
      const double clow = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_lower_bounds, row);
      const double cupp = getSimpleVecFromRowStochVec(*presolve_data.getPresProb().inequality_upper_bounds, row);

      if (!PIPSisZero(iclow) && !PIPSisZero(icupp) && PIPSisEQ(clow, cupp)) {
         if (row.isLinkingRow() && !col.isLinkingCol()) {
            local_singletons = true;

            const int index_row = n_linking_rows_eq + row.getIndex();
            local_linking_column_for_row_in_proc.at(index_row) = my_rank;
            cols.at(index_row) = col;
            return false;
         }
         else
            presolve_data.removeImpliedFreeColumnSingletonEqualityRow(row, col);

      }
      else if (!PIPSisZero(iclow) && !PIPSisZero(icupp)) {
         assert(!PIPSisEQ(clow, cupp));
         return false;
      }
      else {
         assert(!PIPSisEQ(iclow, icupp));

         const double obj_coeff = getSimpleVecFromColStochVec(*presolve_data.getPresProb().objective_gradient, col);
         const double coeff = presolve_data.getRowCoeff(row, col);

         assert(!PIPSisZero(coeff));

         /* convert row to less equal row */
         const double coeff_le_row = PIPSisZero(iclow) ? coeff : -coeff;

         /* both positive */
         if (PIPSisLE(0.0, coeff_le_row) && PIPSisLE(0.0, obj_coeff)) {
            const double ixlow = getSimpleVecFromColStochVec(*presolve_data.getPresProb().primal_lower_bound_indicators, col);
            const double xlow = getSimpleVecFromColStochVec(*presolve_data.getPresProb().primal_lower_bounds, col);
            if (PIPSisZero(ixlow) && PIPSisLT(0.0, obj_coeff))
               PIPS_MPIabortInfeasible("Found unbounded singleton column variable", "StochPresolverSingletonColumns.C", "removeSingletonColumn");

            if (PIPSisZero(ixlow) && PIPSisZero(obj_coeff)) {
               if (row.isLinkingRow() && !col.isLinkingCol()) {
                  local_singletons = true;

                  const int index_row = n_linking_rows_eq + row.getIndex();
                  local_linking_column_for_row_in_proc.at(index_row) = my_rank;
                  cols.at(index_row) = col;
                  coeffs.at(row.getIndex()) = coeff;
                  return false;
               }
               else
                  /* remove variable and associated row from whole system */
                  presolve_data.removeFreeColumnSingletonInequalityRow(row, col, coeff);
            }
            else {
               /* fix variable to lower bound */
               assert(xlow != INF_NEG);
               presolve_data.fixColumnInequalitySingleton(col, row, xlow, coeff);
            }
         }
            /* both negative */
         else if (PIPSisLE(coeff_le_row, 0.0) && PIPSisLE(obj_coeff, 0.0)) {
            const double ixupp = getSimpleVecFromColStochVec(*presolve_data.getPresProb().primal_upper_bound_indicators, col);
            const double xupp = getSimpleVecFromColStochVec(*presolve_data.getPresProb().primal_upper_bounds, col);
            if (PIPSisZero(ixupp) && PIPSisLT(obj_coeff, 0.0))
               PIPS_MPIabortInfeasible("Found unbounded singleton column variable", "StochPresolverSingletonColumns.C", "removeSingletonColumn");

            if (PIPSisZero(ixupp) && PIPSisZero(obj_coeff)) {
               if (row.isLinkingRow() && !col.isLinkingCol()) {
                  local_singletons = true;

                  const int index_row = n_linking_rows_eq + row.getIndex();
                  local_linking_column_for_row_in_proc.at(index_row) = my_rank;
                  cols.at(index_row) = col;
                  coeffs.at(row.getIndex()) = coeff;
                  return false;
               }
               else
                  /* remove variable and associated row from whole system */
                  presolve_data.removeFreeColumnSingletonInequalityRow(row, col, coeff);
            }
            else {
               /* fix variable to upper bound */
               assert(xupp != INF_POS);
               presolve_data.fixColumnInequalitySingleton(col, row, xupp, coeff);
            }
         }
      }
   }
   else
      return false;

   const bool at_root = (row.getNode() == -1 && col.isLinkingCol());
   if (!at_root || my_rank == 0)
      return true;
   else
      return false;

}

INDEX StochPresolverSingletonColumns::findRowForColumnSingleton(const INDEX& col, bool& found) {
   assert(col.isCol());
   assert(-1 <= col.getNode() && col.getNode() < nChildren);

   if (col.getNode() == -1)
      return findRowForLinkingSingleton(col.getIndex(), found);
   else
      return findRowForNonlinkingSingelton(col, found);
}

INDEX StochPresolverSingletonColumns::findRowForLinkingSingleton(int col, bool& found) {
   /* go through all children and check linking variables for the singleton row */

   /* equality part */
   INDEX row = findRowForLinkingSingletonInSystem(col, EQUALITY_SYSTEM, found);
   if (found)
      return row;

   /* inequality part */
   return findRowForLinkingSingletonInSystem(col, INEQUALITY_SYSTEM, found);
}

INDEX StochPresolverSingletonColumns::findRowForNonlinkingSingelton(const INDEX& col, bool& found) {
   assert(col.isCol());
   assert(col.getNode() != -1);

   if (presolve_data.nodeIsDummy(col.getNode())) {
      found = false;
      return INDEX();
   }

   /* check Bmat and Blmat for the singleton row */
   SystemType row_system_type = EQUALITY_SYSTEM;
   bool row_linking = false;
   const int row_node = col.getNode();
   int row_index = -1;

   /* equality part */
   updatePointersForCurrentNode(row_node, row_system_type);

   /* Bmat */
   assert(currBmatTrans);

   found = findRowForSingletonColumnInMatrix(*currBmatTrans, row_index, col.getIndex());
   if (found)
      return INDEX(ROW, row_node, row_index, row_linking, row_system_type);

   /* Blmat */
   assert(currBlmatTrans);
   found = findRowForSingletonColumnInMatrix(*currBlmatTrans, row_index, col.getIndex());
   if (found) {
      row_linking = true;
      return INDEX(ROW, -1, row_index, row_linking, row_system_type);
   }

   /* inequality part */
   row_system_type = INEQUALITY_SYSTEM;
   updatePointersForCurrentNode(row_node, row_system_type);

   /* Bmat */
   assert(currBmatTrans);
   found = findRowForSingletonColumnInMatrix(*currBmatTrans, row_index, col.getIndex());
   if (found)
      return INDEX(ROW, row_node, row_index, row_linking, row_system_type);

   /* Blmat */
   assert(currBlmatTrans);
   found = findRowForSingletonColumnInMatrix(*currBlmatTrans, row_index, col.getIndex());
   if (found) {
      row_linking = true;
      return INDEX(ROW, -1, row_index, row_linking, row_system_type);
   }

   return INDEX();
}

INDEX StochPresolverSingletonColumns::findRowForLinkingSingletonInSystem(int col, SystemType system_type, bool& found) {

   for (int row_node = -1; row_node < nChildren; ++row_node) {
      if (!presolve_data.nodeIsDummy(row_node)) {
         int row_index = -1;
         updatePointersForCurrentNode(row_node, system_type);

         /* check transposed for entry */
         if (row_node == -1) {
            /* Bmat */
            assert(currBmatTrans);
            found = findRowForSingletonColumnInMatrix(*currBmatTrans, row_index, col);
            if (found)
               return INDEX(ROW, row_node, row_index, false, system_type);

            /* Blmat */
            assert(currBlmatTrans);
            found = findRowForSingletonColumnInMatrix(*currBlmatTrans, row_index, col);
            if (found)
               return INDEX(ROW, row_node, row_index, true, system_type);
         }
         else {
            /* Amat */
            assert(currAmatTrans);
            found = findRowForSingletonColumnInMatrix(*currAmatTrans, row_index, col);
            if (found)
               return INDEX(ROW, row_node, row_index, false, system_type);
         }
      }
   }
   assert(found == false);
   return INDEX();
}

bool StochPresolverSingletonColumns::findRowForSingletonColumnInMatrix(const SparseStorageDynamic& mat, int& row, const int& col) {
   assert(col < mat.n_rows());
   if (mat.getRowPtr(col).start != mat.getRowPtr(col).end) {
      assert((mat.getRowPtr(col).end - mat.getRowPtr(col).start) == 1);
      row = mat.getJcolM(mat.getRowPtr(col).start);
      return true;
   }
   return false;
}

void StochPresolverSingletonColumns::checkColImpliedFree(const INDEX& col, const INDEX& row, bool& lb_implied_free, bool& ub_implied_free) {
   assert(col.isCol());
   assert(row.isRow());

   ub_implied_free = lb_implied_free = false;

   if (row.isLinkingRow())
      updatePointersForCurrentNode(col.getNode(), row.getSystemType());
   else
      updatePointersForCurrentNode(row.getNode(), row.getSystemType());

   /* check whether free */
   const double ixupp = getSimpleVecFromColStochVec(*presolve_data.getPresProb().primal_upper_bound_indicators, col);
   const double ixlow = getSimpleVecFromColStochVec(*presolve_data.getPresProb().primal_lower_bound_indicators, col);

   if (PIPSisZero(ixupp))
      ub_implied_free = true;
   if (PIPSisZero(ixlow))
      lb_implied_free = true;

   /* check whether bound tightening found bounds from the variables row that make it implied free */
   ub_implied_free = ub_implied_free || presolve_data.varBoundImpliedFreeBy(true, col, row);
   lb_implied_free = lb_implied_free || presolve_data.varBoundImpliedFreeBy(false, col, row);
}

void StochPresolverSingletonColumns::resetArrays() {
   std::fill(local_linking_column_for_row_in_proc.begin(), local_linking_column_for_row_in_proc.end(), -1);
   std::fill(cols.begin(), cols.end(), INDEX());
   std::fill(coeffs.begin(), coeffs.end(), INF_POS);
}
