/*
 * StochPresolverColumnFixation.h
 *
 *  Created on: 15.05.2019
 *      Author: bzfkempk
 */
#include "StochPresolverColumnFixation.h"

#include "pipsdef.h"
#include "PIPSIPMppOptions.h"
#include <cmath>

StochPresolverColumnFixation::StochPresolverColumnFixation(PresolveData& presolve_data, const DistributedQP& origProb) : StochPresolverBase(
      presolve_data, origProb), limit_fixation_max_fixing_impact(pipsipmpp_options::get_double_parameter("PRESOLVE_COLUMN_FIXATION_MAX_FIXING_IMPACT")),
      fixed_columns(0) {
}

StochPresolverColumnFixation::~StochPresolverColumnFixation() {
}

/* scan through columns and fix those that have tight bounds */
bool StochPresolverColumnFixation::applyPresolving() {
   assert(presolve_data.reductionsEmpty());
   assert(presolve_data.presolve_dataInSync());

#ifndef NDEBUG
   if (my_rank == 0 && verbosity > 1) {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << "\n";
      std::cout << "--- Before column fixation presolving:" << "\n";
   }
   countRowsCols();
#endif

   presolve_data.startColumnFixation();

   int fixed_columns_run = 0;

   /* remove fixed columns from system */
   updatePointersForCurrentNode(-1, EQUALITY_SYSTEM);

   /* linking variables */
   for (int col = 0; col < currgParent->length(); ++col) {
      const INDEX col_index = INDEX(COL, -1, col);

      if (presolve_data.wasColumnRemoved(col_index))
         continue;

      // TODO : make criterion numerically more stable
      // if the column is fixed to a value I want abs(xlow - xupp) to be low enough such that
      // in the whole column for each aij != 0 abs(xlow - xupp) * abs(aij) < epsilon <= feastol
      // this way we ensure that even if the fixing is not quite accurate the column stays valid
      // TODO : include also objective coefficient into this criterion ? - so that fixing the column somehow inaccurately will not majorly impact the objective value
      double absmax_row = 1.0;

      if (!PIPSisZero((*currIxlowParent)[col]) && !PIPSisZero((*currIxuppParent)[col])) {

         assert(PIPSisLE(0.0, (*currxuppParent)[col] - (*currxlowParent)[col]));

         if (PIPSisLT(fabs((*currxuppParent)[col] - (*currxlowParent)[col]) * absmax_row, limit_fixation_max_fixing_impact)) {
            // verify if one of the bounds is integer:
            double intpart = 0.0;
            double value = 0.0;
            if (std::modf((*currxlowParent)[col], &intpart) == 0.0)
               value = (*currxlowParent)[col];
            else if (std::modf((*currxuppParent)[col], &intpart) == 0.0)
               value = (*currxuppParent)[col];
            else  // set the variable to the arithmetic mean:
               value = ((*currxlowParent)[col] + (*currxuppParent)[col]) / 2.0;

            presolve_data.fixColumn(col_index, value);
            if (my_rank == 0)
               ++fixed_columns_run;
         }
      }
   }

   /* child nodes */
   for (int node = 0; node < nChildren; ++node) {
      if (presolve_data.nodeIsDummy(node))
         continue;

      updatePointersForCurrentNode(node, EQUALITY_SYSTEM);

      /* linking variables */
      for (int col = 0; col < currgChild->length(); ++col) {
         // TODO : make criterion numerically more stable
         // if the column is fixed to a value I want abs(xlow - xupp) to be low enough such that
         // in the whole column for each aij != 0 abs(xlow - xupp) * abs(aij) < epsilon <= feastol
         // this way we ensure that even if the fixing is not quite accurate the column stays valid
         // TODO : include also objective coefficient into this criterion ?
         double absmax_row = 1.0;

         if (!PIPSisZero((*currIxlowChild)[col]) && !PIPSisZero((*currIxuppChild)[col])) {
            assert(PIPSisLE(0.0, (*currxuppChild)[col] - (*currxlowChild)[col]));

            if (PIPSisLT(((*currxuppChild)[col] - (*currxlowChild)[col]) * absmax_row, limit_fixation_max_fixing_impact)) {
               // verify if one of the bounds is integer:
               double intpart = 0.0;
               double value = 0.0;
               if (std::modf((*currxlowChild)[col], &intpart) == 0.0)
                  value = (*currxlowChild)[col];
               else if (std::modf((*currxuppChild)[col], &intpart) == 0.0)
                  value = (*currxuppChild)[col];
               else  // set the variable to the arithmetic mean:
                  value = ((*currxlowChild)[col] + (*currxuppChild)[col]) / 2.0;
               presolve_data.fixColumn(INDEX(COL, node, col), value);
               ++fixed_columns_run;
            }
         }
      }
   }

   /* communicate the local changes */
   presolve_data.allreduceAndApplyBoundChanges();
   presolve_data.allreduceAndApplyLinkingRowActivities();
   presolve_data.allreduceAndApplyNnzChanges();
   presolve_data.allreduceObjOffset();

   PIPS_MPIgetSumInPlace(fixed_columns_run, MPI_COMM_WORLD);
   fixed_columns += fixed_columns_run;

#ifndef NDEBUG
   if (my_rank == 0 && verbosity > 1)
      std::cout << "\tFixed columns during column fixation: " << fixed_columns << "\n";
   else if (my_rank == 0 && verbosity == 1)
      std::cout << "Colfix:\t removed " << fixed_columns << " columns" << "\n";

   if (my_rank == 0 && verbosity > 1)
      std::cout << "--- After column fixation presolving:" << "\n";
   countRowsCols();
   if (my_rank == 0 && verbosity > 1)
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << "\n";
#endif

   assert(presolve_data.reductionsEmpty());
   assert(presolve_data.presolve_dataInSync());

   if (fixed_columns_run != 0)
      return true;
   else
      return false;
}
