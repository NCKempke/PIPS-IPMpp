#include "DistributedVariables.h"
#include "Vector.hpp"
#include "Problem.hpp"
#include "DistributedVector.h"
#include "DistributedTree.h"
#include "DistributedTreeCallbacks.h"
#include <iostream>
#include <utility>


DistributedVariables::DistributedVariables(const DistributedTree* tree, std::unique_ptr<Vector<double>> x_in, std::unique_ptr<Vector<double>> s_in,
   std::unique_ptr<Vector<double>> y_in, std::unique_ptr<Vector<double>> z_in, std::unique_ptr<Vector<double>> v_in,
   std::unique_ptr<Vector<double>> gamma_in, std::unique_ptr<Vector<double>> w_in, std::unique_ptr<Vector<double>> phi_in, std::unique_ptr<Vector<double>> t_in,
   std::unique_ptr<Vector<double>> lambda_in, std::unique_ptr<Vector<double>> u_in, std::unique_ptr<Vector<double>> pi_in,
   std::shared_ptr<Vector<double>> ixlow_in, long long nxlowGlobal, std::shared_ptr<Vector<double>> ixupp_in, long long nxuppGlobal,
   std::shared_ptr<Vector<double>> iclow_in, long long mclowGlobal, std::shared_ptr<Vector<double>> icupp_in, long long mcuppGlobal) :
   Variables(std::move(x_in), std::move(s_in), std::move(y_in), std::move(z_in), std::move(v_in), std::move(gamma_in),
      std::move(w_in), std::move(phi_in), std::move(t_in), std::move(lambda_in), std::move(u_in), std::move(pi_in),
      std::move(ixlow_in), std::move(ixupp_in), std::move(iclow_in), std::move(icupp_in)) {

   stochNode = tree;

   /* overwrite local values with distributed globals */
   nxlow = nxlowGlobal;
   nxupp = nxuppGlobal;
   mclow = mclowGlobal;
   mcupp = mcuppGlobal;
}

DistributedVariables::DistributedVariables(const DistributedVariables& vars) : Variables(vars) {
   stochNode = vars.stochNode;
}

std::unique_ptr<Variables> DistributedVariables::clone_full() const {
   return std::make_unique<DistributedVariables>(*this);
}

void
DistributedVariables::collapseHierarchicalStructure(const DistributedProblem& hier_data, const DistributedTree* stochNode_) {
   dynamic_cast<DistributedVector<double>&>(*primals).collapseFromHierarchical(hier_data, *stochNode,
      VectorType::PRIMAL);

   dynamic_cast<DistributedVector<double>&>(*primal_lower_bound_gap).collapseFromHierarchical(hier_data, *stochNode,
      VectorType::PRIMAL);
   dynamic_cast<DistributedVector<double>&>(*primal_upper_bound_gap).collapseFromHierarchical(hier_data, *stochNode,
      VectorType::PRIMAL);
   dynamic_cast<DistributedVector<double>&>(*primal_upper_bound_gap_dual).collapseFromHierarchical(hier_data,
      *stochNode, VectorType::PRIMAL);
   dynamic_cast<DistributedVector<double>&>(*primal_lower_bound_gap_dual).collapseFromHierarchical(hier_data,
      *stochNode, VectorType::PRIMAL);

   dynamic_cast<DistributedVector<double>&>(*equality_duals).collapseFromHierarchical(hier_data, *stochNode,
      VectorType::DUAL_Y);

   dynamic_cast<DistributedVector<double>&>(*slacks).collapseFromHierarchical(hier_data, *stochNode,
      VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*inequality_duals).collapseFromHierarchical(hier_data, *stochNode,
      VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*slack_lower_bound_gap).collapseFromHierarchical(hier_data, *stochNode,
      VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*slack_upper_bound_gap).collapseFromHierarchical(hier_data, *stochNode,
      VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*slack_upper_bound_gap_dual).collapseFromHierarchical(hier_data, *stochNode,
      VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*slack_lower_bound_gap_dual).collapseFromHierarchical(hier_data, *stochNode,
      VectorType::DUAL_Z);

   stochNode = stochNode_;
}


void DistributedVariables::update_indicators(std::shared_ptr<Vector<double>> ixlow_, std::shared_ptr<Vector<double>> ixupp_,
   std::shared_ptr<Vector<double>> iclow_, std::shared_ptr<Vector<double>> icupp_) {
   ixlow = std::move(ixlow_);
   ixupp = std::move(ixupp_);
   iclow = std::move(iclow_);
   icupp = std::move(icupp_);
}

void DistributedVariables::permuteVec0Entries(const std::vector<unsigned int>& perm, bool vars_only) {
   if (!vars_only) {
      dynamic_cast<DistributedVector<double>&>(*ixlow).permuteVec0Entries(perm);
      dynamic_cast<DistributedVector<double>&>(*ixupp).permuteVec0Entries(perm);
   }

   dynamic_cast<DistributedVector<double>&>(*primals).permuteVec0Entries(perm);
   dynamic_cast<DistributedVector<double>&>(*primal_lower_bound_gap).permuteVec0Entries(perm);
   dynamic_cast<DistributedVector<double>&>(*primal_upper_bound_gap).permuteVec0Entries(perm);
   dynamic_cast<DistributedVector<double>&>(*primal_upper_bound_gap_dual).permuteVec0Entries(perm);
   dynamic_cast<DistributedVector<double>&>(*primal_lower_bound_gap_dual).permuteVec0Entries(perm);
}

void DistributedVariables::permuteEqLinkingEntries(const std::vector<unsigned int>& perm) {
   dynamic_cast<DistributedVector<double>&>(*equality_duals).permuteLinkingEntries(perm);
}

void DistributedVariables::permuteIneqLinkingEntries(const std::vector<unsigned int>& perm, bool vars_only) {
   if (!vars_only) {
      dynamic_cast<DistributedVector<double>&>(*iclow).permuteLinkingEntries(perm);
      dynamic_cast<DistributedVector<double>&>(*icupp).permuteLinkingEntries(perm);
   }

   dynamic_cast<DistributedVector<double>&>(*slacks).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*inequality_duals).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*slack_lower_bound_gap).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*slack_upper_bound_gap).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*slack_upper_bound_gap_dual).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*slack_lower_bound_gap_dual).permuteLinkingEntries(perm);
}

bool DistributedVariables::isRootNodeInSync() const {
   bool in_sync = true;
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);
   if (!dynamic_cast<const DistributedVector<double>&>(*primals).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "x not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*slacks).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "s not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*equality_duals).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "y not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*inequality_duals).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "z not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*primal_lower_bound_gap).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "v not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*primal_lower_bound_gap_dual).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "gamma not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*primal_upper_bound_gap).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "w not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*primal_upper_bound_gap_dual).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "phi not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*slack_lower_bound_gap).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "t not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*slack_lower_bound_gap_dual).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "lambda not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*slack_upper_bound_gap).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "u not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*slack_upper_bound_gap_dual).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "pi not in sync" << std::endl;
      in_sync = false;
   }

   MPI_Barrier(MPI_COMM_WORLD);
   return in_sync;
}