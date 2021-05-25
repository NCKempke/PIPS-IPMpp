#include "DistributedVariables.h"
#include "Vector.hpp"
#include "SmartPointer.h"
#include "Problem.h"
#include "DistributedVector.h"
#include "DistributedTree.h"
#include "DistributedTreeCallbacks.h"
#include <iostream>

DistributedVariables::DistributedVariables(const DistributedTree* tree, Vector<double>* x_in, Vector<double>* s_in, Vector<double>* y_in, Vector<double>* z_in,
      Vector<double>* v_in, Vector<double>* gamma_in, Vector<double>* w_in, Vector<double>* phi_in, Vector<double>* t_in, Vector<double>* lambda_in,
      Vector<double>* u_in, Vector<double>* pi_in, Vector<double>* ixlow_in, long long nxlowGlobal, Vector<double>* ixupp_in, long long nxuppGlobal,
      Vector<double>* iclow_in, long long mclowGlobal, Vector<double>* icupp_in, long long mcuppGlobal) : Variables() {

   stochNode = tree;

   SpReferTo(primals, x_in);
   SpReferTo(slacks, s_in);
   SpReferTo(equality_duals, y_in);
   SpReferTo(inequality_duals, z_in);
   SpReferTo(primal_lower_bound_gap, v_in);
   SpReferTo(primal_upper_bound_gap_dual, phi_in);
   SpReferTo(primal_upper_bound_gap, w_in);
   SpReferTo(primal_lower_bound_gap_dual, gamma_in);
   SpReferTo(slack_lower_bound_gap, t_in);
   SpReferTo(slack_lower_bound_gap_dual, lambda_in);
   SpReferTo(slack_upper_bound_gap, u_in);
   SpReferTo(slack_upper_bound_gap_dual, pi_in);
   SpReferTo(ixlow, ixlow_in);
   SpReferTo(ixupp, ixupp_in);
   SpReferTo(iclow, iclow_in);
   SpReferTo(icupp, icupp_in);

   nx = primals->length();
   my = equality_duals->length();
   mz = inequality_duals->length();

   assert(nx == ixlow->length() || 0 == ixlow->length());
   assert(nx == ixupp->length() || 0 == ixupp->length());
   assert(mz == iclow->length() || 0 == iclow->length());
   assert(mz == icupp->length() || 0 == icupp->length());

   nxlow = nxlowGlobal;
   nxupp = nxuppGlobal;
   mclow = mclowGlobal;
   mcupp = mcuppGlobal;
   number_complementarity_pairs = mclow + mcupp + nxlow + nxupp;

   assert(mz == slacks->length());
   assert(nx == primal_lower_bound_gap->length() || (0 == primal_lower_bound_gap->length() && nxlow == 0));
   assert(nx == primal_lower_bound_gap_dual->length() || (0 == primal_lower_bound_gap_dual->length() && nxlow == 0));

   assert(nx == primal_upper_bound_gap->length() || (0 == primal_upper_bound_gap->length() && nxupp == 0));
   assert(nx == primal_upper_bound_gap_dual->length() || (0 == primal_upper_bound_gap_dual->length() && nxupp == 0));

   assert(mz == slack_lower_bound_gap->length() || (0 == slack_lower_bound_gap->length() && mclow == 0));
   assert(mz == slack_lower_bound_gap_dual->length() || (0 == slack_lower_bound_gap_dual->length() && mclow == 0));

   assert(mz == slack_upper_bound_gap->length() || (0 == slack_upper_bound_gap->length() && mcupp == 0));
   assert(mz == slack_upper_bound_gap_dual->length() || (0 == slack_upper_bound_gap_dual->length() && mcupp == 0));

   createChildren();
}

DistributedVariables::DistributedVariables(const DistributedVariables& vars) : Variables(vars) {
   stochNode = vars.stochNode;
   for (auto i : vars.children) {
      children.push_back(new DistributedVariables(*i));
   }
}

DistributedVariables::~DistributedVariables() {
   for (auto & c : children)
      delete c;
}

void DistributedVariables::AddChild(DistributedVariables* child) {
   children.push_back(child);
}


void DistributedVariables::createChildren() {
   DistributedVector<double>& xst = dynamic_cast<DistributedVector<double>&>(*primals);
   DistributedVector<double>& sst = dynamic_cast<DistributedVector<double>&>(*slacks);
   DistributedVector<double>& yst = dynamic_cast<DistributedVector<double>&>(*equality_duals);
   DistributedVector<double>& zst = dynamic_cast<DistributedVector<double>&>(*inequality_duals);
   DistributedVector<double>& vst = dynamic_cast<DistributedVector<double>&>(*primal_lower_bound_gap);
   DistributedVector<double>& gammast = dynamic_cast<DistributedVector<double>&>(*primal_lower_bound_gap_dual);
   DistributedVector<double>& wst = dynamic_cast<DistributedVector<double>&>(*primal_upper_bound_gap);
   DistributedVector<double>& phist = dynamic_cast<DistributedVector<double>&>(*primal_upper_bound_gap_dual);
   DistributedVector<double>& tst = dynamic_cast<DistributedVector<double>&>(*slack_lower_bound_gap);
   DistributedVector<double>& lambdast = dynamic_cast<DistributedVector<double>&>(*slack_lower_bound_gap_dual);
   DistributedVector<double>& ust = dynamic_cast<DistributedVector<double>&>(*slack_upper_bound_gap);
   DistributedVector<double>& pist = dynamic_cast<DistributedVector<double>&>(*slack_upper_bound_gap_dual);
   DistributedVector<double>& ixlowst = dynamic_cast<DistributedVector<double>&>(*ixlow);
   DistributedVector<double>& ixuppst = dynamic_cast<DistributedVector<double>&>(*ixupp);
   DistributedVector<double>& iclowst = dynamic_cast<DistributedVector<double>&>(*iclow);
   DistributedVector<double>& icuppst = dynamic_cast<DistributedVector<double>&>(*icupp);


   for (size_t it = 0; it < xst.children.size(); it++) {
      AddChild(new DistributedVariables(stochNode->getChildren()[it], xst.children[it], sst.children[it], yst.children[it], zst.children[it],
            vst.children[it], gammast.children[it], wst.children[it], phist.children[it], tst.children[it], lambdast.children[it], ust.children[it],
            pist.children[it], ixlowst.children[it], nxlow, ixuppst.children[it], nxupp, iclowst.children[it], mclow, icuppst.children[it], mcupp));
   }

}

void
DistributedVariables::collapseHierarchicalStructure(const DistributedQP& hier_data, const DistributedTree* stochNode_, SmartPointer<Vector<double> > ixlow_,
      SmartPointer<Vector<double> > ixupp_, SmartPointer<Vector<double> > iclow_, SmartPointer<Vector<double> > icupp_) {
   dynamic_cast<DistributedVector<double>&>(*primals).collapseFromHierarchical(hier_data, *stochNode, VectorType::PRIMAL);

   dynamic_cast<DistributedVector<double>&>(*primal_lower_bound_gap).collapseFromHierarchical(hier_data, *stochNode, VectorType::PRIMAL);
   dynamic_cast<DistributedVector<double>&>(*primal_upper_bound_gap).collapseFromHierarchical(hier_data, *stochNode, VectorType::PRIMAL);
   dynamic_cast<DistributedVector<double>&>(*primal_upper_bound_gap_dual).collapseFromHierarchical(hier_data, *stochNode, VectorType::PRIMAL);
   dynamic_cast<DistributedVector<double>&>(*primal_lower_bound_gap_dual).collapseFromHierarchical(hier_data, *stochNode, VectorType::PRIMAL);

   dynamic_cast<DistributedVector<double>&>(*equality_duals).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Y);

   dynamic_cast<DistributedVector<double>&>(*slacks).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*inequality_duals).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*slack_lower_bound_gap).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*slack_upper_bound_gap).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*slack_upper_bound_gap_dual).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*slack_lower_bound_gap_dual).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);

   stochNode = stochNode_;

   ixlow = ixlow_;
   ixupp = ixupp_;
   iclow = iclow_;
   icupp = icupp_;

   for (size_t c = 0; c < children.size(); c++)
      delete children[c];

   children.clear();
   createChildren();
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
