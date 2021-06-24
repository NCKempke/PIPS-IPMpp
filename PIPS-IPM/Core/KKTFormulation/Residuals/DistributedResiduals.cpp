#include "DistributedResiduals.hpp"

#include <utility>
#include "DistributedTree.h"
#include "DistributedVector.h"

DistributedResiduals::DistributedResiduals(std::unique_ptr<Vector<double>> rQ_, std::unique_ptr<Vector<double>> rA_,
   std::unique_ptr<Vector<double>> rC_, std::unique_ptr<Vector<double>> rz_, std::unique_ptr<Vector<double>> rt_,
   std::unique_ptr<Vector<double>> rlambda_, std::unique_ptr<Vector<double>> ru_, std::unique_ptr<Vector<double>> rpi_,
   std::unique_ptr<Vector<double>> rv_, std::unique_ptr<Vector<double>> rgamma_, std::unique_ptr<Vector<double>> rw_,
   std::unique_ptr<Vector<double>> rphi_, std::shared_ptr<Vector<double>> ixlow_,
   std::shared_ptr<Vector<double>> ixupp_, std::shared_ptr<Vector<double>> iclow_,
   std::shared_ptr<Vector<double>> icupp_) : Residuals(
   std::move(rQ_), std::move(rA_), std::move(rC_), std::move(rz_), std::move(rt_), std::move(rlambda_), std::move(ru_),
   std::move(rpi_), std::move(rv_), std::move(rgamma_), std::move(rw_), std::move(rphi_), std::move(ixlow_), std::move(ixupp_), std::move(iclow_),
   std::move(icupp_)) {}

std::unique_ptr<Residuals> DistributedResiduals::cloneFull() const{
   return std::make_unique<DistributedResiduals>(*this);
};

void
DistributedResiduals::collapse_hierarchical_structure(const DistributedProblem& data_hier, const DistributedTree* tree_hier,
   std::shared_ptr<Vector<double>> ixlow_,
   std::shared_ptr<Vector<double>> ixupp_, std::shared_ptr<Vector<double>> iclow_,
   std::shared_ptr<Vector<double>> icupp_) {
   dynamic_cast<DistributedVector<double>&>(*lagrangian_gradient).collapseFromHierarchical(data_hier, *tree_hier,
      VectorType::PRIMAL);

   const bool empty_vec = true;
   if (nxlow > 0) {
      dynamic_cast<DistributedVector<double>&>(*rv).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);
      dynamic_cast<DistributedVector<double>&>(*rgamma).collapseFromHierarchical(data_hier, *tree_hier,
         VectorType::PRIMAL);
   } else {
      dynamic_cast<DistributedVector<double>&>(*rv).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL,
         empty_vec);
      dynamic_cast<DistributedVector<double>&>(*rgamma).collapseFromHierarchical(data_hier, *tree_hier,
         VectorType::PRIMAL, empty_vec);
   }

   if (nxupp > 0) {
      dynamic_cast<DistributedVector<double>&>(*rw).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);
      dynamic_cast<DistributedVector<double>&>(*rphi).collapseFromHierarchical(data_hier, *tree_hier,
         VectorType::PRIMAL);
   } else {
      dynamic_cast<DistributedVector<double>&>(*rw).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL,
         empty_vec);
      dynamic_cast<DistributedVector<double>&>(*rphi).collapseFromHierarchical(data_hier, *tree_hier,
         VectorType::PRIMAL, empty_vec);
   }

   dynamic_cast<DistributedVector<double>&>(*equality_residuals).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Y);

   dynamic_cast<DistributedVector<double>&>(*inequality_residuals).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*inequality_dual_residuals).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);

   if (mcupp > 0) {
      dynamic_cast<DistributedVector<double>&>(*ru).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
      dynamic_cast<DistributedVector<double>&>(*rpi).collapseFromHierarchical(data_hier, *tree_hier,
         VectorType::DUAL_Z);
   } else {
      dynamic_cast<DistributedVector<double>&>(*ru).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z,
         empty_vec);
      dynamic_cast<DistributedVector<double>&>(*rpi).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z,
         empty_vec);
   }

   if (mclow > 0) {
      dynamic_cast<DistributedVector<double>&>(*rt).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
      dynamic_cast<DistributedVector<double>&>(*rlambda).collapseFromHierarchical(data_hier, *tree_hier,
         VectorType::DUAL_Z);
   } else {
      dynamic_cast<DistributedVector<double>&>(*rt).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z,
         empty_vec);
      dynamic_cast<DistributedVector<double>&>(*rlambda).collapseFromHierarchical(data_hier, *tree_hier,
         VectorType::DUAL_Z, empty_vec);
   }


   ixlow = std::move(ixlow_);
   ixupp = std::move(ixupp_);
   iclow = std::move(iclow_);
   icupp = std::move(icupp_);
}

void DistributedResiduals::permute_vec_0_entries(const std::vector<unsigned int>& perm, bool resids_only) {
   if (!resids_only) {
      dynamic_cast<DistributedVector<double>&>(*ixlow).permuteVec0Entries(perm);
      dynamic_cast<DistributedVector<double>&>(*ixupp).permuteVec0Entries(perm);
   }

   dynamic_cast<DistributedVector<double>&>(*lagrangian_gradient).permuteVec0Entries(perm);

   if (nxlow > 0) {
      dynamic_cast<DistributedVector<double>&>(*rv).permuteVec0Entries(perm);
      dynamic_cast<DistributedVector<double>&>(*rgamma).permuteVec0Entries(perm);
   }

   if (nxupp > 0) {
      dynamic_cast<DistributedVector<double>&>(*rw).permuteVec0Entries(perm);
      dynamic_cast<DistributedVector<double>&>(*rphi).permuteVec0Entries(perm);
   }
}

void DistributedResiduals::permute_eq_linking_entries(const std::vector<unsigned int>& perm) {
   dynamic_cast<DistributedVector<double>&>(*equality_residuals).permuteLinkingEntries(perm);
}

void DistributedResiduals::permute_ineq_linking_entries(const std::vector<unsigned int>& perm, bool resids_only) {
   if (!resids_only) {
      dynamic_cast<DistributedVector<double>&>(*iclow).permuteLinkingEntries(perm);
      dynamic_cast<DistributedVector<double>&>(*icupp).permuteLinkingEntries(perm);
   }

   dynamic_cast<DistributedVector<double>&>(*inequality_residuals).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*inequality_dual_residuals).permuteLinkingEntries(perm);

   if (mcupp > 0) {
      dynamic_cast<DistributedVector<double>&>(*ru).permuteLinkingEntries(perm);
      dynamic_cast<DistributedVector<double>&>(*rpi).permuteLinkingEntries(perm);
   }

   if (mclow > 0) {
      dynamic_cast<DistributedVector<double>&>(*rt).permuteLinkingEntries(perm);
      dynamic_cast<DistributedVector<double>&>(*rlambda).permuteLinkingEntries(perm);
   }
}

bool DistributedResiduals::is_root_node_in_sync() const {
   bool in_sync = true;
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);
   if (!dynamic_cast<const DistributedVector<double>&>(*lagrangian_gradient).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rQ not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*inequality_residuals).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rC not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*equality_residuals).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rA not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*inequality_dual_residuals).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rz not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*rt).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rt not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*rlambda).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rlambda not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*ru).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "ru not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*rpi).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rpi not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*rv).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rv not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*rgamma).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rgamma not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*rw).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rw not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*rphi).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rphi not in sync" << std::endl;
      in_sync = false;
   }

   MPI_Barrier(MPI_COMM_WORLD);
   return in_sync;

}
