#include "DistributedResiduals.hpp"
#include "DistributedTree.h"
#include "DistributedVector.h"

DistributedResiduals::DistributedResiduals(Vector<double>* rQ_, Vector<double>* rA_, Vector<double>* rC_, Vector<double>* rz_, Vector<double>* rt_,
      Vector<double>* rlambda_, Vector<double>* ru_, Vector<double>* rpi_, Vector<double>* rv_, Vector<double>* rgamma_, Vector<double>* rw_,
      Vector<double>* rphi_, Vector<double>* ixlow_, double nxlowGlobal, Vector<double>* ixupp_, double nxuppGlobal, Vector<double>* iclow_,
      double mclowGlobal, Vector<double>* icupp_, double mcuppGlobal) {
   SpReferTo(ixlow, ixlow_);
   nxlow = nxlowGlobal;

   SpReferTo(ixupp, ixupp_);
   nxupp = nxuppGlobal;

   SpReferTo(iclow, iclow_);
   mclow = mclowGlobal;

   SpReferTo(icupp, icupp_);
   mcupp = mcuppGlobal;

   SpReferTo(lagrangian_gradient, rQ_);
   SpReferTo(rA, rA_);
   SpReferTo(rC, rC_);
   SpReferTo(rz, rz_);
   SpReferTo(rt, rt_);

   SpReferTo(rlambda, rlambda_);
   SpReferTo(ru, ru_);
   SpReferTo(rpi, rpi_);
   SpReferTo(rv, rv_);
   SpReferTo(rgamma, rgamma_);
   SpReferTo(rw, rw_);
   SpReferTo(rphi, rphi_);

   createChildren();
}


DistributedResiduals::DistributedResiduals(const DistributedTree* tree, Vector<double>* ixlow_, Vector<double>* ixupp_, Vector<double>* iclow_,
      Vector<double>* icupp_) {

   SpReferTo(ixlow, ixlow_);
   nxlow = ixlow->numberOfNonzeros();

   SpReferTo(ixupp, ixupp_);
   nxupp = ixupp->numberOfNonzeros();

   SpReferTo(iclow, iclow_);
   mclow = iclow->numberOfNonzeros();

   SpReferTo(icupp, icupp_);
   mcupp = icupp->numberOfNonzeros();

   const bool empty_vector = true;
   lagrangian_gradient = SmartPointer<Vector<double> >((Vector<double>*) tree->new_primal_vector());
   rA = SmartPointer<Vector<double> >((Vector<double>*) tree->newDualYVector());
   rC = SmartPointer<Vector<double> >((Vector<double>*) tree->newDualZVector());

   rz = SmartPointer<Vector<double> >((Vector<double>*) tree->newDualZVector());
   if (mclow > 0) {
      rt = SmartPointer<Vector<double> >((Vector<double>*) tree->newDualZVector());
      rlambda = SmartPointer<Vector<double> >((Vector<double>*) tree->newDualZVector());
   }
   else {
      rt = SmartPointer<Vector<double> >((Vector<double>*) tree->newDualZVector(empty_vector));
      rlambda = SmartPointer<Vector<double> >((Vector<double>*) tree->newDualZVector(empty_vector));
   }

   if (mcupp > 0) {
      ru = SmartPointer<Vector<double> >((Vector<double>*) tree->newDualZVector());
      rpi = SmartPointer<Vector<double> >((Vector<double>*) tree->newDualZVector());
   }
   else {
      ru = SmartPointer<Vector<double> >((Vector<double>*) tree->newDualZVector(empty_vector));
      rpi = SmartPointer<Vector<double> >((Vector<double>*) tree->newDualZVector(empty_vector));
   }

   if (nxlow > 0) {
      rv = SmartPointer<Vector<double> >((Vector<double>*) tree->new_primal_vector());
      rgamma = SmartPointer<Vector<double> >((Vector<double>*) tree->new_primal_vector());
   }
   else {
      rv = SmartPointer<Vector<double> >((Vector<double>*) tree->new_primal_vector(empty_vector));
      rgamma = SmartPointer<Vector<double> >((Vector<double>*) tree->new_primal_vector(empty_vector));
   }

   if (nxupp > 0) {
      rw = SmartPointer<Vector<double> >((Vector<double>*) tree->new_primal_vector());
      rphi = SmartPointer<Vector<double> >((Vector<double>*) tree->new_primal_vector());
   }
   else {
      rw = SmartPointer<Vector<double> >((Vector<double>*) tree->new_primal_vector(empty_vector));
      rphi = SmartPointer<Vector<double> >((Vector<double>*) tree->new_primal_vector(empty_vector));
   }

   createChildren();
}

DistributedResiduals::DistributedResiduals(const DistributedResiduals& res) : Residuals(res) {
   for (unsigned int i = 0; i < res.children.size(); ++i) {
      children.push_back(new DistributedResiduals(*res.children[i]));
   }
}

DistributedResiduals::~DistributedResiduals() {
   for (DistributedResiduals* child : children)
      delete child;
}

void DistributedResiduals::AddChild(DistributedResiduals* child) {
   children.push_back(child);
}

void DistributedResiduals::createChildren() {
   DistributedVector<double>& rQSt = dynamic_cast<DistributedVector<double>&>(*lagrangian_gradient);

   DistributedVector<double>& rASt = dynamic_cast<DistributedVector<double>&>(*rA);
   DistributedVector<double>& rCSt = dynamic_cast<DistributedVector<double>&>(*rC);
   DistributedVector<double>& rzSt = dynamic_cast<DistributedVector<double>&>(*rz);
   DistributedVector<double>& rtSt = dynamic_cast<DistributedVector<double>&>(*rt);
   DistributedVector<double>& rlambdaSt = dynamic_cast<DistributedVector<double>&>(*rlambda);
   DistributedVector<double>& ruSt = dynamic_cast<DistributedVector<double>&>(*ru);
   DistributedVector<double>& rpiSt = dynamic_cast<DistributedVector<double>&>(*rpi);
   DistributedVector<double>& rvSt = dynamic_cast<DistributedVector<double>&>(*rv);
   DistributedVector<double>& rgammaSt = dynamic_cast<DistributedVector<double>&>(*rgamma);
   DistributedVector<double>& rwSt = dynamic_cast<DistributedVector<double>&>(*rw);
   DistributedVector<double>& rphiSt = dynamic_cast<DistributedVector<double>&>(*rphi);
   DistributedVector<double>& ixlowSt = dynamic_cast<DistributedVector<double>&>(*ixlow);
   DistributedVector<double>& ixuppSt = dynamic_cast<DistributedVector<double>&>(*ixupp);
   DistributedVector<double>& iclowSt = dynamic_cast<DistributedVector<double>&>(*iclow);
   DistributedVector<double>& icuppSt = dynamic_cast<DistributedVector<double>&>(*icupp);

   const size_t nChildren = rASt.children.size();
   for (size_t it = 0; it < nChildren; it++) {

      assert(nChildren == rASt.children.size());
      assert(nChildren == rzSt.children.size());
      assert(nChildren == rCSt.children.size());
      assert(nChildren == rtSt.children.size());
      assert(nChildren == rlambdaSt.children.size());
      assert(nChildren == ruSt.children.size());
      assert(nChildren == rpiSt.children.size());
      assert(nChildren == rvSt.children.size());
      assert(nChildren == rgammaSt.children.size());
      assert(nChildren == rwSt.children.size());
      assert(nChildren == rphiSt.children.size());
      assert(nChildren == ixlowSt.children.size());
      assert(nChildren == ixuppSt.children.size());
      assert(nChildren == iclowSt.children.size());
      assert(nChildren == icuppSt.children.size());

      AddChild(new DistributedResiduals(rQSt.children[it], rASt.children[it], rCSt.children[it], rzSt.children[it], rtSt.children[it],
            rlambdaSt.children[it], ruSt.children[it], rpiSt.children[it], rvSt.children[it], rgammaSt.children[it], rwSt.children[it],
            rphiSt.children[it], ixlowSt.children[it], nxlow, ixuppSt.children[it], nxupp, iclowSt.children[it], mclow, icuppSt.children[it], mcupp));
   }
}

void DistributedResiduals::collapseHierarchicalStructure(const DistributedQP& data_hier, const DistributedTree* tree_hier, SmartPointer<Vector<double> > ixlow_,
      SmartPointer<Vector<double> > ixupp_, SmartPointer<Vector<double> > iclow_, SmartPointer<Vector<double> > icupp_) {
   dynamic_cast<DistributedVector<double>&>(*lagrangian_gradient).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);

   const bool empty_vec = true;
   if (nxlow > 0) {
      dynamic_cast<DistributedVector<double>&>(*rv).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);
      dynamic_cast<DistributedVector<double>&>(*rgamma).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);
   }
   else {
      dynamic_cast<DistributedVector<double>&>(*rv).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL, empty_vec);
      dynamic_cast<DistributedVector<double>&>(*rgamma).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL, empty_vec);
   }

   if (nxupp > 0) {
      dynamic_cast<DistributedVector<double>&>(*rw).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);
      dynamic_cast<DistributedVector<double>&>(*rphi).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);
   }
   else {
      dynamic_cast<DistributedVector<double>&>(*rw).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL, empty_vec);
      dynamic_cast<DistributedVector<double>&>(*rphi).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL, empty_vec);
   }

   dynamic_cast<DistributedVector<double>&>(*rA).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Y);

   dynamic_cast<DistributedVector<double>&>(*rC).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*rz).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);

   if (mcupp > 0) {
      dynamic_cast<DistributedVector<double>&>(*ru).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
      dynamic_cast<DistributedVector<double>&>(*rpi).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
   }
   else {
      dynamic_cast<DistributedVector<double>&>(*ru).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z, empty_vec);
      dynamic_cast<DistributedVector<double>&>(*rpi).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z, empty_vec);
   }

   if (mclow > 0) {
      dynamic_cast<DistributedVector<double>&>(*rt).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
      dynamic_cast<DistributedVector<double>&>(*rlambda).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
   }
   else {
      dynamic_cast<DistributedVector<double>&>(*rt).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z, empty_vec);
      dynamic_cast<DistributedVector<double>&>(*rlambda).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z, empty_vec);
   }


   ixlow = ixlow_;
   ixupp = ixupp_;
   iclow = iclow_;
   icupp = icupp_;

   for (size_t c = 0; c < children.size(); c++)
      delete children[c];

   children.clear();
   createChildren();
}

void DistributedResiduals::permuteVec0Entries(const std::vector<unsigned int>& perm, bool resids_only) {
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

void DistributedResiduals::permuteEqLinkingEntries(const std::vector<unsigned int>& perm) {
   dynamic_cast<DistributedVector<double>&>(*rA).permuteLinkingEntries(perm);
}

void DistributedResiduals::permuteIneqLinkingEntries(const std::vector<unsigned int>& perm, bool resids_only) {
   if (!resids_only) {
      dynamic_cast<DistributedVector<double>&>(*iclow).permuteLinkingEntries(perm);
      dynamic_cast<DistributedVector<double>&>(*icupp).permuteLinkingEntries(perm);
   }

   dynamic_cast<DistributedVector<double>&>(*rC).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*rz).permuteLinkingEntries(perm);

   if (mcupp > 0) {
      dynamic_cast<DistributedVector<double>&>(*ru).permuteLinkingEntries(perm);
      dynamic_cast<DistributedVector<double>&>(*rpi).permuteLinkingEntries(perm);
   }

   if (mclow > 0) {
      dynamic_cast<DistributedVector<double>&>(*rt).permuteLinkingEntries(perm);
      dynamic_cast<DistributedVector<double>&>(*rlambda).permuteLinkingEntries(perm);
   }
}

bool DistributedResiduals::isRootNodeInSync() const {
   bool in_sync = true;
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);
   if (!dynamic_cast<const DistributedVector<double>&>(*lagrangian_gradient).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rQ not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*rC).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rC not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*rA).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "rA not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*rz).isRootNodeInSync()) {
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
