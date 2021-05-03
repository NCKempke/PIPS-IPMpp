#include "DistributedVariables.h"
#include "Vector.hpp"
#include "SmartPointer.h"
#include "Problem.h"
#include "DistributedVector.h"
#include "DistributedTree.h"
#include "DistributedTreeCallbacks.h"
#include "LinearAlgebraPackage.h"
#include <iostream>

DistributedVariables::DistributedVariables(const DistributedTree* tree, Vector<double>* x_in, Vector<double>* s_in, Vector<double>* y_in, Vector<double>* z_in,
      Vector<double>* v_in, Vector<double>* gamma_in, Vector<double>* w_in, Vector<double>* phi_in, Vector<double>* t_in, Vector<double>* lambda_in,
      Vector<double>* u_in, Vector<double>* pi_in, Vector<double>* ixlow_in, long long nxlowGlobal, Vector<double>* ixupp_in, long long nxuppGlobal,
      Vector<double>* iclow_in, long long mclowGlobal, Vector<double>* icupp_in, long long mcuppGlobal) : Variables() {

   stochNode = tree;

   SpReferTo(x, x_in);
   SpReferTo(s, s_in);
   SpReferTo(y, y_in);
   SpReferTo(z, z_in);
   SpReferTo(v, v_in);
   SpReferTo(phi, phi_in);
   SpReferTo(w, w_in);
   SpReferTo(gamma, gamma_in);
   SpReferTo(t, t_in);
   SpReferTo(lambda, lambda_in);
   SpReferTo(u, u_in);
   SpReferTo(pi, pi_in);
   SpReferTo(ixlow, ixlow_in);
   SpReferTo(ixupp, ixupp_in);
   SpReferTo(iclow, iclow_in);
   SpReferTo(icupp, icupp_in);

   nx = x->length();
   my = y->length();
   mz = z->length();

   assert(nx == ixlow->length() || 0 == ixlow->length());
   assert(nx == ixupp->length() || 0 == ixupp->length());
   assert(mz == iclow->length() || 0 == iclow->length());
   assert(mz == icupp->length() || 0 == icupp->length());

   nxlow = nxlowGlobal;
   nxupp = nxuppGlobal;
   mclow = mclowGlobal;
   mcupp = mcuppGlobal;
   nComplementaryVariables = mclow + mcupp + nxlow + nxupp;

   assert(mz == s->length());
   assert(nx == v->length() || (0 == v->length() && nxlow == 0));
   assert(nx == gamma->length() || (0 == gamma->length() && nxlow == 0));

   assert(nx == w->length() || (0 == w->length() && nxupp == 0));
   assert(nx == phi->length() || (0 == phi->length() && nxupp == 0));

   assert(mz == t->length() || (0 == t->length() && mclow == 0));
   assert(mz == lambda->length() || (0 == lambda->length() && mclow == 0));

   assert(mz == u->length() || (0 == u->length() && mcupp == 0));
   assert(mz == pi->length() || (0 == pi->length() && mcupp == 0));

   createChildren();
}

DistributedVariables::DistributedVariables(const DistributedVariables& vars) : Variables(vars) {
   stochNode = vars.stochNode;
   for (unsigned int i = 0; i < vars.children.size(); ++i) {
      children.push_back(new DistributedVariables(*vars.children[i]));
   }
}

DistributedVariables::~DistributedVariables() {
   for (size_t c = 0; c < children.size(); c++)
      delete children[c];
}

void DistributedVariables::AddChild(DistributedVariables* child) {
   children.push_back(child);
}


void DistributedVariables::createChildren() {
   DistributedVector<double>& xst = dynamic_cast<DistributedVector<double>&>(*x);
   DistributedVector<double>& sst = dynamic_cast<DistributedVector<double>&>(*s);
   DistributedVector<double>& yst = dynamic_cast<DistributedVector<double>&>(*y);
   DistributedVector<double>& zst = dynamic_cast<DistributedVector<double>&>(*z);
   DistributedVector<double>& vst = dynamic_cast<DistributedVector<double>&>(*v);
   DistributedVector<double>& gammast = dynamic_cast<DistributedVector<double>&>(*gamma);
   DistributedVector<double>& wst = dynamic_cast<DistributedVector<double>&>(*w);
   DistributedVector<double>& phist = dynamic_cast<DistributedVector<double>&>(*phi);
   DistributedVector<double>& tst = dynamic_cast<DistributedVector<double>&>(*t);
   DistributedVector<double>& lambdast = dynamic_cast<DistributedVector<double>&>(*lambda);
   DistributedVector<double>& ust = dynamic_cast<DistributedVector<double>&>(*u);
   DistributedVector<double>& pist = dynamic_cast<DistributedVector<double>&>(*pi);
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
   dynamic_cast<DistributedVector<double>&>(*x).collapseFromHierarchical(hier_data, *stochNode, VectorType::PRIMAL);

   dynamic_cast<DistributedVector<double>&>(*v).collapseFromHierarchical(hier_data, *stochNode, VectorType::PRIMAL);
   dynamic_cast<DistributedVector<double>&>(*w).collapseFromHierarchical(hier_data, *stochNode, VectorType::PRIMAL);
   dynamic_cast<DistributedVector<double>&>(*phi).collapseFromHierarchical(hier_data, *stochNode, VectorType::PRIMAL);
   dynamic_cast<DistributedVector<double>&>(*gamma).collapseFromHierarchical(hier_data, *stochNode, VectorType::PRIMAL);

   dynamic_cast<DistributedVector<double>&>(*y).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Y);

   dynamic_cast<DistributedVector<double>&>(*s).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*z).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*t).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*u).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*pi).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);
   dynamic_cast<DistributedVector<double>&>(*lambda).collapseFromHierarchical(hier_data, *stochNode, VectorType::DUAL_Z);

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

   dynamic_cast<DistributedVector<double>&>(*x).permuteVec0Entries(perm);
   dynamic_cast<DistributedVector<double>&>(*v).permuteVec0Entries(perm);
   dynamic_cast<DistributedVector<double>&>(*w).permuteVec0Entries(perm);
   dynamic_cast<DistributedVector<double>&>(*phi).permuteVec0Entries(perm);
   dynamic_cast<DistributedVector<double>&>(*gamma).permuteVec0Entries(perm);
}

void DistributedVariables::permuteEqLinkingEntries(const std::vector<unsigned int>& perm) {
   dynamic_cast<DistributedVector<double>&>(*y).permuteLinkingEntries(perm);
}

void DistributedVariables::permuteIneqLinkingEntries(const std::vector<unsigned int>& perm, bool vars_only) {
   if (!vars_only) {
      dynamic_cast<DistributedVector<double>&>(*iclow).permuteLinkingEntries(perm);
      dynamic_cast<DistributedVector<double>&>(*icupp).permuteLinkingEntries(perm);
   }

   dynamic_cast<DistributedVector<double>&>(*s).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*z).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*t).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*u).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*pi).permuteLinkingEntries(perm);
   dynamic_cast<DistributedVector<double>&>(*lambda).permuteLinkingEntries(perm);
}

bool DistributedVariables::isRootNodeInSync() const {
   bool in_sync = true;
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);
   if (!dynamic_cast<const DistributedVector<double>&>(*x).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "x not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*s).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "s not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*y).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "y not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*z).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "z not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*v).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "v not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*gamma).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "gamma not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*w).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "w not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*phi).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "phi not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*t).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "t not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*lambda).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "lambda not in sync" << std::endl;
      in_sync = false;
   }

   if (!dynamic_cast<const DistributedVector<double>&>(*u).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "u not in sync" << std::endl;
      in_sync = false;
   }
   if (!dynamic_cast<const DistributedVector<double>&>(*pi).isRootNodeInSync()) {
      if (my_rank == 0)
         std::cout << "pi not in sync" << std::endl;
      in_sync = false;
   }

   MPI_Barrier(MPI_COMM_WORLD);
   return in_sync;
}
