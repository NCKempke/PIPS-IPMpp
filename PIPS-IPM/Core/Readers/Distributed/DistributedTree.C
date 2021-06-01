/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#include "DistributedTree.h"
#include "DistributedQP.hpp"
#include "DistributedVector.h"
#include "SimpleVector.h"
#include <numeric>
#include <algorithm>
#include <cmath>

Timer DistributedTree::iterMon;
int DistributedTree::rankMe = -1;
int DistributedTree::rankZeroW = 0;
int DistributedTree::rankPrcnd = -1;
int DistributedTree::numProcs = -1;

DistributedTree::DistributedTree(const DistributedTree& other) : commWrkrs{other.commWrkrs}, myProcs(other.myProcs.begin(), other.myProcs.end()),
      myOldProcs(other.myOldProcs.begin(), other.myOldProcs.end()), commP2ZeroW{other.commP2ZeroW}, N{other.N}, MY{other.MY}, MZ{other.MZ},
      MYL{other.MYL}, MZL{other.MZL}, np{other.np}, IPMIterExecTIME{other.IPMIterExecTIME}, is_hierarchical_root{other.is_hierarchical_root},
      is_hierarchical_inner_root{other.is_hierarchical_inner_root}, is_hierarchical_inner_leaf{other.is_hierarchical_inner_leaf} {
   if (other.sub_root)
      sub_root = other.sub_root->clone();
   for (auto& child : other.children)
      children.push_back(child->clone());
}

DistributedTree::~DistributedTree() {
   for (size_t it = 0; it < children.size(); it++)
      delete children[it];
   delete sub_root;
}

int DistributedTree::myl() const {
   return -1;
}

int DistributedTree::mzl() const {
   return -1;
}

bool DistributedTree::distributedPreconditionerActive() const {
   return (rankZeroW != 0) && (rankPrcnd != -1) && (commP2ZeroW != MPI_COMM_NULL) && (rankMe != -1);
}

void DistributedTree::assignProcesses(MPI_Comm comm) {
   assert(comm != MPI_COMM_NULL);
   assert(!is_hierarchical_root);

   const int world_size = PIPS_MPIgetSize(comm);
   /* assign root node */
   commWrkrs = MPI_COMM_WORLD;

   std::vector<int> processes(world_size);
   for (int p = 0; p < world_size; p++)
      processes[p] = p;
   myProcs = processes;

#if 0
   /* inactive but used to calculate loads that children would generate and then assign the processes */
   vector<double> child_nodes_load( children.size() );
   for(size_t i = 0; i < children.size(); i++)
      child_nodes_load[i] = children[i]->processLoad();
#endif

   const unsigned int n_procs = processes.size();
   //**** solve the assignment problem ****
   /* too many MPI processes? */
   std::stringstream ss;
   ss << "too many MPI processes! (max: " << children.size() << ")\n";
   PIPS_MPIabortIf(n_procs > children.size(), ss.str());

   std::vector<unsigned int> map_child_nodes_to_procs;
   mapChildrenToNSubTrees(map_child_nodes_to_procs, children.size(), n_procs);

#ifndef NDEBUG
   for (size_t i = 0; i < children.size(); i++)
      assert(map_child_nodes_to_procs[i] < n_procs);

   for (size_t i = 1; i < children.size(); i++)
      assert(map_child_nodes_to_procs[i - 1] <= map_child_nodes_to_procs[i]);

   std::vector<size_t> load_per_proc(n_procs, 0);
   for (size_t i = 0; i < children.size(); ++i)
      load_per_proc[map_child_nodes_to_procs[i]]++;

   const size_t max_load = *std::max_element(load_per_proc.begin(), load_per_proc.end());
   const size_t min_load = *std::min_element(load_per_proc.begin(), load_per_proc.end());
   assert(max_load == min_load || max_load == min_load + 1);
#endif

   for (size_t i = 0; i < children.size(); i++) {
      children[i]->myProcs.clear();
      children[i]->myProcs.push_back(map_child_nodes_to_procs[i]);
      assert(rankMe >= 0);
      if (static_cast<unsigned int>(rankMe) == map_child_nodes_to_procs[i])
         children[i]->commWrkrs = MPI_COMM_SELF;
      else
         children[i]->commWrkrs = MPI_COMM_NULL;
   }
}

double DistributedTree::processLoad() const {
   assert(!is_hierarchical_root);
   //! need a recursive and also a collective call
   if (IPMIterExecTIME < 0.0)
      //return (NNZQ+NNZA+NNZB+NNZC+NNZD + N+MY+MZ)/1000.0;
      return (N + MY + MZ) / 1000.0;
   return IPMIterExecTIME;
}

void DistributedTree::getGlobalSizes(long long& n, long long& my, long long& mz) const {
   n = N;
   my = MY;
   mz = MZ;
}

void DistributedTree::getGlobalSizes(long long& n, long long& my, long long& myl, long long& mz, long long& mzl) const {
   n = N;
   my = MY;
   mz = MZ;
   myl = MYL, mzl = MZL;
}

DistributedVector<double>* DistributedTree::new_primal_vector(bool empty) const {
   if (commWrkrs == MPI_COMM_NULL) {
      return new DistributedDummyVector<double>();
   }

   if (!sub_root) {
      auto* x = new DistributedVector<double>(empty ? 0 : nx(), commWrkrs);
      for (auto it : children) {
         std::shared_ptr<DistributedVector<double>> child{it->new_primal_vector(empty)};
         x->AddChild(std::move(child));
      }
      return x;
   }
   else {
      assert(children.empty());
      if (sub_root->commWrkrs == MPI_COMM_NULL) {
         DistributedVector<double>* x = new DistributedDummyVector<double>();
         return x;
      }
      else {
         std::unique_ptr<DistributedVector<double>> x_vec{sub_root->new_primal_vector(empty)};
         auto* x = new DistributedVector<double>(std::move(x_vec), nullptr, commWrkrs);
         dynamic_cast<DistributedVector<double>&>(*x->first).parent = x;
         return x;
      }
   }
}

DistributedVector<double>* DistributedTree::newDualYVector(bool empty) const {
   if (commWrkrs == MPI_COMM_NULL)
      return new DistributedDummyVector<double>();

   DistributedVector<double>* y{};
   const int yl = (np == -1) ? myl() : -1;

   if (!sub_root) {
         y = new DistributedVector<double>(empty ? std::min(0, my()) : my(), empty ? std::min(yl, 0) : yl, commWrkrs);

      for (auto it : children) {
         std::shared_ptr<DistributedVector<double>> child{it->newDualYVector(empty)};
         y->AddChild(std::move(child));
      }
   }
   else {
      assert(children.empty());

      if (sub_root->commWrkrs == MPI_COMM_NULL)
         y = new DistributedDummyVector<double>();
      else {
         std::unique_ptr<DistributedVector<double>> y_vec{sub_root->newDualYVector(empty)};
         assert(yl == -1);

         y = new DistributedVector<double>(std::move(y_vec), nullptr, commWrkrs);
         dynamic_cast<DistributedVector<double>&>(*y->first).parent = y;
      }
   }
   return y;
}

DistributedVector<double>* DistributedTree::newDualZVector(bool empty) const {
   if (commWrkrs == MPI_COMM_NULL)
      return new DistributedDummyVector<double>();

   DistributedVector<double>* z{};
   const int zl = (np == -1) ? mzl() : -1;

   if (!sub_root) {

      z = new DistributedVector<double>(empty ? std::min(mz(), 0) : mz(), empty ? std::min(zl, 0) : zl, commWrkrs);
      for (auto it : children) {
         std::shared_ptr<DistributedVector<double>> child{it->newDualZVector(empty)};
         z->AddChild(std::move(child));
      }
   }
   else {
      assert(children.empty());
      if (sub_root->commWrkrs == MPI_COMM_NULL)
         z = new DistributedDummyVector<double>();
      else {
         std::unique_ptr<DistributedVector<double>> z_vec{sub_root->newDualZVector(empty)};
         assert(zl == -1);

         z = new DistributedVector<double>(std::move(z_vec), nullptr, commWrkrs);
         dynamic_cast<DistributedVector<double>&>(*z->first).parent = z;
      }
   }

   return z;
}

DistributedVector<double>* DistributedTree::newRhs() const {
   //is this node a dead-end for this process?
   if (commWrkrs == MPI_COMM_NULL)
      return new DistributedDummyVector<double>();

   DistributedVector<double>* rhs{};
   if (!sub_root) {
      int locmyl = (np == -1) ? myl() : 0;
      int locmzl = (np == -1) ? mzl() : 0;

      locmyl = std::max(locmyl, 0);
      locmzl = std::max(locmzl, 0);

      rhs = new DistributedVector<double>(nx() + std::max(my(), 0) + std::max(mz(), 0) + locmyl + locmzl, commWrkrs);

      for (auto it : children) {
         std::shared_ptr<DistributedVector<double>> child{it->newRhs()};
         rhs->AddChild(std::move(child));
      }
   }
   else
      rhs = sub_root->newRhs();

   return rhs;
}

void DistributedTree::startMonitors() {
   iterMon.start();
   startNodeMonitors();
}

void DistributedTree::stopMonitors() {
   iterMon.stop();
   stopNodeMonitors();
}

void DistributedTree::startNodeMonitors() {
   resMon.reset();
   for (auto & i : children)
      i->startNodeMonitors();
}

void DistributedTree::stopNodeMonitors() {
   for (auto & i : children)
      i->stopNodeMonitors();

   resMon.computeTotal();
}

void DistributedTree::toMonitorsList(std::list<NodeTimer>& lstExecTm) {
   lstExecTm.push_back(resMon.eTotal);

   for (auto & i : children)
      i->toMonitorsList(lstExecTm);
}

void DistributedTree::fromMonitorsList(std::list<NodeTimer>& lstExecTm) {
   resMon.eTotal = lstExecTm.front();
   lstExecTm.pop_front();

   for (auto & i : children)
      i->fromMonitorsList(lstExecTm);
}

bool DistributedTree::balanceLoad() {
   return false; //disabled for now
}

#define maSend 1
#define maRecv 0

void DistributedTree::getSyncInfo(int rank, int& syncNeeded, int& sendOrRecv, int& toFromCPU) {
   // was this node previously assigned to cpu 'rank'?
   int wasAssigned = isInVector(rank, myOldProcs);
   // is currently assigned to cpu 'rank'?
   int isAssigned = isInVector(rank, myProcs);

   syncNeeded = 0;
   sendOrRecv = 0;
   toFromCPU = -1;
   if (isAssigned != wasAssigned) {
      syncNeeded = 1;

      if (wasAssigned) {
         assert(0 == isAssigned);
         assert(!myProcs.empty());
         //where is this node assigned?
         toFromCPU = myProcs[0];
         sendOrRecv = maSend;
      }
      else {
         assert(1 == isAssigned);
         assert(!myOldProcs.empty());
         toFromCPU = myOldProcs[0];
         sendOrRecv = maRecv;
      }
   }
}

void DistributedTree::computeNodeTotal() {
   if (children.empty())
      this->IPMIterExecTIME = resMon.eTotal.children_time;
   else {

      this->IPMIterExecTIME = resMon.eTotal.local_time;
      for (auto & i : children) {
         i->computeNodeTotal();

         this->IPMIterExecTIME += i->IPMIterExecTIME;
      }
   }
}

void DistributedTree::saveCurrentCPUState() {
   myOldProcs = myProcs;

   for (auto & i : children)
      i->saveCurrentCPUState();
}


void DistributedTree::appendPrintTreeLayer(std::vector<std::string>& layer_outputs, unsigned int level) const {
   if (level == layer_outputs.size())
      layer_outputs.emplace_back("");
   if (commWrkrs == MPI_COMM_NULL)
      return;

   std::stringstream curr_level_output;
   const bool I_print_layer = PIPS_MPIgetRank(commWrkrs) == 0;

   if (I_print_layer) {
      if (!children.empty()) {

         std::stringstream level_stream;
         level_stream << "|";
         for (size_t i = 0; i < children.size(); ++i) {
            const auto& child = children[i];

            /* this is a leaf */
            if (child->children.empty() && !child->sub_root) {
               int count = 1;
               while (i + 1 < children.size() && children[i + 1]->myProcs == child->myProcs) {
                  ++count;
                  ++i;
               }

               level_stream << " [ " << child->myProcs.front() << " : " << count << " leafs ]";
            }
            else {
               if (child->myProcs.size() == 1)
                  level_stream << " [ " << child->myProcs.front() << " ]";
               else
                  level_stream << " [ " << child->myProcs.front() << " - " << child->myProcs.back() << " ]";
            }
         }
         level_stream << " |";

         layer_outputs[level].append(level_stream.str());
      }
   }

   for (const auto& child : children) {
      if (child->sub_root)
         child->sub_root->appendPrintTreeLayer(layer_outputs, level + 1);
      else
         child->appendPrintTreeLayer(layer_outputs, level + 1);
   }
}

void DistributedTree::printProcessTree() const {
   assert(commWrkrs != MPI_COMM_NULL);

   unsigned int level = 0;
   bool tree_unbalanced = false;

   if (rankMe == 0)
      std::cout << "Process Tree:\n\n";

   std::vector<std::string> level_outputs;
   if (PIPS_MPIgetRank(commWrkrs) == 0) {
      std::stringstream level_stream;
      if (myProcs.size() == 1)
         level_stream << "| [ " << myProcs.front() << " ] |";
      else
         level_stream << "| [ " << myProcs.front() << " - " << myProcs.back() << " ] |";

      level_outputs.push_back(level_stream.str());
   }
   else
      level_outputs.emplace_back("");

   appendPrintTreeLayer(level_outputs, level + 1);

   const unsigned int levels = PIPS_MPIgetMax(level_outputs.size());

   if (levels > level_outputs.size()) {
      assert(levels == level_outputs.size() + 1);
      tree_unbalanced = true;
      level_outputs.emplace_back("");
   }
   assert(PIPS_MPIisValueEqual(level_outputs.size()));

   for (auto& level : level_outputs) {
      const std::string level_full = PIPS_MPIallgatherString(level);

      if (rankMe == 0)
         std::cout << level_full << "\n\n";
   }

   PIPS_MPIgetLogicOrInPlace(tree_unbalanced);
   if (tree_unbalanced && rankMe == 0)
      std::cout << "WARNING : the split tree is unbalanced and there is branches that are longer than others\n"
                << "\t-> this is bad for performance but that fine of a split usually is bad for performance anyway..\n"
                << "\t-> rethink the amount of layers for the hierarchical approach..\n\n";
}

void DistributedTree::mapChildrenToNSubTrees(std::vector<unsigned int>& map_child_to_sub_tree, unsigned int n_children, unsigned int n_subtrees) {
   assert(n_subtrees <= n_children);
   map_child_to_sub_tree.clear();
   map_child_to_sub_tree.reserve(n_children);

   if (n_subtrees == 0)
      return;

   const unsigned int everyone_gets = std::floor(n_children / n_subtrees);
   const unsigned int n_leftovers = n_children % n_subtrees;

   std::vector<unsigned int> children_per_tree(n_subtrees, everyone_gets);

   if (n_leftovers > 0) {
      const unsigned int free_in_leftover_row = n_subtrees - n_leftovers;

      const unsigned int free_after_each_leftover = std::floor(free_in_leftover_row / n_leftovers);
      const unsigned int additional_frees = free_in_leftover_row % n_leftovers;

      unsigned int additional_frees_assigned = 0;
      for (unsigned int i = 0; i < n_subtrees; i += free_after_each_leftover) {
         ++children_per_tree[i];
         ++i;

         if (additional_frees != additional_frees_assigned) {
            ++i;
            ++additional_frees_assigned;
         }
      }
      assert(additional_frees_assigned == additional_frees);
   }

   assert(std::accumulate(children_per_tree.begin(), children_per_tree.end(), decltype(children_per_tree)::value_type(0)) == n_children);

   for (unsigned int i = 0; i < children_per_tree.size(); ++i) {
      for (unsigned int j = 0; j < children_per_tree[i]; ++j)
         map_child_to_sub_tree.push_back(i);
   }
}
