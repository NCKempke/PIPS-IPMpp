/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef STOCH_TREE_BASE
#define STOCH_TREE_BASE

#include "StochResourcesMonitor.h"
#include "DistributedVector.h"
#include "DistributedMatrix.h"
#include "DistributedSymmetricMatrix.h"

#include <list>
#include <utility>
#include "mpi.h"

class DistributedTreeCallbacks;

class DistributedQP;

class DistributedTree {
   friend DistributedTreeCallbacks;
public:
   StochNodeResourcesMonitor resMon;
   static Timer iterMon;

   long long getN() const { return N; };
   long long getMY() const { return MY; };
   long long getMYL() const { return MYL; };
   long long getMZ() const { return MZ; };
   long long getMZL() const { return MZL; };

   virtual DistributedTree* clone() const = 0;
protected:

   DistributedTree(const DistributedTree& other);

   MPI_Comm commWrkrs{MPI_COMM_NULL};
   std::vector<int> myProcs, myOldProcs;

   MPI_Comm commP2ZeroW{MPI_COMM_NULL};   // preconditioner (rank P+1) and special (rank 0) worker
   static int rankPrcnd;   // rank of preconditioner
   static int rankZeroW;   // rank of the "special" worker (the root, or 0-rank process)
   static int rankMe;      // rank of the running process

   /* global sizes - global meaning on this process - so the sum of all local matrices - MY, MZ do not include linking constraints */
   long long N{0};
   long long MY{0};
   long long MZ{0};
   long long MYL{0};
   long long MZL{0};
   int np{-1}; //n for the parent

   double IPMIterExecTIME{-1.0}; // not used since we currently do not compute loads for nodes and processes...
   std::vector<DistributedTree*> children;
   /* used for hierarchical approach - implies sub structure inside the current Bmat */
   DistributedTree* sub_root{};

   /* global number of all processes available */
   static int numProcs;

   bool is_hierarchical_root{false};
   bool is_hierarchical_inner_root{false};
   bool is_hierarchical_inner_leaf{false};

public:
   // global sizes are still local to each MPI process - they just sum all local data
   virtual void computeGlobalSizes() = 0;
   void getGlobalSizes(long long& n, long long& my, long long& mz) const;
   void getGlobalSizes(long long& n, long long& my, long long& myl, long long& mzlong, long long& mzl) const;

   void assignProcesses(MPI_Comm comm = MPI_COMM_WORLD);

   virtual ~DistributedTree();

   [[nodiscard]] bool distributedPreconditionerActive() const;

   void startMonitors();
   void startNodeMonitors();
   void stopMonitors();
   void stopNodeMonitors();
   bool balanceLoad();

   void getSyncInfo(int myRank, int& syncNeeded, int& sendOrRecv, int& toFromCPU);

   [[nodiscard]] virtual DistributedSymmetricMatrix* createQ() const = 0;
   [[nodiscard]] virtual DistributedVector<double>* createc() const = 0;

   [[nodiscard]] virtual DistributedVector<double>* createxlow() const = 0;
   [[nodiscard]] virtual DistributedVector<double>* createixlow() const = 0;
   [[nodiscard]] virtual DistributedVector<double>* createxupp() const = 0;
   [[nodiscard]] virtual DistributedVector<double>* createixupp() const = 0;

   [[nodiscard]] virtual DistributedMatrix* createA() const = 0;
   [[nodiscard]] virtual DistributedVector<double>* createb() const = 0;

   [[nodiscard]] virtual DistributedMatrix* createC() const = 0;
   [[nodiscard]] virtual DistributedVector<double>* createclow() const = 0;
   [[nodiscard]] virtual DistributedVector<double>* createiclow() const = 0;
   [[nodiscard]] virtual DistributedVector<double>* createcupp() const = 0;
   [[nodiscard]] virtual DistributedVector<double>* createicupp() const = 0;

   [[nodiscard]] DistributedVector<double>* new_primal_vector(bool empty = false) const;
   [[nodiscard]] DistributedVector<double>* newDualYVector(bool empty = false) const;
   [[nodiscard]] DistributedVector<double>* newDualZVector(bool empty = false) const;

   [[nodiscard]] DistributedVector<double>* newRhs() const;

   [[nodiscard]] const DistributedTree* getSubRoot() const { return sub_root; };
   [[nodiscard]] const std::vector<DistributedTree*>& getChildren() const { return children; };
   [[nodiscard]] unsigned int nChildren() const { return children.size(); }
   [[nodiscard]] MPI_Comm getCommWorkers() const { return commWrkrs; };

   [[nodiscard]] virtual int nx() const = 0;
   [[nodiscard]] virtual int my() const = 0;
   [[nodiscard]] virtual int myl() const;
   [[nodiscard]] virtual int mzl() const;
   [[nodiscard]] virtual int mz() const = 0;
   [[nodiscard]] virtual int id() const = 0;

   //returns the global load, i.e. statistic based on the NNZs and
   //dimensions of the node (and subnodes) subproblem before any iteration or the CPU
   //time of this node and its subnodes  after the first iteration.
   [[nodiscard]] double processLoad() const;

   [[nodiscard]] bool isHierarchicalRoot() const { return is_hierarchical_root; };

   void setHierarchicalInnerRoot() { is_hierarchical_inner_root = true; };
   [[nodiscard]] bool isHierarchicalInnerRoot() const { return is_hierarchical_inner_root; };

   void setHierarchicalInnerLeaf() { is_hierarchical_inner_leaf = true; };
   [[nodiscard]] bool isHierarchicalInnerLeaf() const { return is_hierarchical_inner_leaf; };

   /* shave tree and add an additional top layer */
   virtual DistributedTree* shaveDenseBorder(int nx_to_shave, int myl_to_shave, int mzl_to_shave) = 0;
   /* add an additional layer below this one by adding sqrt(nChildren) children each with sqrt(nChildren) of our current children */
   virtual std::pair<int, int> splitTree(int n_layers, DistributedQP* data) = 0;

   // TODO : make sure that none of the not suitable methods get called...
   virtual DistributedTree* switchToHierarchicalTree(DistributedQP*& data) = 0;

   void printProcessTree() const;
protected:
   void appendPrintTreeLayer(std::vector<std::string>& layer_outputs, unsigned int level) const;

   DistributedTree() = default;

   void toMonitorsList(std::list<NodeTimer>&);
   void fromMonitorsList(std::list<NodeTimer>&);

   void computeNodeTotal();

   void saveCurrentCPUState();

   static void mapChildrenToNSubTrees(std::vector<unsigned int>& map_child_to_sub_tree, unsigned int n_children, unsigned int n_subtrees);
};

#endif 
