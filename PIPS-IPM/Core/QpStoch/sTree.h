/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef STOCH_TREE_BASE
#define STOCH_TREE_BASE

#include "StochResourcesMonitor.h"
#include "StochVector.h"
#include "StochGenMatrix.h"
#include "StochSymMatrix.h"

#include "list"

#include "mpi.h"

class sTreeCallbacks;

class sTree
{
  friend sTreeCallbacks;
 public:
  StochNodeResourcesMonitor resMon;
  static StochIterateResourcesMonitor iterMon;

 protected:

  MPI_Comm commWrkrs{ MPI_COMM_NULL };
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
  std::vector<sTree*> children;
  /* used for hierarchical approach - implies sub structure inside the current Bmat */
  sTree* sub_root{};

  /* global number of all processes available */
  static int numProcs;

  bool is_hierarchical_root = false;
  bool is_hierarchical_inner_root = false;
  bool is_hierarchical_inner_leaf = false;

public:
  // global sizes are still local to each MPI process - they just sum all local data
  virtual void computeGlobalSizes() = 0;
  void getGlobalSizes(long long& n, long long& my, long long& mz) const;
  void getGlobalSizes(long long& n, long long& my, long long& myl, long long& mzlong, long long& mzl) const;

  void assignProcesses( MPI_Comm comm = MPI_COMM_WORLD);

  virtual ~sTree();

  bool distributedPreconditionerActive() const;

  void startMonitors(); void startNodeMonitors();
  void stopMonitors();  void stopNodeMonitors();
  void syncMonitoringData(vector<double>& vCPUTotal);
  bool balanceLoad();
  bool balanceLoadPrecond();

  void getSyncInfo(int myRank, int& syncNeeded, int& sendOrRecv, int& toFromCPU );

  virtual StochSymMatrix* createQ() const = 0;
  virtual StochVector* createc() const = 0;

  virtual StochVector* createxlow()  const = 0;
  virtual StochVector* createixlow() const = 0;
  virtual StochVector* createxupp()  const = 0;
  virtual StochVector* createixupp() const = 0;

  virtual StochGenMatrix* createA() const = 0;
  virtual StochVector* createb() const = 0;


  virtual StochGenMatrix* createC() const = 0;
  virtual StochVector* createclow()  const = 0;
  virtual StochVector* createiclow() const = 0;
  virtual StochVector* createcupp()  const = 0;
  virtual StochVector* createicupp() const = 0;

  StochVector* newPrimalVector(bool empty = false) const;
  StochVector* newDualYVector(bool empty = false)  const;
  StochVector* newDualZVector(bool empty = false)  const;

  StochVector* newRhs();

  const sTree* getSubRoot() const { return sub_root; };
  const std::vector<sTree*>& getChildren() const { return children; };
  int nChildren() const { return children.size(); }
  MPI_Comm getCommWorkers() const { return commWrkrs; };

  int innerSize(int which) const;
  virtual int nx() const = 0;
  virtual int my() const = 0;
  virtual int myl() const;
  virtual int mzl() const;
  virtual int mz() const = 0; 
  virtual int id() const = 0;

  //returns the global load, i.e. statistic based on the NNZs and
  //dimensions of the node (and subnodes) subproblem before any iteration or the CPU
  //time of this node and its subnodes  after the first iteration.
  double processLoad() const;

  bool isHierarchicalRoot() const { return is_hierarchical_root; };

  void setHierarchicalInnerRoot() { is_hierarchical_inner_root = true; };
  bool isHierarchicalInnerRoot() const { return is_hierarchical_inner_root; };

  void setHierarchicalInnerLeaf() { is_hierarchical_inner_leaf = true; };
  bool isHierarchicalInnerLeaf() const { return is_hierarchical_inner_leaf; };

  /* shave tree and add an additional top layer */
  virtual sTree* shaveDenseBorder( int nx_to_shave, int myl_to_shave, int mzl_to_shave) = 0;
  /* add an additional layer below this one by adding sqrt(nChildren) children each with sqrt(nChildren) of our current children */
  virtual void splitTreeSquareRoot( const std::vector<int>& twoLinksStartBlockA, const std::vector<int>& twoLinksStartBlockC, bool silent = false ) = 0;

  // TODO : make sure that none of the not suitable methods get called...
  virtual sTree* switchToHierarchicalTree( int nx_to_shave, int myl_to_shave, int mzl_to_shave, const std::vector<int>& twoLinksStartBlockA,
        const std::vector<int>& twoLinksStartBlockC, bool silent = false ) = 0;
  virtual sTree * collapseHierarchicalTree() = 0;

protected:
  void assignProcesses( MPI_Comm, vector<int>&);

  sTree() = default;

  void toMonitorsList( std::list<NodeExecEntry>& );
  void fromMonitorsList( std::list<NodeExecEntry>& );

  void computeNodeTotal();

  void saveCurrentCPUState();

  int isInVector(int elem, const vector<int>& vec) const;

};

#endif 
