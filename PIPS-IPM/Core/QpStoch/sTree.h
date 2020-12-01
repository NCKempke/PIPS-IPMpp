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

class sTree
{
 public:
  virtual ~sTree();

  int NumberOfChildren() const { return children.size(); }

  // global sizes are still local to each MPI process - they just sum all local data
  virtual void computeGlobalSizes() = 0;
  void getGlobalSizes(long long& n, long long& my, long long& mz);

  void assignProcesses( MPI_Comm comm = MPI_COMM_WORLD);

 protected:
  void assignProcesses( MPI_Comm, vector<int>&);

 public:
  MPI_Comm commWrkrs, myOldMpiComm; //workers only
  std::vector<int> myProcs, myOldProcs;

  MPI_Comm commP2ZeroW;   // preconditioner (rank P+1) and special (rank 0) worker
  static int rankPrcnd;   // rank of preconditioner
  static int rankZeroW;   // rank of the "special" worker (the root, or 0-rank process)
  static int rankMe;      // rank of the running process 

  void startMonitors(); void startNodeMonitors();
  void stopMonitors();  void stopNodeMonitors();
  void syncMonitoringData(vector<double>& vCPUTotal);
  bool balanceLoad();
  bool balanceLoadPrecond();

  void getSyncInfo(int myRank, int& syncNeeded, int& sendOrRecv, int& toFromCPU );
  void syncPrimalVector(StochVector& vec) const;
  void syncDualYVector(StochVector& vec) const;
  void syncDualZVector(StochVector& vec) const;
  void syncStochVector(StochVector& vec) const;

  void syncStochGenMatrix(StochGenMatrix& mat) const;
  void syncStochSymMatrix(StochSymMatrix& mat) const;

  virtual StochSymMatrix*   createQ() const = 0;
  virtual StochVector*      createc() const = 0;

  virtual StochVector*      createxlow()  const = 0;
  virtual StochVector*      createixlow() const = 0;
  virtual StochVector*      createxupp()  const = 0;
  virtual StochVector*      createixupp() const = 0;


  virtual StochGenMatrix*   createA() const = 0;
  virtual StochVector*      createb() const = 0;


  virtual StochGenMatrix*   createC() const = 0;
  virtual StochVector*      createclow()  const = 0;
  virtual StochVector*      createiclow() const = 0;
  virtual StochVector*      createcupp()  const = 0;
  virtual StochVector*      createicupp() const = 0;

  StochVector*      newPrimalVector() const;
  StochVector*      newDualYVector()  const;
  StochVector*      newDualZVector()  const;
  StochVector*      newPrimalVectorEmpty() const;
  StochVector*      newDualYVectorEmpty()  const;
  StochVector*      newDualZVectorEmpty()  const;

  StochVector*      newRhs();

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

 protected:
  sTree();

  void toMonitorsList( std::list<NodeExecEntry>& );
  void fromMonitorsList( std::list<NodeExecEntry>& );

  void computeNodeTotal();

  void saveCurrentCPUState();

  int isInVector(int elem, const vector<int>& vec);

  bool is_hierarchical_root = false;
  bool is_hierarchical_inner = false;
  bool is_hierarchical_leaf = false;

 public:
  /* global sizes - global meaning on this process - so the sum of all local matrices - MY, MZ do not include linking constraints */
  long long N, MY, MZ, MYL, MZL;//global sizes
  int np; //n for the parent

  double IPMIterExecTIME; // not used since we currently do not compute loads for nodes and processes...
  std::vector<sTree*> children;

  /* global number of all processes available */
  static int numProcs;

  StochNodeResourcesMonitor resMon;
  static StochIterateResourcesMonitor iterMon;
#ifdef STOCH_TESTING
  void displayProcessInfo(int onWhichRank=0);
  void displayProcessInfo(char* tab);
  void runTestNGP();
  void displayExecTimes(int onWhichRank=0);
  void displayExecTimes(char* szTabbing);

  void displayVectorVsTreeStructure(StochVector& stVec, int rank, char* szTab);
  void displayVectorVsTreeStructure(StochVector& stVec, int rank);

  void displayMatVsTreeStructure(StochGenMatrix& stVec, int myRank, char* tab);
  void displayMatVsTreeStructure(StochGenMatrix& stVec, int myRank);

  void displayMatVsTreeStructure(StochSymMatrix& stVec, int myRank, char* tab);
  void displayMatVsTreeStructure(StochSymMatrix& stVec, int myRank);
#endif
//to be called after assignProcesses
  virtual void loadLocalSizes() = 0;

  bool isHierarchicalRoot() const { return is_hierarchical_root; };
  bool isHierarchicalInner() const { return is_hierarchical_inner; };

  /* shave tree and add an additional top layer */
  virtual sTree* shaveDenseBorder( int nx_to_shave, int myl_to_shave, int mzl_to_shave) = 0;
  /* add an additional layer below this one by adding sqrt(nChildren) children each with sqrt(nChildren) of our current children */
  virtual void splitTreeSquareRoot( const std::vector<int>& twoLinksStartBlockA, const std::vector<int>& twoLinksStartBlockC ) = 0;

  // TODO : make sure that none of the not suitable methods get called...
  virtual sTree* switchToHierarchicalTree( int nx_to_shave, int myl_to_shave, int mzl_to_shave, const std::vector<int>& twoLinksStartBlockA,
        const std::vector<int>& twoLinksStartBlockC ) = 0;
  virtual void collapseHierarchicalTree() = 0;

};

#endif 
