/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#include "sTree.h"
#include "sData.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "SimpleVector.h"
#include "DoubleMatrixTypes.h"
#include "pipsport.h"

#include <numeric>
#include <algorithm>
#include <cmath>

StochIterateResourcesMonitor sTree::iterMon;
int sTree::rankMe    =-1;
int sTree::rankZeroW = 0;
int sTree::rankPrcnd =-1;
int sTree::numProcs  =-1;

sTree::sTree(const sTree& other) : commWrkrs{ other.commWrkrs },
      myProcs( other.myProcs.begin(), other.myProcs.end() ),
      myOldProcs( other.myOldProcs.begin(), other.myOldProcs.end() ),
      commP2ZeroW{ other.commP2ZeroW },
      N{ other.N },
      MY{ other.MY },
      MZ{ other.MZ },
      MYL{ other.MYL },
      MZL{ other.MZL },
      np{ other.np },
      IPMIterExecTIME{ other.IPMIterExecTIME },
      is_hierarchical_root{ other.is_hierarchical_root },
      is_hierarchical_inner_root{ other.is_hierarchical_inner_root },
      is_hierarchical_inner_leaf{ other.is_hierarchical_inner_leaf }
{
   if( other.sub_root )
      sub_root = other.sub_root->clone();
   for( auto& child : other.children )
      children.push_back( child->clone() );
}

sTree::~sTree() 
{
  for(size_t it = 0; it < children.size(); it++)
    delete children[it];
  delete sub_root;
}

int sTree::myl() const
{
   return -1;
}

int sTree::mzl() const
{
   return -1;
}

bool sTree::distributedPreconditionerActive() const
{
   return ( rankZeroW != 0 ) && ( rankPrcnd != -1 ) &&
         ( commP2ZeroW != MPI_COMM_NULL ) && ( rankMe != -1 );
}

void sTree::assignProcesses(MPI_Comm comm)
{
   assert( comm != MPI_COMM_NULL );
   assert( !is_hierarchical_root );

   const int world_size = PIPS_MPIgetSize(comm);
   /* assign root node */
   commWrkrs = MPI_COMM_WORLD;

   std::vector<int> processes(world_size);
   for( int p = 0; p < world_size; p++)
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
   PIPS_MPIabortIf( n_procs > children.size(), ss.str() );

   std::vector<unsigned int> map_child_nodes_to_procs;
   mapChildrenToNSubTrees(map_child_nodes_to_procs, children.size(), n_procs);

#ifndef NDEBUG
   for( size_t i = 0; i < children.size(); i++ )
      assert( map_child_nodes_to_procs[i] < n_procs );

   for( size_t i = 1; i < children.size(); i++ )
      assert( map_child_nodes_to_procs[i - 1] <= map_child_nodes_to_procs[i] );

   std::vector<size_t> load_per_proc(n_procs, 0);
   for( size_t i = 0; i < children.size(); ++i )
      load_per_proc[ map_child_nodes_to_procs[i] ]++;

   const size_t max_load = *std::max_element(load_per_proc.begin(), load_per_proc.end());
   const size_t min_load = *std::min_element(load_per_proc.begin(), load_per_proc.end());
   assert( max_load == min_load || max_load == min_load + 1 );
#endif

   for( size_t i = 0; i < children.size(); i++ )
   {
      children[i]->myProcs.clear();
      children[i]->myProcs.push_back( map_child_nodes_to_procs[i] );
      assert( rankMe >= 0 );
      if( static_cast<unsigned int>(rankMe) == map_child_nodes_to_procs[i] )
         children[i]->commWrkrs = MPI_COMM_SELF;
      else
         children[i]->commWrkrs = MPI_COMM_NULL;
   }
}

double sTree::processLoad() const
{
   assert( !is_hierarchical_root );
  //! need a recursive and also a collective call
  if( IPMIterExecTIME < 0.0 )
     //return (NNZQ+NNZA+NNZB+NNZC+NNZD + N+MY+MZ)/1000.0;
     return (N + MY + MZ) / 1000.0;
  return IPMIterExecTIME;
}

void sTree::getGlobalSizes(long long& n, long long& my, long long& mz) const
{
   n = N; my = MY; mz = MZ;
}

void sTree::getGlobalSizes(long long& n, long long& my, long long& myl, long long& mz, long long& mzl) const
{
   n = N; my = MY; mz = MZ; myl = MYL, mzl = MZL;
}

int sTree::innerSize(int which) const
{
   assert( false );
  if(which==0) return nx();
  if(which==1) return my();
  assert(which==2);
  return mz();
}

StochVector* sTree::newPrimalVector(bool empty) const
{
   if( commWrkrs == MPI_COMM_NULL )
      return new StochDummyVector();

   StochVector* x{};
   if( !sub_root )
   {
      x = new StochVector( empty ? 0 : nx(), commWrkrs);

      for(size_t it = 0; it < children.size(); it++)
      {
         StochVector* child = children[it]->newPrimalVector(empty);
         x->AddChild(child);
      }
   }
   else
   {
      assert( children.size() == 0 );
      if( sub_root->commWrkrs == MPI_COMM_NULL )
         x = new StochDummyVector();
      else
      {
         StochVector* x_vec = sub_root->newPrimalVector(empty);
         x = new StochVector( x_vec, nullptr, commWrkrs );
         x_vec->parent = x;
      }
   }

   return x;
}

StochVector* sTree::newDualYVector(bool empty) const
{
   if( commWrkrs == MPI_COMM_NULL )
      return new StochDummyVector();

   StochVector* y{};
   const int yl = (np == -1) ? myl() : -1;

   if( !sub_root )
   {
      y = new StochVector(empty ? std::min(0, my()) : my(),
            empty ? std::min(yl, 0) : yl, commWrkrs);

      for(size_t it = 0; it < children.size(); it++)
      {
         StochVector* child = children[it]->newDualYVector(empty);
         y->AddChild(child);
      }
   }
   else
   {
      assert( children.size() == 0 );

      if( sub_root->commWrkrs == MPI_COMM_NULL )
         y = new StochDummyVector();
      else
      {
         StochVector* y_vec = sub_root->newDualYVector( empty );
         assert( yl == -1 );

         y = new StochVector( y_vec, nullptr, commWrkrs );
         y_vec->parent = y;
      }
   }
   return y;
}

StochVector* sTree::newDualZVector(bool empty) const
{
   if( commWrkrs == MPI_COMM_NULL )
      return new StochDummyVector();

   StochVector* z{};
   const int zl = (np == -1) ? mzl() : -1;

   if( !sub_root )
   {

      z = new StochVector( empty ? std::min(mz(), 0) : mz(), empty ? std::min(zl, 0) : zl, commWrkrs);
      for( size_t it = 0; it < children.size(); it++)
      {
         StochVector* child = children[it]->newDualZVector(empty);
         z->AddChild(child);
      }
   }
   else
   {
      assert( children.size() == 0 );
      if( sub_root->commWrkrs == MPI_COMM_NULL )
         z = new StochDummyVector();
      else
      {
         StochVector* z_vec = sub_root->newDualZVector(empty);
         assert( zl == -1 );

         z = new StochVector( z_vec, nullptr, commWrkrs );
         z_vec->parent = z;
      }
   }

   return z;
}

StochVector* sTree::newRhs() const
{
  //is this node a dead-end for this process?
  if( commWrkrs == MPI_COMM_NULL )
    return new StochDummyVector();

  StochVector* rhs{};
  if( !sub_root )
  {
     int locmyl = (np == -1) ? myl() : 0;
     int locmzl = (np == -1) ? mzl() : 0;

     locmyl = std::max(locmyl, 0);
     locmzl = std::max(locmzl, 0);

     rhs = new StochVector(nx() + std::max(my(), 0) + std::max(mz(), 0) + locmyl + locmzl, commWrkrs);

     for(size_t it = 0; it < children.size(); it++)
     {
        StochVector* child = children[it]->newRhs();
        rhs->AddChild(child);
     }
  }
  else
     rhs = sub_root->newRhs();

  return rhs;
}

void sTree::startMonitors()
{
  iterMon.recIterateTm_start();
  startNodeMonitors();
}

void sTree::stopMonitors()
{
  iterMon.recIterateTm_stop();
  stopNodeMonitors();
}

void sTree::startNodeMonitors()
{
  resMon.reset();
  for(size_t i=0; i<children.size(); i++)
    children[i]->startNodeMonitors();
}

void sTree::stopNodeMonitors()
{
  for(size_t i=0; i<children.size(); i++)
    children[i]->stopNodeMonitors();

  resMon.computeTotal();
}

void sTree::toMonitorsList(std::list<NodeExecEntry>& lstExecTm)
{
  lstExecTm.push_back(resMon.eTotal);

  for(size_t i=0; i<children.size(); i++)
    children[i]->toMonitorsList(lstExecTm);
}

void sTree::fromMonitorsList(std::list<NodeExecEntry>& lstExecTm)
{
  resMon.eTotal = lstExecTm.front();
  lstExecTm.pop_front();

  for(size_t i=0; i<children.size(); i++)
    children[i]->fromMonitorsList(lstExecTm);
}

void sTree::syncMonitoringData(std::vector<double>& vCPUTotal)
{

  std::list<NodeExecEntry> lstExecTm;
  this->toMonitorsList(lstExecTm);

  int noNodes = lstExecTm.size(); 
  int nCPUs   = vCPUTotal.size();

  double* recvBuf = new double[noNodes+nCPUs];
  double* sendBuf = new double[noNodes+nCPUs];
  
  std::list<NodeExecEntry>::iterator iter = lstExecTm.begin();
  for(int it=0; it<noNodes; it++) { sendBuf[it] = iter->tmChildren; iter++; }
  for(int it=noNodes; it<noNodes+nCPUs; it++) sendBuf[it] = vCPUTotal[it-noNodes];

  
  MPI_Allreduce(sendBuf, recvBuf, noNodes+nCPUs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  iter = lstExecTm.begin();
  for(int it=0; it<noNodes; it++) { iter->tmChildren = recvBuf[it]; iter++; }
  for(int it=noNodes; it<noNodes+nCPUs; it++) vCPUTotal[it-noNodes] = recvBuf[it];

  if(children.size()>0) {
    //local time MPI_MAX, but only for nonleafs
    iter = lstExecTm.begin();
    for(int it=0; it<noNodes; it++) { sendBuf[it] = iter->tmLocal; iter++; }

    MPI_Allreduce(sendBuf, recvBuf, noNodes, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    iter = lstExecTm.begin();
    for(int it=0; it<noNodes; it++) { iter->tmLocal = recvBuf[it]; iter++; }
  }
  delete[] recvBuf; delete[] sendBuf;

  //populate the tree with the global data
  this->fromMonitorsList(lstExecTm);

  //compute syncronized total time for each node of the tree, i.e., local+childs+subchilds
  computeNodeTotal(); //updates this->IPMIterExecTIME
}


bool sTree::balanceLoad()
{
  return false; //disabled for now
  //before synchronization, compute the total time recorded on this CPU
  //updates this->IPMIterExecTIME
  computeNodeTotal();

  int nCPUs; MPI_Comm_size(commWrkrs, &nCPUs);
  std::vector<double> cpuExecTm(nCPUs, 0.0);
  cpuExecTm[rankMe] = this->IPMIterExecTIME;

  //if(sleepFlag && rankMe==1) {cpuExecTm[1] += 4.0;}

  this->syncMonitoringData(cpuExecTm);

  //!log
#if defined(STOCH_TESTING) && !defined(TIMING)
  if(!rankMe) {
    printf("Iteration times per process:\n");
    for(int it=0; it<nCPUs; it++) printf("%8.4f ", cpuExecTm[it]);
    printf("\n\n");
    this->displayExecTimes(0);
  }
#endif
  if(nCPUs==1) return 0;

  double total = 0.0; double maxLoad=0, minLoad=1.e+10;
  for(int it=0; it<nCPUs; it++) {
    total += cpuExecTm[it];
    if(maxLoad<cpuExecTm[it]) maxLoad = cpuExecTm[it];
    if(minLoad>cpuExecTm[it]) minLoad = cpuExecTm[it];
  }

  if(maxLoad<1.0) return 0;

  double balance = std::max(maxLoad/total*nCPUs, total/nCPUs/minLoad);
  if(balance<1.3) { //it is OK, no balancing

    //decide if balancing is needed due to the 'fat' nodes
    total = 0.0; maxLoad=0;
    for(size_t i=0; i<children.size(); i++) {
      total += children[i]->IPMIterExecTIME;
      
      //if(!rankMe) printf("%7.4f ", children[i]->IPMIterExecTIME);
      
      if(maxLoad<children[i]->IPMIterExecTIME) 
	maxLoad = children[i]->IPMIterExecTIME;
    }
    //if(!rankMe) printf("\n");//aaa
    
    balance = maxLoad/total*children.size();

    //if(!rankMe) printf("CPU node balance=%g\n", balance);

    if(balance<2.00) return 0;
    //else if(!rankMe) printf("!!! Balancing NEEDED: node balance=%g\n", balance);
  } else {
    //if(!rankMe) printf("!!! Balancing NEEDED: CPU balance=%g\n", balance);
  }
  return 0;
  //save the current MPI related information
  saveCurrentCPUState();

  cpuExecTm.clear();

  std::vector<int> ranks(nCPUs);
  for(int i=0; i<nCPUs; i++) ranks[i]=i; 
  assignProcesses(MPI_COMM_WORLD, ranks);
  return 1;
}

#define maSend 1
#define maRecv 0

void sTree::getSyncInfo(int rank, int& syncNeeded, int& sendOrRecv, int& toFromCPU )
{
  // was this node previously assigned to cpu 'rank'?
  int wasAssigned = isInVector(rank, myOldProcs);
  // is currently assigned to cpu 'rank'?
  int isAssigned = isInVector(rank, myProcs);

  //if(0==wasAssigned && 0==isAssigned) return;
  //if(1==wasAssigned && 1==isAssigned) return;
  syncNeeded=0;  sendOrRecv = 0; toFromCPU = -1;
  if(isAssigned!=wasAssigned) {
    syncNeeded=1;
    
    if(wasAssigned) {
      assert(0==isAssigned); assert(myProcs.size()>0);
      //where is this node assigned?
      toFromCPU=myProcs[0];
      sendOrRecv = maSend;
    } else {
      assert(1==isAssigned); assert(myOldProcs.size()>0);
      toFromCPU=myOldProcs[0];
      sendOrRecv = maRecv;
    }
  }
}

void sTree::computeNodeTotal()
{
  if(0==children.size())
    this->IPMIterExecTIME = resMon.eTotal.tmChildren;
  else {

    this->IPMIterExecTIME = resMon.eTotal.tmLocal;
    for(size_t i=0; i<children.size(); i++) {
      children[i]->computeNodeTotal();
      
      this->IPMIterExecTIME += children[i]->IPMIterExecTIME;
    }
  }
}

void sTree::saveCurrentCPUState()
{
  myOldProcs = myProcs;

  for(size_t i=0; i<children.size(); i++)
    children[i]->saveCurrentCPUState();
}


void sTree::printProcessTree() const
{
   if( rankMe != 0 )
      return;

   std::cout << "Process Tree:\n\n";

   std::vector<const sTree*> queue;
   queue.insert(queue.end(), this);

   size_t curr_size = queue.size();
   size_t count = 0;
   bool node_layer_reached = false;
   bool node_layer_next = false;
   bool tree_unbalanced = false;
   while( !queue.empty() )
   {
      const sTree* child = queue.front();
      if( !child )
      {
         std::cout << "| ";
         queue.erase(queue.begin());
         ++count;
      }
      else if( node_layer_reached && child->myProcs.size() == 1 )
      {
         static size_t node = 1;
         assert(child);

         const size_t node_begin = node;
         const int curr_proc = child->myProcs[0];
         queue.erase(queue.begin());
         ++count;

         size_t node_end = node;

         while( !queue.empty() && queue.front() != nullptr && queue.front()->myProcs[0] == curr_proc )
         {
            ++node_end;
            queue.erase(queue.begin());
            ++count;
         }

         if( node_begin != node_end )
            std::cout << "[ " << curr_proc << ": " << node_begin << "-" << node_end << " ] ";
         else
            std::cout << "[ " << curr_proc << ": " << node_begin << " ] ";
         node = node_end + 1;
      }
      else
      {
         if( node_layer_reached )
            tree_unbalanced = true;

         assert( child );
         assert( !child->sub_root );
         if( child->myProcs.size() >= 1 )
         {
            if( child->myProcs.size() > 1 )
               std::cout << "[ " << child->myProcs.front() << "-" << child->myProcs.back() << " ] ";
            else
               std::cout << "[ " << child->myProcs.front() << " ] ";
         }

         for( auto& childchild : child->children )
         {
            if( childchild->sub_root )
               queue.insert(queue.end(), childchild->sub_root);
            else
            {
               node_layer_next = true;
               queue.insert(queue.end(), childchild);
            }
         }
         if( !(child->children.empty()) )
            queue.insert(queue.end(), nullptr);
         queue.erase(queue.begin());
         ++count;
      }

      if( count == curr_size )
      {
         if( node_layer_next )
            node_layer_reached = true;
         curr_size = queue.size();
         count = 0;
         std::cout << "\n\n";
      }
   }

   if( tree_unbalanced )
      std::cout << "WARNING : the split tree is unbalanced and there is branches that are longer than others\n" <<
         "\t-> this is bad for performance but that fine of a split usually is bad for performance anyway..\n" <<
         "\t-> rethink the amount of layers for the hierarchical approach..\n\n";
}

void sTree::mapChildrenToNSubTrees( std::vector<unsigned int>& map_child_to_sub_tree, unsigned int n_children, unsigned int n_subtrees )
{
   assert( n_subtrees <= n_children );
   map_child_to_sub_tree.clear();
   map_child_to_sub_tree.reserve(n_children);

   if( n_subtrees == 0 )
      return;

   const unsigned int everyone_gets = std::floor(n_children / n_subtrees);
   const unsigned int n_leftovers = n_children % n_subtrees;

   std::vector<unsigned int> children_per_tree( n_subtrees, everyone_gets );

   if( n_leftovers > 0 )
   {
      const unsigned int free_in_leftover_row = n_subtrees - n_leftovers;

      const unsigned int free_after_each_leftover = std::floor(free_in_leftover_row / n_leftovers);
      const unsigned int additional_frees = free_in_leftover_row % n_leftovers;

      unsigned int additional_frees_assigned = 0;
      for( unsigned int i = 0; i < n_subtrees; i += free_after_each_leftover )
      {
         ++children_per_tree[i];
         ++i;

         if( additional_frees != additional_frees_assigned )
         {
            ++i;
            ++additional_frees_assigned;
         }
      }
      assert( additional_frees_assigned == additional_frees );
   }

   assert( std::accumulate( children_per_tree.begin(), children_per_tree.end(),
         decltype(children_per_tree)::value_type(0) ) == n_children );

   for( unsigned int i = 0; i < children_per_tree.size(); ++i )
   {
      for( unsigned int j = 0; j < children_per_tree[i]; ++j )
         map_child_to_sub_tree.push_back( i );
   }
}
