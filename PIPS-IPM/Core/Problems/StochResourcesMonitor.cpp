#include "StochResourcesMonitor.hpp"

// save diagnostic state
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-override"

#include "mpi.h"
// turn the warnings back on
#pragma GCC diagnostic pop

#include <fstream>
#include <sstream>

#define MAX(a, b) ((a > b) ? a : b)

StochNodeResourcesMonitor::StochNodeResourcesMonitor() : eFact(0, 0), eLsolve(0, 0), eDsolve(0, 0), eLtsolve(0, 0), eTotal(0, 0), eMult(0, 0),
      eReduce(0, 0), eReduceScatter(0, 0), eBcast(0, 0) {
}

StochNodeResourcesMonitor::~StochNodeResourcesMonitor() {}

void StochNodeResourcesMonitor::reset() {
   eFact.clear();
   eLsolve.clear();
   eDsolve.clear();
   eLtsolve.clear();
   eTotal.clear();
   eMult.clear();
   eReduce.clear();
   eReduceScatter.clear();
   eBcast.clear();
   vcSchur.clear();
   vcLsolve.clear();
}


void StochNodeResourcesMonitor::computeTotal() {
   eTotal.local_time = eFact.local_time + eLsolve.local_time + eDsolve.local_time + eLtsolve.local_time;
   eTotal.children_time = eFact.children_time + eLsolve.children_time + eDsolve.children_time + eLtsolve.children_time;
}

void StochNodeResourcesMonitor::recFactTmLocal_start() {
   tmOpStart = MPI_Wtime();
}

void StochNodeResourcesMonitor::recFactTmLocal_stop() {
   double tmOpEnd = MPI_Wtime();
   eFact.local_time += (tmOpEnd - tmOpStart);
}

void StochNodeResourcesMonitor::recFactTmChildren_start() { tmOpStart = MPI_Wtime(); }

void StochNodeResourcesMonitor::recFactTmChildren_stop() {
   double tmOpEnd = MPI_Wtime();
   eFact.children_time += (tmOpEnd - tmOpStart);
}


void StochNodeResourcesMonitor::recLsolveTmChildren_start() { tmOpStart = MPI_Wtime(); }

void StochNodeResourcesMonitor::recLsolveTmChildren_stop() {
   double tmOpEnd = MPI_Wtime();
   eLsolve.children_time += (tmOpEnd - tmOpStart);
}

void StochNodeResourcesMonitor::recDsolveTmLocal_start() { tmOpStart = MPI_Wtime(); }

void StochNodeResourcesMonitor::recDsolveTmLocal_stop() {
   double tmOpEnd = MPI_Wtime();
   eDsolve.local_time += (tmOpEnd - tmOpStart);
}

void StochNodeResourcesMonitor::recDsolveTmChildren_start() { tmOpStart = MPI_Wtime(); }

void StochNodeResourcesMonitor::recDsolveTmChildren_stop() {
   double tmOpEnd = MPI_Wtime();
   eDsolve.children_time += (tmOpEnd - tmOpStart);
}

void StochNodeResourcesMonitor::recLtsolveTmLocal_start() { tmOpStart = MPI_Wtime(); }

void StochNodeResourcesMonitor::recLtsolveTmLocal_stop() {
   double tmOpEnd = MPI_Wtime();
   eLtsolve.local_time += (tmOpEnd - tmOpStart);
}
void StochNodeResourcesMonitor::recLtsolveTmChildren_start() { tmOpStart = MPI_Wtime(); }

void StochNodeResourcesMonitor::recLtsolveTmChildren_stop() {
   double tmOpEnd = MPI_Wtime();
   eLtsolve.children_time += (tmOpEnd - tmOpStart);
}


void StochNodeResourcesMonitor::recSchurCom_start(double size, stCommType type) {
   NodeCommEntry commMon(0.0, size, type);
   vcSchur.push_back(commMon);
   tmOpStart = MPI_Wtime();
}

void StochNodeResourcesMonitor::recSchurCom_stop(double, stCommType) {
   double tmOpEnd = MPI_Wtime();
   vcSchur[vcSchur.size() - 1].time = tmOpEnd - tmOpStart;
}

void StochNodeResourcesMonitor::recSchurMultLocal_start() { eMult.start_local(); }

void StochNodeResourcesMonitor::recSchurMultLocal_stop() { eMult.stop_local(); }

void StochNodeResourcesMonitor::recReduceTmLocal_start() { eReduce.start_local(); }

void StochNodeResourcesMonitor::recReduceTmLocal_stop() { eReduce.stop_local(); }


void StochNodeResourcesMonitor::recReduceScatterTmLocal_start() { eReduceScatter.start_local(); }

void StochNodeResourcesMonitor::recReduceScatterTmLocal_stop() { eReduceScatter.stop_local(); }

//**************************************************
//*************  NodeExecEntry   *******************
//**************************************************
NodeTimer::NodeTimer(double localTime, double childrenTime) : local_time(localTime), children_time(childrenTime) {};

void NodeTimer::clear() {
   local_time = 0.0;
   children_time = 0.0;
}

void NodeTimer::start_local() {
   local_start_time = MPI_Wtime();
}

void NodeTimer::stop_local() {
   local_time += (MPI_Wtime() - local_start_time);
};

void NodeTimer::start_children() {
   children_start_time = MPI_Wtime();
}

void NodeTimer::stop_children() {
   children_time += (MPI_Wtime() - children_start_time);
};

//**************************************************
//*************  ITERATE monitor *******************
//**************************************************

void Timer::start() {
   start_time = MPI_Wtime();
}

void Timer::stop() {
   end_time = MPI_Wtime() - start_time;
}