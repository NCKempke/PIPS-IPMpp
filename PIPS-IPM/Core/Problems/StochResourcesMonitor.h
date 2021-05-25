#ifndef STOCH_RESOURCE_MON
#define STOCH_RESOURCE_MON

#include <vector>
#include <string>

//! not thread safe

enum stCommType { ctAllreduce = 0, ctReduce, ctOther };

class NodeTimer {
public:
   NodeTimer(double localTime, double childrenTime);
   void clear();

   double local_time, children_time;

   void start_local();
   void stop_local();

   void start_children();
   void stop_children();

protected:
   double local_start_time, children_start_time;
};

class NodeCommEntry {
public:
   NodeCommEntry(double time_, double size_, stCommType type_) : time(time_), size(size_), type(type_) {};
   double time, size;
   stCommType type;
};

class StochNodeResourcesMonitor {
public:
   StochNodeResourcesMonitor();
   virtual ~StochNodeResourcesMonitor();

   virtual void reset();
   virtual void computeTotal();

   virtual void recFactTmLocal_start();
   virtual void recFactTmLocal_stop();
   virtual void recFactTmChildren_start();
   virtual void recFactTmChildren_stop();

   virtual void recLsolveTmChildren_start();
   virtual void recLsolveTmChildren_stop();

   virtual void recDsolveTmLocal_start();
   virtual void recDsolveTmLocal_stop();
   virtual void recDsolveTmChildren_start();
   virtual void recDsolveTmChildren_stop();

   virtual void recLtsolveTmLocal_start();
   virtual void recLtsolveTmLocal_stop();
   virtual void recLtsolveTmChildren_start();
   virtual void recLtsolveTmChildren_stop();

   virtual void recSchurCom_start(double size, stCommType type);
   virtual void recSchurCom_stop(double, stCommType);

   //record the dense matrix multiplication within Schur complement computation
   //the recorded time is part of the recFactTmChildren kept in eFact.tmChildren
   virtual void recSchurMultLocal_start();
   virtual void recSchurMultLocal_stop();

   virtual void recReduceTmLocal_start();
   virtual void recReduceTmLocal_stop();

   virtual void recReduceScatterTmLocal_start();
   virtual void recReduceScatterTmLocal_stop();


public:

   NodeTimer eFact, eLsolve, eDsolve, eLtsolve;
   NodeTimer eTotal;
   NodeTimer eMult;
   NodeTimer eReduce, eReduceScatter, eBcast;
   std::vector<NodeCommEntry> vcSchur, vcLsolve;

private:
   double tmOpStart;
};

class Timer {
public:
   Timer() = default;
   virtual ~Timer() = default;

   virtual void start();
   virtual void stop();

public:
   double end_time;
protected:
   double start_time;
};

#endif
