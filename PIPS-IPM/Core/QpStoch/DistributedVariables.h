#ifndef QPGENSTOCHVARS
#define QPGENSTOCHVARS

#include "Variables.h"
#include "Vector.hpp"
#include "SmartPointer.h"
#include "DistributedVector.h"
#include <vector>

class QpGen;

class QP;

class LinearAlgebraPackage;

class sTree;

class DistributedVariables : public Variables {
public:
   /** constructor in which the data and variable pointers are set to point to the given arguments */
   DistributedVariables(const sTree* tree, Vector<double>* x_in, Vector<double>* s_in, Vector<double>* y_in, Vector<double>* z_in, Vector<double>* v_in, Vector<double>* gamma_in,
         Vector<double>* w_in, Vector<double>* phi_in, Vector<double>* t_in, Vector<double>* lambda_in, Vector<double>* u_in, Vector<double>* pi_in, Vector<double>* ixlow_in,
         long long nxlowGlobal, Vector<double>* ixupp_in, long long nxuppGlobal, Vector<double>* iclow_in, long long mclowGlobal, Vector<double>* icupp_in,
         long long mcuppGlobal);

   DistributedVariables(const DistributedVariables& vars);

   virtual ~DistributedVariables();

   bool isRootNodeInSync() const;

   std::vector<DistributedVariables*> children;

   void collapseHierarchicalStructure(const DistributedQP& hier_data, const sTree* stochNode, SmartPointer<Vector<double> > ixlow_, SmartPointer<Vector<double> > ixupp_,
         SmartPointer<Vector<double> > iclow_, SmartPointer<Vector<double> > icupp_);

   void permuteVec0Entries(const std::vector<unsigned int>& perm, bool vars_only = false);
   void permuteEqLinkingEntries(const std::vector<unsigned int>& perm);
   void permuteIneqLinkingEntries(const std::vector<unsigned int>& perm, bool vars_only = false);

   void sync();

protected:
   void createChildren();
   void AddChild(DistributedVariables* child);

   const sTree* stochNode;
};

#endif

