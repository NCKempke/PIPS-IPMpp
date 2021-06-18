#ifndef QPGENSTOCHVARS
#define QPGENSTOCHVARS

#include "Variables.h"
#include "Vector.hpp"
#include "DistributedVector.h"
#include <vector>

class QP;

class DistributedTree;

class DistributedVariables : public Variables {

public:
   /** constructor in which the data and variable pointers are set to point to the given arguments */
   DistributedVariables(const DistributedTree* tree, std::unique_ptr<Vector<double>> x_in, std::unique_ptr<Vector<double>> s_in,
      std::unique_ptr<Vector<double>> y_in, std::unique_ptr<Vector<double>> z_in, std::unique_ptr<Vector<double>> v_in,
      std::unique_ptr<Vector<double>> gamma_in, std::unique_ptr<Vector<double>> w_in, std::unique_ptr<Vector<double>> phi_in, std::unique_ptr<Vector<double>> t_in,
      std::unique_ptr<Vector<double>> lambda_in, std::unique_ptr<Vector<double>> u_in, std::unique_ptr<Vector<double>> pi_in,
      std::shared_ptr<Vector<double>> ixlow_in, long long nxlowGlobal, std::shared_ptr<Vector<double>> ixupp_in,
      long long nxuppGlobal, std::shared_ptr<Vector<double>> iclow_in, long long mclowGlobal, std::shared_ptr<Vector<double>> icupp_in,
         long long mcuppGlobal);

   DistributedVariables(const DistributedVariables& vars);

   [[nodiscard]] std::unique_ptr<Variables> cloneFull() const override;

   ~DistributedVariables() override = default;

   [[nodiscard]] bool isRootNodeInSync() const;

   void collapseHierarchicalStructure(const DistributedQP& hier_data, const DistributedTree* stochNode_,
      std::shared_ptr<Vector<double>> ixlow_, std::shared_ptr<Vector<double>> ixupp_, std::shared_ptr<Vector<double>> iclow_,
      std::shared_ptr<Vector<double>> icupp_);

   void permuteVec0Entries(const std::vector<unsigned int>& perm, bool vars_only = false);
   void permuteEqLinkingEntries(const std::vector<unsigned int>& perm);
   void permuteIneqLinkingEntries(const std::vector<unsigned int>& perm, bool vars_only = false);

protected:
   const DistributedTree* stochNode;
};

#endif

