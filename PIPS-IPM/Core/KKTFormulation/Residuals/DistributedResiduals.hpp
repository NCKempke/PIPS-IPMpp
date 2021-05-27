#ifndef STOCHRESIDUALS
#define STOCHRESIDUALS

#include "Residuals.h"
#include "DistributedVector.h"
#include <vector>

class DistributedTree;

/** 
 * Class added to supply a more generic constructor for its parent, Residuals.
 * The default constructor of Residuals can not be always used since it assumes that vectors
 * are arrays. For stochastic problems, the vectors are trees of arrays.
 *
 * @ingroup QpGen
 */

class DistributedResiduals : public Residuals {
public:
   DistributedResiduals(std::unique_ptr<Vector<double>> rQ_, std::unique_ptr<Vector<double>> rA_,
      std::unique_ptr<Vector<double>> rC_, std::unique_ptr<Vector<double>> rz_, std::unique_ptr<Vector<double>> rt_,
      std::unique_ptr<Vector<double>> rlambda_, std::unique_ptr<Vector<double>> ru_,
      std::unique_ptr<Vector<double>> rpi_, std::unique_ptr<Vector<double>> rv_,
      std::unique_ptr<Vector<double>> rgamma_, std::unique_ptr<Vector<double>> rw_,
      std::unique_ptr<Vector<double>> rphi_, std::shared_ptr<Vector<double>> ixlow_,
      std::shared_ptr<Vector<double>> ixupp_, std::shared_ptr<Vector<double>> iclow_,
      std::shared_ptr<Vector<double>> icupp_);

   DistributedResiduals(const DistributedResiduals& res) = default;

   ~DistributedResiduals() override = default;

   void permuteVec0Entries(const std::vector<unsigned int>& perm, bool resids_only = false);

   void permuteEqLinkingEntries(const std::vector<unsigned int>& perm);

   void permuteIneqLinkingEntries(const std::vector<unsigned int>& perm, bool resids_only = false);

   [[nodiscard]] bool isRootNodeInSync() const;

   void collapseHierarchicalStructure(const DistributedQP& data_hier, const DistributedTree* tree_hier,
      std::shared_ptr<Vector<double>> ixlow_,
      std::shared_ptr<Vector<double>> ixupp_, std::shared_ptr<Vector<double>> iclow_,
      std::shared_ptr<Vector<double>> icupp_);

protected:
};

#endif
