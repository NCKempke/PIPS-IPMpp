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

   [[nodiscard]] std::unique_ptr<Residuals> clone_full() const override;

   ~DistributedResiduals() override = default;

   void permute_vec_0_entries(const std::vector<unsigned int>& perm, bool resids_only = false);

   void permute_eq_linking_entries(const std::vector<unsigned int>& perm);

   void permute_ineq_linking_entries(const std::vector<unsigned int>& perm, bool resids_only = false);

   [[nodiscard]] bool is_root_node_in_sync() const;

   void collapse_hierarchical_structure(const DistributedProblem& data_hier, const DistributedTree* tree_hier);
   void update_indicators(std::shared_ptr<Vector<double>> ixlow_, std::shared_ptr<Vector<double>> ixupp_,
      std::shared_ptr<Vector<double>> iclow_, std::shared_ptr<Vector<double>> icupp_);
};

#endif
