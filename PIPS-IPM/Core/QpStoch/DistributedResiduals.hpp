#ifndef STOCHRESIDUALS
#define STOCHRESIDUALS

#include "Residuals.h"
#include "DistributedVector.h"
#include <vector>

class sTree;

/** 
 * Class added to supply a more generic constructor for its parent, Residuals.
 * The default constructor of Residuals can not be always used since it assumes that vectors
 * are arrays. For stochastic problems, the vectors are trees of arrays.
 *
 * @ingroup QpGen
 */

class DistributedResiduals : public Residuals {
public:
  /**
   * Constructor
   */
  DistributedResiduals( OoqpVector * rQ,
	      OoqpVector * rA, OoqpVector * rC,
	      OoqpVector * rz,
	      OoqpVector * rt, OoqpVector * rlambda,
	      OoqpVector * ru, OoqpVector * rpi,
	      OoqpVector * rv, OoqpVector * rgamma,
	      OoqpVector * rw, OoqpVector * rphi,
	      OoqpVector * ixlow, double nxlowGlobal,
	      OoqpVector * ixupp, double nxuppGlobal,
	      OoqpVector * iclow, double mclowGlobal,
	      OoqpVector * icupp, double mcuppGlobal );

  DistributedResiduals( const sTree* tree, OoqpVector * ixlow_, OoqpVector * ixupp_, OoqpVector * iclow_, OoqpVector * icupp_ );

  DistributedResiduals( const DistributedResiduals& res );

  ~DistributedResiduals() override;

  void permuteVec0Entries( const std::vector<unsigned int>& perm, bool resids_only = false );
  void permuteEqLinkingEntries( const std::vector<unsigned int>& perm );
  void permuteIneqLinkingEntries( const std::vector<unsigned int>& perm, bool resids_only = false );

  bool isRootNodeInSync() const;

  void collapseHierarchicalStructure( const DistributedQP& data, const sTree* tree_hier, OoqpVectorHandle ixlow_, OoqpVectorHandle ixupp_, OoqpVectorHandle iclow_, OoqpVectorHandle icupp_);

  std::vector<DistributedResiduals*> children;
protected:
  void createChildren();
  void AddChild(DistributedResiduals* child);
};

#endif
