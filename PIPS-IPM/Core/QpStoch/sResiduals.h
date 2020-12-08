#ifndef STOCHRESIDUALS
#define STOCHRESIDUALS

#include "QpGenResiduals.h"
#include "StochVector_fwd.h"
#include <vector>

class sTree;

/** 
 * Class added to supply a more generic constructor for its parent, QpGenResiduals.
 * The default constructor of QpGenResiduals can not be always used since it assumes that vectors
 * are arrays. For stochastic problems, the vectors are trees of arrays.
 *
 * @ingroup QpGen
 */

class sResiduals : public QpGenResiduals {
public:
  /**
   * Constructor
   */
  sResiduals( const sTree* tree,OoqpVector * rQ,
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
  
  sResiduals( const sTree* tree,
	      OoqpVector * ixlow_, OoqpVector * ixupp_,
	      OoqpVector * iclow_, OoqpVector * icupp_ );
  
  sResiduals( const sResiduals& res );

  ~sResiduals() override;

  void sync();

  void permuteVec0Entries( const std::vector<unsigned int>& perm, bool resids_only = false );
  void permuteEqLinkingEntries( const std::vector<unsigned int>& perm );
  void permuteIneqLinkingEntries( const std::vector<unsigned int>& perm, bool resids_only = false );

  bool isRootNodeInSync() const;
  void collapseHierarchicalStructure(const sTree* stochNode_, OoqpVectorHandle ixlow_, OoqpVectorHandle ixupp_,
        OoqpVectorHandle iclow_, OoqpVectorHandle icupp_);

  std::vector<sResiduals*> children;
private:
  void createChildren();
  void AddChild(sResiduals* child);
  
 protected:
  const sTree* stochNode;
};

#endif
