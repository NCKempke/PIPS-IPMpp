#ifndef QPGENSTOCHVARS
#define QPGENSTOCHVARS

#include "QpGenVars.h"
#include "OoqpVectorHandle.h"
#include "StochVector.h"
#include <vector>

class QpGen;
class QpGenData;
class LinearAlgebraPackage;
class sTree;

class sVars : public QpGenVars {
public:
  sVars( const sTree* tree,
	 OoqpVector * ixlow_in, OoqpVector * ixupp_in,
	 OoqpVector * iclow_in, OoqpVector * icupp_in);

  /** constructor in which the data and variable pointers are set to
      point to the given arguments */
  sVars( const sTree* tree,
	 OoqpVector * x_in, OoqpVector * s_in,
	 OoqpVector * y_in, OoqpVector * z_in,
	 OoqpVector * v_in, OoqpVector * gamma_in,
	 OoqpVector * w_in, OoqpVector * phi_in,
	 OoqpVector * t_in, OoqpVector * lambda_in,
	 OoqpVector * u_in, OoqpVector * pi_in,
	 OoqpVector * ixlow_in, long long nxlowGlobal,
	 OoqpVector * ixupp_in, long long nxuppGlobal,
	 OoqpVector * iclow_in, long long mclowGlobal,
	 OoqpVector * icupp_in, long long mcuppGlobal);

  sVars(const sVars& vars);

  virtual ~sVars();
  
  bool isRootNodeInSync() const;

  std::vector<sVars*> children;

  void collapseHierarchicalStructure(const sData& hier_data, const sTree* stochNode, OoqpVectorHandle ixlow_, OoqpVectorHandle ixupp_,
        OoqpVectorHandle iclow_, OoqpVectorHandle icupp_);

  void permuteVec0Entries( const std::vector<unsigned int>& perm, bool vars_only = false );
  void permuteEqLinkingEntries( const std::vector<unsigned int>& perm );
  void permuteIneqLinkingEntries( const std::vector<unsigned int>& perm, bool vars_only = false );

  void sync();

protected:
  void createChildren();
  void AddChild(sVars* child);

  const sTree* stochNode;
};
#endif

