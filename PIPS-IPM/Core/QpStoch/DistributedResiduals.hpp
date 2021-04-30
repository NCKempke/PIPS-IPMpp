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
   DistributedResiduals(Vector<double>* rQ, Vector<double>* rA, Vector<double>* rC, Vector<double>* rz, Vector<double>* rt, Vector<double>* rlambda,
         Vector<double>* ru, Vector<double>* rpi, Vector<double>* rv, Vector<double>* rgamma, Vector<double>* rw, Vector<double>* rphi,
         Vector<double>* ixlow, double nxlowGlobal, Vector<double>* ixupp, double nxuppGlobal, Vector<double>* iclow, double mclowGlobal,
         Vector<double>* icupp, double mcuppGlobal);

   DistributedResiduals(const sTree* tree, Vector<double>* ixlow_, Vector<double>* ixupp_, Vector<double>* iclow_, Vector<double>* icupp_);

   DistributedResiduals(const DistributedResiduals& res);

   ~DistributedResiduals() override;

   void permuteVec0Entries(const std::vector<unsigned int>& perm, bool resids_only = false);
   void permuteEqLinkingEntries(const std::vector<unsigned int>& perm);
   void permuteIneqLinkingEntries(const std::vector<unsigned int>& perm, bool resids_only = false);

   bool isRootNodeInSync() const;

   void collapseHierarchicalStructure(const DistributedQP& data, const sTree* tree_hier, SmartPointer<Vector<double> > ixlow_,
         SmartPointer<Vector<double> > ixupp_, SmartPointer<Vector<double> > iclow_, SmartPointer<Vector<double> > icupp_);

   std::vector<DistributedResiduals*> children;
protected:
   void createChildren();
   void AddChild(DistributedResiduals* child);
};

#endif
