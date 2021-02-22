/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef STOCH_TREE_CALLBACKS
#define STOCH_TREE_CALLBACKS

#include "sTree.h"
#include "StochInputTree.h"

/** This class creates objects when  the problem is specified by C callbacks.
 *  Obsolete and present only to ensure compatibility with older versions of the code.
 *  The new sTree implementation, C++-like is sTreeImpl.
 */

class sTreeCallbacks : public sTree
{

 private:
  virtual void initPresolvedData(const StochSymMatrix& Q, const StochGenMatrix& A, const StochGenMatrix& C,
        const StochVector& nxVec, const StochVector& myVec, const StochVector& mzVec, int mylParent, int mzlParent);


 public:
  sTreeCallbacks(StochInputTree* root);
  sTreeCallbacks(const std::vector<StochInputTree::StochInputNode*> &localscens);
  sTreeCallbacks(StochInputTree::StochInputNode* data_);
  ~sTreeCallbacks();

  StochSymMatrix*   createQ() const override;
  StochVector*      createc() const override;

  StochVector*      createxlow()  const override;
  StochVector*      createixlow() const override;
  StochVector*      createxupp()  const override;
  StochVector*      createixupp() const override;


  StochGenMatrix*   createA() const override;
  StochVector*      createb() const override;


  StochGenMatrix*   createC() const override;
  StochVector*      createclow()  const override;
  StochVector*      createiclow() const override;
  StochVector*      createcupp()  const override;
  StochVector*      createicupp() const override;

  int nx() const override;
  int my() const override;
  int myl() const override;
  int mz() const override; 
  int mzl() const override;
  int id() const override; 

  void computeGlobalSizes() override;
  
  int NNZA,NNZQ,NNZB,NNZBl,NNZC,NNZD,NNZDl; //global nnz
  int NNZA_INACTIVE,NNZQ_INACTIVE,NNZB_INACTIVE,NNZBl_INACTIVE,NNZC_INACTIVE,NNZD_INACTIVE,NNZDl_INACTIVE; //global inactive nnz
  long long N_INACTIVE,MY_INACTIVE,MZ_INACTIVE; //global inactive sizes
  int nx_active, my_active, mz_active, myl_active, mzl_active;
  int nx_inactive, my_inactive, mz_inactive, myl_inactive, mzl_inactive;

  void loadLocalSizes() override;

  virtual void switchToPresolvedData();
  virtual void switchToOriginalData();
  virtual bool isPresolved();
  virtual bool hasPresolved();
  virtual void initPresolvedData(const StochSymMatrix& Q, const StochGenMatrix& A, const StochGenMatrix& C, const StochVector& nxVec, const StochVector& myVec, const StochVector& mzVec)
  {
     initPresolvedData(Q, A, C, nxVec, myVec, mzVec, -1, -1);
  }
  virtual void writeSizes(ostream& sout) const;
 protected:
  bool isDataPresolved;
  bool hasPresolvedData;

  sTreeCallbacks();
  StochInputTree::StochInputNode* data; //input data
  // in POOLSCEN case, only root node has non-null data
  StochInputTree* tree;
  std::vector<StochInputTree::StochInputNode*> scens;
  StochInputTree::StochInputNode* fakedata; //convenient struct for holding n,my,mz etc
  // holds stoch trees for each of the scenarios that are combined at this node
  // this is just a convenience to reuse the create* and newVector* functions
  std::vector<sTreeCallbacks*> real_children;
};


#endif
