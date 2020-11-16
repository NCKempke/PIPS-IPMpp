/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef STOCH_TREE_CALLBACKS
#define STOCH_TREE_CALLBACKS

#include "sTree.h"
#include "StochInputTree.h"
#include <functional>
/** This class creates objects when  the problem is specified by C callbacks.
 *  Obsolete and present only to ensure compatibility with older versions of the code.
 *  The new sTree implementation, C++-like is sTreeImpl.
 */

class sData;


class sTreeCallbacks : public sTree
{
   using InputNode = StochInputTree::StochInputNode;
   using DATA_MAT = FMAT InputNode::*;
   using DATA_VEC = FVEC InputNode::*;
   using DATA_NNZ = FNNZ InputNode::*;
   using DATA_INT = int InputNode::*;

private:
  virtual void initPresolvedData(const StochSymMatrix& Q, const StochGenMatrix& A, const StochGenMatrix& C,
        const StochVector& nxVec, const StochVector& myVec, const StochVector& mzVec, int mylParent, int mzlParent);


 public:
  sTreeCallbacks(StochInputTree* root);
  sTreeCallbacks(const std::vector<StochInputTree::StochInputNode*> &localscens);
  sTreeCallbacks(StochInputTree::StochInputNode* data_);
  ~sTreeCallbacks();

  StochSymMatrix*   createQ() const;

 private:
  StochGenMatrix* createMatrix( DATA_INT m_ABmat, DATA_INT n_Mat,
        DATA_INT nnzAmat, DATA_NNZ fnnzAmat, DATA_MAT Amat, DATA_INT nnzBmat,
        DATA_NNZ fnnzBmat, DATA_MAT Bmat, DATA_INT m_Blmat, DATA_INT nnzBlmat,
        DATA_NNZ fnnzBlmat, DATA_MAT Blmat ) const;

  StochVector* createVector( DATA_INT n_vec, DATA_VEC vec, DATA_INT n_linking_vec, DATA_VEC linking_vec ) const;
 public:
  StochVector*      createc() const;

  StochVector*      createxlow()  const;
  StochVector*      createixlow() const;
  StochVector*      createxupp()  const;
  StochVector*      createixupp() const;


  StochGenMatrix*   createA() const;
  StochVector*      createb() const;


  StochGenMatrix*   createC() const;
  StochVector*      createclow()  const;
  StochVector*      createiclow() const;
  StochVector*      createcupp()  const;
  StochVector*      createicupp() const;

  int nx() const;
  int my() const;
  int myl() const;
  int mz() const; 
  int mzl() const;
  int id() const override;

  void computeGlobalSizes();
 public:
  /* inactive sizes store the original state of the tree when switching to the presolved data */
  long long N_INACTIVE, MY_INACTIVE, MZ_INACTIVE, MYL_INACTIVE, MZL_INACTIVE; //global inactive sizes

  int nx_active, my_active, mz_active, myl_active, mzl_active;
  int nx_inactive, my_inactive, mz_inactive, myl_inactive, mzl_inactive;

  void loadLocalSizes();

  virtual void switchToPresolvedData();
  virtual void switchToOriginalData();
  virtual bool isPresolved();
  virtual bool hasPresolved();
  virtual void initPresolvedData(const sData& presolved_data);

  virtual void writeSizes(ostream& sout) const;

  // TODO : make sure that none of the not suitable methods get called...
 private:
  void createSubcommunicatorsAndChildren( std::vector<int>& map_my_procs_to_sub_comm );

  void splitTreeSquareRoot( const std::vector<int>& twoLinksStartBlockA, const std::vector<int>& twoLinksStartBlockC ) override;
  sTree* shaveDenseBorder( int nx_to_shave, int myl_to_shave, int mzl_to_shave) override;
 public:
  sTree* switchToHierarchicalTree( int nx_to_shave, int myl_to_shave, int mzl_to_shave, const std::vector<int>& twoLinksStartBlockA,
        const std::vector<int>& twoLinksStartBlockC ) override;

 protected:
  bool isDataPresolved;
  bool hasPresolvedData;

  sTreeCallbacks();
  InputNode* data; //input data

  // in POOLSCEN case, only root node has non-null data
  StochInputTree* tree;
  std::vector<InputNode*> scens;


  // TODO : remove
  InputNode* fakedata; //convenient struct for holding n,my,mz etc
  // holds stoch trees for each of the scenarios that are combined at this node
  // this is just a convenience to reuse the create* and newVector* functions
  std::vector<sTreeCallbacks*> real_children;
};


#endif
