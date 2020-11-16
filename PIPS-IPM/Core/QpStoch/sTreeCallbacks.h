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

public:
   sTreeCallbacks(StochInputTree* root);
   sTreeCallbacks(StochInputTree::StochInputNode* data_);
   ~sTreeCallbacks();

   StochSymMatrix*   createQ() const;

   StochGenMatrix* createA() const override;
   StochGenMatrix* createC() const override;

   StochVector* createc() const override;

   StochVector* createxlow() const override;
   StochVector* createixlow() const override;
   StochVector* createxupp() const override;
   StochVector* createixupp() const override;

   StochVector* createb() const override;
   StochVector* createclow() const override;
   StochVector* createiclow() const override;
   StochVector* createcupp() const override;
   StochVector* createicupp() const override;

   int nx() const override;
   int my() const override;
   int myl() const override;
   int mz() const override;
   int mzl() const override;
   int id() const override;

   void computeGlobalSizes() override;
   void loadLocalSizes();

   virtual void switchToPresolvedData();
   virtual void switchToOriginalData();
   virtual bool isPresolved();
   virtual bool hasPresolved();
   virtual void initPresolvedData(const sData& presolved_data);

   virtual void writeSizes( std::ostream& sout ) const;

   sTree* switchToHierarchicalTree( int nx_to_shave, int myl_to_shave, int mzl_to_shave, const std::vector<int>& twoLinksStartBlockA,
         const std::vector<int>& twoLinksStartBlockC ) override;

   void splitDataAccordingToTree( sData& data ) const;
   const std::vector<unsigned int>& getMapProcsSubcomms() const
      { assert( is_hierarchical_inner ); return map_proc_subcomm; };
   const std::vector<unsigned int>& getMapBlockSubcomms() const
      { assert( is_hierarchical_inner ); return map_proc_subcomm; };
private:
   virtual void initPresolvedData(const StochSymMatrix& Q, const StochGenMatrix& A, const StochGenMatrix& C,
         const StochVector& nxVec, const StochVector& myVec, const StochVector& mzVec, int mylParent, int mzlParent);

   StochGenMatrix* createMatrix( DATA_INT m_ABmat, DATA_INT n_Mat,
         DATA_INT nnzAmat, DATA_NNZ fnnzAmat, DATA_MAT Amat, DATA_INT nnzBmat,
         DATA_NNZ fnnzBmat, DATA_MAT Bmat, DATA_INT m_Blmat, DATA_INT nnzBlmat,
         DATA_NNZ fnnzBlmat, DATA_MAT Blmat ) const;
   StochVector* createVector( DATA_INT n_vec, DATA_VEC vec, DATA_INT n_linking_vec, DATA_VEC linking_vec ) const;

   void splitMatrixAccordingToTree( StochSymMatrix& mat ) const;
   void splitMatrixAccordingToTree( StochGenMatrix& mat ) const;
   void splitVectorAccordingToTree( StochVector& vec ) const;

   void createSubcommunicatorsAndChildren( std::vector<unsigned int>& map_my_procs_to_sub_comm );
   void splitTreeSquareRoot( const std::vector<int>& twoLinksStartBlockA, const std::vector<int>& twoLinksStartBlockC ) override;
   sTree* shaveDenseBorder( int nx_to_shave, int myl_to_shave, int mzl_to_shave) override;

   /* inactive sizes store the original state of the tree when switching to the presolved data */
   long long N_INACTIVE, MY_INACTIVE, MZ_INACTIVE, MYL_INACTIVE, MZL_INACTIVE; //global inactive sizes

   int nx_active, my_active, mz_active, myl_active, mzl_active;
   int nx_inactive, my_inactive, mz_inactive, myl_inactive, mzl_inactive;

   bool isDataPresolved;
   bool hasPresolvedData;

   /* after this node has been split this will indicate how the procs were split */
   std::vector<unsigned int> map_proc_subcomm;
   std::vector<unsigned int> map_block_subcomm;

   sTreeCallbacks();
   InputNode* data; //input data
};


#endif
