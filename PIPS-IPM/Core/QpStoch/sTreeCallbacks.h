/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */
#ifndef STOCH_TREE_CALLBACKS
#define STOCH_TREE_CALLBACKS

#include "sTree.h"
#include "StochInputTree.h"

#include <memory>
#include <functional>
/** This class creates objects when  the problem is specified by C callbacks.
 *  Obsolete and present only to ensure compatibility with older versions of the code.
 *  The new sTree implementation, C++-like is sTreeImpl.
 */

class sData;

class sTreeCallbacks : public sTree
{
protected:
   sTreeCallbacks( const sTreeCallbacks& other );
public:

   using InputNode = StochInputTree::StochInputNode;
   using DATA_MAT = FMAT InputNode::*;
   using DATA_VEC = FVEC InputNode::*;
   using DATA_NNZ = FNNZ InputNode::*;
   using DATA_INT = int InputNode::*;
   using TREE_SIZE = long long sTree::*;

   sTree* clone() const override;

   sTreeCallbacks(StochInputTree* root);
   sTreeCallbacks(StochInputTree::StochInputNode* data_);
   ~sTreeCallbacks() = default;

   void addChild( sTreeCallbacks* child );

   StochSymMatrix* createQ() const override;

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

   virtual void switchToPresolvedData();
   virtual void switchToOriginalData();
   virtual bool isPresolved();
   virtual bool hasPresolved();
   virtual void initPresolvedData(const sData& presolved_data);

   virtual void writeSizes( std::ostream& sout ) const;

   sTree* switchToHierarchicalTree( sData*& data ) override;

   const std::vector<unsigned int>& getMapBlockSubTrees() const
      { assert( is_hierarchical_inner_root ); return map_node_sub_root; };
   std::vector<MPI_Comm> getChildComms() const;

   void assertTreeStructureCorrect() const;

protected:
   void assertTreeStructureChildren() const;
   void assertSubRoot() const;
   void assertTreeStructureIsNotMyNode() const;
   void assertTreeStructureIsMyNodeChildren() const;
   void assertTreeStructureIsMyNodeSubRoot() const;
   void assertTreeStructureIsMyNode()const;

   static unsigned int getMapChildrenToSqrtNSubTrees( std::vector<unsigned int>& map_child_to_sub_tree, unsigned int n_children );

   void initPresolvedData(const StochSymMatrix& Q, const StochGenMatrix& A, const StochGenMatrix& C,
         const StochVector& nxVec, const StochVector& myVec, const StochVector& mzVec, int mylParent, int mzlParent);

   StochGenMatrix* createMatrix( TREE_SIZE my, TREE_SIZE myl, DATA_INT m_ABmat, DATA_INT n_Mat,
         DATA_INT nnzAmat, DATA_NNZ fnnzAmat, DATA_MAT Amat, DATA_INT nnzBmat,
         DATA_NNZ fnnzBmat, DATA_MAT Bmat, DATA_INT m_Blmat, DATA_INT nnzBlmat,
         DATA_NNZ fnnzBlmat, DATA_MAT Blmat ) const;
   StochVector* createVector( DATA_INT n_vec, DATA_VEC vec, DATA_INT n_linking_vec, DATA_VEC linking_vec ) const;

   void createSubcommunicatorsAndChildren( std::vector<unsigned int>& map_child_to_sub_tree );
   void countTwoLinksForChildTrees(const std::vector<int>& two_links_start_in_child_A, const std::vector<int>& two_links_start_in_child_C,
         std::vector<unsigned int>& two_links_children_eq, std::vector<unsigned int>& two_links_children_ineq,
         unsigned int& two_links_root_eq, unsigned int& two_links_root_ineq ) const;
   void adjustActiveMylBy( int adjustment );
   void adjustActiveMzlBy( int adjustment );
   void adjustSizesAfterSplit( const std::vector<unsigned int>& two_links_children_eq,
         const std::vector<unsigned int>& two_links_children_ineq );

   void splitTree( int n_layers_total, int n_layers_left, sData* data ) override;

   sTree* shaveDenseBorder( int nx_to_shave, int myl_to_shave, int mzl_to_shave) override;

   /* inactive sizes store the original state of the tree when switching to the presolved data */
   long long N_INACTIVE{-1};
   long long MY_INACTIVE{-1};
   long long MZ_INACTIVE{-1};
   long long MYL_INACTIVE{-1};
   long long MZL_INACTIVE{-1};


   int nx_active{0};
   int my_active{0};
   int mz_active{0};
   int myl_active{0};
   int mzl_active{0};

   int nx_inactive{-1};
   int my_inactive{-1};
   int mz_inactive{-1};
   int myl_inactive{-1};
   int mzl_inactive{-1};

   bool isDataPresolved{false};
   bool hasPresolvedData{false};
   bool print_tree_sizes_on_reading{false};

   /* after this node has been split this will indicate how the children were assigned to the (new) sub_roots */
   std::vector<unsigned int> map_node_sub_root;

   sTreeCallbacks();
   InputNode* data{}; //input data
};


#endif
