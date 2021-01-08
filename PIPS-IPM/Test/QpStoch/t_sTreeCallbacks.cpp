#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "mpi.h"
#include "sTreeCallbacks.h"

class HierarchicalMappingParametersTest : public sTreeCallbacks, public ::testing::TestWithParam<std::vector<unsigned int>>
{};

TEST_P(HierarchicalMappingParametersTest, CorrectMappingChildrenToNewRoots )
{
   const std::vector<unsigned int>& param = GetParam();

   std::vector<unsigned int> result;
   const unsigned int n_subtrees = getMapChildrenToSqrtNSubTrees( result, param.size() );

   EXPECT_EQ( n_subtrees, std::round( std::sqrt( param.size() ) ) );
   EXPECT_EQ( param.size(), result.size() );
   EXPECT_THAT( result, ::testing::ContainerEq(param));
}

const std::vector<std::vector<unsigned int>> maps_expected {
   {},
   { 0 },
   { 0, 0},
   { 0, 0, 1},
   { 0, 0, 1, 1 },
   { 0, 0, 0, 1, 1 },
   { 0, 0, 0, 1, 1, 1 },
   { 0, 0, 0, 1, 1, 2, 2},
   { 0, 0, 0, 1, 1, 2, 2, 2},
   { 0, 0, 0, 1, 1, 1, 2, 2, 2},
   { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2},
   { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2},
   { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3},
   { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3},
   { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3},
   { 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3},
   { 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3}
};

INSTANTIATE_TEST_SUITE_P(
      HierarchicalMappingTests,
      HierarchicalMappingParametersTest,
      testing::ValuesIn(maps_expected)
);

class HierarchicalSplittingTest : public sTreeCallbacks, public ::testing::TestWithParam<std::vector<unsigned int>>
{
      void SetUp() override
      {
         sTree::numProcs = 1;
         sTree::rankPrcnd = -1;
         sTree::rankZeroW = 0;
         sTree::rankMe = 0;
         pips_options::setBoolParameter("SILENT", true);
      }

      void TearDown() override
      {
         delete input_tree;
      }

   protected:
      StochInputTree::StochInputNode* root_node{nullptr};
      StochInputTree* input_tree{nullptr};
   public:
      sTreeCallbacks* createTestTree( int n_children, int n_eq_links, int n_ineq_links );

};

sTreeCallbacks* HierarchicalSplittingTest::createTestTree( int nChildren, int n_eq_linkings, int n_ineq_linkings )
{
   const int NX_ROOT = 10;
   const int MY_ROOT = 20;
   const int MZ_ROOT = 30;

   StochInputTree::StochInputNode* root_node = new StochInputTree::StochInputNode( -1, NX_ROOT, MY_ROOT, n_eq_linkings, MZ_ROOT, n_ineq_linkings );

   input_tree = new StochInputTree(root_node);

   for( int i = 0; i < nChildren; ++i )
   {
      const int NX_CHILD = 10 + i;
      const int MY_CHILD = 20 + i;
      const int MZ_CHILD = 30 + i;
      StochInputTree::StochInputNode* leaf_node = new StochInputTree::StochInputNode( i, NX_CHILD, MY_CHILD, n_eq_linkings, MZ_CHILD, n_ineq_linkings );

      input_tree->AddChild( new StochInputTree(leaf_node) );
   }

   sTreeCallbacks* tree = new sTreeCallbacks(input_tree);
   tree->assignProcesses();
   tree->computeGlobalSizes();

   return tree;
}

const std::vector<std::vector<unsigned int>> paramsets_tree_splitting {
   { 2, 3, 10, 12, 1, 3, 1, 10},
   { 1, 1, 10, 12, 1, 1, 1, 15 },
   { 5, 0, 10, 12, 1, 1, 1, 15 },
   { 3, 2, 5, 7, 5, 0, 0, 7 },
   { 3, 2, 5, 7, 5, 7, 2, 5 },
   { 3, 2, 5, 7, 0, 7, 2, 5 }
};

TEST_P(HierarchicalSplittingTest, CorrectTreeSplitAndSizeAdjustment)
{
   std::vector<unsigned int> params = GetParam();
   ASSERT_EQ( params.size(), 8 );
   const int two_links_eq_per_block = params[0];
   const int two_links_ineq_per_block = params[1];

   const unsigned int global_eq_links = params[2];
   const unsigned int global_ineq_links = params[3];

   const int shave_n_eqs = params[4];
   const int shave_n_ineqs = params[5];
   const int shave_n_vars = params[6];

   const int n_children = params[7];
   const int n_eq_links = (n_children - 1) * two_links_eq_per_block + global_eq_links;
   const int n_ineq_links = (n_children - 1) * two_links_ineq_per_block + global_ineq_links;

   sTreeCallbacks* test_tree = createTestTree( n_children, n_eq_links, n_ineq_links);
   test_tree->assertTreeStructureCorrect();

   std::vector<int> twoLinksStartBlockA(n_children, two_links_eq_per_block);
   twoLinksStartBlockA[n_children - 1] = 0;

   std::vector<int> twoLinksStartBlockC(n_children, two_links_ineq_per_block);
   twoLinksStartBlockC[n_children - 1] = 0;

   sTreeCallbacks* test_tree_split =
         dynamic_cast<sTreeCallbacks*>( test_tree->switchToHierarchicalTree(shave_n_vars, shave_n_eqs, shave_n_ineqs, twoLinksStartBlockA, twoLinksStartBlockC ) );

   /// check root for right amount of linking constraints
   long long dummy, MYL, MZL;
   test_tree_split->getGlobalSizes(dummy, dummy, MYL, dummy, MZL);
   EXPECT_EQ( MYL, n_eq_links );
   EXPECT_EQ( MZL, n_ineq_links );

   EXPECT_EQ( test_tree_split->myl(), shave_n_eqs );
   EXPECT_EQ( test_tree_split->mzl(), shave_n_ineqs );

   EXPECT_EQ( test_tree_split->nChildren(), 1 );

   /// check inner root node
   const sTreeCallbacks* inner_root = dynamic_cast<const sTreeCallbacks*>(test_tree_split->getChildren()[0]);
   inner_root->getGlobalSizes(dummy, dummy, MYL, dummy, MZL);
   EXPECT_EQ( MYL, n_eq_links - shave_n_eqs );
   EXPECT_EQ( MZL, n_ineq_links - shave_n_ineqs );

   EXPECT_EQ( inner_root->myl(), global_eq_links - shave_n_eqs + two_links_eq_per_block * ( inner_root->nChildren() - 1) );
   EXPECT_EQ( inner_root->mzl(), global_ineq_links - shave_n_ineqs + two_links_ineq_per_block * (inner_root->nChildren() - 1) );

   int sum_children_exclusive_eq_links = 0;
   int sum_children_exclusive_ineq_links = 0;
   for( auto& child : inner_root->getChildren() )
   {
      EXPECT_EQ( child->myl(), inner_root->myl() );
      EXPECT_EQ( child->mzl(), inner_root->mzl() );

      long long MYL_child, MZL_child, MYL_subroot, MZL_subroot;
      child->getGlobalSizes(dummy, dummy, MYL_child, dummy, MZL_child);
      child->getSubRoot()->getGlobalSizes(dummy, dummy, MYL_subroot, dummy, MZL_subroot);

      EXPECT_EQ( MYL_subroot, child->getSubRoot()->myl() );
      EXPECT_EQ( MZL_subroot, child->getSubRoot()->mzl() );

      /// check sub-roots
      for( auto& childchild : child->getSubRoot()->getChildren() )
      {
         long long MYL_childchild, MZL_childchild;
         childchild->getGlobalSizes(dummy, dummy, MYL_childchild, dummy, MZL_childchild);
         EXPECT_EQ( MYL_childchild, childchild->myl() );
         EXPECT_EQ( MZL_childchild, childchild->mzl() );

         EXPECT_EQ( MYL_childchild, MYL_subroot );
         EXPECT_EQ( MZL_childchild, MZL_subroot );
      }

      EXPECT_EQ( MYL_subroot, MYL_child - child->myl() );
      EXPECT_EQ( MZL_subroot, MZL_child - child->mzl() );

      sum_children_exclusive_eq_links += MYL_child - child->myl();
      sum_children_exclusive_ineq_links += MZL_child - child->mzl();
   }

   EXPECT_EQ( sum_children_exclusive_eq_links + inner_root->myl(), MYL );
   EXPECT_EQ( sum_children_exclusive_ineq_links + inner_root->mzl(), MZL );

   delete test_tree_split;
}

INSTANTIATE_TEST_SUITE_P(
      CorrectTreeSplitAndSizeAdjustment,
      HierarchicalSplittingTest,
      testing::ValuesIn(paramsets_tree_splitting)
);
