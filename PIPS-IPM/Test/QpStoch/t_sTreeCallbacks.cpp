#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "mpi.h"
#include "sTreeCallbacks.h"

class HierarchicalMappingParametersTest : public sTreeCallbacks, public ::testing::TestWithParam<std::vector<unsigned int>>
{};

// Tests factorial of 0.
TEST_F(HierarchicalMappingParametersTest, DummyEquals)
{
   EXPECT_EQ(1, 1);
}

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

class HierarchicalSplittingTest : public sTreeCallbacks, public ::testing::Test
{
      void SetUp() override
      {
         sTree::numProcs = 1;
         sTree::rankPrcnd = -1;
         sTree::rankZeroW = 0;
         sTree::rankMe = 0;
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

TEST_F(HierarchicalSplittingTest, CorrectTreeSplitAndSizeAdjustment)
{
   sTreeCallbacks* test_tree = createTestTree(10, 20, 30);
   test_tree->assertTreeStructureCorrect();

   // say every block has 2 two_links to the next
   std::vector<int> twoLinksStartBlockA(10, 2);
   twoLinksStartBlockA[9] = 0;

   std::vector<int> twoLinksStartBlockC(10, 3);
   twoLinksStartBlockA[9] = 0;

   sTreeCallbacks* test_tree_split = dynamic_cast<sTreeCallbacks*>(test_tree->switchToHierarchicalTree(0, 0, 0, twoLinksStartBlockA, twoLinksStartBlockC ));

   delete test_tree_split;
}

