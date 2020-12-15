#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "sTreeCallbacks.h"

class HierarchicalMappingParametersTest : public sTreeCallbacks, public ::testing::TestWithParam<std::vector<unsigned int>>
{
};

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

