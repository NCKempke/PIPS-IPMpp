#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "mpi.h"
#include "sData.h"
#include <tuple>

class HierarchicalDataSplittingTest : public sData, public ::testing::TestWithParam<std::tuple<std::vector<unsigned int>, std::vector<int>, int, std::vector<unsigned int>>>
{};

TEST_P(HierarchicalDataSplittingTest, TestPermutationOfLinkingConstraintsForSplit )
{
   PERMUTATION result = getChildLinkConsFirstOwnLinkConsLastPermutation( std::get<0>(GetParam()),
         std::get<1>(GetParam()), std::get<2>(GetParam()) );

   EXPECT_EQ( result.size(), std::get<3>(GetParam()).size() );
   EXPECT_THAT( result, ::testing::ContainerEq(std::get<3>(GetParam())) );
}

INSTANTIATE_TEST_CASE_P(
        SplitPermutationTests,
        HierarchicalDataSplittingTest,
        ::testing::Values(
              std::make_tuple(
                    std::vector<unsigned int>{0,0,1,1},
                    std::vector<int>{0,0,0,1,1,2,2,2,2,-1,-1,-1,-1},
                    6,
                    std::vector<unsigned int>{0,1,2,7,8,3,4,5,6,9,10,11,12}
              ),
              std::make_tuple(
                    std::vector<unsigned int>{0,0,1,1,1,2,2,3,3,3},
                    std::vector<int>{0,0,0,1,1,1,1,2,2,2,4,4,5,5,5,6,7,8,8,8,-1,-1,-1},
                    10,
                    std::vector<unsigned int>{0,1,2,13,14,15,16,3,4,5,17,18,6,7,8,19,9,10,11,12,20,21,22}
              ),
              std::make_tuple(
                    std::vector<unsigned int>{0,0,0,0,1,2,2,2},
                    std::vector<int>{0,0,1,1,2,2,3,3,4,4,5,5,5,6,6},
                    4,
                    std::vector<unsigned int>{0,1,2,3,4,5,11,12,13,14,6,7,8,9,10}
              )
        )
);
