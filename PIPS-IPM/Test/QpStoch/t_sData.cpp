#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "mpi.h"
#include "DistributedQP.hpp"
#include <tuple>

class HierarchicalDataSplittingTest
      : public DistributedQP,
        public ::testing::TestWithParam<std::tuple<std::vector<unsigned int>, std::vector<int>, int, std::vector<unsigned int>>> {
};

TEST_P(HierarchicalDataSplittingTest, TestPermutationOfLinkingConstraintsForSplit) {
   const std::vector<unsigned int> map_block_subtree = std::get<0>(GetParam());
   const std::vector<int> link_start_block_id = std::get<1>(GetParam());
   const int n_links_after_split = std::get<2>(GetParam());
   const std::vector<unsigned int>& expected_permutation = std::get<3>(GetParam());

   PERMUTATION result = getChildLinkConsFirstOwnLinkConsLastPermutation(map_block_subtree, link_start_block_id, n_links_after_split);

   EXPECT_EQ(result.size(), expected_permutation.size());
   EXPECT_THAT(result, ::testing::ContainerEq(expected_permutation));
}

INSTANTIATE_TEST_CASE_P
(SplitPermutationTests, HierarchicalDataSplittingTest, ::testing::Values(
      std::make_tuple(std::vector<unsigned int>{0, 0, 1, 1}, std::vector<int>{0, 0, 0, 1, 1, 2, 2, 2, 2, -1, -1, -1, -1}, 6,
            std::vector<unsigned int>{0, 1, 2, 5, 6, 7, 8, 3, 4, 9, 10, 11, 12}),
      std::make_tuple(std::vector<unsigned int>{0, 0, 1, 1, 1, 2, 2, 3, 3, 3},
            std::vector<int>{0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 4, 4, 5, 5, 5, 6, 7, 8, 8, 8, -1, -1, -1}, 10,
            std::vector<unsigned int>{0, 1, 2, 7, 8, 9, 12, 13, 14, 16, 17, 18, 19, 3, 4, 5, 6, 10, 11, 15, 20, 21, 22}),
      std::make_tuple(std::vector<unsigned int>{0, 0, 0, 0, 1, 2, 2, 2}, std::vector<int>{0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6}, 4,
            std::vector<unsigned int>{0, 1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 6, 7, 8, 9}),
      std::make_tuple(std::vector<unsigned int>{0, 0, 0, 1, 1, 2, 2, 2}, std::vector<int>{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}, 14,
            std::vector<unsigned int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13})));
