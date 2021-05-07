#include "gtest/gtest.h"

#include "DistributedTreeCallbacks.h"
#include "StripMatrix.h"
#include "mpi.h"
#include "DistributedMatrix.h"

#include <tuple>
#include <numeric>

class DistributedMatrixSplittingTest
      : public DistributedMatrix, public ::testing::TestWithParam<std::tuple<std::vector<unsigned int>, std::vector<int>, std::vector<unsigned int>>> {
public:
   DistributedMatrix*
   createTestMatrix(unsigned int n_link_vars, unsigned int m_link_vars, unsigned int n_diag, unsigned int m_diag, unsigned int n_blocks,
         unsigned int m_link_cons) const;
};

DistributedMatrix*
DistributedMatrixSplittingTest::createTestMatrix(unsigned int n_link_vars, unsigned int m_link_vars, unsigned int n_diag, unsigned int m_diag,
      unsigned int n_blocks, unsigned int m_link_cons) const {
   DistributedMatrix* root = new DistributedMatrix(0, 0, 0, 0, 0, m_link_vars, n_link_vars, 0, m_link_cons, n_link_vars, 0, MPI_COMM_WORLD);

   for (unsigned int i = 0; i < n_blocks; ++i) {
      DistributedMatrix* child = new DistributedMatrix(0, 0, m_diag, n_link_vars, 0, m_diag, n_diag, 0, m_link_cons, n_diag, 0, MPI_COMM_WORLD);
      root->AddChild(child);
   }

   root->recomputeSize();
   return root;
}

TEST_P(DistributedMatrixSplittingTest, TestSplitDistributedMatrixOnce) {
   const std::vector<unsigned int>& mat_input(std::get<0>(GetParam()));
   const std::vector<int>& twolinks_start_in_block_id(std::get<1>(GetParam()));
   const std::vector<unsigned int>& map_blocks_children(std::get<2>(GetParam()));

   ASSERT_EQ(mat_input.size(), 6);
   const int n_links_vars = mat_input[0];
   const int m_links_vars = mat_input[1];
   const int n_diag = mat_input[2];
   const int m_diag = mat_input[3];
   const int n_blocks = mat_input[4];
   const int m_links_cons = mat_input[5];

   const int n_children = getNDistinctValues(map_blocks_children);
   std::vector<MPI_Comm> child_comms(n_children, MPI_COMM_SELF);

   int n_links_in_root{0};
   for (unsigned int i = 0; i < map_blocks_children.size() - 1; ++i) {
      if (map_blocks_children[i] != map_blocks_children[i + 1])
         n_links_in_root += twolinks_start_in_block_id[i];
   }
   const int sum_twolinks = std::accumulate(twolinks_start_in_block_id.begin(), twolinks_start_in_block_id.end(), 0);

   n_links_in_root += (m_links_cons - sum_twolinks);
   ASSERT_LE(sum_twolinks, m_links_cons);
   ASSERT_LE(n_links_in_root, m_links_cons);
   ASSERT_EQ(twolinks_start_in_block_id.size(), mat_input[4]);

   std::unique_ptr<DistributedMatrix> test_mat(createTestMatrix(n_links_vars, m_links_vars, n_diag, m_diag, n_blocks, m_links_cons));
   test_mat->splitMatrix(twolinks_start_in_block_id, map_blocks_children, n_links_in_root, child_comms);

   EXPECT_EQ(test_mat->children.size(), n_children);

   EXPECT_TRUE(test_mat->Amat->is_a(kSparseGenMatrix));
   EXPECT_TRUE(test_mat->Blmat->is_a(kSparseGenMatrix));

   unsigned int sum_child_children{0};
   unsigned int sum_child_Bmat_children{0};
   unsigned int sum_child_links{0};
   for (auto& child : test_mat->children) {
      EXPECT_TRUE(child->Bmat->is_a(kDistributedMatrix));
      EXPECT_TRUE(child->Blmat->is_a(kStripMatrix));
      EXPECT_TRUE(child->Amat->is_a(kSparseGenMatrix));

      int ma, na;
      child->Amat->getSize(ma, na);
      int mb, nb;
      child->Bmat->getSize(mb, nb);
      int mbl, nbl;
      child->Blmat->getSize(mbl, nbl);

      EXPECT_EQ(na, 0);
      EXPECT_EQ(ma, 0);
      EXPECT_EQ(nb, nbl);

      EXPECT_EQ(mbl, n_links_in_root);

      const DistributedMatrix& bmat = dynamic_cast<const DistributedMatrix&>(*child->Bmat);
      const StripMatrix& blmat = dynamic_cast<const StripMatrix&>(*child->Blmat);
      EXPECT_EQ(blmat.children.size(), bmat.children.size());
      EXPECT_EQ(0, child->children.size());

      bmat.Blmat->getSize(mbl, nbl);

      sum_child_links += mbl;
      sum_child_children += child->children.size();
      sum_child_Bmat_children += bmat.children.size();
   }

   EXPECT_EQ(n_links_in_root + sum_child_links, m_links_cons);
   EXPECT_EQ(sum_child_children, 0);
   EXPECT_EQ(sum_child_Bmat_children, map_blocks_children.size());
};

INSTANTIATE_TEST_CASE_P
(TestSplitDistributedMatrixOnce, DistributedMatrixSplittingTest, ::testing::Values(
      std::make_tuple(std::vector<unsigned int>{5, 5, 10, 3, 10, 25}, std::vector<int>{2, 2, 2, 2, 2, 2, 2, 2, 2, 0},
            std::vector<unsigned int>{0, 0, 0, 1, 1, 2, 2, 2, 2, 3}),
      std::make_tuple(std::vector<unsigned int>{10, 2, 2, 2, 10, 18}, std::vector<int>{2, 2, 2, 2, 2, 2, 2, 2, 2, 0},
            std::vector<unsigned int>{0, 0, 0, 1, 1, 2, 2, 2, 2, 3}),
      std::make_tuple(std::vector<unsigned int>{4, 6, 1, 5, 10, 25}, std::vector<int>{2, 1, 0, 0, 3, 2, 2, 2, 2, 0},
            std::vector<unsigned int>{0, 0, 0, 1, 1, 2, 2, 2, 3, 3}),
      std::make_tuple(std::vector<unsigned int>{5, 5, 10, 3, 10, 25}, std::vector<int>{2, 3, 4, 2, 2, 1, 0, 3, 3, 0},
            std::vector<unsigned int>{0, 1, 2, 3, 3, 4, 5, 5, 5, 5})));
