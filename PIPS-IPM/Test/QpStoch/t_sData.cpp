#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "mpi.h"
#include "sData.h"

class HierarchicalMappingParametersTest : public sData, public ::testing::TestWithParam<std::vector<unsigned int>>
{};

TEST_F(HierarchicalMappingParametersTest, DummyEquals)
{
   EXPECT_EQ(1, 1);
}
