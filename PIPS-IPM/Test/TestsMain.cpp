#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "mpi.h"

class MPIDummyTestingEnvironment : public testing::Environment {
 public:
  ~MPIDummyTestingEnvironment() override = default;

  void SetUp() override
  {
     char** argv = nullptr;
     int argc = 0;
     const int mpiError = MPI_Init(&argc, &argv);
     ASSERT_FALSE(mpiError);
  }

  void TearDown() override
  {
     const int mpiError = MPI_Finalize();
     ASSERT_FALSE(mpiError);
  }
};


int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
//  testing::FLAGS_gtest_death_test_style = "fast";

  testing::AddGlobalTestEnvironment(new MPIDummyTestingEnvironment);

  return RUN_ALL_TESTS();
}
