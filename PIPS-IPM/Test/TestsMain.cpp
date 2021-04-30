#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <memory>

#include "mpi.h"
#include "utilities.hpp"

class MPITestingEnvironment : public testing::Environment {
public:
   ~MPITestingEnvironment() override = default;

   void SetUp() override {
   }

   void TearDown() override {
   }
};


int main(int argc, char** argv) {
#ifndef NDEBUG
   const int mpi_init_error = MPI_Init(&argc, &argv);
   assert(!mpi_init_error);
#else
   MPI_Init(&argc, &argv);
#endif

   testing::InitGoogleTest(&argc, argv);
//  testing::FLAGS_gtest_death_test_style = "fast";

   /* supress all test output except for rank 0 */
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   int size;
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   if (myrank != 0) {
      ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
      delete listeners.Release(listeners.default_result_printer());
   }
   else
      std::cout << "Running tests with " << size << " MPI processes\n";

   testing::AddGlobalTestEnvironment(new MPITestingEnvironment);

   const int test_result = RUN_ALL_TESTS();

   if (myrank == 0)
      std::cout << "Google Test exited with " << test_result << "\n";

#ifndef NDEBUG
   const int mpi_finalize_error = MPI_Finalize();
   assert(!mpi_finalize_error);
#else
   MPI_Finalize();
#endif

   return 0;
}
