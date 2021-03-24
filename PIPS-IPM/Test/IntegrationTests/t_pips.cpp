/*
 * t_pips.cpp
 *
 *  Created on: 19.03.2021
 *      Author: bzfkempk
 */
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "sFactoryAug.h"
#include "GondzioStochLpSolver.h"
#include "PIPSIpmInterface.h"
#include "gmspips_reader.hpp"

#include "utilities.hpp"

#include <memory>
#include <vector>
#include <tuple>
#include <string>

#ifndef __ROOT_DIR__
#define __ROOT_DIR__ "define_this_variable_for_the_tests_to_succeed"
#endif

// TODO : criterion should be mathematical and adaptive to scaling ..
const double solution_tol = 1e-4;

class ScenarioTests : public ::testing::TestWithParam<Instance>
{
   protected:
      const std::string root{ __ROOT_DIR__ };
      const int world_size{PIPS_MPIgetSize()};

      void SetUp() override
      {
         ASSERT_GE(world_size, 1);
      }

      void TearDown() override
      {
      }

   public:

      ScenarioTests()
      {
         if(const char* gams_env = std::getenv("GAMSSYSDIR"))
         {
            if( !gams_env )
               std::cout << "For this test suite to run please set your environment variable GAMSSYSDIR pointing to the gams directory\n";
            gams_path = std::string(gams_env);
         }
      }

      std::string gams_path;

      double solveInstance( const std::string& path_instance, size_t n_blocks, PresolverType presolver, ScalerType scaler );
};

std::vector<Instance> getInstances()
{
   const std::string root{ __ROOT_DIR__ };
   const std::string file_name{ root + "/PIPS-IPM/Test/IntegrationTests/gamssmall_instance_data.txt" };
   return readInstancesFromFile(file_name);
}


double ScenarioTests::solveInstance( const std::string& path_instance, size_t n_blocks, PresolverType presolver, ScalerType scaler )
{
   testing::internal::CaptureStdout();

   gmspips_reader reader( path_instance, gams_path, n_blocks );
   std::unique_ptr<StochInputTree> tree( reader.read_problem() );

   pips_options::setBoolParameter("GONDZIO_ADAPTIVE_LINESEARCH", false);

   PIPSIpmInterface<sFactoryAug, GondzioStochLpSolver> pipsIpm(tree.get(), MPI_COMM_WORLD,
         scaler, presolver );

   double result = std::numeric_limits<double>::infinity();
   try
   {
      pipsIpm.go();
      result = pipsIpm.getObjective();
   }
   catch( ... )
   {
      EXPECT_TRUE(false);
   }

   testing::internal::GetCapturedStdout();
   return result;
};

TEST_P( ScenarioTests, TestAllGamssmallScaleGeoPresolveDisabled )
{
   const std::string& problem_paths(GetParam().name);
   const size_t& n_blocks(GetParam().n_blocks);
   const double& result(GetParam().result);

   ASSERT_GE( world_size, 1 );

   if( static_cast<size_t>(world_size) >= n_blocks )
      GTEST_SKIP();

   const double result_solve = solveInstance( root + problem_paths, n_blocks, PRESOLVER_NONE, SCALER_GEO_STOCH );

   EXPECT_NEAR( result, result_solve, solution_tol ) << " while solving " << problem_paths;
};

TEST_P( ScenarioTests, TestAllGamssmallScaleGeoPresolveDefault )
{
   const std::string& problem_paths(GetParam().name);
   const size_t& n_blocks(GetParam().n_blocks);
   const double& result(GetParam().result);

   ASSERT_GE( world_size, 1 );

   /* skip test if too many mpi procs were run for this example */
   if( static_cast<size_t>(world_size) >= n_blocks )
      GTEST_SKIP();

   const double result_solve = solveInstance( root + problem_paths, n_blocks, PRESOLVER_STOCH, SCALER_GEO_STOCH );

   EXPECT_NEAR( result, result_solve, solution_tol ) << " while solving " << problem_paths;
};

TEST_P( ScenarioTests, TestAllGamssmallNoScalePresolveDisabled )
{
   const std::string& problem_paths(GetParam().name);
   const size_t& n_blocks(GetParam().n_blocks);
   const double& result(GetParam().result);

   ASSERT_GE( world_size, 1 );

   if( static_cast<size_t>(world_size) >= n_blocks )
      GTEST_SKIP();

   const double result_solve = solveInstance( root + problem_paths, n_blocks, PRESOLVER_NONE, SCALER_NONE );

   EXPECT_NEAR( result, result_solve, solution_tol ) << " while solving " << problem_paths;
};

TEST_P( ScenarioTests, TestAllGamssmallNoScalePresolveDefault )
{
   const std::string& problem_paths(GetParam().name);
   const size_t& n_blocks(GetParam().n_blocks);
   const double& result(GetParam().result);

   ASSERT_GE( world_size, 1 );

   if( static_cast<size_t>(world_size) >= n_blocks )
      GTEST_SKIP();

   const double result_solve = solveInstance( root + problem_paths, n_blocks, PRESOLVER_STOCH, SCALER_NONE );

   EXPECT_NEAR( result, result_solve, solution_tol ) << " while solving " << problem_paths;
};

INSTANTIATE_TEST_SUITE_P(
      InstantiateTestsWithAllGamssmallInstances,
      ScenarioTests,
      ::testing::ValuesIn(getInstances())
);
