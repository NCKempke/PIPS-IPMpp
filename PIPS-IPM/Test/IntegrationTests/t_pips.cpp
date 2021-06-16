/*
 * t_pips.cpp
 *
 *  Created on: 19.03.2021
 *      Author: bzfkempk
 */
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "../Verbosity.hpp"

#include "InteriorPointMethod.hpp"
#include "PIPSIPMppInterface.hpp"
#include "gmspips_reader.hpp"
#include "PIPSIPMppOptions.h"
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

class ScenarioTests : public ::testing::TestWithParam<Instance> {
protected:
   const std::string root{__ROOT_DIR__};
   const int world_size{PIPS_MPIgetSize()};

   void SetUp() override {
      ASSERT_GE(world_size, 1);
   }

   void TearDown() override {
   }

public:

   ScenarioTests() {
      if (const char* gams_env = std::getenv("GAMSSYSDIR")) {
         gams_path = std::string(gams_env);
      } else {
        std::cout << "For this test suite to run please set your environment variable GAMSSYSDIR pointing to the gams directory\n";
      }
   }

   std::string gams_path;

   void solveInstanceAndCheckResult(double expected_objective, int expected_iterations, const std::string& path, size_t n_blocks,
      PresolverType presolver_type, ScalerType scaler_type, MehrotraStrategyType primal_dual_type);

   [[nodiscard]] std::tuple<TerminationStatus, double, int, std::string>
   solveInstance(const std::string& path_instance, size_t n_blocks, PresolverType presolver, ScalerType scaler,
      MehrotraStrategyType primal_dual_type) const;
};

std::vector<Instance> getInstances() {
   const std::string root{__ROOT_DIR__};
   const std::string file_name{root + "/PIPS-IPM/Test/IntegrationTests/gamssmall_instance_data.txt"};
   return readInstancesFromFile(file_name);
}


std::tuple<TerminationStatus, double, int, std::string>
ScenarioTests::solveInstance(const std::string& path_instance, size_t n_blocks, PresolverType presolver,
   ScalerType scaler, MehrotraStrategyType primal_dual_type) const {

   if (!verbose) {
      testing::internal::CaptureStdout();
   }

   gmspips_reader reader(path_instance, gams_path, n_blocks);
   std::unique_ptr<DistributedInputTree> tree(reader.read_problem());

   pipsipmpp_options::set_bool_parameter("GONDZIO_ADAPTIVE_LINESEARCH", false);

   double objective = std::numeric_limits<double>::infinity();
   int n_iterations = -1;

   PIPSIPMppInterface pipsIpm(tree.get(), primal_dual_type, MPI_COMM_WORLD, scaler, presolver);

   TerminationStatus result = TerminationStatus::DID_NOT_RUN;
   try {
      result = pipsIpm.run();
      objective = pipsIpm.getObjective();
      n_iterations = pipsIpm.n_iterations();
   }
   catch (...) {
      EXPECT_TRUE(false) << " PIPS threw while solving " << path_instance;
   }

   if (!verbose) {
      const std::string output = testing::internal::GetCapturedStdout();
      return {result, objective, n_iterations, output};
   } else {
      const std::string output = "";
      return {result, objective, n_iterations, output};
   }
};

void ScenarioTests::solveInstanceAndCheckResult(double expected_objective, int expected_iterations, const std::string& path, size_t n_blocks,
   PresolverType presolver_type, ScalerType scaler_type, MehrotraStrategyType primal_dual_type) {

   ASSERT_GE(world_size, 1);

   if (static_cast<size_t>(world_size) >= n_blocks)
      GTEST_SKIP();

   const auto[result, objective_solve, n_iterations, output_solve] = solveInstance(path, n_blocks, presolver_type, scaler_type, primal_dual_type);

   EXPECT_EQ(result, TerminationStatus::SUCCESSFUL_TERMINATION);
   EXPECT_NEAR(expected_objective, objective_solve, solution_tol) << " while solving " << path << "\nOutput_run: " << output_solve << "\n";

   /* allow a 10% increase in the number of iterations */
   EXPECT_LE(n_iterations, std::ceil(expected_iterations * 1.1)) << " solving took too may iterations - performance might be affected\n";
}

TEST_P(ScenarioTests, TestGamssmallPrimalDualStepScaleGeo) {
   const std::string& problem_paths(GetParam().name);
   const size_t n_blocks(GetParam().n_blocks);
   const double result(GetParam().result);

   const int n_expected_iterations(GetParam().n_iterations);

   solveInstanceAndCheckResult(result, n_expected_iterations, root + problem_paths, n_blocks, PresolverType::PRESOLVER_NONE, ScalerType::SCALER_GEO_STOCH, MehrotraStrategyType::PRIMAL_DUAL);
};

TEST_P(ScenarioTests, TestGamssmallPrimalDualStepScaleGeoPresolve) {
   const std::string& problem_paths(GetParam().name);
   const size_t n_blocks(GetParam().n_blocks);
   const double result(GetParam().result);
   const int n_expected_iterations(GetParam().n_iterations);

   solveInstanceAndCheckResult(result, n_expected_iterations, root + problem_paths, n_blocks, PresolverType::PRESOLVER_STOCH, ScalerType::SCALER_GEO_STOCH, MehrotraStrategyType::PRIMAL_DUAL);
};

TEST_P(ScenarioTests, TestGamssmallPrimalDualStep) {
   const std::string& problem_paths(GetParam().name);
   const size_t n_blocks(GetParam().n_blocks);
   const double result(GetParam().result);
   const int n_expected_iterations(GetParam().n_iterations);

   solveInstanceAndCheckResult(result, n_expected_iterations, root + problem_paths, n_blocks, PresolverType::PRESOLVER_NONE, ScalerType::SCALER_NONE, MehrotraStrategyType::PRIMAL_DUAL);
};

TEST_P(ScenarioTests, TestGamssmallPrimalDualStepPresolve) {
   const std::string& problem_paths(GetParam().name);
   const size_t n_blocks(GetParam().n_blocks);
   const double result(GetParam().result);
   const int n_expected_iterations(GetParam().n_iterations);

   solveInstanceAndCheckResult(result, n_expected_iterations, root + problem_paths, n_blocks, PresolverType::PRESOLVER_STOCH, ScalerType::SCALER_NONE, MehrotraStrategyType::PRIMAL_DUAL);
};

TEST_P(ScenarioTests, TestGamssmallNoSettings) {
   const std::string& problem_paths(GetParam().name);
   const size_t n_blocks(GetParam().n_blocks);
   const double result(GetParam().result);
   const int n_expected_iterations(GetParam().n_iterations);

   solveInstanceAndCheckResult(result, n_expected_iterations, root + problem_paths, n_blocks, PresolverType::PRESOLVER_NONE, ScalerType::SCALER_NONE, MehrotraStrategyType::PRIMAL);
};

TEST_P(ScenarioTests, TestGamssmallPresolve) {
   const std::string& problem_paths(GetParam().name);
   const size_t n_blocks(GetParam().n_blocks);
   const double result(GetParam().result);
   const int n_expected_iterations(GetParam().n_iterations);

   solveInstanceAndCheckResult(result, n_expected_iterations, root + problem_paths, n_blocks, PresolverType::PRESOLVER_STOCH, ScalerType::SCALER_NONE, MehrotraStrategyType::PRIMAL);
};

TEST_P(ScenarioTests, TestGamssmallScaleGeoPresolve) {
   const std::string& problem_paths(GetParam().name);
   const size_t n_blocks(GetParam().n_blocks);
   const double result(GetParam().result);
   const int n_expected_iterations(GetParam().n_iterations);

   solveInstanceAndCheckResult(result, n_expected_iterations, root + problem_paths, n_blocks, PresolverType::PRESOLVER_STOCH, ScalerType::SCALER_GEO_STOCH, MehrotraStrategyType::PRIMAL);
};

INSTANTIATE_TEST_SUITE_P(InstantiateTestsWithAllGamssmallInstances, ScenarioTests, ::testing::ValuesIn(getInstances()));
