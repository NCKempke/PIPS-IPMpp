/*
 * t_gmspips_reader.cpp
 *
 *  Created on: 19.03.2021
 *      Author: bzfkempk
 */

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "gmspips_reader.hpp"
#include "DistributedFactory.hpp"
#include <cstdlib>
#include <tuple>
#include <memory>

#ifndef __ROOT_DIR__
#define __ROOT_DIR__ "define_this_variable_for_the_tests_to_succeed"
#endif

class GmspipsReaderTest : public gmspips_reader, public ::testing::TestWithParam<std::tuple<const std::string, size_t>> {
public:
   std::string gams_path;

   GmspipsReaderTest() {
      if (const char* gams_env = std::getenv("GAMSSYSDIR")) {
         if (!gams_env)
            std::cout << "For this test suite to run please set your environment variable GAMSSYSDIR pointing to the gams directory\n";
         gams_path = std::string(gams_env);
      }
   }

   ~GmspipsReaderTest() override = default;
};


TEST_P(GmspipsReaderTest, TestReadGamsSmallProblems) {
   const std::string& problem_paths(std::get<0>(GetParam()));
   const size_t& n_blocks_per_problem(std::get<1>(GetParam()));

   gmspips_reader reader(problem_paths, gams_path, n_blocks_per_problem);

   std::unique_ptr<DistributedInputTree> tree(reader.read_problem());
   ASSERT_TRUE(tree);
   DistributedFactory factory(tree.get(), MPI_COMM_WORLD);
};

// TODO ... this is hardcoded and needs to change...
const std::string root{__ROOT_DIR__};
const std::string examples_path{root + "/PIPS-IPM/Drivers/simple/GAMSsmall/"};
INSTANTIATE_TEST_CASE_P
(TestReadGamsSmallProblems, GmspipsReaderTest,
      ::testing::Values(std::make_tuple(examples_path + "examples_boundTightening/run_exampleAC_boundStrength/exampleAC_boundStrength", 3),
//            std::make_tuple( examples_path + "examples_dualfixing/run_dual_fixing_A2/exampleAC_boundStrength", 3 ),
            std::make_tuple(examples_path + "examples_hierarchical_approach/run_hier_approach_2blocks_2by2/hier_approach_2blocks_2by2", 4),
            std::make_tuple(examples_path + "examples_hierarchical_approach/run_hier_approach_4blocks_2by2/hier_approach_4blocks_2by2", 5),
            std::make_tuple(examples_path + "examples_hierarchical_approach/run_hier_approach_4blocks_2by3/hier_approach_4blocks_2by3", 5),
            std::make_tuple(examples_path + "examples_hierarchical_approach/run_hier_approach_8blocks_2by3/hier_approach_8blocks_2by3", 9),

            std::make_tuple(examples_path +
                            "examples_nearlyParallelRows/run_nearlyParallelEqualityAndInequalityRows_B0A2/nearlyParallelEqualityAndInequalityRows_B0A2",
                  4),
//            std::make_tuple( examples_path + "examples_nearlyParallelRows/run_nearlyParallelEqualityInequalityRowsMixed/nearlyParallelEqualityInequalityRowsMixed", 8 ),
            std::make_tuple(examples_path +
                            "examples_nearlyParallelRows/run_nearlyParallelEqualityRowsBothSingletons_B0A2/nearlyParallelEqualityRowsBothSingletons_B0A2",
                  4), std::make_tuple(examples_path +
                                      "examples_nearlyParallelRows/run_nearlyParallelEqualityRowsOneRowNoSingleton_B0A2/nearlyParallelEqualityRowsOneRowNoSingleton_B0A2",
                  4),
            std::make_tuple(examples_path + "examples_nearlyParallelRows/run_nearlyParallelInequalityRows_B0A2/nearlyParallelInequalityRows_B0A2", 4),

            std::make_tuple(examples_path + "examples_parallelRows/run_parallelEqualityAndInequalityRow_B0A2/parallelEqualityAndInequalityRow_B0A2",
                  4),
            std::make_tuple(examples_path + "examples_parallelRows/run_parallelEqualityInequalityRowsMixed/parallelEqualityInequalityRowsMixed", 4),
            std::make_tuple(examples_path + "examples_parallelRows/run_parallelEqualityRows_B0A2/parallelEqualityRows_B0A2", 4),

            std::make_tuple(examples_path + "examples_parallelRows/run_parallelInequalityRows_B0A2/parallelInequalityRows_B0A2", 4),

            std::make_tuple(examples_path + "examples_singletonEqualityColumn/run_singletonEqualityColumn_A2/singletonEqualityColumn_A2", 4),
            std::make_tuple(examples_path + "examples_singletonEqualityColumn/run_singletonEqualityColumn_B0B1A2/singletonEqualityColumn_B0B1A2", 4),
            std::make_tuple(examples_path + "examples_singletonEqualityColumn/run_singletonEqualityColumn_B0Bl0/singletonEqualityColumn_B0Bl0", 4),
            std::make_tuple(examples_path + "examples_singletonEqualityColumn/run_singletonEqualityColumn_B0Bl2/singletonEqualityColumn_B0Bl2", 4),
            std::make_tuple(examples_path + "examples_singletonEqualityColumn/run_singletonEqualityColumn_B0/singletonEqualityColumn_B0", 4),
            std::make_tuple(examples_path + "examples_singletonEqualityColumn/run_singletonEqualityColumn_B1/singletonEqualityColumn_B1", 4),
            std::make_tuple(examples_path +
                            "examples_singletonEqualityColumn/run_singletonEqualityColumn_multiple_noLink/singletonEqualityColumn_multiple_noLink",
                  4), std::make_tuple(examples_path +
                                      "examples_singletonEqualityColumn/run_singletonEqualityColumn_multiple_resulting/singletonEqualityColumn_multiple_resulting",
                  4), std::make_tuple(examples_path +
                                      "examples_singletonEqualityColumn/run_singletonEqualityColumn_multiple_resulting_noLink/singletonEqualityColumn_multiple_resulting_noLink",
                  4),

            std::make_tuple(examples_path + "examples_singletonInequalityColumn/run_singletonInequalityColumn_A2/singletonInequalityColumn_A2", 4),
            std::make_tuple(
                  examples_path + "examples_singletonInequalityColumn/run_singletonInequalityColumn_B0B2A1/singletonInequalityColumn_B0B2A1", 4),
            std::make_tuple(examples_path + "examples_singletonInequalityColumn/run_singletonInequalityColumn_B0Bl0/singletonInequalityColumn_B0Bl0",
                  4),
            std::make_tuple(examples_path + "examples_singletonInequalityColumn/run_singletonInequalityColumn_B0Bl2/singletonInequalityColumn_B0Bl2",
                  4),
            std::make_tuple(examples_path + "examples_singletonInequalityColumn/run_singletonInequalityColumn_B0/singletonInequalityColumn_B0", 4),
            std::make_tuple(examples_path + "examples_singletonInequalityColumn/run_singletonInequalityColumn_B1/singletonInequalityColumn_B1", 4),
            std::make_tuple(examples_path +
                            "examples_singletonInequalityColumn/run_singletonInequalityColumn_multipleLinkingFixationsAndSingletonCols/singletonInequalityColumn_multipleLinkingFixationsAndSingletonCols",
                  4), std::make_tuple(examples_path +
                                      "examples_singletonInequalityColumn/run_singletonInequalityColumn_multiple_resulting/singletonInequalityColumn_multiple_resulting",
                  4),

            std::make_tuple(examples_path + "examples_singletonRow/run_exampleAC_singletonRow2/exampleAC_singletonRow2", 3), std::make_tuple(
                  examples_path + "examples_singletonRow/run_exampleAC_singletonRow3_singletonLinkingRow/exampleAC_singletonRow3_singletonLinkingRow",
                  3), std::make_tuple(examples_path + "examples_singletonRow/run_exampleAC_singletonRow/exampleAC_singletonRow", 3),
            std::make_tuple(examples_path + "examples_singletonRow/run_example_breakSingletonRows/example_breakSingletonRows", 4)));
