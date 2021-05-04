message(STATUS "Running google test")

add_subdirectory("${PROJECT_SOURCE_DIR}/ThirdPartyLibs/googletest" "ThirdPartyLibs/googltest")

mark_as_advanced(
        BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
        gmock_build_tests gtest_build_samples gtest_build_tests
        gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

add_subdirectory(Utilities)
include_directories(Utilities)

add_library(customMain TestsMain.cpp)
target_link_libraries(customMain MPI::MPI_CXX gtest gmock)

macro(package_add_test TESTNAME)
    # create an exectuable in which the tests will be stored
    add_executable(${TESTNAME} ${ARGN})
    # link the Google test infrastructure, mocking library, and a default main fuction to
    # the test executable.  Remove g_test_main if writing your own main function.
    target_link_libraries(${TESTNAME}
            PRIVATE gmspips_reader
            PRIVATE ooqpstoch
            PRIVATE ooqpstochla
            PRIVATE interior_point_method
            PRIVATE ooqpequistochscaler
            PRIVATE ooqpgensparse
            PRIVATE ooqpsparse
            PRIVATE ooqpdense
            PRIVATE ooqpbase
            PRIVATE ooqpqpstochpresolver
            PRIVATE OpenMP::OpenMP_CXX
            PRIVATE customMain
            PRIVATE utilities
            )
    # the root PIPS directory (used for finding test related data)
    target_compile_definitions(${TESTNAME} PRIVATE -D__ROOT_DIR__="${CMAKE_SOURCE_DIR}")

    # gtest_discover_tests replaces gtest_add_tests,
    # see https://cmake.org/cmake/help/v3.10/module/GoogleTest.html for more options to pass to it
    gtest_discover_tests(${TESTNAME}
            # set a working directory so your project root so that you can find test data via paths relative to the project root
            WORKING_DIRECTORY ${PROJECT_DIR}
            DISCOVERY_TIMEOUT 10
            )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

add_subdirectory(QpStoch)
add_subdirectory(StochLinearAlgebra)
add_subdirectory(Preprocessing)
add_subdirectory(Drivers)
add_subdirectory(IntegrationTests)
