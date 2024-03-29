add_library(pips-ipmpp STATIC
        Base/Observer.C

        Globalization/Filter.cpp
        Globalization/FilterLineSearch.cpp
        Globalization/FilterStrategy.cpp

        Interface/PIPSIPMppInterface.cpp

        Options/AbstractOptions.C
        Options/Options.C
        Options/PIPSIPMppOptions.C

        Problems/DistributedFactory.cpp
        Problems/DistributedProblem.cpp
        Problems/Problem.cpp
        Problems/StochResourcesMonitor.cpp

        Readers/Distributed/DistributedInputTree.C
        Readers/Distributed/DistributedTree.C
        Readers/Distributed/DistributedTreeCallbacks.C
        Readers/MpsReader.C

        Utilities/pipschecks.C
        Utilities/hash.C
        Utilities/sort.cpp
        )

target_include_directories(pips-ipmpp
    PRIVATE
        Base
        Globalization
        Interface
        Options
        Problems
        Readers
        Utilities
        ${BOOST_INCLUDE_DIRS}
    PUBLIC
        Readers/Distributed
    )

add_subdirectory(InteriorPointMethod)
add_subdirectory(KKTFormulation)
add_subdirectory(LinearAlgebra)

add_subdirectory(LinearSolvers)

MESSAGE(STATUS "Solver libraries: ${lin_solvers_list}")
foreach(linear_solver ${lin_solvers_list})
    target_link_libraries(pips-ipmpp PRIVATE ${linear_solver})
endforeach(linear_solver)

add_subdirectory(Preprocessing)

target_link_libraries(pips-ipmpp
        PRIVATE Boost::boost
        PRIVATE ${MATH_LIBS}
)
