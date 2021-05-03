add_subdirectory(InteriorPointMethod)
add_subdirectory(QpPreprocess)
add_subdirectory(LinearAlgebra/Sparse)
add_subdirectory(LinearAlgebra/Dense)
add_subdirectory(LinearAlgebra/Distributed)
add_subdirectory(LinearSolvers)

add_library(ooqpbase STATIC Abstract/OoqpVersion.C Abstract/Options.C Abstract/Observer.C Abstract/Variables.C Problems/Problem.C Abstract/Residuals.cpp Abstract/Status.C Abstract/Scaler.C
        QpStoch/Monitor.C Abstract/IotrRefCount.C Abstract/DoubleLinearSolver.C InteriorPointMethod/Solver.C Vector/Vector.cpp Vector/SimpleVector.C Utilities/drand.C Utilities/sort.C
        Utilities/pipschecks.C QpGen/QpGenOptions.C QpStoch/StochOptions.C Problems/ProblemFactory.cpp)

add_library(ooqpgensparse STATIC
        Abstract/Variables.C Problems/QP.cpp LinearSystems/LinearSystem.C #QpGen
        Problems/QPFactory.C LinearSystems/SparseLinearSystem.C #QpGenSparse
        Readers/MpsReader.C Readers/hash.C #Readers
        LinearSolvers/Preconditioners/SCsparsifier.C
        )

add_library(ooqpstoch STATIC
        QpStoch/DistributedFactory.C
        Problems/DistributedQP.cpp
        LinearSystems/DistributedLinearSystem.C
        QpStoch/sLinsysRoot.C
        QpStoch/sLinsysRootAug.C
        QpStoch/sLinsysRootAugHierInner.C
        QpStoch/sLinsysLeaf.C
        QpStoch/sLinsysLeafSchurSlv.C
        QpStoch/DistributedVariables.C
        QpStoch/DistributedResiduals.cpp
        QpStoch/sTree.C
        QpStoch/sLinsysRootBordered.C
        QpStoch/sTreeCallbacks.C
        LinearAlgebra/Distributed/StringGenMatrix.C # TODO : bad hack...
        QpStoch/RegularizationStrategy.h QpStoch/RegularizationStrategy.cpp)

if (HAVE_MA27)
    target_compile_definitions(ooqpstoch PUBLIC -DWITH_MA27)
    target_link_libraries(ooqpstoch PRIVATE ma27_solver)
endif (HAVE_MA27)

if (HAVE_MA57)
    target_compile_definitions(ooqpstoch PUBLIC -DWITH_MA57)
    target_link_libraries(ooqpstoch PRIVATE ma57_solver)
endif (HAVE_MA57)

if (HAVE_MUMPS)
    target_sources(ooqpstoch QpStoch/sLinsysLeafMumps.C)
    target_compile_definitions(ooqpstoch PUBLIC -DWITH_MUMPS)
    target_link_libraries(ooqpstoch PRIVATE mumps_solver)
endif (HAVE_MUMPS)

if (HAVE_PARDISO)
    target_compile_definitions(ooqpstoch PUBLIC -DWITH_PARDISO)
    target_link_libraries(ooqpstoch PRIVATE pardiso_solver)
endif (HAVE_PARDISO)

if (HAVE_MKL_PARDISO)
    target_compile_definitions(ooqpstoch PUBLIC -DWITH_MKL_PARDISO)
    target_link_libraries(ooqpstoch PRIVATE pardiso_mkl_solver)
endif (HAVE_MKL_PARDISO)

if (HAVE_PARDISO OR HAVE_MKL_PARDISO)
    target_link_libraries(ooqpstoch PUBLIC pardiso_abstract)
endif (HAVE_PARDISO OR HAVE_MKL_PARDISO)

target_link_libraries(ooqpstoch
        PRIVATE MPI::MPI_CXX
        PRIVATE ${MATH_LIBS}
        )
