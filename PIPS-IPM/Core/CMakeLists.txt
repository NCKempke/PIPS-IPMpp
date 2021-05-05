add_subdirectory(LinearSolvers)

add_library(pips-ipmpp STATIC
        Abstract/Options.C Abstract/Observer.C Abstract/Variables.C Problems/Problem.C Abstract/Residuals.cpp Abstract/Status.C Abstract/Scaler.C
        #Filter.cpp FilterStrategy.cpp FilterLineSearch.cpp
        InteriorPointMethod/Solver.C InteriorPointMethod/GondzioStochSolver.C InteriorPointMethod/GondzioStochLpSolver.C InteriorPointMethod/InteriorPointMethod.cpp
        QpStoch/Monitor.C Abstract/IotrRefCount.C Abstract/DoubleLinearSolver.C Vector/Vector.cpp Vector/SimpleVector.C Utilities/drand.C Utilities/sort.C
        Utilities/pipschecks.C Abstract/QpGenOptions.C QpStoch/DistributedOptions.C
        Abstract/Variables.C Problems/QP.cpp LinearSystems/LinearSystem.C #QpGen
        Readers/MpsReader.C Readers/hash.C #Readers
        LinearSolvers/Preconditioners/SCsparsifier.C
        Problems/DistributedQP.cpp
        Problems/DistributedFactory.C 
        InteriorPointMethod/Solver.C
        LinearSystems/DistributedLinearSystem.C
        LinearSystems/DistributedRootLinearSystem.C
        LinearSystems/sLinsysRootAug.C
        LinearSystems/sLinsysRootAugHierInner.C
        LinearSystems/DistributedLeafLinearSystem.C
        LinearSystems/sLinsysLeafSchurSlv.C
        QpStoch/DistributedVariables.C
        QpStoch/DistributedResiduals.cpp
        QpStoch/DistributedTree.C
        LinearSystems/sLinsysRootBordered.C
        QpStoch/DistributedTreeCallbacks.C
        LinearAlgebra/Distributed/StripMatrix.C # TODO : bad hack...
        LinearSystems/RegularizationStrategy.h LinearSystems/RegularizationStrategy.cpp
        
        LinearAlgebra/Dense/DenseStorage.C LinearAlgebra/Dense/DenseSymMatrix.C
        LinearAlgebra/Dense/DeSymIndefSolver.C LinearAlgebra/Dense/DeSymIndefSolver2.C
        LinearAlgebra/Dense/DenseMatrix.C LinearAlgebra/Dense/DeSymPSDSolver.C
        
        LinearAlgebra/Distributed/DistributedVector.C LinearAlgebra/Distributed/DistributedSymmetricMatrix.C LinearAlgebra/Distributed/DistributedMatrix.C
        LinearAlgebra/Distributed/BorderedMatrix.C LinearAlgebra/Distributed/BorderedSymmetricMatrix.C
        LinearAlgebra/Distributed/StochResourcesMonitor.C LinearAlgebra/Distributed/DistributedInputTree.C
        
        LinearAlgebra/Sparse/SparseStorage.C LinearAlgebra/Sparse/SparseStorageDynamic.C
        LinearAlgebra/Sparse/SparseMatrix.C LinearAlgebra/Sparse/SparseSymmetricMatrix.C
        
        Preprocessing/EquiStochScaler.C Preprocessing/GeoStochScaler.C Preprocessing/StochScaler.C Preprocessing/QpScaler.C
        Preprocessing/StochRowStorage.C Preprocessing/StochColumnStorage.C Preprocessing/QpPostsolver.C Preprocessing/StochPostsolver.C Preprocessing/StochPresolver.C Preprocessing/QpPresolver.C Preprocessing/StochPresolverBase.C
        Preprocessing/StochPresolverModelCleanup.C Preprocessing/PresolveData.C Preprocessing/StochPresolverBoundStrengthening.C Preprocessing/StochPresolverColumnFixation.C Preprocessing/StochPresolverParallelRows.C
        Preprocessing/StochPresolverSingletonRows.C Preprocessing/StochPresolverSingletonColumns.C Preprocessing/QpPresolver.C Preprocessing/QpPostsolver.C
)

if (HAVE_MA27)
    target_compile_definitions(pips-ipmpp PUBLIC -DWITH_MA27)
    target_link_libraries(pips-ipmpp PRIVATE ma27_solver)
endif (HAVE_MA27)

if (HAVE_MA57)
    target_compile_definitions(pips-ipmpp PUBLIC -DWITH_MA57)
    target_link_libraries(pips-ipmpp PRIVATE ma57_solver)
endif (HAVE_MA57)

if (HAVE_MUMPS)
    target_sources(pips-ipmpp QpStoch/sLinsysLeafMumps.C)
    target_compile_definitions(pips-ipmpp PUBLIC -DWITH_MUMPS)
    target_link_libraries(pips-ipmpp PRIVATE mumps_solver)
endif (HAVE_MUMPS)

if (HAVE_PARDISO)
    target_compile_definitions(pips-ipmpp PUBLIC -DWITH_PARDISO)
    target_link_libraries(pips-ipmpp PRIVATE pardiso_solver)
endif (HAVE_PARDISO)

if (HAVE_MKL_PARDISO)
    target_compile_definitions(pips-ipmpp PUBLIC -DWITH_MKL_PARDISO)
    target_link_libraries(pips-ipmpp PRIVATE pardiso_mkl_solver)
endif (HAVE_MKL_PARDISO)

if (HAVE_PARDISO OR HAVE_MKL_PARDISO)
    target_link_libraries(pips-ipmpp PUBLIC pardiso_abstract)
endif (HAVE_PARDISO OR HAVE_MKL_PARDISO)

target_link_libraries(pips-ipmpp
        PRIVATE MPI::MPI_CXX
        PRIVATE ${MATH_LIBS}
)
