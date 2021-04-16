add_subdirectory(QpSolvers)
add_subdirectory(QpPreprocess)
add_subdirectory(SparseLinearAlgebra)
add_subdirectory(DenseLinearAlgebra)
add_subdirectory(StochLinearAlgebra)
add_subdirectory(LinearSolvers)

add_library(ooqpbase STATIC Abstract/OoqpVersion.C Abstract/Options.C Abstract/Observer.C Abstract/Variables.C Abstract/Problem.C Abstract/Residuals.cpp Abstract/Solver.C Abstract/Status.C Abstract/Scaler.C
  Abstract/OoqpMonitor.C Abstract/IotrRefCount.C Abstract/DoubleLinearSolver.C 
  Vector/OoqpVector.C Vector/SimpleVector.C
Utilities/drand.C Utilities/sort.C Utilities/pipschecks.C QpGen/QpGenOptions.C QpStoch/StochOptions.C
        Abstract/Residuals.cpp)

add_library(ooqpgensparse STATIC
  QpGen/QpGenVars.C QpGen/QuadraticProblem.C QpGen/QpGen.C QpGen/QpGenLinsys.C #QpGen
  QpGen/QpGenSparseSeq.C QpGen/QpGenSparseLinsys.C #QpGenSparse
  Readers/MpsReader.C Readers/hash.C #Readers
  LinearSolvers/Preconditioners/SCsparsifier.C
)

add_library(ooqpstoch STATIC 
   QpStoch/sFactory.C
   QpStoch/sFactoryAug.C
   QpStoch/sData.C
   QpStoch/sLinsys.C
   QpStoch/sLinsysRoot.C
   QpStoch/sLinsysRootAug.C
   QpStoch/sLinsysRootAugHierInner.C
   QpStoch/sLinsysLeaf.C
   QpStoch/sLinsysLeafSchurSlv.C
   QpStoch/sVars.C
   QpStoch/StochMonitor.C
   QpStoch/sResiduals.C
   QpStoch/sTree.C
   QpStoch/sLinsysRootBordered.C
   QpStoch/sTreeCallbacks.C
   StochLinearAlgebra/StringGenMatrix.C # TODO : bad hack...
)

if( HAVE_MA27 )
   target_compile_definitions( ooqpstoch PUBLIC -DWITH_MA27 )
   target_link_libraries( ooqpstoch PRIVATE ma27_solver )
endif( HAVE_MA27 )

if( HAVE_MA57 )
   target_compile_definitions( ooqpstoch PUBLIC -DWITH_MA57 )
   target_link_libraries( ooqpstoch PRIVATE ma57_solver )
endif( HAVE_MA57 )

if(HAVE_MUMPS)
   target_sources( ooqpstoch QpStoch/sLinsysLeafMumps.C)
   target_compile_definitions( ooqpstoch PUBLIC -DWITH_MUMPS )
   target_link_libraries(ooqpstoch PRIVATE mumps_solver)
endif(HAVE_MUMPS)

if(HAVE_PARDISO)
   target_compile_definitions( ooqpstoch PUBLIC -DWITH_PARDISO )
   target_link_libraries( ooqpstoch PRIVATE pardiso_solver )
endif(HAVE_PARDISO)

if( HAVE_MKL_PARDISO )
   target_compile_definitions( ooqpstoch PUBLIC -DWITH_MKL_PARDISO )
   target_link_libraries(ooqpstoch PRIVATE pardiso_mkl_solver)
endif( HAVE_MKL_PARDISO)

if( HAVE_PARDISO OR HAVE_MKL_PARDISO )
   target_link_libraries( ooqpstoch PUBLIC pardiso_abstract )
endif( HAVE_PARDISO OR HAVE_MKL_PARDISO )

target_link_libraries( ooqpstoch 
   PRIVATE MPI::MPI_CXX
   PRIVATE ${MATH_LIBS}
)









