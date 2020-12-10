add_subdirectory(QpSolvers)
add_subdirectory(QpPreprocess)
add_subdirectory(SparseLinearAlgebra)
add_subdirectory(DenseLinearAlgebra)
add_subdirectory(StochLinearAlgebra)
add_subdirectory(QpStoch)
add_subdirectory(LinearSolvers)

add_library(ooqpbase STATIC Abstract/OoqpVersion.C Abstract/Options.C Abstract/Observer.C Abstract/Variables.C Abstract/Data.C Abstract/Solver.C Abstract/Status.C Abstract/Scaler.C
  Abstract/OoqpMonitor.C Abstract/IotrRefCount.C Abstract/DoubleLinearSolver.C 
  Vector/OoqpVector.C Vector/SimpleVector.C Vector/VectorUtilities.C
  Utilities/drand.C Utilities/sort.C Utilities/pipschecks.C)

add_library(ooqpgensparse 
  QpGen/QpGenVars.C QpGen/QpGenData.C QpGen/QpGenResiduals.C QpGen/QpGen.C QpGen/QpGenLinsys.C QpGen/QpGenOptions.C #QpGen
  QpGen/QpGenSparseSeq.C QpGen/QpGenSparseLinsys.C #QpGenSparse
  Readers/MpsReader.C Readers/hash.C #Readers
  LinearSolvers/Preconditioners/SCsparsifier.C
)

## TODO : does not make much sense to link the solvers here I guess - later (when sfactory ist there) would be better

if( HAVE_MA27 )
   target_compile_definitions( ooqpgensparse PUBLIC -DWITH_MA27 )
   target_link_libraries( ooqpgensparse PRIVATE ma27_solver )
endif( HAVE_MA27 )

if( HAVE_MA57 )
   target_compile_definitions( ooqpgensparse PUBLIC -DWITH_MA57 )
   target_link_libraries( ooqpgensparse PRIVATE ma57_solver )
endif( HAVE_MA57 )

#if(HAVE_MUMPS)
#   target_compile_definitions( ooqpgensparse PUBLIC -DWITH_MUMPS )
#   target_link_libraries(ooqpgensparse PRIVATE mumps_solver)
#endif(HAVE_MUMPS)

if(HAVE_PARDISO)
   target_compile_definitions( ooqpgensparse PUBLIC -DWITH_PARDISO )
   target_link_libraries(ooqpgensparse PRIVATE pardiso_solver)
endif(HAVE_PARDISO)

if( HAVE_MKL_PARDISO )
   target_compile_definitions( ooqpgensparse PUBLIC -DWITH_MKL_PARDISO )
   target_link_libraries(ooqpgensparse PRIVATE pardiso_mkl_solver)
endif( HAVE_MKL_PARDISO)
   