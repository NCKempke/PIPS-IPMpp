add_subdirectory(QpSolvers)
add_subdirectory(QpPreprocess)
add_subdirectory(SparseLinearAlgebra)
add_subdirectory(DenseLinearAlgebra)
add_subdirectory(StochLinearAlgebra)
add_subdirectory(QpStoch)

add_library(ooqpbase STATIC Abstract/OoqpVersion.C Abstract/Options.C Abstract/Observer.C Abstract/Variables.C Abstract/Data.C Abstract/Solver.C Abstract/Status.C Abstract/Scaler.C
  Abstract/OoqpMonitor.C Abstract/IotrRefCount.C Abstract/DoubleLinearSolver.C 
  Vector/OoqpVector.C Vector/SimpleVector.C Vector/VectorUtilities.C
  Utilities/drand.C Utilities/sort.C Utilities/pipschecks.C)

set(SOLVERS_LIBRARIES "")
if(HAVE_MA27)
	set(solvers LinearSolvers/Mc30Scaler/Mc30Scaler.C LinearSolvers/Ma27Solver/Ma27Solver.C LinearSolvers/Ma27Solver/Ma27SolverRoot.C QpGen/QpGenSparseMa27.C )
   set(SOLVERS_LIBRARIES "${SOLVERS_LIBRARIES};${MA27_LIBRARY}")
endif(HAVE_MA27)
if(HAVE_MA57)
	set(solvers ${solvers} LinearSolvers/Ma57Solver/Ma57Solver.C LinearSolvers/Ma57Solver/Ma57SolverRoot.C QpGen/QpGenSparseMa57.C)
   set(SOLVERS_LIBRARIES "${SOLVERS_LIBRARIES};${MA57_LIBRARY}")
endif(HAVE_MA57)
#if(HAVE_MA86)
#	set(solvers ${solvers} LinearSolvers/Ma86Solver/Ma86Solver.C QpGen/QpGenSparseMa86.C)
#endif(HAVE_MA86)

if(HAVE_MUMPS)
	set(mumps_files ${mumps_files} LinearSolvers/MumpsSolver/MumpsSolverBase.C)
   
   if(WITH_MUMPS_LEAF)
      set(mumps_files ${mumps_files} LinearSolvers/MumpsSolver/MumpsSolverLeaf.C)
   endif(WITH_MUMPS_LEAF)
   
   if(WITH_MUMPS_ROOT)
      set(mumps_files ${mumps_files} LinearSolvers/MumpsSolver/MumpsSolverRoot.C)
   endif(WITH_MUMPS_ROOT)
   
   add_library( mumps_solver STATIC ${mumps_files} )
   target_link_libraries( mumps_solver PRIVATE ${MUMPS_LIBRARY} ${PARMETIS_LIBRARY} ${MATH_LIBS} blas lapack esmumps scotch scotcherr gfortran pthread m )
endif(HAVE_MUMPS)

if(HAVE_PARDISO)
   add_library( pardiso_solver STATIC LinearSolvers/PardisoSolver/PardisoSolver.C LinearSolvers/PardisoSolver/PardisoSchurSolver.C LinearSolvers/PardisoSolver/PardisoIndefSolver.C)

   if(HAVE_MKL_PARDISO AND WITH_MKL_PARDISO)
      message(STATUS "Linking MKL PARDISO")
      target_link_libraries( pardiso_solver PRIVATE ${MKL_PARDISO_LIBS} )
   else(HAVE_MKL_PARDISO AND WITH_MKL_PARDISO)
      message(STATUS "Linking Schenk PARDISO")
      target_link_libraries( pardiso_solver PRIVATE ${PARDISO_LIBRARY} ${MATH_LIBS} gfortran pthread m )
   endif(HAVE_MKL_PARDISO AND WITH_MKL_PARDISO)
endif(HAVE_PARDISO)

set(solvers ${solvers} LinearSolvers/Preconditioners/SCsparsifier.C)

add_library(ooqpgensparse 
  QpGen/QpGenVars.C QpGen/QpGenData.C QpGen/QpGenResiduals.C QpGen/QpGen.C QpGen/QpGenLinsys.C QpGen/QpGenOptions.C #QpGen
  QpGen/QpGenSparseSeq.C QpGen/QpGenSparseLinsys.C #QpGenSparse
  Readers/MpsReader.C Readers/hash.C #Readers
  ${solvers})

target_link_libraries( ooqpgensparse PRIVATE ${METIS_4_LIBRARY})
target_link_libraries( ooqpgensparse PRIVATE ${SOLVERS_LIBRARIES})

if(HAVE_MUMPS)
   target_link_libraries(ooqpgensparse PRIVATE mumps_solver)
endif(HAVE_MUMPS)

if(HAVE_PARDISO)
   target_link_libraries(ooqpgensparse PRIVATE pardiso_solver)
endif(HAVE_PARDISO)
