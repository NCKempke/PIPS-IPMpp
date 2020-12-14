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

add_library(ooqpgensparse STATIC
  QpGen/QpGenVars.C QpGen/QpGenData.C QpGen/QpGenResiduals.C QpGen/QpGen.C QpGen/QpGenLinsys.C QpGen/QpGenOptions.C #QpGen
  QpGen/QpGenSparseSeq.C QpGen/QpGenSparseLinsys.C #QpGenSparse
  Readers/MpsReader.C Readers/hash.C #Readers
  LinearSolvers/Preconditioners/SCsparsifier.C
)

add_library(ooqpstoch STATIC 
   QpStoch/sFactory.C
   QpStoch/sFactoryAug.C
   #   sFactoryAugSchurLeaf.C
   #   sFactoryAugComm2SchurLeaf.C
   QpStoch/sData.C
   QpStoch/sLinsys.C
   QpStoch/sLinsysRoot.C
   QpStoch/sLinsysRootAug.C
   #   sLinsysRootComm2.C
   #   sLinsysRootAugComm2.C
   QpStoch/sLinsysLeaf.C
   QpStoch/sLinsysLeafSchurSlv.C
   QpStoch/sVars.C
   QpStoch/StochMonitor.C
   QpStoch/StochOptions.C
   QpStoch/sResiduals.C
   QpStoch/sTree.C
   QpStoch/sLinsysRootBordered.C
   QpStoch/sTreeCallbacks.C
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

target_link_libraries( ooqpstoch PRIVATE ${MATH_LIBS} )









