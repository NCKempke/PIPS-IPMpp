if(HAVE_MUMPS)
	set(mumps_files ${mumps_files} MumpsSolver/MumpsSolverBase.C)
   
   if(WITH_MUMPS_LEAF)
      set(mumps_files ${mumps_files} MumpsSolver/MumpsSolverLeaf.C)
   endif(WITH_MUMPS_LEAF)
   
   if(WITH_MUMPS_ROOT)
      set(mumps_files ${mumps_files} MumpsSolver/MumpsSolverRoot.C)
   endif(WITH_MUMPS_ROOT)
   
   add_library( mumps_solver STATIC ${mumps_files} )
   target_include_directories( mumps_solver PRIVATE ${MUMPS_INCLUDE_DIR})
   target_link_libraries( mumps_solver PRIVATE ${MUMPS_LIBRARY} ${PARMETIS_LIBRARY} ${MATH_LIBS} blas lapack esmumps scotch scotcherr gfortran pthread m )
endif(HAVE_MUMPS)


if(HAVE_MA27)
   add_library( ma27_solver STATIC Mc30Scaler/Mc30Scaler.C Ma27Solver/Ma27Solver.C Ma27Solver/Ma27SolverRoot.C ../QpGen/QpGenSparseMa27.C )
   target_link_libraries( ma27_solver PRIVATE ${MA27_LIBRARY} )
endif(HAVE_MA27)


if(HAVE_MA57)
   add_library( ma57_solver STATIC Ma57Solver/Ma57Solver.C Ma57Solver/Ma57SolverRoot.C ../QpGen/QpGenSparseMa57.C )
   target_link_libraries( ma57_solver ${MA57_LIBRARY} )
endif(HAVE_MA57)


if( HAVE_PARDISO )
   add_library( pardiso_solver STATIC PardisoSolver/PardisoSolver.C PardisoSolver/PardisoProjectSolver.C PardisoSolver/PardisoSchurSolver.C PardisoSolver/PardisoIndefSolver.C)
   target_link_libraries( pardiso_solver PRIVATE ${PARDISO_LIBRARY} ${MATH_LIBS} gfortran pthread m )
endif(HAVE_PARDISO)

if( HAVE_MKL_PARDISO )
   add_library( pardiso_mkl_solver STATIC PardisoSolver/PardisoSolver.C PardisoSolver/PardisoMKLSolver.C )
   target_link_libraries( pardiso_mkl_solver PRIVATE ${MKL_PARDISO_LIBS} )
endif( HAVE_MKL_PARDISO )