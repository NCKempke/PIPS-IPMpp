if(HAVE_MUMPS)
   add_library( mumps_solver SHARED MumpsSolver/MumpsSolverBase.C MumpsSolver/MumpsSolverLeaf.C MumpsSolver/MumpsSolverRoot.C )
   target_link_libraries( mumps_solver 
      PRIVATE ${MUMPS_D_LIBRARY}
      PRIVATE ${MUMPS_COMMON_LIBRARY} 
      PRIVATE ${MUMPS_FORT_LIB} 
      PRIVATE ${MUMPS_PORD_LIBRARY} 
      PRIVATE ${SCALAPACK_LIBRARIES} 
      #PRIVATE ${PARMETIS_LIBRARY}
      PRIVATE ${MATH_LIBS} 
      PRIVATE blas
   #PRIVATE lapack
PRIVATE esmumps PRIVATE scotch PRIVATE scotcherr PRIVATE gfortran PRIVATE pthread PRIVATE m )
   target_link_libraries( mumps_solver PRIVATE ${METIS_5_LIBRARY} )
endif(HAVE_MUMPS)

if(HAVE_MA27)
   add_library( ma27_solver SHARED Mc30Scaler/Mc30Scaler.C Ma27Solver/Ma27Solver.C Ma27Solver/Ma27SolverRoot.C )# ../QpGen/QpGenSparseMa27.C )
   target_link_libraries( ma27_solver PRIVATE ${MA27_LIBRARY} )
endif(HAVE_MA27)

if(HAVE_MA57)
   add_library( ma57_solver SHARED Ma57Solver/Ma57Solver.C Ma57Solver/Ma57SolverRoot.C )#../QpGen/QpGenSparseMa57.C )
   target_link_libraries( ma57_solver 
      PRIVATE ${MA57_LIBRARY}
      PRIVATE ${METIS_4_LIBRARY} )
endif(HAVE_MA57)

if( HAVE_PARDISO OR HAVE_MKL_PARDISO )
   add_library( pardiso_abstract SHARED
      PardisoSolver/PardisoSolver/PardisoSolver.C
      PardisoSolver/PardisoSchurSolver/PardisoSchurSolver.C
      PardisoSolver/PardisoIndefSolver/PardisoIndefSolver.C   
)
endif( HAVE_PARDISO OR HAVE_MKL_PARDISO)

if( HAVE_PARDISO )
   add_library( pardiso_solver SHARED 
      PardisoSolver/PardisoSolver/PardisoProjectSolver.C
      PardisoSolver/PardisoIndefSolver/PardisoProjectIndefSolver.C
      PardisoSolver/PardisoSchurSolver/PardisoProjectSchurSolver.C
   )
   target_link_libraries( pardiso_solver 
      PRIVATE pardiso_abstract 
      PRIVATE ${PARDISO_LIBRARY}
      PRIVATE ${MATH_LIBS} 
      PRIVATE gfortran
      PRIVATE pthread
      PRIVATE m
   )
endif(HAVE_PARDISO)

if( HAVE_MKL_PARDISO )
   add_library( pardiso_mkl_solver SHARED 
      PardisoSolver/PardisoSolver/PardisoMKLSolver.C
      PardisoSolver/PardisoSchurSolver/PardisoMKLSchurSolver.C
      PardisoSolver/PardisoIndefSolver/PardisoMKLIndefSolver.C
   )
   target_link_libraries( pardiso_mkl_solver PUBLIC pardiso_abstract PRIVATE ${MKL_PARDISO_LIBS} )
endif( HAVE_MKL_PARDISO )
