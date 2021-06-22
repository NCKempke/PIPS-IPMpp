target_include_directories(pips-ipmpp PRIVATE ./)
target_include_directories(pips-ipmpp PRIVATE
        ./
        BiCGStabSolver
        CGSolver
        DensePSDSolver
        DenseSymmetricIndefinitSolver
        PCGSolver
        Preconditioners
    )

target_sources(pips-ipmpp
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/DoubleLinearSolver.C
            ${CMAKE_CURRENT_SOURCE_DIR}/BiCGStabSolver/BiCGStabSolver.C
            ${CMAKE_CURRENT_SOURCE_DIR}/CGSolver/CGSolver.C
            ${CMAKE_CURRENT_SOURCE_DIR}/DensePSDSolver/DeSymPSDSolver.C
            ${CMAKE_CURRENT_SOURCE_DIR}/DenseSymmetricIndefinitSolver/DeSymIndefSolver.C
            ${CMAKE_CURRENT_SOURCE_DIR}/DenseSymmetricIndefinitSolver/DeSymIndefSolver2.C
            ${CMAKE_CURRENT_SOURCE_DIR}/PCGSolver/PCGSolver.C
            ${CMAKE_CURRENT_SOURCE_DIR}/Preconditioners/SCsparsifier.C
        )

if (HAVE_MUMPS)
    add_library(mumps_solver SHARED MumpsSolver/MumpsSolverBase.C MumpsSolver/MumpsSolverLeaf.C MumpsSolver/MumpsSolverRoot.C)
    target_include_directories(mumps_solver PRIVATE ${MUMPS_INCLUDE_DIR})
    target_link_libraries(mumps_solver
            PRIVATE ${MUMPS_D_LIBRARY}
            PRIVATE ${MUMPS_COMMON_LIBRARY}
            PRIVATE ${MUMPS_FORT_LIB}
            PRIVATE ${MUMPS_PORD_LIBRARY}
            PRIVATE ${SCALAPACK_LIBRARIES}
            PRIVATE ${PARMETIS_LIBRARY}
            PRIVATE ${MATH_LIBS}
            PRIVATE esmumps PRIVATE scotch PRIVATE scotcherr PRIVATE gfortran PRIVATE pthread PRIVATE m)
    target_link_libraries(mumps_solver PRIVATE ${METIS_5_LIBRARY})

    target_include_directories(pips-ipmpp PRIVATE ${MUMPS_INCLUDE_DIR})
    target_sources(pips-ipmpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../KKTFormulation/LinearSystems/sLinsysLeafMumps.C)
    target_compile_definitions(pips-ipmpp PUBLIC -DWITH_MUMPS)

    list(APPEND lin_solvers_list "mumps_solver")
    set(lin_solvers_list ${lin_solvers_list} PARENT_SCOPE)
endif (HAVE_MUMPS)

if (HAVE_MA27)
    add_library(ma27_solver STATIC Mc30Scaler/Mc30Scaler.C Ma27Solver/Ma27Solver.C Ma27Solver/Ma27SolverRoot.C)
    target_include_directories(ma27_solver PRIVATE ./)
    target_link_libraries(ma27_solver
            PRIVATE ${MA27_LIBRARY}
            PRIVATE ${GFORTRAN_LIBRARY}
            )

    target_include_directories(pips-ipmpp PRIVATE Mc30Scaler Ma27Solver)
    target_compile_definitions(pips-ipmpp PUBLIC -DWITH_MA27)
    list(APPEND lin_solvers_list "ma27_solver")
    set(lin_solvers_list ${lin_solvers_list} PARENT_SCOPE)

endif (HAVE_MA27)

if (HAVE_MA57)
    add_library(ma57_solver STATIC Ma57Solver/Ma57Solver.C Ma57Solver/Ma57SolverRoot.C)
    target_link_libraries(ma57_solver
            PRIVATE ${MA57_LIBRARY}
            PRIVATE ${METIS_4_LIBRARY}
            PRIVATE ${GFORTRAN_LIBRARY}
            )

    target_include_directories(pips-ipmpp PRIVATE Ma57Solver)
    target_compile_definitions(pips-ipmpp PUBLIC -DWITH_MA57)

    list(APPEND lin_solvers_list "ma57_solver")
    set(lin_solvers_list ${lin_solvers_list} PARENT_SCOPE)
endif (HAVE_MA57)

if (HAVE_PARDISO OR HAVE_MKL_PARDISO)
    add_library(pardiso_abstract STATIC
            PardisoSolver/PardisoSolver/PardisoSolver.C
            PardisoSolver/PardisoSchurSolver/PardisoSchurSolver.C
            PardisoSolver/PardisoIndefSolver/PardisoIndefSolver.C
            )

    target_include_directories(pips-ipmpp PRIVATE
            PardisoSolver/PardisoSolver
            PardisoSolver/PardisoSchurSolver
            PardisoSolver/PardisoIndefSolver
            )

    list(APPEND lin_solvers_list "pardiso_abstract")
    set(lin_solvers_list ${lin_solvers_list} PARENT_SCOPE)

endif (HAVE_PARDISO OR HAVE_MKL_PARDISO)

if (HAVE_PARDISO)
    add_library(pardiso_solver STATIC
            PardisoSolver/PardisoSolver/PardisoProjectSolver.C
            PardisoSolver/PardisoIndefSolver/PardisoProjectIndefSolver.C
            PardisoSolver/PardisoSchurSolver/PardisoProjectSchurSolver.C
            )
    set_target_properties(pardiso_solver PROPERTIES CXX_VISIBILITY_PRESET hidden)
    target_link_libraries(pardiso_solver
            PRIVATE pardiso_abstract
            PRIVATE ${PARDISO_LIBRARY}
            PRIVATE ${MATH_LIBS}
            PRIVATE gfortran
            PRIVATE pthread
            PRIVATE m
            )
    target_compile_definitions(pips-ipmpp PUBLIC -DWITH_PARDISO)
    target_sources(pips-ipmpp PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/../KKTFormulation/LinearSystems/sLinsysLeafSchurSlv.C
            )
    list(APPEND lin_solvers_list "pardiso_solver")
    set(lin_solvers_list ${lin_solvers_list} PARENT_SCOPE)
endif (HAVE_PARDISO)

if (HAVE_MKL_PARDISO)
    add_library(pardiso_mkl_solver STATIC
            PardisoSolver/PardisoSolver/PardisoMKLSolver.C
            PardisoSolver/PardisoSchurSolver/PardisoMKLSchurSolver.C
            PardisoSolver/PardisoIndefSolver/PardisoMKLIndefSolver.C
            )
    target_link_libraries(pardiso_mkl_solver PRIVATE pardiso_abstract PRIVATE ${MKL_PARDISO_LIBS})
    target_compile_definitions(pips-ipmpp PUBLIC -DWITH_MKL_PARDISO)

    list(APPEND lin_solvers_list "pardiso_mkl_solver")
    set(lin_solvers_list ${lin_solvers_list} PARENT_SCOPE)
endif (HAVE_MKL_PARDISO)