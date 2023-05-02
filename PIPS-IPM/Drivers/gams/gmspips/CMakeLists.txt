
set(GAMSDIR ..)
set(GAMSSTATDIR ../statgdx)

if (BUILD_GDX_SOURCE)
    add_library(statgdx
            ${GAMSSTATDIR}/clibtypes.cpp
            ${GAMSSTATDIR}/datastorage.cpp
            ${GAMSSTATDIR}/dtoaLoc.cpp
            ${GAMSSTATDIR}/gamsxcc.cpp
            ${GAMSSTATDIR}/gcmt.cpp
            ${GAMSSTATDIR}/gdlaudit.cpp
            ${GAMSSTATDIR}/gdxdclib.cpp
            ${GAMSSTATDIR}/gdxstatic.cpp
            ${GAMSSTATDIR}/gmsdata.cpp
            ${GAMSSTATDIR}/gmsgen.cpp
            ${GAMSSTATDIR}/gmsglob.cpp
            ${GAMSSTATDIR}/gmsglobx.cpp
            ${GAMSSTATDIR}/gmsheapnew.cpp
            ${GAMSSTATDIR}/gmslibname.cpp
            ${GAMSSTATDIR}/gmsobj.cpp
            ${GAMSSTATDIR}/gmsspecs.cpp
            ${GAMSSTATDIR}/gmsstrm.cpp
            ${GAMSSTATDIR}/gxdefs.cpp
            ${GAMSSTATDIR}/gxfile.cpp
            ${GAMSSTATDIR}/idglobal_p3.cpp
            ${GAMSSTATDIR}/math_p3.cpp
            ${GAMSSTATDIR}/optcc.cpp
            ${GAMSSTATDIR}/p3Custom2.cpp
            ${GAMSSTATDIR}/p3ieeefp.cpp
            ${GAMSSTATDIR}/p3io.cpp
            ${GAMSSTATDIR}/p3library.cpp
            ${GAMSSTATDIR}/p3platform.cpp
            ${GAMSSTATDIR}/p3private.cpp
            ${GAMSSTATDIR}/p3process.cpp
            ${GAMSSTATDIR}/p3threads.cpp
            ${GAMSSTATDIR}/p3utils.cpp
            ${GAMSSTATDIR}/paldoorg.cpp
            ${GAMSSTATDIR}/palmdcon.cpp
            ${GAMSSTATDIR}/pchutil.cpp
            ${GAMSSTATDIR}/runner.cpp
            ${GAMSSTATDIR}/strhash.cpp
            ${GAMSSTATDIR}/strutilx.cpp
            ${GAMSSTATDIR}/system_p3.cpp
            ${GAMSSTATDIR}/sysutils_p3.cpp
            ${GAMSSTATDIR}/xcompress.cpp
            )

    target_compile_definitions(statgdx
            PRIVATE _P3_DllExport=
            PRIVATE STDCALL=
            PRIVATE GDX_CALLCONV=unwanted
            PRIVATE LEG
            PRIVATE GC_NO_MUTEX
            )

    target_compile_options(statgdx
            PRIVATE -Wno-register
            PRIVATE -Wno-stringop-overflow
            PRIVATE -Wno-unused-result
            PRIVATE -Wno-unused-value
            PRIVATE -Wno-parentheses
            PRIVATE -Wno-unused-local-typedefs
            PRIVATE -Wno-dangling-else
            PRIVATE -Wno-unused-function
            PRIVATE -Wno-unused-but-set-variable
            PRIVATE -Wno-unused-variable
            PRIVATE -Wno-sign-compare
            PRIVATE -Wno-format-truncation
            PRIVATE -Wno-narrowing
            PRIVATE "$<$<COMPILE_LANGUAGE:C>:-Wno-implicit-function-declaration>"
            )

    target_link_libraries(statgdx
            PRIVATE Threads::Threads
            )

    add_library(gmspipsio gmspipsio.c)

    target_link_libraries(gmspipsio
            PRIVATE statgdx
            )
    target_compile_definitions(gmspipsio
            PRIVATE GDXSOURCE
            )
    target_compile_options(gmspipsio
            PRIVATE "$<$<COMPILE_LANGUAGE:C>:-Wno-implicit-function-declaration>"
            )
    target_include_directories(gmspipsio
            PUBLIC ${GAMSSTATDIR}
            )
else (BUILD_GDX_SOURCE)
    add_library(gmspipsio
            gmspipsio.c
            ${GAMSDIR}/apifiles/C/api/gmomcc.c
            ${GAMSDIR}/apifiles/C/api/gevmcc.c
            ${GAMSDIR}/apifiles/C/api/gdxcc.c
            )
    target_include_directories(gmspipsio
            PUBLIC ${GAMSDIR}/apifiles/C/api/
            )
endif (BUILD_GDX_SOURCE)

target_include_directories(gmspipsio
        PUBLIC .
        )
target_compile_definitions(gmspipsio
        PRIVATE GC_NO_MUTEX
        )
target_link_libraries(gmspipsio
        PRIVATE ${DL_LIBRARY}
        )

########### gmpipsreader ###########
add_library(gmspips_read_write gmspips_reader.cpp)

target_include_directories(gmspips_read_write
        PRIVATE ../PIPS-IPM/Core/Readers/Distributed
        PRIVATE ../PIPS-IPM/Core/Utilities)

target_link_libraries(gmspips_read_write
        PRIVATE gmspipsio
        PRIVATE pips-ipmpp
        PRIVATE MPI::MPI_CXX
        )

############## gmschk ##############
add_executable(gmschk
        gmspipschk.cpp
        )

target_link_libraries(gmschk
        PRIVATE gmspipsio
        PRIVATE OpenMP::OpenMP_CXX
        )

############## gmspips ##############
add_executable(gmspips gmspips.cpp)

set(DIRECTORIES
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/Interface
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/KKTFormulation/Variables
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/KKTFormulation/Residuals
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/LinearAlgebra/Abstract
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/LinearAlgebra/Dense
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/LinearAlgebra/Distributed
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/LinearAlgebra/Sparse
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/Preprocessing
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/Problems
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/Readers/Distributed
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/LinearSolvers
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/LinearSolvers/Ma27Solver
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/LinearSolvers/Ma57Solver
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/LinearSolvers/PardisoSolver/PardisoIndefSolver
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/LinearSolvers/PardisoSolver/PardisoSchurSolver
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/LinearSolvers/PardisoSolver/PardisoSolver
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Core/Options
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Drivers
        ${PROJECT_SOURCE_DIR}/PIPS-IPM/Drivers/gams/gmspips
)


target_include_directories(gmspips PUBLIC ${DIRECTORIES})
target_include_directories(gmspips PRIVATE ${includes_for_interface})
target_compile_options(gmspips
        PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-Wextra>"
        PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-Wsuggest-override>"
        )

#target_compile_definitions(gmspips
#        PRIVATE GMS_LOG
#        )

target_link_libraries(gmspips
        PRIVATE pips-ipmpp
        PRIVATE gmspips_read_write
        PRIVATE OpenMP::OpenMP_CXX
        PRIVATE MPI::MPI_CXX
        )
