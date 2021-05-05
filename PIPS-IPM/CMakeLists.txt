#add_definitions(-DTIMING -DSTOCH_TESTING) # timing output
include_directories(
        Core/Abstract
        Core/InteriorPointMethod
        Core/LinearAlgebra/Sparse
        Core/LinearAlgebra/Dense
        Core/LinearAlgebra/Distributed
        Core/LinearSolvers/Ma27Solver
        Core/LinearSolvers/Ma57Solver
        Core/LinearSolvers/Mc30Scaler
        Core/LinearSolvers/Preconditioners
        Core/LinearSolvers/PardisoSolver/PardisoSolver
        Core/LinearSolvers/PardisoSolver/PardisoIndefSolver
        Core/LinearSolvers/PardisoSolver/PardisoSchurSolver
        Core/LinearSolvers/BiCGStabSolver
        Core/LinearSystems
        Core/Problems
        Core/Preprocessing
        Core/QpStoch
        Core/Readers
        Core/Utilities
        Core/Vector
)

if (HAVE_MUMPS)
    include_directories(Core/LinearSolvers/MumpsSolver)
endif (HAVE_MUMPS)

add_subdirectory(Core)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_executable(pipsipmCallbackExample Drivers/callbackExample.cpp)
target_compile_options(pipsipmCallbackExample PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(pipsipmCallbackExample
        PRIVATE pips-ipmpp
        PRIVATE OpenMP::OpenMP_CXX
        PRIVATE MPI::MPI_CXX
        )

set(GAMSDIR Drivers/gams)
set(GAMSSTATDIR Drivers/statgdx)

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

    add_library(gmspipsio Drivers/gmspips/gmspipsio.c)

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
            Drivers/gmspips/gmspipsio.c
            ${GAMSDIR}/apifiles/C/api/gmomcc.c
            ${GAMSDIR}/apifiles/C/api/gevmcc.c
            ${GAMSDIR}/apifiles/C/api/gdxcc.c
            )
    target_include_directories(gmspipsio
            PUBLIC ${GAMSDIR}/apifiles/C/api/
            )
endif (BUILD_GDX_SOURCE)

target_include_directories(gmspipsio
        PUBLIC Drivers/gmspips/
        )
target_compile_definitions(gmspipsio
        PRIVATE GC_NO_MUTEX
        )
target_link_libraries(gmspipsio
        PRIVATE ${DL_LIBRARY}
        )

############## gmspips ##############
add_executable(gmspips Drivers/gmspips/gmspips.cpp)

target_compile_options(gmspips
        PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-Wextra>"
        PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:-Wsuggest-override>"
        )

target_compile_definitions(gmspips
        PRIVATE GMS_PIPS
        PRIVATE GMS_MPI
        PRIVATE LINKCONSTR
        ) # -DGMS_LOG

target_link_libraries(gmspips
        PRIVATE gmspipsio
        PRIVATE pips-ipmpp
        PRIVATE OpenMP::OpenMP_CXX
        PRIVATE MPI::MPI_CXX
        )

############## gmschk ##############
add_executable(gmschk
        Drivers/gmspips/gmspipschk.cpp
        )
target_link_libraries(gmschk
        PRIVATE gmspipsio
        )

########### gmpipsreader ###########
add_library(gmspips_reader Drivers/gmspips/gmspips_reader.cpp)

target_link_libraries(gmspips_reader
        PRIVATE gmspipsio
        PRIVATE pips-ipmpp
        PRIVATE MPI::MPI_CXX
        )
