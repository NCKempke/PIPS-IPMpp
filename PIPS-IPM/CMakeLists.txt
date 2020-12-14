#add_definitions(-DTIMING -DSTOCH_TESTING) # timing output
include_directories(
   Core/Abstract 
   Core/Vector 
   Core/Utilities 
   Core/QpSolvers 
   Core/QpPreprocess 
   Core/QpGen
   Core/SparseLinearAlgebra 
   Core/DenseLinearAlgebra 
   Core/Readers
   Core/LinearSolvers/Ma27Solver 
   Core/LinearSolvers/Ma57Solver 
   Core/LinearSolvers/Mc30Scaler 
   Core/LinearSolvers/Preconditioners
   #Core/LinearSolvers/Ma86Solver 
   Core/LinearSolvers/PardisoSolver/PardisoSolver
   Core/LinearSolvers/PardisoSolver/PardisoIndefSolver
   Core/LinearSolvers/PardisoSolver/PardisoSchurSolver
   Core/LinearSolvers/BiCGStabSolver 
   Core/StochLinearAlgebra 
   Core/QpStoch
)

if( HAVE_MUMPS )
   include_directories(Core/LinearSolvers/MumpsSolver)
endif( HAVE_MUMPS )

add_subdirectory(Core)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_executable(pipsipmCallbackExample Drivers/callbackExample.cpp)
target_compile_options(pipsipmCallbackExample PRIVATE -Wall ${OpenMP_CXX_FLAGS})
target_link_libraries(pipsipmCallbackExample 
   PRIVATE ooqpstoch 
   PRIVATE ooqpstochla 
   PRIVATE ooqpgondziostoch 
   PRIVATE ooqpequistochscaler
   PRIVATE ooqpgensparse 
   PRIVATE ooqpbase 
   PRIVATE ooqpsparse 
   PRIVATE ooqpdense 
   PRIVATE ooqpqpstochpresolver
   PRIVATE OpenMP::OpenMP_CXX
   PRIVATE MPI::MPI_CXX
)
  
if (BUILD_GDX_SOURCE)
   add_definitions(-D_P3_DllExport= -DSTDCALL= -DGDX_CALLCONV=unwanted -DLEG)
   add_library(statgdx
     Drivers/statgdx/clibtypes.cpp
     Drivers/statgdx/datastorage.cpp
     Drivers/statgdx/gdlaudit.cpp
     Drivers/statgdx/gdxdclib.cpp
     Drivers/statgdx/gdxstatic.cpp
     Drivers/statgdx/gmsdata.cpp
     Drivers/statgdx/gmsgen.cpp
     Drivers/statgdx/gmsglob.cpp
     Drivers/statgdx/gmsglobx.cpp
     Drivers/statgdx/gmsheapnew.cpp
     Drivers/statgdx/gmslibname.cpp
     Drivers/statgdx/gmsobj.cpp
     Drivers/statgdx/gmsspecs.cpp
     Drivers/statgdx/gmsstrm.cpp
     Drivers/statgdx/gxdefs.cpp
     Drivers/statgdx/gxfile.cpp
     Drivers/statgdx/idglobal_p3.cpp
     Drivers/statgdx/math_p3.cpp
     Drivers/statgdx/p3Custom2.cpp
     Drivers/statgdx/p3ieeefp.cpp
     Drivers/statgdx/p3io.cpp
     Drivers/statgdx/p3library.cpp
     Drivers/statgdx/p3platform.cpp
     Drivers/statgdx/p3private.cpp
     Drivers/statgdx/p3process.cpp
     Drivers/statgdx/p3threads.cpp
     Drivers/statgdx/p3utils.cpp
     Drivers/statgdx/paldoorg.cpp
     Drivers/statgdx/pchutil.cpp
     Drivers/statgdx/runner.cpp
     Drivers/statgdx/strhash.cpp
     Drivers/statgdx/strutilx.cpp
     Drivers/statgdx/system_p3.cpp
     Drivers/statgdx/sysutils_p3.cpp
     Drivers/statgdx/utils2.cpp
     Drivers/statgdx/xcompress.cpp
     Drivers/statgdx/gclgms.c )
     
   add_executable(gmspips Drivers/gmspips/gmspips.cpp Drivers/gmspips/gmspipsio.c)
   add_definitions(-DGMS_PIPS -DGMS_MPI -DLINKCONSTR -DGDXSOURCE -fno-strict-aliasing -fwrapv ) # -DGMS_LOG
   target_include_directories(gmspips PUBLIC Drivers/statgdx Drivers/gmspips/)
   target_link_libraries(gmspips
       PRIVATE statgdx
       PRIVATE ooqpstoch 
       PRIVATE ooqpstochla
       PRIVATE ooqpgondziostoch
       PRIVATE ooqpequistochscaler
       PRIVATE ooqpgensparse
       PRIVATE ooqpbase
       PRIVATE ooqpsparse
       PRIVATE ooqpdense
       PRIVATE ooqpqpstochpresolver
       PRIVATE OpenMP::OpenMP_CXX 
       PRIVATE MPI::MPI_CXX
   )

   add_executable(gmschk Drivers/gmspips/gmspipschk.cpp Drivers/gmspips/gmspipsio.c)
   target_include_directories(gmschk 
      PUBLIC Drivers/statgdx
      PRIVATE Drivers/gmspips/
   )
   target_link_libraries(gmschk 
      PRIVATE dl
      PRIVATE statgdx
   )
else()
   set(GAMSDIR Drivers/gams)
   add_executable(gmspips Drivers/gmspips/gmspips.cpp Drivers/gmspips/gmspipsio.c ${GAMSDIR}/apifiles/C/api/gmomcc.c ${GAMSDIR}/apifiles/C/api/gevmcc.c ${GAMSDIR}/apifiles/C/api/gdxcc.c)
   add_definitions(-DGMS_PIPS -DGMS_MPI -DLINKCONSTR) # -DGMS_LOG
   target_compile_options(gmspips PRIVATE -Wall ${OpenMP_CXX_FLAGS})
   target_include_directories(gmspips PUBLIC ${GAMSDIR}/apifiles/C/api/ Drivers/gmspips/)
   target_link_libraries(gmspips
       PRIVATE ooqpstoch
       PRIVATE ooqpstochla
       PRIVATE ooqpgondziostoch
       PRIVATE ooqpequistochscaler
       PRIVATE ooqpgensparse
       PRIVATE ooqpbase
       PRIVATE ooqpsparse
       PRIVATE ooqpdense
       PRIVATE ooqpqpstochpresolver
       PRIVATE OpenMP::OpenMP_CXX
       PRIVATE MPI::MPI_CXX
   )
   
   add_executable(gmschk Drivers/gmspips/gmspipschk.cpp Drivers/gmspips/gmspipsio.c ${GAMSDIR}/apifiles/C/api/gmomcc.c ${GAMSDIR}/apifiles/C/api/gevmcc.c ${GAMSDIR}/apifiles/C/api/gdxcc.c)
   target_include_directories(gmschk PUBLIC ${GAMSDIR}/apifiles/C/api/ Drivers/gmspips/)
   target_link_libraries(gmschk PRIVATE dl)
endif(BUILD_GDX_SOURCE)
