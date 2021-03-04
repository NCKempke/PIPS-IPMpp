project(PIPSIPMpp LANGUAGES CXX C)

cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

find_package(MPI REQUIRED)
message(STATUS "MPI_C_COMPILER=${MPI_C_COMPILER}")
message(STATUS "MPI_C_WORKS=${MPI_C_WORKS}")
message(STATUS "MPI_C_LIBRARIES=${MPI_C_LIBRARIES}")
message(STATUS "MPI_C_INCLUDE_DIRS=${MPI_C_INCLUDE_DIRS}")
message(STATUS "MPI_C_LINK_FLAGS=${MPI_C_LINK_FLAGS}")
message(STATUS "MPI_CXX_COMPILER=${MPI_CXX_COMPILER}")
message(STATUS "MPI_CXX_WORKS=${MPI_CXX_WORKS}")
message(STATUS "MPI_CXX_LIBRARIES=${MPI_CXX_LIBRARIES}")
message(STATUS "MPI_CXX_INCLUDE_DIRS=${MPI_CXX_INCLUDE_DIRS}")
message(STATUS "MPI_CXX_LINK_FLAGS=${MPI_CXX_LINK_FLAGS}")
message(STATUS "MPIEXEC=${MPIEXEC}")
message(STATUS "MPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}")


include_directories(SYSTEM ${MPI_CXX_INCLUDE_DIRS})

# WATCH OUT : these iwll be added to all executables - also the GAMS part - this might or might not lead to many warnings
# add gmspips/pips specific WARNINGS via target_compile_options in PIPS-IPM/CMakeLists.txt
 
### TODO : introduce Wextra/pedantic -Wconversion -Wno-zero-as-null-pointer-const maybe
#add_compile_options( -Wall -Wno-long-long -Wnon-virtual-dtor -Woverloaded-virtual )
set(c_flags -Wall -g -Wno-long-long -fPIC)
set(cxx_flags -Wall -g -Wno-long-long -Wnon-virtual-dtor -Wno-zero-as-null-pointer-constant -Woverloaded-virtual -fopenmp -fPIC)
add_compile_options(
  "$<$<COMPILE_LANGUAGE:C>:${c_flags}>"
  "$<$<COMPILE_LANGUAGE:CXX>:${cxx_flags}>"
)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_BUILD_TYPE)
  message("CMake build type is ${CMAKE_BUILD_TYPE}")
else()
  message("CMAKE_BUILD_TYPE not set - using DEBUG as default!")
  set( CMAKE_BUILD_TYPE "Debug" )
endif(CMAKE_BUILD_TYPE)

#with dumping matrices 1st stage (and global soon)
if(DUMP)
  add_definitions(-DDUMP)
endif()

#with testing of (selected) executables
option(WITH_MAKETEST "Enable 'make test'" OFF)
message(STATUS "WITH_MAKETEST=${WITH_MAKETEST}")

option(BUILD_GDX_SOURCE "Build GDX Source" OFF)
message(STATUS "BUILD_GDX_SOURCE=${BUILD_GDX_SOURCE}")

option(WITH_MKL_PARDISO "Building with MKL PARDISO (if available) instead of Schenk PARDISO" OFF)
message(STATUS "WITH_MKL_PARDISO=${WITH_MKL_PARDISO}")

option(WITH_TIMING "Build with timing recording and reporting for PIPS-IPM" OFF)
message(STATUS "WITH_TIMING=${WITH_TIMING}")

option(WITH_VERBOSE "Build with extra verbosity level" OFF)
message(STATUS "WITH_VERBOSE=${WITH_VERBOSE}")

if(WITH_TIMING)
	add_definitions(-DTIMING -DSTOCH_TESTING -DNLPTIMING)
endif(WITH_TIMING)

if(WITH_VERBOSE)
	add_definitions(-DVERBOSE)
endif(WITH_VERBOSE)

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   #set(CMAKE_CXX_COMPILER "mpicxx")
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
else(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   message( FATAL_ERROR "Operating system not supported")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

## required packages 
find_library(GFORTRAN_LIBRARY NAMES libgfortran.so REQUIRED)
find_library(DL_LIBRARY NAMES dl REQUIRED)

if( BUILD_GDX_SOURCE )
   find_package(Threads REQUIRED)
endif( BUILD_GDX_SOURCE )

find_package(OpenMP)
if(NOT OpenMP_FOUND)
   message(FATAL_ERROR "Configuration aborted due to missing OpenMP libraries!")
endif(NOT OpenMP_FOUND)
message(STATUS "OpenMP ... found")

if(MATH_LIBS)
	message(STATUS "Using user-defined MATH_LIBS=${MATH_LIBS}")
elseif(IS_DIRECTORY $ENV{EBROOTGCC})
    message(STATUS "Using GNU Compiler MKL")
    #try mkl_avx2 /512
    set(MATH_LIBS "-Wl,--start-group -Wl,--no-as-needed -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_sequential -lmkl_core -lgfortran -lpthread -lm -ldl -Wl,--end-group -std=c++11")
    ## mkl advisor                                               -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

    ####set(MATH_LIBS "-Wl,--start-group -Wl,--no-as-needed -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_avx -lmkl_def -lgomp -lgfortran -lpthread -lm -ldl -Wl,--end-group")
    message(STATUS "Using MKLROOT:$ENV{MKLROOT} MATH_LIBS=${MATH_LIBS}")
elseif(IS_DIRECTORY $ENV{EBROOTICCIFORT})
    message(STATUS "Using Intel Compiler MKL")
    #try mkl_avx2 / 512
    set(MATH_LIBS "-Wl,--start-group -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_sequential -lmkl_core -liomp5 -lgfortran -lpthread -lm -ldl -Wl,--end-group -std=c++11")
    ###set(MATH_LIBS "-Wl,--start-group -Wl,--no-as-needed -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_sequential -lmkl_core -liomp5 -lgfortran -fopenmp -lpthread -lm -Wl,--end-group -std=c++11")
    ###set(MATH_LIBS "-Wl,--start-group  -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_avx -lmkl_def -liomp5 -lpthread -lgfortran -lm -ldl -Wl,--end-group")
    message(STATUS "Using MKLROOT:$ENV{MKLROOT} MATH_LIBS=${MATH_LIBS}")
elseif(IS_DIRECTORY $ENV{MKLROOT})
   set( MATH_LIBS  "-Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_gnu_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -m64" )
   # set(MATH_LIBS "-Wl,--start-group  -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group -lgomp -lgfortran -lm -ldl")
   message(STATUS "Using MKLROOT:$ENV{MKLROOT} MATH_LIBS=${MATH_LIBS}")
else(MATH_LIBS)
    message(STATUS "Searching system path for MATH_LIBS")
    enable_language(Fortran)
    message(STATUS "FC = ${FC}")
    message(STATUS "Fortran compiler name is ${CMAKE_Fortran_COMPILER}")
    message(STATUS "Fortran compiler ID is ${CMAKE_Fortran_COMPILER_ID}")
    find_package(LAPACK REQUIRED)
    message(STATUS " LAPACK_LIBRARIES:  ${LAPACK_LIBRARIES}")
    message(STATUS " CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES:  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
    set(MATH_LIBS ${LAPACK_LIBRARIES} -ldl -l${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ${OpenMP_CXX_FLAGS})
    message(STATUS " MATH_LIBS:  ${MATH_LIBS}")
    find_package(BLAS REQUIRED)
endif(MATH_LIBS)

#Boost
find_package(Boost COMPONENTS log log_setup regex system filesystem thread chrono date_time)
message(STATUS "Boost ... found")

if(NOT Boost_FOUND)
   message("! A Boost version containing BOOST.LOG has NOT been found.")
   message("! Please specify the location of the Boost library containing 'log' package.")
   message("! For this create a 'Toolchain.cmake' file or edit the one from the root ")
   message("! directory and use 'set(BOOST_ROOT \"/path/to/boost\")'.")
   message("! Then run  'cmake -DCMAKE_TOOLCHAIN_FILE=../Toolchain.cmake ..' in the build directory.")
   message(FATAL_ERROR "Configuration aborted due to missing BOOST libraries!")
endif(NOT Boost_FOUND)

########################################
# ThirdPartyLibs
########################################
set(SHARED_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLibs)

set(SOLVERS_FOUND "")
########################################
# MA27
########################################
message(STATUS "Looking for MA27")
if(DEFINED ENV{MA27_DIR})
    set(MA27_DIR "$ENV{MA27_DIR}")
else(DEFINED ENV{MA27_DIR})
    set(MA27_DIR ${SHARED_DIR}/MA27/src)
endif(DEFINED ENV{MA27_DIR})

########################################
# MC30 (used with MA27)
########################################
if(DEFINED ENV{MC30_DIR})
    set(MC30_DIR "$ENV{MC30_DIR}")
else(DEFINED ENV{MC30_DIR})
    set(MC30_DIR ${SHARED_DIR}/MC30/src)
endif(DEFINED ENV{MC30_DIR})

find_library(MA27_LIBRARY libma27.a PATHS ${MA27_DIR}/lib)

if(MA27_LIBRARY)
   find_library(MC30_LIBRARY libmc30.a PATHS ${MC30_DIR}/lib)
    
   if(NOT MC30_LIBRARY)
      message(STATUS "Found MA27 but did not find MC30 which is required - will not use MA27")
      set(HAVE_MA27 FALSE)
   else(NOT MC30_LIBRARY)
      set(HAVE_MA27 TRUE)
      set(SOLVERS_FOUND "${SOLVERS_FOUND} MA27")
      message(STATUS "MA27_LIBRARY = ${MA27_LIBRARY}")
      message(STATUS "MC30_LIBRARY = ${MC30_LIBRARY}")
      set(MA27_LIBRARY ${MA27_LIBRARY} ${MC30_LIBRARY})
   endif(NOT MC30_LIBRARY)
else(MA27_LIBRARY)
   set(HAVE_MA27 FALSE)
   message(STATUS "Did not find MA27")
endif(MA27_LIBRARY)


########################################
# MA57
########################################
message(STATUS "Looking for MA57")
if(DEFINED ENV{MA57_DIR})
    set(MA57_DIR "$ENV{MA57_DIR}")
else(DEFINED ENV{MA57_DIR})
    set(MA57_DIR ${SHARED_DIR}/MA57/src)
endif(DEFINED ENV{MA57_DIR})

########################################
# METIS 4 (used with MA57)
########################################
if(DEFINED ENV{METIS_DIR})
    set(METIS_DIR "$ENV{METIS_DIR}")
else(DEFINED ENV{METIS_DIR})
    set(METIS_DIR ${SHARED_DIR}/METIS_4)
endif(DEFINED ENV{METIS_DIR})

find_library(MA57_LIBRARY libma57.a PATHS ${MA57_DIR}/lib)

if(MA57_LIBRARY)
   find_library(METIS_4_LIBRARY NAME libmetis.a HINTS ${METIS_DIR}/src PATHS ${METIS_DIR})
    
   if(NOT METIS_4_LIBRARY)
      message(STATUS "Found MA57 but did not find METIS 4 which is required - will not use MA57")
      set(HAVE_MA57 FALSE)
   else(NOT METIS_4_LIBRARY)
      set(HAVE_MA57 TRUE)
      set(SOLVERS_FOUND "${SOLVERS_FOUND} MA57")
      message(STATUS "MA57_LIBRARY = ${MA57_LIBRARY}")
      message(STATUS "METIS_4_LIBRARY = ${METIS_4_LIBRARY}")
   endif(NOT METIS_4_LIBRARY)
else(MA57_LIBRARY)
   set(HAVE_MA57 FALSE)
   message(STATUS "Did not find MA57")
endif(MA57_LIBRARY)

########################################
# MUMPS
########################################
if( FALSE )
message(STATUS "Looking for MUMPS")
if(DEFINED ENV{MUMPS_DIR})
   set(MUMPS_DIR "$ENV{MUMPS_DIR}")
else(DEFINED ENV{MUMPS_DIR})
   set(MUMPS_DIR ${SHARED_DIR}/MUMPS)
endif(DEFINED ENV{MUMPS_DIR})

########################################
# PARMETIS and METIS 5 (used with MUMPS)
########################################

if(DEFINED ENV{PARMETIS_DIR})
   set(PARMETIS_DIR "$ENV{PARMETIS_DIR}")
else(DEFINED ENV{PARMETIS_DIR})
   set(PARMETIS_DIR ${SHARED_DIR}/METIS)
endif(DEFINED ENV{PARMETIS_DIR})
     
find_library(MUMPS_D_LIBRARY libdmumps.so PATHS ${MUMPS_DIR}/lib NO_DEFAULT_PATH)
find_library(MUMPS_COMMON_LIBRARY libmumps_common.so PATHS ${MUMPS_DIR}/lib NO_DEFAULT_PATH)
find_library(MUMPS_PORD_LIBRARY libpord.so PATHS ${MUMPS_DIR}/lib NO_DEFAULT_PATH)

find_path(MUMPS_INCLUDE_DIR NAME dmumps_c.h HINTS ${MUMPS_DIR}/include NO_DEFAULT_PATH)

if( MUMPS_D_LIBRARY AND MUMPS_COMMON_LIBRARY AND MUMPS_PORD_LIBRARY AND MUMPS_INCLUDE_DIR )
   
   find_library(PARMETIS_LIBRARY NAME libparmetis.a HINTS ${PARMETIS_DIR}/lib PATHS ${PARMETIS_DIR}/lib)
   find_library(METIS_5_LIBRARY NAME libmetis.a HINTS ${PARMETIS_DIR}/lib PATHS ${PARMETIS_DIR}/lib)
   if(NOT SCALAPACK_LIBRARIES)
      MESSAGE(STATUS "SCALAPACK_LIBRARIES for MUMPS not set. Try default linking to scalapack-openmpi.")
      find_library(SCALAPACK_OMPI_LIBRARY NAMES scalapack-openmpi)
      set(SCALAPACK_LIBRARIES ${SCALAPACK_OMPI_LIBRARY})
   endif(NOT SCALAPACK_LIBRARIES)
   
   if( NOT MUMPS_FORT_LIB )
      MESSAGE(STATUS "MUMPS_FORT_LIB not set for linking the Fortran MPI interface library required by MUMPS. Try using OpenMPI Fortran interface -lmpi_mpifh")
      find_library(MUMPS_FORT_LIB mpi_mpifh)
   endif(NOT MUMPS_FORT_LIB)
   
   if( NOT PARMETIS_LIBRARY OR NOT METIS_5_LIBRARY OR NOT SCALAPACK_LIBRARIES OR NOT MUMPS_FORT_LIB )
      message(STATUS "Found MUMPS but did not find METIS 5 and or PARMETIS or SCALAPACK or MUMPS_FORT_LIB which are required - will not use MUMPS")
      set(HAVE_MUMPS FALSE)
   else( NOT PARMETIS_LIBRARY OR NOT METIS_5_LIBRARY OR NOT SCALAPACK_LIBRARIES OR NOT MUMPS_FORT_LIB )
      set(HAVE_MUMPS TRUE)
      set(SOLVERS_FOUND "${SOLVERS_FOUND} MUMPS")
      #set(MUMPS_LIBRARY "-Wl,--start-group ${MUMPS_D_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${MUMPS_FORT_LIB} ${MUMPS_PORD_LIBRARY} ${SCALAPACK_LIBRARIES} ${METIS_5_LIBRARY} \
      #${PARMETIS_LIBRARY} -Wl,--end-group")
      message(STATUS "MUMPS_D_LIBRARY = ${MUMPS_D_LIBRARY}")
      message(STATUS "MUMPS_COMMON_LIBRARY = ${MUMPS_COMMON_LIBRARY}")
      message(STATUS "MUMPS_PORD_LIBRARY = ${MUMPS_PORD_LIBRARY}")
      message(STATUS "MUMPS_INCLUDE_DIR = ${MUMPS_INCLUDE_DIR}")
   endif( NOT PARMETIS_LIBRARY OR NOT METIS_5_LIBRARY OR NOT SCALAPACK_LIBRARIES OR NOT MUMPS_FORT_LIB )

   message(STATUS "MUMPS_FORT_LIB = ${MUMPS_FORT_LIB}")
   message(STATUS "PARMETIS_LIBRARY = ${PARMETIS_LIBRARY}")
   message(STATUS "METIS_5_LIBRARY = ${METIS_5_LIBRARY}")
   message(STATUS "SCALAPACK_LIBRARIES = ${SCALAPACK_LIBRARIES}")
   message(STATUS "MUMPS_FORT_LIB:=${MUMPS_FORT_LIB}")
else( MUMPS_D_LIBRARY AND MUMPS_COMMON_LIBRARY AND MUMPS_PORD_LIBRARY AND MUMPS_INCLUDE_DIR )
   set(HAVE_MUMPS FALSE)
   message(STATUS "Did not find MUMPS")
endif( MUMPS_D_LIBRARY AND MUMPS_COMMON_LIBRARY AND MUMPS_PORD_LIBRARY AND MUMPS_INCLUDE_DIR )
else( FALSE )
   set(HAVE_MUMPS FALSE)
endif( FALSE )
########################################
# MKL_PARDISO Solver
########################################
if( WITH_MKL_PARDISO ) 
   message(STATUS "Looking for MKL_PARDISO")

   if( IS_DIRECTORY $ENV{MKLROOT} )
      find_library(MKL_INTEL_LP64_LIBRARY NAME libmkl_intel_lp64.a HINTS $ENV{MKLROOT}/lib/intel64)
      find_library(MKL_CORE_LIBRARY NAME libmkl_core.a HINTS $ENV{MKLROOT}/lib/intel64)
      find_library(MKL_SEQUENTIAL_LIBRARY NAME libmkl_sequential.a HINTS $ENV{MKLROOT}/lib/intel64)
      find_path(MKL_INCLUDE_DIR $ENV{MKLROOT}/include)
      
      if( !MKL_INTEL_LP64_LIBRARY OR !MKL_CORE_LIBRARY OR !MKL_SEQUENTIAL_LIBRARY OR !MKL_INCLUDE_DIR )
         message(STATUS "Did not find MKL_PARDISO" )
         set( HAVE_MKL_PARDISO FALSE )
      else( !MKL_INTEL_LP64_LIBRARY OR !MKL_CORE_LIBRARY OR !MKL_SEQUENTIAL_LIBRARY OR !MKL_INCLUDE_DIR )
         set( HAVE_MKL_PARDISO TRUE )
         set(SOLVERS_FOUND "${SOLVERS_FOUND} MKL_PARDISO")
         set(MKL_PARDISO_LIBS "-Wl,--start-group ${MKL_INTEL_LP64_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_SEQUENTIAL_LIBRARY} -Wl,--end-group -lgfortran -lpthread -lm")
         message(STATUS "MKL_INTEL_LP64_LIBRARY = ${MKL_INTEL_LP64_LIBRARY}")
         message(STATUS "MKL_CORE_LIBRARY = ${MKL_CORE_LIBRARY}")
         message(STATUS "MKL_SEQUENTIAL_LIBRARY = ${MKL_SEQUENTIAL_LIBRARY}")
         message(STATUS "MKL_INCLUDE_DIR = ${MKL_INCLUDE_DIR}")
      endif( !MKL_INTEL_LP64_LIBRARY OR !MKL_CORE_LIBRARY OR !MKL_SEQUENTIAL_LIBRARY OR !MKL_INCLUDE_DIR )
   else( IS_DIRECTORY $ENV{MKLROOT} )
      message( STATUS "Did not find MKL_PARDISO")
      set( HAVE_MKL_PARDISO FALSE )
   endif( IS_DIRECTORY $ENV{MKLROOT} )
else( WITH_MKL_PARDISO )
   message( STATUS "NOT looking for MKL_PARDISO")
   set( HAVE_MKL_PARDISO FALSE )
endif( WITH_MKL_PARDISO )

########################################
# Schenk PARDISO
########################################
if( NOT WITH_MKL_PARDISO )
   message(STATUS "Looking for Schenk PARDISO")

   if(DEFINED ENV{PARDISO_DIR})
      set(PARDISO_DIR "$ENV{PARDISO_DIR}")
   else(DEFINED ENV{PARDISO_DIR})
      set(PARDISO_DIR ${SHARED_DIR}/PARDISO/src)
   endif(DEFINED ENV{PARDISO_DIR})

   find_library(PARDISO_LIBRARY NAME libpardiso.so HINTS ${PARDISO_DIR})

   if( PARDISO_LIBRARY )
      set(HAVE_PARDISO TRUE)
      set(SOLVERS_FOUND "${SOLVERS_FOUND} SCHENK_PARDISO")
      message( STATUS "PARDISO_LIBRARY = ${PARDISO_LIBRARY}")
   else( PARDISO_LIBRARY )
      set(HAVE_PARDISO FALSE)
      message(STATUS "Did not find Schenk PARDSIO")
   endif( PARDISO_LIBRARY )
else( NOT WITH_MKL_PARDISO )
   message(STATUS "NOT looking for Schenk PARDISO")
   set(HAVE_PARDISO FALSE)
endif( NOT WITH_MKL_PARDISO )
   

########################################
# print solvers found or abort if none
########################################
if( !HAVE_MA27 AND !HAVE_MA57 AND !HAVE_MUMPS AND !HAVE_MKL_PARDISO AND !HAVE_PARDISO )
   message( FATAL_ERROR "Did not find any of the supported linear solvers - cannot build PIPSIPMpp..")
else( !HAVE_MA27 AND !HAVE_MA57 AND !HAVE_MUMPS AND !HAVE_MKL_PARDISO AND !HAVE_PARDISO )
   message( STATUS "Found the following solvers: ${SOLVERS_FOUND}")
endif( !HAVE_MA27 AND !HAVE_MA57 AND !HAVE_MUMPS AND !HAVE_MKL_PARDISO AND !HAVE_PARDISO )

if( HAVE_MA27 )
   add_definitions( -DWITH_MA27 )
endif( HAVE_MA27 )

if( HAVE_MA57 )
   add_definitions( -DWITH_MA57 )
endif( HAVE_MA57 ) 

if( HAVE_MUMPS )
   add_definitions( -DWITH_MUMPS )
endif( HAVE_MUMPS )

if( HAVE_MKL_PARDISO )
   add_definitions( -DWITH_MKL_PARDISO )
endif( HAVE_MKL_PARDISO )

if( HAVE_PARDISO )
   add_definitions( -DWITH_PARDISO )
endif( HAVE_PARDISO )
###############################################################################
###############################################################################

add_subdirectory(Utilities)
add_subdirectory(PIPS-IPM)

##########################################################
# CMake Tests
##########################################################
if(WITH_MAKETEST)
   enable_testing()

   include(GoogleTest)
   add_subdirectory(PIPS-IPM/Test/)
   
   #add_test(NAME PIPS-IPM-linkingConsTest COMMAND sh ${PROJECT_SOURCE_DIR}/PIPS-IPM/Test/pipsipmLinkConsTest.sh $<TARGET_FILE:pipsipmCallbackExample> )
endif(WITH_MAKETEST)


get_directory_property( DirDefs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )

foreach( d ${DirDefs} )
    message( STATUS "Found Define: " ${d} )
endforeach()

message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_STANDARD_LIBRIES: ${CMAKE_CXX_STANDARD_LIBRARIES}")
message(STATUS "CMAKE_LIBRARY_PATH: ${CMAKE_LIBRARY_PATH}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message(STATUS "${inc_dirs}")
get_property(link_dirs DIRECTORY PROPERTY LINK_DIRECTORIES)
message(STATUS "${link_dirs}")

