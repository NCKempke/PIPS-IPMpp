cmake_minimum_required(VERSION 3.9)

project(PIPSIPMpp VERSION 1.0 LANGUAGES CXX C)

message("Default C++ compiler flags: ${CMAKE_CXX_FLAGS}")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# For debugging build, set CMAKE_BUILD_TYPE to DEBUG
# For optimized build, set CMAKE_BUILD_TYPE to RELEASE
# For optimized build with debug info, set CMAKE_BUILD_TYPE to RELWITHDEBINFO
# CMAKE_BUILD_TYPE should be set at the command line or in a toolchain file

### TODO : this will not work cross compiler / cross platform.. G++ only
### TODO : introduce Wextra/pedantic -Wno-zero-as-null-pointer-const maybe
#add_compile_options( -Wall -Wno-long-long -Wnon-virtual-dtor -Woverloaded-virtual )
set(c_flags -Wall -Wno-long-long)
set(cxx_flags -Wall -Wno-long-long -Wnon-virtual-dtor -Woverloaded-virtual -fopenmp)
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
option(WITH_MAKETEST "Enable 'make test'" ON)
message(STATUS "WITH_MAKETEST=${WITH_MAKETEST}")

#use mkl_pardiso solver instead of schenk solver
option(WITH_MKL_PARDISO "Use mkl_pardiso from Intel for PIPS-IPM" OFF)
message(STATUS "WITH_MKL_PARDISO=${WITH_MKL_PARDISO}")

#use mkl_pardiso solver instead of schenk solver
option(WITH_HIERARCHICAL "Use hierarchical approach for solving" OFF)
message(STATUS "WITH_HIERARCHICAL=${WITH_HIERARCHICAL}")

#use (sparse) MUMPS solver for local Schur complement
option(WITH_MUMPS_LEAF "Use MUMPS (if available) for (local) Schur complement at each leaf" OFF)
message(STATUS "WITH_MUMPS_LEAF=${WITH_MUMPS_LEAF}")

#use (sparse) MUMPS solver for global Schur complement
option(WITH_MUMPS_ROOT "Use MUMPS (if available) for (global) Schur complement at root" OFF)
message(STATUS "WITH_MUMPS_ROOT=${WITH_MUMPS_ROOT}")

#compute local schur complement blockwise with pardiso
option(WITH_PARDISO_BLOCKSC "Use Pardiso (if available) for comuting the Schur complement in blocks" OFF)
message(STATUS "WITH_PARDISO_BLOCKSC=${WITH_PARDISO_BLOCKSC}")

option(BUILD_GDX_SOURCE "Build GDX Source" ON)
message(STATUS "BUILD_GDX_SOURCE=${BUILD_GDX_SOURCE}")

option(USE_MA27 "Use MA27 solver for global and local Schur complement" OFF)
message(STATUS "USE_MA27=${USE_MA27}")

option(USE_MA57 "Use MA57 solver for global and local Schur complement" OFF)
message(STATUS "USE_MA57=${USE_MA57}")

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
   set(CMAKE_CXX_COMPILER "mpicxx")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
else(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   message( FATAL_ERROR "Operating system not supported")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

## required packages 
find_library(GFORTRAN_LIBRARY NAMES libgfortran.so REQUIRED)
find_library(DL_LIBRARY NAMES dl REQUIRED)

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
    set(MATH_LIBS "-Wl,--start-group -Wl,--no-as-needed -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_sequential -lmkl_core -lgfortran -fopenmp -lpthread -lm -Wl,--end-group -std=c++11")
    ####set(MATH_LIBS "-Wl,--start-group -Wl,--no-as-needed -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_avx -lmkl_def -lgomp -lgfortran -lpthread -lm -ldl -Wl,--end-group")
    message(STATUS "Using MKLROOT:$ENV{MKLROOT} MATH_LIBS=${MATH_LIBS}")
elseif(IS_DIRECTORY $ENV{EBROOTICCIFORT})
    message(STATUS "Using Intel Compiler MKL")
    #try mkl_avx2 / 512
    set(MATH_LIBS "-Wl,--start-group  -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_sequential -lmkl_core -liomp5 -lgfortran -qopenmp -lpthread -lm -Wl,--end-group -std=c++11")
    ###set(MATH_LIBS "-Wl,--start-group -Wl,--no-as-needed -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_lapack95_lp64 -lmkl_sequential -lmkl_core -liomp5 -lgfortran -fopenmp -lpthread -lm -Wl,--end-group -std=c++11")
    ###set(MATH_LIBS "-Wl,--start-group  -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_avx -lmkl_def -liomp5 -lpthread -lgfortran -lm -ldl -Wl,--end-group")
    message(STATUS "Using MKLROOT:$ENV{MKLROOT} MATH_LIBS=${MATH_LIBS}")
elseif(IS_DIRECTORY $ENV{MKLROOT})
   set(MATH_LIBS "-Wl,--start-group  -L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group -lgomp -lgfortran -lm -ldl")
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
      message(STATUS "MA57_LIBRARY = ${MA57_LIBRARY}")
      message(STATUS "METIS_4_LIBRARY = ${METIS_4_LIBRARY}")
      set(MA57_LIBRARY ${MA57_LIBRARY} ${METIS_4_LIBRARY})
   endif(NOT METIS_4_LIBRARY)
else(MA57_LIBRARY)
   set(HAVE_MA57 FALSE)
   message(STATUS "Did not find MA57")
endif(MA57_LIBRARY)

########################################
# MUMPS
########################################
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
     
find_library(MUMPS_D_LIBRARY libdmumps.a PATHS ${MUMPS_DIR}/lib NO_DEFAULT_PATH)
find_library(MUMPS_COMMON_LIBRARY libmumps_common.a PATHS ${MUMPS_DIR}/lib NO_DEFAULT_PATH)
find_library(MUMPS_PORD_LIBRARY libpord.a PATHS ${MUMPS_DIR}/lib NO_DEFAULT_PATH)

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
      set(MUMPS_LIBRARY "-Wl,--start-group ${MUMPS_D_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${MUMPS_FORT_LIB} ${MUMPS_PORD_LIBRARY} ${SCALAPACK_LIBRARIES} ${METIS_5_LIBRARY} \
      ${PARMETIS_LIBRARY} -Wl,--end-group")

      message(STATUS "MUMPS_D_LIBRARY = ${MUMPS_D_LIBRARY}")
      message(STATUS "MUMPS_COMMON_LIBRARY = ${MUMPS_COMMON_LIBRARY}")
      message(STATUS "MUMPS_PORD_LIBRARY = ${MUMPS_PORD_LIBRARY}")
      message(STATUS "MUMPS_INCLUDE_DIR = ${MUMPS_INCLUDE_DIR}")
   endif( NOT PARMETIS_LIBRARY OR NOT METIS_5_LIBRARY OR NOT SCALAPACK_LIBRARIES OR NOT MUMPS_FORT_LIB )

   message(STATUS "PARMETIS_LIBRARY = ${PARMETIS_LIBRARY}")
   message(STATUS "METIS_5_LIBRARY = ${METIS_5_LIBRARY}")
   message(STATUS "SCALAPACK_LIBRARIES = ${SCALAPACK_LIBRARIES}")
   message(STATUS "MUMPS_FORT_LIB:=${MUMPS_FORT_LIB}")
else( MUMPS_D_LIBRARY AND MUMPS_COMMON_LIBRARY AND MUMPS_PORD_LIBRARY AND MUMPS_INCLUDE_DIR )
   set(HAVE_MUMPS FALSE)
   message(STATUS "Did not find MUMPS")
endif( MUMPS_D_LIBRARY AND MUMPS_COMMON_LIBRARY AND MUMPS_PORD_LIBRARY AND MUMPS_INCLUDE_DIR )

########################################
# MKL_PARDISO Solver
########################################
message(STATUS "Looking for MKL_PARDISO")

if(IS_DIRECTORY $ENV{MKLROOT})
   find_library(MKL_INTEL_LP64_LIBRARY NAME libmkl_intel_lp64.a HINTS $ENV{MKLROOT}/lib/intel64)
   find_library(MKL_CORE_LIBRARY NAME libmkl_core.a HINTS $ENV{MKLROOT}/lib/intel64)
   find_library(MKL_SEQUENTIAL_LIBRARY NAME libmkl_sequential.a HINTS $ENV{MKLROOT}/lib/intel64)
   find_path(MKL_INCLUDE_DIR $ENV{MKLROOT}/include)
   
   if( !MKL_INTEL_LP64_LIBRARY OR !MKL_CORE_LIBRARY OR !MKL_SEQUENTIAL_LIBRARY OR !MKL_INCLUDE_DIR )
      message(STATUS "Did not find MKL_PARDISO" )
      set( HAVE_MKL_PARDISO FALSE )
   else( !MKL_INTEL_LP64_LIBRARY OR !MKL_CORE_LIBRARY OR !MKL_SEQUENTIAL_LIBRARY OR !MKL_INCLUDE_DIR )
      set( HAVE_MKL_PARDISO TRUE )
      set(MKL_PARDISO_LIBS "-Wl,--start-group ${MKL_INTEL_LP64_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_SEQUENTIAL_LIBRARY} -Wl,--end-group -lgfortran -lpthread -lm")
      message(STATUS "MKL_INTEL_LP64_LIBRARY = ${MKL_INTEL_LP64_LIBRARY}")
      message(STATUS "MKL_CORE_LIBRARY = ${MKL_CORE_LIBRARY}")
      message(STATUS "MKL_SEQUENTIAL_LIBRARY = ${MKL_SEQUENTIAL_LIBRARY}")
      message(STATUS "MKL_INCLUDE_DIR = ${MKL_INCLUDE_DIR}")
   endif( !MKL_INTEL_LP64_LIBRARY OR !MKL_CORE_LIBRARY OR !MKL_SEQUENTIAL_LIBRARY OR !MKL_INCLUDE_DIR )
else(IS_DIRECTORY $ENV{MKLROOT})
   message( STATUS "Did not find MKL_PARDISO")
   set( HAVE_MKL_PARDISO FALSE )
endif(IS_DIRECTORY $ENV{MKLROOT})
   
########################################
# Schenk PARDISO
########################################
message(STATUS "Looking for Schenk PARDISO")

if(DEFINED ENV{PARDISO_DIR})
   set(PARDISO_DIR "$ENV{PARDISO_DIR}")
else(DEFINED ENV{PARDISO_DIR})
   set(PARDISO_DIR ${SHARED_DIR}/PARDISO/src)
endif(DEFINED ENV{PARDISO_DIR})

find_library(PARDISO_LIBRARY NAME libpardiso.so HINTS ${PARDISO_DIR})

if( PARDISO_LIBRARY )
   set(HAVE_PARDISO TRUE)
   message( STATUS "PARDISO_LIBRARY = ${PARDISO_LIBRARY}")
else( PARDISO_LIBRARY )
   set(HAVE_PARDISO FALSE)
   message(STATUS "Did not find Schenk PARDSIO")
endif( PARDISO_LIBRARY )

if( !HAVE_MA27 AND !HAVE_MA57 AND !HAVE_MUMPS AND !HAVE_MKL_PARDISO AND !HAVE_PARDISO )
   message( FATAL_ERROR "Did not find any of the supported linear solvers - cannot build PIPSIPMpp..")
endif( !HAVE_MA27 AND !HAVE_MA57 AND !HAVE_MUMPS AND !HAVE_MKL_PARDISO AND !HAVE_PARDISO )

########################################
# check solver settings
########################################
if( WITH_MKL_PARDISO AND !HAVE_MKL_PARDISO )
   message( WARNING "Specified WITH_MKL_PARDISO but did not find it - using different solver")
   set( WITH_MKL_PARDISO FALSE )
endif( WITH_MKL_PARDISO AND !HAVE_MKL_PARDISO )

if( WITH_MUMPS_LEAF OR WITH_MUMPS_ROOT ) 
   if( NOT HAVE_MUMPS )  
      message( WARNING "Specified MUMPS as linear solver but did not find it - using different solver")
      set( WITH_MUMPS_LEAF FALSE )
      set( WITH_MUMPS_ROOT FALSE )
   endif( NOT HAVE_MUMPS )  
endif( WITH_MUMPS_LEAF OR WITH_MUMPS_ROOT )   
   
if( USE_MA27 AND !HAVE_MA27 )   
   message( WARNING "Specified USE_MA27 but did not find it - using different solver")
   set( USE_MA27 FALSE )
endif( USE_MA27 AND !HAVE_MA27 )   

if( USE_MA57 AND !HAVE_MA57 )   
   message( WARNING "Specified USE_MA57 but did not find it - using different solver")
   set( USE_MA57 FALSE )
endif( USE_MA57 AND !HAVE_MA57 )   
   
########################################
# user specified MA27
########################################
if( USE_MA27 )
   message(STATUS "Using MA27 as linear solver")
   set(WITH_PARDISO_BLOCKSC TRUE)
   
   set( USE_MA57 false )
   set( HAVE_MA57 false )
   
   set( HAVE_MUMPS false )
   set( WITH_MUMPS_LEAF false )
   set( WITH_MUMPS_ROOT false )
   
   set( HAVE_MKL_PARDISO false )
   set( WITH_MKL_PARDISO false )
   
   set( HAVE_PARDISO false )
endif( USE_MA27 )
if( HAVE_MA27 )
   add_definitions(-DWITH_MA27)
endif( HAVE_MA27)

########################################
# user specified MA57
########################################
if( USE_MA57 )
   message(STATUS "Using MA57 as linear solver")
   set(WITH_PARDISO_BLOCKSC TRUE)
   
   set( HAVE_MUMPS false )
   set( WITH_MUMPS_LEAF false )
   set( WITH_MUMPS_ROOT false )
   
   set( HAVE_MKL_PARDISO false )
   set( WITH_MKL_PARDISO false )
   
   set( HAVE_PARDISO false )
endif( USE_MA57 )
if( HAVE_MA57 )
   add_definitions(-DWITH_MA57)
endif( HAVE_MA57)
########################################
# user specified MUMPS
########################################
if( HAVE_MUMPS )
   if( WITH_MUMPS_LEAF )
      message( STATUS "Using MUMPS for leaf computations")
      add_definitions(-DWITH_MUMPS_LEAF)
   endif(WITH_MUMPS_LEAF)
   if(WITH_MUMPS_ROOT)
      message( STATUS "Using MUMPS for root computations")
      add_definitions(-DWITH_MUMPS_ROOT)
   endif(WITH_MUMPS_ROOT)
   
   if( WITH_MUMPS_LEAF AND WITH_MUMPS_ROOT )
      set( HAVE_MKL_PARDISO false )
      set( WITH_MKL_PARDISO false )
      
      set( HAVE_PARDISO false )
   endif( WITH_MUMPS_LEAF AND WITH_MUMPS_ROOT )
endif(HAVE_MUMPS)

if( WITH_MKL_PARDISO )
   message( STATUS "Using MKL_PARDISO as linear solver")
   add_definitions(-DWITH_PARDISO)
   
   set( HAVE_PARDISO false )
endif( WITH_MKL_PARDISO )

if( HAVE_PARDISO )
   message( STATUS "Using Schenk PARDISO as default linear solver")
   add_definitions(-DWITH_PARDISO)
endif(HAVE_PARDISO)

if( WITH_HIERARCHICAL )
   message(STATUS "Using Hierarchical approach")
   add_definitions(-DHIERARCHICAL)
   SET( WITH_PARDISO_BLOCKSC TRUE )
endif(WITH_HIERARCHICAL)

if(WITH_PARDISO_BLOCKSC)
   message(STATUS "Blockwise Schur Complement computation enabled")
   add_definitions( -DPARDISO_BLOCKSC)
else(WITH_PARDISO_BLOCKSC)
   message(STATUS "Augmented approach for Schur Complement computation")
endif(WITH_PARDISO_BLOCKSC)

###############################################################################
###############################################################################

add_subdirectory(Utilities)
add_subdirectory(PIPS-IPM)

##########################################################
# CMake Tests
##########################################################
if(WITH_MAKETEST)
   enable_testing()

   add_test(NAME PIPS-IPM-linkingConsTest COMMAND sh ${PROJECT_SOURCE_DIR}/PIPS-IPM/Test/pipsipmLinkConsTest.sh $<TARGET_FILE:pipsipmCallbackExample> )
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

