project(PIPSAll)

cmake_minimum_required(VERSION 3.4)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# For debugging build, set CMAKE_BUILD_TYPE to DEBUG
# For optimized build, set CMAKE_BUILD_TYPE to RELEASE
# For optimized build with debug info, set CMAKE_BUILD_TYPE to RELWITHDEBINFO
# CMAKE_BUILD_TYPE should be set at the command line or in a toolchain file
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -lgfortran -Wno-long-long -Wnon-virtual-dtor -Woverloaded-virtual")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_BUILD_TYPE  RELEASE)
if(CMAKE_BUILD_TYPE)
  message("CMake build type is ${CMAKE_BUILD_TYPE}")
  message("Build type default flags are: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
else()
  message("CMAKE_BUILD_TYPE not set!")
endif(CMAKE_BUILD_TYPE)

message(STATUS "Building shared libs: ${BUILD_SHARED_LIBS}")

#with dumping matrices 1st stage (and global soon)
if(DUMP)
  add_definitions(-DDUMP)
endif()

#with testing of (selected) executables
option(WITH_MAKETEST "Enable 'make test'" ON)

#use mkl_pardiso solver instead of schenk solver
option(WITH_MKL_PARDISO "Use mkl_pardiso from Intel for PIPS-IPM" OFF)
message(STATUS "WITH_MKL_PARDISO=${WITH_MKL_PARDISO}")

#use mkl_pardiso solver instead of schenk solver
option(WITH_HIERARCHICAL "Use hierarchical approach for solving" OFF)
message(STATUS "WITH_HIERARCHICAL=${WITH_HIERARCHICAL}")

#use (sparse) Pardiso solver for summed Schur complement
option(WITH_PARDISOINDEF "Use Pardiso (if available) for summed Schur complement" ON)
message(STATUS "WITH_PARDISOINDEF=${WITH_PARDISOINDEF}")

#use (sparse) MUMPS solver for local Schur complement
option(WITH_MUMPS_LEAF "Use MUMPS (if available) for (local) Schur complement at each leaf" OFF)
message(STATUS "WITH_MUMPS_LEAF=${WITH_MUMPS_LEAF}")

#use (sparse) MUMPS solver for global Schur complement
option(WITH_MUMPS_ROOT "Use MUMPS (if available) for (global) Schur complement at root" OFF)
message(STATUS "WITH_MUMPS_ROOT=${WITH_MUMPS_ROOT}")

#compute local schur complement blockwise with pardiso
option(WITH_PARDISO_BLOCKSC "Use Pardiso (if available) for comuting the Schur complement in blocks" OFF)
message(STATUS "WITH_PARDISO_BLOCKSC=${WITH_PARDISO_BLOCKSC}")

option(BUILD_GDX_SOURCE "Build GDX Source" OFF)
message(STATUS "BUILD_GDX_SOURCE=${BUILD_GDX_SOURCE}")

option(USE_MA27 "Use MA27 solver for global and local Schur complement" OFF)
message(STATUS "USE_MA27=${USE_MA27}")

option(USE_MA57 "Use MA57 solver for global and local Schur complement" OFF)
message(STATUS "USE_MA57=${USE_MA57}")

#with timing for PIPS-IPM
option(WITH_TIMING "Build with timing recording and reporting for PIPS-IPM" OFF)
#with PIPS IPM AND NLP profile metric
option(WITH_VERBOSE "Build with extra verbosity level" OFF)

if(WITH_TIMING)
	add_definitions(-DTIMING -DSTOCH_TESTING -DNLPTIMING)
endif(WITH_TIMING)

if(WITH_VERBOSE)
	add_definitions(-DVERBOSE)
endif(WITH_VERBOSE)

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   set(CMAKE_CXX_COMPILER "mpicxx")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -fpermissive")
else(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
   message( FATAL_ERROR "Operating system not supported")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

## required packages 
find_library(GFORTRAN_LIBRARY NAMES libgfortran.so REQUIRED)
find_library(DL_LIBRARY NAMES dl REQUIRED)

find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

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

# include different "whole archive" linking options depending on compiler
#if (CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
#  set(WHOLE_ARCHIVE "-Wl,-all_load")
#  set(NO_WHOLE_ARCHIVE "-Wl,-noall_load")
#  message(STATUS "SETTING HAVE_CLANG")
#  set(HAVE_CLANG TRUE)

  # Clang doesn't have a Fortran compiler in its suite (yet),
  # so detect libraries for gfortran; we need equivalents to
  # libgfortran and libquadmath, which are implicitly
  # linked by flags in CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES
#  find_package(GFortranLibs REQUIRED)

  # Add directory containing libgfortran and libquadmath to
  # linker. Should also contain libgomp, if not using
  # Intel OpenMP runtime
#  link_directories(${GFORTRAN_LIBRARIES_DIR})
#  if((BUILD_PIPS_IPM OR BUILD_PIPS_NLP) AND NOT OpenMP_FOUND)
#    message(STATUS "OpenMP not found!")

    # Try to find Intel OpenMP library first, since this
    # OpenMP runtime is the default in Clang 3.7
#    find_package(LibIomp)
#    if(LIBIOMP_FOUND)
#      if(APPLE)
	# Apple clang's -fopenmp flag only works with -lgomp
#	set(OpenMP_CXX_FLAGS "-liomp5")
#      else()
	# LLVM clang's trunk works with libiomp
#	set(OpenMP_CXX_FLAGS "-fopenmp -liomp5")
#      endif()
#      set(OpenMP_C_FLAGS "${OpenMP_CXX_FLAGS}")
#      include_directories(${LIBIOMP_INCLUDE_DIR})
#      link_directories(${LIBIOMP_LIBRARIES_DIR})
#      message(STATUS "Intel OpenMP libraries directory is ${LIBIOMP_LIBRARIES_DIR}")

#    else() # Otherwise, fall back to GOMP (gcc OpenMP) runtime
#      set(OpenMP_CXX_FLAGS "-fopenmp -lgomp")
#      set(OpenMP_C_FLAGS "${OpenMP_CXX_FLAGS}")
#      include_directories(${LIBGOMP_INCLUDE_DIR})
#    endif()
# endif()
#else ()
 # set(HAVE_CLANG FALSE)
 # set(WHOLE_ARCHIVE "-Wl,--whole-archive")
  #set(NO_WHOLE_ARCHIVE "-Wl,--no-whole-archive")
#endif (CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")

#find_library(FL_LIBRARY NAMES fl PATHS /usr/lib /usr/local/lib /usr/local/opt/flex)

### ThirdPartyLibs
set(SHARED_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLibs)

# HSL Solver and METIS_4
# MA27
if(DEFINED ENV{MA27_DIR})
    set(MA27_DIR "$ENV{MA27_DIR}")
else(DEFINED ENV{MA27_DIR})
    set(MA27_DIR ${SHARED_DIR}/MA27/src)
endif(DEFINED ENV{MA27_DIR})

# MC30
if(DEFINED ENV{MC30_DIR})
    set(MC30_DIR "$ENV{MC30_DIR}")
else(DEFINED ENV{MC30_DIR})
    set(MC30_DIR ${SHARED_DIR}/MC30/src)
endif(DEFINED ENV{MC30_DIR})

# MA60
#set(MA60_DIR ${SHARED_DIR}/MA60/src)
  
# MA57
if(DEFINED ENV{MA57_DIR})
    set(MA57_DIR "$ENV{MA57_DIR}")
else(DEFINED ENV{MA57_DIR})
    set(MA57_DIR ${SHARED_DIR}/MA57/src)
endif(DEFINED ENV{MA57_DIR})

find_library(MA27_LIBRARY libma27.a PATHS ${MA27_DIR}/lib)
find_library(MA57_LIBRARY libma57.a PATHS ${MA57_DIR}/lib)

if(MA27_LIBRARY)
   find_library(MC30_LIBRARY libmc30.a PATHS ${MC30_DIR}/lib)
   #find_library(MA60_LIBRARY libma60.a PATHS ${MA60_DIR}/lib)
    
   if(NOT MC30_LIBRARY)
      message(STATUS "Found MA27 but not MC30 which is required - will not use MA27")
      set(HAVE_MA27 FALSE)
   #elseif(NOT MA60_LIBRARY)
   #   message(STATUS "Found MA27 but not MA60 which is required - will not use MA27")
   #   set(HAVE_MA27 FALSE)
   else(NOT MC30_LIBRARY)
      set(HAVE_MA27 TRUE)
      message(STATUS "MA27_LIBRARY = ${MA27_LIBRARY}")
   #   message(STATUS "MA60_LIBRARY = ${MA60_LIBRARY}")
      message(STATUS "MC30_LIBRARY = ${MC30_LIBRARY}")
      set(MA27_LIBRARY ${MA27_LIBRARY} ${MC30_LIBRARY})# ${MA60_LIBRARY})
   endif(NOT MC30_LIBRARY)
endif(MA27_LIBRARY)
  
if(NOT MA27_LIBRARY AND NOT MA57_LIBRARY)
  message(FATAL_ERROR "Could not find MA27 nor MA57. We need one HSL solver")
endif(NOT MA27_LIBRARY AND NOT MA57_LIBRARY)

# METIS
if(DEFINED ENV{METIS_DIR})
    set(METIS_DIR "$ENV{METIS_DIR}")
else(DEFINED ENV{METIS_DIR})
    set(METIS_DIR ${SHARED_DIR}/METIS_4)
endif(DEFINED ENV{METIS_DIR})

# PARMETIS
set(PARMETIS_DIR ${SHARED_DIR}/METIS)
  
find_library(PARMETIS_LIBRARY NAME libparmetis.a HINTS ${PARMETIS_DIR}/lib PATHS ${PARMETIS_DIR}/lib)
find_library(METIS_5_LIBRARY NAME libmetis.a HINTS ${PARMETIS_DIR}/lib PATHS ${PARMETIS_DIR}/lib)
  
find_library(METIS_4_LIBRARY NAME libmetis.a HINTS ${METIS_DIR}/src PATHS ${METIS_DIR})
if( METIS_4_LIBRARY )
   message(STATUS "METIS_4_LIBRARY = ${METIS_4_LIBRARY}")
   set(HAVE_METIS TRUE)
else( METIS_4_LIBRARY )
   message(FATAL_ERROR "METIS not found! METIS is required for PIPS-IPM. Please run ./wgetMETIS.sh in ${METIS_DIR}")
endif( METIS_4_LIBRARY )
  
if( PARMETIS_LIBRARY AND METIS_5_LIBRARY )
   set( HAVE_PARMETIS TRUE )
else( PARMETIS_LIBRARY AND METIS_5_LIBRARY )
   message(STATUS "PARMETIS not found! Parmetis is required for build PIPS-IPM with MUMPS. Please run ./wgetMETIS.sh in ${PARMETIS_DIR}")
endif( PARMETIS_LIBRARY AND METIS_5_LIBRARY )
  
if(MA57_LIBRARY)
   message(STATUS "MA57_LIBRARY = ${MA57_LIBRARY}")
   include_directories(${METIS_DIR}/Lib)
   set(HAVE_MA57 TRUE)
else(MA57_LIBRARY)
   set(HAVE_MA57 FALSE)
endif(MA57_LIBRARY)

  #MA86 is not needed at this time
  #find_library(MA86_LIBRARY hsl_ma86 PATHS ${MA86_DIR}/lib)

if(WITH_MUMPS_LEAF OR WITH_MUMPS_ROOT)
   set(MUMPS_DIR ${SHARED_DIR}/MUMPS)
    
   find_library(MUMPS_D_LIBRARY libdmumps.a PATHS ${MUMPS_DIR}/lib NO_DEFAULT_PATH)
   find_library(MUMPS_COMMON_LIBRARY libmumps_common.a PATHS ${MUMPS_DIR}/lib NO_DEFAULT_PATH)
   find_library(MUMPS_PORD_LIBRARY libpord.a PATHS ${MUMPS_DIR}/lib NO_DEFAULT_PATH)
    
   find_path(MUMPS_INCLUDE_DIR NAME dmumps_c.h HINTS ${MUMPS_DIR}/include NO_DEFAULT_PATH)
    
   if(NOT MUMPS_D_LIBRARY OR NOT MUMPS_COMMON_LIBRARY OR NOT MUMPS_PORD_LIBRARY OR NOT MUMPS_INCLUDE_DIR)
      message(WARNING "even though WITH_MUMPS_LEAF/WITH_MUMPS_ROOT was specified : no MUMPS was found under ${MUMPS_DIR}/lib")
   endif(NOT MUMPS_D_LIBRARY OR NOT MUMPS_COMMON_LIBRARY OR NOT MUMPS_PORD_LIBRARY OR NOT MUMPS_INCLUDE_DIR)
endif(WITH_MUMPS_LEAF OR WITH_MUMPS_ROOT)

if( MUMPS_INCLUDE_DIR AND MUMPS_D_LIBRARY AND HAVE_PARMETIS )
   # we assume the other libraries are in there
   set(HAVE_MUMPS TRUE)
   message(STATUS "FOUND MUMPS")
   # you may need parmetis if mumps was build with parmetis (we recommend using parmetis). For more info
   # see the readme file in ThirdPartyLibs/MUMPS and the sample build file in PIPS-NLP/Core/LinearSolvers/MumpsSolver
   # In short: insert something like the following line to MUMPS_LIBRARY below
   # -L/home/petra1/work/installs/parmetis-4.0.3/lib -lparmetis
   # hack(ish) way of specifying the many mumps libraries
   # - put all the  mumps libraries: -ldmumps -lmumps_common -lpord
   # - scalapack should be set in ${SCALAPACK_LIBRARIES} using Toolchain.cmake (or hardcoded here)
   # - mpifort also needed
   # - parmetis 4.3.0 or later required by mumps but no need to specify it here: it should be detected by the build system at this point
   # NOTE: METIS 4 is no longer sufficient to link MUMPS against - needs METIS 5...
   if(NOT SCALAPACK_LIBRARIES)
      MESSAGE(STATUS "SCALAPACK_LIBRARIES for MUMPS not set. Default linking for Ubuntu.")
      find_library(SCALAPACK_OMPI_LIBRARY NAMES scalapack-openmpi REQUIRED)
      set(SCALAPACK_LIBRARIES "${SCALAPACK_OMPI_LIBRARY} ")
   endif(NOT SCALAPACK_LIBRARIES)
   message(STATUS "SCALAPACK_LIBRARIES:=${SCALAPACK_LIBRARIES}")
   find_library(BLACS_LIBRARY NAMES blacs-openmpi)
   message(STATUS "BLACS_LIBRARY:=${BLACS_LIBRARY}")
   # -lmpi_mpifh includes the OpenMPI Fortran MPI interface used when
   # building MUMPS. MPICH requires -lmpifort
   if(NOT MUMPS_FORT_LIB)
      MESSAGE(STATUS "MUMPS_FORT_LIB not set for linking the Fortran MPI interface library required by MUMPS. Using OpenMPI Fortran interface -lmpi_mpifh")
      find_library(MUMPS_FORT_LIB mpi_mpifh)
   endif(NOT MUMPS_FORT_LIB)
   set(MUMPS_LIBRARY "-Wl,--start-group ${MUMPS_D_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${MUMPS_FORT_LIB} ${MUMPS_PORD_LIBRARY} ${SCALAPACK_LIBRARIES} ${METIS_5_LIBRARY} \
   ${PARMETIS_LIBRARY} -Wl,--end-group")
   message(STATUS "MUMPS_LIBRARIES:=${MUMPS_LIBRARY}")
   include_directories(${MUMPS_INCLUDE_DIR})
else( MUMPS_INCLUDE_DIR AND MUMPS_D_LIBRARY AND HAVE_PARMETIS )
   set(HAVE_MUMPS FALSE)
   message(STATUS "MUMPS/Parmets library and/or headers are not present. Will build PIPS without MUMPS")
   set(WITH_MUMPS_LEAF FALSE)
   set(WITH_MUMPS_ROOT FALSE)
   set(MUMPS_LIBRARY "")
endif( MUMPS_INCLUDE_DIR AND MUMPS_D_LIBRARY AND HAVE_PARMETIS )

  #MKL_PARDISO Solver
set(HAVE_MKL_PARDISO FALSE)
if(WITH_MKL_PARDISO)
   if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
      if(IS_DIRECTORY $ENV{MKLROOT})
         find_library(MKL_INTEL_LP64_LIBRARY NAME libmkl_intel_lp64.a HINTS $ENV{MKLROOT}/lib/intel64_lin)
         find_library(MKL_CORE_LIBRARY NAME libmkl_core.a HINTS $ENV{MKLROOT}/lib/intel64_lin)
         find_library(MKL_SEQUENTIAL_LIBRARY NAME libmkl_sequential.a HINTS $ENV{MKLROOT}/lib/intel64_lin)
         set(MKL_PARDISO_LIBS "-Wl,--start-group ${MKL_INTEL_LP64_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_SEQUENTIAL_LIBRARY} -Wl,--end-group -lgfortran -lpthread -lm")
         include_directories($ENV{MKLROOT}/include)
         set(HAVE_MKL_PARDISO TRUE)
      endif(IS_DIRECTORY $ENV{MKLROOT})
   endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
message(STATUS "MKL_pardiso libs: \n\t${MKL_INTEL_LP64_LIBRARY} \n\t${MKL_CORE_LIBRARY} \n\t${MKL_SEQUENTIAL_LIBRARY}")
message(STATUS "HAVE_MKL_PARDISO is ${HAVE_MKL_PARDISO}")
endif(WITH_MKL_PARDISO)

### Looking for PARDISO Solver
set(PARDISO_DIR ${SHARED_DIR}/PARDISO/src)
set(HAVE_PARDISO FALSE)
find_library(PARDISO_LIBRARY NAME libpardiso.so HINTS ${PARDISO_DIR})
if(NOT PARDISO_LIBRARY)
   set(PARDISO_LIBRARY "")
   message(STATUS "Could not find Schenk PARDSIO - compiling without it")
else( NOT PARDISO_LIBRARY )
   set(HAVE_PARDISO TRUE)
   message(STATUS "Using PARDISO from ${PARDISO_DIR}/libpardiso.so")
endif( NOT PARDISO_LIBRARY )

#### switch for using MA27 over all other available solvers - since MA27 has no augmented schur complement appoach 
#### implemented we have to set -DPARDISO_BLOCKSC
if( USE_MA27 )
   if( !HAVE_MA27 )
      message(FATAL_ERROR "specified WITH_MA27 but could not find MA27..")
   endif( !HAVE_MA27 )
      
   if( USE_MA57 )
      message(WARNING "specified WITH_MA27 AND WITH_MA57 - only using MA57...")
      set( HAVE_MA27 false )
      set( MA27_LIBRARY  "")
   endif( USE_MA57 )
      
   if( HAVE_MA27 )
      set( HAVE_MA57 false )
      set( HAVE_MUMPS false )
      set( HAVE_PARDISO false )
      set( WITH_MKL_PARDISO false )
      set( HAVE_MKL_PARDISO false )
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPARDISO_BLOCKSC")
   endif( HAVE_MA27 )
endif( USE_MA27 )
   
   
#### switch for using MA57 over all other available solvers - since MA27 has no augmented schur complement appoach 
#### implemented we have to set -DPARDISO_BLOCKSC
if( USE_MA57 )
   if( !HAVE_MA57 )
      message(FATAL_ERROR "specified WITH_MA57 but could not find MA57..")
   endif( !HAVE_MA57 )
   
   set( HAVE_MA27 false )
   set( HAVE_MUMPS false )
   set( HAVE_PARDISO false )
   set( WITH_MKL_PARDISO false )
   set( HAVE_MKL_PARDISO false )
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPARDISO_BLOCKSC")
endif( USE_MA57 )     
   
if(HAVE_MA27)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MA27")
else(HAVE_MA27)
   set(MA27_LIBRARY "")
endif(HAVE_MA27)

if(HAVE_MA57)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MA57")
else(HAVE_MA57)
   set(MA57_LIBRARY "")
endif(HAVE_MA57)

  #if(HAVE_MA86)
  #  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MA86")
  #else(HAVE_MA57)
  #  set(MA86_LIBRARY "")
  #endif(HAVE_MA86)

if(WITH_HIERARCHICAL)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHIERARCHICAL -DPARDISO_BLOCKSC")
endif(WITH_HIERARCHICAL)

if(HAVE_PARDISO)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_PARDISO")
   if(WITH_PARDISO_BLOCKSC)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPARDISO_BLOCKSC")
   endif(WITH_PARDISO_BLOCKSC)    
else(HAVE_PARDISO)
   set(PARDISO_LIBRARY "")
endif(HAVE_PARDISO)

if(HAVE_MUMPS)
   if(WITH_MUMPS_LEAF)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MUMPS_LEAF")
   endif(WITH_MUMPS_LEAF)
   if(WITH_MUMPS_ROOT)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MUMPS_ROOT")
   endif(WITH_MUMPS_ROOT)
endif(HAVE_MUMPS)


if(HAVE_MKL_PARDISO AND WITH_MKL_PARDISO)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MKL_PARDISO")
endif(HAVE_MKL_PARDISO AND WITH_MKL_PARDISO)

#Boost
find_package(Boost COMPONENTS log log_setup regex system filesystem thread chrono date_time)
include_directories(${Boost_INCLUDE_DIR})
#message("---- ${Boost_INCLUDE_DIR}") #message("---- ${Boost_LIBRARIES}")
if (NOT Boost_FOUND)
   message("! A Boost version containing BOOST.LOG has NOT been found.")
   message("! Please specify the location of the Boost library containing 'log' package.")
   message("! For this create a 'Toolchain.cmake' file or edit the one from the root ")
   message("! directory and use 'set(BOOST_ROOT \"/path/to/boost\")'.")
   message("! Then run  'cmake -DCMAKE_TOOLCHAIN_FILE=../Toolchain.cmake ..' in the build directory.")
   message(FATAL_ERROR "Configuration aborted due to missing BOOST libraries!")
endif(NOT Boost_FOUND)

###############################################################################
###############################################################################

include_directories(Input)
include_directories(Input/OPF_Matpower)
include_directories(SolverInterface)

add_subdirectory(Utilities)

add_subdirectory(PIPS-IPM)

##########################################################
# CMake Tests
##########################################################
if(WITH_MAKETEST)
   enable_testing()

   add_test(NAME PIPS-IPM-linkingConsTest COMMAND sh ${PROJECT_SOURCE_DIR}/PIPS-IPM/Test/pipsipmLinkConsTest.sh $<TARGET_FILE:pipsipmCallbackExample> )
endif(WITH_MAKETEST)

message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_STANDARD_LIBRIES: ${CMAKE_CXX_STANDARD_LIBRARIES}")
message(STATUS "CMAKE_LIBRARY_PATH: ${CMAKE_LIBRARY_PATH}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message(STATUS "${inc_dirs}")
get_property(link_dirs DIRECTORY PROPERTY LINK_DIRECTORIES)
message(STATUS "${link_dirs}")

