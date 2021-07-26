#add_definitions(-DTIMING -DSTOCH_TESTING) # timing output

add_subdirectory(Core)
add_subdirectory(Drivers)

set(includes_for_interface
		Core/Preprocessing
		Core/Interface
		Core/InteriorPointMethod
		Core/Utilities
        Core/Base
		Core/Options
		)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

############## callbackExample ##############
add_executable(pipsipmCallbackExample Drivers/CallbackExample/callbackExample.cpp)

target_include_directories(pipsipmCallbackExample PRIVATE ${includes_for_interface})
target_compile_options(pipsipmCallbackExample PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(pipsipmCallbackExample
        PRIVATE pips-ipmpp
        PRIVATE OpenMP::OpenMP_CXX
        PRIVATE MPI::MPI_CXX
        )

############## gmspips ##############
add_executable(gmspips Drivers/gams/gmspips/gmspips.cpp)

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
