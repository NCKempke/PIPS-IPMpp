############## callbackExample ##############
add_executable(pipsipmCallbackExample callbackExample.cpp)

target_include_directories(pipsipmCallbackExample PRIVATE ${includes_for_interface})
target_compile_options(pipsipmCallbackExample PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(pipsipmCallbackExample
        PRIVATE pips-ipmpp
        PRIVATE OpenMP::OpenMP_CXX
        PRIVATE MPI::MPI_CXX
        )