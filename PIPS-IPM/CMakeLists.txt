#add_definitions(-DTIMING -DSTOCH_TESTING) # timing output
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

add_subdirectory(Core)
add_subdirectory(Drivers)