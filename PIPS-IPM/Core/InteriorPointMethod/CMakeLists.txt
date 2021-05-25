target_include_directories(pips-ipmpp PRIVATE ./)

target_sources(pips-ipmpp
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/InteriorPointMethod.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MehrotraStrategy.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Solver.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Statistics.cpp
)