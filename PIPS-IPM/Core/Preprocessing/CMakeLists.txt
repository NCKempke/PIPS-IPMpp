target_include_directories(pips-ipmpp PRIVATE ./)

target_sources(pips-ipmpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/CurtisReidScaler.cpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/EquilibriumScaler.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/GeoStochScaler.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/PresolveData.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Postsolver.cpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Presolver.cpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Scaler.cpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StochColumnStorage.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StochPostsolver.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StochPresolver.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StochPresolverBase.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StochPresolverBoundStrengthening.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StochPresolverColumnFixation.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StochPresolverModelCleanup.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StochPresolverParallelRows.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StochPresolverSingletonColumns.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StochPresolverSingletonRows.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/StochRowStorage.C
        )