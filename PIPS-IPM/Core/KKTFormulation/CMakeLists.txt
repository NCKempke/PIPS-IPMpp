target_include_directories(pips-ipmpp
    PRIVATE
        LinearSystems
        Residuals
        Variables
    )

target_sources(pips-ipmpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/BorderMod_Block.cpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/DistributedLeafLinearSystem.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/DistributedLinearSystem.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/DistributedRootLinearSystem.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/FriedlanderOrbanRegularization.cpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/IpoptRegularization.cpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/LinearSystem.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/RACFG_BLOCK.cpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/RegularizationStrategy.cpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/sLinsysLeafSchurSlv.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/sLinsysRootAug.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/sLinsysRootAugHierInner.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinearSystems/sLinsysRootBordered.C

        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Residuals/DistributedResiduals.cpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Residuals/Residuals.cpp

        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Variables/DistributedVariables.cpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Variables/Variables.C
    )