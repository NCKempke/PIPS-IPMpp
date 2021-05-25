target_include_directories(pips-ipmpp
        PRIVATE
            Abstract
            Dense
            Distributed
            Sparse
        )

target_sources(pips-ipmpp
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Abstract/Vector.cpp

        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Dense/DenseMatrix.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Dense/DenseStorage.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Dense/DenseSymmetricMatrix.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Dense/SimpleVector.C

        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Distributed/BorderedMatrix.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Distributed/BorderedSymmetricMatrix.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Distributed/DistributedMatrix.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Distributed/DistributedSymmetricMatrix.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Distributed/DistributedVector.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Distributed/StripMatrix.C

        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Sparse/SparseMatrix.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Sparse/SparseStorage.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Sparse/SparseStorageDynamic.C
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Sparse/SparseSymmetricMatrix.C
)