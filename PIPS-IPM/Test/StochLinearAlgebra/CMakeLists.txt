include_directories(../../Core/Options)
include_directories(../../Core/Interface)
include_directories(../../Core/LinearAlgebra/Distributed)
include_directories(../../Core/Vector)
include_directories(../../Core/Base)
include_directories(../../Core/Utilities)
include_directories(../../Core/LinearAlgebra/Sparse)

package_add_test(DistributedMatrixTest t_DistributedMatrix.cpp)
