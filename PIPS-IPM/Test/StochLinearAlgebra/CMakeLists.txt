include_directories(../../Core/QpGen)
include_directories(../../Core/QpStoch)
include_directories(../../Core/LinearAlgebra/Distributed)
include_directories(../../Core/Vector)
include_directories(../../Core/Abstract)
include_directories(../../Core/Utilities)
include_directories(../../Core/LinearAlgebra/Sparse)

package_add_test(DistributedMatrixTest t_DistributedMatrix.cpp)
