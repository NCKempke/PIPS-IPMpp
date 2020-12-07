add_library(ooqpsparse STATIC SparseStorage.C SparseStorageDynamic.C SparseLinearAlgebraPackage.C 
  SparseGenMatrix.C SparseSymMatrix.C)
target_link_libraries(ooqpsparse PRIVATE ${MATH_LIBS})