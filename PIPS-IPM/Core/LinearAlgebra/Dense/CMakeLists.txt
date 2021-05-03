add_library(ooqpdense STATIC DenseStorage.C DenseSymMatrix.C
        DeSymIndefSolver.C DeSymIndefSolver2.C
        DenseGenMatrix.C DeSymPSDSolver.C
    )
target_link_libraries(ooqpdense PRIVATE ${MATH_LIBS})
