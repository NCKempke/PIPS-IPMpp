# add_library(ooqpmehrotra STATIC MehrotraSolver.C)
# add_library(ooqpmehrotrastoch STATIC MehrotraStochSolver.C MehrotraSolver.C)
# add_library(ooqpgondzio STATIC GondzioSolver.C)
add_library(ooqpgondziostoch STATIC Solver.C GondzioStochSolver.C GondzioStochLpSolver.C InteriorPointMethod.cpp)



