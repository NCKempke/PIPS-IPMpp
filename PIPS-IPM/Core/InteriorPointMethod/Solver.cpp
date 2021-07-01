#include "Solver.hpp"
#include "Problem.hpp"
#include "Variables.h"
#include "Residuals.h"
#include "AbstractLinearSystem.h"
#include "AbstractOptions.h"
#include "DistributedFactory.hpp"
#include <cmath>

double g_iterNumber = 0.;
int print_level = 1000;

int gOuterBiCGIter = 0;
double gOuterBiCGIterAvg = 0.;

Solver::Solver(DistributedFactory& factory, Problem& problem) : factory(factory), step(factory.make_variables(problem)) {
}

void Solver::solve_linear_system(Variables& iterate, Problem& problem, Residuals& residuals, Variables& step, AbstractLinearSystem& linear_system) {
   residuals.evaluate(problem, iterate);
   residuals.set_complementarity_residual(iterate, 0.);

   linear_system.factorize(iterate);
   linear_system.solve(iterate, residuals, step);
   step.negate();

   // take the full affine scaling step
   iterate.add(step, 1.);
   double shift = 1e3 + 2 * iterate.violation();
   iterate.shift_bound_variables(shift, shift);
}