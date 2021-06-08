#ifndef FILTERLINESEARCH_H
#define FILTERLINESEARCH_H

#include "FilterStrategy.hpp"

class Problem;

class Residuals;

class Variables;

/*! \class LineSearch
 * \brief Line-search
 *
 *  Line-search strategy
 */
class FilterLineSearch {
public:
   /*!
    *  Constructor
    */
   FilterLineSearch(FilterStrategyParameters& filter_strategy_parameters, FilterParameters& filter_parameters, int max_iterations = 30,
         double backtracking_ratio = 0.5, double min_step_length = 1e-9);
   FilterLineSearch(int max_iterations = 30, double backtracking_ratio = 0.5, double min_step_length = 1e-9);

   FilterStrategy filter_strategy;
   /* ratio of step length update in ]0, 1[ */
   double backtracking_ratio;
   int number_iterations;
   bool verbose{false};

   void initialize(Residuals& initial_residuals);
   double predicted_reduction(Problem& problem, Variables& direction, Variables& trial_iterate, double step_length);
   void compute_acceptable_iterate(Problem& problem, Variables& current_iterate, Variables& direction, Residuals& current_residuals,
         double primal_step_length, double dual_step_length);

private:
   double min_step_length;
   int max_iterations;
   bool termination_(bool is_accepted);
};

#endif // FILTERLINESEARCH_H
