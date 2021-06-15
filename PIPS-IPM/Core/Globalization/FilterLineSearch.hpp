#ifndef FILTERLINESEARCH_H
#define FILTERLINESEARCH_H

#include "FilterStrategy.hpp"
#include "Scaler.hpp"

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
   explicit FilterLineSearch(const Scaler* scaler = nullptr, int max_iterations = 30, double backtracking_ratio = 0.5, double min_step_length = 1e-9);

   FilterStrategy filter_strategy;
   /* ratio of step length update in ]0, 1[ */
   double backtracking_ratio;
   int number_iterations;
   bool verbose{false};

   void initialize(Residuals& initial_residuals);
   void compute_acceptable_iterate(Problem& problem, Variables& current_iterate, Variables& direction, Residuals& current_residuals,
         double& primal_step_length, double& dual_step_length);

private:
   const Scaler* scaler{};
   double min_step_length;
   int max_iterations;
   bool termination_(bool is_accepted);
};

#endif // FILTERLINESEARCH_H
