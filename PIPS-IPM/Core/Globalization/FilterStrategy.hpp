#ifndef FILTERSTRATEGY_H
#define FILTERSTRATEGY_H

#include <iostream>
#include <memory>
#include "Filter.hpp"

class FilterParameters;

class Residuals;

class Variables;

/*! \class FilterStrategyParameters
 * \brief Constants for filter and tube strategies
 *
 *  Set of constants to control the filter and tube strategies
 */
struct FilterStrategyParameters {
   double Sigma; /*!< Sufficient reduction constant */
   double Delta; /*!< Switching constant */
   double ubd;
   double fact;
};

/*! \class FilterStrategy
 * \brief Step acceptance strategy based on a filter
 *
 *  Strategy that accepts or declines a trial step
 */
class FilterStrategy {
public:
   FilterStrategy(FilterStrategyParameters& filter_strategy_parameters, FilterParameters& filter_parameters);
   FilterStrategy();

   Filter filter;
   FilterStrategyParameters parameters; /*!< Set of parameters */
   double verbose{false};

   void initialize(Residuals& initial_residuals);
   bool check_acceptance(Variables& current_iterate, Residuals& current_residuals, Variables& trial_iterate, Residuals& trial_residuals,
         double predicted_reduction, double step_length);
};

#endif // FILTERSTRATEGY_H
