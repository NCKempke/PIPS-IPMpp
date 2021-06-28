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
   void initialize(Residuals& initial_residuals);

   bool check_acceptance(Residuals& current_residuals, Residuals& trial_residuals, double predicted_reduction);

private:
   bool verbose{false};
};

#endif // FILTERSTRATEGY_H
