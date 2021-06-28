#ifndef FILTERLINESEARCH_H
#define FILTERLINESEARCH_H

#include "FilterStrategy.hpp"
#include "Scaler.hpp"
#include "InteriorPointMethod.hpp"

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
   explicit FilterLineSearch(DistributedFactory& factory, Problem& problem, double dnorm, InteriorPointMethodType interior_point_method_type,
         const Scaler* scaler = nullptr);

   FilterStrategy filter_strategy;
   /* ratio of step length update in ]0, 1[ */
   const double backtracking_ratio{0.5};
   int number_iterations{0};
   void initialize(Residuals& initial_residuals);
   void register_observer(AbstractLinearSystem* linear_system);
   void compute_acceptable_iterate(Problem& problem, Variables& iterate, Residuals& residuals, Variables& step, AbstractLinearSystem& linear_system,
         int iteration);
   void compute_acceptable_iterate(Problem& problem, Variables& current_iterate, Variables& direction, Residuals& current_residuals);

private:

   std::unique_ptr<InteriorPointMethod> interior_point_method;
   const Scaler* scaler{};
   double min_step_length{1e-9};
   const int max_iterations{30};
   const bool verbose{false};

   [[nodiscard]] bool termination_(bool is_accepted) const;
};

#endif // FILTERLINESEARCH_H
