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
   void compute_acceptable_iterate(Problem& problem, Variables& current_iterate, Residuals& current_residuals, Variables& step, AbstractLinearSystem& linear_system,
         int iteration);
   void print_statistics(const Problem& problem, const Variables& iterate, const Residuals& residuals, int i, double mu,
         TerminationStatus stop_code, int level);

private:

   std::unique_ptr<InteriorPointMethod> interior_point_method;
   const Scaler* scaler{};
   /** norm of problem data */
   double dnorm{0.};
   double min_step_length{1e-9};
   const int max_iterations{20};
   const bool verbose{false};

   [[nodiscard]] bool termination_(bool is_accepted) const;
};

#endif // FILTERLINESEARCH_H
