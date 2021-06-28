#ifndef FILTER_H
#define FILTER_H

#include <ostream>
#include <vector>
#include <list>
#include <map>
#include <memory>

struct FilterParameters {
   double Beta; /*!< Margin around filter */
   double Gamma; /*!< Margin around filter (sloping margin) */
};

struct FilterEntry {
   double infeasibility_measure;
   double optimality_measure;
};

/*! \class Filter
 * \brief Filter
 *
 *  Filter
 */
class Filter {
public:
   Filter(FilterParameters& constants);
   Filter();
   ~Filter() = default;

   double upper_bound; /*!< Upper bound on constraint violation */
   unsigned int max_size; /*!< Max filter size */
   FilterParameters constants; /*!< Set of constants */

   void reset();
   void add(double infeasibility_measure, double optimality_measure);
   bool accept(double infeasibility_measure, double optimality_measure);
   bool improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
         double trial_optimality_measure) const;

   friend std::ostream& operator<<(std::ostream& stream, Filter& filter);

protected:
   std::list<FilterEntry> entries_;
};

#endif // FILTER_H
