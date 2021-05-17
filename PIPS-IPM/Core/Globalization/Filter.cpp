#include <iostream>
#include <cmath>
#include "Filter.hpp"

Filter::Filter(FilterParameters& constants) : constants(constants) {
   this->reset();
}

Filter::Filter() : constants({0.999, 0.001}) {
   this->reset();
}

void Filter::reset() {
   /* initialize the maximum filter size (not critical) */
   this->max_size = 50;
   this->upper_bound = INFINITY;
   this->entries_.clear();
   return;
}

/*  add (infeasibility_measure, optimality_measure) to the filter */
void Filter::add(double infeasibility_measure, double optimality_measure) {
   /* remove dominated filter entries */
   std::list<FilterEntry>::iterator entry = this->entries_.begin();
   while (entry != this->entries_.end()) {
      if (infeasibility_measure < entry->infeasibility_measure && optimality_measure <= entry->optimality_measure) {
         entry = this->entries_.erase(entry);
      }
      else {
         entry++;
      }
   }

   /* check sufficient space available for new entry (remove last entry, if not) */
   if (this->max_size <= this->entries_.size()) {
      FilterEntry& last_element = this->entries_.back();
      this->upper_bound = this->constants.Beta * std::max(this->upper_bound, last_element.infeasibility_measure);
      this->entries_.pop_back();
   }

   /* add new entry to the filter */
   std::list<FilterEntry>::iterator position = this->entries_.begin();
   while (position != this->entries_.end() && infeasibility_measure >= this->constants.Beta * position->infeasibility_measure) {
      position++;
   }
   FilterEntry new_entry{infeasibility_measure, optimality_measure};
   this->entries_.insert(position, new_entry);

   return;
}

// filter must be nonempty

double Filter::eta_min() {
   FilterEntry& last_element = this->entries_.back();
   return last_element.infeasibility_measure;
}

// filter must be nonempty

double Filter::omega_min() {
   FilterEntry& last_element = this->entries_.back();
   return last_element.optimality_measure;
}

/* query: return true if (infeasibility_measure, optimality_measure) acceptable, false otherwise */
bool Filter::accept(double infeasibility_measure, double optimality_measure) {
   /* check upper bound first */
   if (this->constants.Beta * this->upper_bound <= infeasibility_measure) {
      return false;
   }

   std::list<FilterEntry>::iterator position = this->entries_.begin();
   while (position != this->entries_.end() && infeasibility_measure >= this->constants.Beta * position->infeasibility_measure) {
      position++;
   }

   /* check acceptability */
   if (position == this->entries_.begin()) {
      return true; // acceptable as left-most entry
   }
   else if (optimality_measure <= std::prev(position)->optimality_measure - this->constants.Gamma * infeasibility_measure) {
      return true; // point acceptable
   }
   else {
      return false; // point rejected
   }
}

//! improves_current_iterate: check acceptable wrt current point 

bool Filter::improves_current_iterate(double current_infeasibility_measure, double current_optimality_measure, double trial_infeasibility_measure,
      double trial_optimality_measure) {
   return (trial_optimality_measure <= current_optimality_measure - this->constants.Gamma * trial_infeasibility_measure) ||
          (trial_infeasibility_measure < this->constants.Beta * current_infeasibility_measure);
}

double Filter::compute_actual_reduction(double current_objective, double /*current_residual*/, double trial_objective) {
   return current_objective - trial_objective;
}

//! print: print the content of the filter

std::ostream& operator<<(std::ostream& stream, Filter& filter) {
   stream << "************\n";
   stream << "  Current filter (constraint residual, objective):\n";
   for (FilterEntry const& entry: filter.entries_) {
      stream << "\t" << entry.infeasibility_measure << "\t" << entry.optimality_measure << "\n";
   }
   stream << "************\n";
   return stream;
}