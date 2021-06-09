//
// Created by bzfkempk on 27.04.21.
//

#include "RegularizationStrategy.h"
#include "PIPSIPMppOptions.h"
#include "pipsdef.h"
#include<cmath>

RegularizationStrategy::RegularizationStrategy(unsigned int positive_eigenvalues_expected_, unsigned int negative_eigenvalues_expected_)
      : positive_eigenvalues_expected{positive_eigenvalues_expected_}, negative_eigenvalues_expected{negative_eigenvalues_expected_}
{}

bool RegularizationStrategy::is_inertia_correct(const Inertia& inertia) const {
   auto[positive_eigenvalues, negative_eigenvalues, zero_eigenvalues] = inertia;

   std::cout << "comparing inertia of (" << positive_eigenvalues << "," << negative_eigenvalues << "," << zero_eigenvalues << ") against expected: (" <<
      positive_eigenvalues_expected << "," << negative_eigenvalues_expected << "," << "0)\n";
   return positive_eigenvalues_expected == positive_eigenvalues && negative_eigenvalues_expected == negative_eigenvalues && zero_eigenvalues == 0;
}

void RegularizationStrategy::notify_new_step() {
   new_factorization = true;
};