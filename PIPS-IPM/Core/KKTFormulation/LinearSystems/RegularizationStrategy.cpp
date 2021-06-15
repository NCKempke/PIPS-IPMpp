//
// Created by bzfkempk on 27.04.21.
//

#include "RegularizationStrategy.h"
#include "PIPSIPMppOptions.h"
#include "pipsdef.h"

RegularizationStrategy::RegularizationStrategy(unsigned int positive_eigenvalues_expected_, unsigned int negative_eigenvalues_expected_, MPI_Comm mpi_comm_)
      : mpi_comm{mpi_comm_}, positive_eigenvalues_expected{positive_eigenvalues_expected_}, negative_eigenvalues_expected{negative_eigenvalues_expected_}
{}

bool RegularizationStrategy::is_inertia_correct(const Inertia& inertia) const {
   auto[positive_eigenvalues, negative_eigenvalues, zero_eigenvalues] = inertia;

   if (pipsipmpp_options::get_bool_parameter("REGULARIZATION_VERBOSE") && PIPS_MPIgetRank(mpi_comm) == 0) {
      std::cout << "comparing inertia of (" << positive_eigenvalues << "," << negative_eigenvalues << "," << zero_eigenvalues << ") against expected: (" <<
                positive_eigenvalues_expected << "," << negative_eigenvalues_expected << "," << "0)\n";
   }
   return positive_eigenvalues_expected == positive_eigenvalues && negative_eigenvalues_expected == negative_eigenvalues && zero_eigenvalues == 0;
}

void RegularizationStrategy::notify_new_step() {
   new_factorization = true;
};