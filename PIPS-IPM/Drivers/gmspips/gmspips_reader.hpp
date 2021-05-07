/*
 * gmspips_reader.h
 *
 *  Created on: 19.03.2021
 *      Author: bzfkempk
 */


#ifndef PIPS_IPM_DRIVERS_GMSPIPS_GMSPIPSREADER_HPP_
#define PIPS_IPM_DRIVERS_GMSPIPS_GMSPIPSREADER_HPP_

#include <string>

#include "pipsdef.h"
#include "gmspipsio.h"
#include "DistributedInputTree.h"

class gmspips_reader {
public:
   gmspips_reader(const std::string& path_to_problem, const std::string& path_to_gams, size_t n_blocks);
   virtual ~gmspips_reader();

   DistributedInputTree* read_problem();

protected:
   gmspips_reader() {};


   std::vector<GMSPIPSBlockData_t*> blocks;

   size_t n_blocks{};
   bool log_reading{false};
};

#endif /* PIPS_IPM_DRIVERS_GMSPIPS_GMSPIPSREADER_HPP_ */
