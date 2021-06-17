/*
 * gmspips_reader.h
 *
 *  Created on: 19.03.2021
 *      Author: bzfkempk
 */


#ifndef PIPS_IPM_DRIVERS_GMSPIPS_GMSPIPSREADER_HPP_
#define PIPS_IPM_DRIVERS_GMSPIPS_GMSPIPSREADER_HPP_

#include <string>
#include <cassert>
#include "pipsdef.h"

#include "../pips_reader.h"
#include "gmspipsio.h"
#include "DistributedInputTree.h"

class gmspips_reader : pips_reader {
public:
   gmspips_reader(const std::string& path_to_problem, const std::string& path_to_gams, size_t n_blocks);
   ~gmspips_reader() override;

   std::unique_ptr<DistributedInputTree> read_problem() override;

protected:
   std::vector<GMSPIPSBlockData_t*> blocks;
   bool log_reading{false};
};

#endif /* PIPS_IPM_DRIVERS_GMSPIPS_GMSPIPSREADER_HPP_ */
