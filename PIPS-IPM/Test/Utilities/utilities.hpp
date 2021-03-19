/*
 * utilities.hpp
 *
 *  Created on: 19.03.2021
 *      Author: bzfkempk
 */

#include <iostream>
#include <fstream>
#include <memory>

#include <string>
#include <vector>
#include <tuple>


struct Instance{
      const std::string name;
      const size_t n_blocks;
      const double result;
};

std::vector<Instance> readInstancesFromFile( const std::string& file_name );
