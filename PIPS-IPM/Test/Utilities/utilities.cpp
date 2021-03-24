/*
 * utilities.cpp
 *
 *  Created on: 19.03.2021
 *      Author: bzfkempk
 */
#include "utilities.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <tuple>

bool isComment(std::string& line)
{
   if( line.size() > 1 && line[0] == '#' )
      return true;
   else if( line.size() > 2 && line[0] == '/' && line[1] == '/')
      return true;
   else
      return false;
}

/* expects the file to be like [name] [n_blocks] [expected_objective] - "#" and "//" is ignored as comments */
std::vector<Instance> readInstancesFromFile( const std::string& file_name )
{
   std::vector<Instance> instances;

   std::ifstream params;
   params.open(file_name.c_str(), std::ios::in);

   if( !params.good() )
   {
      std::cout << "Failed to open provided instance file \"" << file_name << "\"\n";
      return instances;
   }

   std::string line;
   while( std::getline(params, line) )
   {
      if( line.compare("") == 0 )
         continue;

      if( isComment(line) )
         continue;

      std::istringstream iss(line);

      std::string instance_path;
      size_t n_blocks;
      double objective;
      if( !(iss >> instance_path >> n_blocks >> objective) )
      {
         /* some error while reading that line occured */
         std::cout << "Error while reading line \"" << line << "\" from options file - skipping that line\n";
         continue;
      }

      instances.push_back( { instance_path, n_blocks, objective } );
   }
   return instances;
}
