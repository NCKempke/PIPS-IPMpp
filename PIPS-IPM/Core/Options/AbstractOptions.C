#include "AbstractOptions.h"
#include "pipsdef.h"

#include <fstream>
#include <sstream>

#ifdef PRE_CPP11
#include <cstdlib>
#endif

namespace abstract_options {
   std::map<std::string, double> AbstractOptions::double_options;
   std::map<std::string, int> AbstractOptions::int_options;
   std::map<std::string, bool> AbstractOptions::bool_options;

   AbstractOptions::AbstractOptions() {
      /// INTERIOR-POINT ALGORITHM
      bool_options["IP_ACCURACY_REDUCED"] = false;
      bool_options["IP_PRINT_TIMESTAMP"] = false;
      bool_options["IP_STEPLENGTH_CONSERVATIVE"] = false;
   }

   int AbstractOptions::get_int_param(const std::string& identifier) const {
      const std::map<std::string, int>::const_iterator& it = int_options.find(identifier);

      if (it != int_options.end())
         return it->second;
      else {
         std::cout << "No element \"" << identifier << "\" of type int in options - using 0 instead" << std::endl;
         return 0;
      }
   }

   double AbstractOptions::get_double_param(const std::string& identifier) const {
      const std::map<std::string, double>::const_iterator it = double_options.find(identifier);

      if (it != double_options.end())
         return it->second;
      else {
         std::cout << "No element \"" << identifier << "\" of type double in options - using 0 instead" << std::endl;
         return 0;
      }
   }

   bool AbstractOptions::get_bool_param(const std::string& identifier) const {
      const std::map<std::string, bool>::const_iterator it = bool_options.find(identifier);

      if (it != bool_options.end())
         return it->second;
      else {
         std::cout << "No element \"" << identifier << "\" of type bool in options - using false instead" << std::endl;
         return false;
      }
   }

   void AbstractOptions::set_int_param(const std::string& param, int value) {
      int_options[param] = value;
   }

   void AbstractOptions::set_bool_param(const std::string& param, bool value) {
      bool_options[param] = value;
   }

   void AbstractOptions::set_double_param(const std::string& param, double value) {
      double_options[param] = value;
   }

   bool isComment(std::string& line) {
      if (line.size() > 1 && line[0] == '#')
         return true;
      else if (line.size() > 2 && line[0] == '/' && line[1] == '/')
         return true;
      else
         return false;
   }

   void AbstractOptions::load_options_from_file(const std::string& filename) {
      std::ifstream params;
      params.open(filename.c_str(), std::ios::in);
      const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);

      if (!params.good()) {
         if (my_rank == 0) {
            std::cout << "Failed to open provided options file \"" << filename << "\"" << std::endl;
            std::cout << "Using default configuration" << std::endl;
         }
         return;
      }

      std::string line;
      while (std::getline(params, line)) {
         if (line.compare("") == 0)
            continue;

         if (isComment(line))
            continue;

         std::istringstream iss(line);

         std::string identifier;
         std::string value;
         std::string type;
         if (!(iss >> identifier >> value >> type)) {
            /* some error while reading that line occured */
            if (my_rank == 0)
               std::cout << "Error while reading line \"" << line << "\" from options file - skipping that line" << std::endl;
            continue;
         }

         if (identifier.length() == 0 || identifier[0] == '#')
            continue;

         try {
            if (!identifier_exists(identifier)) {
               if (my_rank == 0)
                  std::cout << "Warning - unknown identifier - skipping it: " << identifier << std::endl;
            }
            else if (type.compare("int") == 0 || type.compare("integer") == 0) {
#ifdef PRE_CPP11
               int_options[identifier] = atoi(value.c_str());
#else
               int_options[identifier] = std::stoi(value);
#endif
            }
            else if (type.compare("double") == 0) {
#ifdef PRE_CPP11
               double_options[identifier] = atof(value.c_str());
#else
               double_options[identifier] = std::stod(value);
#endif
            }
            else if (type.compare("bool") == 0 || type.compare("boolean") == 0) {
               if (value.compare("true") == 0 || value.compare("TRUE") == 0 || value.compare("True") == 0)
                  bool_options[identifier] = true;
               else if (value.compare("false") == 0 || value.compare("FALSE") == 0 || value.compare("False") == 0)
                  bool_options[identifier] = false;
               else if (my_rank == 0)
                  std::cout << "Unknown value \"" << value << "\" for bool parameter \"" << identifier << "\" in options file - skipping that line"
                            << std::endl;
            }
            else {
               if (my_rank == 0)
                  std::cout << "Unknown type \"" << type << "\" in options file (supported: bool/boolean, double, int/integer)- skipping that line"
                            << std::endl;
            }
         }
         catch (std::exception& e) {
            if (my_rank == 0)
               std::cout << "Error while reading line \"" << line << "\" from options file - could not parse value of " << identifier << " " << value
                         << " - skipping that line" << std::endl;
         }
      }
   }

   bool AbstractOptions::identifier_exists(const std::string& identifier) const {
      return int_options.find(identifier) != int_options.end() || bool_options.find(identifier) != bool_options.end() ||
             double_options.find(identifier) != double_options.end();
   }

   bool AbstractOptions::is_identifier_unique(const std::string& identifier) const {
      int n_found = 0;

      if (int_options.find(identifier) != int_options.end())
         n_found++;
      if (bool_options.find(identifier) != bool_options.end())
         n_found++;
      if (double_options.find(identifier) != double_options.end())
         n_found++;

      return n_found <= 1;
   }
}
