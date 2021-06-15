#include "AbstractOptions.h"
#include "pipsdef.h"

#include <fstream>
#include <sstream>

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

   int AbstractOptions::get_int_param(const std::string& identifier) {
      const std::map<std::string, int>::const_iterator& it = int_options.find(identifier);

      if (it != int_options.end())
         return it->second;
      else {
         std::cout << "No element \"" << identifier << "\" of type int in options - using 0 instead\n";
         return 0;
      }
   }

   double AbstractOptions::get_double_param(const std::string& identifier) {
      const auto it = double_options.find(identifier);

      if (it != double_options.end())
         return it->second;
      else {
         std::cout << "No element \"" << identifier << "\" of type double in options - using 0 instead\n";
         return 0;
      }
   }

   bool AbstractOptions::get_bool_param(const std::string& identifier) {
      const auto it = bool_options.find(identifier);

      if (it != bool_options.end())
         return it->second;
      else {
         std::cout << "No element \"" << identifier << "\" of type bool in options - using false instead\n";
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
            std::cout << "Failed to open provided options file \"" << filename << "\"\n";
            std::cout << "Using default configuration\n";
         }
         return;
      }

      std::string line;
      while (std::getline(params, line)) {
         if (line.empty())
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
               std::cout << "Error while reading line \"" << line << "\" from options file - skipping that line\n";
            continue;
         }

         if (identifier.length() == 0 || identifier[0] == '#')
            continue;

         try {
            if (!identifier_exists(identifier)) {
               if (my_rank == 0)
                  std::cout << "Warning - unknown identifier - skipping it: " << identifier << "\n";
            }
            else if (type == "int" || type == "integer") {
               int_options[identifier] = std::stoi(value);
            }
            else if (type == "double") {
               double_options[identifier] = std::stod(value);
            }
            else if (type == "bool" || type == "boolean") {
               if (value == "true" || value == "TRUE" || value == "True")
                  bool_options[identifier] = true;
               else if (value == "false" || value == "FALSE" || value == "False")
                  bool_options[identifier] = false;
               else if (my_rank == 0)
                  std::cout << "Unknown value \"" << value << "\" for bool parameter \"" << identifier << "\" in options file - skipping that line\n";
            }
            else {
               if (my_rank == 0)
                  std::cout << "Unknown type \"" << type << "\" in options file (supported: bool/boolean, double, int/integer)- skipping that line\n";
            }
         }
         catch (std::exception& e) {
            if (my_rank == 0)
               std::cout << "Error while reading line \"" << line << "\" from options file - could not parse value of " << identifier << " " << value
                         << " - skipping that line\n";
         }
      }
   }

   bool AbstractOptions::identifier_exists(const std::string& identifier) {
      return int_options.find(identifier) != int_options.end() || bool_options.find(identifier) != bool_options.end() ||
             double_options.find(identifier) != double_options.end();
   }

   bool AbstractOptions::is_identifier_unique(const std::string& identifier) {
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
