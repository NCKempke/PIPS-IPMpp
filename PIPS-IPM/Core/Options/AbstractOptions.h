/*
 * Options.h
 *
 *  Created on: 03.04.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_ABSTRACT_OPTIONS_H_
#define PIPS_IPM_CORE_ABSTRACT_OPTIONS_H_

#include <string>
#include <map>
#include "Singleton.h"

/**
 * base class for options class.
 */

namespace abstract_options {
   int get_int_parameter(const std::string& identifier);
   double get_double_parameter(const std::string& identifier);
   bool get_bool_parameter(const std::string& identifier);

   class AbstractOptions : public Singleton {
   private:
      virtual void setDefaults() {};

   protected:
      // not thread safe when modified..
      // TODO : there is no hash_map in C++03 I think
      static std::map<std::string, double> double_options;
      static std::map<std::string, int> int_options;
      static std::map<std::string, bool> bool_options;

      friend int get_int_parameter(const std::string& identifier);
      friend double get_double_parameter(const std::string& identifier);
      friend bool get_bool_parameter(const std::string& identifier);

      AbstractOptions();
      ~AbstractOptions() override = default;

      static AbstractOptions& get_instance() {
         static AbstractOptions opt;
         return opt;
      }

      static bool is_identifier_unique(const std::string& identifier);
      static bool identifier_exists(const std::string& identifier);
      static void load_options_from_file(const std::string& filename);

      static int get_int_param(const std::string& identifier) ;
      static double get_double_param(const std::string& identifier) ;
      static bool get_bool_param(const std::string& identifier) ;

      static void set_int_param(const std::string& param, int value);
      static void set_bool_param(const std::string& param, bool value);
      static void set_double_param(const std::string& param, double value);
   };

   inline int get_int_parameter(const std::string& identifier) {
      return AbstractOptions::get_instance().get_int_param(identifier);
   }

   inline bool get_bool_parameter(const std::string& identifier) {
      return AbstractOptions::get_instance().get_bool_param(identifier);
   }

   inline double get_double_parameter(const std::string& identifier) {
      return AbstractOptions::get_instance().get_double_param(identifier);
   }
}

#endif /* PIPS_IPM_CORE_ABSTRACT_OPTIONS_H_ */
