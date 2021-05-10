/*
 * QpGenOptions.h
 *
 *  Created on: 01.07.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPGEN_QPGENOPTIONS_H_
#define PIPS_IPM_CORE_QPGEN_QPGENOPTIONS_H_

#include "AbstractOptions.h"
#include "pipsport.h"

#include <cassert>

/**
 * The getInstanceMethod must be specified in this BaseClass.
 * Defines default options.
 */

namespace options {
   void setOptions(const std::string& opt_file);
   int getIntParameter(const std::string& identifier);
   double getDoubleParameter(const std::string& identifier);
   bool getBoolParameter(const std::string& identifier);

   class Options : public abstract_options::AbstractOptions {

   protected :
      friend void setOptions(const std::string& opt_file);
      friend int getIntParameter(const std::string& identifier);
      friend double getDoubleParameter(const std::string& identifier);
      friend bool getBoolParameter(const std::string& identifier);

      static Options& getInstance() {
         static Options opt;
         return opt;
      }

      void setDefaults() override;
      Options();

      virtual ~Options() {};
   };

   inline void setOptions(const std::string& opt_file) {
      return Options::getInstance().load_options_from_file(opt_file);
   }

   inline int getIntParameter(const std::string& identifier) {
      return Options::getInstance().get_int_param(identifier);
   }

   inline bool getBoolParameter(const std::string& identifier) {
      return Options::getInstance().get_bool_param(identifier);
   }

   inline double getDoubleParameter(const std::string& identifier) {
      return Options::getInstance().get_double_param(identifier);
   }
}

#endif /* PIPS_IPM_CORE_QPGEN_QPGENOPTIONS_H_ */
