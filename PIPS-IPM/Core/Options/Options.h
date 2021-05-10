/*
 * QpGenOptions.h
 *
 *  Created on: 01.07.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPGEN_QPGENOPTIONS_H_
#define PIPS_IPM_CORE_QPGEN_QPGENOPTIONS_H_

#include "Options.h"
#include "pipsport.h"

#include <cassert>

/**
 * The getInstanceMethod must be specified in this BaseClass.
 * Defines default options.
 */

namespace qpgen_options {
   void setOptions(const std::string& opt_file);
   int getIntParameter(const std::string& identifier);
   double getDoubleParameter(const std::string& identifier);
   bool getBoolParameter(const std::string& identifier);

   class QpGenOptions : public base_options::Options {

   protected :
      friend void setOptions(const std::string& opt_file);
      friend int getIntParameter(const std::string& identifier);
      friend double getDoubleParameter(const std::string& identifier);
      friend bool getBoolParameter(const std::string& identifier);

      static QpGenOptions& getInstance() {
         static QpGenOptions opt;
         return opt;
      }

      void setDefaults() override;
      QpGenOptions();

      virtual ~QpGenOptions() {};
   };

   inline void setOptions(const std::string& opt_file) {
      return QpGenOptions::getInstance().fillOptionsFromFile(opt_file);
   }

   inline int getIntParameter(const std::string& identifier) {
      return QpGenOptions::getInstance().getIntParam(identifier);
   }

   inline bool getBoolParameter(const std::string& identifier) {
      return QpGenOptions::getInstance().getBoolParam(identifier);
   }

   inline double getDoubleParameter(const std::string& identifier) {
      return QpGenOptions::getInstance().getDoubleParam(identifier);
   }
}

#endif /* PIPS_IPM_CORE_QPGEN_QPGENOPTIONS_H_ */