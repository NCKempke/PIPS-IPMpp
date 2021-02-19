/*
 * StochOptions.h
 *
 *  Created on: 03.04.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_
#define PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_

#include "QpGenOptions.h"
#include "pipsport.h"

#include <cassert>
#include <vector>
#include <ostream>
/**
 * The getInstanceMethod must be specified in this BaseClass.
 * Defines default options for StochPIPS.
 */


enum SolverType
{
   SOLVER_NONE = 0, SOLVER_MA27 = 1, SOLVER_MA57 = 2, SOLVER_PARDISO = 3, SOLVER_MKL_PARDISO = 4, SOLVER_MUMPS = 5
};

enum SolverTypeDense
{
   SOLVER_DENSE_SYM_INDEF = 0, SOLVER_DENSE_SYM_INDEF_SADDLE_POINT = 1, SOLVER_DENSE_SYM_PSD = 2
};

std::ostream& operator<<(std::ostream& os, const SolverType solver);
std::ostream& operator<<(std::ostream& os, const SolverTypeDense solver);


namespace pips_options
{
   bool isSolverAvailable(SolverType solver);
   SolverType getSolverRoot();
   SolverType getSolverSubRoot();
   SolverType getSolverLeaf();
   SolverTypeDense getSolverDense();

   void activateHierarchialApproach();

   void setOptions(const std::string& opt_file);
   void setIntParameter(const std::string& identifier, int value);
   void setDoubleParameter(const std::string& identifier, double value);
   void setBoolParameter(const std::string& identifier, bool value);

   int getIntParameter(const std::string& identifier);
   double getDoubleParameter(const std::string& identifier);
   bool getBoolParameter(const std::string& identifier);

   class StochOptions : public qpgen_options::QpGenOptions
   {

   private:
      friend void setIntParameter(const std::string& identifier, int value);
      friend void setDoubleParameter(const std::string& identifier, double value);
      friend void setBoolParameter(const std::string& identifier, bool value);
      friend void setOptions(const std::string& opt_file);
      friend void activateHierarchialApproach();
      friend int getIntParameter(const std::string& identifier);
      friend double getDoubleParameter(const std::string& identifier);
      friend bool getBoolParameter(const std::string& identifier);

      static StochOptions& getInstance()
      {
         static StochOptions opt;
         return opt;
      }

      void setHierarchical();
      void setDefaults() override;
      void setPresolveDefaults();
      StochOptions();

      ~StochOptions() override = default;
   };
}

#endif /* PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_ */
