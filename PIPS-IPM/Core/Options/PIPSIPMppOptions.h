/*
 * DistributedOptions.h
 *
 *  Created on: 03.04.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_
#define PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_

#include "Options.h"
#include "pipsport.h"

#include <cassert>
#include <vector>
#include <ostream>

/**
 * The getInstanceMethod must be specified in this BaseClass.
 * Defines default options for StochPIPS.
 */


enum SolverType {
   SOLVER_NONE = 0, SOLVER_MA27 = 1, SOLVER_MA57 = 2, SOLVER_PARDISO = 3, SOLVER_MKL_PARDISO = 4, SOLVER_MUMPS = 5
};

enum SolverTypeDense {
   SOLVER_DENSE_SYM_INDEF = 0, SOLVER_DENSE_SYM_INDEF_SADDLE_POINT = 1, SOLVER_DENSE_SYM_PSD = 2
};

std::ostream& operator<<(std::ostream& os, const SolverType solver);
std::ostream& operator<<(std::ostream& os, const SolverTypeDense solver);


namespace pipsipmpp_options {
   bool is_solver_available(SolverType solver);
   SolverType get_solver_root();
   SolverType get_solver_sub_root();
   SolverType get_solver_leaf();
   SolverTypeDense get_solver_dense();

   void activate_hierarchial_approach();

   void set_options(const std::string& opt_file);
   void set_int_parameter(const std::string& identifier, int value);
   void set_double_parameter(const std::string& identifier, double value);
   void set_bool_parameter(const std::string& identifier, bool value);

   int get_int_parameter(const std::string& identifier);
   double get_double_parameter(const std::string& identifier);
   bool get_bool_parameter(const std::string& identifier);

   class PIPSIPMppOptions : public options::Options {

   private:
      friend void set_int_parameter(const std::string& identifier, int value);
      friend void set_double_parameter(const std::string& identifier, double value);
      friend void set_bool_parameter(const std::string& identifier, bool value);
      friend void set_options(const std::string& opt_file);
      friend void activate_hierarchial_approach();
      friend int get_int_parameter(const std::string& identifier);
      friend double get_double_parameter(const std::string& identifier);
      friend bool get_bool_parameter(const std::string& identifier);

      static PIPSIPMppOptions& getInstance() {
         static PIPSIPMppOptions opt;
         return opt;
      }

      static void setHierarchical();
      void setDefaults() override;
      void setPresolveDefaults();
      PIPSIPMppOptions();

      ~PIPSIPMppOptions() override = default;
   };
}

#endif /* PIPS_IPM_CORE_QPSTOCH_STOCHOPTIONS_H_ */
