/*
 * QpGenOptions.h
 *
 *  Created on: 01.07.2020
 *      Author: bzfkempk
 */

#include "QpGenOptions.h"

namespace qpgen_options {
   QpGenOptions::QpGenOptions() {
      /* initializes base class options first ny calling base class default constructor */
      Options::getInstance();
      /* override with own set of options */
      setDefaults();
   }

   void QpGenOptions::setDefaults() {
      /// SCALER
      bool_options["SCALER_OUTPUT"] = true;

      /// INTERIOR-POINT ALGORITHM
      bool_options["IP_ACCURACY_REDUCED"] = false;
      bool_options["IP_PRINT_TIMESTAMP"] = false;
      bool_options["IP_STEPLENGTH_CONSERVATIVE"] = false;

      /// GONDZIO SOLVER
      /** maximum of Gondzio correctors computed */
      int_options["GONDZIO_MAX_CORRECTORS"] = 3;

      /// SOLVER CONTROLS

      /// ERROR ABSORBTION / ITERATIVE REFINEMENT
      // controls the type of error absorbtion at the outer level of the linear system
      // - 0:no error absortion (OOQP works just fine)
      // - 1:iterative refinement (used when error absortion is
      // also done at a lower level, for example in the solve with
      // the dense Schur complement
      // - 2:BiCGStab with the factorization as preconditioner
      int_options["OUTER_SOLVE"] = 0;
      // controls the type of error absortion/correction at the inner level when solving
      //with the dense Schur complement
      // - 0: no error correction
      // - 1: iter. refin.
      // - 2: BiCGStab
      int_options["INNER_SC_SOLVE"] = 0;

      /// OUTER BIGCSTAB
      double_options["OUTER_BICG_TOL"] = 1e-10;
      double_options["OUTER_BICG_EPSILON"] = 1e-15;

      bool_options["OUTER_BICG_DYNAMIC_TOL"] = true;
      bool_options["OUTER_BICG_PRINT_STATISTICS"] = false;

      int_options["OUTER_BICG_MAX_ITER"] = 75;
      int_options["OUTER_BICG_MAX_NORMR_DIVERGENCES"] = 4;
      int_options["OUTER_BICG_MAX_STAGNATIONS"] = 4;

      bool_options["XYZS_SOLVE_PRINT_RESISDUAL"] = false;

      /// REGULARIZATION FOR LINEAR SYSTEM
      bool_options["REGULARIZATION"] = true;

      double_options["REGULARIZATION_INITIAL_PRIMAL"] = 1e-10;
      double_options["REGULARIZATION_INITIAL_DUAL_Y"] = -1e-10;
      double_options["REGULARIZATION_INITIAL_DUAL_Z"] = -1e-10;

      double_options["REGULARIZATION_MIN_PRIMAL"] = 1e-20;
      double_options["REGULARIZATION_MIN_DUAL"] = -1e-20;
   }
}
