/*
 * StochOptions.C
 *
 *  Created on: 03.04.2020
 *      Author: bzfkempk
 */

#include "PIPSIPMppOptions.h"

#include <limits>
#include <vector>
#include <algorithm>

#include "pipsdef.h"

std::ostream& operator<<(std::ostream& os, const SolverType solver) {
   switch (solver) {
      case SolverType::SOLVER_NONE:
         os << "SOLVER_NONE";
         break;
      case SolverType::SOLVER_MA27:
         os << "SOLVER_MA27";
         break;
      case SolverType::SOLVER_MA57:
         os << "SOLVER_MA57";
         break;
      case SolverType::SOLVER_PARDISO:
         os << "SOLVER_PARDISO";
         break;
      case SolverType::SOLVER_MKL_PARDISO:
         os << "SOLVER_MKL_PARDISO";
         break;
      case SolverType::SOLVER_MUMPS:
         os << "SOLVER_MUMPS";
         break;
   }
   return os;
}

std::ostream& operator<<(std::ostream& os, const SolverTypeDense solver) {
   switch (solver) {
      case SolverTypeDense::SOLVER_DENSE_SYM_INDEF:
         os << "SOLVER_DENSE_SYM_INDEF";
         break;
      case SolverTypeDense::SOLVER_DENSE_SYM_INDEF_SADDLE_POINT:
         os << "SOLVER_DENSE_SYM_INDEF_SADDLE_POINT";
         break;
      case SolverTypeDense::SOLVER_DENSE_SYM_PSD:
         os << "SOLVER_DENSE_SYM_PSD";
         break;
   }
   return os;
}

namespace pipsipmpp_options {
   const std::vector<SolverType> solvers_available{
#ifdef WITH_PARDISO
         SolverType::SOLVER_PARDISO,
#endif
#ifdef WITH_MKL_PARDISO
         SolverType::SOLVER_MKL_PARDISO,
#endif
#ifdef WITH_MA57
         SolverType::SOLVER_MA57,
#endif
#ifdef WITH_MA27
         SolverType::SOLVER_MA27,
#endif
#ifdef WITH_MUMPS
         SolverType::SOLVER_MUMPS,
#endif
         SolverType::SOLVER_NONE};

   bool is_solver_available(SolverType solver) {
      return std::find(solvers_available.begin(), solvers_available.end(), solver) != solvers_available.end();
   }

   void printAvailableSolvers() {
      std::cout << "Available solvers are ";
      for (SolverType s : solvers_available) {
         if (s != SolverType::SOLVER_NONE)
            std::cout << s << " = " << static_cast<int>(s) << "\t";
      }
      std::cout << "\n";
   }

   SolverType get_solver_root() {
      const int solver_int = get_int_parameter("LINEAR_ROOT_SOLVER");
      if (solver_int < 1 || solver_int > 5) {
         if (PIPS_MPIgetRank() == 0)
            std::cout << "Error: unknown solver type LINEAR_ROOT_SOLVER: " << solver_int << "\n";
         printAvailableSolvers();
         MPI_Barrier(MPI_COMM_WORLD);
         MPI_Abort(MPI_COMM_WORLD, -1);
      }

      SolverType solver_root = static_cast<SolverType>(solver_int);
      if (!is_solver_available(solver_root)) {
         if (PIPS_MPIgetRank() == 0) {
            std::cout << "Error: sprecified root solver \"" << solver_root << "\" is not available\n";
            printAvailableSolvers();
         }
         MPI_Barrier(MPI_COMM_WORLD);
         MPI_Abort(MPI_COMM_WORLD, -1);
      }

      return solver_root;
   }

   SolverType get_solver_sub_root() {
      const int solver_int = get_int_parameter("LINEAR_SUB_ROOT_SOLVER");
      if (solver_int < 1 || solver_int > 5) {
         if (PIPS_MPIgetRank() == 0)
            std::cout << "Error: unknown solver type LINEAR_SUB_ROOT_SOLVER: " << solver_int << "\n";
         printAvailableSolvers();
         MPI_Barrier(MPI_COMM_WORLD);
         MPI_Abort(MPI_COMM_WORLD, -1);
      }

      SolverType solver_sub_root = static_cast<SolverType>(solver_int);
      if (!is_solver_available(solver_sub_root)) {
         if (PIPS_MPIgetRank() == 0) {
            std::cout << "Error: sprecified sub-root solver \"" << solver_sub_root << "\" is not available\n";
            printAvailableSolvers();
         }
         MPI_Barrier(MPI_COMM_WORLD);
         MPI_Abort(MPI_COMM_WORLD, -1);
      }

      return solver_sub_root;
   }

   SolverTypeDense get_solver_dense() {
      const int solver_int = get_int_parameter("LINEAR_DENSE_SOLVER");
      if (solver_int < 0 || solver_int > 2) {
         if (PIPS_MPIgetRank() == 0)
            std::cout << "Error: unknown solver type LINEAR_DENSE_SOLVER: " << solver_int << "\n";
         MPI_Barrier(MPI_COMM_WORLD);
         MPI_Abort(MPI_COMM_WORLD, -1);
      }

      SolverTypeDense solver_dense = static_cast<SolverTypeDense>(solver_int);
      return solver_dense;
   }

   SolverType get_solver_leaf() {
      const int solver_int = get_int_parameter("LINEAR_LEAF_SOLVER");
      if (solver_int < 1 || solver_int > 5) {
         if (PIPS_MPIgetRank() == 0)
            std::cout << "Error: unknown solver type LINEAR_LEAF_SOLVER: " << solver_int << "\n";
         printAvailableSolvers();
         MPI_Barrier(MPI_COMM_WORLD);
         MPI_Abort(MPI_COMM_WORLD, -1);
      }

      SolverType solver_leaf = static_cast<SolverType>(solver_int);
      if (!is_solver_available(solver_leaf)) {
         if (PIPS_MPIgetRank() == 0) {
            std::cout << "Error: sprecified leaf solver \"" << solver_leaf << "\" is not available\n";
            printAvailableSolvers();
         }
         MPI_Barrier(MPI_COMM_WORLD);
         MPI_Abort(MPI_COMM_WORLD, -1);
      }

      return solver_leaf;
   }

   void PIPSIPMppOptions::setHierarchical() {
      bool_options["HIERARCHICAL"] = true;

      bool_options["SC_COMPUTE_BLOCKWISE"] = true;
      bool_options["ALLREDUCE_SCHUR_COMPLEMENT"] = true;
      bool_options["PRECONDITION_DISTRIBUTED"] = false;

      int_options["GONDZIO_MAX_CORRECTORS"] = 5;
      int_options["GONDZIO_STOCH_ADDITIONAL_CORRECTORS_MAX"] = 3;

      if (is_solver_available(SolverType::SOLVER_MA57)) {
         int_options["LINEAR_LEAF_SOLVER"] = SolverType::SOLVER_MA57;
         int_options["LINEAR_ROOT_SOLVER"] = SolverType::SOLVER_MA57;
         int_options["LINEAR_SUB_ROOT_SOLVER"] = SolverType::SOLVER_MA57;
      }
   }

   PIPSIPMppOptions::PIPSIPMppOptions() {
      /* initialize base class options first (QpGenOptions) */
      Options::getInstance();

      /* now override with own set of options */
      setDefaults();
   }

   void PIPSIPMppOptions::setDefaults() {
      /// GENERAL
      bool_options["PRINT_TREESIZES_ON_READ"] = false;
      /* surpresses some of the output */
      bool_options["SILENT"] = false;

      /// SCALER
      bool_options["SCALER_OUTPUT"] = true;

      /// LINEAR SOLVERS
      assert(solvers_available.size() > 1);

      const SolverType default_solver = solvers_available[0];

      int_options["LINEAR_LEAF_SOLVER"] = default_solver;
      int_options["LINEAR_ROOT_SOLVER"] = default_solver;
      int_options["LINEAR_SUB_ROOT_SOLVER"] = default_solver;

      int_options["LINEAR_DENSE_SOLVER"] = SolverTypeDense::SOLVER_DENSE_SYM_INDEF;

      bool_options["PARDISO_FOR_GLOBAL_SC"] = true;
      bool_options["PARDISO_SPARSE_RHS_LEAF"] = false;
      /** -1 is choose default */
      int_options["PARDISO_SYMB_INTERVAL"] = -1;
      int_options["PARDISO_PIVOT_PERTURBATION"] = -1;
      int_options["PARDISO_NITERATIVE_REFINS"] = -1;
      int_options["PARDISO_PIVOT_PERTURBATION_ROOT"] = -1;
      int_options["PARDISO_NITERATIVE_REFINS_ROOT"] = -1;

      /// Schur Complement Computation
      /// PRECONDITIONERS

      if (default_solver != SolverType::SOLVER_PARDISO && default_solver != SolverType::SOLVER_MUMPS) {
         bool_options["SC_COMPUTE_BLOCKWISE"] = true;
         bool_options["PRECONDITION_DISTRIBUTED"] = false;
      }
      else {
         bool_options["SC_COMPUTE_BLOCKWISE"] = false;
         bool_options["PRECONDITION_DISTRIBUTED"] = true;
      }
      bool_options["PRECONDITION_SPARSE"] = true;

      int_options["SC_BLOCKWISE_BLOCKSIZE_MAX"] = 20;

      /// HIERARCHICAL APPROACH
      bool_options["HIERARCHICAL"] = false;
      bool_options["HIERARCHICAL_APPLY_SPLIT"] = true;
      bool_options["HIERARCHICAL_TESTING"] = false;
      bool_options["HIERARCHICAL_PRINT_HIER_DATA"] = false;
      bool_options["HIERARCHICAL_MOVE_A0_TO_DENSE_LAYER"] = false;

      /** should all columns be computed at once when computing the hierarchical schur complements ?
       * this will render SC_BLOCKSIZE_HIERARCHICAL obsolete but might consume big amounts of memory and
       * might fail if columns * length columns is too big (stroage allocation might fail)
       */
      bool_options["SC_HIERARCHICAL_COMPUTE_BLOCKWISE"] = true;

      /** columns to solve for at once when computing the hierarchical schur complements*/
      int_options["SC_BLOCKSIZE_HIERARCHICAL"] = 20;

      /** 1 -> only dense border, 2 -> dense border + 1 additional layer ... */
      int_options["HIERARCHICAL_APPROACH_N_LAYERS"] = 2;

      /// SCHUR COMPLEMENT
      /** should the schur complement be allreduced to all processes or to a single one */
      bool_options["ALLREDUCE_SCHUR_COMPLEMENT"] = false;
      /// GONDZIO SOLVERS
      /** should adaptive linesearch be applied in the GondzioStoch solvers */
      bool_options["GONDZIO_STOCH_ADAPTIVE_LINESEARCH"] = false;
      /** if GONDZIO adaptive linesearch is true determines number of linesearch points */
      int_options["GONDZIO_STOCH_N_LINESEARCH"] = 10;
      /** should additional corrector steps for small complementarity pairs be applied */
      bool_options["GONDZIO_STOCH_ADDITIONAL_CORRECTORS_SMALL_VARS"] = true;
      /** how many additional steps should be applied at most (in addition to the still existing gondzio corrector limit) */
      int_options["GONDZIO_STOCH_ADDITIONAL_CORRECTORS_MAX"] = 1;
      /** first iteration at which to look for small corrector steps */
      int_options["GONDZIO_STOCH_FIRST_ITER_SMALL_CORRECTORS"] = 10;
      /** alpha must be lower equal to this value for the IPM to try and apply small corrector steps */
      double_options["GONDZIO_STOCH_MAX_ALPHA_SMALL_CORRECTORS"] = 0.95;
      /** should the amount of gondzio correctors be scheduled dynamically - invalidates the max correctors setting */
      bool_options["GONDZIO_STOCH_USE_DYNAMIC_CORRECTOR_SCHEDULE"] = false;

      /** should relatively early converged variables be pushed away artificially from their bounds */
      bool_options["GONDZIO_STOCH_PUSH_CONVERGED_VARS_FROM_BOUND"] = false;
      /** at which frequency should converged variables be pushed away from their bounds */
      int_options["GONDZIO_STOCH_FREQUENCY_PUSH_CONVERGED_VARS"] = 4;
      /** starting with which mu should the pushing be done */
      double_options["GONDZIO_STOCH_MU_LIMIT_PUSH_CONVERGED_VARS"] = 1e-3;

      /// SOLVER CONTROLS

      /// ERROR ABSORBTION / ITERATIVE REFINEMENT
      // controls the type of error absorbtion at the outer level of the linear system
      // - 0:no error absortion (OOQP works just fine)
      // - 1:iterative refinement (used when error absortion is
      // also done at a lower level, for example in the solve with
      // the dense Schur complement
      // - 2:BiCGStab with the factorization as preconditioner
      int_options["OUTER_SOLVE"] = 2;
      bool_options["OUTER_SOLVE_REFINE_ORIGINAL_SYSTEM"] = false;

      // controls the type of error absortion/correction at the inner level when solving
      //with the dense Schur complement
      // - 0: no error correction
      // - 1: iter. refin.
      // - 2: BiCGStab
      int_options["INNER_SC_SOLVE"] = 0;

      /// OUTER BIGCSTAB
      double_options["OUTER_BICG_TOL"] = 1e-10;
      double_options["OUTER_BICG_EPSILON"] = 1e-15;

      bool_options["OUTER_BICG_PRINT_STATISTICS"] = false;

      int_options["OUTER_BICG_MAX_ITER"] = 75;
      int_options["OUTER_BICG_MAX_NORMR_DIVERGENCES"] = 4;
      int_options["OUTER_BICG_MAX_STAGNATIONS"] = 4;

      bool_options["XYZS_SOLVE_PRINT_RESISDUAL"] = false;

      /// REGULARIZATION FOR LINEAR SYSTEM
      bool_options["REGULARIZATION"] = false;
      bool_options["REGULARIZATION_VERBOSE"] = false;
      int_options["REGULARIZATION_STRATEGY"] = 1; // 0 -> Ipopt, 1 -> FriedlanderOrban

      double_options["FRIEDLANDER_ORBAN_REGULARIZATION_INITIAL_PRIMAL"] = 1;
      double_options["FRIEDLANDER_ORBAN_REGULARIZATION_INITIAL_DUAL_Y"] = 1;
      double_options["FRIEDLANDER_ORBAN_REGULARIZATION_INITIAL_DUAL_Z"] = 1;
      double_options["FRIEDLANDER_ORBAN_REGULARIZATION_PRIMAL_MIN"] = 1e-10;
      double_options["FRIEDLANDER_ORBAN_REGULARIZATION_DUAL_MIN"] = 1e-10;

      double_options["IPOPT_REGULARIZATION_MIN_PRIMAL"] = 1e-20;
      double_options["IPOPT_REGULARIZATION_MAX_PRIMAL"] = 1e40;

      bool_options["SCHUR_COMPLEMENT_FORCE_SPARSE_COMPUTATIONS"] = false;
      setPresolveDefaults();
   }

   void PIPSIPMppOptions::setPresolveDefaults() {
      /** all presolve/postsolve constants and settings */
      // TODO : many of these need adjustments/ have to be thought about
      double_options["PRESOLVE_INFINITY"] = std::numeric_limits<double>::infinity();

      /// STOCH PRESOLVER
      /** limit for max rounds to apply all presolvers */
      int_options["PRESOLVE_MAX_ROUNDS"] = 2;
      /** should the problem be written to std::cout before and after presolve */
      bool_options["PRESOLVE_PRINT_PROBLEM"] = false;
      /** should the presolved problem be written out in MPS format */
      bool_options["PRESOLVE_WRITE_PRESOLVED_PROBLEM_MPS"] = false;
      /** should free variables' bounds be reset after presolve (given the row implying these bounds was not removed */
      bool_options["PRESOLVE_RESET_FREE_VARIABLES"] = false;
      /** should inequalities be transformed into euqalities at the end of presolve */
      bool_options["PRESOLVE_TRANSFROM_INEQUALITIES_INTO_EQUALITIES"] = false;

      /** verbosity */
      int_options["PRESOLVE_VERBOSITY"] = 1;

      /** turn respective presolvers on/off */
      bool_options["PRESOLVE_BOUND_STRENGTHENING"] = true;
      bool_options["PRESOLVE_PARALLEL_ROWS"] = true;
      bool_options["PRESOLVE_COLUMN_FIXATION"] = true;
      bool_options["PRESOLVE_SINGLETON_ROWS"] = true;
      bool_options["PRESOLVE_SINGLETON_COLUMNS"] = true;

      /// BOUND STRENGTHENING
      /** limit for rounds of bound strengthening per call of presolver */
      int_options["PRESOLVE_BOUND_STR_MAX_ITER"] = 2;
      /** min entry to devide by in order to derive a bound */
      double_options["PRESOLVE_BOUND_STR_NUMERIC_LIMIT_ENTRY"] = 1e-7;
      /** max activity to be devided */
      double_options["PRESOLVE_BOUND_STR_MAX_PARTIAL_ACTIVITY"] = std::numeric_limits<double>::max();
      /** max bounds proposed from bounds strengthening presolver */
      double_options["PRESOLVE_BOUND_STR_NUMERIC_LIMIT_BOUNDS"] = 1e12;

      /// COLUMN FIXATION
      /** limit on the possible impact a column can have on the problem */
      double_options["PRESOLVE_COLUMN_FIXATION_MAX_FIXING_IMPACT"] = 1.0e-12; // for variable fixing

      /// MODEL CLEANUP
      /** limit for the size of a matrix entry below which it will be removed from the problem */
      double_options["PRESOLVE_MODEL_CLEANUP_MIN_MATRIX_ENTRY"] = 1.0e-10;//1.0e-10; // for model cleanup // was 1.0e-10
      /** max for the matrix entry when the impact of entry times (bux-blx) is considered */
      double_options["PRESOLVE_MODEL_CLEANUP_MAX_MATRIX_ENTRY_IMPACT"] = 1.0e-3; // was 1.0e-3
      /** difference in orders between feastol and the impact of entry times (bux-blx) for an entry to get removed */
      double_options["PRESOLVE_MODEL_CLEANUP_MATRIX_ENTRY_IMPACT_FEASDIST"] = 1.0e-2;  // for model cleanup // was 1.0e-2

      /// PARALLEL ROWS
      /** tolerance for comparing two double values in two different rows and for them being considered equal */
      double_options["PRESOLVE_PARALLEL_ROWS_TOL_COMPARE_ENTRIES"] = 1.0e-8;

      /// PRESOLVE DATA
      /** absolute maximum for newly found bounds accepted by presolve */
      double_options["PRESOLVE_MAX_BOUND_ACCEPTED"] = 1e10;

      /** track row/column through presolve */
      bool_options["PRESOLVE_TRACK_ROW"] = false;
      bool_options["PRESOLVE_TRACK_COL"] = false;
      /** the row */
      int_options["PRESOLVE_TRACK_ROW_INDEX"] = 0;
      int_options["PRESOLVE_TRACK_ROW_NODE"] = -1;
      int_options["PRESOLVE_TRACK_ROW_SYSTEM"] = 0; // 0 -> EQ, 1 -> INEQ
      bool_options["PRESOLVE_TRACK_ROW_LINKING"] = false;
      /** the column */
      int_options["PRESOLVE_TRACK_COL_INDEX"] = 0;
      int_options["PRESOLVE_TRACK_COL_NODE"] = -1;

      /// POSTSOLVE
      /** should postsolve be applied */
      bool_options["POSTSOLVE"] = true;

      /** tolerance used for checking residuals after postsolve */
      double_options["POSTSOLVE_TOLERANCE"] = feastol * 1e2;

      /** should the residuals before unscaling, after unscaling before postsolve, after postsolve be printed */
      bool_options["POSTSOLVE_PRINT_RESIDS"] = true;

      /// FILTER
      bool_options["FILTER_VERBOSE"] = false;
   }


   void activate_hierarchial_approach() {
      PIPSIPMppOptions::getInstance().setHierarchical();
   }

   void set_options(const std::string& opt_file) {
      return PIPSIPMppOptions::getInstance().load_options_from_file(opt_file);
   }

   void set_int_parameter(const std::string& identifier, int value) {
      PIPSIPMppOptions::getInstance().set_int_param(identifier, value);
   }

   void set_double_parameter(const std::string& identifier, double value) {
      PIPSIPMppOptions::getInstance().set_double_param(identifier, value);
   }

   void set_bool_parameter(const std::string& identifier, bool value) {
      PIPSIPMppOptions::getInstance().set_bool_param(identifier, value);
   }

   int get_int_parameter(const std::string& identifier) {
      return PIPSIPMppOptions::getInstance().get_int_param(identifier);
   }

   bool get_bool_parameter(const std::string& identifier) {
      return PIPSIPMppOptions::getInstance().get_bool_param(identifier);
   }

   double get_double_parameter(const std::string& identifier) {
      return PIPSIPMppOptions::getInstance().get_double_param(identifier);
   }

}
