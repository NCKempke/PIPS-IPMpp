//
// Created by charlie on 26.04.21.
//

#ifndef MEHROTRAHEURISTIC_H
#define MEHROTRAHEURISTIC_H

#include <memory>
#include "Statistics.hpp"
#include "TerminationStatus.h"

class Problem;

class Variables;

class Residuals;

class AbstractLinearSystem;

class DistributedFactory;

class Scaler;

enum MehrotraHeuristic { PRIMAL, PRIMAL_DUAL };

class MehrotraStrategy {
public:
   MehrotraStrategy(DistributedFactory& factory, Problem& problem, MehrotraHeuristic mehrotra_heuristic, const Scaler* scaler);
   TerminationStatus corrector_predictor(DistributedFactory& factory, Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
         AbstractLinearSystem& linear_system);
   void set_BiCGStab_tolerance(int iteration) const;
   ~MehrotraStrategy();

protected:
   MehrotraHeuristic mehrotra_heuristic;
   const Scaler* scaler{};
   Variables* corrector_step;
   /** storage for residual vectors */
   Residuals* corrector_residuals;
   unsigned int n_linesearch_points;
   Variables* temp_step;
   Statistics statistics;

   std::unique_ptr<Residuals> residuals_unscaled{};
   /** termination parameters */
   double mutol{1.e-6};
   double artol{1.e-4};

   /** norm of problem data */
   double dnorm{0.};
   /** norm of original unscaled problem */
   double dnorm_orig{0.};

   /* observer stuff for checking convergence of BiCGStab */
   bool bicgstab_skipped;
   bool bicgstab_converged;
   double bigcstab_norm_res_rel;
   int bicg_iterations;

   /** should a BiCGDependent Schedule for the number of Gondzio Correctors be used */
   const bool dynamic_corrector_schedule;
   /** should additional corrector steps for small complementarity pairs be applied */
   const bool additional_correctors_small_comp_pairs;
   /** should additional corrector steps for small complementarity pairs be applied */
   const int max_additional_correctors;
   /** first iteration at which to look for small corrector steps */
   const int first_iter_small_correctors;
   /** alpha must be lower equal to this value for the IPM to try and apply small corrector steps */
   const double max_alpha_small_correctors;
   int NumberSmallCorrectors;

   /** maximum number of Gondzio corrector steps */
   int maximum_correctors;
   /** actual number of Gondzio corrections needed */
   int number_gondzio_corrections;

   /** various parameters associated with Gondzio correction */
   double step_factor0, step_factor1, acceptance_tolerance, beta_min, beta_max;

   // controls whether setBiCGTol applies an dynamic schedule for the BiCGStab tolerance or just uses the user defined input (OUTER_BICG_TOL)
   bool dynamic_bicg_tol;

   /** exponent in Mehrotra's centering parameter, which is usually chosen to me (muaff/mu)^tsig, where muaff is the predicted
 *  complementarity gap obtained from an affine-scaling step, while mu is the current complementarity gap */
   double tsig;

   /** parameters associated with the step length heuristic */
   double gamma_f{0.99};
   double gamma_a{1. / (1. - 0.99)};
   /** number in (0,1) with which the step length is multiplied */
   double steplength_factor{0.99999999};

   /** maximum number of iterations allowed */
   int max_iterations{0};

   /** history of values of mu obtained on all iterations to date */
   double* mu_history{};
   /** history of values of residual norm obtained on all iterations to date */
   double* residual_norm_history{};

   /** history of values of phi obtained on all iterations to date */
   double* phi_history{};
   /** the i-th entry of this array contains the minimum value of phi encountered by the algorithm on or before iteration i */
   double* phi_min_history{};

   bool print_timestamp{true};
   double start_time{-1.};

   TerminationStatus
   corrector_predictor_primal(DistributedFactory& factory, Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
         AbstractLinearSystem& linear_system);
   TerminationStatus
   corrector_predictor_primal_dual(DistributedFactory& factory, Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
         AbstractLinearSystem& linear_system);
   void gondzio_correction_loop_primal(Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
         AbstractLinearSystem& linear_system, int iteration, double& alpha, double sigma, double mu, bool& small_corr, bool& numerical_troubles);
   void compute_predictor_step(Problem& problem, Variables& iterate, Residuals& residuals, AbstractLinearSystem& linear_system, Variables& step);
   void compute_corrector_step(Problem& problem, Variables& iterate, AbstractLinearSystem& linear_system, Variables& step, double sigma, double mu);
   void
   compute_gondzio_corrector(Problem& problem, Variables& iterate, AbstractLinearSystem& linear_system, double rmin, double rmax, bool small_corr);
   void calculate_alpha_weight_candidate(Variables* iterate, Variables* predictor_step, Variables* corrector_step, double alpha_predictor,
         double& alpha_candidate, double& weight_candidate);
   void
   calculate_alpha_pd_weight_candidate(Variables* iterate, Variables* predictor_step, Variables* corrector_step, double alpha_primal, double alpha_dual,
         double& alpha_primal_candidate, double& alpha_dual_candidate, double& weight_primal_candidate, double& weight_dual_candidate);
   void do_probing(Problem* problem, Variables* iterate, Residuals* residuals, Variables* step, double& alpha);
   void do_probing(Problem* problem, Variables* iterate, Residuals* residuals, Variables* step, double& alpha_primal, double& alpha_dual);
   bool is_poor_step(bool& pure_centering_step, bool precond_decreased, double alpha_max) const;
   void compute_probing_step(Variables* probing_step, const Variables* iterate, const Variables* step, double alpha) const;
   void compute_probing_step(Variables* probing_step, const Variables* iterate, const Variables* step, double alpha_primal, double alpha_dual) const;
   double compute_step_factor_probing(double resids_norm_last, double resids_norm_probing, double mu_last, double mu_probing) const;
   bool decrease_preconditioner_impact(AbstractLinearSystem* sys) const;
   void adjust_limit_gondzio_correctors();
   void check_numerical_troubles(Residuals* residuals, bool& numerical_troubles, bool& small_corr) const;
   void
   print_statistics(const Problem* problem, const Variables* iterate, const Residuals* residuals, double dnorm, double alpha, double sigma, int i,
         double mu, int stop_code, int level);
   void print_statistics(const Problem* problem, const Variables* iterate, const Residuals* residuals, double dnorm, double alpha_primal,
         double alpha_dual, double sigma, int i, double mu, int stop_code, int level);
   double mehrotra_step_length(Variables* iterate, Variables* step);
   void mehrotra_step_length(Variables* iterate, Variables* step, double& alpha_primal, double& alpha_dual);
   TerminationStatus compute_status(const Problem* data, const Variables* iterate /* iterate */, const Residuals* residuals, int iteration, double mu);
   void set_problem_norm(const Problem& problem);
   std::pair<double, double> compute_unscaled_gap_and_residual_norm(const Residuals& residuals);
   void default_monitor(const Problem* problem /* problem */, const Variables* iterate /* iterate */, const Residuals* residuals, double alpha,
         double sigma, int i, double mu, int status_code, int level) const;
};


#endif //MEHROTRAHEURISTIC_H
