#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <memory>
#include <Status.h>
#include "Scaler.h"
#include "Residuals.h"

class Problem;

class Variables;

class AbstractLinearSystem;

class Status;

class OoqpMonitor;

class ProblemFactory;

/**  * @defgroup QpSolvers
 *
 * Interior-point QP solvers
 * @{
 */

const unsigned int max_linesearch_points = 50;

/** 
 * Abstract base class for QP solvers.
 */
class Solver {
public:
   Solver(ProblemFactory& problem_formulation, Problem& problem, const Scaler* = nullptr);

   virtual ~Solver();

   /** solve the IPM system */
   virtual void solve_linear_system(Variables& iterate, Problem& problem, Residuals& residuals, Variables& step);

   /** implements the interior-point method for solving the subproblem */
   virtual TerminationCode solve(Problem& problem, Variables& iterate, Residuals& residuals) = 0;

   /** Mehrotra's heuristic to calculate the final step length */
   virtual double mehrotra_step_length(Variables* iterate, Variables* step);

   /** Mehrotra's heuristic to calculate the final step length in primal and dual direction */
   virtual void mehrotra_step_length_PD(Variables* iterate, Variables* step, double& alpha_primal, double& alpha_dual);

   /** perform monitor operation at each interior-point iteration */
   virtual void
   do_monitor(const Problem* data, const Variables* iterate, const Residuals* residuals, double alpha, double sigma, int i, double mu, int stop_code,
         int level);

   /** perform monitor operation at each interior-point iteration */
   virtual void
   do_monitor_Pd(const Problem* data, const Variables* iterate, const Residuals* residuals, double alpha_primal, double alpha_dual, double sigma, int i,
         double mu, int stop_code, int level);

   /** default monitor: prints out one line of information on each interior-point iteration */
   void default_monitor(const Problem* problem, const Variables* iterate, const Residuals* residuals, double alpha, double sigma, int i, double mu,
         int status_code, int level) const;

   /** this method called to test for convergence status at the end of each interior-point iteration */
   virtual TerminationCode
   do_status(const Problem* problem, const Variables* iterate, const Residuals* residuals, int i, double mu, TerminationCode level);

   /** default method for checking status. May be replaced by a user-defined method */
   virtual TerminationCode
   default_status(const Problem* data, const Variables* iterate, const Residuals* residuals, int iteration, double mu, TerminationCode level);

   /** method to add user-defined monitors to the monitor operations performed at each iteration */
   void add_monitor(OoqpMonitor* m);

   /** method to replace the default_status method with a user-defined status checking method */
   void use_status(Status* s) { status = s; }

   /** enables default_monitor as one of the monitors */
   void monitorSelf();

   void setMuTol(double m) { mutol = m; }

   double getMuTol() const { return mutol; }

   void setArTol(double ar) { artol = ar; }

   double getArTol() const { return artol; }

   double dataNorm() const { return dnorm; }

   double dataNormOrig() const { return dnorm_orig; }

   /** returns a pointed to the linear system object stored in this
    *  class */
   AbstractLinearSystem* getLinearSystem() const { return linear_system; };

   /** reset parameters to their default values */
   virtual void reset_parameters() {};

protected:
   OoqpMonitor* itsMonitors{};
   Status* status{};
   const Scaler* scaler{};
   ProblemFactory& factory;
   /**  storage for step vectors */
   Variables* step, * corrector_step;
   /** storage for residual vectors */
   Residuals* corrector_residuals;

   std::unique_ptr<Residuals> residuals_unscaled{};

   /** norm of problem data */
   double dnorm{0.};

   /** norm of original unscaled problem */
   double dnorm_orig{0.};

   /** termination parameters */
   double mutol{1.e-6};
   double artol{1.e-4};

   /** number in (0,1) with which the step length is multiplied */
   double steplength_factor{0.99999999};

   /** parameters associated with the step length heuristic */
   double gamma_f{0.99};
   double gamma_a{1. / (1. - 0.99)};

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

   AbstractLinearSystem* linear_system{};

   /** iteration counter */
   int iteration{0};

   /** initialize dnorm and dnorm_orig */
   void set_problem_norm(const Problem& problem);

   /** parameter in range [0,100] determines verbosity. (Higher value
    *  => more verbose.) */
   int print_level;

   /** exponent in Mehrotra's centering parameter, which is usually chosen to me (muaff/mu)^tsig, where muaff is the predicted
    *  complementarity gap obtained from an affine-scaling step, while mu is the current complementarity gap */
   double tsig;

   /** maximum number of Gondzio corrector steps */
   int maximum_correctors;

   /** actual number of Gondzio corrections needed */
   int number_gondzio_corrections;

   /** various parameters associated with Gondzio correction */
   double step_factor0, step_factor1, acceptance_tolerance, beta_min, beta_max;

private:
   std::pair<double, double> compute_unscaled_gap_and_residual_norm(const Residuals& residuals);
};

//@}

#endif

