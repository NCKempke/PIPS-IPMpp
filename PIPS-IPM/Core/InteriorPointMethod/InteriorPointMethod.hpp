#ifndef GONDZIO_INTERIOR_POINT_METHOD_H
#define GONDZIO_INTERIOR_POINT_METHOD_H

#include "pipsdef.h"
#include "Solver.h"
#include "Observer.h"
#include "Variables.h"
#include "Problem.h"

enum StepLengthType { PRIMAL, PRIMAL_DUAL };

class InteriorPointMethod : public Solver, public Observer {
public:
   InteriorPointMethod(ProblemFactory& problem_formulation, Problem& problem, const Scaler* scaler);

   virtual ~InteriorPointMethod();

   TerminationCode solve(Problem& problem, Variables& iterate, Residuals& residuals) override;

protected:
   StepLengthType step_length_type;
   unsigned int n_linesearch_points;
   Variables* temp_step;

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

   /* should early converged variables push away a bit from their respecitve bounds */
   const bool push_converged_vars_from_bound;
   /* at which frequnecy should small vars be checked */
   const int fequency_push_converged_vars_from_bound;
   /* starting with which mu should checking start */
   const double mu_limit_push_converged_vars_from_bound;

   /* observer stuff for checking convergence of BiCGStab */
   bool bicgstab_skipped;
   bool bicgstab_converged;
   double bigcstab_norm_res_rel;
   int bicg_iterations;
   // controls whether setBiCGTol applies an dynamic schedule for the BiCGStab tolerance or just uses the user defined input (OUTER_BICG_TOL)
   bool dynamic_bicg_tol;

   void compute_predictor_step(Problem& problem, Variables& iterate, Residuals& residuals);
   void compute_corrector_step(Problem& problem, Variables& iterate, double sigma, double mu);
   void compute_gondzio_corrector(Problem* problem, Variables* iterate, double rmin, double rmax, bool small_corr);
   void check_linear_system_solve_numerical_troubles(Residuals* residuals, bool& numerical_troubles, bool& small_corr) const;
   void registerBiCGStabOvserver(LinearSystem* sys);
   void setBiCGStabTol(int iteration) const;
   void adjustLimitGondzioCorrectors();
   bool decreasePreconditionerImpact(LinearSystem* sys) const;
   double computeStepFactorProbing(double resids_norm_last, double resids_norm_probing, double mu_last, double mu_probing) const;
   void computeProbingStep_pd(Variables* probing_step, const Variables* iterate, const Variables* step, double alpha_primal, double alpha_dual) const;
   void doProbing_pd(Problem* prob, Variables* iterate, Residuals* resid, double& alpha_pri, double& alpha_dual);
   bool restartIterateBecauseOfPoorStep(bool& pure_centering_step, bool precond_limit, double alpha_max) const;
   void notifyFromSubject() override;

private:
   // returns Gondzio weight for corrector step
   virtual void
   calculate_alpha_PD_weight_candidate(Variables* iterate, Variables* predictor_step, Variables* corrector_step, double alpha_primal, double alpha_dual,
         double& step_length_primal_candidate, double& step_length_dual_candidate, double& weight_primal_candidate, double& weight_dual_candidate);
};

#endif