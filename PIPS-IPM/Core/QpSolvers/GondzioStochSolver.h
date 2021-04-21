/*
 * GondzioStochSolver.h
 *
 *  Created on: Dec 7, 2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_GONDZIOSTOCHSOLVER_H
#define PIPS_IPM_GONDZIOSTOCHSOLVER_H

#include "Solver.h"
#include "Observer.h"

class Problem;

class Variables;

class ProblemFormulation;

/**
 * Derived class of Solver implementing Gondzio-correction version of
 * Mehrotra's original predictor-corrector algorithm.
 * @ingroup QpSolvers
 */
class GondzioStochSolver : public Solver, public Observer {
private:
   // returns Gondzio weight for corrector step
   virtual void calculateAlphaWeightCandidate(Variables* iterate, Variables* predictor_step, Variables* corrector_step, double predictor_alpha,
         double& alpha_candidate, double& weight_candidate);

protected:
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

   void computePredictorStep(Problem* prob, Variables* iterate, Residuals* residuals);

   void computeCorrectorStep(Problem* prob, Variables* iterate, double sigma, double mu);

   void computeGondzioCorrector(Problem* prob, Variables* iterate, double rmin, double rmax, bool small_corr);

   void checkLinsysSolveNumericalTroublesAndReact(Residuals* residuals, bool& numerical_troubles, bool& small_corr) const;

   void registerBiCGStabOvserver(LinearSystem* sys);

   void setBiCGStabTol(int iteration) const;

   // controls whether setBiCGTol applies an dynamic schedule for the BiCGStab tolerance or just uses the user defined input (OUTER_BICG_TOL)
   bool dynamic_bicg_tol;

   void adjustLimitGondzioCorrectors();

   bool decreasePreconditionerImpact(LinearSystem* sys) const;

   double computeStepFactorProbing(double resids_norm_last, double resids_norm_probing, double mu_last, double mu_probing) const;

   void computeProbingStep(Variables* probing_step, const Variables* iterate, const Variables* step, double alpha) const;

   void doProbing(Problem* prob, Variables* iterate, Residuals* residuals, double& alpha);

   bool restartIterateBecauseOfPoorStep(bool& pure_centering_step, bool precond_limit, double alpha_max) const;

   void pushConvergedVarsAwayFromBounds(Problem& problem, Variables& vars) const;

   void pushSmallComplementarityProducts(const Problem& problem, Variables& iterate, Residuals& residuals) const;

public:

   GondzioStochSolver(ProblemFormulation& problem_formulation, Problem& problem, const Scaler* scaler = nullptr);

   virtual ~GondzioStochSolver();

   TerminationCode solve(Problem& problem, Variables* iterate, Residuals* residuals) override;

   void notifyFromSubject() override;
};

#endif /* PIPS_IPM_GONDZIOSTOCHSOLVER_H */
