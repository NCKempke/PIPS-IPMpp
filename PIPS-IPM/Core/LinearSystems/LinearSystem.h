/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENLINSYS
#define QPGENLINSYS

#include "AbstractLinearSystem.h"
#include "DoubleMatrixHandle.h"
#include "OoqpVector.h"
#include "Observer.h"
#include "RegularizationStrategy.h"

#include <functional>
#include <memory>

class Problem;

class ProblemFactory;

class Variables;

class Residuals;

/** 
 * Linear System solvers. This class
 * contains definitions of methods and data common to the sparse and
 * dense special cases of the general formulation. The derived classes
 * QpGenSparseLinsys and QpGenDenseLinsys contain the aspects that are
 * specific to the sparse and dense forms.
 *
 * @see QpGenSparseLinsys
 * @see QpGenDenseLinsys
 *
 * @ingroup QpGen 
 */

class LinearSystem : public AbstractLinearSystem, public Subject {

public:
      enum class IterativeSolverSolutionStatus : int
      {
         DID_NOT_RUN = -2,
         NOT_CONVERGED_MAX_ITERATIONS = -1,
         CONVERGED = 0,
         SKIPPED = 1,
         STAGNATION = 3,
         BREAKDOWN = 4,
         DIVERGED = 5
      };

      friend std::ostream& operator<<(std::ostream& os, IterativeSolverSolutionStatus status)
      {
          switch (status)
          {
              case IterativeSolverSolutionStatus::DID_NOT_RUN : return os << "did not run" ;
              case IterativeSolverSolutionStatus::NOT_CONVERGED_MAX_ITERATIONS: return os << "not converged int max iterations";
              case IterativeSolverSolutionStatus::CONVERGED: return os << "converged";
              case IterativeSolverSolutionStatus::SKIPPED: return os << "skipped";
              case IterativeSolverSolutionStatus::STAGNATION: return os << "stagnation occurred";
              case IterativeSolverSolutionStatus::BREAKDOWN: return os << "breakdown occurred";
              case IterativeSolverSolutionStatus::DIVERGED: return os << "diverged";
              // omit default case to trigger compiler warning for missing cases
          };
          return os << static_cast<std::uint16_t>(status);
      }

protected:
   /** observer pattern for convergence status of BiCGStab when calling solve */
   IterativeSolverSolutionStatus bicg_conv_flag{IterativeSolverSolutionStatus::DID_NOT_RUN};
   int bicg_niterations{-1};

   double bicg_resnorm{0.0};
   double bicg_relresnorm{0.0};

   int getIntValue(const std::string& s) const override;

   double getDoubleValue(const std::string& s) const override;

   bool getBoolValue(const std::string& s) const override;

   /** stores a critical diagonal matrix as a vector */
   OoqpVector* nomegaInv{};

   ProblemFactory* factory{};

   /** right-hand side of the system */
   OoqpVector* rhs{};

   // TODO : add parameters
   /** regularization parameters */
   const bool apply_regularization{};
   std::unique_ptr<RegularizationStrategy> regularization_strategy;

   OoqpVector* primal_regularization_diagonal{};
   OoqpVector* dual_equality_regularization_diagonal{};
   OoqpVector* dual_inequality_regularization_diagonal{};

   LinearSystem(ProblemFactory* factory_, Problem* problem, bool create_iter_ref_vecs);

   /** dimensions of the vectors */
   long long nx{0};
   long long my{0};
   long long mz{0};

   /** dq = diag(Q); dd = dq - gamma/ v + phi/w */
   OoqpVector* primal_diagonal{};
   OoqpVector* dq{};

   /** index matrices for the upper and lower bounds on x and Cx */
   OoqpVector* ixupp{};
   OoqpVector* icupp{};
   OoqpVector* ixlow{};
   OoqpVector* iclow{};

   /** dimensions of the upper and lower bound vectors */
   long long nxupp{0};
   long long nxlow{0};
   long long mcupp{0};
   long long mclow{0};

   bool useRefs{false};

   /** Work vectors for iterative refinement of the XYZ linear system */
   OoqpVector* sol{};
   OoqpVector* sol2{};
   OoqpVector* res{};
   OoqpVector* resx{};
   OoqpVector* resy{};
   OoqpVector* resz{};

   /** Work vectors for BiCGStab */
   OoqpVector* sol3{};
   OoqpVector* res2{};
   OoqpVector* res3{};
   OoqpVector* res4{};
   OoqpVector* res5{};

   /// error absorbtion in linear system outer level
   const int outerSolve;
   const int innerSCSolve;

   /// parameters for the bicg solve
   const bool outer_bicg_print_statistics;

   const double outer_bicg_eps;

   const int outer_bicg_max_iter;
   const int outer_bicg_max_normr_divergences;
   const int outer_bicg_max_stagnations;

   const bool xyzs_solve_print_residuals;

   double barrier_parameter_current_iterate{std::numeric_limits<double>::infinity()};
public:
   LinearSystem(ProblemFactory* factory, Problem* problem);

   LinearSystem(ProblemFactory* factory_, Problem* problem, OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_, OoqpVector* primal_regularization_,
        OoqpVector* dual_equality_regularization, OoqpVector* dual_inequality_regularization_, OoqpVector* rhs_, bool create_iter_ref_vecs);

   ~LinearSystem() override;


   /** sets up the matrix for the main linear system in "augmented
    * system" form. The actual factorization is performed by a routine
    * specific to either the sparse or dense case.
    *
    * @see QpGenSparseLinsys::factorize
    * @see QpGenDenseLinsys::factorize
    */
   void factorize(Problem* problem, Variables* iterate) override;

   /** solves the system for a given set of residuals. Assembles the
    * right-hand side appropriate to the matrix factored in factor,
    * solves the system using the factorization produced there,
    * partitions the solution vector into step components, then
    * recovers the step components eliminated during the block
    * elimination that produced the augmented system form
    *
    * @see QpGenSparseLinsys::solveCompressed
    * @see QpGenDenseLinsys::solveCompressed
 */
   void solve(Problem* problem, Variables* variables, Residuals* residuals, Variables* step) override;

   /** assembles a single vector object from three given vectors
    *
    * @param rhs (output) final joined vector
    * @param rhs1 (input) first part of rhs
    * @param rhs2 (input) middle part of rhs
    * @param rhs3 (input) last part of rhs
    */
   virtual void joinRHS(OoqpVector& rhs, const OoqpVector& rhs1, const OoqpVector& rhs2, const OoqpVector& rhs3) const;

   /** extracts three component vectors from a given aggregated vector.
    *
    * @param vars (input) aggregated vector
    * @param vars1 (output) first part of vars
    * @param vars2 (output) middle part of vars
    * @param vars3 (output) last part of vars
    */
   virtual void separateVars(OoqpVector& vars1, OoqpVector& vars2, OoqpVector& vars3, const OoqpVector& vars) const;

   /** assemble right-hand side of augmented system and call
       solveCompressed to solve it */
   virtual void solveXYZS(OoqpVector& stepx, OoqpVector& stepy, OoqpVector& stepz, OoqpVector& steps, OoqpVector& ztemp, Problem* problem);

   /** perform the actual solve using the factors produced in factor.
    *
    * @param rhs on input contains the aggregated right-hand side of
    * the augmented system; on output contains the solution in
    * aggregated form
    *
    * @see QpGenSparseLinsys::solveCompressed
    * @see QpGenDenseLinsys::solveCompressed
    */
   virtual void solveCompressed(OoqpVector& rhs) = 0;

   /** places the diagonal resulting from the bounds on x into the
    * augmented system matrix */
   virtual void put_primal_diagonal() = 0;

   /** places the diagonal resulting from the bounds on Cx into the
    * augmented system matrix */
   virtual void put_dual_inequalites_diagonal() = 0;

   /** computes the diagonal matrices in the augmented system from the
       current set of variables */
   virtual void computeDiagonals(OoqpVector& t, OoqpVector& lambda, OoqpVector& u, OoqpVector& pi, OoqpVector& v,
         OoqpVector& gamma, OoqpVector& w, OoqpVector& phi);

   /** will factorize and regularize kkt until inertia criterion is met */
   virtual void factorize_with_correct_inertia() = 0;

   void print_regularization_statistics() const;

protected:
   void compute_regularized_system_residuals(const OoqpVector& sol, OoqpVector& res, OoqpVector& solx, OoqpVector& soly, OoqpVector& solz, const Problem& problem);
    void compute_system_residuals(const OoqpVector& sol, OoqpVector& res, OoqpVector& solx, OoqpVector& soly, OoqpVector& solz, const Problem& problem);
   //void computeResidualsReducedSlacks(const QP& data);

   //void computeResidualsFull(const QP& data);

   void system_mult(double beta, OoqpVector& res, double alpha, const OoqpVector& sol, const Problem& problem, OoqpVector& solx,
      OoqpVector& soly, OoqpVector& solz, bool use_regularized_system);

   double matXYZinfnorm(const Problem& problem, OoqpVector& solx, OoqpVector& soly, OoqpVector& solz, bool use_regularized_system);

   //void matReducedInfnorm(const QP& data);

   //void matFullInfnorm(const QP& data);
   void printDiagonalNorms() const;

   // TODO : move to LinearSystem level
   void solveCompressedBiCGStab(const std::function<void(double, OoqpVector&, double, OoqpVector&)>& matMult, const std::function<double()>& matInfnorm);

   void solveCompressedIterRefin(const std::function<void(OoqpVector& sol, OoqpVector& res)>& computeResidual);

};

#endif
