/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENLINSYS
#define QPGENLINSYS

#include "AbstractLinearSystem.h"
#include "Vector.hpp"
#include "Observer.h"
#include "RegularizationStrategy.h"

#include <functional>
#include <memory>

class Problem;

class DistributedFactory;

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
   enum class IterativeSolverSolutionStatus : int {
      DID_NOT_RUN = -2, NOT_CONVERGED_MAX_ITERATIONS = -1, CONVERGED = 0, SKIPPED = 1, STAGNATION = 3, BREAKDOWN = 4, DIVERGED = 5
   };

   friend std::ostream& operator<<(std::ostream& os, IterativeSolverSolutionStatus status) {
      switch (status) {
         case IterativeSolverSolutionStatus::DID_NOT_RUN :
            return os << "did not run";
         case IterativeSolverSolutionStatus::NOT_CONVERGED_MAX_ITERATIONS:
            return os << "not converged int max iterations";
         case IterativeSolverSolutionStatus::CONVERGED:
            return os << "converged";
         case IterativeSolverSolutionStatus::SKIPPED:
            return os << "skipped";
         case IterativeSolverSolutionStatus::STAGNATION:
            return os << "stagnation occurred";
         case IterativeSolverSolutionStatus::BREAKDOWN:
            return os << "breakdown occurred";
         case IterativeSolverSolutionStatus::DIVERGED:
            return os << "diverged";
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

   [[nodiscard]] int getIntValue(const std::string& s) const override;

   [[nodiscard]] double getDoubleValue(const std::string& s) const override;

   [[nodiscard]] bool getBoolValue(const std::string& s) const override;


   /** - (T_n^-1 \Lambda_n + U_n^-1 \Pi_n)^-1*/
   std::shared_ptr<Vector<double>> nomegaInv{};

   DistributedFactory* factory;

   /** right-hand side of the system */
   std::shared_ptr<Vector<double>> rhs{};

   // TODO : add parameters
   /** regularization parameters */
   const bool apply_regularization{};
   std::unique_ptr<RegularizationStrategy> regularization_strategy;

   std::shared_ptr<Vector<double>> primal_regularization_diagonal{};
   std::shared_ptr<Vector<double>> dual_equality_regularization_diagonal{};
   std::shared_ptr<Vector<double>> dual_inequality_regularization_diagonal{};

   LinearSystem(DistributedFactory* factory_, const Problem& problem, bool create_iter_ref_vecs);

   /** dq = diag(Q); dd = dq - gamma/ v + phi/w */
   std::shared_ptr<Vector<double>> primal_diagonal{};
   std::shared_ptr<Vector<double>> dq{};

   /** index matrices for the upper and lower bounds on x and Cx */
   std::shared_ptr<Vector<double>> ixupp{};
   std::shared_ptr<Vector<double>> icupp{};
   std::shared_ptr<Vector<double>> ixlow{};
   std::shared_ptr<Vector<double>> iclow{};

   /** dimensions of the upper and lower bound vectors */
   long long nxupp{0};
   long long nxlow{0};
   long long mcupp{0};
   long long mclow{0};

   bool useRefs{false};

   /** Work vectors for iterative refinement of the XYZ linear system */
   std::unique_ptr<Vector<double>> sol{};
   std::unique_ptr<Vector<double>> sol2{};
   std::unique_ptr<Vector<double>> res{};
   std::unique_ptr<Vector<double>> resx{};
   std::unique_ptr<Vector<double>> resy{};
   std::unique_ptr<Vector<double>> resz{};

   /** Work vectors for BiCGStab */
   std::unique_ptr<Vector<double>> sol3{};
   std::unique_ptr<Vector<double>> res2{};
   std::unique_ptr<Vector<double>> res3{};
   std::unique_ptr<Vector<double>> res4{};
   std::unique_ptr<Vector<double>> res5{};

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

   const Problem& problem;

public:
   LinearSystem(DistributedFactory* factory, const Problem& problem);

   LinearSystem(DistributedFactory* factory_, const Problem& problem, std::shared_ptr<Vector<double>> dd_, std::shared_ptr<Vector<double>> dq_, std::shared_ptr<Vector<double>> nomegaInv_,
         std::shared_ptr<Vector<double>> primal_regularization_, std::shared_ptr<Vector<double>> dual_equality_regularization, std::shared_ptr<Vector<double>> dual_inequality_regularization_,
         std::shared_ptr<Vector<double>> rhs_, bool create_iter_ref_vecs);

   ~LinearSystem() override = default;


   /** sets up the matrix for the main linear system in "augmented
    * system" form. The actual factorization is performed by a routine
    * specific to either the sparse or dense case.
    *
    * @see QpGenSparseLinsys::factorize
    * @see QpGenDenseLinsys::factorize
    */
   void factorize(Variables& iterate) override;

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
   void solve(Variables& variables, Residuals& residuals, Variables& step) override;

   /** assembles a single vector object from three given vectors
    *
    * @param rhs (output) final joined vector
    * @param rhs1 (input) first part of rhs
    * @param rhs2 (input) middle part of rhs
    * @param rhs3 (input) last part of rhs
    */
   static void joinRHS(Vector<double>& rhs, const Vector<double>& rhs1, const Vector<double>& rhs2, const Vector<double>& rhs3);

   /** extracts three component vectors from a given aggregated vector.
    *
    * @param vars (input) aggregated vector
    * @param vars1 (output) first part of vars
    * @param vars2 (output) middle part of vars
    * @param vars3 (output) last part of vars
    */
   static void separateVars(Vector<double>& vars1, Vector<double>& vars2, Vector<double>& vars3, const Vector<double>& vars);

   /** assemble right-hand side of augmented system and call
       solveCompressed to solve it */
   void solveXYZS(Vector<double>& stepx, Vector<double>& stepy, Vector<double>& stepz, Vector<double>& steps);

   /** perform the actual solve using the factors produced in factor.
    *
    * @param rhs on input contains the aggregated right-hand side of
    * the augmented system; on output contains the solution in
    * aggregated form
    *
    * @see QpGenSparseLinsys::solveCompressed
    * @see QpGenDenseLinsys::solveCompressed
    */
   virtual void solveCompressed(Vector<double>& rhs) = 0;

   /** places the diagonal resulting from the bounds on x into the
    * augmented system matrix */
   virtual void put_primal_diagonal() = 0;

   /** set diagonal in kkt system associated with equalities to zero (necessary after regularization) */
   virtual void clear_dual_equality_diagonal() = 0;

   /** places the diagonal resulting from the bounds on Cx into the
    * augmented system matrix */
   virtual void put_dual_inequalites_diagonal() = 0;

   /** communicate barrier down the linear system hierarchy */
   virtual void put_barrier_parameter(double barrier) = 0;

   /** computes the diagonal matrices in the augmented system from the
       current set of variables */
   virtual void
   computeDiagonals(Vector<double>& t, Vector<double>& lambda, Vector<double>& u, Vector<double>& pi, Vector<double>& v, Vector<double>& gamma,
         Vector<double>& w, Vector<double>& phi);

   /** will factorize and regularize kkt until inertia criterion is met */
   virtual void factorize_with_correct_inertia() = 0;

   void print_regularization_statistics() const;
protected:
   void compute_regularized_system_residuals(const Vector<double>& sol, Vector<double>& res, Vector<double>& solx, Vector<double>& soly,
         Vector<double>& solz, const Problem& problem);

   void compute_system_residuals(const Vector<double>& sol, Vector<double>& res, Vector<double>& solx, Vector<double>& soly, Vector<double>& solz,
         const Problem& problem);
   //void computeResidualsReducedSlacks(const QP& data);

   //void computeResidualsFull(const QP& data);

   void system_mult(double beta, Vector<double>& res, double alpha, const Vector<double>& sol, const Problem& problem, Vector<double>& solx,
         Vector<double>& soly, Vector<double>& solz, bool use_regularized_system);

   double matXYZinfnorm(const Problem& problem, Vector<double>& solx, Vector<double>& soly, Vector<double>& solz, bool use_regularized_system);

   //void matReducedInfnorm(const QP& data);

   //void matFullInfnorm(const QP& data);
   void printDiagonalNorms() const;

   // TODO : move to LinearSystem level
   void solveCompressedBiCGStab(const std::function<void(double, Vector<double>&, double, Vector<double>&)>& matMult,
         const std::function<double()>& matInfnorm);

   void solveCompressedIterRefin(const std::function<void(Vector<double>& sol, Vector<double>& res)>& computeResidual);

};

#endif
