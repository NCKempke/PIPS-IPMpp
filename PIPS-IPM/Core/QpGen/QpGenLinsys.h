/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENLINSYS
#define QPGENLINSYS

#include "LinearSystem.h"
#include "DoubleMatrixHandle.h"
#include "OoqpVector.h"
#include "Observer.h"

#include <functional>

class Data;
class QpGenData;
class QpGen;
class Variables;
class Residuals;
class DoubleLinearSolver;

/** 
 * Linear System solvers for the general QP formulation. This class
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

class QpGenLinsys : public LinearSystem, public Subject {
protected:
  /** observer pattern for convergence status of BiCGStab when calling solve */
  int bicg_conv_flag{-2};
  int bicg_niterations{-1};

  double bicg_resnorm{0.0};
  double bicg_relresnorm{0.0};

  int getIntValue(const std::string& s) const override;
  double getDoubleValue(const std::string& s) const override;
  bool getBoolValue(const std::string& s) const override;

  /** stores a critical diagonal matrix as a vector */
  OoqpVector* nomegaInv{};

  QpGen* factory{};

  /** right-hand side of the system */
  OoqpVector* rhs{};

  QpGenLinsys( QpGen* factory_, QpGenData* prob, bool create_iter_ref_vecs );

  /** dimensions of the vectors in the general QP formulation */
  long long nx{0};
  long long my{0};
  long long mz{0};

  /** temporary storage vectors */
  OoqpVector* dd{};
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

  int useRefs{0};

  /** Work vectors for iterative refinement of the XYZ linear system */
  OoqpVector* sol{};
  OoqpVector* res{};
  OoqpVector* resx{};
  OoqpVector* resy{};
  OoqpVector* resz{};

  /** Work vectors for BiCGStab */
  OoqpVector* sol2{};
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

public:
  QpGenLinsys( QpGen* factory, QpGenData* data );
  QpGenLinsys( QpGen* factory_, QpGenData* prob, OoqpVector* dd_, OoqpVector* dq_,
        OoqpVector* nomegaInv_, OoqpVector* rhs_, bool create_iter_ref_vecs );

  ~QpGenLinsys() override;


  /** sets up the matrix for the main linear system in "augmented
   * system" form. The actual factorization is performed by a routine
   * specific to either the sparse or dense case.
   *
   * @see QpGenSparseLinsys::factor
   * @see QpGenDenseLinsys::factor
   */
  virtual void factor(Data *prob, Variables *vars);

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
  virtual void solve(Data *prob, Variables *vars, Residuals *res,
		     Variables *step);

  /** assembles a single vector object from three given vectors
   *
   * @param rhs (output) final joined vector
   * @param rhs1 (input) first part of rhs
   * @param rhs2 (input) middle part of rhs
   * @param rhs3 (input) last part of rhs
   */
  virtual void joinRHS( OoqpVector& rhs, const OoqpVector& rhs1,
			const OoqpVector& rhs2, const OoqpVector& rhs3 ) const;

  /** extracts three component vectors from a given aggregated vector.
   *
   * @param vars (input) aggregated vector
   * @param vars1 (output) first part of vars
   * @param vars2 (output) middle part of vars
   * @param vars3 (output) last part of vars
   */
  virtual void separateVars( OoqpVector& vars1, OoqpVector& vars2,
			     OoqpVector& vars3, const OoqpVector& vars ) const;

  /** assemble right-hand side of augmented system and call
      solveCompressed to solve it */
  virtual void solveXYZS( OoqpVector& stepx, OoqpVector& stepy,
			  OoqpVector& stepz, OoqpVector& steps,
			  OoqpVector& ztemp, QpGenData * data );

  /** perform the actual solve using the factors produced in factor.
   *
   * @param rhs on input contains the aggregated right-hand side of
   * the augmented system; on output contains the solution in
   * aggregated form 
   *
   * @see QpGenSparseLinsys::solveCompressed
   * @see QpGenDenseLinsys::solveCompressed
   */
  virtual void solveCompressed( OoqpVector& rhs ) = 0;

  /** places the diagonal resulting from the bounds on x into the
   * augmented system matrix */
  virtual void putXDiagonal( OoqpVector& xdiag ) = 0;

  /** places the diagonal resulting from the bounds on Cx into the
   * augmented system matrix */
  virtual void putZDiagonal( OoqpVector& zdiag ) = 0;

  /** computes the diagonal matrices in the augmented system from the
      current set of variables */
  virtual void computeDiagonals( OoqpVector& dd, OoqpVector& omega,
				 OoqpVector& t,  OoqpVector& lambda,
				 OoqpVector& u,  OoqpVector& pi,
				 OoqpVector& v,  OoqpVector& gamma,
				 OoqpVector& w,  OoqpVector& phi );
   protected:
      void computeResidualXYZ(const OoqpVector& sol, OoqpVector& res, OoqpVector& solx,
            OoqpVector& soly, OoqpVector& solz, const QpGenData& data);
      void computeResidualsReducedSlacks( const QpGenData& data );
      void computeResidualsFull( const QpGenData& data );

      void matXYZMult(double beta, OoqpVector& res, double alpha, const OoqpVector& sol,
            const QpGenData& data, OoqpVector& solx, OoqpVector& soly,
            OoqpVector& solz);
      void matReducedSlacksMult( const QpGenData& data );
      void matFullMult( const QpGenData& data );

      double matXYZinfnorm(const QpGenData& data, OoqpVector &solx, OoqpVector &soly,
            OoqpVector &solz);
      void matReducedInfnorm( const QpGenData& data );
      void matFullInfnorm( const QpGenData& data );


  // TODO : move to LinearSystem level
  void solveCompressedBiCGStab( const std::function<void(double, OoqpVector&, double, OoqpVector&)>& matMult, const std::function<double()>& matInfnorm );

  void solveCompressedIterRefin( const std::function<void(OoqpVector& sol, OoqpVector& res)>& computeResidual );

};

#endif
