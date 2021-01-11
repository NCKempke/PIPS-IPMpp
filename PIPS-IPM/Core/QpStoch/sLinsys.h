/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#ifndef SLINSYS
#define SLINSYS

#include "QpGenLinsys.h"
#include "DoubleLinearSolver.h"
#include "OoqpVectorHandle.h"
#include "DenseSymMatrix.h"
#include "SparseSymMatrix.h"
#include "DenseGenMatrix.h"
#include "SimpleVector.h"
#include "StochVector.h"

#include <vector>
#include <memory>

#include "mpi.h"

class sTree;
class sFactory;
class sData;
class StringGenMatrix;
class StringSymMatrix;

class sLinsys : public QpGenLinsys
{
 public:
  sLinsys(sFactory* factory, sData* prob, bool is_hierarchy_root = false);
  sLinsys(sFactory* factory,
		   sData* prob, 
		   OoqpVector* dd, 
		   OoqpVector* dq,
		   OoqpVector* nomegaInv,
		   OoqpVector* rhs,
		   OoqpVector* primal_reg,
		   OoqpVector* dual_y_reg,
		   OoqpVector* dual_z_reg,
		   bool create_iter_ref_vecs
  );

  ~sLinsys() override;

  void factor (Data *prob, Variables *vars) override;
  virtual void factor2(sData *prob, Variables *vars) = 0;


  virtual void Lsolve( sData *prob, OoqpVector& x ) = 0;
  virtual void Dsolve( sData *prob, OoqpVector& x ) = 0;
  virtual void Ltsolve( sData *prob, OoqpVector& x ) = 0;
  virtual void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp) = 0;

  virtual void putZDiagonal( OoqpVector& zdiag )=0;
  virtual void solveCompressed( OoqpVector& rhs );
  virtual void putXDiagonal( OoqpVector& xdiag_ )=0;

  void joinRHS( OoqpVector& rhs_in, const OoqpVector& rhs1_in,
		const OoqpVector& rhs2_in, const OoqpVector& rhs3_in ) const;

  void separateVars( OoqpVector& x_in, OoqpVector& y_in,
		     OoqpVector& z_in, const OoqpVector& vars_in ) const;
  
  virtual void deleteChildren() = 0;

  virtual bool isDummy() const { return false; };

 protected:
  int locnx, locmy, locmyl, locmz, locmzl;
  sData* data;
  
  int iAmDistrib;

  /* members for blockwise schur complement computation */
  bool computeBlockwiseSC{false};

  int blocksizemax;
  double* colsBlockDense{};
  int* colId{};
  int* colSparsity{};

  /* is this linsys the overall root */
  const bool is_hierarchy_root{false};

  int n_solvers{-1};
  int n_threads_solvers{-1};
  std::vector<std::unique_ptr<DoubleLinearSolver>> solvers_blocked{1};
  std::vector<std::unique_ptr<SymMatrix>> problems_blocked{1};

  SymMatrix* kkt{};
  DoubleLinearSolver* solver{};

 public:
  MPI_Comm mpiComm;
  sTree* stochNode;

  template<typename T>
  struct RFGAC_BLOCK
  {
     /* represents a block like
      * [ R_i F_i^T G_i^T ]             [ R_i^T A_i^T C_i^T ]
      * [ A_i   0     0   ] or possibly [  F_i    0     0   ]
      * [ C_i   0     0   ]             [  G_i    0     0   ]
      */
     T& R;
     T& A;
     T& C;
     T& F;
     T& G;

     RFGAC_BLOCK( T& R, T& A, T& C, T& F, T& G ) :
        R(R), A(A), C(C), F(F), G(G) {};
  };

  typedef RFGAC_BLOCK<StringGenMatrix> BorderLinsys;
  typedef RFGAC_BLOCK<SparseGenMatrix> BorderBiBlock;

  virtual void addLnizi(sData *prob, OoqpVector& z0, OoqpVector& zi);

  virtual void addLniziLinkCons(sData *prob, OoqpVector& z0, OoqpVector& zi, int parentmy, int parentmz);

  virtual void addLniZiHierarchyBorder( DenseGenMatrix& result, BorderLinsys& border );

  /* adds mat to res starting at row_0 col_0 */
  void addMatAt( DenseGenMatrix& res, const SparseGenMatrix& mat, int row_0, int col_0 ) const;

  /* add Bi_{outer}^T to res */
  virtual void addBiTBorder( DenseGenMatrix& res, const BorderBiBlock& BiT) const;

  /* compute Bi_{outer}^T X_i = Bi_{outer}^T Ki^-1 (Bi_{outer} - Bi_{inner} X0) and add it to SC */
  virtual void LniTransMultHierarchyBorder( DenseSymMatrix& SC, const DenseGenMatrix& X0, BorderLinsys& border,
        int parent_nx, int parent_my, int parent_mz );

  /** y += alpha * Lni^T * x */
  void LniTransMult(sData *prob, 
		    SimpleVector& y, 
		    double alpha, SimpleVector& x);

  /** Methods that use dense matrices U and V to compute the 
   *  terms from the Schur complement.
   */
  virtual void allocU(DenseGenMatrix ** Ut, int np);
  virtual void allocV (DenseGenMatrix ** V, int np);
  virtual void computeU_V(sData *prob, DenseGenMatrix* U, DenseGenMatrix* V);

  /** Method(s) that use a memory-friendly mechanism for computing
   *  the terms from the Schur Complement
   */
  virtual void addTermToDenseSchurCompl(sData *prob, DenseSymMatrix& SC);

  virtual void addTermToSchurComplBlocked(sData* /*prob*/, bool /*sparseSC*/, SymMatrix& /*SC*/) { assert( 0 && "not implemented here" ); };
 protected:
//  virtual void addBiTLeftKiBiRightToResBlocked( bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
//        /* const */ BorderBiBlock &border_right, DoubleMatrix& result);

 public:
  /* add you part of the border times rhs to b0 */
  virtual void addBorderTimesRhsToB0( StochVector& /*rhs*/, SimpleVector& /*b0*/, BorderLinsys& /*border*/ )
  { assert( false && "not implemented here" ); };

  /* add you part of the border times rhs to b0 */
  virtual void addBorderX0ToRhs( StochVector& /*rhs*/, const SimpleVector& /*x0*/, BorderLinsys& /*border*/ )
  { assert( false && "not implemented here" ); };

  virtual void addBiTLeftKiBiRightToResBlockedParallelSolvers( bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
        /* const */ BorderBiBlock& border_right, DoubleMatrix& result);

  void addBiTLeftKiDenseToResBlockedParallelSolvers( bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
        /* const */ DenseGenMatrix& BT, DoubleMatrix& result);

  virtual void addTermToSparseSchurCompl(sData* /*prob*/, SparseSymMatrix& /*SC*/ ) { assert(0 && "not implemented here"); };
					
  /** Used in the iterative refinement for the dense Schur complement systems
   * Computes res += [0 A^T C^T ]*inv(KKT)*[0;A;C] x
   */
  virtual void addTermToSchurResidual(sData* prob, 
				      SimpleVector& res, 
				      SimpleVector& x);

  /* solve for all border rhs -> calculate SC = Bborder^T K^-1 Bborder */
  virtual void addInnerToHierarchicalSchurComplement( DenseSymMatrix& /*schur_comp*/, sData* /*prob*/ )
  { assert( false && "not implemented here"); }

  /* compute B_{inner}^T K^{-1} B_{outer} and add it up in result */
  virtual void LsolveHierarchyBorder( DenseGenMatrix& /*result*/, BorderLinsys& /*border*/ )
  { assert( false && "not implemented here" ); };

  /* solve with SC and comput X_0 = SC^-1 B_0 */
  virtual void DsolveHierarchyBorder( DenseGenMatrix& /*buffer_b0*/ )
  { assert( false && "not implemented here" ); };

  /* compute SUM_i Bi_{outer}^T X_i = Bi_{outer}^T Ki^-1 (Bi_{outer} - Bi_{inner} X0) */
  virtual void LtsolveHierarchyBorder( DenseSymMatrix& /*SC*/, const DenseGenMatrix& /*X0*/, BorderLinsys& /*border_outer*/ )
  { assert( false && "not implemented here" ); };

 protected:
  void addLeftBorderTimesDenseColsToResTransp( const BorderBiBlock& border_left, const double* cols,
        const int* cols_id, int length_col, int n_cols, bool sparse_res, bool sym_res, DoubleMatrix& res) const;

  void addLeftBorderTimesDenseColsToResTranspSparse( const BorderBiBlock& border_left, const double* cols,
        const int* cols_id, int length_col, int n_cols, SparseSymMatrix& res) const;

  void addLeftBorderTimesDenseColsToResTranspDense( const BorderBiBlock& border_left, const double* cols,
        const int* cols_id, int length_col, int n_cols, int n_cols_res, double** res) const;

  /* calculate res += X_i * B_i^T */
  void multRightDenseSchurComplBlocked( /* const */ sData* prob, const DenseGenMatrix& X, DenseGenMatrix& result, int parent_nx, int parent_my, int parent_mz );
};

#endif
