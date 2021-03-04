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
#include "StringGenMatrix.h"

#include <vector>
#include <memory>
#include <tuple>

#include "mpi.h"

class sTree;
class sFactory;
class sData;
class StringSymMatrix;

class sLinsys : public QpGenLinsys
{
 public:
      template<typename T>
      struct RACFG_BLOCK
      {
         private:
         std::unique_ptr<T> dummy{ new T() };

         public:
         const bool use_local_RAC{};
         const bool has_RAC{};
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

         RACFG_BLOCK( T& R, T& A, T& C, T& F, T& G ) :
            has_RAC{true}, R{R}, A{A}, C{C}, F{F}, G{G} {};

         RACFG_BLOCK( T& F, T& G, bool use_local_RAC ) :
            use_local_RAC{ use_local_RAC }, has_RAC{false}, R{*dummy}, A{*dummy}, C{*dummy}, F{F}, G{G} {};

         RACFG_BLOCK( const RACFG_BLOCK<T>& block ) :
            use_local_RAC{ block.use_local_RAC }, has_RAC{ block.has_RAC }, R{ block.R }, A{ block.A }, C{ block.C },
             F{ block.F }, G{ block.G } {};
      };

      using BorderLinsys = RACFG_BLOCK<StringGenMatrix>;
      using BorderBiBlock = RACFG_BLOCK<SparseGenMatrix>;

      static BorderLinsys getChild( BorderLinsys& border, unsigned int i )
      {
         const bool dummy = border.F.children[i]->isKindOf(kStringGenDummyMatrix);
         assert( i < border.F.children.size() );
         if( border.has_RAC )
         {
            if( !dummy && border.F.children[i]->mat->isKindOf(kStringGenMatrix) )
               return BorderLinsys( dynamic_cast<StringGenMatrix&>(*border.R.children[i]->mat),
                     dynamic_cast<StringGenMatrix&>(*border.A.children[i]->mat),
                     dynamic_cast<StringGenMatrix&>(*border.C.children[i]->mat),
                     dynamic_cast<StringGenMatrix&>(*border.F.children[i]->mat),
                     dynamic_cast<StringGenMatrix&>(*border.G.children[i]->mat)
                  );
            else
               return BorderLinsys( *border.R.children[i], *border.A.children[i], *border.C.children[i],
                  *border.F.children[i], *border.G.children[i] );
         }
         else
         {
            if( !dummy && border.F.children[i]->mat->isKindOf(kStringGenMatrix) )
               return BorderLinsys( dynamic_cast<StringGenMatrix&>(*border.F.children[i]->mat),
                     dynamic_cast<StringGenMatrix&>(*border.G.children[i]->mat),
                     border.use_local_RAC
                  );
            else
               return BorderLinsys( *border.F.children[i], *border.G.children[i], border.use_local_RAC );
         }
      }

      template<typename T>
      struct BorderMod_Block
      {
         public:
            BorderLinsys border;
            const T& multiplier;

            BorderMod_Block( BorderLinsys& border_, const T& multiplier ) :
               border{ border_ }, multiplier{ multiplier } {};
      };

      using BorderMod = BorderMod_Block<DenseGenMatrix>;
      using BorderModVector = BorderMod_Block<StochVector>;

      template<typename T>
      static BorderMod_Block<T> getChild( BorderMod_Block<T>& bordermod, unsigned int i )
      {
         BorderLinsys child = getChild( bordermod.border, i );
         return BorderMod_Block<T>( child, bordermod.multiplier );
      }

  sLinsys(sFactory* factory, sData* prob, bool is_hierarchy_root = false);
  sLinsys(sFactory* factory,
		   sData* prob, 
		   OoqpVector* dd, 
		   OoqpVector* dq,
		   OoqpVector* nomegaInv,
		   OoqpVector* rhs,
		   bool create_iter_ref_vecs
  );

  ~sLinsys() override;

  void factor (Data *prob, Variables *vars) override;

  virtual void factor2(sData *prob, Variables *vars) = 0;
  virtual void assembleKKT(sData *prob, Variables *vars) = 0;
  virtual void allreduceAndFactorKKT(sData *prob, Variables *vars) = 0;

  virtual void Lsolve( sData *prob, OoqpVector& x ) = 0;
  virtual void Dsolve( sData *prob, OoqpVector& x ) = 0;
  virtual void Ltsolve( sData *prob, OoqpVector& x ) = 0;
  virtual void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp, bool use_local_RAC) = 0;

  void solveCompressed( OoqpVector& rhs ) override;

  void joinRHS( OoqpVector& rhs_in, const OoqpVector& rhs1_in,
		const OoqpVector& rhs2_in, const OoqpVector& rhs3_in ) const override;

  void separateVars( OoqpVector& x_in, OoqpVector& y_in,
		     OoqpVector& z_in, const OoqpVector& vars_in ) const override;

  virtual void deleteChildren() = 0;

  virtual bool isDummy() const { return false; };

 protected:
  int locnx, locmy, locmyl, locmz, locmzl;
  sData* data{};
  
  int iAmDistrib;

  /* members for blockwise schur complement computation */
  bool computeBlockwiseSC{false};

  int blocksizemax{0};
  double* colsBlockDense{};
  int* colId{};
  int* colSparsity{};

  /* is this linsys the overall root */
  const bool is_hierarchy_root{false};

  int n_solvers{-1};
  int n_threads_solvers{-1};

  std::unique_ptr<SymMatrix> kkt{};
  std::unique_ptr<DoubleLinearSolver> solver{};

 public:
  MPI_Comm mpiComm{MPI_COMM_NULL};
  sTree* stochNode{};

  virtual void addLnizi(sData *prob, OoqpVector& z0, OoqpVector& zi);
  virtual void addLniziLinkCons( sData */*prob*/, OoqpVector& /*z0*/, OoqpVector& /*zi*/, bool /*use_local_RAC*/ ) { assert( false && "not implemented here"); };


  /* adds mat to res starting at row_0 col_0 */
  void addMatAt( DenseGenMatrix& res, const SparseGenMatrix& mat, int row_0, int col_0 ) const;

  /* add BiT to res */
  virtual void addBiTBorder( DenseGenMatrix& res, const BorderBiBlock& BiT) const;

  /* compute Bli^T X_i = Bli^T Ki^-1 (Bri - Bi_{inner} X0) and add it to SC */
  virtual void LniTransMultHierarchyBorder( DoubleMatrix& /*SC*/, const DenseGenMatrix& /*X0*/, BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/,
        std::vector<BorderMod>& /*Br_mod_border*/, bool /*sparse_res*/, bool /*sym_res*/, bool /*use_local_RAC*/) { assert( false && "not implemented here"); };

  /** y += alpha * Lni^T * x */
  virtual void LniTransMult(sData *prob, 
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

  virtual void addTermToSchurComplBlocked(sData* /*prob*/, bool /*sparseSC*/, SymMatrix& /*SC*/, bool /*use_local_RAC*/ ) { assert( 0 && "not implemented here" ); };

  virtual void computeInnerSystemRightHandSide( StochVector& /*rhs_inner*/, const SimpleVector& /*b0*/, bool /*use_local_RAC*/ ) { assert( false && "not implemented here" ); };
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

  /* compute result += Bl^T K^-1 Br where K is our own linear system */
  virtual void addBTKiInvBToSC( DoubleMatrix& /*result*/, BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/,
        bool /*sym_res*/, bool /*sparse_res*/ )
  { assert( false && "not implemented here"); }

  /* compute Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij ) and add it up in result */
  virtual void LsolveHierarchyBorder( DenseGenMatrix& /*result*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/ )
  { assert( false && "not implemented here" ); };

  virtual void LsolveHierarchyBorder( DenseGenMatrix& /*result*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/, bool /*use_local_RAC*/ )
  { assert( false && "not implemented here" ); };

  /* solve with SC and comput X_0 = SC^-1 B_0 */
  virtual void DsolveHierarchyBorder( DenseGenMatrix& /*buffer_b0*/ )
  { assert( false && "not implemented here" ); };

  /* compute RES += SUM_i Bli_^T X_i = Bli^T Ki^-1 ( ( Bri - sum_j Bmodij Xij ) - Bi_{inner} X0) */
  virtual void LtsolveHierarchyBorder( DoubleMatrix& /*res*/, const DenseGenMatrix& /*X0*/,
        BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/, bool /*sym_res*/, bool /*sparse_res*/)
  { assert( false && "not implemented here" ); };

  virtual void LtsolveHierarchyBorder( DoubleMatrix& /*res*/, const DenseGenMatrix& /*X0*/,
        BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/, bool /*sym_res*/, bool /*sparse_res*/, bool /*use_local_RAC*/)
  { assert( false && "not implemented here" ); };

  /* compute Bi_{inner}^T Ki^{-1} ( Bri - sum_j Brmod_ij Xj )and add it to result */
  virtual void addInnerBorderKiInvBrToRes( DenseGenMatrix& /*result*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/, bool /*has_RAC*/ )
  { assert( false && "not implemented here" ); };

 protected:
  void addLeftBorderTimesDenseColsToResTransp( const BorderBiBlock& Bl, const double* cols,
        const int* cols_id, int length_col, int n_cols, bool sparse_res, bool sym_res, DoubleMatrix& res ) const;

  void addLeftBorderTimesDenseColsToResTranspSparse( const BorderBiBlock& Bl, const double* cols,
        const int* cols_id, int length_col, int n_cols, SparseSymMatrix& res ) const;

  void addLeftBorderTimesDenseColsToResTranspDense( const BorderBiBlock& Bl, const double* cols,
        const int* cols_id, int length_col, int n_cols, int n_cols_res, double** res) const;

  /* calculate res -= BT * X */
  void finalizeDenseBorderBlocked( BorderLinsys& B, const DenseGenMatrix& X, DenseGenMatrix& result );

  /* calculate res -= X * BT */
  void multRightDenseBorderBlocked( BorderBiBlock& BT, const DenseGenMatrix& X, DenseGenMatrix& result );

  /* calculate res -= (sum_j XjT * BjT ) */
  void multRightDenseBorderModBlocked( std::vector<BorderMod>& border_mod, DenseGenMatrix& result );

  /* calculate res -= (sum_j X0jT * B0JT ) */
  void finalizeDenseBorderModBlocked( std::vector<BorderMod>& border_mod, DenseGenMatrix& result );

};

#endif
