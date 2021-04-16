/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHLEAFLINSYS
#define STOCHLEAFLINSYS

#include "sLinsys.h"
#include "sTree.h"
#include "sFactory.h"
#include "DistributedQP.hpp"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"

#include "pipsport.h"

#include "omp.h"

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to QpGenSparseLinsys.
 */
class sLinsysLeaf : public sLinsys
{
 public:
    sLinsysLeaf(sFactory* factory,
		DistributedQP* prob_,				    
		OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
		OoqpVector* rhs_);

  ~sLinsysLeaf() override = default;

  void factor2( DistributedQP *prob, Variables *vars) override;
  void assembleKKT(DistributedQP*, Variables*) override {};
  void allreduceAndFactorKKT(DistributedQP* prob, Variables* vars) override { factor2( prob, vars ); };

  void Lsolve( DistributedQP*, OoqpVector& ) override {};
  void Dsolve( DistributedQP*, OoqpVector& x ) override;
  void Ltsolve( DistributedQP*, OoqpVector& ) override {};

  //void Lsolve2 ( OoqpVector& x ) override;
  //void Dsolve2 ( OoqpVector& x ) override;
  void Ltsolve2( DistributedQP *prob, StochVector& x, SimpleVector& xp, bool) override;

  void putZDiagonal( OoqpVector& zdiag ) override;
  //void solveCompressed( OoqpVector& rhs ) override;
  void putXDiagonal( OoqpVector& xdiag_ ) override;

  //void Ltsolve_internal(  DistributedQP *prob, StochVector& x, SimpleVector& xp);
  void deleteChildren() override;

  void addTermToSchurComplBlocked( sData *prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC, int ) override;

  void addLniziLinkCons(DistributedQP *prob, OoqpVector& z0_, OoqpVector& zi_, bool ) override;

  void addInnerBorderKiInvBrToRes( DoubleMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool, bool sparse_res, bool sym_res, int begin_cols, int end_cols, int ) override;
  void LniTransMultHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
        bool sparse_res, bool sym_res, bool, int begin_cols, int end_cols, int n_empty_rows_inner_border ) override;

 protected:

  static void mySymAtPutSubmatrix(SymMatrix& kkt, GenMatrix& B, GenMatrix&, int locnx, int locmy, int);

  void addBorderTimesRhsToB0( StochVector& rhs, SimpleVector& b0, BorderLinsys& border ) override;
  void addBorderX0ToRhs( StochVector& rhs, const SimpleVector& x0, BorderLinsys& border ) override;
 private:
  void addBorderTimesRhsToB0( SimpleVector& rhs, SimpleVector& b0, BorderBiBlock& border );
  void addBorderX0ToRhs( SimpleVector& rhs, const SimpleVector& x0, BorderBiBlock& border );

protected:

  /* compute result += B_inner^T K^-1 Br */
  void addInnerBorderKiInvBrToRes( DenseGenMatrix& result, BorderLinsys& Br, int begin_cols, int end_cols );

  /* compute result += B_inner^T K^-1 ( Br - Br_mod_border ) */
  void addLeftBorderKiInvBrToRes( DoubleMatrix& result, BorderBiBlock& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sparse_res, bool sym_res, int begin_cols_br, int end_cols_br,
              int begin_cols_res, int end_cols_res );
}; 

#endif
