/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHLEAFLINSYS
#define STOCHLEAFLINSYS

#include "sLinsys.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
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
		sData* prob_,				    
		OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
      OoqpVector* primal_reg_,
      OoqpVector* dual_y_reg_,
      OoqpVector* dual_z_reg_,
		OoqpVector* rhs_);

  ~sLinsysLeaf() override = default;

  void factor2( sData *prob, Variables *vars) override;
  void assembleKKT(sData*, Variables*) override {};
  void allreduceAndFactorKKT(sData* prob, Variables* vars) override { factor2( prob, vars ); };

  void Lsolve( sData*, OoqpVector& ) override {};
  void Dsolve( sData*, OoqpVector& x ) override;
  void Ltsolve( sData*, OoqpVector& ) override {};

  //void Lsolve2 ( OoqpVector& x ) override;
  //void Dsolve2 ( OoqpVector& x ) override;
  void Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp, bool) override;

  void putXDiagonal( const OoqpVector& xdiag_ ) override;
  void putZDiagonal( const OoqpVector& zdiag_ ) override;

  void addRegularization( OoqpVector& regP_, OoqpVector& regDy_, OoqpVector& regDz_ ) const override;
  void addRegularizationsToKKTs( const OoqpVector& regP_, const OoqpVector& regDy_, const OoqpVector& regDz_ ) override;

  //void Ltsolve_internal(  sData *prob, StochVector& x, SimpleVector& xp);
  void deleteChildren() override;

  void addTermToSchurComplBlocked( sData *prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC ) override;

  void addLniziLinkCons(sData *prob, OoqpVector& z0_, OoqpVector& zi_, bool ) override;

  void addInnerBorderKiInvBrToRes( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool ) override;
  void LniTransMultHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
        bool sparse_res, bool sym_res, bool ) override;

 protected:

  static void mySymAtPutSubmatrix(SymMatrix& kkt, GenMatrix& B, GenMatrix&, int locnx, int locmy, int);

  template<class LINSOLVER>
  void initBlockedSolvers();

  void addBorderTimesRhsToB0( StochVector& rhs, SimpleVector& b0, BorderLinsys& border ) override;
  void addBorderX0ToRhs( StochVector& rhs, const SimpleVector& x0, BorderLinsys& border ) override;
 private:
  void addBorderTimesRhsToB0( SimpleVector& rhs, SimpleVector& b0, BorderBiBlock& border );
  void addBorderX0ToRhs( SimpleVector& rhs, const SimpleVector& x0, BorderBiBlock& border );

  /* compute result += B_inner^T K^-1 Br */
  void addInnerBorderKiInvBrToRes( DenseGenMatrix& result, BorderLinsys& Br );

  /* compute result += B_inner^T K^-1 ( Br - Br_mod_border ) */
  void addInnerBorderKiInvBrToResDense( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border );
}; 

#endif
