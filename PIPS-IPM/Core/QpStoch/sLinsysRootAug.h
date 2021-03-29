/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYS
#define SAUGLINSYS

#ifdef WITH_PARDISO
#include "PardisoProjectIndefSolver.h"
#endif

#ifdef WITH_MKL_PARDISO
#include "PardisoMKLIndefSolver.h"
#endif

#ifdef WITH_MA57
#include "Ma57SolverRoot.h"
#endif

#ifdef WITH_MA27
#include "Ma27SolverRoot.h"
#endif

#ifdef WITH_MUMPS
#include "MumpsSolverRoot.h"
#endif

#include "sLinsysRoot.h"
#include "pipsport.h"
#include <memory>

class sData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAug : public sLinsysRoot {

 public:
  sLinsysRootAug(sFactory * factory_, sData * prob_);
  sLinsysRootAug(sFactory* factory,
		 sData* prob_,
		 OoqpVector* dd_, OoqpVector* dq_,
		 OoqpVector* nomegaInv_,
		 OoqpVector* rhs_, bool creat_solvers);
  ~sLinsysRootAug() override = default;

  void finalizeKKT( sData* prob, Variables* vars) override;
  void finalizeKKTdist(sData* prob) override;

  void Lsolve(sData *prob, OoqpVector& x) override;
  void Dsolve(sData *prob, OoqpVector& x) override;
  void Ltsolve(sData *prob, OoqpVector& x) override;

  using sLinsys::LsolveHierarchyBorder;
  void LsolveHierarchyBorder( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool two_link_border ) override;

  using sLinsys::LtsolveHierarchyBorder;
  void LtsolveHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br,
        std::vector<BorderMod>& br_mod_border, bool sym_res, bool sparse_res, bool two_link_border ) override;

 protected:
  SymMatrix* createKKT (sData* prob) const;
  void createSolversSparse( SolverType solver );
  void createSolversDense();

  void assembleLocalKKT( sData* prob ) override;
//  void solveReduced( sData *prob, SimpleVector& b);
  void solveReducedLinkCons( sData *prob, SimpleVector& b);
  void solveReducedLinkConsBlocked( sData* data, DenseGenMatrix& rhs_mat_transp, int rhs_start, int n_rhs );
  void addBTKiInvBToSC( DoubleMatrix& result, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
        bool sym_res, bool sparse_res ) override;

 private:
  void createSolversAndKKts(sData* prob);
  void finalizeKKTdense( sData* prob, Variables* vars);
  void finalizeKKTsparse( sData* prob, Variables* vars);
  void solveWithIterRef( sData *prob, SimpleVector& b);
  void solveWithBiCGStab( sData *prob, SimpleVector& b);

  void DsolveHierarchyBorder( DenseGenMatrix& b ) override;

  // add specified columns of given matrix Ht (either Ft or Gt) to Schur complement
  void addLinkConsBlock0Matrix( sData *prob, SparseGenMatrix& Ht, int nHtOffsetCols, int nKktOffsetCols, int startCol, int endCol);

  /** y = beta*y - alpha* SC * x */
  void SCmult ( double beta, SimpleVector& y, double alpha, SimpleVector& x, sData* prob);

  SymMatrix* CtDC{};

  std::vector<double> reduced_rhss_blocked;
  std::unique_ptr<SimpleVector> redRhs;
};

#endif

