/*
 * sLinsysLeafMumps.h
 *
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPSTOCH_SLINSYSLEAFMUMPS_H_
#define PIPS_IPM_CORE_QPSTOCH_SLINSYSLEAFMUMPS_H_

#include "sLinsysLeaf.h"
#include "pipsport.h"

class sFactory;
class DistributedQP;

/** This class solves the linear system corresponding to a leaf node.
 */
class sLinsysLeafMumps : public sLinsysLeaf
{
 public:
  static constexpr int bufferMaxSize = (1024*1024*64);

  sLinsysLeafMumps(sFactory* factory,
        DistributedQP* prob_,
        OoqpVector* dd_, OoqpVector* dq_,
        OoqpVector* nomegaInv_,
        OoqpVector* rhs_)
        : sLinsysLeaf(factory, prob_, dd_, dq_, nomegaInv_, rhs_) {};

  ~sLinsysLeafMumps();

  void addTermToSparseSchurCompl(DistributedQP *prob,
            SparseSymMatrix& SC) override;

  void addTermToDenseSchurCompl(DistributedQP *prob,
            DenseSymMatrix& SC) override;

 private:
  void addTermToSchurComplMumps(DistributedQP *prob, bool sparseSC,
            SymMatrix& SC);

  /* build right matrix for Schur complement; Fortran indexed, CSC, and without empty columns */
  void buildSchurRightMatrix(DistributedQP *prob, SymMatrix& SC);

  SparseGenMatrix* schurRightMatrix_csc{};
  int* schurRightNzColId{};
  int nSC{-1};
  int mSchurRight{-1};
  int nSchurRight{-1};

  double* buffer{};
  int bufferSize{-1};
};

#endif /* PIPS_IPM_CORE_QPSTOCH_SLINSYSLEAFMUMPS_H_ */
