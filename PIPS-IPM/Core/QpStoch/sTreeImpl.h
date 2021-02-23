/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef STOCH_TREE_IMPL
#define STOCH_TREE_IMPL

#include "sTree.h"
#include "stochasticInput.hpp"

/** A full implementation of sTree that is currently used, the other implementation
 *  being sTreeCallbacks (obsolete)
 * 
 */

class sTreeImpl : public sTree
{
 public:
  sTreeImpl(stochasticInput &in, MPI_Comm comm=MPI_COMM_WORLD);
 private: 
  sTreeImpl(int idx, stochasticInput &in);
 public:
  virtual ~sTreeImpl();

  StochSymMatrix*   createQ() const override;
  StochVector*      createc() const override;

  StochVector*      createxlow()  const override;
  StochVector*      createixlow() const override;
  StochVector*      createxupp()  const override;
  StochVector*      createixupp() const override;


  StochGenMatrix*   createA() const override;
  StochVector*      createb() const override;


  StochGenMatrix*   createC() const override;
  StochVector*      createclow()  const override;
  StochVector*      createiclow() const override;
  StochVector*      createcupp()  const override;
  StochVector*      createicupp() const override;

  int nx() const override;
  int my() const override; 
  int mz() const override; 
  int id() const override; 

  void computeGlobalSizes() override;
  void loadLocalSizes() override;
 private:
  int m_id;
  stochasticInput& in;

  sTreeImpl* parent;

  size_t m_nx, m_my, m_mz;

  //int compute_nFirstStageEq();
  void splitConstraints_stage1();
  void splitConstraints_stage2(int scen);
  //int compute_nSecondStageEq(int scen);
};



#endif
