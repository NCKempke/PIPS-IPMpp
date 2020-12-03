/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#ifndef PIPSIPM_INTERFACE
#define PIPSIPM_INTERFACE

#include <algorithm>
#include <functional>

#include "stochasticInput.hpp"
#include "sTreeImpl.h"

#include "sTree.h"
#include "sData.h"
#include "sResiduals.h"
#include "sVars.h"
#include "StochMonitor.h"
#include <cstdlib>
#include <stdexcept>
#include <algorithm>


#include "PreprocessFactory.h"
#include "Scaler.h"
#include "Presolver.h"
#include "Postsolver.h"

#include "sTreeCallbacks.h"
#include "pipsport.h"

#include "StochOptions.h"

//#define PRESOLVE_POSTSOLVE_ONLY // will not call solve routine an just presolve and then postsolve the problem - for debugging presolve and postsolve operations

template<class FORMULATION, class IPMSOLVER> 
class PIPSIpmInterface 
{
 public:
  PIPSIpmInterface(stochasticInput &in, MPI_Comm = MPI_COMM_WORLD);
  PIPSIpmInterface(StochInputTree* in, MPI_Comm = MPI_COMM_WORLD,
        ScalerType scaler_type = SCALER_NONE, PresolverType presolver_type = PRESOLVER_NONE, std::string settings = "PIPSIPMpp.opt");
  ~PIPSIpmInterface();

  void go();
  double getObjective();
  double getFirstStageObjective() const;

  void setPrimalTolerance(double val);
  void setDualTolerance(double val);

  std::vector<double> gatherPrimalSolution();
  std::vector<double> gatherDualSolutionEq();
  std::vector<double> gatherDualSolutionIneq();
  std::vector<double> gatherDualSolutionIneqUpp();
  std::vector<double> gatherDualSolutionIneqLow();
  std::vector<double> gatherDualSolutionVarBounds();
  std::vector<double> gatherDualSolutionVarBoundsUpp();
  std::vector<double> gatherDualSolutionVarBoundsLow();

  std::vector<double> getFirstStagePrimalColSolution() const;
  std::vector<double> getSecondStagePrimalColSolution(int scen) const;
  std::vector<double> getFirstStageDualRowSolution() const;
  std::vector<double> getSecondStageDualRowSolution(int scen) const;

  void postsolveComputedSolution();

private:
  void printComplementarityResiduals(const sVars& vars) const;
  std::vector<double> gatherFromSolution( OoqpVectorHandle sVars::* member_to_gather );

public:
  std::vector<double> gatherEqualityConsValues();
  std::vector<double> gatherInequalityConsValues();

  void getVarsUnscaledUnperm();
  void getResidsUnscaledUnperm();
  //more get methods to follow here

  static bool isDistributed() { return true; }

 protected:
  std::unique_ptr<FORMULATION> factory{};
  std::unique_ptr<PreprocessFactory> prefactory{};

  std::unique_ptr<sData> data{};       // possibly presolved data
  std::unique_ptr<sData> dataUnpermNotHier{}; // data after presolve before permutation, scaling and hierarchical data
  std::unique_ptr<sData> origData{};   // original data
  std::unique_ptr<sVars> vars{};
  std::unique_ptr<sVars> unscaleUnpermNotHierVars{};
  std::unique_ptr<sVars> postsolvedVars{};

  std::unique_ptr<sResiduals> resids{};
  std::unique_ptr<sResiduals> unscaleUnpermNotHierResids{};
  std::unique_ptr<sResiduals> postsolvedResids{};

  std::unique_ptr<Presolver> presolver{};
  std::unique_ptr<Postsolver> postsolver{};
  std::unique_ptr<Scaler> scaler{};
  std::unique_ptr<IPMSOLVER> solver{};

  MPI_Comm comm = MPI_COMM_NULL;
  const int my_rank = -1;
  bool ran_solver = false;
};

//----------------------------------------------------------------------
// IMPLEMENTATION
//----------------------------------------------------------------------


template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::PIPSIpmInterface(stochasticInput &in, MPI_Comm comm) :comm(comm), my_rank( PIPS_MPIgetRank(comm) )
{
  factory.reset( new FORMULATION( in, comm) );
#ifdef TIMING
  if( my_rank == 0 ) printf("factory created\n");
  if( my_rank == 0 ) printf("prefactory created\n");
#endif

  data.reset( dynamic_cast<sData*>( factory->makeData() ) );
#ifdef TIMING
  if( my_rank == 0 ) printf("data created\n");
#endif

  vars.reset( dynamic_cast<sVars*>( factory->makeVariables( data.get() ) ) );
#ifdef TIMING
  if( my_rank == 0 ) printf("variables created\n");
#endif

  resids.reset( dynamic_cast<sResiduals*>( factory->makeResiduals( data.get() ) ) );
#ifdef TIMING
  if( my_rank == 0 ) printf("resids created\n");
#endif

  solver.reset( new IPMSOLVER( factory.get(), data.get(), scaler.get() ) );
  solver->addMonitor(new StochMonitor( factory.get() ));
#ifdef TIMING
  if( my_rank == 0 ) printf("solver created\n");
  //solver->monitorSelf();
#endif

}

template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::PIPSIpmInterface(StochInputTree* in, MPI_Comm comm, ScalerType scaler_type,
      PresolverType presolver_type, std::string settings) : comm(comm), my_rank( PIPS_MPIgetRank() )
{
  pips_options::setOptions(settings);
  const bool postsolve = pips_options::getBoolParameter("POSTSOLVE");

  MPI_Barrier(comm);
  const double t0 = MPI_Wtime();

  factory.reset( new FORMULATION( in, comm) );
#ifdef TIMING
  if( my_rank == 0 ) printf("factory created\n");
#endif

  prefactory.reset( new PreprocessFactory() );
#ifdef TIMING
  if( my_rank == 0 ) printf("prefactory created\n");
#endif

  // presolving activated?
  if( presolver_type != PRESOLVER_NONE )
  {
     origData.reset( dynamic_cast<sData*>(factory->makeData()) );

     MPI_Barrier(comm);
     const double t0_presolve = MPI_Wtime();

     if( postsolve )
        postsolver.reset( prefactory->makePostsolver( origData.get() ) );

     presolver.reset( prefactory->makePresolver(dynamic_cast<sFactory*>(factory.get())->tree, origData.get(), presolver_type, postsolver.get()) );

     data.reset( dynamic_cast<sData*>(presolver->presolve()) );

     factory->data = data.get(); // todo update also sTree* of factory

     MPI_Barrier(comm);
     const double t_presolve = MPI_Wtime();
     if( my_rank == 0 )
        std::cout << "---presolve time (in sec.): " << t_presolve - t0_presolve << "\n";
  }
  else
  {
     data.reset( dynamic_cast<sData*>(factory->makeData()) );
     assert( data );
  }

//  data->writeToStreamDense(std::cout);

#if 0
  ofstream myfile;
  myfile.open ("PipsToMPS_prslv.mps");
  data->writeMPSformat(myfile);
  myfile.close();
#endif

#ifdef TIMING
  if( my_rank == 0 ) printf("data created\n");
#endif

  dataUnpermNotHier.reset( data->cloneFull() );

  // after identifying the linking structure switch to hierarchical data structure -> will this do anything to the scaler?
  if( pips_options::getBoolParameter("PARDISO_FOR_GLOBAL_SC") )
     data->activateLinkStructureExploitation();

  // TODO : save "old" data somewhere?
#ifdef HIERARCHICAL
  data.reset( dynamic_cast<sData*>(factory->switchToHierarchicalData( data.release()) ) );
//  data->writeToStreamDense(std::cout);
#endif

  vars.reset( dynamic_cast<sVars*>( factory->makeVariables( data.get() ) ) );
#ifdef TIMING
  if( my_rank == 0 ) printf("variables created\n");
#endif

  resids.reset( dynamic_cast<sResiduals*>( factory->makeResiduals( data.get() ) ) );
#ifdef TIMING
  if( my_rank == 0 ) printf("resids created\n");
#endif

  scaler.reset( prefactory->makeScaler(data.get(), scaler_type) );

#ifdef TIMING
  if( my_rank == 0 ) printf("scaler created\n");
#endif

  if( scaler )
  {
     MPI_Barrier(comm);
     const double t0_scaling = MPI_Wtime();

     scaler->scale();

     MPI_Barrier(comm);
     const double t_scaling = MPI_Wtime();
     if( my_rank == 0 )
        std::cout << "---scaling time (in sec.): " << t_scaling - t0_scaling << "\n";
  }

  solver.reset( new IPMSOLVER( factory.get(), data.get(), scaler.get() ) );
  solver->addMonitor(new StochMonitor( factory.get(), scaler.get() ));
#ifdef TIMING
  if( my_rank == 0 ) printf("solver created\n");
  //solver->monitorSelf();
#endif

  MPI_Barrier(comm);
  const double t1 = MPI_Wtime();
  if( my_rank == 0 )
     std::cout << "---reading time (in sec.): " << t1 - t0 << "\n";
}


template<typename FORMULATION, typename IPMSOLVER>
void PIPSIpmInterface<FORMULATION,IPMSOLVER>::go()
{
  if( my_rank == 0 )
     std::cout << "solving ...\n";

  // TODO : use unlifted data....
  if( my_rank == 0 )
  {
#ifndef HIERARCHICAL
     std::cout << "1st stage " << data->getLocalnx() << " variables, " << data->getLocalmy()
	       << " equality constraints, " << data->getLocalmz() << " inequality constraints.\n";

    const int nscens = data->children.size();
    if( nscens )
    {
       std::cout << "2nd stage " << data->children[0]->getLocalnx() << " variables, "
             << data->children[0]->getLocalmy() << " equality constraints, "
             << data->children[0]->getLocalmz() << " inequality constraints.\n";

       std::cout << nscens << " scenarios." << "\n";
       std::cout << "Total " << data->getLocalnx() + nscens * data->children[0]->getLocalnx() << " variables, "
             << data->getLocalmy() + nscens * data->children[0]->getLocalmy()  << " equality constraints, "
             << data->getLocalmz() + nscens * data->children[0]->getLocalmz() << " inequality constraints.\n";
    }
#endif
  }
#ifdef TIMING
  double tmElapsed=MPI_Wtime();
#endif

#if defined(PRESOLVE_POSTSOLVE_ONLY) && !defined(NDEBUG)
  const int result = 0;
#else
  //---------------------------------------------
  const int result = solver->solve(data.get(), vars.get(), resids.get());
  //---------------------------------------------
#endif

  if( result != 0 && my_rank == 0 )
     std::cout << "failed to solve instance, result code: " << result << "\n";
  
  ran_solver = true;

#ifdef TIMING
   if ( 0 != result )
      return;

   tmElapsed = MPI_Wtime()-tmElapsed;

   const double objective = getObjective();

   if( my_rank == 0 ) {
    //std::cout << " " << data->nx << " variables, " << data->my
    // << " equality constraints, " << data->mz << " inequality constraints.\n";

    std::cout << " Iterates: " << solver->iter <<",    Optimal Solution:  "
          << objective << "\n";

    std::cout << "Solve time: " << tmElapsed << " seconds." << endl;

    char *var = getenv("OMP_NUM_THREADS");
    if(var != nullptr) {
      int num_threads;
      sscanf( var, "%d", &num_threads );
      std::cout << "Num threads: " << num_threads << "\n";
    }
  }
#endif

#if !defined(NDEBUG) && defined(PRESOLVE_POSTSOLVE_ONLY)
  postsolveComputedSolution();
#endif
}

template<typename FORMULATION, typename SOLVER>
double PIPSIpmInterface<FORMULATION,SOLVER>::getObjective() {

  if(!ran_solver)
    throw std::logic_error("Must call go() and start solution process before trying to retrieve original solution");

  if( postsolver != nullptr && postsolvedVars == nullptr)
    this->postsolveComputedSolution();

  double obj;
  if(postsolvedVars != nullptr)
    obj = origData->objectiveValue( postsolvedVars.get() );
  else
  {
    obj = data->objectiveValue( vars.get() );
    if( scaler )
       obj = scaler->getObjUnscaled(obj);
  }

  return obj;
}


template<typename FORMULATION, typename SOLVER>
double PIPSIpmInterface<FORMULATION,SOLVER>::getFirstStageObjective() const
{
  OoqpVector& x = *(dynamic_cast<StochVector&>(*vars->x).vec);
  OoqpVector& c = *(dynamic_cast<StochVector&>(*data->g).vec);
  return c.dotProductWith(x);
}

template<class FORMULATION, class IPMSOLVER>
PIPSIpmInterface<FORMULATION, IPMSOLVER>::~PIPSIpmInterface()
{
}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::getVarsUnscaledUnperm()
{
  assert(unscaleUnpermNotHierVars == nullptr);
  assert(dataUnpermNotHier);

  if(!ran_solver)
    throw std::logic_error("Must call go() and start solution process before trying to retrieve unscaled unpermutated solution");
  if( scaler )
  {
    std::unique_ptr<sVars> unscaled_vars{ dynamic_cast<sVars*>(scaler->getVariablesUnscaled(*vars)) };
    unscaleUnpermNotHierVars.reset( data->getVarsUnperm(*unscaled_vars, *dataUnpermNotHier) );
  }
  else
    unscaleUnpermNotHierVars.reset( data->getVarsUnperm(*vars, *dataUnpermNotHier) );

}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::getResidsUnscaledUnperm()
{
  assert(unscaleUnpermNotHierResids == nullptr);
  assert(dataUnpermNotHier);

  if(!ran_solver)
    throw std::logic_error("Must call go() and start solution process before trying to retrieve unscaled unpermutated residuals");
  if( scaler )
  {
    std::unique_ptr<sResiduals> unscaled_resids{ dynamic_cast<sResiduals*>(scaler->getResidualsUnscaled(*resids)) };
    unscaleUnpermNotHierResids.reset( data->getResidsUnperm(*unscaled_resids, *dataUnpermNotHier) );
  }
  else
    unscaleUnpermNotHierResids.reset( data->getResidsUnperm(*resids, *dataUnpermNotHier) );
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherFromSolution( OoqpVectorHandle sVars::* member_to_gather )
{
  if( unscaleUnpermNotHierVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr)
    this->postsolveComputedSolution();

  std::vector<double> vec;
  if( postsolver == nullptr)
    vec = dynamic_cast<const StochVector&>(*(unscaleUnpermNotHierVars.get()->*member_to_gather)).gatherStochVector();
  else
    vec = dynamic_cast<const StochVector&>(*(postsolvedVars.get()->*member_to_gather)).gatherStochVector();

  return vec;

}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherPrimalSolution()
{
  return gatherFromSolution( &sVars::x );
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionEq()
{
  return gatherFromSolution( &sVars::y );
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneq()
{
  return gatherFromSolution( &sVars::z );
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneqUpp()
{
  return gatherFromSolution( &sVars::pi );
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionIneqLow()
{
  return gatherFromSolution( &sVars::lambda );
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBounds()
{
  std::vector<double> duals_varbounds_upp = gatherDualSolutionVarBoundsUpp();
  std::vector<double> duals_varbounds_low = gatherDualSolutionVarBoundsLow();

  assert(duals_varbounds_low.size() == duals_varbounds_upp.size());

  std::vector<double> duals_varbounds;
  duals_varbounds.reserve(duals_varbounds_low.size());

  std::transform(duals_varbounds_low.begin(), duals_varbounds_low.end(), duals_varbounds_upp.begin(), std::back_inserter(duals_varbounds), std::minus<double>());

  return duals_varbounds;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBoundsUpp()
{
  return gatherFromSolution( &sVars::phi );
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherDualSolutionVarBoundsLow()
{
  return gatherFromSolution( &sVars::gamma );
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherEqualityConsValues()
{
  if( unscaleUnpermNotHierResids == nullptr)
    this->getResidsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr )
    this->postsolveComputedSolution();

  StochVector* eq_vals = (postsolvedVars == nullptr) ? dynamic_cast<StochVector*>(unscaleUnpermNotHierResids->rA->cloneFull()) :
    dynamic_cast<StochVector*>(postsolvedResids->rA->cloneFull());

  if( origData == nullptr || postsolvedVars == nullptr )
    eq_vals->axpy(1.0, *data->bA);
  else
    eq_vals->axpy(1.0, *origData->bA);

  std::vector<double> eq_vals_vec = eq_vals->gatherStochVector();

  delete eq_vals;

  return eq_vals_vec;
}


template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::gatherInequalityConsValues()
{
  if( unscaleUnpermNotHierVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( unscaleUnpermNotHierResids == nullptr)
    this->getResidsUnscaledUnperm();

  if( postsolver != nullptr && postsolvedVars == nullptr )
    this->postsolveComputedSolution();

  StochVector* ineq_vals = (postsolvedVars == nullptr) ? dynamic_cast<StochVector*>(unscaleUnpermNotHierResids->rC->cloneFull()) :
    dynamic_cast<StochVector*>(postsolvedResids->rC->cloneFull());

  if( postsolvedVars == nullptr )
    ineq_vals->axpy(1.0, *unscaleUnpermNotHierVars->s);
  else
    ineq_vals->axpy(1.0, *postsolvedVars->s);

  std::vector<double> ineq_vals_vec = ineq_vals->gatherStochVector();

  delete ineq_vals;

  return ineq_vals_vec;
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getFirstStagePrimalColSolution() const
{
	SimpleVector const &v = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->x).vec);
	return std::vector<double>(&v[0],&v[0]+v.length());
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getSecondStagePrimalColSolution(int scen) const {
	SimpleVector const &v = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->x).children[scen]->vec);
	if(!v.length())
	  return std::vector<double>(); //this vector is not on this processor
	else
	  return std::vector<double>( &v[0], &v[0] + v.length() );
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getFirstStageDualRowSolution() const
{
  assert( false && "not in use - can only be called when using sTreeImpl - we use sTreeCallbacks ...");
  SimpleVector const &y =
        *dynamic_cast<SimpleVector const*>((dynamic_cast<StochVector const&>(*vars->y)).vec);
  SimpleVector const &z =
        *dynamic_cast<SimpleVector const*>((dynamic_cast<StochVector const&>(*vars->z)).vec);

  if( !y.length() && !z.length() )
     return std::vector<double>(); //this vector is not on this processor
  else
  {
     std::vector<int> const &map = dynamic_cast<sTreeImpl*>(factory->tree)->idx_EqIneq_Map;

     std::vector<double> multipliers(map.size());
     for( size_t i = 0; i < map.size(); i++ )
     {
        int idx = map[i];
        if( idx < 0 )
        {
           //equality
           idx = -idx - 1;
           assert(idx >= 0);
           multipliers[i] = y[idx];
        }
        else
        {
           //inequality - since, we have z-\lambda+\pi=0, where \lambda is the multiplier for low and
           //\pi is the multiplier for upp, therefore z containts the right multiplier for this row.
#ifndef NDEBUG
           SimpleVector const &iclow = *dynamic_cast<SimpleVector const*>((dynamic_cast<StochVector const&>(*vars->iclow)).vec);
           SimpleVector const &icupp = *dynamic_cast<SimpleVector const*>((dynamic_cast<StochVector const&>(*vars->icupp)).vec);
           assert(iclow[idx] > 0 || icupp[idx] > 0);
#endif
           multipliers[i] = z[idx];
        }
     }
     return multipliers;
  }
}

template<class FORMULATION, class IPMSOLVER>
std::vector<double> PIPSIpmInterface<FORMULATION, IPMSOLVER>::getSecondStageDualRowSolution(int scen) const 
{
  assert( false && "not in use - can only be called when using sTreeImpl - we use sTreeCallbacks ...");
  SimpleVector const &y = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->y).children[scen]->vec);
  SimpleVector const &z = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->z).children[scen]->vec);
  SimpleVector const &iclow = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->iclow).children[scen]->vec);
  SimpleVector const &icupp = *dynamic_cast<SimpleVector const*>(dynamic_cast<StochVector const&>(*vars->icupp).children[scen]->vec);
  //assert(v.length());
  if( !y.length() && !z.length() )
    return std::vector<double>(); //this vector is not on this processor
  else 
  {
    std::vector<int> const &map= dynamic_cast<sTreeImpl*>(factory->tree->children[scen])->idx_EqIneq_Map;

    std::vector<double> multipliers(map.size());
    for(size_t i = 0; i < map.size(); i++)
    {
      int idx = map[i];
      if(idx<0) 
      {
	      //equality
	      idx =- idx - 1;
        assert( idx >= 0 );
        multipliers[i] = y[idx];
      } 
      else
      {
	      //inequality - since, we have z-\lambda+\pi=0, where \lambda is the multiplier for low and
	      //\pi is the multiplier for upp, therefore z containts the right multiplier for this row.
#ifndef NDEBUG
	      assert(iclow[idx] > 0 || icupp[idx] > 0);
	      multipliers[i] = z[idx];
#endif
      }
    }
    return multipliers;
  }
  //return std::vector<double>(&v[0],&v[0]+v.length());
}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::printComplementarityResiduals(const sVars& svars) const
{
  const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);

  /* complementarity residuals before postsolve */
  std::unique_ptr<OoqpVectorBase<double>> t_clone{ svars.t->cloneFull() };
  std::unique_ptr<OoqpVectorBase<double>> u_clone{ svars.u->cloneFull() };
  std::unique_ptr<OoqpVectorBase<double>> v_clone{ svars.v->cloneFull() };
  std::unique_ptr<OoqpVectorBase<double>> w_clone{ svars.w->cloneFull() };

  t_clone->componentMult(*svars.lambda);
  t_clone->selectNonZeros(*svars.iclow);

  u_clone->componentMult(*svars.pi);
  u_clone->selectNonZeros(*svars.icupp);

  v_clone->componentMult(*svars.gamma);
  v_clone->selectNonZeros(*svars.ixlow);

  w_clone->componentMult(*svars.phi);
  w_clone->selectNonZeros(*svars.ixupp);

  const double rlambda_infnorm = t_clone->infnorm();
  const double rpi_infnorm = u_clone->infnorm();
  const double rgamma_infnorm = v_clone->infnorm();
  const double rphi_infnorm = w_clone->infnorm();

  if( my_rank == 0 )
  {
     std::cout << " rl norm = " << rlambda_infnorm << "\n";
     std::cout << " rp norm = " << rpi_infnorm << "\n";
     std::cout << " rg norm = " << rgamma_infnorm << "\n";
     std::cout << " rf norm = " << rphi_infnorm << "\n";
     std::cout << "\n";
  }
}

template<class FORMULATION, class IPMSOLVER>
void PIPSIpmInterface<FORMULATION, IPMSOLVER>::postsolveComputedSolution()
{
  const bool print_residuals = pips_options::getBoolParameter("POSTSOLVE_PRINT_RESIDS");
  const int my_rank = PIPS_MPIgetRank(comm);

  assert(origData);
  assert(data);

#if !defined(NDEBUG) && defined(PRESOLVE_POSTSOLVE_ONLY) // todo : resids for C also need recomputation.. - s variable
  /* todo: randomize all vectors x since it has not actually been set to anything */
  vars->x->setToConstant(0.1);
  resids->calcresids(data, vars);
#endif

  if( unscaleUnpermNotHierVars == nullptr)
    this->getVarsUnscaledUnperm();

  if( unscaleUnpermNotHierResids == nullptr)
    this->getResidsUnscaledUnperm();

  if( postsolvedVars != nullptr || postsolvedResids != nullptr )
    return;

  if( postsolver == nullptr )
  {
    assert( "no postsolver available" && 0 );
    return;
  }

  if( print_residuals )
  {
     if( my_rank == 0 )
        std::cout << "\n" << "Residuals before postsolve:" << "\n";
     resids->calcresids(data.get(), vars.get(), print_residuals);
     printComplementarityResiduals(*vars);

     if( my_rank == 0 )
        std::cout << "Residuals after unscaling/permuting:" << "\n";
     unscaleUnpermNotHierResids->calcresids(dataUnpermNotHier.get(), unscaleUnpermNotHierVars.get(), print_residuals);
     printComplementarityResiduals(*unscaleUnpermNotHierVars);
  }

  MPI_Barrier(comm);
  const double t0_postsolve = MPI_Wtime();

#ifdef HIERARCHICAL
  factory->collapseHierarchicalTree();
#endif
  dynamic_cast<sTreeCallbacks*>(factory->tree)->switchToOriginalData();
  factory->data = origData.get();

  postsolvedVars.reset( dynamic_cast<sVars*>( factory->makeVariables( origData.get() ) ) );

  postsolvedResids.reset( dynamic_cast<sResiduals*>( factory->makeResiduals( origData.get() ) ) );
  postsolver->postsolve(*unscaleUnpermNotHierVars, *postsolvedVars);

  double obj_postsolved = origData->objectiveValue(postsolvedVars.get());

  MPI_Barrier(comm);
  const double t_postsolve = MPI_Wtime();

  if( my_rank == 0 )
  {
     std::cout << "---postsolve time (in sec.): " << t_postsolve - t0_postsolve << "\n";
     std::cout << "Objective value after postsolve: " << obj_postsolved << "\n";
  }

  /* compute residuals for postprocessed solution and check for feasibility */
  if( print_residuals )
  {
     if( my_rank == 0 )
        std::cout << "\n" << "Residuals after postsolve:" << "\n";
     postsolvedResids->calcresids(origData.get(), postsolvedVars.get(), print_residuals);

     printComplementarityResiduals(*postsolvedVars);
  }
}

#endif
