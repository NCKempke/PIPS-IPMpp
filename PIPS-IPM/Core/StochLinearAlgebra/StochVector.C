#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "VectorUtilities.h"
#include "StochVector.h"
#include "pipsport.h"

#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#include <math.h>
#include "StochVector_fwd.h"

template<typename T>
StochVectorBase<T>::StochVectorBase( SimpleVectorBase<T>* vec, SimpleVectorBase<T>* vecl, MPI_Comm mpi_comm)
   : OoqpVectorBase<T>(0), vec(vec), vecl(vecl), parent(nullptr), mpiComm( mpi_comm ), iAmDistrib( PIPS_MPIgetDistributed(mpiComm) ),
     iAmSpecial( PIPS_MPIiAmSpecial( iAmDistrib, mpiComm) )
{
   assert( vec || vecl );

   if( vec )
      this->n += vec->length();
   if( vecl )
      this->n += vecl->length();
}

template<typename T>
StochVectorBase<T>::StochVectorBase(int n_, MPI_Comm mpiComm_ )
  : OoqpVectorBase<T>(n_), vecl(nullptr), parent(nullptr), mpiComm(mpiComm_), iAmDistrib( PIPS_MPIgetDistributed(mpiComm) ),
    iAmSpecial( PIPS_MPIiAmSpecial( iAmDistrib, mpiComm) )
{
   assert( n_ >= 0 );

   vec = new SimpleVectorBase<T>(n_);
   vecl = nullptr;
}

template<typename T>
StochVectorBase<T>::StochVectorBase(int n_, int nl_, MPI_Comm mpiComm_ )
  : OoqpVectorBase<T>(0), parent(nullptr), mpiComm(mpiComm_), iAmDistrib( PIPS_MPIgetDistributed(mpiComm) ),
    iAmSpecial( PIPS_MPIiAmSpecial( iAmDistrib, mpiComm) )
{
   this->n = 0;

   if( n_ >= 0)
   {
      vec = new SimpleVectorBase<T>(n_);
      this->n += n_;
   }
   else
      vec = nullptr;

   if( nl_ >= 0 )
   {
      vecl = new SimpleVectorBase<T>(nl_);
      this->n += nl_;
   }
   else
      vecl = nullptr;
}

template<typename T>
void StochVectorBase<T>::AddChild(StochVectorBase<T>* child)
{
  child->parent = this;
  children.push_back(child);

  this->n += child->n;
}

template<typename T>
void StochVectorBase<T>::AddChild(OoqpVectorBase<T>* child_)
{
  StochVectorBase<T>* child = reinterpret_cast<StochVectorBase<T>*>(child_);
  AddChild(child);
}

template<typename T>
StochVectorBase<T>::~StochVectorBase()
{
  for (size_t it = 0; it < children.size(); it++)
    delete children[it];

  if( vec )
    delete vec;

  if( vecl )
	 delete vecl;
}

template<typename T>
OoqpVectorBase<T>*
StochVectorBase<T>::clone() const
{
   assert( vec || vecl );
   StochVectorBase<T> *clone;
   clone = new StochVectorBase<T>(vec ? vec->length() : -1, vecl ? vecl->length() : -1, mpiComm);

   for( size_t it = 0; it < children.size(); it++ )
      clone->AddChild(children[it]->clone());

   return clone;
}

template<typename T>
OoqpVectorBase<T>* StochVectorBase<T>::cloneFull() const
{
   assert( vec || vecl );
   StochVectorBase<T>* clone = new StochVectorBase<T>(vec ? vec->length() : -1, vecl ? vecl->length() : -1, mpiComm);

   if( vec )
      clone->vec->copyFrom(*vec);
   if( vecl )
      clone->vecl->copyFrom(*vecl);

   for( size_t it = 0; it < children.size(); it++ )
      clone->AddChild(children[it]->cloneFull());

   return clone;
}

template<typename T>
void StochVectorBase<T>::setNotIndicatedEntriesToVal(T val, const OoqpVectorBase<T>& ind )
{
   const StochVectorBase<T>& ind_vec = dynamic_cast<const StochVectorBase<T>&>(ind);

   assert(this->children.size() == ind_vec.children.size());
   assert( (this->vec && ind_vec.vec) || (this->vec == nullptr && ind_vec.vec == nullptr) );
   assert( (this->vecl && ind_vec.vecl) || (this->vecl == nullptr && ind_vec.vecl == nullptr) );

   if( this->vec )
      this->vec->setNotIndicatedEntriesToVal(val, *ind_vec.vec);

   if( this->vecl )
      this->vecl->setNotIndicatedEntriesToVal(val, *ind_vec.vecl);

   for( size_t node = 0; node < children.size(); ++node )
      this->children[node]->setNotIndicatedEntriesToVal(val, *ind_vec.children[node] );
}

template<typename T>
void StochVectorBase<T>::jointCopyFrom(const StochVectorBase<T>& vx, const StochVectorBase<T>& vy, const StochVectorBase<T>& vz)
{
   assert( this->vec );
   SimpleVectorBase<T>& sv  = dynamic_cast<SimpleVectorBase<T>&>(*this->vec);
   assert( sizeof(T) == sizeof(sv[0]) );

   assert( this->children.size() == vx.children.size() );
   assert( this->children.size() == vy.children.size() );
   assert( this->children.size() == vz.children.size() );
   const int N = sv.length();
   int n1 = 0;
   int n2 = 0;
   int n3 = 0;
   int n4 = 0;
   int n5 = 0;
   int n6 = 0;

   if( vx.vec )
   {
      const SimpleVectorBase<T>& svx = dynamic_cast<const SimpleVectorBase<T>&>(*vx.vec);
      n1 = svx.length();
      assert( n1 >= 0 );

      assert( n1 <= N );
      if( n1 > 0 )
         memcpy(&sv[0], &svx[0], n1 * sizeof(T));
   }

   if( vy.vec )
   {
      const SimpleVectorBase<T>& svy = dynamic_cast<const SimpleVectorBase<T>&>(*vy.vec);
      n2 = svy.length();
      assert( n2 >= 0 );

      assert( n1 + n2 <= N );
      if( n2 > 0 )
         memcpy(&sv[n1], &svy[0], n2 * sizeof(T));
   }

   if( vz.vec )
   {
      const SimpleVectorBase<T>& svz = dynamic_cast<const SimpleVectorBase<T>&>(*vz.vec);
      n3 = svz.length();
      assert( n3 >= 0 );

      assert( n1 + n2 + n3 <= N );
      if( n3 > 0 )
         memcpy(&sv[n1 + n2], &svz[0], n3 * sizeof(T));
   }

   if( vx.vecl )
   {
      const SimpleVectorBase<T>& svxl = dynamic_cast<const SimpleVectorBase<T>&>(*vx.vecl);
      n4 = svxl.length();
      assert( n4 >= 0 );

      assert( n1 + n2 + n3 + n4 <= N );
      if( n4 > 0 )
         memcpy(&sv[n1 + n2 + n3], &svxl[0], n4 * sizeof(T) );
   }

   if( vy.vecl )
   {
      const SimpleVectorBase<T>& svyl = dynamic_cast<const SimpleVectorBase<T>&>(*vy.vecl);
      n5 = svyl.length();
      assert( n5 >= 0 );

      assert( n1 + n2 + n3 + n4 + n5 <= N );
      if( n5 > 0 )
         memcpy(&sv[n1 + n2 + n3 + n4], &svyl[0], n5 * sizeof(T));
   }

   if( vz.vecl )
   {
      const SimpleVectorBase<T>& svzl = dynamic_cast<const SimpleVectorBase<T>&>(*vz.vecl);
      n6 = svzl.length();
      assert( n6 >= 0 );

      assert( n1 + n2 + n3 + n4 + n5 + n6 <= N );
      if( n6 > 0 )
         memcpy(&sv[n1 + n2 + n3 + n4 + n5], &svzl[0], n6 * sizeof(T));
   }

   assert( n1 + n2 + n3 + n4 + n5 + n6 == N );

   for(size_t it = 0; it < children.size(); it++)
      children[it]->jointCopyFrom(*vx.children[it], *vy.children[it], *vz.children[it]);
}

template<typename T>
void StochVectorBase<T>::jointCopyTo(StochVectorBase<T>& vx, StochVectorBase<T>& vy, StochVectorBase<T>& vz) const
{
   assert( this->vec );
   const SimpleVectorBase<T>& sv  = dynamic_cast<const SimpleVectorBase<T>&>(*this->vec);
   assert( sizeof(T) == sizeof(sv[0]) );

   const int N = sv.length();
   int n1 = 0;
   int n2 = 0;
   int n3 = 0;
   int n4 = 0;
   int n5 = 0;
   int n6 = 0;

   if( vx.vec )
   {
      SimpleVectorBase<T> &svx = dynamic_cast<SimpleVectorBase<T>&>(*vx.vec);
      n1 = svx.length();
      assert(n1 >= 0);

      assert(n1 <= N);
      if( n1 > 0 )
         memcpy(&svx[0], &sv[0], n1 * sizeof(T));
   }

   if( vy.vec )
   {
      SimpleVectorBase<T>& svy = dynamic_cast<SimpleVectorBase<T>&>(*vy.vec);
      n2 = svy.length();
      assert( n2 >= 0 );

      assert( n1 + n2 <= N );
      if( n2 > 0 )
         memcpy(&svy[0], &sv[n1], n2 * sizeof(T));
   }

   if( vz.vec )
   {
      SimpleVectorBase<T>& svz = dynamic_cast<SimpleVectorBase<T>&>(*vz.vec);
      n3 = svz.length();
      assert( n3 >= 0 );

      assert( n1 + n2 + n3 <= N );
      if( n3 > 0 )
         memcpy(&svz[0], &sv[n1 + n2], n3 * sizeof(T));
   }

   if( vx.vecl )
   {
      SimpleVectorBase<T>& svxl = dynamic_cast<SimpleVectorBase<T>&>(*vx.vecl);
      n4 = svxl.length();
      assert( n4 >= 0 );

      assert( n1 + n2 + n3 + n4 <= N );
      if( n4 > 0 )
         memcpy(&svxl[0], &sv[n1 + n2 + n3], n4 * sizeof(T) );
   }

   if( vy.vecl )
   {
      SimpleVectorBase<T>& svyl = dynamic_cast<SimpleVectorBase<T>&>(*vy.vecl);
      n5 = svyl.length();
      assert(n5 >= 0);

      assert( n1 + n2 + n3 + n4 + n5 <= N );
      if( n5 > 0 )
         memcpy(&svyl[0], &sv[n1 + n2 + n3 + n4], n5 * sizeof(T));
   }

   if( vz.vecl )
   {
      SimpleVectorBase<T>& svzl = dynamic_cast<SimpleVectorBase<T>&>(*vz.vecl);
      n6 = svzl.length();
      assert(n6 >= 0);

      assert( n1 + n2 + n3 + n4 + n5 +n6 <= N );
      if( n6 > 0 )
         memcpy(&svzl[0], &sv[n1 + n2 + n3 + n4 + n5], n6 * sizeof(T));
   }      assert( n1 + n2 + n3 + n4 + n5 <= N );

   assert( n1 + n2 + n3 + n4 + n5 + n6 == N );

   for(size_t it = 0; it < children.size(); it++)
      children[it]->jointCopyTo(*vx.children[it], *vy.children[it], *vz.children[it]);
}


template<typename T>
bool StochVectorBase<T>::isKindOf( int kind ) const
{
  return (kind == kStochVector);
}

template<typename T>
void StochVectorBase<T>::scale( T alpha )
{
   if( vec )
      vec->scale(alpha);

   if( vecl )
      vecl->scale(alpha);

   for(size_t it = 0; it < children.size(); it++)
      children[it]->scale(alpha);
}

template<typename T>
bool StochVectorBase<T>::isZero() const
{
	bool is_zero = true;

	for( size_t node = 0; node < children.size(); ++node )
	{
	   const bool is_zero_tmp = children[node]->isZero();
		is_zero = ( is_zero && is_zero_tmp );
	}

	PIPS_MPIgetLogicAndInPlace(is_zero, mpiComm);

	if( vec )
	{
	   const bool is_zero_tmp = dynamic_cast<SimpleVectorBase<T>&>(*vec).isZero();
	   is_zero = is_zero && is_zero_tmp;
	}

	if( vecl )
   {
	   const bool is_zero_tmp = dynamic_cast<SimpleVectorBase<T>&>(*vecl).isZero();
	   is_zero = is_zero && is_zero_tmp;
   }

	return is_zero;
}

template<typename T>
void StochVectorBase<T>::setToZero()
{
   if( vec )
      vec->setToZero();

   if( vecl )
      vecl->setToZero();

   for(size_t it = 0; it < children.size(); it++)
      children[it]->setToZero();
}

template<typename T>
void StochVectorBase<T>::setToConstant(T c)
{
   if( vec )
      vec->setToConstant(c);

   if( vecl )
      vecl->setToConstant(c);

   for(size_t it = 0; it < children.size(); it++)
      children[it]->setToConstant(c);
}

template<typename T>
void StochVectorBase<T>::randomize( T alpha, T beta, T *ix )
{
  assert( "Not implemented" && 0 );
}


template<typename T>
void StochVectorBase<T>::copyFrom( const OoqpVectorBase<T>& v_ )
{
   const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

   if( this->vec )
   {
      assert( v.vec );
      this->vec->copyFrom(*v.vec);
   }
   else
      assert( v.vec == nullptr );

   if( this->vecl )
   {
      assert( v.vecl );
      this->vecl->copyFrom(*v.vecl);
   }
   else
      assert( v.vecl == nullptr );

   assert( children.size() == v.children.size() );

   for(size_t it = 0; it < children.size(); it++)
      children[it]->copyFrom( *v.children[it] );
}

template<typename T>
void StochVectorBase<T>::copyFromAbs(const OoqpVectorBase<T>& v_ )
{
  const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

  if( this->vec )
  {
     assert( v.vec );
     this->vec->copyFromAbs(*v.vec);
  }
  else
     assert( v.vec == nullptr );

  if( this->vecl )
  {
     assert(v.vecl);
     this->vecl->copyFromAbs(*v.vecl);
  }
  else
     assert( v.vecl == nullptr );

  assert(children.size() == v.children.size());
  for(size_t it = 0; it < children.size(); it++)
    children[it]->copyFromAbs(*v.children[it]);
}

template<typename T>
T StochVectorBase<T>::infnorm() const
{
   T infnrm = 0.0;

   for(size_t it = 0; it < children.size(); it++)
      infnrm = std::max(infnrm, children[it]->infnorm());

   if( iAmDistrib )
      PIPS_MPIgetMaxInPlace(infnrm, mpiComm);

   if( vec )
      infnrm = std::max(vec->infnorm(), infnrm);

   if( vecl )
      infnrm = std::max(vecl->infnorm(), infnrm);

   return infnrm;
}

template<typename T>
double StochVectorBase<T>::twonorm() const
{
  const T scale = this->infnorm();
  assert(scale >= 0.0);

  if( PIPSisZero(scale) )
     return 0.0;

  return scale * std::sqrt( this->dotProductSelf( 1.0 / scale ) );
}

template<typename T>
T StochVectorBase<T>::onenorm() const
{
  T onenorm = 0.0;

  for( size_t it = 0; it < children.size(); it++ )
     onenorm += children[it]->onenorm();

  if( iAmSpecial && vec )
     onenorm += vec->onenorm();

  if( iAmSpecial && vecl )
     onenorm += vecl->onenorm();

  if( iAmDistrib && parent == nullptr )
     PIPS_MPIgetSumInPlace(onenorm, mpiComm);

  return onenorm;
}


template<typename T>
void StochVectorBase<T>::min( T& m, int& index ) const
{
   // index is broken for StochVector
   index = -1;

   m = std::numeric_limits<T>::max();

   if( vec )
   {
      T min_tmp;
      vec->min( min_tmp, index );

      m = std::min( min_tmp, m );
   }

   if( vecl )
   {
      T min_tmp;
      vecl->min(min_tmp, index);

      m = std::min( min_tmp, m );
   }

   for(size_t it = 0; it < children.size(); it++)
   {
      T min_tmp;
      children[it]->min(min_tmp, index);

      m = std::min( min_tmp, m );
   }

   if( iAmDistrib )
      PIPS_MPIgetMinInPlace(m, mpiComm);
}

template<typename T>
void StochVectorBase<T>::max( T& m, int& index ) const
{
   // index is broken for StochVector
   index = -1;

   m = -std::numeric_limits<T>::max();
   if( vec )
   {
      T max_tmp;
      vec->max( max_tmp, index );

      m = std::max( m, max_tmp );
   }

   if( vecl )
   {
      T max_tmp;
      vecl->max(max_tmp, index);

      m = std::max( m, max_tmp );
   }

   for( size_t it = 0; it < children.size(); it++ )
   {
      T max_tmp;
      children[it]->max(max_tmp, index);

      m = std::max( m, max_tmp );
   }

   if( iAmDistrib )
      PIPS_MPIgetMaxInPlace(m, mpiComm);
}

template<typename T>
void StochVectorBase<T>::absminVecUpdate(OoqpVectorBase<T>& absminvec) const
{
   StochVectorBase<T>& absminvecStoch = dynamic_cast<StochVectorBase<T>&>(absminvec);

   if( vec )
   {
      assert( absminvecStoch.vec );
      vec->absminVecUpdate(*(absminvecStoch.vec));
   }
   else
      assert( absminvecStoch.vec == nullptr );

   if( vecl )
   {
      assert( absminvecStoch.vecl );
      vecl->absminVecUpdate( *absminvecStoch.vecl );
   }
   else
      assert( absminvecStoch.vecl == nullptr );

   assert( absminvecStoch.children.size() == children.size() );
   for( size_t it = 0; it < children.size(); it++ )
      children[it]->absminVecUpdate( *absminvecStoch.children[it] );
}

template<typename T>
void StochVectorBase<T>::absmaxVecUpdate(OoqpVectorBase<T>& absmaxvec) const
{
   StochVectorBase<T>& absmaxvecStoch = dynamic_cast<StochVectorBase<T>&>(absmaxvec);

   if( vec )
   {
      assert( absmaxvecStoch.vec );
      vec->absmaxVecUpdate( *absmaxvecStoch.vec);
   }
   else
      assert( absmaxvecStoch.vec == nullptr );

   if( vecl )
   {
      assert( absmaxvecStoch.vecl );
      vecl->absmaxVecUpdate( *absmaxvecStoch.vecl );
   }
   else
      assert( absmaxvecStoch.vecl == nullptr );

   assert( absmaxvecStoch.children.size() == children.size() );
   for( size_t it = 0; it < children.size(); it++ )
      children[it]->absmaxVecUpdate( *absmaxvecStoch.children[it] );
}

template<typename T>
void StochVectorBase<T>::absmin(T& m) const
{
   T min_tmp = std::numeric_limits<T>::infinity();

   if( vec )
   {
      vec->absmin( min_tmp );
      m = std::min( m, min_tmp );
   }

   if( vecl )
   {
      vecl->absmin( min_tmp );
      m = std::min( m, min_tmp );
   }

   for( size_t it = 0; it < children.size(); ++it )
   {
      children[it]->absmin( min_tmp );
      m = std::min( m, min_tmp );
   }

   if( iAmDistrib )
      PIPS_MPIgetMinInPlace(m, mpiComm);

   assert( m >= 0.0 );
}

template<typename T>
void StochVectorBase<T>::absminNonZero(T& m, T zero_eps) const
{
   assert(zero_eps >= 0.0);
   m = std::numeric_limits<T>::infinity();
   T min_tmp = std::numeric_limits<T>::infinity();

   if( vec )
   {
      vec->absminNonZero(min_tmp, zero_eps);
      m = std::min( min_tmp, m );
   }

   if( vecl )
   {
      vecl->absminNonZero(min_tmp, zero_eps);
      m = std::min( min_tmp, m );
   }

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->absminNonZero(min_tmp, zero_eps);
      m = std::min( min_tmp, m );
   }

   if( iAmDistrib )
      PIPS_MPIgetMinInPlace(m, mpiComm);

   assert( m >= zero_eps );
}


template<typename T>
T StochVectorBase<T>::stepbound(const OoqpVectorBase<T> & v_, T maxStep ) const
{
  const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

  T step = 1.0;

  if( vec )
  {
     assert( v.vec );
     T stepvec = this->vec->stepbound(*v.vec, maxStep);
     if( stepvec < step )
        step = stepvec;
  }

  if( vecl )
  {
     assert(v.vecl);
     T stepl = vecl->stepbound(*v.vecl, maxStep);
     if( stepl < step )
        step = stepl;
  }

  //check tree compatibility
  assert( children.size() == v.children.size() );

  for(size_t it = 0; it < children.size(); it++)
    step = children[it]->stepbound(*v.children[it], step);

  if(iAmDistrib )
     PIPS_MPIgetMinInPlace(step, mpiComm);

  return step;
}

template<typename T>
T StochVectorBase<T>::findBlocking(const OoqpVectorBase<T> & wstep_vec,
			      const OoqpVectorBase<T> & u_vec,
			      const OoqpVectorBase<T> & ustep_vec,
			      T maxStep,
			      T *w_elt,
			      T *wstep_elt,
			      T *u_elt,
			      T *ustep_elt,
			      int& first_or_second) const
{
  const StochVectorBase<T>& w = *this;
  const StochVectorBase<T>& u = dynamic_cast<const StochVectorBase<T>&>(u_vec);

  const StochVectorBase<T>& wstep = dynamic_cast<const StochVectorBase<T>&>(wstep_vec);
  const StochVectorBase<T>& ustep = dynamic_cast<const StochVectorBase<T>&>(ustep_vec);
  const T local_eps = 1e-14;

  T step = maxStep;

  // todo only if i am special?
  if( w.vecl )
  {
    assert(wstep.vecl);
    assert(u.vecl);
    assert(ustep.vecl);

    step = w.vecl->findBlocking(*wstep.vecl, *u.vecl, *ustep.vecl, step,
                 w_elt, wstep_elt, u_elt, ustep_elt,
                 first_or_second);
  }

  if( w.vec )
  {
     assert(wstep.vec);
     assert(u.vec);
     assert(ustep.vec);

     step = w.vec->findBlocking(*wstep.vec, *u.vec, *ustep.vec, step,
                  w_elt, wstep_elt, u_elt, ustep_elt,
                  first_or_second);
  }

  const int nChildren = w.children.size();
  //check tree compatibility
  assert( nChildren - u.children.size() == 0);
  assert( wstep.children.size() == ustep.children.size() );
  assert( nChildren - ustep.children.size() == 0);

  for(int it = 0; it < nChildren; it++)
  {
     step = w.children[it]->findBlocking(*wstep.children[it], *u.children[it], *ustep.children[it], step,
           w_elt, wstep_elt, u_elt,ustep_elt, first_or_second);
  }

  if(iAmDistrib == 1) {
    T stepG;
    assert(PIPSisLE(step, 1.0));
    assert(PIPSisLE(0.0, step));

    stepG = PIPS_MPIgetMin(step, mpiComm);
    const bool iHaveMinStep = PIPSisEQ(step, stepG, local_eps);

    //we prefer a AllReduce instead of a bcast, since the step==stepG m
    //may occur for two different processes and a deadlock may occur.
    T buffer[5]; //0-primal val, 1-primal step, 2-dual value, 3-step, 4-1st or 2nd

    int count;
    if( iHaveMinStep ) {
      buffer[0]=*w_elt; buffer[1]=*wstep_elt;
      buffer[2]=*u_elt; buffer[3]=*ustep_elt;
      buffer[4]=first_or_second;

      count = 1;
    } else {

      count = 0;
      buffer[0]=buffer[1]=buffer[2]=buffer[3]=buffer[4]= -std::numeric_limits<T>::max();
    }

    count = PIPS_MPIgetSum(count, mpiComm);
    // MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, mpiComm); // not working properly in templated version
    assert(count >= 1);

    // is there more than one process with step==stepG?
    if( count > 1 )
    {
       int myrank;
       int mineqrank;

       MPI_Comm_rank(mpiComm, &myrank);

       if( iHaveMinStep )
          mineqrank = myrank;
       else
          mineqrank = std::numeric_limits<int>::max();

       mineqrank = PIPS_MPIgetMin(mineqrank, mpiComm);
       // MPI_Allreduce(MPI_IN_PLACE, &mineqrank, 1, MPI_INT, MPI_MIN, mpiComm); // not working properly in templated version

       // step==stepG and not smallest rank?
      if( iHaveMinStep && mineqrank != myrank )
         buffer[0]=buffer[1]=buffer[2]=buffer[3]=buffer[4]= -std::numeric_limits<T>::max();
    }

    T bufferOut[5];
    PIPS_MPImaxArray(buffer, bufferOut, 5, mpiComm);

    *w_elt = bufferOut[0]; *wstep_elt=bufferOut[1];
    *u_elt = bufferOut[2]; *ustep_elt=bufferOut[3];

    // negative or 0 means no blocking, so set first_or_second to 0.
    if( bufferOut[4] <= 0.5 )
       first_or_second = 0;
    else if( bufferOut[4] <= 1.5 )
       first_or_second = 1;
    else
       first_or_second = 2;

    step=stepG;
  }
  return step;
}

template<typename T>
void StochVectorBase<T>::findBlocking_pd(const OoqpVectorBase<T> & wstep_vec,
			      const OoqpVectorBase<T> & u_vec,
			      const OoqpVectorBase<T> & ustep_vec,
			      T& maxStepPri, T& maxStepDual,
			      T& w_elt_p, T& wstep_elt_p, T& u_elt_p, T& ustep_elt_p,
				   T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d,
				   bool& primalBlocking, bool& dualBlocking) const
{
  const StochVectorBase<T>& w = *this;
  const StochVectorBase<T>& u = dynamic_cast<const StochVectorBase<T>&>(u_vec);
  const T local_eps = 1e-14;

  const StochVectorBase<T>& wstep = dynamic_cast<const StochVectorBase<T>&>(wstep_vec);
  const StochVectorBase<T>& ustep = dynamic_cast<const StochVectorBase<T>&>(ustep_vec);

  // todo only if i am special?
  if( w.vecl )
  {
    assert(wstep.vecl);
    assert(u.vecl);
    assert(ustep.vecl);

    w.vecl->findBlocking_pd(*wstep.vecl, *u.vecl, *ustep.vecl, maxStepPri, maxStepDual,
                 w_elt_p, wstep_elt_p, u_elt_p, ustep_elt_p,
				     w_elt_d, wstep_elt_d, u_elt_d, ustep_elt_d,
                 primalBlocking, dualBlocking);
  }

  if( w.vec )
  {
     assert( wstep.vec );
     assert( u.vec );
     assert( ustep.vec );

     w.vec->findBlocking_pd(*wstep.vec, *u.vec, *ustep.vec, maxStepPri, maxStepDual,
		  	  	  w_elt_p, wstep_elt_p, u_elt_p, ustep_elt_p,
				  w_elt_d, wstep_elt_d, u_elt_d, ustep_elt_d,
				  primalBlocking, dualBlocking);
  }

  int nChildren=w.children.size();
  //check tree compatibility
  assert( nChildren             - u.children.size() == 0);
  assert( wstep.children.size() == ustep.children.size() );
  assert( nChildren             - ustep.children.size() == 0);

  for(int it=0; it<nChildren; it++) {
    w.children[it]->findBlocking_pd(*wstep.children[it],
			       *u.children[it],
			       *ustep.children[it],
			       maxStepPri, maxStepDual,
				    w_elt_p, wstep_elt_p, u_elt_p, ustep_elt_p,
				    w_elt_d, wstep_elt_d, u_elt_d, ustep_elt_d,
				    primalBlocking, dualBlocking);
   }

   if( iAmDistrib == 1 )
   {
      T maxStepGlobalPri, maxStepGlobalDual;
      assert(PIPSisLE(maxStepPri, 1.0) && PIPSisLE(maxStepDual, 1.0));
      assert(PIPSisLE(0.0, maxStepPri) && PIPSisLE(0.0, maxStepDual));

      maxStepGlobalPri = PIPS_MPIgetMin(maxStepPri, mpiComm);
      maxStepGlobalDual = PIPS_MPIgetMin(maxStepDual, mpiComm);
      const bool iHaveMinStepPri = PIPSisEQ(maxStepPri, maxStepGlobalPri, local_eps);
      const bool iHaveMinStepDual = PIPSisEQ(maxStepDual, maxStepGlobalDual, local_eps);

      //we prefer a AllReduce instead of a bcast, since the step==stepG
      //may occur for two different processes and a deadlock may occur.
      T buffer[10];
      int count[2];
      //values for computation of the primal steplength:
      //0-primal val, 1-primal step, 2-dual value, 3-dual step, 4-primalBlocking
      if( iHaveMinStepPri )
      {
         buffer[0] = w_elt_p;
         buffer[1] = wstep_elt_p;
         buffer[2] = u_elt_p;
         buffer[3] = ustep_elt_p;
         buffer[4] = primalBlocking ? 1.0 : 0.0;

         count[0] = 1;
      }
      else
      {
         buffer[0] = buffer[1] = buffer[2] = buffer[3] = buffer[4] =
               -std::numeric_limits<T>::max();
         count[0] = 0;
      }

      //values for computation of the dual steplength:
      //5-primal val, 6-primal step, 7-dual value, 8-dual step, 9-dualBlocking
      if( iHaveMinStepDual )
      {
         buffer[5] = w_elt_d;
         buffer[6] = wstep_elt_d;
         buffer[7] = u_elt_d;
         buffer[8] = ustep_elt_d;
         buffer[9] = dualBlocking ? 1.0 : 0.0;

         count[1] = 1;
      }
      else
      {
         buffer[5] = buffer[6] = buffer[7] = buffer[8] = buffer[9] =
               -std::numeric_limits<T>::max();
         count[1] = 0;
      }

      PIPS_MPIsumArrayInPlace(count, 2, mpiComm);
      // MPI_Allreduce(MPI_IN_PLACE, count, 2, MPI_INT, MPI_SUM, mpiComm);

      assert(count[0] >= 1 && count[1] >= 1);

      int myrank;
      MPI_Comm_rank(mpiComm, &myrank);

      // is there more than one process with maxStepPri==stepG?
      if( count[0] > 1 )
      {
         int mineqrank;

         if( iHaveMinStepPri )
            mineqrank = myrank;
         else
            mineqrank = std::numeric_limits<int>::max();

          mineqrank = PIPS_MPIgetMin(mineqrank, mpiComm);
         // MPI_Allreduce(MPI_IN_PLACE, &mineqrank, 1, MPI_INT, MPI_MIN, mpiComm);

         // step==stepG and not smallest rank?
         if( iHaveMinStepPri && mineqrank != myrank )
            buffer[0] = buffer[1] = buffer[2] = buffer[3]=buffer[4]= -std::numeric_limits<T>::max();
      }

      // is there more than one process with maxStepDual==stepF?
      if( count[1] > 1 )
      {
         int mineqrank;

         if( iHaveMinStepDual )
            mineqrank = myrank;
         else
            mineqrank = std::numeric_limits<int>::max();

          mineqrank = PIPS_MPIgetMin(mineqrank, mpiComm);
         // MPI_Allreduce(MPI_IN_PLACE, &mineqrank, 1, MPI_INT, MPI_MIN, mpiComm);

         // stepDual==stepF and not smallest rank?
         if( iHaveMinStepDual && mineqrank != myrank )
            buffer[5]=buffer[6]=buffer[7]=buffer[8]=buffer[9]= -std::numeric_limits<T>::max();
      }

      T bufferOut[10];
      PIPS_MPImaxArray(buffer, bufferOut, 10, mpiComm);

      w_elt_p = bufferOut[0];
      wstep_elt_p = bufferOut[1];
      u_elt_p = bufferOut[2];
      ustep_elt_p = bufferOut[3];

      w_elt_d = bufferOut[5];
      wstep_elt_d = bufferOut[6];
      u_elt_d = bufferOut[7];
      ustep_elt_d = bufferOut[8];

      primalBlocking = bufferOut[4] <= 0.5 ? false : true;
      maxStepPri = maxStepGlobalPri;

      dualBlocking = bufferOut[9] <= 0.5 ? false : true;
      maxStepDual = maxStepGlobalDual;
   }
}

template<typename T>
void StochVectorBase<T>::componentMult( const OoqpVectorBase<T>& v_ )
{
   const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

   if( vec )
   {
      assert( v.vec );
      vec->componentMult( *v.vec );
   }
   else
      assert( v.vec == nullptr );

   if( vecl )
   {
      assert( v.vecl );
      vecl->componentMult( *v.vecl );
   }
   else
      assert( v.vecl == nullptr );

   assert( v.children.size() == children.size() );
   for( size_t it = 0; it < children.size(); it++ )
      children[it]->componentMult( *v.children[it] );
}

template<typename T>
void StochVectorBase<T>::componentDiv ( const OoqpVectorBase<T>& v_ )
{
   const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

   if( vec )
   {
      assert( v.vec );
      vec->componentDiv( *v.vec );
   }
   else
      assert( v.vec == nullptr );

   if( vecl )
   {
      assert( v.vecl );
      vecl->componentDiv( *v.vecl );
   }
   else
      assert( v.vecl == nullptr );

   assert( v.children.size() == children.size() );
   for( size_t it = 0; it < children.size(); ++it )
      children[it]->componentDiv( *v.children[it] );
}

template<typename T>
bool StochVectorBase<T>::componentEqual( const OoqpVectorBase<T>& v_, T tol) const
{
   const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);
   assert(v.children.size() == children.size());

   bool component_equal = true;

   if( vec )
   {
      assert( v.vec );
      const bool component_equal_tmp = vec->componentEqual(*v.vec, tol);
      component_equal = component_equal && component_equal_tmp;
   }
   else
      assert( v.vec == nullptr );

   if( !component_equal )
      if( parent == nullptr )
        std::cout << "not equal in root node non-link" << std::endl;

   if( vecl )
   {
      assert( v.vecl );
      const bool component_equal_tmp = vecl->componentEqual(*v.vecl, tol);
      component_equal = component_equal && component_equal_tmp;
   }
   else
      assert( v.vecl == nullptr );

   if( !component_equal )
      if( parent == nullptr )
         std::cout << "not equal in root node link" << std::endl;

   for(size_t child = 0; child < children.size(); child++)
   {
      const bool component_equal_tmp = children[child]->componentEqual(*v.children[child], tol);
      component_equal = component_equal_tmp && component_equal;

      if( !component_equal )
         std::cout << "not equal in root child node " << child << std::endl;
   }

   if( iAmDistrib )
      PIPS_MPIgetLogicAndInPlace(component_equal, mpiComm);

   return component_equal;
}

template<typename T>
bool StochVectorBase<T>::componentNotEqual( const T val, const T tol ) const
{
   bool not_equal = true;

   if( vec )
   {
      const bool not_equal_tmp = vec->componentNotEqual( val, tol );
      not_equal = not_equal_tmp && not_equal;
   }

   if( !not_equal )
      if( parent == nullptr )
        std::cout << "equal in root node non-link" << std::endl;

   if( vecl )
   {
      const bool not_equal_tmp = vecl->componentNotEqual(val, tol);
      not_equal = not_equal_tmp && not_equal;
   }

   if( !not_equal )
      std::cout << "equal in root node link" << std::endl;

   for(size_t child = 0; child < children.size(); child++)
   {
      const bool not_equal_tmp = children[child]->componentNotEqual(val, tol);
      not_equal = not_equal && not_equal_tmp;

      if( !not_equal )
         std::cout << "equal in root child node " << child << std::endl;
   }

   if( iAmDistrib )
      PIPS_MPIgetLogicAndInPlace(not_equal, mpiComm);

   return not_equal;
}

template<typename T>
void StochVectorBase<T>::scalarMult( T num )
{
   if( vec )
      vec->scalarMult(num);
   if( vecl )
      vecl->scalarMult(num);

  for( size_t it = 0; it < children.size(); ++it )
     children[it]->scalarMult( num );
}

template<typename T>
void StochVectorBase<T>::writeToStream( std::ostream& out, int offset ) const
{
   // TODO modify for hierarchical approach
   const int rank = PIPS_MPIgetRank(mpiComm);
   const int world_size = PIPS_MPIgetSize(mpiComm);

   MPI_Status status;
   int l;
   std::stringstream sout;
   if( rank == 0 )
   {
      for( int i = 0; i < offset; ++i )
         sout << "\t";
      sout << "--vec--" << std::endl;
      if( vec )
         vec->writeToStream( sout, offset );
      for( int i = 0; i < offset; ++i )
         sout << "\t";
      sout << "-------" << std::endl;
   }

   if( rank == 0 )
   {
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->writeToStream( sout, offset + 1 );

      out << sout.str();
      sout.str( std::string() );

      for( int p = 1; p < world_size; p++ )
      {
         MPI_Probe(p, p, mpiComm, &status);
         MPI_Get_count(&status, MPI_CHAR, &l);
         char *buf = new char[l];
         MPI_Recv(buf, l, MPI_CHAR, p, p, mpiComm, &status);
         std::string rowPartFromP(buf, l);
         out << rowPartFromP;
         delete[] buf;
      }
   }
   else if( iAmDistrib == 1 )
   {
      // rank != 0
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->writeToStream( sout, offset + 1 );

      std::string str = sout.str();
      MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, rank, mpiComm);
   }

   if( rank == 0 )
   {
      for( int i = 0; i < offset; ++i )
         sout << "\t";
      sout << "--vecl-" << std::endl;

      if( vecl )
         vecl->writeToStream( sout, offset );

      for( int i = 0; i < offset; ++i )
         sout << "\t";
      sout << "-------" << std::endl;

      out << sout.str();
   }

   if( iAmDistrib == 1 )
      MPI_Barrier( mpiComm );
}

template<typename T>
void StochVectorBase<T>::pushAwayFromZero( double tol, double amount, const OoqpVectorBase<T>* select )
{
   const StochVectorBase<T>* selects = select ? dynamic_cast<const StochVectorBase<T>*>(select) : nullptr;

   if( vec )
      vec->pushAwayFromZero( tol, amount, selects ? selects->vec : nullptr );

   if( vecl )
      vecl->pushAwayFromZero( tol, amount, selects ? selects->vecl : nullptr );

   for( size_t i = 0; i < this->children.size(); ++i )
      this->children[i]->pushAwayFromZero( tol, amount, selects ? selects->children[i] : nullptr );
}

template<typename T>
void StochVectorBase<T>::getSumCountIfSmall( double tol, double& sum_small, int& n_close, const OoqpVectorBase<T>* select ) const
{
   const StochVectorBase<T>* selects = dynamic_cast<const StochVectorBase<T>*>(select);

   for( size_t i = 0; i < this->children.size(); ++i )
      this->children[i]->getSumCountIfSmall( tol, sum_small, n_close, selects ? selects->children[i] : nullptr );

   if( iAmSpecial && vec )
   {
      if( selects )
         assert(selects->vec);

      vec->getSumCountIfSmall( tol, sum_small, n_close, selects ? selects->vec : nullptr );
   }

   if( iAmSpecial && vecl )
   {
      if( selects )
         assert(selects->vecl);
      vecl->getSumCountIfSmall( tol, sum_small, n_close, selects ? selects->vecl : nullptr );
   }

   if( iAmDistrib && parent == nullptr )
   {
      PIPS_MPIgetSumInPlace(sum_small, mpiComm);
      PIPS_MPIgetSumInPlace(n_close, mpiComm);
   }
}

template<typename T>
void StochVectorBase<T>::writefToStream( std::ostream& out,
				  const char format[] ) const
{
  vec->writefToStream(out, format);
  if( vecl ) vecl->writefToStream(out, format);

  for(size_t it=0; it<children.size(); it++)
    children[it]->writefToStream(out, format);
}

template<typename T>
void StochVectorBase<T>::writeMPSformatRhs( std::ostream& out, int rowType, const OoqpVectorBase<T>* irhs) const
{
   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);
   std::string rt;
   if( rowType == 0 )
      rt = "E";
   else if( rowType == 1 )
      rt = "L";
   else if( rowType == 2 )
      rt = "G";
   else
      assert(0);

   const StochVectorBase<T>* ic = nullptr;
   if( irhs )
      ic = dynamic_cast<const StochVectorBase<T>*>(irhs);

   if( myRank == 0 )
   {
      std::string rowNameStub = " B row_";
      rowNameStub += rt;
      rowNameStub += "_R_";
      if( irhs && ic )
         vec->writeMPSformatOnlyRhs( out, rowNameStub, dynamic_cast<const SimpleVectorBase<T>*>(ic->vec));
      else
         vec->writeMPSformatOnlyRhs( out, rowNameStub, nullptr);
      if(vecl)
      {
         rowNameStub = " B row_";
         rowNameStub += rt;
         rowNameStub += "_L_";
         if( irhs )
            vecl->writeMPSformatOnlyRhs( out, rowNameStub, dynamic_cast<const SimpleVectorBase<T>*>(ic->vecl));
         else
            vecl->writeMPSformatOnlyRhs( out, rowNameStub, nullptr);
      }
   }
   for(int it=0; it<(int)children.size(); it++)
   {
      std::stringstream sstm;
      sstm << " B row_" << rt << "_" << it << "_";
      std::string rowNameStub = sstm.str();
      if( irhs )
         children[it]->vec->writeMPSformatOnlyRhs( out, rowNameStub, dynamic_cast<const SimpleVectorBase<T>*>(ic->children[it]->vec));
      else
         children[it]->vec->writeMPSformatOnlyRhs( out, rowNameStub, nullptr);
   }
}

template<typename T>
void StochVectorBase<T>::writeMPSformatBounds(std::ostream& out, const OoqpVectorBase<T>* ix, bool upperBound) const
{
   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);

   const StochVectorBase<T>* ixStoch = dynamic_cast<const StochVectorBase<T>*>(ix);

   if( myRank==0 )
   {
      std::string varNameStub = "var_L_";
      vec->writeMPSformatBoundsWithVar(out, varNameStub, (ixStoch->vec), upperBound);
   }
   for(int it=0; it<(int)children.size(); it++)
   {
      std::stringstream sstm2;
      sstm2 << "var_" << it << "_";
      std::string varNameStub = sstm2.str();
      children[it]->vec->writeMPSformatBoundsWithVar(out, varNameStub, (ixStoch->children[it]->vec), upperBound);
   }
}

/** this += alpha * x */
template<typename T>
void StochVectorBase<T>::axpy ( T alpha, const OoqpVectorBase<T>& x_ )
{
   const StochVectorBase<T>& x = dynamic_cast<const StochVectorBase<T>&>(x_);

   if( alpha == 0.0 )
      return;

   if( vec )
   {
      assert( x.vec );
      vec->axpy( alpha, *x.vec );
   }
   else
      assert( x.vec == nullptr );

   if( vecl )
   {
      assert( x.vecl );
      vecl->axpy( alpha, *x.vecl );
   }
   else
      assert( x.vecl == nullptr );

   assert( x.children.size() == children.size() );
   for( size_t it = 0; it < children.size(); ++it )
      children[it]->axpy( alpha, *x.children[it] );
}

/** this += alpha * x * z */
template<typename T>
void StochVectorBase<T>::axzpy( T alpha, const OoqpVectorBase<T>& x_, const OoqpVectorBase<T>& z_ )
{
  const StochVectorBase<T>& x = dynamic_cast<const StochVectorBase<T>&>(x_);
  const StochVectorBase<T>& z = dynamic_cast<const StochVectorBase<T>&>(z_);

  if( vec )
  {
     assert( x.vec );
     assert( z.vec );
     vec->axzpy( alpha, *x.vec, *z.vec );
  }
  else
  {
     assert( x.vec == nullptr );
     assert( z.vec == nullptr );
  }

  if( vecl )
  {
     assert( x.vecl );
     assert( z.vecl );
     vecl->axzpy( alpha, *x.vecl, *z.vecl );
  }
  else
  {
     assert( x.vecl == nullptr );
     assert( z.vecl == nullptr );
  }

  assert( x.children.size() == children.size() );
  assert( z.children.size() == children.size() );
  for( size_t it = 0; it < children.size(); ++it )
     children[it]->axzpy( alpha, *x.children[it], *z.children[it] );
}

/** this += alpha * x / z */
template<typename T>
void StochVectorBase<T>::axdzpy( T alpha, const OoqpVectorBase<T>& x_, const OoqpVectorBase<T>& z_ )
{
   const StochVectorBase<T>& x = dynamic_cast<const StochVectorBase<T>&>(x_);
   const StochVectorBase<T>& z = dynamic_cast<const StochVectorBase<T>&>(z_);

   if( vec )
   {
      assert( x.vec );
      assert( z.vec );
      vec->axdzpy( alpha, *x.vec, *z.vec );
   }
   else
   {
      assert( x.vec == nullptr );
      assert( z.vec == nullptr );
   }

   if( vecl )
   {
      assert(x.vecl);
      assert(z.vecl);
      vecl->axdzpy(alpha, *x.vecl, *z.vecl );
   }
   else
   {
      assert( x.vecl == nullptr );
      assert( z.vecl == nullptr );
   }

   assert( x.children.size() == children.size() );
   assert( z.children.size() == children.size() );
   for( size_t it = 0; it < children.size(); it++ )
      children[it]->axdzpy( alpha, *x.children[it], *z.children[it] );
}


template<typename T>
void StochVectorBase<T>::addConstant( T c )
{
   if( vec )
      vec->addConstant(c);

   if( vecl )
      vecl->addConstant(c);

  for( size_t it = 0; it < children.size(); it++ )
     children[it]->addConstant(c);
}


template<typename T>
void StochVectorBase<T>::gondzioProjection( T rmin, T rmax )
{
   if( vec )
      vec->gondzioProjection( rmin, rmax );

   if( vecl )
      vecl->gondzioProjection( rmin, rmax );

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->gondzioProjection( rmin, rmax );
}

template<typename T>
T StochVectorBase<T>::dotProductWith( const OoqpVectorBase<T>& v_ ) const
{
   const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

   T dot_product = 0.0;

   assert( v.children.size() == children.size() );
   for(size_t it = 0; it < children.size(); it++)
      dot_product += children[it]->dotProductWith( *v.children[it] );

   if( iAmSpecial && vec )
   {
      assert( v.vec );
      dot_product += vec->dotProductWith( *v.vec );
   }
   else if( !vec )
      assert( v.vec == nullptr );

   if( iAmSpecial && vecl )
   {
      assert( v.vecl );
      dot_product += vecl->dotProductWith( *v.vecl );
   }
   else if( !vecl )
      assert( v.vecl == nullptr );

   if( iAmDistrib && parent == nullptr )
      PIPS_MPIgetSumInPlace(dot_product, mpiComm);

   return dot_product;
}

template<typename T>
T StochVectorBase<T>::dotProductSelf(T scaleFactor) const
{
   T dot_product = 0.0;

   for( size_t it = 0; it < children.size(); it++ )
      dot_product += children[it]->dotProductSelf( scaleFactor );

   if( iAmSpecial && vec )
      dot_product += vec->dotProductSelf( scaleFactor );

   if( iAmSpecial && vecl )
      dot_product += vecl->dotProductSelf( scaleFactor );

   if( iAmDistrib && parent == nullptr )
      PIPS_MPIgetSumInPlace(dot_product, mpiComm);

   return dot_product;
}

/** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
 */
template<typename T>
T StochVectorBase<T>::shiftedDotProductWith( T alpha, const OoqpVectorBase<T>& mystep_,
					const OoqpVectorBase<T>& yvec_,
					T beta,  const OoqpVectorBase<T>& ystep_ ) const
{
   const StochVectorBase<T>& mystep = dynamic_cast<const StochVectorBase<T>&>(mystep_);
   const StochVectorBase<T>& yvec   = dynamic_cast<const StochVectorBase<T>&>(yvec_);
   const StochVectorBase<T>& ystep  = dynamic_cast<const StochVectorBase<T>&>(ystep_);

   T dot_product = 0.0;

   assert( this->children.size() == mystep.children.size() );
   assert( this->children.size() == yvec.children.size() );
   assert( this->children.size() == ystep.children.size() );

   for( size_t it = 0; it < children.size(); it++ )
      dot_product += children[it]->shiftedDotProductWith( alpha, *mystep.children[it], *yvec.children[it], beta, *ystep.children[it] );

   if( iAmSpecial && vec )
   {
      assert( mystep.vec );
      assert( yvec.vec );
      assert( ystep.vec );
      dot_product += vec->shiftedDotProductWith( alpha, *mystep.vec, *yvec.vec, beta, *ystep.vec );
   }
   else if( !vec )
   {
      assert( mystep.vec == nullptr );
      assert( yvec.vec == nullptr );
      assert( ystep.vec == nullptr );
   }

   if( iAmSpecial && vecl )
   {
      assert( mystep.vecl );
      assert( yvec.vecl );
      assert( ystep.vecl );
      dot_product += vecl->shiftedDotProductWith( alpha, *mystep.vecl, *yvec.vecl, beta, *ystep.vecl );
   }
   else if( !vecl )
   {
      assert( mystep.vecl == nullptr );
      assert( yvec.vecl == nullptr );
      assert( ystep.vecl == nullptr );
   }

   if( iAmDistrib && parent == nullptr )
      PIPS_MPIgetSumInPlace(dot_product, mpiComm);

   return dot_product;
}

template<typename T>
void StochVectorBase<T>::negate()
{
   if( vec )
      vec->negate();
   if( vecl )
      vecl->negate();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->negate();
}

template<typename T>
void StochVectorBase<T>::invert()
{
   if( vec )
      vec->invert();

   if( vecl )
      vecl->invert();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->invert();
}

template<typename T>
void StochVectorBase<T>::invertSave(T zeroReplacementVal)
{
   if( vec )
      vec->invertSave( zeroReplacementVal );

   if( vecl )
      vecl->invertSave( zeroReplacementVal );

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->invertSave( zeroReplacementVal );
}

template<typename T>
void StochVectorBase<T>::applySqrt()
{
   if( vec )
      vec->applySqrt();

   if( vecl )
      vecl->applySqrt();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->applySqrt();
}

template<typename T>
void StochVectorBase<T>::roundToPow2()
{
   if( vec )
      vec->roundToPow2();

   if( vecl )
      vecl->roundToPow2();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->roundToPow2();
}

template<typename T>
bool StochVectorBase<T>::allPositive() const
{
   bool all_pos = true;

   if( vec )
   {
      const bool all_pos_tmp = vec->allPositive();
      all_pos = all_pos && all_pos_tmp;
   }

   if( vecl )
   {
      const bool all_pos_tmp = vecl->allPositive();
      all_pos = all_pos && all_pos_tmp;
   }

   for(size_t it = 0; it < children.size(); it++)
   {
      const bool all_pos_tmp = children[it]->allPositive();
      all_pos = all_pos && all_pos_tmp;
   }

   if( iAmDistrib )
      PIPS_MPIgetLogicAndInPlace( all_pos, mpiComm );

   return all_pos;
}

template<typename T>
bool StochVectorBase<T>::matchesNonZeroPattern( const OoqpVectorBase<T>& select_ ) const
{
  const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);

  bool match = true;

  if( vec )
  {
     assert( select.vec );
     const bool match_tmp = vec->matchesNonZeroPattern( *select.vec );
     match = match && match_tmp;
  }
  else
     assert( select.vec == nullptr );

  if( vecl )
  {
     assert(select.vecl);
     const bool match_tmp = vecl->matchesNonZeroPattern( *select.vecl );
     match = match && match_tmp;
  }
  else
     assert( select.vecl == nullptr );

  assert( children.size() == select.children.size() );
  for( size_t it = 0; it < children.size() && match; it++ )
  {
     const bool match_tmp = children[it]->matchesNonZeroPattern( *select.children[it] );
     match = match && match_tmp;
  }

  if( iAmDistrib )
     PIPS_MPIgetLogicAndInPlace( match, mpiComm );

  return match;
}

template<typename T>
void StochVectorBase<T>::selectNonZeros( const OoqpVectorBase<T>& select_ )
{
   const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);

   if( vec )
   {
      assert( select.vec );
      vec->selectNonZeros( *select.vec );
   }
   else
      assert( select.vec == nullptr );
   if( vecl )
   {
      assert( select.vecl );
      vecl->selectNonZeros( *select.vecl );
   }
   else
      assert( select.vecl == nullptr );

   assert(children.size() == select.children.size());
   for( size_t it = 0; it < children.size(); it++ )
      children[it]->selectNonZeros( *select.children[it] );
}

template<typename T>
void StochVectorBase<T>::selectPositive()
{
   if( vec )
      vec->selectPositive();

   if( vecl )
      vecl->selectPositive();

   for(size_t it = 0; it < children.size(); it++)
     children[it]->selectPositive();
}

template<typename T>
void StochVectorBase<T>::selectNegative()
{
   if( vec )
      vec->selectNegative();

   if( vecl )
      vecl->selectNegative();

   for(size_t it = 0; it < children.size(); it++)
     children[it]->selectNegative();
}

template<typename T>
long long StochVectorBase<T>::numberOfNonzeros() const
{
   //!opt - store the number of nnz to avoid communication
  long long nnz = 0;

  for(size_t it = 0; it < children.size(); it++)
     nnz += children[it]->numberOfNonzeros();

  if( iAmSpecial && vec )
     nnz += vec->numberOfNonzeros();

  if( iAmSpecial && vecl )
     nnz += vecl->numberOfNonzeros();

  if( iAmDistrib && parent == nullptr )
     PIPS_MPIgetSumInPlace( nnz, mpiComm );

  return nnz;
}

template<typename T>
void
StochVectorBase<T>::addSomeConstants(T c, const OoqpVectorBase<T> &select_)
{
   const StochVectorBase<T> &select = dynamic_cast<const StochVectorBase<T>&>(select_);
   assert(children.size() == select.children.size());

   if( vec )
   {
      assert( select.vec );
      vec->addSomeConstants(c, *select.vec);
   }
   else
      assert( select.vec == nullptr );

   if( vecl )
   {
      assert(select.vecl);
      vecl->addSomeConstants(c, *select.vecl);
   }
   else
      assert( select.vecl == nullptr );

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->addSomeConstants(c, *select.children[it]);
}

template<typename T>
void StochVectorBase<T>::writefSomeToStream( std::ostream& out,
			 const char format[],
			 const OoqpVectorBase<T>& select_ ) const
{
   assert( "Not yet implemented" && 0 );
}

template<typename T>
void StochVectorBase<T>::axdzpy( T alpha, const OoqpVectorBase<T>& x_,
		       const OoqpVectorBase<T>& z_, const OoqpVectorBase<T>& select_ )
{
   const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);
   const StochVectorBase<T>& x = dynamic_cast<const StochVectorBase<T>&>(x_);
   const StochVectorBase<T>& z = dynamic_cast<const StochVectorBase<T>&>(z_);

   assert(children.size() == select.children.size());
   assert(children.size() == x.children.size());
   assert(children.size() == z.children.size());

   if( vec )
   {
      assert( x.vec );
      assert( z.vec );
      assert( select.vec );
      vec->axdzpy(alpha, *x.vec, *z.vec, *select.vec);
   }
   else
   {
      assert( x.vec == nullptr );
      assert( z.vec == nullptr );
      assert( select.vec == nullptr );
   }

   if( vecl )
   {
      assert( x.vecl );
      assert( z.vecl );
      assert( select.vecl );
      vecl->axdzpy(alpha, *x.vecl, *z.vecl, *select.vecl);
   }
   else
   {
      assert( x.vecl == nullptr );
      assert( z.vecl == nullptr );
      assert( select.vecl == nullptr );
   }

   for(size_t it = 0; it < children.size(); it++)
      children[it]->axdzpy(alpha, *x.children[it], *z.children[it], *select.children[it]);
}

template<typename T>
bool StochVectorBase<T>::somePositive( const OoqpVectorBase<T>& select_ ) const
{
   const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);
   assert( children.size() == select.children.size() );

   bool some_positive = true;;

   for(size_t it = 0; it < children.size(); it++)
   {
      const bool some_pos_tmp = children[it]->somePositive( *select.children[it] );
      some_positive = some_positive && some_pos_tmp;
   }

   if( iAmDistrib )
      PIPS_MPIgetLogicAndInPlace( some_positive, mpiComm );

   if( vec )
   {
      assert( select.vec );
      const bool some_pos_tmp = vec->somePositive(*select.vec);
      some_positive = some_positive && some_pos_tmp;
   }
   else
      assert( select.vec == nullptr );

   if( vecl )
   {
      assert(select.vecl);
      const bool some_pos_tmp = vecl->somePositive(*select.vecl);
      some_positive = some_positive && some_pos_tmp;
   }
   else
      assert( select.vecl == nullptr );

   return some_positive;
}

template<typename T>
void StochVectorBase<T>::divideSome( const OoqpVectorBase<T>& div_, const OoqpVectorBase<T>& select_ )
{
   const StochVectorBase<T> &div = dynamic_cast<const StochVectorBase<T>&>(div_);
   const StochVectorBase<T> &select = dynamic_cast<const StochVectorBase<T>&>(select_);

   assert(children.size() == div.children.size());
   assert(children.size() == select.children.size());

   if( vec )
   {
      assert(div.vec);
      assert(select.vec);
      vec->divideSome(*div.vec, *select.vec);
   }
   else
   {
      assert( div.vec == nullptr );
      assert( select.vec == nullptr );
   }

   if( vecl )
   {
      assert(div.vecl);
      assert(select.vecl);
      vecl->divideSome(*div.vecl, *select.vecl);
   }
   else
   {
      assert( div.vecl == nullptr );
      assert( select.vecl == nullptr );
   }

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->divideSome(*div.children[it], *select.children[it]);
}

template<typename T>
void StochVectorBase<T>::copyIntoArray( T v[] ) const
{
  assert( "Not supported" && 0 );
}

template<typename T>
void StochVectorBase<T>::copyFromArray( const T v[] )
{
  assert( "Not supported" && 0 );
}

template<typename T>
void StochVectorBase<T>::copyFromArray( const char v[] )
{
  assert( "Not supported" && 0 );
}

template<typename T>
void StochVectorBase<T>::removeEntries( const OoqpVectorBase<int>& select_ )
{
   const StochVectorBase<int>& select = dynamic_cast<const StochVectorBase<int>&>(select_);
   assert(children.size() == select.children.size());

   this->n = 0;

   if( vec )
   {
      assert( select.vec );
      vec->removeEntries(*select.vec);
      this->n = vec->n;
   }
   else
      assert( select.vec == nullptr );

   if( vecl )
   {
      assert( select.vecl );
      vecl->removeEntries(*select.vecl);
      this->n += vecl->n;
   }
   else
      assert( select.vecl == nullptr );

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->removeEntries(*select.children[it]);
      this->n += children[it]->n;
   }
}

/* uses PIPSisZero() as check instead of == 0.0 */
template<typename T>
int StochVectorBase<T>::getNnzs() const
{
   int non_zeros = 0;

   for( size_t it = 0; it < children.size(); it++ )
      non_zeros += children[it]->getNnzs();

   if( iAmSpecial && vec )
      non_zeros += vec->getNnzs();

   if( iAmSpecial && vecl )
      non_zeros += vecl->getNnzs();

   if( iAmDistrib && parent == nullptr )
      PIPS_MPIgetSumInPlace(non_zeros, mpiComm);

   return non_zeros;
}

template<typename T>
void StochVectorBase<T>::permuteVec0Entries(const std::vector<unsigned int>& permvec)
{
   if( vec )
      dynamic_cast<SimpleVectorBase<T>*>(vec)->permuteEntries(permvec);
}

template<typename T>
void StochVectorBase<T>::permuteLinkingEntries(const std::vector<unsigned int>& permvec)
{
   if( vecl )
      dynamic_cast<SimpleVectorBase<T>*>(vecl)->permuteEntries(permvec);
}

template<typename T>
std::vector<T> StochVectorBase<T>::gatherStochVector() const
{
#ifdef HIERARCHICAL
   // TODO adapt for hier approach
//   assert( false && "TODO : implement" );
#endif
   const SimpleVectorBase<T>& firstvec = dynamic_cast<const SimpleVectorBase<T>&>(*vec);
   const size_t nChildren = children.size();

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   const int my_size = PIPS_MPIgetSize(mpiComm);

   std::vector<T> gatheredVecLocal(0);

   for( size_t i = 0; i < nChildren; ++i )
   {
      const SimpleVectorBase<T>& vec = dynamic_cast<const SimpleVectorBase<T>&>(*children[i]->vec);

      if( vec.length() > 0 )
         gatheredVecLocal.insert(gatheredVecLocal.end(), &vec[0], &vec[0] + vec.length());
   }

   size_t solLength = firstvec.length();

   // final vector
   std::vector<T> gatheredVec(0);

   if( my_size > 0 )
   {
      // get all lengths
      std::vector<int> recvcounts(my_size);
      std::vector<int> recvoffsets(my_size);

      int mylength = int(gatheredVecLocal.size());

      PIPS_MPIallgather(&mylength, 1, &recvcounts[0], 1, mpiComm);

      // all-gather local components
      recvoffsets[0] = 0;
      for( size_t i = 1; i < size_t(my_size); ++i )
         recvoffsets[i] = recvoffsets[i - 1] + recvcounts[i - 1];

      if( my_rank == 0 )
      {
         solLength += recvoffsets[my_size - 1] + recvcounts[my_size - 1];
         gatheredVec = std::vector<T>(solLength);

         PIPS_MPIgatherv(&gatheredVecLocal[0], mylength, &gatheredVec[0] + firstvec.length(),
            &recvcounts[0], &recvoffsets[0], 0, mpiComm);
      }
      else
      {
        T dummy;
        PIPS_MPIgatherv(&gatheredVecLocal[0], mylength, &dummy, &recvcounts[0], &recvoffsets[0], 0, mpiComm);
      }
   }
   else
   {
      solLength += gatheredVecLocal.size();

      gatheredVec = std::vector<T>(solLength);

      std::copy(gatheredVecLocal.begin(), gatheredVecLocal.end(), gatheredVec.begin() + firstvec.length());
   }

   if( my_rank == 0 )
   {
      std::copy(&firstvec[0], &firstvec[0] + firstvec.length(), &gatheredVec[0]);

      if( vecl && vecl->length() > 0 )
      {
         const SimpleVectorBase<T>& linkvec = dynamic_cast<const SimpleVectorBase<T>&>(*vecl);
         gatheredVec.insert(gatheredVec.end(), &linkvec[0], &linkvec[0] + linkvec.length());
      }
   }
   return gatheredVec;
}

// is root node data of StochVector same on all procs?
template<typename T>
bool StochVectorBase<T>::isRootNodeInSync() const
{
#ifdef HIERARCHICAL
   // TODO adapt for hier approach
//   assert( false && "TODO : implement" );
#endif

   assert(vec);
   assert(mpiComm);

   bool in_sync = true;
   const SimpleVectorBase<T>& vec_simple = dynamic_cast<const SimpleVectorBase<T>&>(*vec);

   /* no need to check if not distributed or not at root node */
   if( !iAmDistrib || parent != nullptr)
      return in_sync;

   int my_rank, world_size;
   MPI_Comm_rank(mpiComm, &my_rank);
   MPI_Comm_size(mpiComm, &world_size);

   /* if there is a linking part we have to chekc it as well */
   const int vec_length = vec_simple.length();
   const int vecl_length = (vecl) ? dynamic_cast<const SimpleVectorBase<T>&>(*vecl).length() : 0;

   const long long count = vec_length + vecl_length;

   assert( count < std::numeric_limits<int>::max());

   /* mpi reduce on vector */
   std::vector<T> sendbuf(count, 0.0);
   std::vector<T> recvbuf(count, 0.0);
   std::copy(vec_simple.elements(), vec_simple.elements() + vec_simple.length(), sendbuf.begin());

   if( vecl )
   {
      const SimpleVectorBase<T>& vecl_simple = dynamic_cast<const SimpleVectorBase<T>&>(*vecl);
      std::copy(vecl_simple.elements(), vecl_simple.elements() + vecl_simple.length(),
            sendbuf.begin() + vec_simple.length());
   }
   PIPS_MPImaxArray(&sendbuf[0], &recvbuf[0], count, mpiComm);

   for( int i = 0; i < count; ++i )
   {
      if( !PIPSisEQ(sendbuf[i], recvbuf[i]) && (sendbuf[i] != recvbuf[i]) && !(std::isnan(sendbuf[i]) && std::isnan(recvbuf[i])) )
      {
//         std::cout << std::setprecision(16);
//         std::cout << sendbuf[i] << " != " << recvbuf[i] << std::endl;

         /* someone else had a higher value here */
         in_sync = false;
      }
   }

   return in_sync;
}

template<typename T>
StochVectorBase<T>* StochVectorBase<T>::raiseBorder( int n_vars, bool linking_part, bool shave_top )
{
   assert( parent == nullptr );
   assert( vec || vecl );

   SimpleVectorBase<T>* vecs;
   if( linking_part )
      vecs = dynamic_cast<SimpleVectorBase<T>*>(vecl);
   else
      vecs = dynamic_cast<SimpleVectorBase<T>*>(vec);

   assert( vecs );
   assert( n_vars <= vecs->length() );

   SimpleVectorBase<T>* border = vecs->shaveBorder(n_vars, shave_top);

   StochVectorBase<T>* top_layer;
   if( shave_top )
       top_layer = new StochVectorBase<T>( border, nullptr, mpiComm );
   else
       top_layer = new StochVectorBase<T>( nullptr, border, mpiComm );

   this->n -= n_vars;

   this->parent = top_layer;
   top_layer->AddChild(this);

   return top_layer;
}

template class StochVectorBase<int>;
template class StochVectorBase<double>;
