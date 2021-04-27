#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "VectorUtilities.h"
#include "StochVector.h"
#include "pipsport.h"
#include "sTree.h"
#include "DistributedQP.hpp"

#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <memory>

#include "StochVector_fwd.h"

template<typename T>
StochVectorBase<T>::StochVectorBase(OoqpVectorBase<T>* first, OoqpVectorBase<T>* last, MPI_Comm mpi_comm)
   : first(first), last(last), mpiComm(mpi_comm ), iAmDistrib(PIPS_MPIgetDistributed(mpiComm) ),
     iAmSpecial( PIPS_MPIiAmSpecial( iAmDistrib, mpiComm) )
{
   assert(first || last );

   if( first )
      this->n += first->length();
   if( last )
      this->n += last->length();
}

template<typename T>
StochVectorBase<T>::StochVectorBase(int n_, MPI_Comm mpiComm_ )
  : OoqpVectorBase<T>(n_), mpiComm(mpiComm_), iAmDistrib( PIPS_MPIgetDistributed(mpiComm) ),
    iAmSpecial( PIPS_MPIiAmSpecial( iAmDistrib, mpiComm) )
{
    first = new SimpleVectorBase<T>(n_);
    last = nullptr;
}

template<typename T>
StochVectorBase<T>::StochVectorBase(int n_, int nl_, MPI_Comm mpiComm_ )
  : mpiComm(mpiComm_), iAmDistrib( PIPS_MPIgetDistributed(mpiComm) ),
    iAmSpecial( PIPS_MPIiAmSpecial( iAmDistrib, mpiComm ) )
{
   if( n_ >= 0)
   {
       first = new SimpleVectorBase<T>(n_);
      this->n += n_;
   }

   if( nl_ >= 0 )
   {
       last = new SimpleVectorBase<T>(nl_);
      this->n += nl_;
   }
}

template<typename T>
void StochVectorBase<T>::AddChild(StochVectorBase<T>* child)
{
   child->parent = this;
   if( child->first->isKindOf(kStochVector) )
      dynamic_cast<StochVectorBase<T>*>(child->first)->parent = this;

   children.push_back(child);

   this->n += child->n;
}

template<typename T>
StochVectorBase<T>::~StochVectorBase()
{
  for (size_t it = 0; it < children.size(); it++)
    delete children[it];

  if( first )
    delete first;

  if( last )
	 delete last;
}

template<typename T>
OoqpVectorBase<T>*
StochVectorBase<T>::clone() const
{
   assert(first || last );
   StochVectorBase<T>* clone = new StochVectorBase<T>(first ? first->clone() : nullptr, last ? last->clone() : nullptr, mpiComm);

   for( size_t it = 0; it < children.size(); it++ )
      clone->AddChild( dynamic_cast<StochVectorBase<T>*>(children[it]->clone()) );

   assert( this->n == clone->n );
   return clone;
}

template<typename T>
OoqpVectorBase<T>* StochVectorBase<T>::cloneFull() const
{
   assert(first || last );
   StochVectorBase<T>* clone = new StochVectorBase<T>(first ? first->cloneFull() : nullptr, last ? last->cloneFull() : nullptr, mpiComm);

   for( size_t it = 0; it < children.size(); it++ )
      clone->AddChild(dynamic_cast<StochVectorBase<T>*>(children[it]->cloneFull()));

   assert( this->n == clone->n );
   return clone;
}

template<typename T>
void StochVectorBase<T>::setNotIndicatedEntriesToVal(T val, const OoqpVectorBase<T>& ind )
{
   const StochVectorBase<T>& ind_vec = dynamic_cast<const StochVectorBase<T>&>(ind);

   assert(this->children.size() == ind_vec.children.size());
   assert((this->first && ind_vec.first) || (this->first == nullptr && ind_vec.first == nullptr) );
   assert((this->last && ind_vec.last) || (this->last == nullptr && ind_vec.last == nullptr) );

   if( this->first )
      this->first->setNotIndicatedEntriesToVal(val, *ind_vec.first);

   if( this->last )
      this->last->setNotIndicatedEntriesToVal(val, *ind_vec.last);

   for( size_t node = 0; node < children.size(); ++node )
      this->children[node]->setNotIndicatedEntriesToVal(val, *ind_vec.children[node] );
}

template<typename T>
void StochVectorBase<T>::jointCopyFrom(const StochVectorBase<T>& vx, const StochVectorBase<T>& vy, const StochVectorBase<T>& vz)
{
   if( vx.first->isKindOf(kStochVector) )
   {
      assert(vy.first->isKindOf(kStochVector) && vz.first->isKindOf(kStochVector) );
      assert( !last );
      assert( children.size() == dynamic_cast<const StochVectorBase&>(*vy.first).children.size() );

      this->jointCopyFrom(dynamic_cast<const StochVectorBase&>(*vx.first), dynamic_cast<const StochVectorBase&>(*vy.first),
                          dynamic_cast<const StochVectorBase&>(*vz.first) );
   }
   else
   {
      assert(first );
      SimpleVectorBase<T>& sv = dynamic_cast<SimpleVectorBase<T>&>(*this->first);

      int n1 = 0;
      int n2 = 0;
      int n3 = 0;
      int n4 = 0;
      int n5 = 0;
      int n6 = 0;

      if( vx.first )
      {
         const SimpleVectorBase<T>& svx = dynamic_cast<const SimpleVectorBase<T>&>(*vx.first);
         n1 = svx.length();
         assert( n1 >= 0 );

         assert( n1 <= sv.length() );
         if( n1 > 0 )
            memcpy(&sv[0], &svx[0], n1 * sizeof(T));
      }

      if( vy.first )
      {
         const SimpleVectorBase<T>& svy = dynamic_cast<const SimpleVectorBase<T>&>(*vy.first);
         n2 = svy.length();
         assert( n2 >= 0 );

         assert( n1 + n2 <= sv.length() );
         if( n2 > 0 )
            memcpy(&sv[n1], &svy[0], n2 * sizeof(T));
      }

      if( vz.first )
      {
         const SimpleVectorBase<T>& svz = dynamic_cast<const SimpleVectorBase<T>&>(*vz.first);
         n3 = svz.length();
         assert( n3 >= 0 );

         assert( n1 + n2 + n3 <= sv.length() );
         if( n3 > 0 )
            memcpy(&sv[n1 + n2], &svz[0], n3 * sizeof(T));
      }

      if( vx.last )
      {
         const SimpleVectorBase<T>& svxl = dynamic_cast<const SimpleVectorBase<T>&>(*vx.last);
         n4 = svxl.length();
         assert( n4 >= 0 );

         assert( n1 + n2 + n3 + n4 <= sv.length() );
         if( n4 > 0 )
            memcpy(&sv[n1 + n2 + n3], &svxl[0], n4 * sizeof(T) );
      }

      if( vy.last )
      {
         const SimpleVectorBase<T>& svyl = dynamic_cast<const SimpleVectorBase<T>&>(*vy.last);
         n5 = svyl.length();
         assert( n5 >= 0 );

         assert( n1 + n2 + n3 + n4 + n5 <= sv.length() );
         if( n5 > 0 )
            memcpy(&sv[n1 + n2 + n3 + n4], &svyl[0], n5 * sizeof(T));
      }

      if( vz.last )
      {
         const SimpleVectorBase<T>& svzl = dynamic_cast<const SimpleVectorBase<T>&>(*vz.last);
         n6 = svzl.length();
         assert( n6 >= 0 );

         assert( n1 + n2 + n3 + n4 + n5 + n6 <= sv.length() );
         if( n6 > 0 )
            memcpy(&sv[n1 + n2 + n3 + n4 + n5], &svzl[0], n6 * sizeof(T));
      }

      assert( n1 + n2 + n3 + n4 + n5 + n6 == sv.length() );

      assert( children.size() == vx.children.size() );
      assert( children.size() == vy.children.size() );
      assert( children.size() == vz.children.size() );

      for(size_t it = 0; it < children.size(); it++)
         children[it]->jointCopyFrom(*vx.children[it], *vy.children[it], *vz.children[it]);
   }
}

template<typename T>
void StochVectorBase<T>::jointCopyTo(StochVectorBase<T>& vx, StochVectorBase<T>& vy, StochVectorBase<T>& vz) const
{
   if( vx.first->isKindOf(kStochVector) )
   {
      assert(vy.first->isKindOf(kStochVector) && vz.first->isKindOf(kStochVector) );
      assert( !last );
      assert( children.size() == dynamic_cast<StochVectorBase&>(*vy.first).children.size() );

      this->jointCopyTo(dynamic_cast<StochVectorBase&>(*vx.first), dynamic_cast<StochVectorBase&>(*vy.first),
                        dynamic_cast<StochVectorBase&>(*vz.first) );
   }
   else
   {
      assert( this->first );
      const SimpleVectorBase<T>& sv  = dynamic_cast<const SimpleVectorBase<T>&>(*this->first);
      int n1 = 0;
      int n2 = 0;
      int n3 = 0;
      int n4 = 0;
      int n5 = 0;
      int n6 = 0;

      if( vx.first )
      {
         SimpleVectorBase<T> &svx = dynamic_cast<SimpleVectorBase<T>&>(*vx.first);
         n1 = svx.length();
         assert(n1 >= 0);

         assert(n1 <= sv.length());
         if( n1 > 0 )
            memcpy(&svx[0], &sv[0], n1 * sizeof(T));
      }

      if( vy.first )
      {
         SimpleVectorBase<T>& svy = dynamic_cast<SimpleVectorBase<T>&>(*vy.first);
         n2 = svy.length();
         assert( n2 >= 0 );

         assert( n1 + n2 <= sv.length() );
         if( n2 > 0 )
            memcpy(&svy[0], &sv[n1], n2 * sizeof(T));
      }

      if( vz.first )
      {
         SimpleVectorBase<T>& svz = dynamic_cast<SimpleVectorBase<T>&>(*vz.first);
         n3 = svz.length();
         assert( n3 >= 0 );

         assert( n1 + n2 + n3 <= sv.length() );
         if( n3 > 0 )
            memcpy(&svz[0], &sv[n1 + n2], n3 * sizeof(T));
      }

      if( vx.last )
      {
         SimpleVectorBase<T>& svxl = dynamic_cast<SimpleVectorBase<T>&>(*vx.last);
         n4 = svxl.length();
         assert( n4 >= 0 );

         assert( n1 + n2 + n3 + n4 <= sv.length() );
         if( n4 > 0 )
            memcpy(&svxl[0], &sv[n1 + n2 + n3], n4 * sizeof(T) );
      }

      if( vy.last )
      {
         SimpleVectorBase<T>& svyl = dynamic_cast<SimpleVectorBase<T>&>(*vy.last);
         n5 = svyl.length();
         assert(n5 >= 0);

         assert( n1 + n2 + n3 + n4 + n5 <= sv.length() );
         if( n5 > 0 )
            memcpy(&svyl[0], &sv[n1 + n2 + n3 + n4], n5 * sizeof(T));
      }

      if( vz.last )
      {
         SimpleVectorBase<T>& svzl = dynamic_cast<SimpleVectorBase<T>&>(*vz.last);
         n6 = svzl.length();
         assert(n6 >= 0);

         assert( n1 + n2 + n3 + n4 + n5 +n6 <= sv.length() );
         if( n6 > 0 )
            memcpy(&svzl[0], &sv[n1 + n2 + n3 + n4 + n5], n6 * sizeof(T));
      }      assert( n1 + n2 + n3 + n4 + n5 <= sv.length() );

      assert( n1 + n2 + n3 + n4 + n5 + n6 == sv.length() );

      assert( children.size() == vx.children.size() );
      assert( children.size() == vy.children.size() );
      assert( children.size() == vz.children.size() );

      for(size_t it = 0; it < children.size(); it++)
         children[it]->jointCopyTo(*vx.children[it], *vy.children[it], *vz.children[it]);
   }
}


template<typename T>
bool StochVectorBase<T>::isKindOf( int kind ) const
{
  return (kind == kStochVector);
}

template<typename T>
void StochVectorBase<T>::scale( T alpha )
{
   if( first )
      first->scale(alpha );

   if( last )
      last->scale(alpha );

   for(size_t it = 0; it < children.size(); it++)
      children[it]->scale( alpha );
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

	if( first )
	{
	   const bool is_zero_tmp = dynamic_cast<SimpleVectorBase<T>&>(*first).isZero();
	   is_zero = is_zero && is_zero_tmp;
	}

	if( last )
   {
	   const bool is_zero_tmp = dynamic_cast<SimpleVectorBase<T>&>(*last).isZero();
	   is_zero = is_zero && is_zero_tmp;
   }

	return is_zero;
}

template<typename T>
void StochVectorBase<T>::setToZero()
{
   setToConstant( T{} );
}

template<typename T>
void StochVectorBase<T>::setToConstant(T c)
{
   if( first )
      first->setToConstant(c);

   if( last )
      last->setToConstant(c);

   for(size_t it = 0; it < children.size(); it++)
      children[it]->setToConstant(c);
}

template<typename T>
void StochVectorBase<T>::copyFrom( const OoqpVectorBase<T>& v_ )
{
   const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

   if( first )
   {
      assert( v.first );
      first->copyFrom(*v.first);
   }
   else
      assert(v.first == nullptr );

   if( last )
   {
      assert( v.last );
      last->copyFrom(*v.last);
   }
   else
      assert(v.last == nullptr );

   assert( children.size() == v.children.size() );

   for(size_t it = 0; it < children.size(); it++)
      children[it]->copyFrom(*v.children[it]);
}

template<typename T>
void StochVectorBase<T>::copyFromAbs(const OoqpVectorBase<T>& v_ )
{
  const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

  if( first )
  {
     assert( v.first );
     first->copyFromAbs(*v.first);
  }
  else
     assert(v.first == nullptr );

  if( last )
  {
     assert(v.last);
     last->copyFromAbs(*v.last);
  }
  else
     assert(v.last == nullptr );

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

   if( first )
      infnrm = std::max(first->infnorm(), infnrm);

   if( last )
      infnrm = std::max(last->infnorm(), infnrm);

   if( iAmDistrib )
      PIPS_MPIgetMaxInPlace(infnrm, mpiComm);

   return infnrm;
}

template<typename T>
double StochVectorBase<T>::twonorm() const
{
  const T scale = this->infnorm();
#ifndef NDEBUG
  if( ! (scale > 0.0) ){
     std::cout << "ERROR : infnorm smaller 0 .. : " << scale << std::endl;
  }
   assert(scale >= 0.0);
#endif

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

  if(first && (iAmSpecial || first->isKindOf(kStochVector)) )
     onenorm += first->onenorm();

  if(iAmSpecial && last )
     onenorm += last->onenorm();

  if( iAmDistrib && !parent )
     PIPS_MPIgetSumInPlace(onenorm, mpiComm);

  return onenorm;
}


template<typename T>
void StochVectorBase<T>::min( T& m, int& index ) const
{
   // index is broken for StochVector
   index = -1;

   m = std::numeric_limits<T>::max();

   if( first )
   {
      T min_tmp = std::numeric_limits<T>::max();
      first->min(min_tmp, index );

      m = std::min( min_tmp, m );
   }

   if( last )
   {
      T min_tmp = std::numeric_limits<T>::max();
      last->min(min_tmp, index);

      m = std::min( min_tmp, m );
   }

   for(size_t it = 0; it < children.size(); it++)
   {
      T min_tmp = std::numeric_limits<T>::max();
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
   if( first )
   {
      T max_tmp = -std::numeric_limits<T>::max();
      first->max(max_tmp, index );

      m = std::max( m, max_tmp );
   }

   if( last )
   {
      T max_tmp = -std::numeric_limits<T>::max();
      last->max(max_tmp, index);

      m = std::max( m, max_tmp );
   }

   for( size_t it = 0; it < children.size(); it++ )
   {
      T max_tmp = -std::numeric_limits<T>::max();
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

   if( first )
   {
      assert( absminvecStoch.first );
      first->absminVecUpdate(*(absminvecStoch.first));
   }
   else
      assert(absminvecStoch.first == nullptr );

   if( last )
   {
      assert( absminvecStoch.last );
      last->absminVecUpdate(*absminvecStoch.last );
   }
   else
      assert(absminvecStoch.last == nullptr );

   assert( absminvecStoch.children.size() == children.size() );

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->absminVecUpdate( *absminvecStoch.children[it] );
}

template<typename T>
void StochVectorBase<T>::absmaxVecUpdate(OoqpVectorBase<T>& absmaxvec) const
{
   StochVectorBase<T>& absmaxvecStoch = dynamic_cast<StochVectorBase<T>&>(absmaxvec);

   if( first )
   {
      assert( absmaxvecStoch.first );
      first->absmaxVecUpdate(*absmaxvecStoch.first);
   }
   else
      assert(absmaxvecStoch.first == nullptr );

   if( last )
   {
      assert( absmaxvecStoch.last );
      last->absmaxVecUpdate(*absmaxvecStoch.last );
   }
   else
      assert(absmaxvecStoch.last == nullptr );

   assert( absmaxvecStoch.children.size() == children.size() );
   for( size_t it = 0; it < children.size(); it++ )
      children[it]->absmaxVecUpdate( *absmaxvecStoch.children[it] );
}

template<typename T>
void StochVectorBase<T>::absmin(T& m) const
{
   T min_tmp = std::numeric_limits<T>::infinity();

   if( first )
   {
      first->absmin(min_tmp );
      m = std::min( m, min_tmp );
   }

   if( last )
   {
      last->absmin(min_tmp );
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

   if( first )
   {
      first->absminNonZero(min_tmp, zero_eps);
      assert( min_tmp > zero_eps );
      m = std::min( min_tmp, m );
   }

   if( last )
   {
      last->absminNonZero(min_tmp, zero_eps);
      assert( min_tmp > zero_eps );
      m = std::min( min_tmp, m );
   }

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->absminNonZero(min_tmp, zero_eps);
      assert( min_tmp > zero_eps );
      m = std::min( min_tmp, m );
   }

   if( iAmDistrib )
      PIPS_MPIgetMinInPlace(m, mpiComm);

   assert(m > zero_eps );
}


template<typename T>
T StochVectorBase<T>::stepbound(const OoqpVectorBase<T> & v_, T maxStep ) const
{
  const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

  T step = 1.0;

  if( first )
  {
     assert( v.first );
     T stepvec = first->stepbound(*v.first, maxStep);
     if( stepvec < step )
        step = stepvec;
  }

  if( last )
  {
     assert(v.last);
     T stepl = last->stepbound(*v.last, maxStep);
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
  if( w.last )
  {
    assert(wstep.last);
    assert(u.last);
    assert(ustep.last);

    step = w.last->findBlocking(*wstep.last, *u.last, *ustep.last, step,
                                w_elt, wstep_elt, u_elt, ustep_elt,
                                first_or_second);
  }

  if( w.first )
  {
     assert(wstep.first);
     assert(u.first);
     assert(ustep.first);

     step = w.first->findBlocking(*wstep.first, *u.first, *ustep.first, step,
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
  if( w.last )
  {
    assert(wstep.last);
    assert(u.last);
    assert(ustep.last);

    w.last->findBlocking_pd(*wstep.last, *u.last, *ustep.last, maxStepPri, maxStepDual,
                            w_elt_p, wstep_elt_p, u_elt_p, ustep_elt_p,
                            w_elt_d, wstep_elt_d, u_elt_d, ustep_elt_d,
                            primalBlocking, dualBlocking);
  }

  if( w.first )
  {
     assert( wstep.first );
     assert( u.first );
     assert( ustep.first );

     w.first->findBlocking_pd(*wstep.first, *u.first, *ustep.first, maxStepPri, maxStepDual,
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

   if( first )
   {
      assert( v.first );
      first->componentMult(*v.first );
   }
   else
      assert(v.first == nullptr );

   if( last )
   {
      assert( v.last );
      last->componentMult(*v.last );
   }
   else
      assert(v.last == nullptr );

   assert( v.children.size() == children.size() );
   for( size_t it = 0; it < children.size(); it++ )
      children[it]->componentMult( *v.children[it] );
}

template<typename T>
void StochVectorBase<T>::componentDiv ( const OoqpVectorBase<T>& v_ )
{
   const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

   if( first )
   {
      assert( v.first );
      first->componentDiv(*v.first );
   }
   else
      assert(v.first == nullptr );

   if( last )
   {
      assert( v.last );
      last->componentDiv(*v.last );
   }
   else
      assert(v.last == nullptr );

   assert( v.children.size() == children.size() );

   for( size_t it = 0; it < children.size(); ++it )
      children[it]->componentDiv( *v.children[it] );
}

template<typename T>
bool StochVectorBase<T>::componentEqual( const OoqpVectorBase<T>& v_, T tol ) const
{
   const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

   bool component_equal = true;

   if( first )
   {
      assert( v.first );
      const bool component_equal_tmp = first->componentEqual(*v.first, tol);
      component_equal = component_equal && component_equal_tmp;
   }
   else
      assert(v.first == nullptr );

   if( !component_equal )
      if( parent == nullptr )
        std::cout << "not equal in root node non-link" << std::endl;

   if( last )
   {
      assert( v.last );
      const bool component_equal_tmp = last->componentEqual(*v.last, tol);
      component_equal = component_equal && component_equal_tmp;
   }
   else
      assert(v.last == nullptr );

   if( !component_equal )
      if( parent == nullptr )
         std::cout << "not equal in root node link" << std::endl;

   assert(v.children.size() == children.size());

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

   if( first )
   {
      const bool not_equal_tmp = first->componentNotEqual(val, tol );
      not_equal = not_equal_tmp && not_equal;
   }

   if( !not_equal )
      if( parent == nullptr )
        std::cout << "equal in root node non-link" << std::endl;

   if( last )
   {
      const bool not_equal_tmp = last->componentNotEqual(val, tol);
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
   if( first )
      first->scalarMult(num );
   if( last )
      last->scalarMult(num );

  for( size_t it = 0; it < children.size(); ++it )
     children[it]->scalarMult( num );
}

template<typename T>
void StochVectorBase<T>::writeToStream( std::ostream& out, int offset ) const
{
   const int my_rank = PIPS_MPIgetRank(mpiComm);
   const int world_size = PIPS_MPIgetSize(mpiComm);

   std::stringstream sout;

   /* print first and keep only for rank 0 */
   for( int i = 0; i < offset; ++i )
      sout << "\t";
   sout << "--first--\n";
   if( first )
      first->writeToStream(sout, offset + 1 );
   for( int i = 0; i < offset; ++i )
      sout << "\t";
   sout << "-------\n";

   if( my_rank != 0 )
   {
      sout.str("");
      sout.clear();
   }

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->writeToStream( sout, offset + 1 );

   const std::string my_row_part = sout.str();
   const std::string full_row = PIPS_MPIallgatherString( my_row_part, mpiComm );

   if( my_rank == 0 )
      out << full_row;

   if( my_rank == 0 && last )
   {
      sout.str("");
      sout.clear();

      for( int i = 0; i < offset; ++i )
         sout << "\t";
      sout << "--last--\n";

      if( last )
         last->writeToStream(sout, offset + 1 );

      for( int i = 0; i < offset; ++i )
         sout << "\t";
      sout << "-------\n";

      out << sout.str();
   }

   if( world_size > 1 )
      MPI_Barrier( mpiComm );
}

template<typename T>
void StochVectorBase<T>::pushAwayFromZero( double tol, double amount, const OoqpVectorBase<T>* select )
{
   const StochVectorBase<T>* selects = select ? dynamic_cast<const StochVectorBase<T>*>(select) : nullptr;

   if( first )
      first->pushAwayFromZero(tol, amount, selects ? selects->first : nullptr );

   if( last )
      last->pushAwayFromZero(tol, amount, selects ? selects->last : nullptr );

   if( selects )
      assert( children.size() == selects->children.size() );

   for( size_t i = 0; i < this->children.size(); ++i )
      this->children[i]->pushAwayFromZero( tol, amount, selects ? selects->children[i] : nullptr );
}

template<typename T>
void StochVectorBase<T>::getSumCountIfSmall( double tol, double& sum_small, int& n_close, const OoqpVectorBase<T>* select ) const
{
   const StochVectorBase<T>* selects = dynamic_cast<const StochVectorBase<T>*>(select);

   if( selects )
      assert( children.size() == selects->children.size() );

   for( size_t i = 0; i < children.size(); ++i )
      this->children[i]->getSumCountIfSmall( tol, sum_small, n_close, selects ? selects->children[i] : nullptr );

   if(first && (iAmSpecial || first->isKindOf(kStochVector)) )
   {
      if( selects )
         assert(selects->first);

      first->getSumCountIfSmall(tol, sum_small, n_close, selects ? selects->first : nullptr );
   }

   if(iAmSpecial && last )
   {
      if( selects )
         assert(selects->last);
      last->getSumCountIfSmall(tol, sum_small, n_close, selects ? selects->last : nullptr );
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
   if( first )
      first->writefToStream(out, format);

   if( last )
      last->writefToStream(out, format);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->writefToStream(out, format);
}

template<typename T>
void StochVectorBase<T>::writeMPSformatRhs( std::ostream& out, int rowType, const OoqpVectorBase<T>* irhs) const
{
   // TODO : will not work with hierarchical data
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
         first->writeMPSformatOnlyRhs(out, rowNameStub, dynamic_cast<const SimpleVectorBase<T>*>(ic->first));
      else
         first->writeMPSformatOnlyRhs(out, rowNameStub, nullptr);
      if(last)
      {
         rowNameStub = " B row_";
         rowNameStub += rt;
         rowNameStub += "_L_";
         if( irhs )
            last->writeMPSformatOnlyRhs(out, rowNameStub, dynamic_cast<const SimpleVectorBase<T>*>(ic->last));
         else
            last->writeMPSformatOnlyRhs(out, rowNameStub, nullptr);
      }
   }
   for(int it=0; it<(int)children.size(); it++)
   {
      std::stringstream sstm;
      sstm << " B row_" << rt << "_" << it << "_";
      std::string rowNameStub = sstm.str();
      if( irhs )
         children[it]->first->writeMPSformatOnlyRhs(out, rowNameStub, dynamic_cast<const SimpleVectorBase<T>*>(ic->children[it]->first));
      else
         children[it]->first->writeMPSformatOnlyRhs(out, rowNameStub, nullptr);
   }
}

template<typename T>
void StochVectorBase<T>::writeMPSformatBounds(std::ostream& out, const OoqpVectorBase<T>* ix, bool upperBound) const
{
   // TODO : will not work with hierarchical data
   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);

   const StochVectorBase<T>* ixStoch = dynamic_cast<const StochVectorBase<T>*>(ix);

   if( myRank==0 )
   {
      std::string varNameStub = "var_L_";
      first->writeMPSformatBoundsWithVar(out, varNameStub, (ixStoch->first), upperBound);
   }
   for(int it=0; it<(int)children.size(); it++)
   {
      std::stringstream sstm2;
      sstm2 << "var_" << it << "_";
      std::string varNameStub = sstm2.str();
      children[it]->first->writeMPSformatBoundsWithVar(out, varNameStub, (ixStoch->children[it]->first), upperBound);
   }
}

/** this += alpha * x */
template<typename T>
void StochVectorBase<T>::axpy ( T alpha, const OoqpVectorBase<T>& x_ )
{
   const StochVectorBase<T>& x = dynamic_cast<const StochVectorBase<T>&>(x_);

   if( alpha == 0.0 )
      return;

   if( first )
   {
      assert( x.first );
      first->axpy(alpha, *x.first );
   }
   else
      assert(x.first == nullptr );

   if( last )
   {
      assert( x.last );
      last->axpy(alpha, *x.last );
   }
   else
      assert(x.last == nullptr );

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

  if( first )
  {
     assert( x.first );
     assert( z.first );
     first->axzpy(alpha, *x.first, *z.first );
  }
  else
  {
     assert(x.first == nullptr );
     assert(z.first == nullptr );
  }

  if( last )
  {
     assert( x.last );
     assert( z.last );
     last->axzpy(alpha, *x.last, *z.last );
  }
  else
  {
     assert(x.last == nullptr );
     assert(z.last == nullptr );
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

   if( first )
   {
      assert( x.first );
      assert( z.first );
      first->axdzpy(alpha, *x.first, *z.first );
   }
   else
   {
      assert(x.first == nullptr );
      assert(z.first == nullptr );
   }

   if( last )
   {
      assert(x.last);
      assert(z.last);
      last->axdzpy(alpha, *x.last, *z.last );
   }
   else
   {
      assert(x.last == nullptr );
      assert(z.last == nullptr );
   }

   assert( x.children.size() == children.size() );
   assert( z.children.size() == children.size() );

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->axdzpy( alpha, *x.children[it], *z.children[it] );
}


template<typename T>
void StochVectorBase<T>::addConstant( T c )
{
   if( first )
      first->addConstant(c);

   if( last )
      last->addConstant(c);

  for( size_t it = 0; it < children.size(); it++ )
     children[it]->addConstant(c);
}


template<typename T>
void StochVectorBase<T>::gondzioProjection( T rmin, T rmax )
{
   if( first )
      first->gondzioProjection(rmin, rmax );

   if( last )
      last->gondzioProjection(rmin, rmax );

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

   if(first && (iAmSpecial || first->isKindOf(kStochVector)) )
   {
      assert( v.first );
      dot_product += first->dotProductWith(*v.first );
   }
   else if( !first )
      assert(v.first == nullptr );

   if(iAmSpecial && last )
   {
      assert( v.last );
      dot_product += last->dotProductWith(*v.last );
   }
   else if( !last )
      assert(v.last == nullptr );

   if( iAmDistrib && parent == nullptr )
      PIPS_MPIgetSumInPlace(dot_product, mpiComm);

   return dot_product;
}

template<typename T>
T StochVectorBase<T>::dotProductSelf(T scaleFactor) const
{
#ifndef NDEBUG
   if( scaleFactor < 0.0 ){
      std::cout << "ERROR : infnorm smaller 0 .. : " << scaleFactor << std::endl;
   }
   assert(scaleFactor >= 0.0);
#endif

   T dot_product = 0.0;

   for( size_t it = 0; it < children.size(); it++ )
      dot_product += children[it]->dotProductSelf( scaleFactor );

   if(first && (iAmSpecial || first->isKindOf(kStochVector)) )
      dot_product += first->dotProductSelf(scaleFactor );

   if(iAmSpecial && last )
      dot_product += last->dotProductSelf(scaleFactor );

   if( iAmDistrib && !parent )
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

   assert( children.size() == mystep.children.size() );
   assert( children.size() == yvec.children.size() );
   assert( children.size() == ystep.children.size() );

   for( size_t it = 0; it < children.size(); it++ )
      dot_product += children[it]->shiftedDotProductWith( alpha, *mystep.children[it], *yvec.children[it], beta, *ystep.children[it] );

   if(first && (iAmSpecial || first->isKindOf(kStochVector)) )
   {
      assert( mystep.first );
      assert( yvec.first );
      assert( ystep.first );
      dot_product += first->shiftedDotProductWith(alpha, *mystep.first, *yvec.first, beta, *ystep.first );
   }
   else if( !first )
   {
      assert(mystep.first == nullptr );
      assert(yvec.first == nullptr );
      assert(ystep.first == nullptr );
   }

   if(iAmSpecial && last )
   {
      assert( mystep.last );
      assert( yvec.last );
      assert( ystep.last );
      dot_product += last->shiftedDotProductWith(alpha, *mystep.last, *yvec.last, beta, *ystep.last );
   }
   else if( !last )
   {
      assert(mystep.last == nullptr );
      assert(yvec.last == nullptr );
      assert(ystep.last == nullptr );
   }

   if( iAmDistrib && parent == nullptr )
      PIPS_MPIgetSumInPlace(dot_product, mpiComm);

   return dot_product;
}

template<typename T>
void StochVectorBase<T>::negate()
{
   if( first )
      first->negate();
   if( last )
      last->negate();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->negate();
}

template<typename T>
void StochVectorBase<T>::invert()
{
   if( first )
      first->invert();

   if( last )
      last->invert();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->invert();
}

template<typename T>
void StochVectorBase<T>::invertSave(T zeroReplacementVal)
{
   if( first )
      first->invertSave(zeroReplacementVal );

   if( last )
      last->invertSave(zeroReplacementVal );

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->invertSave( zeroReplacementVal );
}

template<typename T>
void StochVectorBase<T>::applySqrt()
{
   if( first )
      first->applySqrt();

   if( last )
      last->applySqrt();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->applySqrt();
}

template<typename T>
void StochVectorBase<T>::roundToPow2()
{
   if( first )
      first->roundToPow2();

   if( last )
      last->roundToPow2();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->roundToPow2();
}

template<typename T>
bool StochVectorBase<T>::allPositive() const
{
   bool all_pos = true;

   if( first )
   {
      const bool all_pos_tmp = first->allPositive();
      all_pos = all_pos && all_pos_tmp;
   }

   if( last )
   {
      const bool all_pos_tmp = last->allPositive();
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
bool StochVectorBase<T>::allOf( const std::function<bool(const T&)>& pred ) const
{
   bool all = true;

   if( first )
   {
      const bool all_vec = first->allOf(pred );
      all = all && all_vec;
   }

   if( last )
   {
      const bool all_vecl = last->allOf(pred );
      all = all && all_vecl;
   }

   for( const auto& child : children )
   {
      const bool all_child = child->allOf( pred );
      all = all && all_child;
   }

   return all;
}


template<typename T>
bool StochVectorBase<T>::matchesNonZeroPattern( const OoqpVectorBase<T>& select_ ) const
{
  const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);

  bool match = true;

  if( first )
  {
     assert( select.first );
     const bool match_tmp = first->matchesNonZeroPattern(*select.first );
     match = match && match_tmp;
  }
  else
     assert(select.first == nullptr );

  if( last )
  {
     assert(select.last);
     const bool match_tmp = last->matchesNonZeroPattern(*select.last );
     match = match && match_tmp;
  }
  else
     assert(select.last == nullptr );

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

   if( first )
   {
      assert( select.first );
      first->selectNonZeros(*select.first );
   }
   else
      assert(select.first == nullptr );
   if( last )
   {
      assert( select.last );
      last->selectNonZeros(*select.last );
   }
   else
      assert(select.last == nullptr );

   assert(children.size() == select.children.size());

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->selectNonZeros( *select.children[it] );
}

template<typename T>
void StochVectorBase<T>::selectPositive()
{
   if( first )
      first->selectPositive();

   if( last )
      last->selectPositive();

   for(size_t it = 0; it < children.size(); it++)
     children[it]->selectPositive();
}

template<typename T>
void StochVectorBase<T>::selectNegative()
{
   if( first )
      first->selectNegative();

   if( last )
      last->selectNegative();

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

  if(first && (iAmSpecial || first->isKindOf(kStochVector)) )
     nnz += first->numberOfNonzeros();

  if(iAmSpecial && last )
     nnz += last->numberOfNonzeros();

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

   if( first )
   {
      assert( select.first );
      first->addSomeConstants(c, *select.first);
   }
   else
      assert(select.first == nullptr );

   if( last )
   {
      assert(select.last);
      last->addSomeConstants(c, *select.last);
   }
   else
      assert(select.last == nullptr );

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->addSomeConstants(c, *select.children[it]);
}

template<typename T>
void StochVectorBase<T>::axdzpy( T alpha, const OoqpVectorBase<T>& x_,
		       const OoqpVectorBase<T>& z_, const OoqpVectorBase<T>& select_ )
{
   const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);
   const StochVectorBase<T>& x = dynamic_cast<const StochVectorBase<T>&>(x_);
   const StochVectorBase<T>& z = dynamic_cast<const StochVectorBase<T>&>(z_);


   if( first )
   {
      assert( x.first );
      assert( z.first );
      assert( select.first );
      first->axdzpy(alpha, *x.first, *z.first, *select.first);
   }
   else
   {
      assert(x.first == nullptr );
      assert(z.first == nullptr );
      assert(select.first == nullptr );
   }

   if( last )
   {
      assert( x.last );
      assert( z.last );
      assert( select.last );
      last->axdzpy(alpha, *x.last, *z.last, *select.last);
   }
   else
   {
      assert(x.last == nullptr );
      assert(z.last == nullptr );
      assert(select.last == nullptr );
   }

   assert(children.size() == select.children.size());
   assert(children.size() == x.children.size());
   assert(children.size() == z.children.size());

   for(size_t it = 0; it < children.size(); it++)
      children[it]->axdzpy(alpha, *x.children[it], *z.children[it], *select.children[it]);
}

template<typename T>
bool StochVectorBase<T>::somePositive( const OoqpVectorBase<T>& select_ ) const
{
   const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);

   bool some_positive = true;;

   assert( children.size() == select.children.size() );

   for(size_t it = 0; it < children.size(); it++)
   {
      const bool some_pos_tmp = children[it]->somePositive( *select.children[it] );
      some_positive = some_positive && some_pos_tmp;
   }

   if( iAmDistrib )
      PIPS_MPIgetLogicAndInPlace( some_positive, mpiComm );

   if( first )
   {
      assert( select.first );
      const bool some_pos_tmp = first->somePositive(*select.first);
      some_positive = some_positive && some_pos_tmp;
   }
   else
      assert(select.first == nullptr );

   if( last )
   {
      assert(select.last);
      const bool some_pos_tmp = last->somePositive(*select.last);
      some_positive = some_positive && some_pos_tmp;
   }
   else
      assert(select.last == nullptr );

   return some_positive;
}

template<typename T>
void StochVectorBase<T>::divideSome( const OoqpVectorBase<T>& div_, const OoqpVectorBase<T>& select_ )
{
   const StochVectorBase<T> &div = dynamic_cast<const StochVectorBase<T>&>(div_);
   const StochVectorBase<T> &select = dynamic_cast<const StochVectorBase<T>&>(select_);

   if( first )
   {
      assert(div.first);
      assert(select.first);
      first->divideSome(*div.first, *select.first);
   }
   else
   {
      assert(div.first == nullptr );
      assert(select.first == nullptr );
   }

   if( last )
   {
      assert(div.last);
      assert(select.last);
      last->divideSome(*div.last, *select.last);
   }
   else
   {
      assert(div.last == nullptr );
      assert(select.last == nullptr );
   }

   assert(children.size() == div.children.size());
   assert(children.size() == select.children.size());

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->divideSome(*div.children[it], *select.children[it]);
}

template<typename T>
void StochVectorBase<T>::removeEntries( const OoqpVectorBase<int>& select_ )
{
   const StochVectorBase<int>& select = dynamic_cast<const StochVectorBase<int>&>(select_);

   this->n = 0;

   if( first )
   {
      assert( select.first );
      first->removeEntries(*select.first);
      this->n = first->length();
   }
   else
      assert(select.first == nullptr );

   if( last )
   {
      assert( select.last );
      last->removeEntries(*select.last);
      this->n += last->length();
   }
   else
      assert(select.last == nullptr );

   assert(children.size() == select.children.size());
   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->removeEntries(*select.children[it]);
      this->n += children[it]->length();
   }
}

/* uses PIPSisZero() as check instead of == 0.0 */
template<typename T>
int StochVectorBase<T>::getNnzs() const
{
   int non_zeros = 0;

   for( size_t it = 0; it < children.size(); it++ )
      non_zeros += children[it]->getNnzs();

   if(first && (iAmSpecial || first->isKindOf(kStochVector)) )
      non_zeros += first->getNnzs();

   if(iAmSpecial && last )
      non_zeros += last->getNnzs();

   if( iAmDistrib && parent == nullptr )
      PIPS_MPIgetSumInPlace(non_zeros, mpiComm);

   return non_zeros;
}

template<typename T>
void StochVectorBase<T>::recomputeSize()
{
   if(first && first->isKindOf(kStochVector) )
      dynamic_cast<StochVector*>(first)->recomputeSize();
   if(last && last->isKindOf(kStochVector) )
      dynamic_cast<StochVector*>(last)->recomputeSize();

   this->n = 0;
   if( first )
      this->n += first->length();
   if( last )
      this->n += last->length();

   for( auto& child : children )
   {
      child->recomputeSize();
      this->n += child->length();
   }
}

template<typename T>
void StochVectorBase<T>::permuteVec0Entries(const std::vector<unsigned int>& permvec)
{
   if( first )
      dynamic_cast<SimpleVectorBase<T>*>(first)->permuteEntries(permvec);
}

template<typename T>
void StochVectorBase<T>::permuteLinkingEntries(const std::vector<unsigned int>& permvec)
{
   if( last )
      dynamic_cast<SimpleVectorBase<T>*>(last)->permuteEntries(permvec);
}

template<typename T>
std::vector<T> StochVectorBase<T>::gatherStochVector() const
{
   if( pips_options::getBoolParameter( "HIERARCHICAL" ) )
   {
      // TODO adapt for hier approach
      // assert( false && "TODO : implement" );
   }

   const SimpleVectorBase<T>& firstvec = dynamic_cast<const SimpleVectorBase<T>&>(*first);
   const size_t nChildren = children.size();

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   const int my_size = PIPS_MPIgetSize(mpiComm);

   std::vector<T> gatheredVecLocal(0);

   for( size_t i = 0; i < nChildren; ++i )
   {
      const SimpleVectorBase<T>& vec = dynamic_cast<const SimpleVectorBase<T>&>(*children[i]->first);

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

      if(last && last->length() > 0 )
      {
         const SimpleVectorBase<T>& linkvec = dynamic_cast<const SimpleVectorBase<T>&>(*last);
         gatheredVec.insert(gatheredVec.end(), &linkvec[0], &linkvec[0] + linkvec.length());
      }
   }
   return gatheredVec;
}

// is root node data of StochVector same on all procs?
template<typename T>
bool StochVectorBase<T>::isRootNodeInSync() const
{
   if( pips_options::getBoolParameter( "HIERARCHICAL" ) )
   {
      // TODO adapt for hier approach
      // assert( false && "TODO : implement" );
   }

   assert(first);
   assert(mpiComm);

   bool in_sync = true;
   const SimpleVectorBase<T>& vec_simple = dynamic_cast<const SimpleVectorBase<T>&>(*first);

   /* no need to check if not distributed or not at root node */
   if( !iAmDistrib || parent != nullptr)
      return in_sync;

   int my_rank, world_size;
   MPI_Comm_rank(mpiComm, &my_rank);
   MPI_Comm_size(mpiComm, &world_size);

   /* if there is a linking part we have to chekc it as well */
   const int vec_length = vec_simple.length();
   const int vecl_length = (last) ? dynamic_cast<const SimpleVectorBase<T>&>(*last).length() : 0;

   const long long count = vec_length + vecl_length;

   assert( count < std::numeric_limits<int>::max());

   /* mpi reduce on vector */
   std::vector<T> sendbuf(count, 0.0);
   std::vector<T> recvbuf(count, 0.0);
   std::copy(vec_simple.elements(), vec_simple.elements() + vec_simple.length(), sendbuf.begin());

   if( last )
   {
      const SimpleVectorBase<T>& vecl_simple = dynamic_cast<const SimpleVectorBase<T>&>(*last);
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
void StochVectorBase<T>::split( const std::vector<unsigned int>& map_blocks_children, const std::vector<MPI_Comm>& child_comms,
      const std::vector<int>& twolinks_start_in_block, int
#ifndef NDEBUG
      n_links_in_root
#endif
      )
{
   const unsigned int n_curr_children = children.size();
   assert( n_curr_children == map_blocks_children.size() );

   if( last )
   {
      assert( n_curr_children == twolinks_start_in_block.size() );
      assert(std::accumulate( twolinks_start_in_block.begin(), twolinks_start_in_block.end(), 0 ) <= last->length() );
      assert( twolinks_start_in_block.back() == 0 );
   }
   else
      assert( twolinks_start_in_block.empty() );

   const unsigned int n_new_children = getNDistinctValues(map_blocks_children);
   std::vector<StochVectorBase<T>*> new_children(n_new_children);

   SimpleVectorBase<T>* vecl_leftover = dynamic_cast<SimpleVectorBase<T>*>(last);

   unsigned int begin_curr_child_blocks{0};
   unsigned int end_curr_child_blocks{0};

   for( unsigned int i = 0; i < n_new_children; ++i )
   {
      while( end_curr_child_blocks != (n_curr_children - 1) &&
            map_blocks_children[end_curr_child_blocks] == map_blocks_children[end_curr_child_blocks + 1] )
         ++end_curr_child_blocks;

      const int n_links_for_child = last ? std::accumulate(twolinks_start_in_block.begin() + begin_curr_child_blocks,
            twolinks_start_in_block.begin() + end_curr_child_blocks, 0 ) : -1;
      const unsigned int n_blocks_for_child = end_curr_child_blocks - begin_curr_child_blocks + 1;

      SimpleVectorBase<T>* vecl_child = last ? vecl_leftover->shaveBorder(n_links_for_child, true) : nullptr;

      if( child_comms[i] == MPI_COMM_NULL )
      {
         delete vecl_child;
         vecl_child = nullptr;
      }

      StochVectorBase<T>* vec = (child_comms[i] == MPI_COMM_NULL) ? nullptr :
            new StochVectorBase<T>( new SimpleVectorBase<T>(0), vecl_child, child_comms[i] );

      for( unsigned int j = 0; j < n_blocks_for_child; ++j )
      {
         StochVectorBase<T>* child = children.front();
         assert( !child->last );
         children.erase(children.begin());

         if( child_comms[i] == MPI_COMM_NULL )
            assert( child->mpiComm == MPI_COMM_NULL );

         if( child_comms[i] != MPI_COMM_NULL )
            vec->AddChild(child);
//         else
//            delete child;
      }

      /* create child holding the new StochVector as it's first part */
      new_children[i] = (child_comms[i] == MPI_COMM_NULL) ? new StochDummyVectorBase<T>() : new StochVectorBase<T>( vec, nullptr, child_comms[i] );
      if( vec )
         vec->parent = new_children[i];

      ++end_curr_child_blocks;
      begin_curr_child_blocks = end_curr_child_blocks;
   }

   if( last )
      assert( vecl_leftover->length() == n_links_in_root );
   assert( children.empty() );

   this->n = 0;
   if( first )
      this->n += first->length();
   if( last )
      this->n += last->length();

   for( auto& child : new_children )
      AddChild( child );
}

template<typename T>
StochVectorBase<T>* StochVectorBase<T>::raiseBorder( int n_vars, bool linking_part, bool shave_top )
{
   assert( parent == nullptr );
   assert(first || last );

   SimpleVectorBase<T>* vecs;
   if( linking_part )
      vecs = dynamic_cast<SimpleVectorBase<T>*>(last);
   else
      vecs = dynamic_cast<SimpleVectorBase<T>*>(first);

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

   StochVectorBase<T>* me = this;
   IotrAddRef(&me);

   return top_layer;
}

template<typename T>
void StochVectorBase<T>::collapseFromHierarchical( const DistributedQP& data_hier, const sTree& tree_hier, VectorType type, bool empty_vec )
{
   SimpleVectorBase<T>* new_vec = new SimpleVectorBase<T>();
   SimpleVectorBase<T>* new_vecl{};

   if( (tree_hier.getMYL() > 0 && type == VectorType::DUAL_Y) ||
         (tree_hier.getMYL() > 0 && type == VectorType::DUAL_Z) )
      new_vecl = new SimpleVectorBase<T>();

   assert( tree_hier.nChildren() == 1 );
   assert( children.size() == 1 );
   assert( data_hier.children.size() == 1 );

   /* hierarchical top */
   std::vector<StochVectorBase<T>*> new_children;
   children[0]->appendHierarchicalToThis( new_vec, new_vecl, new_children, *tree_hier.getChildren()[0], *data_hier.children[0], type, empty_vec );

   delete children[0];
   children.clear();

   if(first && !empty_vec )
      new_vec->appendToFront( dynamic_cast<SimpleVectorBase<T>&>(*first) );
   if(last && !empty_vec )
      new_vecl->appendToBack( dynamic_cast<SimpleVectorBase<T>&>(*last) );

   this->n = new_vec->length();
   if( new_vecl )
      this->n += new_vecl->length();
   for( auto child : new_children )
   {
      assert( child->children.empty() );
      this->AddChild( child );
   }

   delete first;
   delete last;
    first = new_vec;
    last = new_vecl;

   if( first )
      PIPS_MPImaxArrayInPlace(dynamic_cast<SimpleVector&>(*first).elements(), first->length() );

   if( last )
      PIPS_MPImaxArrayInPlace(dynamic_cast<SimpleVector&>(*last).elements(), last->length() );
}

template<typename T>
void StochVectorBase<T>::appendHierarchicalToThis( SimpleVectorBase<T>* new_vec, SimpleVectorBase<T>* new_vecl,
      std::vector<StochVectorBase<T>*>& new_children, const sTree& tree_hier, const DistributedQP& data_hier, VectorType type, bool empty_vec )
{
   assert( children.size() == tree_hier.nChildren() );
   assert( children.size() == data_hier.children.size() );

   for( size_t i = 0; i < children.size(); ++i )
   {
      StochVectorBase<T>& child = *children[i];
      const sTree* sub_root = tree_hier.getChildren()[i]->getSubRoot();

      // not a leaf
      if( sub_root )
      {
         if( child.isKindOf(kStochDummy) )
            child.appendHierarchicalToThis( new_vec, new_vecl, new_children, *sub_root, *data_hier.children[i], type, empty_vec );
         else
         {
            assert( child.first->isKindOf(kStochVector) );
            assert( !child.last);
            dynamic_cast<StochVectorBase<T>&>(*child.first).appendHierarchicalToThis(new_vec, new_vecl, new_children, *sub_root, *data_hier.children[i], type, empty_vec);
         }
      }
      else
      {
         // a leaf gets added to new_children
         assert( child.children.empty() );
         new_children.insert(new_children.end(), &child);
         children[i] = nullptr;
      }
   }

   if(first && !empty_vec)
      new_vec->appendToFront( dynamic_cast<SimpleVectorBase<T>&>(*first) );

   if(last && !empty_vec )
      new_vecl->appendToBack( dynamic_cast<SimpleVectorBase<T>&>(*last) );

   if( !empty_vec && !data_hier.isHierarchyInnerRoot() && data_hier.isHierarchyInnerLeaf() && type != VectorType::PRIMAL )
   {
      PERMUTATION link_vec_perm;

      if( type == VectorType::DUAL_Y )
         link_vec_perm = data_hier.getLinkConsEqPermInv();
      else if( type == VectorType::DUAL_Z )
         link_vec_perm = data_hier.getLinkConsIneqPermInv();

      assert( static_cast<unsigned int>(new_vecl->length()) >= link_vec_perm.size() );
      SimpleVectorBase<T> new_vecl_curr_part = SimpleVectorBase<T>(new_vecl->elements() + new_vecl->length() - link_vec_perm.size(), link_vec_perm.size() );
      new_vecl_curr_part.permuteEntries(link_vec_perm);
   }
}

template<typename T>
void StochDummyVectorBase<T>::appendHierarchicalToThis( SimpleVectorBase<T>*, SimpleVectorBase<T>* new_vecl,
      std::vector<StochVectorBase<T>*>& new_children, const sTree& tree_hier, const DistributedQP&, VectorType type, bool empty_vec )
{
   assert( tree_hier.getCommWorkers() == MPI_COMM_NULL );
   const int n_dummies = tree_hier.nChildren();
   /* insert the children this dummy is representing */
   for ( int i = 0; i < n_dummies; ++i )
      new_children.insert(new_children.end(), new StochDummyVectorBase<T>() );

   if ( type == VectorType::DUAL_Y && !empty_vec )
      new_vecl->appendToBack(tree_hier.getMYL(), -std::numeric_limits<T>::infinity() );
   else if ( type == VectorType::DUAL_Z && !empty_vec )
      new_vecl->appendToBack(tree_hier.getMZL(), -std::numeric_limits<T>::infinity() );
}


template<typename T>
OoqpVectorBase<T>* StochVectorBase<T>::getLinkingVecNotHierarchicalTop() const
{
   const StochVectorBase<T>* curr_par = parent;
   if( curr_par == nullptr )
   {
      /* we are the top */
      assert(first);
      return first;
   }

   while( curr_par->parent != nullptr )
      curr_par = curr_par->parent;

   if( curr_par->children.size() == 1 && curr_par->children[0]->children.size() != 0 )
   {
      /* the current parent is the hierarchical top */
      assert( curr_par->children[0]->first );
      return curr_par->children[0]->first;
   }
   else
   {
      /* the current parent is a normal top */
      assert( curr_par->first );
      return curr_par->first;
   }
}

template<typename T>
void StochVectorBase<T>::pushSmallComplementarityPairs( OoqpVectorBase<T>& other_vec_in, const OoqpVectorBase<T>& select_in, double tol_this, double tol_other, double tol_pairs )
{
   StochVectorBase<T>& other_vec = dynamic_cast<StochVectorBase<T>&>(other_vec_in);
   const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_in);

   assert(children.size() == other_vec.children.size());
   assert(children.size() == select.children.size());

   if( first )
   {
      assert( other_vec.first );
      assert( select.first );

      first->pushSmallComplementarityPairs(*other_vec.first, *select.first, tol_this, tol_other, tol_pairs );
   }

   if( last )
   {
      assert( other_vec.last );
      assert( select.last );

      last->pushSmallComplementarityPairs(*other_vec.last, *select.last, tol_this, tol_other, tol_pairs );
   }

   for( size_t i = 0; i < children.size(); ++i )
   {
      assert( children[i] );
      assert( other_vec.children[i] );
      assert( select.children[i] );

      children[i]->pushSmallComplementarityPairs( *other_vec.children[i], *select.children[i], tol_this, tol_other, tol_pairs );
   }

}

template class StochVectorBase<int>;
template class StochVectorBase<double>;
