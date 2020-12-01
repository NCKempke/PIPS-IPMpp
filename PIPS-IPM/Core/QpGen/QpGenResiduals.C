/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGenResiduals.h"
#include "QpGenVars.h"
#include "QpGenData.h"

#include "OoqpVector.h"
#include "LinearAlgebraPackage.h"

#include "pipsdef.h"

#include <iostream>
#include <fstream>

#include "mpi.h"

QpGenResiduals::QpGenResiduals( LinearAlgebraPackage * la,
				long long nx_, long long my_, long long mz_,
				OoqpVector * ixlow_in, OoqpVector * ixupp_in,
				OoqpVector * iclow_in, OoqpVector * icupp_in )
{
  assert(false && "Cannot be used with StochLinearAlgebra");
  nx = nx_;
  my = my_;
  mz = mz_;

  SpReferTo( ixlow, ixlow_in );
  nxlow = ixlow->numberOfNonzeros();

  SpReferTo( ixupp, ixupp_in );
  nxupp = ixupp->numberOfNonzeros();

  SpReferTo( iclow, iclow_in );
  mclow = iclow->numberOfNonzeros();

  SpReferTo( icupp, icupp_in );
  mcupp = icupp->numberOfNonzeros();

  rQ = OoqpVectorHandle( la->newVector( nx ) );
  rA = OoqpVectorHandle( la->newVector( my ) );
  rC = OoqpVectorHandle( la->newVector( mz ) );

  rz = OoqpVectorHandle( la->newVector( mz ) );
  if ( mclow > 0 ) {
    rt      = OoqpVectorHandle( la->newVector( mz ) );
    rlambda = OoqpVectorHandle( la->newVector( mz ) );
  }
  if ( mcupp > 0 ) {
    ru     = OoqpVectorHandle( la->newVector( mz ) );
    rpi    = OoqpVectorHandle( la->newVector( mz ) );
  }
  if( nxlow > 0 ) {
    rv     = OoqpVectorHandle( la->newVector( nx ) );
    rgamma = OoqpVectorHandle( la->newVector( nx ) );
  }
  if( nxupp > 0 ) {
    rw   = OoqpVectorHandle( la->newVector( nx ) );
    rphi = OoqpVectorHandle( la->newVector( nx ) );
  }
}

QpGenResiduals::QpGenResiduals( const QpGenResiduals& res) : Residuals(res)
{
  nx = res.nx;
  my = res.my;
  mz = res.mz;

  ixlow = OoqpVectorHandle( res.ixlow->cloneFull() );
  nxlow = res.nxlow;

  ixupp = OoqpVectorHandle( res.ixupp->cloneFull() );
  nxupp = res.nxupp;

  iclow = OoqpVectorHandle( res.iclow->cloneFull() );
  mclow = res.mclow;

  icupp = OoqpVectorHandle( res.icupp->cloneFull() );
  mcupp = res.mcupp;

  rQ = OoqpVectorHandle( res.rQ->cloneFull() );
  rA = OoqpVectorHandle( res.rA->cloneFull() );
  rC = OoqpVectorHandle( res.rC->cloneFull() );

  rz = OoqpVectorHandle( res.rz->cloneFull() );

  rt = OoqpVectorHandle( res.rt->cloneFull() );
  rlambda = OoqpVectorHandle( res.rlambda->cloneFull() );

  ru = OoqpVectorHandle( res.ru->cloneFull() );
  rpi = OoqpVectorHandle( res.rpi->cloneFull() );

  rv = OoqpVectorHandle( res.rv->cloneFull() );
  rgamma = OoqpVectorHandle( res.rgamma->cloneFull() );
  
  rw = OoqpVectorHandle( res.rw->cloneFull() );
  rphi = OoqpVectorHandle( res.rphi->cloneFull() );
}

double updateNormAndPrint( double norm, const OoqpVector& vec, bool print, std::string&& name )
{
   const double infnorm = vec.infnorm();

   if( print )
   {
      const double twonorm = vec.twonorm();

      if( 0 == PIPS_MPIgetRank() )
         std::cout << name << " infnorm = " << infnorm << " | twonorm = " << twonorm << "\n";
   }

   return std::max( norm, infnorm );
}

void QpGenResiduals::calcresids(Data *prob_in, Variables *vars_in, bool print_resids)
{
#ifdef TIMING
  print_resids = true;
#endif
  const int myRank = PIPS_MPIgetRank();
  QpGenVars * vars = (QpGenVars *) vars_in;
  QpGenData * prob = (QpGenData *) prob_in;

  double norm = 0.0, gap = 0.0;
 
  /*** rQ = Qx + g - A^T y - C^T z - gamma + phi ***/
  prob->getg( *rQ );
  prob->Qmult( 1.0, *rQ,  1.0, *vars->x );

  // contribution calculate x^T (g + Qx) to duality gap */
  gap += rQ->dotProductWith(*vars->x); 

  prob->ATransmult( 1.0, *rQ, -1.0, *vars->y );
  prob->CTransmult( 1.0, *rQ, -1.0, *vars->z );

  vars->gamma->selectNonZeros( *ixlow );
  vars->phi->selectNonZeros( *ixupp );
  if( nxlow > 0 ) rQ->axpy( -1.0, *vars->gamma );
  if( nxupp > 0 ) rQ->axpy(  1.0, *vars->phi );

  norm = updateNormAndPrint( norm, *rQ, print_resids, "rQ");

  /*** rA = Ax - b ***/
  prob->getbA( *rA );
  prob->Amult( -1.0, *rA, 1.0, *vars->x );

  norm = updateNormAndPrint( norm, *rA, print_resids, "rA");

  // contribution -b^T y to duality gap
  gap -= prob->bA->dotProductWith(*vars->y);

  /*** rC = Cx - s ***/
  rC->copyFrom( *vars->s );
  prob->Cmult( -1.0, *rC, 1.0, *vars->x );

  norm = updateNormAndPrint( norm, *rC, print_resids, "rC");

  /*** rz = z - lambda + pi ***/
  rz->copyFrom( *vars->z );

  if( mclow > 0 )
    rz->axpy( -1.0, *vars->lambda );
  if( mcupp > 0 )
    rz->axpy(  1.0, *vars->pi );


  norm = updateNormAndPrint( norm, *rz, print_resids, "rz");

  /*** rt = s - d - t ***/
  if( mclow > 0 )
  {
     rt->copyFrom( *vars->s );
     rt->axpy( -1.0, prob->slowerBound() );
     rt->selectNonZeros( *iclow );
     rt->axpy( -1.0, *vars->t );

     norm = updateNormAndPrint( norm, *rt, print_resids, "rt");

     // contribution - d^T lambda to duality gap
     gap -= prob->bl->dotProductWith(*vars->lambda);
  }


  if( mcupp > 0 )
  {
     /*** ru = s - f + u ***/
     ru->copyFrom( *vars->s );
     ru->axpy( -1.0, prob->supperBound() );
     ru->selectNonZeros( *icupp );
     ru->axpy( 1.0, *vars->u );

     norm = updateNormAndPrint( norm, *ru, print_resids, "ru");

     // contribution - f^T pi to duality gap
     gap += prob->bu->dotProductWith(*vars->pi);
  }

  if( nxlow > 0 )
  {
    /*** rv = x - lx - v ***/
    rv->copyFrom( *vars->x );
    rv->axpy( -1.0, prob->xlowerBound() );
    rv->selectNonZeros( *ixlow );
    rv->axpy( -1.0, *vars->v );

    norm = updateNormAndPrint( norm, *rv, print_resids, "rv");

    // contribution - lx^T gamma to duality gap
    gap -= prob->blx->dotProductWith(*vars->gamma);
  }

  if( nxupp > 0 )
  {
    /*** rw = x - ux + w ***/
    rw->copyFrom( *vars->x );
    rw->axpy( -1.0, prob->xupperBound() );
    rw->selectNonZeros( *ixupp );
    rw->axpy(  1.0, *vars->w );

    norm = updateNormAndPrint( norm, *rw, print_resids, "rw");

    // contribution + bu^T phi to duality gap
    gap += prob->bux->dotProductWith(*vars->phi);
  }
   
  mDualityGap = gap;
  mResidualNorm = norm;

  if( print_resids && myRank == 0 )
    std::cout << "Norm residuals: " << mResidualNorm << "\tduality gap: " << mDualityGap << "\n";

}

double QpGenResiduals::recomputeResidualNorm()
{
   mResidualNorm = 0.0;

   double componentNorm = 0.0;
   componentNorm = rQ->infnorm();

   if( componentNorm > mResidualNorm )
      mResidualNorm = componentNorm;

   componentNorm = rA->infnorm();
   if( componentNorm > mResidualNorm )
      mResidualNorm = componentNorm;

   componentNorm = rC->infnorm();
   if( componentNorm > mResidualNorm )
      mResidualNorm= componentNorm;

   if( mclow > 0 )
   {
     componentNorm = rt->infnorm();
     if( componentNorm > mResidualNorm )
        mResidualNorm = componentNorm;
   }

   if( mcupp > 0 )
   {
     componentNorm = ru->infnorm();
     if( componentNorm > mResidualNorm )
        mResidualNorm = componentNorm;
   }

   componentNorm = rz->infnorm();
   if( componentNorm > mResidualNorm )
      mResidualNorm = componentNorm;

   if( nxlow > 0 )
   {
     componentNorm = rv->infnorm();
     if( componentNorm > mResidualNorm )
        mResidualNorm = componentNorm;
   }

   if( nxupp > 0 )
   {
     componentNorm = rw->infnorm();
     if( componentNorm > mResidualNorm )
        mResidualNorm = componentNorm;
   }
   return mResidualNorm;
}

void QpGenResiduals::add_r3_xz_alpha(const Variables *vars_in, double alpha)
{
  QpGenVars * vars = (QpGenVars *) vars_in;

  if( mclow > 0 ) rlambda->axzpy( 1.0, *vars->t, *vars->lambda );
  if( mcupp > 0 ) rpi    ->axzpy( 1.0, *vars->u, *vars->pi );
  if( nxlow > 0 ) rgamma ->axzpy( 1.0, *vars->v, *vars->gamma );
  if( nxupp > 0 ) rphi   ->axzpy( 1.0, *vars->w, *vars->phi );

  if( alpha != 0.0 )
  {
    if( mclow > 0 ) rlambda->addSomeConstants( alpha, *iclow );
    if( mcupp > 0 ) rpi    ->addSomeConstants( alpha, *icupp );
    if( nxlow > 0 ) rgamma ->addSomeConstants( alpha, *ixlow );
    if( nxupp > 0 ) rphi   ->addSomeConstants( alpha, *ixupp );
  }
}

void QpGenResiduals::set_r3_xz_alpha(const Variables *vars, double alpha)
{
  this->clear_r3();
  this->add_r3_xz_alpha( vars, alpha );
}
  
void QpGenResiduals::clear_r3()
{
  if( mclow > 0 ) rlambda->setToZero();
  if( mcupp > 0 ) rpi    ->setToZero();
  if( nxlow > 0 ) rgamma ->setToZero();
  if( nxupp > 0 ) rphi   ->setToZero();
}
  
void QpGenResiduals::clear_r1r2()
{
  rQ->setToZero();
  rA->setToZero();
  rC->setToZero();
  rz->setToZero();
  if( nxlow > 0 ) rv->setToZero();
  if( nxupp > 0 ) rw->setToZero();
  if( mclow > 0 ) rt->setToZero();
  if( mcupp > 0 ) ru->setToZero();
}

void QpGenResiduals::project_r3(double rmin, double rmax)
{
  if( mclow > 0 ) {
    rlambda->gondzioProjection( rmin, rmax );
    rlambda->selectNonZeros( *iclow );
  }
  if( mcupp > 0 ) {
    rpi    ->gondzioProjection( rmin, rmax );
    rpi    ->selectNonZeros( *icupp );
  }
  if( nxlow > 0 ) {
    rgamma ->gondzioProjection( rmin, rmax );
    rgamma ->selectNonZeros( *ixlow );
  }
  if( nxupp > 0 ) {
    rphi   ->gondzioProjection( rmin, rmax );
    rphi   ->selectNonZeros( *ixupp );
  }

}
  

int QpGenResiduals::validNonZeroPattern()
{
  if( nxlow > 0 && 
      ( !rv    ->matchesNonZeroPattern( *ixlow ) ||
	!rgamma->matchesNonZeroPattern( *ixlow ) ) ) {
    return 0;
  }

  if( nxupp > 0 &&
      ( !rw  ->matchesNonZeroPattern( *ixupp ) ||
	!rphi->matchesNonZeroPattern( *ixupp ) ) ) {
    return 0;
  }
  if( mclow > 0 &&
      ( !rt     ->matchesNonZeroPattern( *iclow ) ||
	!rlambda->matchesNonZeroPattern( *iclow ) ) ) {
    return 0;
  }

  if( mcupp > 0 &&
      ( !ru ->matchesNonZeroPattern( *icupp ) ||
	!rpi->matchesNonZeroPattern( *icupp ) ) ) {
    return 0;
  }
  
  return 1;
}

void QpGenResiduals::writeToStream(std::ostream& out)
{
  /*
  printf("--------------rQ\n");
  rQ->writeToStream(out);printf("---------------------------\n");
  
  printf("rA\n");
  rA->writeToStream(out);printf("---------------------------\n");
  */
  printf("rC\n");
  rC->writeToStream(out);printf("---------------------------\n");


  printf("rz\n");
  rz->writeToStream(out);printf("---------------------------\n");
  /*  
  if ( mclow > 0 ) {
    printf("rt\n");
    rt->writeToStream(out);printf("---------------------------\n");
    printf("rlambda\n");
    rlambda->writeToStream(out);printf("---------------------------\n");
  }
  if ( mcupp > 0 ) {
    printf("ru\n");
    ru->writeToStream(out);printf("---------------------------\n");
    printf("rpi\n");
    rpi->writeToStream(out);printf("---------------------------\n");
  }

  
  if( nxlow > 0 ) {
    printf("rv\n");
    rv->writeToStream(out);printf("---------------------------\n");
    printf("rgamma\n");
    rgamma->writeToStream(out);printf("---------------------------\n");
  }
  if( nxupp > 0 ) {
    printf("rw\n");
    rw->writeToStream(out);printf("---------------------------\n");
    printf("rphi\n");
    rphi->writeToStream(out);printf("---------------------------\n");
    }
  */
}
