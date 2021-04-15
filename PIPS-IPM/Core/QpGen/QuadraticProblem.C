/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QuadraticProblem.h"
#include "QpGenVars.h"
#include "DoubleMatrix.h"
#include "OoqpVector.h"
#include <cmath>

#include "SimpleVector.h"
#include "LinearAlgebraPackage.h"
#include "MpsReader.h"

QuadraticProblem::QuadraticProblem( LinearAlgebraPackage * la_in, OoqpVector * c_in, SymMatrix * Q_in, OoqpVector * xlow_in, OoqpVector * ixlow_in,
	 OoqpVector * xupp_in, OoqpVector * ixupp_in, GenMatrix  * A_in, OoqpVector * bA_in, GenMatrix  * C_in, OoqpVector * clow_in, OoqpVector * iclow_in, OoqpVector * cupp_in, OoqpVector * icupp_in ) :
        la{la_in},
        nxlow{ ixlow_in->numberOfNonzeros() },
        nxupp{ ixupp_in->numberOfNonzeros() },
        mclow{ iclow_in->numberOfNonzeros() },
        mcupp{ icupp_in->numberOfNonzeros() }
{
  SpReferTo( g,     c_in  );
  SpReferTo( bA,    bA_in );
  SpReferTo( blx,   xlow_in  );
  SpReferTo( ixlow, ixlow_in );
  SpReferTo( bux,   xupp_in  );
  SpReferTo( ixupp, ixupp_in );
  SpReferTo( bl,    clow_in  );
  SpReferTo( iclow, iclow_in );
  SpReferTo( bu,    cupp_in  );
  SpReferTo( icupp, icupp_in );

  long long dummy;

  nx = g->length();
  SpReferTo( Q, Q_in );

  SpReferTo( A, A_in );
  A->getSize( my, dummy );
  
  SpReferTo( C, C_in );
  C->getSize( mz, dummy );
}

void QuadraticProblem::Qmult( double beta,  OoqpVector& y,
		       double alpha, const OoqpVector& x ) const
{
  Q->mult( beta, y, alpha, x );
}

void QuadraticProblem::Amult( double beta,  OoqpVector& y,
		       double alpha, const OoqpVector& x) const
{
  A->mult( beta, y, alpha, x );
}

void QuadraticProblem::Cmult( double beta,  OoqpVector& y,
		       double alpha, const OoqpVector& x ) const
{
  C->mult( beta, y, alpha, x );
}

void QuadraticProblem::ATransmult( double beta,  OoqpVector& y,
			    double alpha, const OoqpVector& x ) const
{
  A->transMult( beta, y, alpha, x );
}

void QuadraticProblem::CTransmult( double beta,  OoqpVector& y,
			    double alpha, const OoqpVector& x ) const
{
  C->transMult( beta, y, alpha, x );
}

void QuadraticProblem::getg( OoqpVector& myG ) const
{
  myG.copyFrom( *g );
}

void QuadraticProblem::getbA( OoqpVector& bout ) const
{
  bout.copyFrom( *bA );
}

double QuadraticProblem::datanorm() const
{
  double norm = 0.0;
  double componentNorm;

  componentNorm = g->infnorm();
  if( componentNorm > norm ) norm = componentNorm;
  
  componentNorm = Q->abmaxnorm();
  if( componentNorm > norm ) norm = componentNorm;

  componentNorm = bA->infnorm();
  if( componentNorm > norm ) norm = componentNorm;

  componentNorm = A->abmaxnorm();
  if( componentNorm > norm ) norm = componentNorm;

  componentNorm = C->abmaxnorm();
  if( componentNorm > norm ) norm = componentNorm;

  assert( blx->matchesNonZeroPattern( *ixlow ) );
  componentNorm = blx->infnorm();
  if( componentNorm > norm ) norm = componentNorm;

  assert( bux->matchesNonZeroPattern( *ixupp ) );
  componentNorm = bux->infnorm();
  if( componentNorm > norm ) norm = componentNorm;

  assert( bl->matchesNonZeroPattern( *iclow ) );
  componentNorm = bl->infnorm();
  if( componentNorm > norm ) norm = componentNorm;

  assert( bu->matchesNonZeroPattern( *icupp ) );
  componentNorm = bu->infnorm();
  if( componentNorm > norm ) norm = componentNorm;

  return norm;
}

void QuadraticProblem::datainput( MpsReader * reader, int& iErr ) 
{
    reader->readQpGen( *g, *Q, *blx, *ixlow, *bux, *ixupp,
		     *A, *bA,
		     *C, *bl, *iclow, *bu, *icupp, iErr );

    if( reader->scalingOption == 1){
        // Create the scaling vector
        this->createScaleFromQ();

        //Scale the variables
        this->scaleQ();
        this->scaleA();
        this->scaleC();
        this->scaleg();
        this->scalexlow();
        this->scalexupp();
        }

    /* If objective sense is "MAX", flip the C and Q matrices */
    if( !strncmp( reader->objectiveSense, "MAX", 3)){
        this->flipg();
        this->flipQ();
        }  
}

void 
QuadraticProblem::randomlyChooseBoundedVariables( OoqpVector& x,
					   OoqpVector& dualx,
					   OoqpVector& xlow_,
					   OoqpVector& ixlow_,
					   OoqpVector& xupp_,
					   OoqpVector& ixupp_,
					   double * ix,
					   double percentLowerOnly,
					   double percentUpperOnly,
					   double percentBound )
{
  int i;
  double drand( double * );
  // Initialize the upper and lower bounds on x
  int n = x.length();

  
  double * sxlow  = new double[n];
  double * sixlow = new double[n];
  double * sxupp  = new double[n];
  double * sixupp = new double[n];
  double * sx     = new double[n];
  double * sdualx = new double[n];
  
  for( i = 0; i < n; i++ ) {
    double r = drand(ix);
    //cout << " r: " << r << "   ";
	
    if( r < percentLowerOnly ) {
      //cout << "i= " << i << " Lower bound " << endl;
      sixlow[i]  = 1.0;
      sxlow[i]    = (drand(ix) - 0.5) * 3.0;
      sixupp[i]  = 0.0;
      sxupp[i]    = 0.0;
    } else if ( r < percentLowerOnly + percentUpperOnly ) { 
      //cout << "i= " << i << " Upper bound " << endl;
      sixlow[i]  = 0.0;
      sxlow[i]    = 0.0;
      sixupp[i]  = 1.0;
      sxupp[i]    = (drand(ix) - 0.5) * 3.0;
    } else if ( r < percentLowerOnly + percentUpperOnly
		+ percentBound ) {
      //cout << "i= " << i << " Two-sided bound " << endl;
      sixlow[i]  = 1.0;
      sxlow[i]    = (drand(ix) - 0.5) * 3.0;
      sixupp[i]  = 1.0;
      sxupp[i]    = sxlow[i] +  drand(ix) * 10.0;
    } else {
      // it is free
      //cout << "i= " << i << " Free " << endl;
      sixlow[i]  = 0.0;
      sxlow[i]    = 0.0;
      sixupp[i]  = 0.0;
      sxupp[i]    = 0.0;
    }
  }
  xlow_. copyFromArray( sxlow );
  ixlow_.copyFromArray( sixlow );
  xupp_. copyFromArray( sxupp );
  ixupp_.copyFromArray( sixupp );

  for ( i = 0; i < n; i++ ) {
    if( sixlow[i] == 0.0 && sixupp[i] == 0.0 ) {
      // x[i] not bounded
      sx[i] = 20.0 * drand(ix) - 10.0;
      sdualx[i] = 0.0;
    } else if ( sixlow[i] != 0.0 && sixupp[i] != 0.0 ) {
      // x[i] is bounded above and below
      double r = drand( ix );
      if( r < 0.33 ) {
	// x[i] is on its lower bound
	sx[i]     =  sxlow[i];
	sdualx[i] =  10.0 * drand( ix );
      } else if ( r > .66 ) {
	// x[i] is on its upper bound
	sx[i]     =  sxupp[i];
	sdualx[i] = -10.0 * drand( ix );
      } else {
	// x[i] is somewhere in between
	double theta = .99 * drand( ix ) + .005;
	sx[i] = (1 - theta) * sxlow[i] + theta * sxupp[i];
	sdualx[i] = 0.0;
      }
    } else if ( sixlow[i] != 0.0 ) {
      // x[i] is only bounded below
      if( drand( ix ) < .33 ) {
	// x[i] is on its lower bound
	sx[i]     =  sxlow[i];
	sdualx[i] =  10.0 * drand( ix );
      } else {
	// x[i] is somewhere above its lower bound
	sx[i]     = sxlow[i] + 0.005 + 10.0 * drand(ix);
	sdualx[i] = 0.0;
      }
    } else { // x[i] only has an upper bound
      if( drand(ix) > .66 ) {
	// x[i] is on its upper bound
	sx[i]     =  sxupp[i];
	sdualx[i] = -10.0 * drand( ix );
      } else {
	// x[i] is somewhere below its upper bound
	sx[i]     =  sxupp[i] - 0.005 - 10.0 * drand(ix);
	sdualx[i] = 0.0;
      }
    } // end else x[i] only has an upper bound
  } // end for ( i = 0; i < n; i++ )
  x.copyFromArray( sx );
  dualx.copyFromArray( sdualx );

  delete [] sxlow;
  delete [] sxupp;
  delete [] sixlow;
  delete [] sixupp;
  delete [] sx;
  delete [] sdualx;
}

void QuadraticProblem::print()
{
  std::cout << "begin Q\n";
  Q->writeToStream( std::cout );
  std::cout << "end Q\n";
  std::cout << "begin c\n";
  g->writeToStream( std::cout );
  std::cout << "end c\n";

  std::cout << "begin xlow\n";
  blx->writeToStream( std::cout );
  std::cout << "end xlow\n";
  std::cout << "begin ixlow\n";
  ixlow->writeToStream( std::cout );
  std::cout << "end ixlow\n";

  std::cout << "begin xupp\n";
  bux->writeToStream( std::cout );
  std::cout << "end xupp\n";
  std::cout << "begin ixupp\n";
  ixupp->writeToStream( std::cout );
  std::cout << "end ixupp\n";
  std::cout << "begin A\n";

  A->writeToStream( std::cout );
  std::cout << "end A\n";
  std::cout << "begin b\n";
  bA->writeToStream( std::cout );
  std::cout << "end b\n";
  std::cout << "begin C\n";
  C->writeToStream( std::cout );
  std::cout << "end C\n";
  
  std::cout << "begin clow\n";
  bl->writeToStream( std::cout );
  std::cout << "end clow\n";
  std::cout << "begin iclow\n";
  iclow->writeToStream( std::cout );
  std::cout << "end iclow\n";

  std::cout << "begin cupp\n";
  bu->writeToStream( std::cout );
  std::cout << "end cupp\n";
  std::cout << "begin icupp\n";
  icupp->writeToStream( std::cout );
  std::cout << "end icupp\n";

}

#include <fstream>

void QuadraticProblem::datarandom( OoqpVector & x, OoqpVector & y,
			    OoqpVector & z, OoqpVector & s )
{
  double drand( double * );
  double ix = 3074.20374;
  
  OoqpVectorHandle xdual(la->newVector( nx ));
  this->randomlyChooseBoundedVariables( x, *xdual,
  					*blx, *ixlow, *bux, *ixupp,
  					&ix, .25, .25, .25 );

  {
//      ofstream x_vec( "x" );
//      x->writeToStream( x_vec );
//      ofstream eta_vec( "xdual" );
//      eta_vec.precision(16);
//      xdual->writeToStream( eta_vec );
//      ofstream blx_vec( "blx" );
//      blx->writeToStream( blx_vec );
//      ofstream bux_vec( "bux" );
//      bux->writeToStream( bux_vec );
  } 

  OoqpVectorHandle sprime(la->newVector( mz ));
  this->randomlyChooseBoundedVariables( *sprime, z,
  					*bl, *iclow, *bu, *icupp,
  					&ix, .25, .25, .5 );
  
  {
    //      	ofstream z_vec( "z" );
    //      	z_vec << z << endl;
  }
  {
    Q->randomizePSD( &ix );
//      ofstream Q_mat( "Q" );
//      Q_mat << Q << endl;
  }
  { 
    A->randomize( -10.0, 10.0, &ix );
    //  	ofstream A_mat( "A" );
    //  	A_mat << A << endl;
  }
  {
    C->randomize( -10.0, 10.0, &ix );
    //  	ofstream C_mat( "C" );
    //  	C_mat << C << endl;
  }

  y.randomize( -10.0, 10.0, &ix );
  {
    //  	ofstream y_vec( "y" );
    //  	y_vec << y << endl;
  }
  // g = - Q x + A\T y + C\T z + xdual 
  g->copyFrom( *xdual );
  Q->mult( 1.0, *g, -1.0, x );
  A->transMult( 1.0, *g, 1.0, y );
  C->transMult( 1.0, *g, 1.0, z );
  // bA = A x
  A->mult( 0.0, *bA, 1.0, x );
  {
    //  	ofstream bA_vec( "bA" );
    //  	bA_vec << bA << endl;
  }
  // Have a randomly generated sprime.
  // C x - s = 0. Let q + sprime = s, i.e. q = s - sprime.
  // Compute s and temporarily store in q
  C->mult( 0.0, s, 1.0, x );
  // Now compute the real q = s - sprime
  OoqpVectorHandle q(la->newVector( mz ));
  q->copyFrom( s );
  q->axpy( -1.0, *sprime );
  // Adjust bl and bu appropriately
  bl->axpy( 1.0, *q );
  bu->axpy( 1.0, *q );
  
  bl->selectNonZeros( *iclow );
  bu->selectNonZeros( *icupp );
  
  {
    //   	ofstream bl_vec( "bl" );
    //      	bl_vec << bl << endl;
    //      	ofstream bu_vec( "bu" );
    //      	bu_vec << bu << endl;
  }
}


void QuadraticProblem::putQIntoAt( GenMatrix& M, int row, int col )
{
  M.atPutSubmatrix( row, col, *Q, 0, 0, nx, nx );
}

void QuadraticProblem::putQIntoAt( SymMatrix& M, int row, int col )
{
  M.symAtPutSubmatrix( row, col, *Q, 0, 0, nx, nx );
}

void QuadraticProblem::putAIntoAt( GenMatrix& M, int row, int col )
{
  M.atPutSubmatrix( row, col, *A, 0, 0, my, nx );
}

void QuadraticProblem::putAIntoAt( SymMatrix& M, int row, int col )
{
  M.symAtPutSubmatrix( row, col, *A, 0, 0, my, nx );
}

void QuadraticProblem::putCIntoAt( GenMatrix& M, int row, int col )
{
  M.atPutSubmatrix( row, col, *C, 0, 0, mz, nx );
}

void QuadraticProblem::putCIntoAt( SymMatrix& M, int row, int col )
{
  M.symAtPutSubmatrix( row, col, *C, 0, 0, mz, nx );
}

void QuadraticProblem::getDiagonalOfQ( OoqpVector& dq )
{
  Q->fromGetDiagonal(0, dq);
}

void QuadraticProblem::objective_gradient(const QpGenVars* vars, OoqpVector& gradient) {
   this->getg(gradient);
   this->Qmult(1., gradient, 1., *vars->x);
   return;
}

double QuadraticProblem::objective_value(const QpGenVars * vars ) const {
  OoqpVectorHandle gradient(la->newVector(nx));
  this->getg( *gradient );
  this->Qmult(1., *gradient, 0.5, *vars->x );

  return gradient->dotProductWith(*vars->x );
}

void QuadraticProblem::createScaleFromQ()
{
  // Stuff the diagonal elements of Q into the vector "sc"
  this->getDiagonalOfQ( *sc);

  // Modifying scVector is equivalent to modifying sc
  SimpleVector & scVector = dynamic_cast<SimpleVector &>(*sc);

  int scLength = scVector.length();

  for( int i = 0; i < scLength; i++){
    if( scVector[i] > 1)
        scVector[i] = 1.0/sqrt( scVector[i]);
    else
        scVector[i] = 1.0;
    }
}

void QuadraticProblem::scaleQ()
{
  Q->symmetricScale( *sc);
}


void QuadraticProblem::scaleA()
{
  A->columnScale( *sc);
}

void QuadraticProblem::scaleC()
{
  C->columnScale( *sc);
}

void QuadraticProblem::scaleg()
{
  SimpleVector & scVector = dynamic_cast<SimpleVector &>(*sc);

  assert ( scVector.length() == g->length());

  // D * g
  g->componentMult( scVector);
}

void QuadraticProblem::scalexupp()
{
  SimpleVector & scVector = dynamic_cast<SimpleVector &>(*sc);

  assert ( scVector.length() == bux->length());

  // inverse(D) * bux
  bux->componentDiv( scVector);

}


void QuadraticProblem::scalexlow()
{
  SimpleVector & scVector = dynamic_cast<SimpleVector &>(*sc);

  assert ( scVector.length() == blx->length());

  // inverse(D) * blx
  blx->componentDiv( scVector);

}

void QuadraticProblem::flipg()
{
  // Multiply C matrix by -1
  g->scalarMult( -1.0);
}

void QuadraticProblem::flipQ()
{
  // Multiply Q matrix by -1
  Q->scalarMult( -1.0);
}

