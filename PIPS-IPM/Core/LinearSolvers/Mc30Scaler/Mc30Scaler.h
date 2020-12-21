/*
 * Mc30Scaler.h
 *
 *  Created on: Nov 6, 2020
 *      Author: bzfkempk
 */

#ifndef PIPSIPM_CORE_LINEARSOLVERS_MC30SCALER_H
#define PIPSIPM_CORE_LINEARSOLVERS_MC30SCALER_H

#include "DoubleLinearSolver.h"
#include "pipsport.h"

#include <vector>

#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif

extern "C"
{
   void FNAME(mc30ad)( int* n, int* ne, const double a[], const int irn[], const int icn[],
         double s[], double w[], int* lp, int* ifail);
}

class Mc30Scaler : public SymmetricLinearScaler
{
   public:
      Mc30Scaler() = default;

      ~Mc30Scaler() override = default;

      const double* getScaling() const override { return scaling_factors.data(); }

      /** compute scaling for a symmetric indefinite matrix and scale it */
      void scaleMatrixTripletFormat( int n, int nnz, double* M, const int* rowM, const int* colM, bool fortran_indexed ) override;

      void scaleMatrixCSRFormat( int, int, double*, const int*, const int*, bool) override
      { assert( false && "TODO : implement" ); };

      /* scale a vector */
      void scaleVector( OoqpVector& vec_in ) const override;

      /* unscale a vector */
      void unscaleVector( OoqpVector& vec_in ) const override;

   private:
      void getFortranIndex(const int* rowM, const int* colM, int length, int max_index);

      std::vector<int> rowM_ft_indexed;
      std::vector<int> colM_ft_indexed;

      std::vector<double> scaling_factors;
      std::vector<double> scaling_workspace;

      int scaling_output_control{1};
      int scaling_error{0};
};

#endif /* PIPSIPM_CORE_LINEARSOLVERS_MC30SCALER_H */
