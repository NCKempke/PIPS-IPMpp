//
// Created by nils-christian on 18.06.21.
//

#ifndef PIPSIPMPP_CURTISREIDSCALER_H
#define PIPSIPMPP_CURTISREIDSCALER_H

#include "Scaler.hpp"
#include "DistributedVector.h"

class Problem;

class CurtisReidScaler : Scaler {
private:
   std::tuple<std::unique_ptr<Vector<int>>, std::unique_ptr<Vector<int>>, std::unique_ptr<Vector<int>>> get_nonzero_vectors() const;
   std::tuple<std::unique_ptr<Vector<double>>, std::unique_ptr<Vector<double>>, std::unique_ptr<Vector<double>>> get_log_sum_vectors() const;
protected:
   void doObjScaling() const override;

   void applyGeoMean(Vector<double>& maxvec, const Vector<double>& minvec);
public:

   CurtisReidScaler(const ProblemFactory& problem_factory, const Problem& problem, bool bitshifting = false);
   ~CurtisReidScaler() override = default;

   /** scale */
   void scale() override;
};

#endif //PIPSIPMPP_CURTISREIDSCALER_H
