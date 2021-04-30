#ifndef STOCHVECTORHANDLE
#define STOCHVECTORHANDLE

#include "IotrRefCount.h"
#include "SmartPointer.h"
#include "pipsport.h"

template<typename T>
class DistributedVector;

template<typename T>
class StochDummyVectorBase;

template<typename T> using StochVectorBaseHandle = SmartPointer<DistributedVector<T> >;
typedef SmartPointer<DistributedVector<double>> StochVectorHandle;

#endif
