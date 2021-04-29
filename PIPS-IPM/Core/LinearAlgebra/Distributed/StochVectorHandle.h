#ifndef STOCHVECTORHANDLE
#define STOCHVECTORHANDLE

#include "IotrRefCount.h"
#include "SmartPointer.h"
#include "StochVector_fwd.h"
#include "pipsport.h"

template<typename T> using StochVectorBaseHandle = SmartPointer<StochVectorBase<T> >;
typedef SmartPointer<StochVector> StochVectorHandle;

#endif
