/*
 * pipsdef.h
 *
 *  Created on: 29.01.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_UTILITIES_PIPSDEF_H_
#define PIPS_IPM_CORE_UTILITIES_PIPSDEF_H_

// save diagnostic state
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-override"

#include "mpi.h"
// turn the warnings back on
#pragma GCC diagnostic pop

#include <iostream>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <limits>

#include "omp.h"
#include <algorithm>

//#define PIPS_DEBUG

using Permutation = std::vector<unsigned int>;

const double pips_eps = 1e-13;
const double pips_eps0 = 1e-40;

static const double feastol = 1.0e-6; // was 1.0e-6
static const double infinity = 1.0e30;

static inline double relativeDiff(double val1, double val2) {
   const double val1Abs = std::fabs(val1);
   const double val2Abs = std::fabs(val2);
   const double div = std::max(1.0, std::max(val1Abs, val2Abs));

   return (val1 - val2) / div;
}

#ifdef PIPS_DEBUG
#define PIPSdebugMessage                printf("[%s:%d] debug: ", __FILE__, __LINE__), printf
#else
#define PIPSdebugMessage      while( 0 ) /*lint -e{530}*/ printf
#endif

inline bool PIPSisEQFeas(double val1, double val2) {
   return (std::fabs(val1 - val2) <= feastol);
}

inline bool PIPSisEQ(double val1, double val2, double eps = pips_eps) {
   return (std::fabs(val1 - val2) <= eps);
}

inline bool PIPSisRelEQ(double val1, double val2, double eps = pips_eps) {
   const double reldiff = relativeDiff(val1, val2);

   return (std::fabs(reldiff) <= eps);
}

inline bool PIPSisRelEQFeas(double val1, double val2) {
   const double reldiff = relativeDiff(val1, val2);

   return std::fabs(reldiff) <= feastol;
}

inline bool PIPSisLE(double val1, double val2, double eps = pips_eps) {
   return (val1 <= val2 + eps);
}

inline bool PIPSisLEFeas(double val1, double val2) {
   return (val1 <= val2 + feastol);
}

inline bool PIPSisRelLEFeas(double val1, double val2) {
   const double reldiff = relativeDiff(val1, val2);

   return reldiff <= feastol;
}

inline bool PIPSisLT(double val1, double val2, double eps = pips_eps) {
   return (val1 < val2 - eps);
}

inline bool PIPSisLTFeas(double val1, double val2) {
   return (val1 < val2 - feastol);
}

inline bool PIPSisRelLT(double val1, double val2, double eps = pips_eps) {
   const double reldiff = relativeDiff(val1, val2);

   return (reldiff < -eps);
}

inline bool PIPSisRelLTFeas(double val1, double val2) {
   const double reldiff = relativeDiff(val1, val2);

   return (reldiff < -feastol);
}

inline bool PIPSisZero(double val, double eps0 = pips_eps0) {
   return (std::fabs(val) < eps0);
}

inline bool PIPSisZeroFeas(double val) {
   return (std::fabs(val) < feastol);
}

template<typename T>
inline bool isInVector(const T& elem, const std::vector<T>& vec) {
   return std::find(vec.begin(), vec.end(), elem) != vec.end();
}

template<typename T>
inline bool containsSorted(const std::vector<T>& subset, const std::vector<T>& vec) {
   assert(std::is_sorted(subset.begin(), subset.end()));
   assert(std::is_sorted(vec.begin(), vec.end()));

   return std::includes(vec.begin(), vec.end(), subset.begin(), subset.end());
}

template<typename T>
inline bool contains(const std::vector<T>& subset, const std::vector<T>& vec) {
   std::vector<T> subset_cpy(subset);
   std::vector<T> vec_cpy(vec);

   std::sort(subset_cpy.begin(), subset_cpy.end());
   std::sort(vec_cpy.begin(), vec_cpy.end());

   return containsSorted(subset_cpy, vec_cpy);
}

void doubleLexSort(int first[], int n, int second[], double data[]);


template<typename T>
inline void permuteVector(const Permutation& perm, std::vector<T>& vec) {
   assert(perm.size() == vec.size());

   std::vector<T> tmp(vec.size());

   for (size_t i = 0; i < vec.size(); ++i)
      tmp[i] = vec[perm[i]];
   vec = tmp;
}

inline Permutation getInversePermutation(const Permutation& perm) {
   size_t size = perm.size();
   Permutation perm_inv(size, 0);

   for (size_t i = 0; i < size; i++)
      perm_inv[perm[i]] = i;

   return perm_inv;
}

template<typename T>
inline unsigned int getNDistinctValues(const std::vector<T>& values) {
   return std::set<T>(values.begin(), values.end()).size();
}

inline int PIPSgetnOMPthreads() {
   return omp_get_max_threads();
}

inline MPI_Comm PIPS_MPIcreateGroupFromRanks(const int* chosen_ranks, unsigned int n_chosen_ranks, MPI_Comm mpi_comm_all = MPI_COMM_WORLD) {
   MPI_Group all_group;
   MPI_Comm_group(MPI_COMM_WORLD, &all_group);

   MPI_Group sub_group;
   MPI_Group_incl(all_group, n_chosen_ranks, chosen_ranks, &sub_group);

   MPI_Comm sub_comm;
   MPI_Comm_create(mpi_comm_all, sub_group, &sub_comm);

   MPI_Group_free(&sub_group);
   return sub_comm;
}

inline MPI_Comm PIPS_MPIcreateGroupFromRanks(const std::vector<int> chosen_ranks, MPI_Comm mpi_comm_all = MPI_COMM_WORLD) {
   return PIPS_MPIcreateGroupFromRanks(chosen_ranks.data(), chosen_ranks.size(), mpi_comm_all);
}

inline bool PIPS_MPIiAmSpecial(int iAmDistrib, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   bool iAmSpecial = true;

   if (iAmDistrib) {
      int rank;

      MPI_Comm_rank(mpiComm, &rank);
      if (rank > 0)
         iAmSpecial = false;
   }

   return iAmSpecial;
}

inline bool iAmSpecial(int iAmDistrib, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   return PIPS_MPIiAmSpecial(iAmDistrib, mpiComm);
}

inline int PIPS_MPIgetSize(MPI_Comm comm = MPI_COMM_WORLD) {
   if (comm == MPI_COMM_NULL)
      return -1;
   int mysize;
   MPI_Comm_size(comm, &mysize);
   return mysize;
}

inline bool PIPS_MPIgetDistributed(MPI_Comm comm = MPI_COMM_WORLD) {
   return PIPS_MPIgetSize(comm) > 1;
}

void inline PIPS_MPIabortInfeasible(std::string message, std::string file, std::string function, MPI_Comm comm = MPI_COMM_WORLD) {
   std::cerr << "Infesibility detected in " << file << " function " << function << "!\n";
   std::cerr << "Message: " << message << "\n";
   std::cerr << "Aborting now.\n";
   MPI_Abort(comm, 1);
}

void inline _PIPS_MPIabortIf(bool cond, const std::string& message, const char* file, int line) {
   if (cond) {
      std::cerr << "Called MPI_Abort in " << file << " at " << line << " saying " << message << "\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
}
#define PIPS_MPIabortIf(cond, message) _PIPS_MPIabortIf( cond, message, __FILE__, __LINE__ )

inline std::vector<int> PIPSallgathervInt(const std::vector<int>& vecLocal, MPI_Comm mpiComm) {
   int myrank;
   MPI_Comm_rank(mpiComm, &myrank);
   int mysize;
   MPI_Comm_size(mpiComm, &mysize);

   std::vector<int> vecGathered;

   if (mysize > 0) {
      // get all lengths
      std::vector<int> recvcounts(mysize);
      std::vector<int> recvoffsets(mysize);

      int lengthLocal = int(vecLocal.size());

      MPI_Allgather(&lengthLocal, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, mpiComm);

      // all-gather local components
      recvoffsets[0] = 0;
      for (int i = 1; i < mysize; ++i)
         recvoffsets[i] = recvoffsets[i - 1] + recvcounts[i - 1];

      const int lengthGathered = recvoffsets[mysize - 1] + recvcounts[mysize - 1];

      vecGathered = std::vector<int>(lengthGathered);

      MPI_Allgatherv(&vecLocal[0], lengthLocal, MPI_INT, &vecGathered[0], &recvcounts[0], &recvoffsets[0], MPI_INT, mpiComm);
   }
   else {
      vecGathered = vecLocal;
   }

   return vecGathered;
}

inline std::vector<int> PIPSallgathervInt(const std::vector<int>& vecLocal, MPI_Comm mpiComm, int& startMy, int& endMy) {
   int myrank;
   MPI_Comm_rank(mpiComm, &myrank);
   int mysize;
   MPI_Comm_size(mpiComm, &mysize);

   std::vector<int> vecGathered;

   if (mysize > 0) {
      // get all lengths
      std::vector<int> recvcounts(mysize);
      std::vector<int> recvoffsets(mysize);

      int lengthLocal = int(vecLocal.size());

      MPI_Allgather(&lengthLocal, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, mpiComm);

      // all-gather local components
      recvoffsets[0] = 0;
      for (int i = 1; i < mysize; ++i)
         recvoffsets[i] = recvoffsets[i - 1] + recvcounts[i - 1];

      const int lengthGathered = recvoffsets[mysize - 1] + recvcounts[mysize - 1];

      startMy = recvoffsets[myrank];

      // not the last rank?
      if (myrank < mysize - 1) {
         endMy = recvoffsets[myrank + 1];
      }
      else {
         assert(myrank == mysize - 1);
         endMy = lengthGathered;
      }

      vecGathered = std::vector<int>(lengthGathered);

      MPI_Allgatherv(&vecLocal[0], lengthLocal, MPI_INT, &vecGathered[0], &recvcounts[0], &recvoffsets[0], MPI_INT, mpiComm);
   }
   else {
      startMy = 0;
      endMy = int(vecLocal.size());
      vecGathered = vecLocal;
   }

   return vecGathered;
}

template<typename T>
struct get_mpi_datatype_t {
   static inline MPI_Datatype get();
};

// specialization for particular types:
template<>
inline MPI_Datatype get_mpi_datatype_t<int>::get() {
   return MPI_INT;
};

template<>
inline MPI_Datatype get_mpi_datatype_t<long>::get() {
   return MPI_LONG;
};

template<>
inline MPI_Datatype get_mpi_datatype_t<long long>::get() {
   return MPI_LONG_LONG;
};

template<>
inline MPI_Datatype get_mpi_datatype_t<double>::get() {
   return MPI_DOUBLE;
};

template<>
inline MPI_Datatype get_mpi_datatype_t<char>::get() {
   return MPI_CHAR;
};

template<>
inline MPI_Datatype get_mpi_datatype_t<bool>::get() {
   return MPI_CXX_BOOL;
};

template<>
inline MPI_Datatype get_mpi_datatype_t<unsigned int>::get() {
   return MPI_UNSIGNED;
};

template<>
inline MPI_Datatype get_mpi_datatype_t<long unsigned int>::get() {
   return MPI_UNSIGNED_LONG;
};

template<typename T>
inline MPI_Datatype get_mpi_datatype(const T&) {
   return get_mpi_datatype_t<T>::get();
};

template<typename T>
inline MPI_Datatype get_mpi_datatype(T*) {
   return get_mpi_datatype_t<T>::get();
};

template<typename T>
inline MPI_Datatype get_mpi_datatype(const T* const) {
   return get_mpi_datatype_t<T>::get();
};

template<typename T>
struct get_mpi_locdatatype_t {
   static inline MPI_Datatype get();
};

template<>
inline MPI_Datatype get_mpi_locdatatype_t<float>::get() {
   return MPI_FLOAT_INT;
};

template<>
inline MPI_Datatype get_mpi_locdatatype_t<double>::get() {
   return MPI_DOUBLE_INT;
};

template<>
inline MPI_Datatype get_mpi_locdatatype_t<long>::get() {
   return MPI_LONG_INT;
};

template<>
inline MPI_Datatype get_mpi_locdatatype_t<int>::get() {
   return MPI_2INT;
};

template<>
inline MPI_Datatype get_mpi_locdatatype_t<short>::get() {
   return MPI_SHORT_INT;
};

template<>
inline MPI_Datatype get_mpi_locdatatype_t<long double>::get() {
   return MPI_LONG_DOUBLE_INT;
};

template<typename T>
inline MPI_Datatype get_mpi_locdatatype(const T&) {
   return get_mpi_locdatatype_t<T>::get();
};

template<typename T>
inline MPI_Datatype get_mpi_locdatatype(const T* const) {
   return get_mpi_locdatatype_t<T>::get();
};

template<typename T>
inline MPI_Datatype get_mpi_locdatatype(T*) {
   return get_mpi_locdatatype_t<T>::get();
};


inline int PIPS_MPIgetRank(MPI_Comm mpiComm = MPI_COMM_WORLD) {
   if (mpiComm == MPI_COMM_NULL)
      return 0;
   int myrank;
   MPI_Comm_rank(mpiComm, &myrank);
   return myrank;
}

template<typename T>
inline std::vector<std::pair<T, int> > PIPS_MPIminlocArray(const T* localmin, int length, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   if (length <= 0)
      return std::vector<std::pair<T, int> >();

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   std::vector<std::pair<T, int> > pairs(length);

   for (unsigned int i = 0; i < pairs.size(); ++i) {
      pairs[i].first = localmin[i];
      pairs[i].second = my_rank;
   }

   MPI_Allreduce(MPI_IN_PLACE, &pairs[0], length, get_mpi_locdatatype(localmin), MPI_MINLOC, mpiComm);

   return pairs;
}

template<typename T>
inline std::vector<std::pair<T, int> > PIPS_MPIminlocArray(const std::vector<T>& localmin, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   if (localmin.size() == 0)
      return std::vector<std::pair<T, int> >();

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   std::vector<std::pair<T, int> > pairs(localmin.size());

   for (unsigned int i = 0; i < pairs.size(); ++i) {
      pairs[i].first = localmin[i];
      pairs[i].second = my_rank;
   }

   MPI_Allreduce(MPI_IN_PLACE, &pairs[0], localmin.size(), get_mpi_locdatatype(localmin[0]), MPI_MINLOC, mpiComm);

   return pairs;
}

template<typename T>
inline std::vector<std::pair<T, int> > PIPS_MPImaxlocArray(const T* localmax, int length, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   if (length <= 0)
      return std::vector<std::pair<T, int> >();

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   std::vector<std::pair<T, int> > pairs(length);

   for (unsigned int i = 0; i < pairs.size(); ++i) {
      pairs[i].first = localmax[i];
      pairs[i].second = my_rank;
   }

   MPI_Allreduce(MPI_IN_PLACE, &pairs[0], length, get_mpi_locdatatype(localmax), MPI_MAXLOC, mpiComm);

   return pairs;
}

template<typename T>
inline std::vector<std::pair<T, int> > PIPS_MPImaxlocArray(const std::vector<T>& localmax, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   if (localmax.size() == 0)
      return std::vector<std::pair<T, int> >();

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   std::vector<std::pair<T, int> > pairs(localmax.size());

   for (unsigned int i = 0; i < pairs.size(); ++i) {
      pairs[i].first = localmax[i];
      pairs[i].second = my_rank;
   }

   MPI_Allreduce(MPI_IN_PLACE, &pairs[0], localmax.size(), get_mpi_locdatatype(localmax[0]), MPI_MAXLOC, mpiComm);

   return pairs;
}

template<typename T>
inline T PIPS_MPIgetMin(const T& localmin, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   T globalmin = 0.0;
   MPI_Allreduce(&localmin, &globalmin, 1, get_mpi_datatype(localmin), MPI_MIN, mpiComm);

   return globalmin;
}

template<typename T>
inline T PIPS_MPIgetMax(const T& localmax, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   T globalmax = 0.0;
   MPI_Allreduce(&localmax, &globalmax, 1, get_mpi_datatype(localmax), MPI_MAX, mpiComm);

   return globalmax;
}

template<typename T>
void PIPS_MPIgetMaxInPlace(T& max, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   MPI_Allreduce(MPI_IN_PLACE, &max, 1, get_mpi_datatype(max), MPI_MAX, mpiComm);
}

template<typename T>
void PIPS_MPIgetMinInPlace(T& min, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   MPI_Allreduce(MPI_IN_PLACE, &min, 1, get_mpi_datatype(min), MPI_MIN, mpiComm);
}

template<typename T>
inline bool PIPS_MPIisValueEqual(const T& val, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   std::vector<T> minmax{val, -val};
   PIPS_MPImaxArrayInPlace(minmax, mpiComm);

   return (minmax[0] == -minmax[1]);
}

template<typename T>
inline void PIPS_MPImaxArrayInPlace(T* localmax, int length, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   assert(length >= 0);
   if (length == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, localmax, length, get_mpi_datatype(localmax), MPI_MAX, mpiComm);
}

template<typename T>
inline void PIPS_MPImaxArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   if (elements.size() == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_MAX, mpiComm);
}

template<typename T>
inline void PIPS_MPIminArrayInPlace(T* elements, int length, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   assert(length >= 0);
   if (length == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, elements, length, get_mpi_datatype(elements), MPI_MIN, mpiComm);
}

template<typename T>
inline void PIPS_MPIminArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   if (elements.size() == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_MIN, mpiComm);
}

template<typename T>
inline T PIPS_MPIgetSum(const T& localsummand, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   T sum;
   MPI_Allreduce(&localsummand, &sum, 1, get_mpi_datatype(localsummand), MPI_SUM, mpiComm);

   return sum;
}

template<typename T>
inline void PIPS_MPIgetLogicOrInPlace(T& localval, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   MPI_Allreduce(MPI_IN_PLACE, &localval, 1, get_mpi_datatype(localval), MPI_LOR, mpiComm);
}

template<typename T>
inline T PIPS_MPIgetLogicOr(const T& localval, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   T lor;
   MPI_Allreduce(&localval, &lor, 1, get_mpi_datatype(localval), MPI_LOR, mpiComm);

   return lor;
}

template<typename T>
inline void PIPS_MPIgetLogicAndInPlace(T& localval, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   MPI_Allreduce(MPI_IN_PLACE, &localval, 1, get_mpi_datatype(localval), MPI_LAND, mpiComm);
}

template<typename T>
inline T PIPS_MPIgetLogicAnd(const T& localval, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   T land;
   MPI_Allreduce(&localval, &land, 1, get_mpi_datatype(localval), MPI_LAND, mpiComm);

   return land;
}

template<typename T>
inline void PIPS_MPIlogicOrArrayInPlace(T* elements, int length, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   assert(length >= 0);
   if (length == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, elements, length, get_mpi_datatype(elements), MPI_LOR, mpiComm);
}

template<typename T>
inline void PIPS_MPIlogicOrArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   if (elements.size() == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_LOR, mpiComm);
}

template<typename T>
inline void PIPS_MPIgetSumInPlace(T& sum, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   MPI_Allreduce(MPI_IN_PLACE, &sum, 1, get_mpi_datatype(sum), MPI_SUM, mpiComm);
}

template<typename T>

inline void PIPS_MPIsumArrayInPlace(T* elements, int length, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   assert(length >= 0);

   if (length == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, elements, length, get_mpi_datatype(elements), MPI_SUM, mpiComm);
}

template<typename T>
inline void PIPS_MPIsumArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   if (elements.size() == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_SUM, mpiComm);
}

template<typename T>
inline void PIPS_MPIsumArray(const T* source, T* dest, int length, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   assert(length >= 0);

   if (length == 0)
      return;

   MPI_Allreduce(source, dest, length, get_mpi_datatype(source), MPI_SUM, mpiComm);
}

template<typename T>
inline void PIPS_MPImaxArray(const T* source, T* dest, int length, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   assert(length >= 0);
   if (length == 0)
      return;
   MPI_Allreduce(source, dest, length, get_mpi_datatype(source), MPI_MAX, mpiComm);
}

template<typename T>
inline void PIPS_MPImaxArray(const std::vector<T>& source, std::vector<T>& dest, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   assert(source.size() <= dest.size());
   if (source.size() == 0)
      return;

   MPI_Allreduce(&source[0], &dest[0], source.size(), get_mpi_datatype(&source[0]), MPI_MAX, mpiComm);
}

template<typename T>
inline void PIPS_MPIminArray(const std::vector<T>& source, std::vector<T>& dest, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   assert(source.size() <= dest.size());
   if (source.size() == 0)
      return;

   MPI_Allreduce(&source[0], &dest[0], source.size(), get_mpi_datatype(&source[0]), MPI_MIN, mpiComm = MPI_COMM_WORLD);
}

template<typename T>
inline void PIPS_MPIallgather(const T* vec, int n_vec, int& n_res, T*& res, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   const int size = PIPS_MPIgetSize(mpiComm);
   if (size > 1) {
      std::vector<int> recieve_sizes(size);
      std::vector<int> recieve_offsets(size);

      MPI_Allgather(&n_vec, 1, MPI_INT, recieve_sizes.data(), 1, MPI_INT, mpiComm);

      recieve_offsets[0] = 0;
      for (int i = 1; i < size; ++i)
         recieve_offsets[i] = recieve_offsets[i - 1] + recieve_sizes[i - 1];

      n_res = recieve_offsets.back() + recieve_sizes.back();

      res = new T[n_res];

      MPI_Allgatherv(vec, n_vec, get_mpi_datatype(vec), res, recieve_sizes.data(), recieve_offsets.data(), get_mpi_datatype(res), mpiComm);
   }
   else {
      n_res = n_vec;
      res = new T[n_res];
      std::uninitialized_copy(vec, vec + n_vec, res);
   }
}

inline std::string PIPS_MPIallgatherString(const std::string& str, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   const char* send = str.c_str();
   char* res{};
   int n_res{0};

   PIPS_MPIallgather(send, str.size(), n_res, res, mpiComm);

   std::string str_gathered(res, n_res);
   delete[] res;
   return str_gathered;
}

template<typename T>
inline void
PIPS_MPIgatherv(const T* sendbuf, int sendcnt, T* recvbuf, int* recvcnts, const int* recvoffsets, int root, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   assert(sendcnt >= 0);
   MPI_Gatherv(sendbuf, sendcnt, get_mpi_datatype(sendbuf), recvbuf, recvcnts, recvoffsets, get_mpi_datatype(recvbuf), root,
         mpiComm = MPI_COMM_WORLD);
}

template<typename T>
inline void PIPS_MPIallgather(const T* sendbuf, int sendcnt, T* recvbuf, int recvcnt, MPI_Comm mpiComm = MPI_COMM_WORLD) {
   assert(sendcnt >= 0);
   assert(recvcnt >= 0);
   MPI_Allgather(sendbuf, sendcnt, get_mpi_datatype(recvbuf), recvbuf, recvcnt, get_mpi_datatype(recvbuf), mpiComm);
}

#endif
