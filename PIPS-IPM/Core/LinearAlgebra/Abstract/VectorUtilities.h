/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef VECTORUTILITIES
#define VECTORUTILITIES

#include <iostream>
#include <fstream>
#include <cstring>

template<typename T>
T
find_blocking(const T w[], int n, int incw, const T wstep[], int incwstep, const T u[], int incu, const T ustep[], int incustep, T maxStep, T* w_elt,
      T* wstep_elt, T* u_elt, T* ustep_elt, int& first_or_second) {
   T bound = maxStep;

   int i = n - 1, lastBlocking = -1;

   // Search backward so that we find the blocking constraint of lowest
   // index. We do this to make things consistent with MPI's MPI_MINLOC,
   // which returns the processor with smallest rank where a min occurs.
   //
   // Still, going backward is ugly!
   const T* pw = w + (n - 1) * incw;
   const T* pwstep = wstep + (n - 1) * incwstep;
   const T* pu = u + (n - 1) * incu;
   const T* pustep = ustep + (n - 1) * incustep;

   while (i >= 0) {
      T temp = *pwstep;
      if (*pw > 0 && temp < 0) {
         temp = -*pw / temp;
         if (temp <= bound) {
            bound = temp;
            lastBlocking = i;
            first_or_second = 1;
         }
      }
      temp = *pustep;
      if (*pu > 0 && temp < 0) {
         temp = -*pu / temp;
         if (temp <= bound) {
            bound = temp;
            lastBlocking = i;
            first_or_second = 2;
         }
      }

      i--;
      if (i >= 0) {
         // It is safe to decrement the pointers
         pw -= incw;
         pwstep -= incwstep;
         pu -= incu;
         pustep -= incustep;
      }
   }

   if (lastBlocking > -1) {
      // fill out the elements
      *w_elt = w[lastBlocking];
      *wstep_elt = wstep[lastBlocking];
      *u_elt = u[lastBlocking];
      *ustep_elt = ustep[lastBlocking];
   }
   return bound;
}

template<typename T>
void
find_blocking_pd(const T w[], const int n, const T wstep[], const T u[], const T ustep[], T& maxStep_primal, T& maxStep_dual, T& w_elt, T& wstep_elt,
      T& u_elt, T& ustep_elt, T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d, bool& primalBlocking, bool& dualBlocking) {
   int lastBlockingPrimal = -1, lastBlockingDual = -1;

   for (int i = 0; i < n; i++) {
      T temp = wstep[i];
      if (w[i] > 0 && temp < 0) {
         temp = -w[i] / temp;
         if (temp < maxStep_primal) {
            maxStep_primal = temp;
            lastBlockingPrimal = i;
            primalBlocking = true;
         }
      }
      temp = ustep[i];
      if (u[i] > 0 && temp < 0) {
         temp = -u[i] / temp;
         if (temp < maxStep_dual) {
            maxStep_dual = temp;
            lastBlockingDual = i;
            dualBlocking = true;
         }
      }
   }

   if (lastBlockingPrimal > -1) {
      // fill out the elements
      w_elt = w[lastBlockingPrimal];
      wstep_elt = wstep[lastBlockingPrimal];
      u_elt = u[lastBlockingPrimal];
      ustep_elt = ustep[lastBlockingPrimal];
   }
   if (lastBlockingDual > -1) {
      // fill out the elements
      w_elt_d = w[lastBlockingDual];
      wstep_elt_d = wstep[lastBlockingDual];
      u_elt_d = u[lastBlockingDual];
      ustep_elt_d = ustep[lastBlockingDual];
   }
}

#endif
