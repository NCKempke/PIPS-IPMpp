/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */
#include "Vector.hpp"

template<typename T>
Vector<T>::Vector(int n_) : n{n_} {
   assert(n_ >= 0);
}

template
class Vector<int>;

template
class Vector<double>;
