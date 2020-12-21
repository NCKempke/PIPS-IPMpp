/*
 * QpPresolver.C
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#include "QpPresolver.h"


QpPresolver::QpPresolver(const Data* prob, Postsolver* postsolver)
 : origprob(prob), postsolver(postsolver)
{
}
