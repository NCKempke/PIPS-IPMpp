/*
 * QpPresolver.C
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#include "Presolver.hpp"


Presolver::Presolver(const Problem& problem, Postsolver* postsolver) : original_problem{problem}, postsolver{postsolver} {}
