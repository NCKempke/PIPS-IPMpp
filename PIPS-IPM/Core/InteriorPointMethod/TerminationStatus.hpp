/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include <iostream>

#ifndef TERMINATIONSTATUS_H
#define TERMINATIONSTATUS_H

enum class TerminationStatus {
   SUCCESSFUL_TERMINATION = 0, NOT_FINISHED, MAX_ITS_EXCEEDED, INFEASIBLE, UNKNOWN, DID_NOT_RUN
};

std::ostream& operator<<(std::ostream& os, TerminationStatus status);
#endif










