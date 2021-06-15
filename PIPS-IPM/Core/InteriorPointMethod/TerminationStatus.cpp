//
// Created by nils-christian on 09.06.21.
//
#include "TerminationStatus.hpp"

std::ostream& operator<<(std::ostream& os, TerminationStatus status) {
   switch (status) {
      case TerminationStatus::SUCCESSFUL_TERMINATION:
         os << "SUCCESSFUL_TERMINATION";
         break;
      case TerminationStatus::NOT_FINISHED:
         os << "NOT_FINISHED";
         break;
      case TerminationStatus::MAX_ITS_EXCEEDED:
         os << "MAX_ITS_EXCEEDED";
         break;
      case TerminationStatus::INFEASIBLE:
         os << "INFEASIBLE";
         break;
      case TerminationStatus::UNKNOWN:
         os << "UNKNOWN";
         break;
      case TerminationStatus::DID_NOT_RUN:
         os << "DID_NOT_RUN";
         break;
   }
   return os;
}
