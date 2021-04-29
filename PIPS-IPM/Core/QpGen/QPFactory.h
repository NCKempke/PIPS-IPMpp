/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SPSEQQPGENFACTORY
#define SPSEQQPGENFACTORY

#include "ProblemFactory.h"

class QP;

class Variables;

class QPFactory : public ProblemFactory {
protected:
   int nnzQ;
   int nnzA;
   int nnzC;
public:
   QPFactory(int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_) : ProblemFactory(nx_, my_, mz_), nnzQ(nnzQ_), nnzA(nnzA_), nnzC(nnzC_) {}

   Problem* create_problem(double c[], int krowQ[], int jcolQ[], double dQ[], double xlow[], char ixlow[], double xupp[], char ixupp[], int krowA[],
         int jcolA[], double dA[], double b[], int krowC[], int jcolC[], double dC[], double clow[], char iclow[], double cupp[], char icupp[]);

   void join_right_hand_side(OoqpVector& rhs_in, const OoqpVector& rhs1_in, const OoqpVector& rhs2_in, const OoqpVector& rhs3_in) const override;

   void separate_variables(OoqpVector& x_in, OoqpVector& y_in, OoqpVector& z_in, const OoqpVector& vars_in) const override;
};

#endif
