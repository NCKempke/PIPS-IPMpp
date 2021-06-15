/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENDATA
#define QPGENDATA

#include "Problem.hpp"
#include "Vector.hpp"

class MpsReader;

class Variables;

#ifdef TESTING
class QpGenDataTester;
#endif

/**
 * Data for the general QP formulation.
 *
 * @ingroup QpGen
 */

class QP : public Problem {
#ifdef TESTING
   friend QpGenDataTester;
#endif

protected:

   QP() = default;

public:
   std::shared_ptr<SymmetricMatrix> Q;

   /** constructor that sets up pointers to the data objects that are
       passed as arguments */
   QP(std::shared_ptr<Vector<double>> c, std::shared_ptr<SymmetricMatrix> Q, std::shared_ptr<Vector<double>> xlow,
      std::shared_ptr<Vector<double>> ixlow, std::shared_ptr<Vector<double>> xupp,
      std::shared_ptr<Vector<double>> ixupp, std::shared_ptr<GeneralMatrix> A, std::shared_ptr<Vector<double>> bA,
      std::shared_ptr<GeneralMatrix> C, std::shared_ptr<Vector<double>> clow, std::shared_ptr<Vector<double>> iclow,
      std::shared_ptr<Vector<double>> cupp,
      std::shared_ptr<Vector<double>> icupp);

   /** insert the Hessian Q into the matrix M for the fundamental linear system, where M is stored as a SymMatrix */
   virtual void putQIntoAt(SymmetricMatrix& M, int row, int col);

   /** insert the Hessian Q into the matrix M for the fundamental linear system, where M is stored as a GenMatrix */
   virtual void putQIntoAt(GeneralMatrix& M, int row, int col);

   /** y = beta * y + alpha * Q * x */
   void
   hessian_multiplication(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   /** extract the diagonal of the Hessian and put it in the Vector<double> hessian_diagonal */
   void hessian_diagonal(Vector<double>& hessian_diagonal) const override;

   void createScaleFromQ();

   void scaleQ();

   void flipQ();

   [[nodiscard]] double datanorm() const override;

   virtual void datainput() {};

   virtual void datainput(MpsReader* reader, int& iErr);

   void print() override;

   virtual void objective_gradient(const Variables& variables, Vector<double>& gradient) const override;

   virtual double objective_value(const Variables& variables) const override;

   ~QP() override = default;
};

#endif
