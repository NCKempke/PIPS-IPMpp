//
// Created by nils-christian on 18.06.21.
//

#ifndef PIPSIPMPP_PROBLEMFACTORY_H
#define PIPSIPMPP_PROBLEMFACTORY_H

#include <memory>

class QP;

class Variables;

class Residuals;

class Variables;

class DoubleLinearSolver;

class AbstractMatrix;

class Problem;

class Residuals;

class AbstractLinearSystem;

class ProblemFactory {
public:
   ProblemFactory() = default;

   virtual ~ProblemFactory() = default;

   [[nodiscard]] virtual std::unique_ptr<Problem> make_problem() const = 0;

   [[nodiscard]] virtual std::unique_ptr<Residuals> make_residuals(const Problem& problem) const = 0;

   [[nodiscard]] virtual std::unique_ptr<Variables> make_variables(const Problem& problem) const = 0;

   [[nodiscard]] virtual std::unique_ptr<AbstractLinearSystem> make_linear_system(Problem& problem) const = 0;

   /** create x vector */
   [[nodiscard]] virtual std::unique_ptr<Vector<double>> make_primal_vector() const = 0;

   /** create dual A vector */
   [[nodiscard]] virtual std::unique_ptr<Vector<double>> make_equalities_dual_vector() const = 0;

   /** create dual C vector */
   [[nodiscard]] virtual std::unique_ptr<Vector<double>> make_inequalities_dual_vector() const = 0;

   /** create rhs for whole system */
   [[nodiscard]] virtual std::unique_ptr<Vector<double>> make_right_hand_side() const = 0;
};

#endif //PIPSIPMPP_PROBLEMFACTORY_H
