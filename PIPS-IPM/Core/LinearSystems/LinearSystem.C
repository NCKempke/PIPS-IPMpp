/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "LinearSystem.h"
#include "Residuals.h"
#include "Problem.h"
#include "Variables.h"
#include "Vector.hpp"
#include "SmartPointer.h"
#include "mpi.h"
#include "QpGenOptions.h"
#include "DistributedOptions.h"
#include "DistributedFactory.h"
#include <vector>
#include <functional>
#include <type_traits>
#include <memory>

extern int gOuterBiCGIter;
extern double gOuterBiCGIterAvg;
extern double g_iterNumber;
extern int gOuterBiCGFails;

static std::vector<int> bicgIters;

LinearSystem::LinearSystem(DistributedFactory* factory_, Problem* problem, bool create_iter_ref_vecs) : factory(factory_),
      apply_regularization(qpgen_options::getBoolParameter("REGULARIZATION")), outerSolve(qpgen_options::getIntParameter("OUTER_SOLVE")),
      innerSCSolve(qpgen_options::getIntParameter("INNER_SC_SOLVE")),
      outer_bicg_print_statistics(qpgen_options::getBoolParameter("OUTER_BICG_PRINT_STATISTICS")),
      outer_bicg_eps(qpgen_options::getDoubleParameter("OUTER_BICG_EPSILON")),
      outer_bicg_max_iter(qpgen_options::getIntParameter("OUTER_BICG_MAX_ITER")),
      outer_bicg_max_normr_divergences(qpgen_options::getIntParameter("OUTER_BICG_MAX_NORMR_DIVERGENCES")),
      outer_bicg_max_stagnations(qpgen_options::getIntParameter("OUTER_BICG_MAX_STAGNATIONS")),
      xyzs_solve_print_residuals(qpgen_options::getBoolParameter("XYZS_SOLVE_PRINT_RESISDUAL")) {

   assert(factory_);
   assert(problem);

   nx = problem->nx;
   my = problem->my;
   mz = problem->mz;
   ixlow = problem->ixlow;
   ixupp = problem->ixupp;
   iclow = problem->iclow;
   icupp = problem->icupp;

   nxlow = problem->nxlow;
   nxupp = problem->nxupp;
   mclow = problem->mclow;
   mcupp = problem->mcupp;

   if (create_iter_ref_vecs) {
      if (outerSolve || xyzs_solve_print_residuals) {
         //for iterative refinement or BICGStab
         sol = factory->make_right_hand_side();
         sol2 = factory->make_right_hand_side();
         res = factory->make_right_hand_side();
         resx = factory->make_primal_vector();
         resy = factory->make_equalities_dual_vector();
         resz = factory->make_inequalities_dual_vector();

         if (outerSolve == 2) {
            //BiCGStab; additional vectors needed
            sol3 = factory->make_right_hand_side();
            res2 = factory->make_right_hand_side();
            res3 = factory->make_right_hand_side();
            res4 = factory->make_right_hand_side();
            res5 = factory->make_right_hand_side();
         }
      }
   }
}

LinearSystem::LinearSystem(DistributedFactory* factory_, Problem* problem, Vector<double>* primal_diagonal_, Vector<double>* dq_,
      Vector<double>* nomegaInv_, Vector<double>* primal_regularization_, Vector<double>* dual_equality_regularization_,
      Vector<double>* dual_inequality_regularization_, Vector<double>* rhs_, bool create_iter_ref_vecs) : LinearSystem(factory_, problem,
      create_iter_ref_vecs) {
   primal_diagonal = primal_diagonal_;
   dq = dq_;
   nomegaInv = nomegaInv_;
   primal_regularization_diagonal = primal_regularization_;
   dual_equality_regularization_diagonal = dual_equality_regularization_;
   dual_inequality_regularization_diagonal = dual_inequality_regularization_;
   rhs = rhs_;
}

LinearSystem::LinearSystem(DistributedFactory* factory_, Problem* problem) : LinearSystem(factory_, problem, true) {
   if (nxupp + nxlow > 0) {
      primal_diagonal = factory->make_primal_vector();
      dq = factory->make_primal_vector();
      problem->hessian_diagonal(*dq);
   }

   nomegaInv = factory->make_inequalities_dual_vector();
   rhs = factory->make_right_hand_side();


   primal_regularization_diagonal = factory->make_primal_vector();
   dual_equality_regularization_diagonal = factory->make_equalities_dual_vector();
   dual_inequality_regularization_diagonal = factory->make_inequalities_dual_vector();
}

LinearSystem::~LinearSystem() {
   if (!useRefs) {
      delete dual_inequality_regularization_diagonal;
      delete dual_equality_regularization_diagonal;
      delete primal_regularization_diagonal;
      delete primal_diagonal;
      delete dq;
      delete rhs;
      delete nomegaInv;
   }

   delete sol;
   delete res;
   delete resx;
   delete resy;
   delete resz;
   delete sol2;
   delete sol3;
   delete res2;
   delete res3;
   delete res4;
   delete res5;
}

int LinearSystem::getIntValue(const std::string& s) const {
   if (s == "BICG_NITERATIONS")
      return bicg_niterations;
   else if (s == "BICG_CONV_FLAG")
      return static_cast<std::underlying_type<IterativeSolverSolutionStatus>::type>(bicg_conv_flag);
   else {
      std::cout << "Unknown observer int request in LinearSystem.C: " << s << "\n";
      return -1;
   }
}

bool LinearSystem::getBoolValue(const std::string& s) const {
   if (s == "BICG_CONVERGED")
      return bicg_conv_flag == IterativeSolverSolutionStatus::CONVERGED;
   else if (s == "BICG_SKIPPED")
      return bicg_conv_flag == IterativeSolverSolutionStatus::SKIPPED;
   else if (s == "BICG_DIVERGED")
      return bicg_conv_flag == IterativeSolverSolutionStatus::DIVERGED;
   else if (s == "BICG_BREAKDOWN")
      return bicg_conv_flag == IterativeSolverSolutionStatus::BREAKDOWN;
   else if (s == "BICG_STAGNATION")
      return bicg_conv_flag == IterativeSolverSolutionStatus::STAGNATION;
   else if (s == "BICG_EXCEED_MAX_ITER")
      return bicg_conv_flag == IterativeSolverSolutionStatus::NOT_CONVERGED_MAX_ITERATIONS;
   else {
      std::cout << "Unknown observer bool request in LinearSystem.C: " << s << "\n";
      return false;
   }
}

double LinearSystem::getDoubleValue(const std::string& s) const {
   if (s == "BICG_RESNORM")
      return bicg_resnorm;
   else if (s == "BICG_RELRESNORM")
      return bicg_relresnorm;
   else {
      std::cout << "Unknown observer double request in LinearSystem.C: " << s << "\n";
      return 0.0;
   }
}

static void biCGStabCommunicateStatus(int flag, int it) {
   double iterAvg = 0.0;

   /* IP algorithm started? */
   if (g_iterNumber >= 0.5) {
      bicgIters.push_back(it);

      for (size_t i = 0; i < bicgIters.size(); i++)
         iterAvg += double(bicgIters[i]);

      iterAvg /= bicgIters.size();
   }
   else {
      iterAvg = it;
   }

   gOuterBiCGIterAvg = iterAvg;
   gOuterBiCGIter = it;

   if (flag != 0 && flag != 1)
      gOuterBiCGFails++;
}

static bool isZero(double val, LinearSystem::IterativeSolverSolutionStatus& status) {
   if (PIPSisZero(val)) {
      status = LinearSystem::IterativeSolverSolutionStatus::BREAKDOWN;
      return true;
   }

   return false;
}

void LinearSystem::factorize(Problem* /* problem */, Variables* vars) {

   assert(vars->validNonZeroPattern());
   assert(vars->validNonZeroPattern());

   put_barrier_parameter(vars->mu());

   computeDiagonals(*vars->t, *vars->lambda, *vars->u, *vars->pi, *vars->v, *vars->gamma, *vars->w, *vars->phi);

   if (pips_options::get_bool_parameter("HIERARCHICAL_TESTING")) {
      std::cout << "Setting diags to 1.0 for Hierarchical debugging\n";
      primal_diagonal->setToConstant(1.0);
      nomegaInv->setToConstant(1.0);
   }

   if (nxlow + nxupp > 0) {
      put_primal_diagonal();
   }

   clear_dual_equality_diagonal();

   if (mclow + mcupp > 0) {
      put_dual_inequalites_diagonal();
   }

   primal_regularization_diagonal->setToZero();
   dual_equality_regularization_diagonal->setToZero();
   dual_inequality_regularization_diagonal->setToZero();

   printDiagonalNorms();
}

void LinearSystem::printDiagonalNorms() const {
   assert(primal_diagonal);
   assert(nomegaInv);
   const double infnorm_primal_diagonal = primal_diagonal->infnorm();
   const double twonorm_primal_diagonal = primal_diagonal->twonorm();
   const double infnorm_inequalities_diagonal = nomegaInv->infnorm();
   const double twonorm_inequalities_diagonal = nomegaInv->twonorm();

   double min_primal_diagonal;
   int dummy;
   primal_diagonal->min(min_primal_diagonal, dummy);

   double min_inequalities_diagonal;
   primal_diagonal->min(min_inequalities_diagonal, dummy);

   if (PIPS_MPIgetRank() == 0) {
      std::cout << "Primal diagonal: ||.||_inf=" << infnorm_primal_diagonal << ", ||.||_2=" << twonorm_primal_diagonal << ", min="
                << min_primal_diagonal << "\n";
      std::cout << "Inequalities diagonal: ||.||_inf=" << infnorm_inequalities_diagonal << ", ||.||_2=" << twonorm_inequalities_diagonal << ", min "
                << min_inequalities_diagonal << "\n";
   }
}

void LinearSystem::print_regularization_statistics() const {
   assert(primal_regularization_diagonal);
   assert(dual_equality_regularization_diagonal);
   assert(dual_inequality_regularization_diagonal);

   const double infnorm_primal_regularization = primal_regularization_diagonal->infnorm();
   const double infnorm_dual_equality_regularization = dual_equality_regularization_diagonal->infnorm();
   const double infnorm_dual_inequality_regularization = dual_inequality_regularization_diagonal->infnorm();

   if (PIPS_MPIgetRank() == 0) {
      std::cout << "Regularized system with (||dP||_inf, ||dDy||_inf, ||dDz||_inf) = (" << infnorm_primal_regularization << ", "
                << infnorm_dual_equality_regularization << ", " << infnorm_dual_inequality_regularization << ")\n";
   }
}

void LinearSystem::computeDiagonals(Vector<double>& t, Vector<double>& lambda, Vector<double>& u, Vector<double>& pi, Vector<double>& v,
      Vector<double>& gamma, Vector<double>& w, Vector<double>& phi) {

   /*** dd = dQ + Gamma/V + Phi/W ***/
   if (nxlow + nxupp > 0) {
      primal_diagonal->copyFrom(*dq);
   }
   if (nxupp + nxlow > 0) {
      if (nxlow > 0)
         primal_diagonal->axdzpy(1.0, gamma, v, *ixlow);
      if (nxupp > 0)
         primal_diagonal->axdzpy(1.0, phi, w, *ixupp);
   }
   assert(primal_diagonal->allOf([](const double& d) {
      return d >= 0;
   }));

   nomegaInv->setToZero();
   /*** omega = Lambda/T + Pi/U ***/
   if (mclow > 0)
      nomegaInv->axdzpy(1.0, lambda, t, *iclow);
   if (mcupp > 0)
      nomegaInv->axdzpy(1.0, pi, u, *icupp);

   assert(nomegaInv->allOf([](const double& d) {
      return d >= 0;
   }));

   /*** omega = -omega^-1 ***/
   nomegaInv->invert();
   nomegaInv->negate();
}

void LinearSystem::solve(Problem* problem, Variables* variables, Residuals* residuals, Variables* step) {
   assert(variables->validNonZeroPattern());
   assert(residuals->valid_non_zero_pattern());

   /*** compute rX ***/
   /* rx = rQ */
   step->x->copyFrom(*residuals->lagrangian_gradient);
   if (nxlow > 0) {
      Vector<double>& gamma_by_v = *step->v;
      gamma_by_v.copyFrom(*variables->gamma);
      gamma_by_v.divideSome(*variables->v, *ixlow);

      /* rx = rQ + Gamma/V rv */
      step->x->axzpy(1.0, gamma_by_v, *residuals->rv);
      /* rx = rQ + Gamma/V rv + rGamma/V */
      step->x->axdzpy(1.0, *residuals->rgamma, *variables->v, *ixlow);
   }

   if (nxupp > 0) {
      Vector<double>& phi_by_w = *step->w;
      phi_by_w.copyFrom(*variables->phi);
      phi_by_w.divideSome(*variables->w, *ixupp);

      /* rx = rQ + Gamma/V * rv + rGamma/V + Phi/W * rw */
      step->x->axzpy(1.0, phi_by_w, *residuals->rw);
      /* rx = rQ + Gamma/V * rv + rGamma/V + Phi/W * rw - rphi/W */
      step->x->axdzpy(-1.0, *residuals->rphi, *variables->w, *ixupp);
   }

   // start by partially computing step->s
   /*** compute rs ***/
   /* step->s = rz */
   step->s->copyFrom(*residuals->rz);
   if (mclow > 0) {
      Vector<double>& lambda_by_t = *step->t;
      lambda_by_t.copyFrom(*variables->lambda);
      lambda_by_t.divideSome(*variables->t, *iclow);

      /* step->s = rz + Lambda/T * rt */
      step->s->axzpy(1.0, lambda_by_t, *residuals->rt);
      /* step->s = rz + Lambda/T * rt + rlambda/T */
      step->s->axdzpy(1.0, *residuals->rlambda, *variables->t, *iclow);
   }

   if (mcupp > 0) {
      Vector<double>& pi_by_u = *step->u;
      pi_by_u.copyFrom(*variables->pi);
      pi_by_u.divideSome(*variables->u, *icupp);

      /* step->s = rz + Lambda/T * rt + rlambda/T + Pi/U *ru */
      step->s->axzpy(1.0, pi_by_u, *residuals->ru);
      /* step->s = rz + Lambda/T * rt + rlambda/T + Pi/U *ru - rpi/U */
      step->s->axdzpy(-1.0, *residuals->rpi, *variables->u, *icupp);
   }

   /*** ry = rA ***/
   step->y->copyFrom(*residuals->rA);
   /*** rz = rC ***/
   step->z->copyFrom(*residuals->rC);

   {
      // Unfortunately, we need a temporary Vector<double> for the solve,
      // Use step->lambda or step->pi
      SmartPointer<Vector<double> > ztemp;
      if (mclow > 0) {
         ztemp = step->lambda;
      }
      else {
         ztemp = step->pi;
      }

      solveXYZS(*step->x, *step->y, *step->z, *step->s, *ztemp, problem);
   }

   if (mclow > 0) {
      /* Dt = Ds - rt */
      step->t->copyFrom(*step->s);
      step->t->axpy(-1.0, *residuals->rt);
      step->t->selectNonZeros(*iclow);

      /* Dlambda = T^-1 (rlambda - Lambda * Dt ) */
      step->lambda->copyFrom(*residuals->rlambda);
      step->lambda->axzpy(-1.0, *variables->lambda, *step->t);
      step->lambda->divideSome(*variables->t, *iclow);
      //!
      step->lambda->selectNonZeros(*iclow);
   }

   if (mcupp > 0) {
      /* Du = ru - Ds */
      step->u->copyFrom(*residuals->ru);
      step->u->axpy(-1.0, *step->s);
      step->u->selectNonZeros(*icupp);

      /* Dpi = U^-1 ( rpi - Pi * Du ) */
      step->pi->copyFrom(*residuals->rpi);
      step->pi->axzpy(-1.0, *variables->pi, *step->u);
      step->pi->divideSome(*variables->u, *icupp);
      //!
      step->pi->selectNonZeros(*icupp);
   }

   if (nxlow > 0) {
      /* Dv = Dx - rv */
      step->v->copyFrom(*step->x);
      step->v->axpy(-1.0, *residuals->rv);
      step->v->selectNonZeros(*ixlow);

      /* Dgamma = V^-1 ( rgamma - Gamma * Dv ) */
      step->gamma->copyFrom(*residuals->rgamma);
      step->gamma->axzpy(-1.0, *variables->gamma, *step->v);
      step->gamma->divideSome(*variables->v, *ixlow);
      //!
      step->gamma->selectNonZeros(*ixlow);
   }

   if (nxupp > 0) {
      /* Dw = rw - Dx */
      step->w->copyFrom(*residuals->rw);
      step->w->axpy(-1.0, *step->x);
      step->w->selectNonZeros(*ixupp);

      /* Dphi = W^-1 ( rphi - Phi * Dw ) */
      step->phi->copyFrom(*residuals->rphi);
      step->phi->axzpy(-1.0, *variables->phi, *step->w);
      step->phi->divideSome(*variables->w, *ixupp);
      //!
      step->phi->selectNonZeros(*ixupp);
   }
   assert(step->validNonZeroPattern());

}

void LinearSystem::solveXYZS(Vector<double>& stepx, Vector<double>& stepy, Vector<double>& stepz, Vector<double>& steps, Vector<double>& /* ztemp */,
      Problem* problem) {
   /* step->z = rC */
   /* step->s = rz + Lambda/T * rt + rlambda/T + Pi/U *ru - rpi/U */

   /* rx = rQ + Gamma/V * rv + rGamma/V + Phi/W * rw - rphi/W */
   /* ry = rA */
   /* rz = rC + Omega^-1 ( rz + Lambda/T * rt + rlambda/T + Pi/U *ru - rpi/U ) */
   stepz.axzpy(-1.0, *nomegaInv, steps);

   std::unique_ptr<Vector<double>> residual;
   if (xyzs_solve_print_residuals) {
      residual.reset(rhs->cloneFull());
      joinRHS(*residual, stepx, stepy, stepz);

      const double xinf = stepx.infnorm();
      const double yinf = stepy.infnorm();
      const double zinf = stepz.infnorm();

      if (PIPS_MPIgetRank() == 0)
         std::cout << "rhsx norm : " << xinf << ",\trhsy norm : " << yinf << ",\trhsz norm : " << zinf << "\n";
   }

   assert(rhs);
   this->joinRHS(*rhs, stepx, stepy, stepz);

   if (outerSolve == 1) {
      ///////////////////////////////////////////////////////////////
      // Iterative refinement
      ///////////////////////////////////////////////////////////////
      auto computeResiduals = [this, &stepx, &stepy, &stepz, &capture0 = *problem](auto&& PH1, auto&& PH2) {
         compute_regularized_system_residuals(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2), stepx, stepy, stepz, capture0);
      };

      solveCompressedIterRefin(computeResiduals);

      this->separateVars(stepx, stepy, stepz, *sol);

   }
   else if (outerSolve == 0) {
      ///////////////////////////////////////////////////////////////
      // Default solve - Schur complement based decomposition
      ///////////////////////////////////////////////////////////////
      solveCompressed(*rhs);
      separateVars(stepx, stepy, stepz, *rhs);

   }
   else {
      assert(outerSolve == 2);
      ///////////////////////////////////////////////////////////////
      // BiCGStab
      ///////////////////////////////////////////////////////////////

      const bool use_regularized_system = true;
      auto matMult = [this, &capture0 = *problem, &stepx, &stepy, &stepz, use_regularized_system](auto&& PH1, auto&& PH2, auto&& PH3, auto&& PH4) {
         system_mult(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2), std::forward<decltype(PH3)>(PH3),
               std::forward<decltype(PH4)>(PH4), capture0, stepx, stepy, stepz, use_regularized_system);
      };

      auto matInfnorm = [this, &capture0 = *problem, &stepx, &stepy, &stepz] {
         return matXYZinfnorm(capture0, stepx, stepy, stepz, use_regularized_system);
      };

      solveCompressedBiCGStab(matMult, matInfnorm);

      this->separateVars(stepx, stepy, stepz, *sol);

      /* notify observers about result of BiCGStab */
      notifyObservers();
   }

   if (xyzs_solve_print_residuals) {
      assert(sol);
      const double bnorm = residual->infnorm();
      this->joinRHS(*sol, stepx, stepy, stepz);
      this->system_mult(1.0, *residual, -1.0, *sol, *problem, stepx, stepy, stepz, true);

      this->separateVars(*resx, *resy, *resz, *residual);
      const double resxnorm = resx->infnorm();
      const double resynorm = resy->infnorm();
      const double resznorm = resz->infnorm();

      if (PIPS_MPIgetRank() == 0) {
         std::cout << "bnorm " << bnorm << "\n";
         std::cout << "resx norm: " << resxnorm << "\tnorm/bnorm " << resxnorm / bnorm << "\n";
         std::cout << "resy norm: " << resynorm << "\tnorm/bnorm " << resynorm / bnorm << "\n";
         std::cout << "resz norm: " << resznorm << "\tnorm/bnorm " << resznorm / bnorm << "\n";
      }
   }

   stepy.negate();
   stepz.negate();

   /* Ds = Omega^-1 (rz + Lambda/T * rt + rlambda/T + Pi/U *ru - rpi/U - Dz ) */
   steps.axpy(-1.0, stepz);
   steps.componentMult(*nomegaInv);
   steps.negate();
}

void LinearSystem::solveCompressedBiCGStab(const std::function<void(double, Vector<double>&, double, Vector<double>&)>& matMult,
      const std::function<double()>& matInfnorm) {
   //aliases
   Vector<double>& r0 = *res2, & dx = *sol2, & best_x = *sol3, & v = *res3, & t = *res4, & p = *res5;
   Vector<double>& x = *sol, & r = *res, & b = *rhs;

   const double tol = qpgen_options::getDoubleParameter("OUTER_BICG_TOL");
   const double n2b = b.twonorm();
   const double tolb = std::max(n2b * tol, outer_bicg_eps);

   gOuterBiCGIter = 0;
   bicg_niterations = 0;

   assert(n2b >= 0);

   const int myRank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   //starting guess/point
   x.copyFrom(b);

   //solution to the approx. system
   solveCompressed(x);
   //initial residual: res = b - Ax
   r.copyFrom(b);
   matMult(1.0, r, -1.0, x);

   double normr = r.twonorm(), normr_min = normr;
   best_x.copyFrom(x);

   bicg_resnorm = normr;
   bicg_relresnorm = bicg_resnorm / n2b;

   //quick return if solve is accurate enough
   if (normr <= tolb) {

      bicg_conv_flag = IterativeSolverSolutionStatus::SKIPPED;
      biCGStabCommunicateStatus(static_cast<std::underlying_type<IterativeSolverSolutionStatus>::type>(bicg_conv_flag), bicg_niterations);

      if (myRank == 0)
         std::cout << "BiCGStab (it=" << bicg_niterations << ", rel.res.norm=" << bicg_relresnorm << ", rel.r.norm=" << normr / n2b << ", avg.iter="
                   << gOuterBiCGIterAvg << ") " << bicg_conv_flag << "\n";
      return;
   }

   if (outer_bicg_print_statistics) {
      const double infb = b.infnorm();
      const double glbinfnorm = matInfnorm();
      const double xonenorm = x.onenorm();

      if (myRank == 0) {
         std::cout << "global system infnorm=" << glbinfnorm << " x 1norm=" << xonenorm << " tolb/tolnew: " << tolb << " "
                   << (tol * xonenorm * glbinfnorm) << "\n";
         std::cout << "outer BiCGStab starts: " << normr << " > " << tolb << " normb2=" << n2b << " normbinf=" << infb << " (tolerance=" << tol << ")"
                   << "\n";
      }
   }

   r0.copyFrom(r);

   //normalize
   r0.scale(1 / normr);

   bicg_conv_flag = IterativeSolverSolutionStatus::NOT_CONVERGED_MAX_ITERATIONS;
   int normrNDiv = 0;
   int nstags = 0;
   double rho = 1., omega = 1., alpha = 1.;

   //main loop
   for (bicg_niterations = 0; bicg_niterations < outer_bicg_max_iter; bicg_niterations++) {
      const double rho1 = rho;

      rho = r0.dotProductWith(r);

      if (isZero(rho, bicg_conv_flag))
         break;

      //first half of the iterate
      {
         if (bicg_niterations == 0)
            p.copyFrom(r);
         else {
            const double beta = (rho / rho1) * (alpha / omega);

            if (isZero(beta, bicg_conv_flag))
               break;

            //-------- p = r + beta*(p - omega*v) --------
            p.axpy(-omega, v);
            p.scale(beta);
            p.axpy(1.0, r);
         }

         //precond: ph = \tilde{K}^{-1} p
         dx.copyFrom(p);
         solveCompressed(dx);

         //mat-first: v = K*ph
         matMult(0.0, v, 1.0, dx);

         const double rtv = r0.dotProductWith(v);

         if (isZero(rtv, bicg_conv_flag))
            break;

         alpha = rho / rtv;

         if ((std::fabs(alpha) * dx.twonorm()) <= outer_bicg_eps * x.twonorm())
            nstags++;
         else
            nstags = 0;

         // x = x + alpha * dx ( x = x + alpha * ph)
         x.axpy(alpha, dx);
         // r = r - alpha * v ( s = r - alpha * v)
         r.axpy(-alpha, v);

         //check for convergence
         normr = r.twonorm();

         if (normr <= tolb || nstags >= outer_bicg_max_stagnations) {
            //compute the actual residual
            Vector<double>& res = dx; //use dx
            res.copyFrom(b);
            matMult(1.0, res, -1.0, x);

            bicg_resnorm = res.twonorm();
            if (bicg_resnorm <= tolb) {
               //converged
               bicg_conv_flag = IterativeSolverSolutionStatus::CONVERGED;
               break;
            }
         } //~end of convergence test
      }

      //second half of the iterate now
      {
         //preconditioner
         dx.copyFrom(r);
         solveCompressed(dx);

         //mat-first
         matMult(0.0, t, 1.0, dx);

         const double tt = t.dotProductSelf(1.0);

         if (isZero(tt, bicg_conv_flag))
            break;

         omega = t.dotProductWith(r) / tt;

         if ((std::fabs(omega) * dx.twonorm()) <= outer_bicg_eps * x.twonorm())
            nstags++;
         else
            nstags = 0;

         // x=x+omega*dx  (x=x+omega*sh)
         x.axpy(omega, dx);
         // r = r-omega*t (r=s-omega*sh)
         r.axpy(-omega, t);
         //check for convergence
         normr = r.twonorm();

         if (normr <= tolb || nstags >= outer_bicg_max_stagnations) {
            //compute the actual residual
            Vector<double>& res = dx; //use dx
            res.copyFrom(b);
            matMult(1.0, res, -1.0, x);

            bicg_resnorm = res.twonorm();

            if (bicg_resnorm <= tolb) {
               //converged
               bicg_conv_flag = IterativeSolverSolutionStatus::CONVERGED;
               break;
            }
         }
         else {
            if (normr >= normr_min)
               normrNDiv++;
            else
               normrNDiv = 0;

            if (normrNDiv > outer_bicg_max_normr_divergences) {
               // rollback to best iterate
               x.copyFrom(best_x);
               normr = normr_min;

               bicg_conv_flag = IterativeSolverSolutionStatus::DIVERGED;
               break;
            }
         } //~end of convergence test

         if (normr < normr_min) {
            // update best for rollback
            normr_min = normr;
            best_x.copyFrom(x);
         }
      } //~end of scoping

      if (nstags >= outer_bicg_max_stagnations) {
         // rollback to best iterate
         if (normr_min < normr) {
            normr = normr_min;
            x.copyFrom(best_x);
         }

         bicg_conv_flag = IterativeSolverSolutionStatus::STAGNATION;
         break;
      }

      if (isZero(omega, bicg_conv_flag))
         break;

   } //~ end of BiCGStab loop

   bicg_relresnorm = bicg_resnorm / n2b;

   biCGStabCommunicateStatus(static_cast<std::underlying_type<IterativeSolverSolutionStatus>::type>(bicg_conv_flag), std::max(bicg_niterations, 1));
   if (myRank == 0)
      std::cout << "BiCGStab (it=" << std::max(1, bicg_niterations) << ", rel.res.norm=" << bicg_relresnorm << ", rel.r.norm=" << normr / n2b
                << ", avg.iter=" << gOuterBiCGIterAvg << ") " << bicg_conv_flag << "\n";
}

/**
 * res = beta * res + alpha * mat * sol
 *       [ Q + dq + gamma/ v + phi/w + regP                   AT                                 CT               ]
 * mat = [            A                     -dual_equality_regularization_diagonal               0                ]
 *       [            C                                       0                     -(lambda/V + pi/u)^-1 + regDz ]
 * stepx, stepy, stepz are used as temporary buffers
 * if use_regularized_sysyem == false primal_regularization_diagonal, dual_equality_regularization_diagonal and dual_inequality_regularization_diagonal are not used
 */
void
LinearSystem::system_mult(double beta, Vector<double>& res, double alpha, const Vector<double>& sol, const Problem& problem, Vector<double>& solx,
      Vector<double>& soly, Vector<double>& solz, bool use_regularized_system) {
   assert(resx);
   assert(resy);
   assert(resz);
   assert(nomegaInv);
   assert(primal_diagonal);
   assert((primal_regularization_diagonal && dual_equality_regularization_diagonal && dual_inequality_regularization_diagonal) ||
          (!primal_regularization_diagonal && !dual_equality_regularization_diagonal && !dual_inequality_regularization_diagonal));

   separateVars(solx, soly, solz, sol);
   separateVars(*resx, *resy, *resz, res);

   /* resx = beta resx + alpha Q solx + alpha dd solx + alpha primal_regularization_diagonal solx */
   problem.hessian_multiplication(beta, *resx, alpha, solx);
   resx->axzpy(alpha, *primal_diagonal, solx);
   if (use_regularized_system && primal_regularization_diagonal)
      resx->axzpy(alpha, *primal_regularization_diagonal, solx);

   /* resx = beta resx + alpha Q solx + alpha dd solx + alpha AT soly + alpha CT solz */
   problem.ATransmult(1.0, *resx, alpha, soly);
   problem.CTransmult(1.0, *resx, alpha, solz);

   /* resy = beta resy + alpha A solx */
   problem.Amult(beta, *resy, alpha, solx);
   if (use_regularized_system && dual_equality_regularization_diagonal)
      resy->axzpy(alpha, *dual_equality_regularization_diagonal, soly);

   /* resz = beta resz + alpha C solx + alpha nomegaInv solz */
   problem.Cmult(beta, *resz, alpha, solx);
   resz->axzpy(alpha, *nomegaInv, solz);
   if (use_regularized_system && dual_inequality_regularization_diagonal)
      resz->axzpy(alpha, *dual_inequality_regularization_diagonal, solz);

   this->joinRHS(res, *resx, *resy, *resz);
}

/* computes infinity norm of entire system; solx, soly, solz are used as temporary buffers */
double
LinearSystem::matXYZinfnorm(const Problem& problem, Vector<double>& solx, Vector<double>& soly, Vector<double>& solz, bool use_regularized_system) {
   solx.copyFromAbs(*primal_diagonal);
   if (use_regularized_system && primal_regularization_diagonal)
      solx.axpy(1.0, *primal_regularization_diagonal);

   problem.A->addColSums(solx);
   problem.C->addColSums(solx);
   double infnorm = solx.infnorm();

   soly.setToZero();
   if (use_regularized_system && dual_equality_regularization_diagonal)
      soly.axpy(1.0, *dual_equality_regularization_diagonal);

   problem.A->addRowSums(soly);
   infnorm = std::max(infnorm, soly.infnorm());

   solz.copyFromAbs(*nomegaInv);
   solz.negate();
   if (use_regularized_system && dual_inequality_regularization_diagonal)
      solz.axpy(1.0, *dual_inequality_regularization_diagonal);

   problem.C->addRowSums(solz);
   infnorm = std::max(infnorm, solz.infnorm());

   return infnorm;
}

void LinearSystem::solveCompressedIterRefin(const std::function<void(Vector<double>& solution, Vector<double>& residual)>& computeResidual) {
   assert(this->rhs);
   assert(this->res);
   assert(this->sol);
   assert(this->sol2);

   /* aliases */
   Vector<double>& residual = *this->res;
   Vector<double>& solution = *this->sol;
   Vector<double>& best_solution = *this->sol2;
   const Vector<double>& right_hand_side = *this->rhs;

   residual.copyFrom(right_hand_side);
   solution.setToZero();

   const double rhs_norm = std::max(1.0, residual.twonorm());
   const double tolerance_iterative_refinement = 1e-8;
   const int max_iterative_refinement_steps = 50;
   const int max_steps_stalling = 5;
   double residual_norm_best_solution = std::numeric_limits<double>::max();
   double residual_norm_current_solution{};
   double residual_last_iterate = std::numeric_limits<double>::max();

   const double min_improvement_factor = 1e-1;
   const double max_divergence_factor = 1e4;
   int n_steps_no_improvement = 0;

   IterativeSolverSolutionStatus convergence_status = IterativeSolverSolutionStatus::DID_NOT_RUN;

   int current_iterative_refinement_step = -1;

   do {
      current_iterative_refinement_step++;

      /* solve A'x = b */
      this->solveCompressed(residual);

      /* x = x + dx */
      solution.axpy(1.0, residual);

      residual.copyFrom(right_hand_side);

      computeResidual(solution, residual);
      residual_norm_current_solution = residual.twonorm();

      if (residual_norm_current_solution > min_improvement_factor * residual_last_iterate) {
         ++n_steps_no_improvement;
      }
      else {
         n_steps_no_improvement = 0;
      }

      if (residual_norm_current_solution < tolerance_iterative_refinement * rhs_norm) {
         convergence_status = IterativeSolverSolutionStatus::CONVERGED;
         break;
      }
      else if (current_iterative_refinement_step == max_iterative_refinement_steps) {
         convergence_status = IterativeSolverSolutionStatus::NOT_CONVERGED_MAX_ITERATIONS;
         break;
      }
      else if (max_divergence_factor * residual_norm_best_solution < residual_norm_current_solution) {
         convergence_status = IterativeSolverSolutionStatus::DIVERGED;
         break;
      }
      else if (n_steps_no_improvement > max_steps_stalling) {
         convergence_status = IterativeSolverSolutionStatus::STAGNATION;
         break;
      }

      /* store best solution if we compute one more in the next round */
      if (residual_norm_current_solution < residual_norm_best_solution) {
         residual_norm_best_solution = residual_norm_current_solution;
         best_solution.copyFrom(solution);
      }
      residual_last_iterate = residual_norm_current_solution;

   } while (true);

   assert(convergence_status != IterativeSolverSolutionStatus::DID_NOT_RUN);

   /* backtrack to best solution if it was not the last one */
   if (residual_norm_best_solution < residual_norm_current_solution)
      solution.copyFrom(best_solution);

   if (PIPS_MPIgetRank() == 0) {
      const double norm_final = std::min(residual_norm_best_solution, residual_norm_current_solution);
      std::cout << "IterRef (it=" << current_iterative_refinement_step << ", residual norm=" << norm_final << ", relative residual norm="
                << norm_final / rhs_norm << ") " << convergence_status << "\n";
   }
}

/**
 * res = res - mat * sol
 *
 * [ resx ]   [ resx ]   [ Q + dd + regularization           AT                    CT               ] [ solx ]
 * [ resy ] = [ resy ] - [       A                     regularization              0                ] [ soly ]
 * [ resz ]   [ resz ]   [       C                            0          nOmegaInv + regularization ] [ solz ]
 *
 * stepx, stepy, stepz are used as temporary buffers
 */
void LinearSystem::compute_regularized_system_residuals(const Vector<double>& sol, Vector<double>& res, Vector<double>& solx, Vector<double>& soly,
      Vector<double>& solz, const Problem& problem) {
   const double alpha = -1.0;
   const double beta = 1.0;

   system_mult(beta, res, alpha, sol, problem, solx, soly, solz, true);
}

void LinearSystem::compute_system_residuals(const Vector<double>& sol, Vector<double>& res, Vector<double>& solx, Vector<double>& soly,
      Vector<double>& solz, const Problem& problem) {
   const double alpha = -1.0;
   const double beta = 1.0;

   system_mult(beta, res, alpha, sol, problem, solx, soly, solz, false);
}

void
LinearSystem::joinRHS(Vector<double>& rhs_in, const Vector<double>& rhs1_in, const Vector<double>& rhs2_in, const Vector<double>& rhs3_in) const {
   rhs_in.jointCopyFrom(rhs1_in, rhs2_in, rhs3_in);
}

void LinearSystem::separateVars(Vector<double>& x_in, Vector<double>& y_in, Vector<double>& z_in, const Vector<double>& vars_in) const {
   vars_in.jointCopyTo(x_in, y_in, z_in);
}
