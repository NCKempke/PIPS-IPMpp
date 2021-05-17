#include <iostream>
#include <SimpleVector.h>
#include "Variables.h"
#include "Vector.hpp"
#include "SmartPointer.h"

#include "Problem.h"
#include "MpsReader.h"

Variables::Variables(Vector<double>* x_in, Vector<double>* s_in, Vector<double>* y_in, Vector<double>* z_in, Vector<double>* v_in,
      Vector<double>* gamma_in, Vector<double>* w_in, Vector<double>* phi_in, Vector<double>* t_in, Vector<double>* lambda_in, Vector<double>* u_in,
      Vector<double>* pi_in, Vector<double>* ixlow_in, Vector<double>* ixupp_in, Vector<double>* iclow_in, Vector<double>* icupp_in) {
   SpReferTo(x, x_in);
   SpReferTo(s, s_in);
   SpReferTo(y, y_in);
   SpReferTo(z, z_in);
   SpReferTo(v, v_in);
   SpReferTo(phi, phi_in);
   SpReferTo(w, w_in);
   SpReferTo(gamma, gamma_in);
   SpReferTo(t, t_in);
   SpReferTo(lambda, lambda_in);
   SpReferTo(u, u_in);
   SpReferTo(pi, pi_in);
   SpReferTo(ixlow, ixlow_in);
   SpReferTo(ixupp, ixupp_in);
   SpReferTo(iclow, iclow_in);
   SpReferTo(icupp, icupp_in);

   nx = x->length();
   my = y->length();
   mz = z->length();

   assert(nx == ixlow->length() || 0 == ixlow->length());
   assert(nx == ixlow->length() || 0 == ixlow->length());
   assert(mz == iclow->length() || 0 == iclow->length());
   assert(mz == icupp->length() || 0 == icupp->length());

   nxlow = ixlow->numberOfNonzeros();
   nxupp = ixupp->numberOfNonzeros();
   mclow = iclow->numberOfNonzeros();
   mcupp = icupp->numberOfNonzeros();
   nComplementaryVariables = mclow + mcupp + nxlow + nxupp;

   assert(mz == s->length());
   assert(nx == v->length() || (0 == v->length() && nxlow == 0));
   assert(nx == gamma->length() || (0 == gamma->length() && nxlow == 0));

   assert(nx == w->length() || (0 == w->length() && nxupp == 0));
   assert(nx == phi->length() || (0 == phi->length() && nxupp == 0));

   assert(mz == t->length() || (0 == t->length() && mclow == 0));
   assert(mz == lambda->length() || (0 == lambda->length() && mclow == 0));

   assert(mz == u->length() || (0 == u->length() && mcupp == 0));
   assert(mz == pi->length() || (0 == pi->length() && mcupp == 0));
}

Variables::Variables(const Variables& vars) {
   ixlow = SmartPointer<Vector<double> >(vars.ixlow->cloneFull());
   ixupp = SmartPointer<Vector<double> >(vars.ixupp->cloneFull());
   iclow = SmartPointer<Vector<double> >(vars.iclow->cloneFull());
   icupp = SmartPointer<Vector<double> >(vars.icupp->cloneFull());

   nx = vars.nx;
   my = vars.my;
   mz = vars.mz;

   nxlow = ixlow->numberOfNonzeros();
   nxupp = ixupp->numberOfNonzeros();
   mclow = iclow->numberOfNonzeros();
   mcupp = icupp->numberOfNonzeros();

   s = SmartPointer<Vector<double> >(vars.s->cloneFull());

   t = SmartPointer<Vector<double> >(vars.t->cloneFull());
   lambda = SmartPointer<Vector<double> >(vars.lambda->cloneFull());

   u = SmartPointer<Vector<double> >(vars.u->cloneFull());
   pi = SmartPointer<Vector<double> >(vars.pi->cloneFull());

   v = SmartPointer<Vector<double> >(vars.v->cloneFull());
   gamma = SmartPointer<Vector<double> >(vars.gamma->cloneFull());

   w = SmartPointer<Vector<double> >(vars.w->cloneFull());
   phi = SmartPointer<Vector<double> >(vars.phi->cloneFull());

   x = SmartPointer<Vector<double> >(vars.x->cloneFull());
   y = SmartPointer<Vector<double> >(vars.y->cloneFull());
   z = SmartPointer<Vector<double> >(vars.z->cloneFull());
   nComplementaryVariables = mclow + mcupp + nxlow + nxupp;
}

double Variables::getAverageDistanceToBoundForConvergedVars(const Problem&, double tol) const {
   assert(0 < tol);

   double sum_small_distance = 0.0;
   int n_close = 0;
   v->getSumCountIfSmall(tol, sum_small_distance, n_close, &*ixlow);
   w->getSumCountIfSmall(tol, sum_small_distance, n_close, &*ixupp);
   u->getSumCountIfSmall(tol, sum_small_distance, n_close, &*icupp);
   t->getSumCountIfSmall(tol, sum_small_distance, n_close, &*iclow);

   if (n_close == 0)
      return std::numeric_limits<double>::infinity();
   else
      return sum_small_distance / (double) n_close;
}


void Variables::pushSlacksFromBound(double tol, double amount) {
   if (nxlow > 0)
      v->pushAwayFromZero(tol, amount, &*ixlow);
   if (nxupp > 0)
      w->pushAwayFromZero(tol, amount, &*ixupp);
   if (mclow > 0)
      t->pushAwayFromZero(tol, amount, &*iclow);
   if (mcupp > 0)
      u->pushAwayFromZero(tol, amount, &*icupp);
}


double Variables::mu() {
   double mu = 0.;
   if (nComplementaryVariables == 0) {
      return 0.;
   }
   else {

      if (mclow > 0)
         mu += t->dotProductWith(*lambda);
      if (mcupp > 0)
         mu += u->dotProductWith(*pi);
      if (nxlow > 0)
         mu += v->dotProductWith(*gamma);
      if (nxupp > 0)
         mu += w->dotProductWith(*phi);

      mu /= nComplementaryVariables;
      return mu;
   }
}

double Variables::mustep_pd(const Variables* iterate, double alpha_primal, double alpha_dual) {
   double mu = 0.;
   if (nComplementaryVariables == 0) {
      return 0.;
   }
   else {
      if (mclow > 0) {
         mu += t->shiftedDotProductWith(alpha_primal, *iterate->t, *lambda, alpha_dual, *iterate->lambda);
      }
      if (mcupp > 0) {
         mu += u->shiftedDotProductWith(alpha_primal, *iterate->u, *pi, alpha_dual, *iterate->pi);
      }
      if (nxlow > 0) {
         mu += v->shiftedDotProductWith(alpha_primal, *iterate->v, *gamma, alpha_dual, *iterate->gamma);
      }
      if (nxupp > 0) {
         mu += w->shiftedDotProductWith(alpha_primal, *iterate->w, *phi, alpha_dual, *iterate->phi);
      }
      mu /= nComplementaryVariables;
      return mu;
   }
}

void Variables::saxpy(const Variables* iterate, double alpha) {
   x->axpy(alpha, *iterate->x);
   y->axpy(alpha, *iterate->y);
   z->axpy(alpha, *iterate->z);
   s->axpy(alpha, *iterate->s);
   if (mclow > 0) {
      assert(iterate->t->matchesNonZeroPattern(*iclow) && iterate->lambda->matchesNonZeroPattern(*iclow));

      t->axpy(alpha, *iterate->t);
      lambda->axpy(alpha, *iterate->lambda);
   }
   if (mcupp > 0) {
      assert(iterate->u->matchesNonZeroPattern(*icupp) && iterate->pi->matchesNonZeroPattern(*icupp));

      u->axpy(alpha, *iterate->u);
      pi->axpy(alpha, *iterate->pi);
   }
   if (nxlow > 0) {
      assert(iterate->v->matchesNonZeroPattern(*ixlow) && iterate->gamma->matchesNonZeroPattern(*ixlow));

      v->axpy(alpha, *iterate->v);
      gamma->axpy(alpha, *iterate->gamma);
   }
   if (nxupp > 0) {
      assert(iterate->w->matchesNonZeroPattern(*ixupp) && iterate->phi->matchesNonZeroPattern(*ixupp));

      w->axpy(alpha, *iterate->w);
      phi->axpy(alpha, *iterate->phi);
   }
}

void Variables::saxpy_pd(const Variables* iterate, double alpha_primal, double alpha_dual) {
   x->axpy(alpha_primal, *iterate->x);
   y->axpy(alpha_dual, *iterate->y);
   z->axpy(alpha_dual, *iterate->z);
   s->axpy(alpha_primal, *iterate->s);
   if (mclow > 0) {
      assert(iterate->t->matchesNonZeroPattern(*iclow) && iterate->lambda->matchesNonZeroPattern(*iclow));

      t->axpy(alpha_primal, *iterate->t);
      lambda->axpy(alpha_dual, *iterate->lambda);
   }
   if (mcupp > 0) {
      assert(iterate->u->matchesNonZeroPattern(*icupp) && iterate->pi->matchesNonZeroPattern(*icupp));

      u->axpy(alpha_primal, *iterate->u);
      pi->axpy(alpha_dual, *iterate->pi);
   }
   if (nxlow > 0) {
      assert(iterate->v->matchesNonZeroPattern(*ixlow) && iterate->gamma->matchesNonZeroPattern(*ixlow));

      v->axpy(alpha_primal, *iterate->v);
      gamma->axpy(alpha_dual, *iterate->gamma);
   }
   if (nxupp > 0) {
      assert(iterate->w->matchesNonZeroPattern(*ixupp) && iterate->phi->matchesNonZeroPattern(*ixupp));

      w->axpy(alpha_primal, *iterate->w);
      phi->axpy(alpha_dual, *iterate->phi);
   }
}

void Variables::negate() {
   s->negate();
   x->negate();
   y->negate();
   z->negate();
   if (mclow > 0) {
      t->negate();
      lambda->negate();
   }
   if (mcupp > 0) {
      u->negate();
      pi->negate();
   }
   if (nxlow > 0) {
      v->negate();
      gamma->negate();
   }
   if (nxupp > 0) {
      w->negate();
      phi->negate();
   }
}

double Variables::stepbound(const Variables* iterate) {
   double maxStep = 1.0;
   if (mclow > 0) {
      assert(t->somePositive(*iclow));
      assert(lambda->somePositive(*iclow));

      maxStep = t->stepbound(*iterate->t, maxStep);
      maxStep = lambda->stepbound(*iterate->lambda, maxStep);
   }

   if (mcupp > 0) {
      assert(u->somePositive(*icupp));
      assert(pi->somePositive(*icupp));

      maxStep = u->stepbound(*iterate->u, maxStep);
      maxStep = pi->stepbound(*iterate->pi, maxStep);
   }

   if (nxlow > 0) {
      assert(v->somePositive(*ixlow));
      assert(gamma->somePositive(*ixlow));

      maxStep = v->stepbound(*iterate->v, maxStep);
      maxStep = gamma->stepbound(*iterate->gamma, maxStep);
   }

   if (nxupp > 0) {
      assert(w->somePositive(*ixupp));
      assert(phi->somePositive(*ixupp));

      maxStep = w->stepbound(*iterate->w, maxStep);
      maxStep = phi->stepbound(*iterate->phi, maxStep);
   }

   assert(maxStep <= 1.0);
   return maxStep;
}

void Variables::stepbound_pd(const Variables* iterate, double& alpha_primal, double& alpha_dual) {
   double maxStep_primal = 1.0;
   double maxStep_dual = 1.0;

   if (mclow > 0) {
      assert(t->somePositive(*iclow));
      assert(lambda->somePositive(*iclow));

      maxStep_primal = t->stepbound(*iterate->t, maxStep_primal);
      maxStep_dual = lambda->stepbound(*iterate->lambda, maxStep_dual);
   }

   if (mcupp > 0) {
      assert(u->somePositive(*icupp));
      assert(pi->somePositive(*icupp));

      maxStep_primal = u->stepbound(*iterate->u, maxStep_primal);
      maxStep_dual = pi->stepbound(*iterate->pi, maxStep_dual);
   }

   if (nxlow > 0) {
      assert(v->somePositive(*ixlow));
      assert(gamma->somePositive(*ixlow));

      maxStep_primal = v->stepbound(*iterate->v, maxStep_primal);
      maxStep_dual = gamma->stepbound(*iterate->gamma, maxStep_dual);
   }

   if (nxupp > 0) {
      assert(w->somePositive(*ixupp));
      assert(phi->somePositive(*ixupp));

      maxStep_primal = w->stepbound(*iterate->w, maxStep_primal);
      maxStep_dual = phi->stepbound(*iterate->phi, maxStep_dual);
   }

   assert(maxStep_primal <= 1.0);
   assert(maxStep_dual <= 1.0);

   alpha_primal = maxStep_primal;
   alpha_dual = maxStep_dual;
}

double
Variables::findBlocking(const Variables* step, double& primalValue, double& primalStep, double& dualValue, double& dualStep, int& firstOrSecond) {
   double alpha = 1.;
   firstOrSecond = 0;

   if (mclow > 0) {
      alpha = t->find_blocking(*step->t, *lambda, *step->lambda, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   if (mcupp > 0) {
      alpha = u->find_blocking(*step->u, *pi, *step->pi, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   if (nxlow > 0) {
      alpha = v->find_blocking(*step->v, *gamma, *step->gamma, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   if (nxupp > 0) {
      alpha = w->find_blocking(*step->w, *phi, *step->phi, alpha, &primalValue, &primalStep, &dualValue, &dualStep, firstOrSecond);
   }
   return alpha;
}

void
Variables::findBlocking_pd(const Variables* step, double& primalValue, double& primalStep, double& dualValue, double& dualStep, double& primalValue_d,
      double& primalStep_d, double& dualValue_d, double& dualStep_d, double& alphaPrimal, double& alphaDual, bool& primalBlocking,
      bool& dualBlocking) {
   alphaPrimal = 1.0, alphaDual = 1.0;
   primalBlocking = false, dualBlocking = false;

   if (mclow > 0) {
      t->find_blocking_pd(*step->t, *lambda, *step->lambda, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d,
            primalStep_d, dualValue_d, dualStep_d, primalBlocking, dualBlocking);
   }

   if (mcupp > 0) {
      u->find_blocking_pd(*step->u, *pi, *step->pi, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d, primalStep_d,
            dualValue_d, dualStep_d, primalBlocking, dualBlocking);
   }

   if (nxlow > 0) {
      v->find_blocking_pd(*step->v, *gamma, *step->gamma, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d,
            primalStep_d, dualValue_d, dualStep_d, primalBlocking, dualBlocking);
   }

   if (nxupp > 0) {
      w->find_blocking_pd(*step->w, *phi, *step->phi, alphaPrimal, alphaDual, primalValue, primalStep, dualValue, dualStep, primalValue_d,
            primalStep_d, dualValue_d, dualStep_d, primalBlocking, dualBlocking);
   }
}

void Variables::push_to_interior(double alpha, double beta) {
   s->setToZero();
   x->setToZero();
   y->setToZero();
   z->setToZero();

   if (nxlow > 0) {
      v->setToConstant(alpha);
      v->selectNonZeros(*ixlow);
      gamma->setToConstant(beta);
      gamma->selectNonZeros(*ixlow);
   }
   if (nxupp > 0) {
      w->setToConstant(alpha);
      w->selectNonZeros(*ixupp);
      phi->setToConstant(beta);
      phi->selectNonZeros(*ixupp);
   }

   if (mclow > 0) {
      t->setToConstant(alpha);
      t->selectNonZeros(*iclow);
      lambda->setToConstant(beta);
      lambda->selectNonZeros(*iclow);
   }
   if (mcupp > 0) {
      u->setToConstant(alpha);
      u->selectNonZeros(*icupp);
      pi->setToConstant(beta);
      pi->selectNonZeros(*icupp);
   }
}

double Variables::violation() {
   double viol = 0.0, cmin = 0.0;
   int iblock;

   if (nxlow > 0) {
      v->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;

      gamma->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;
   }
   if (nxupp > 0) {
      w->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;

      phi->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;
   }
   if (mclow > 0) {
      t->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;

      lambda->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;
   }
   if (mcupp > 0) {
      u->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;

      pi->min(cmin, iblock);
      if (cmin < viol)
         viol = cmin;
   }
   return -viol;
}

void Variables::shiftBoundVariables(double alpha, double beta) {
   if (nxlow > 0) {
      v->add_constant(alpha, *ixlow);
      gamma->add_constant(beta, *ixlow);
   }
   if (nxupp > 0) {
      w->add_constant(alpha, *ixupp);
      phi->add_constant(beta, *ixupp);
   }
   if (mclow > 0) {
      t->add_constant(alpha, *iclow);
      lambda->add_constant(beta, *iclow);
   }
   if (mcupp > 0) {
      u->add_constant(alpha, *icupp);
      pi->add_constant(beta, *icupp);
   }
}

void Variables::copy(const Variables* b_in) {
   const Variables* b = (const Variables*) b_in;

   s->copyFrom(*b->s);
   if (nxlow > 0) {
      v->copyFrom(*b->v);
      gamma->copyFrom(*b->gamma);
   }
   if (nxupp > 0) {
      w->copyFrom(*b->w);
      phi->copyFrom(*b->phi);
   }
   if (mclow > 0) {
      t->copyFrom(*b->t);
      lambda->copyFrom(*b->lambda);
   }
   if (mcupp > 0) {
      u->copyFrom(*b->u);
      pi->copyFrom(*b->pi);
   }
   x->copyFrom(*b->x);
   y->copyFrom(*b->y);
   z->copyFrom(*b->z);

}

double Variables::onenorm() const {
   double norm = x->onenorm();
   norm += s->onenorm();
   norm += y->onenorm();
   norm += z->onenorm();

   norm += v->onenorm();
   norm += phi->onenorm();
   norm += w->onenorm();
   norm += gamma->onenorm();
   norm += t->onenorm();
   norm += lambda->onenorm();
   norm += u->onenorm();
   norm += pi->onenorm();

   return norm;
}

double Variables::infnorm() const {
   double norm = 0.0;
   double temp = x->infnorm();
   if (temp > norm)
      norm = temp;
   temp = s->infnorm();
   if (temp > norm)
      norm = temp;
   temp = y->infnorm();
   if (temp > norm)
      norm = temp;
   temp = z->infnorm();
   if (temp > norm)
      norm = temp;

   temp = v->infnorm();
   if (temp > norm)
      norm = temp;
   temp = phi->infnorm();
   if (temp > norm)
      norm = temp;

   temp = w->infnorm();
   if (temp > norm)
      norm = temp;
   temp = gamma->infnorm();
   if (temp > norm)
      norm = temp;

   temp = t->infnorm();
   if (temp > norm)
      norm = temp;
   temp = lambda->infnorm();
   if (temp > norm)
      norm = temp;

   temp = u->infnorm();
   if (temp > norm)
      norm = temp;
   temp = pi->infnorm();
   if (temp > norm)
      norm = temp;

   return norm;
}

void Variables::setToZero() {
   x->setToZero();
   s->setToZero();
   y->setToZero();
   z->setToZero();

   v->setToZero();
   gamma->setToZero();

   w->setToZero();
   phi->setToZero();

   t->setToZero();
   lambda->setToZero();

   u->setToZero();
   pi->setToZero();
}

int Variables::validNonZeroPattern() {
   if (nxlow > 0 && (!v->matchesNonZeroPattern(*ixlow) || !gamma->matchesNonZeroPattern(*ixlow))) {

      if (!v->matchesNonZeroPattern(*ixlow))
         printf("invalidNonZeroPattern v\n");
      if (!gamma->matchesNonZeroPattern(*ixlow))
         printf("invalidNonZeroPattern gamma\n");
      return 0;
   }

   if (nxupp > 0 && (!w->matchesNonZeroPattern(*ixupp) || !phi->matchesNonZeroPattern(*ixupp))) {
      if (!w->matchesNonZeroPattern(*ixupp))
         printf("invalidNonZeroPattern w\n");
      if (!phi->matchesNonZeroPattern(*ixupp))
         printf("invalidNonZeroPattern phi\n");
      return 0;
   }
   if (mclow > 0 && (!t->matchesNonZeroPattern(*iclow) || !lambda->matchesNonZeroPattern(*iclow))) {
      if (!t->matchesNonZeroPattern(*iclow))
         printf("invalidNonZeroPattern t\n");
      if (!lambda->matchesNonZeroPattern(*iclow))
         printf("invalidNonZeroPattern lambda\n");
      return 0;
   }

   if (mcupp > 0 && (!u->matchesNonZeroPattern(*icupp) || !pi->matchesNonZeroPattern(*icupp))) {
      if (!u->matchesNonZeroPattern(*icupp))
         printf("invalidNonZeroPattern u\n");
      if (!phi->matchesNonZeroPattern(*icupp))
         printf("invalidNonZeroPattern phi\n");
      return 0;
   }

   return 1;
}

void Variables::unscaleSolution(Problem* problem) {
// Modifying sx is equivalent to modifying x
   SimpleVector<double>& sx = (SimpleVector<double>&) *this->x;

// x = D * x'
   sx.componentMult(problem->scale());
}

void Variables::unscaleBounds(Problem* problem) {
   SimpleVector<double>& sxlow = (SimpleVector<double>&) problem->xlowerBound();
   SimpleVector<double>& sxupp = (SimpleVector<double>&) problem->xupperBound();

// l = D * l'
   sxlow.componentMult(problem->scale());

// u = D * u'
   sxupp.componentMult(problem->scale());
}

void Variables::printSolution(MpsReader* reader, Problem* problem, int& iErr) {
   assert(x->isKindOf(kSimpleVector)); // Otherwise this routine

   SimpleVector<double> g(nx);
   problem->getg(g);
   problem->hessian_multiplication(1.0, g, 0.5, *x);
   double objective = g.dotProductWith(*x);

   SimpleVector<double>& sx = (SimpleVector<double>&) *this->x;
   SimpleVector<double>& sxlow = (SimpleVector<double>&) problem->xlowerBound();
   SimpleVector<double>& sixlow = (SimpleVector<double>&) problem->ixlowerBound();
   SimpleVector<double>& sxupp = (SimpleVector<double>&) problem->xupperBound();
   SimpleVector<double>& sixupp = (SimpleVector<double>&) problem->ixupperBound();
   SimpleVector<double>& sgamma = (SimpleVector<double>&) *this->gamma;
   SimpleVector<double>& sphi = (SimpleVector<double>&) *this->phi;
   SimpleVector<double>& sy = (SimpleVector<double>&) *this->y;
   SimpleVector<double>& ss = (SimpleVector<double>&) *this->s;
   SimpleVector<double>& slambda = (SimpleVector<double>&) *this->lambda;
   SimpleVector<double>& spi = (SimpleVector<double>&) *this->pi;
   SimpleVector<double>& sz = (SimpleVector<double>&) *this->z;
   SimpleVector<double>& sclow = (SimpleVector<double>&) problem->slowerBound();
   SimpleVector<double>& siclow = (SimpleVector<double>&) problem->islowerBound();
   SimpleVector<double>& scupp = (SimpleVector<double>&) problem->supperBound();
   SimpleVector<double>& sicupp = (SimpleVector<double>&) problem->isupperBound();

   char* cxupp = new char[nx];
   char* cxlow = new char[nx];
   for (int j = 0; j < nx; j++) {
      if (nxupp > 0 && sixupp[j] != 0) {
         cxupp[j] = 1;
      }
      else {
         cxupp[j] = 0;
      }
      if (nxlow > 0 && sixlow[j] != 0) {
         cxlow[j] = 1;
      }
      else {
         cxlow[j] = 0;
      }
   }
   char* cclow, * ccupp;
   if (mz <= 0) {
      cclow = 0;
      ccupp = 0;
   }
   else {
      cclow = new char[mz];
      ccupp = new char[mz];
      for (int i = 0; i < mz; i++) {
         if (mclow > 0 && siclow[i] != 0.0) {
            cclow[i] = 1;
         }
         else {
            cclow[i] = 0;
         }
         if (mcupp > 0 && sicupp[i] != 0.0) {
            ccupp[i] = 1;
         }
         else {
            ccupp[i] = 0;
         }
      }
   }

   if (reader->scalingOption == 1) {
      // Unscale the solution and bounds before printing
      this->unscaleSolution(problem);
      this->unscaleBounds(problem);
   }

   reader->printSolution(sx.elements(), nx, sxlow.elements(), cxlow, sxupp.elements(), cxupp, sgamma.elements(), sphi.elements(), sy.elements(), my,
         ss.elements(), mz, sclow.elements(), cclow, scupp.elements(), ccupp, slambda.elements(), spi.elements(), sz.elements(), objective, iErr);
   delete[] cclow;
   delete[] ccupp;
   delete[] cxlow;
   delete[] cxupp;
}

void Variables::printNorms() const {
   const int my_rank = PIPS_MPIgetRank();

   const double infnorm = this->infnorm();

   if (my_rank == 0)
      std::cout << "||vars||_INF = " << infnorm << std::endl;

   double temp_inf = x->infnorm();
   double temp_2 = x->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_x||_INF = " << temp_inf << "\t||vars_x||_2 = " << temp_2 << std::endl;

   temp_inf = s->infnorm();
   temp_2 = s->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_s||_INF = " << temp_inf << "\t||vars_s||_2 = " << temp_2 << std::endl;

   temp_inf = y->infnorm();
   temp_2 = y->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_y||_INF = " << temp_inf << "\t||vars_y||_2 = " << temp_2 << std::endl;

   temp_inf = z->infnorm();
   temp_2 = z->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_z||_INF = " << temp_inf << "\t||vars_z||_2 = " << temp_2 << std::endl;

   temp_inf = v->infnorm();
   temp_2 = v->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_v||_INF = " << temp_inf << "\t||vars_v||_2 = " << temp_2 << std::endl;
   temp_inf = phi->infnorm();
   temp_2 = phi->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_phi||_INF = " << temp_inf << "\t||vars_phi||_2 = " << temp_2 << std::endl;

   temp_inf = w->infnorm();
   temp_2 = w->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_w||_INF = " << temp_inf << "\t||vars_w||_2 = " << temp_2 << std::endl;
   temp_inf = gamma->infnorm();
   temp_2 = gamma->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_gamma||_INF = " << temp_inf << "\t||vars_gamma||_2 = " << temp_2 << std::endl;

   temp_inf = t->infnorm();
   temp_2 = t->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_t||_INF = " << temp_inf << "\t||vars_t||_2 = " << temp_2 << std::endl;
   temp_inf = lambda->infnorm();
   temp_2 = lambda->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_lambda||_INF = " << temp_inf << "\t||vars_lambda||_2 = " << temp_2 << std::endl;

   temp_inf = u->infnorm();
   temp_2 = u->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_u||_INF = " << temp_inf << "\t||vars_u||_2 = " << temp_2 << std::endl;
   temp_inf = pi->infnorm();
   temp_2 = pi->twonorm();
   if (my_rank == 0)
      std::cout << "||vars_pi||_INF = " << temp_inf << "\t||vars_pi||_2 = " << temp_2 << std::endl;
}

// default implementation for Variables::print() prints abusive
// message. Since we don't have any knowledge of how the variables are
// stored at this top level, we can't do much except print their
// dimensions.

void Variables::print() {
   std::cout << " Complementary Variables = " << nComplementaryVariables << std::endl;
   std::cout << "(Cannot tell you more at this level)" << std::endl;
}
