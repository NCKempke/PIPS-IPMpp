#ifndef PROXLINFTRUSTMANGAGER_HPP
#define PROXLINFTRUSTMANGAGER_HPP

#include "bundleManager.hpp"

#include "proximalBAQP.hpp"

// like levelManager, but only accept a step if it's in the right direction
template<typename BAQPSolver, typename LagrangeSolver, typename RecourseSolver>
class proxQPManager : public bundleManager<BAQPSolver,LagrangeSolver,RecourseSolver> {
public:
	proxQPManager(stochasticInput &input, BAContext & ctx) : 
		bundleManager<BAQPSolver,LagrangeSolver,RecourseSolver>(input,ctx) {
		int nscen = input.nScenarios();
		trialSolution.resize(nscen,std::vector<double>(input.nFirstStageVars(),0.));
		lastModelObj.resize(nscen);
		proxCenterModelObj.resize(nscen);
		u = 1.;
		t = MPI_Wtime();
		t2 = 0.;
		mL = 0.1;
		mR = 0.1;
		eps_sol = 10.; // absolute scale, need to adjust
	}


protected: 
	virtual void doStep() {
		using namespace std;
		
		int nvar1 = this->input.nFirstStageVars();
		int nscen = this->input.nScenarios();
		vector<int> const &localScen = this->ctx.localScenarios();

		if (this->nIter == 1) {
			for (int scen = 0; scen < nscen; scen++) {
			//for (unsigned r = 1; r < localScen.size(); r++) {
			//	int scen = localScen[r];
				proxCenterModelObj[scen] = this->bundle[scen].at(0).objmax;
			}
		}
		
		if (this->ctx.mype() == 0) printf("Iter %d Current Objective: %f Best Primal: %f, Relerr: %g Elapsed: %f (%f in QP solve)\n",this->nIter-1,this->currentObj,this->bestPrimalObj,fabs(lastModelObjSum-this->currentObj)/fabs(this->currentObj),MPI_Wtime()-t,t2);
		if (this->terminated_) return;	
		

		proximalQPModel lm(nvar1,this->bundle,this->currentSolution,u);
		BAQPSolver solver(lm);
		
		double tstart = MPI_Wtime();
		solver.go();
		t2 += MPI_Wtime() - tstart;
	
		std::vector<double> const& y = solver.getFirstStagePrimalColSolution();
		double lastModelObjSum_this = 0.;
		double v_this = 0.;
		//for (unsigned r = 1; r < localScen.size(); r++) { // for distributed solver
		for (int scen = 0; scen < nscen; scen++) { // for serial solver
			//int scen = localScen[r];
			std::vector<double> const& z = solver.getSecondStagePrimalColSolution(scen);
			for (int k = 0; k < nvar1; k++) {
				double sum = y[k]/sqrt((double)nscen);
				for (unsigned i = 0; i < this->bundle[scen].size(); i++) {
					sum -= this->bundle[scen][i].subgradient[k]*z[i];
				}
				trialSolution[scen][k] = this->currentSolution[scen][k]+sum/u;
			}
			lastModelObj[scen] = -solver.getSecondStageDualRowSolution(scen)[0];
			lastModelObjSum_this += -lastModelObj[scen];
			v_this += lastModelObj[scen] - proxCenterModelObj[scen] - eps_sol/nscen;
		}
		double v;
		v = v_this;
		lastModelObjSum = lastModelObjSum_this;
		/* if distributed solver:
		MPI_Allreduce(&lastModelObjSum_this,&lastModelObjSum,1,MPI_DOUBLE,MPI_SUM,this->ctx.comm());
		MPI_Allreduce(&v_this,&v,1,MPI_DOUBLE,MPI_SUM,this->ctx.comm());*/
	
		eps_sol = -(mR-mL)*v; 
		if (this->ctx.mype() == 0) printf("v^k = %g, eps_sol = %g\n",v, eps_sol);
		if (v > -1e-2) {
			this->terminated_ = true;
			/*std::vector<double> y = solver.getFirstStagePrimalColSolution();
			for (unsigned k = 0; k < y.size(); k++) y[k] /= -sqrt((double)nscen);
			double val = this->testPrimal(y);
			printf("Primal Obj: %.10g\n",val);*/
			//checkLastPrimals();
			return;
		}
		

		double newObj = this->evaluateSolution(trialSolution, eps_sol/nscen);
		
		// update tau according to p.536 of kiwiel's 1995 paper
		double u_int = 2.*u*(1-(this->currentObj-newObj)/v);
		u = min(max(max(u_int,u/10.),1e-4),10.*u);
		if (this->ctx.mype() == 0) printf("updated u to %g\n",u);

		if (this->ctx.mype() == 0) cout << "Trial solution has obj = " << newObj;
		if (newObj - this->currentObj > -mL*v) {
			if (this->ctx.mype() == 0) cout << ", accepting\n";
			swap(this->currentSolution,trialSolution);
			for (int scen = 0; scen < nscen; scen++) {
			//for (unsigned r = 1; r < localScen.size(); r++) {
			//	int scen = localScen[r];
				proxCenterModelObj[scen] = this->bundle[scen][this->bundle[scen].size()-1].objmax;
			}
			this->currentObj = newObj;
			

		} else {
			if (this->ctx.mype() == 0) cout << ", null step\n";
		}
		


	}

private:
	std::vector<double> lastModelObj;
	std::vector<double> proxCenterModelObj;
	double lastModelObjSum;
	double u;
	double t, t2;
	double mL, mR;
	double eps_sol; // allowable imprecision in the solution 

	std::vector<std::vector<double> > trialSolution; 


};


#endif