#ifndef RAWINPUT_HPP
#define RAWINPUT_HPP

#include "stochasticInput.hpp"

// save diagnostic state
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wsuggest-override"
#include "mpi.h"
// turn the warnings back on
#pragma GCC diagnostic pop

// reads format written by dumpSmlModel.cpp
class rawInput : public stochasticInput {
public:

  rawInput(const std::string &datarootname, int overrideScenarioNumber = 0, MPI_Comm comm = MPI_COMM_WORLD);
  rawInput(const std::string &datarootname, const std::string& zerodata, int overrideScenarioNumber = 0, MPI_Comm comm = MPI_COMM_SELF);

	int nScenarios() override { return nScenarios_; }
	int nFirstStageVars() override { return nFirstStageVars_; }
	int nFirstStageCons() override  { return nFirstStageCons_; }
	int nSecondStageVars(int scen) override { return nSecondStageVars_; }
	int nSecondStageCons(int scen) override { return nSecondStageCons_; }

	std::vector<double> getFirstStageColLB() override { return firstStageData.collb; }
	std::vector<double> getFirstStageColUB() override { return firstStageData.colub; }
	std::vector<double> getFirstStageObj() override { return firstStageData.obj; }
	std::vector<std::string> getFirstStageColNames() override { return firstStageData.colnames; }
	std::vector<double> getFirstStageRowLB() override { return firstStageData.rowlb; }
	std::vector<double> getFirstStageRowUB() override { return firstStageData.rowub; }
	std::vector<std::string> getFirstStageRowNames() override { return firstStageData.rownames; }
	bool isFirstStageColInteger(int col) override { return false; }

	std::vector<double> getSecondStageColLB(int scen) override;
	std::vector<double> getSecondStageColUB(int scen) override;
	std::vector<double> getSecondStageObj(int scen) override;
	std::vector<std::string> getSecondStageColNames(int scen) override;
	std::vector<double> getSecondStageRowUB(int scen) override;
	std::vector<double> getSecondStageRowLB(int scen) override;
	std::vector<std::string> getSecondStageRowNames(int scen) override;
	double scenarioProbability(int scen) override { return 1.0/nScenarios_; }
	bool isSecondStageColInteger(int scen, int col) override { return false; }

	CoinPackedMatrix getFirstStageConstraints() override { return Amat; }
	CoinPackedMatrix getSecondStageConstraints(int scen) override { return Wmat; }
	CoinPackedMatrix getLinkingConstraints(int scen) override { return Tmat; }

	

	bool scenarioDimensionsEqual() override { return true; }
	bool onlyBoundsVary() override { return true; }
	bool allProbabilitiesEqual() override { return true; }
	bool continuousRecourse() override { return true; }
	

protected:
	struct scenData {
		std::vector<double> collb, colub, rowlb, rowub, obj;
		std::vector<std::string> rownames, colnames;
		bool didLoad;
		scenData() : didLoad(false) {}
		void initialize(int nvar, int ncons);
	};
	void loadLocalScenData(int scen);
	CoinPackedMatrix Amat, Tmat, Wmat;
	std::vector<scenData> localData;
	const std::string datarootname;
	scenData firstStageData;
	int nScenariosTrue;
	int mype_;

	int nScenarios_, nFirstStageVars_, nFirstStageCons_, nSecondStageVars_, nSecondStageCons_;
private:
  void parseZeroData(const std::string &zerodata, int overrideScenarioNumber);
};


#endif
