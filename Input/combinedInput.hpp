#ifndef COMBINEDINPUTHPP
#define COMBINEDINPUTHPP

#include "stochasticInput.hpp"

// wrapper for combining scenarios for lagrangian subproblems

class combinedInput : public stochasticInput {
public:
	combinedInput(stochasticInput &inner, std::vector<std::vector<int> > const& scenarioMap);
	int nScenarios() override { return scenarioMap.size(); }
	int nFirstStageVars() override { return inner.nFirstStageVars(); }
	int nFirstStageCons() override { return inner.nFirstStageCons(); }
	int nSecondStageVars(int scen) override;
	int nSecondStageCons(int scen) override;

	std::vector<double> getFirstStageColLB() override { return inner.getFirstStageColLB(); }
	std::vector<double> getFirstStageColUB() override { return inner.getFirstStageColUB(); }
	std::vector<double> getFirstStageObj() override { return inner.getFirstStageObj(); }
	std::vector<std::string> getFirstStageColNames() override { return inner.getFirstStageColNames(); }
	std::vector<double> getFirstStageRowLB() override { return inner.getFirstStageRowLB(); }
	std::vector<double> getFirstStageRowUB() override { return inner.getFirstStageRowUB(); }
	std::vector<std::string> getFirstStageRowNames() override { return inner.getFirstStageRowNames(); }
	bool isFirstStageColInteger(int col) override { return inner.isFirstStageColInteger(col); }

	std::vector<double> getSecondStageColLB(int scen) override;
	std::vector<double> getSecondStageColUB(int scen) override;
	std::vector<double> getSecondStageObj(int scen) override;
	std::vector<std::string> getSecondStageColNames(int scen) override;
	std::vector<double> getSecondStageRowUB(int scen) override;
	std::vector<double> getSecondStageRowLB(int scen) override;
	std::vector<std::string> getSecondStageRowNames(int scen) override;
	double scenarioProbability(int scen) override;
	bool isSecondStageColInteger(int scen, int col) override;

	CoinPackedMatrix getFirstStageConstraints() override { return inner.getFirstStageConstraints(); }
	CoinPackedMatrix getSecondStageConstraints(int scen) override;
	CoinPackedMatrix getLinkingConstraints(int scen) override;

	

	
	bool scenarioDimensionsEqual() override { return inner.scenarioDimensionsEqual() && equalScenarios; }
	bool onlyBoundsVary() override { return inner.onlyBoundsVary(); }
	bool allProbabilitiesEqual() override { return equalScenarios && inner.allProbabilitiesEqual(); }
	bool continuousRecourse() override { return inner.continuousRecourse(); }

private:
	// map from "fake" scenario index to group of scenarios it represents
	std::vector<std::vector<int> > scenarioMap;
	bool equalScenarios; // equal number of scenarios in each combined scenario
	stochasticInput &inner;

};



#endif
