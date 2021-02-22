#ifndef SMPSINPUT_HPP
#define SMPSINPUT_HPP

#include "stochasticInput.hpp"
#include "CoinMpsIO.hpp"
#include <cmath>

// Only works with SMPS files with SCENARIO format!
class SMPSInput : public stochasticInput {
public:
	virtual ~SMPSInput() {}
	SMPSInput(std::string const& cor, std::string const& tim, std::string const& sto);
	int nScenarios() override { return nscen; }
	int nFirstStageVars() override { return nvar1; }
	int nFirstStageCons() override { return ncons1; }
	int nSecondStageVars(int scen) override { return nvar2; }
	int nSecondStageCons(int scen) override { return ncons2; }

	std::vector<double> getFirstStageColLB() override { return firstStageData.collb; }
	std::vector<double> getFirstStageColUB() override { return firstStageData.colub; }
	std::vector<double> getFirstStageObj() override { return firstStageData.obj; }
	std::vector<std::string> getFirstStageColNames() override { return firstStageData.colname; }
	std::vector<double> getFirstStageRowLB() override { return firstStageData.rowlb; }
	std::vector<double> getFirstStageRowUB() override { return firstStageData.rowub; }
	std::vector<std::string> getFirstStageRowNames() override { return firstStageData.rowname; }
	bool isFirstStageColInteger(int col) override { return firstStageData.isColInteger.at(col); } 
   virtual bool isFirstStageColBinary(int col) {
	  bool isInteger = this->isFirstStageColInteger(col);
	  // CoinMpsIO has no isBinary member function, but some preprocessing features require
	  // knowledge of binary variables, so kludge in an "isBinary" member function by
	  // relying on CoinMpsIO setting lower and upper bounds to zero and one, respectively.
	  // Also note: CoinMpsIO uses a default tolerance of 1.0e-8 on integrality comparisons.
	  const double intTol = 1.0e-8;
	  bool isLBzero = (fabs(this->getFirstStageColLB().at(col)) < intTol);
	  bool isUBone = (fabs(this->getFirstStageColUB().at(col) - 1.0) < intTol);
	  return (isInteger && isLBzero && isUBone);
	}

	std::vector<double> getSecondStageColLB(int scen) override;
	std::vector<double> getSecondStageColUB(int scen) override;
	// objective vector, already multiplied by probability
	std::vector<double> getSecondStageObj(int scen) override;
	std::vector<std::string> getSecondStageColNames(int scen) override;
	std::vector<double> getSecondStageRowUB(int scen) override;
	std::vector<double> getSecondStageRowLB(int scen) override;
	std::vector<std::string> getSecondStageRowNames(int scen) override;
	double scenarioProbability(int scen) override { return (probabilitiesequal) ? 1.0/nscen : probabilities.at(scen); }
	bool isSecondStageColInteger(int scen, int col) override { return secondStageTemplate.isColInteger.at(col); }
   virtual bool isSecondStageColBinary(int scen, int col) {
	  bool isInteger = this->isSecondStageColInteger(scen, col);
	  // CoinMpsIO has no isBinary member function, but some preprocessing features require
	  // knowledge of binary variables, so kludge in an "isBinary" member function by
	  // relying on CoinMpsIO setting lower and upper bounds to zero and one, respectively.
	  // Also note: CoinMpsIO uses a default tolerance of 1.0e-8 on integrality comparisons.
	  const double intTol = 1.0e-8;
	  bool isLBzero = (fabs(this->getSecondStageColLB(scen).at(col)) < intTol);
	  bool isUBone = (fabs(this->getSecondStageColUB(scen).at(col) - 1.0) < intTol);
	  return (isInteger && isLBzero && isUBone);
	}

	// returns the column-oriented first-stage constraint matrix (A matrix)
	CoinPackedMatrix getFirstStageConstraints() override { return firstStageData.mat; }
	// returns the column-oriented second-stage constraint matrix (W matrix)
	CoinPackedMatrix getSecondStageConstraints(int scen) override;
	// returns the column-oriented matrix linking the first-stage to the second (T matrix)
	CoinPackedMatrix getLinkingConstraints(int scen) override;



	bool scenarioDimensionsEqual() override { return true; }
	bool onlyBoundsVary() override { return onlyboundsvary; }
	bool allProbabilitiesEqual() override { return probabilitiesequal; }
	bool continuousRecourse() override { return continuousrecourse; }


private:
	struct problemData {

		problemData() : ncol(-1), nrow(-1) {}
		bool isInitialized() { return (ncol >= 0 && nrow >= 0); }

		int ncol, nrow;
		std::vector<double> collb, colub, rowlb, rowub, obj;
		std::vector<bool> isColInteger;
		CoinPackedMatrix mat;
		std::vector<std::string> colname, rowname;

		void initialize(int ncol, int nrow) {
			this->ncol = ncol;
			this->nrow = nrow;
			colub.resize(ncol);
			collb.resize(ncol);
			colname.resize(ncol);
			obj.resize(ncol);
			rowub.resize(nrow);
			rowlb.resize(nrow);
			rowname.resize(nrow);
			isColInteger.resize(ncol);
		}

	};

	void cacheScenario(int scen);

	int nscen, nvar1, ncons1, nvar2, ncons2;
	int nvar, ncons; // total variables
	std::vector<problemData> scenarioData;
	problemData firstStageData, secondStageTemplate;
	CoinPackedMatrix TmatTemplate;
	std::vector<CoinPackedMatrix> Tmats;
	std::vector<double> probabilities;
	std::string const corfile, timfile, stofile;
	CoinMpsIO reader;
	bool onlyboundsvary;
	bool probabilitiesequal;
	bool continuousrecourse;

	// save locations in sto file for later reading
	std::vector<std::streampos> scenarioStarts;
	std::vector<int> scenarioLens; // number of lines in sto file for each scenario, not including "SC" line



};

#endif
