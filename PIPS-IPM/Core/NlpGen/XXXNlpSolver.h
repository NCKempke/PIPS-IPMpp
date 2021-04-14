

class XXXNlpSolver : public Solver
{
protected: 
  NlpGenFactory * factory;

  /**  storage for step vectors */
  NlpGenVars * step;

public:
  XXXNlpSolver( NlpGenFactory * opt, NlpDataFull * prob );

  ~XXXNlpSolver();

  virtual int solve( Problem *prob, Variables *iterate, Residuals * resids );


  virtual void defaultMonitor( Problem * data, Variables * vars,
			       Residuals * resids,
			       double alpha, double sigma,
			       int i, double mu,
			       int status_code,
			       int level );

};
