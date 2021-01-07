#ifndef STOCH_INPUT_TREE
#define STOCH_INPUT_TREE

#include <vector>

/**
 * The following types define callback functions passed by user to pass 
 * data for each node to the solver.
 */
extern "C" 
typedef int (*FNNZ)(void* user_data, int id, int* nnz);

extern "C" 
typedef int (*FMAT)(void* user_data, int id, 
		    int* krowM, int* jcolM, double* M);
extern "C"
typedef int (*FVEC)(void* user_data, int id, double* vec, int len);

extern "C"
typedef int (*FLEN)(void* user_data, int id, int* len);


class StochInputTree
{
      friend class sTreeCallbacks;
   public:

      /**
       * Inner class that contains the node related data.
       */
      class StochInputNode
      {
         friend class StochInputTree; friend class StochTree; friend class sTreeCallbacks;
  public:
    StochInputNode(int id_ = -1);

	StochInputNode(void* user_data, int id, 
	   int n, int my, int mz,
	   FMAT fQ, FNNZ fnnzQ, FVEC fc,  
	   FMAT fA, FNNZ fnnzA, FMAT fB, FNNZ fnnzB, 
	   FVEC fb, 
	   FMAT fC, FNNZ fnnzC, FMAT fD, FNNZ fnnzD,
	   FVEC fclow, FVEC ficlow, FVEC fcupp, FVEC ficupp,
	   FVEC fxlow, FVEC fixlow, FVEC fxupp, FVEC fixupp,
	   bool deleteUserData=false);

	// full callback constructor without constraints
   StochInputNode(void* user_data, int id,
      FNNZ n, FNNZ my, FNNZ mz,
      FMAT fQ, FNNZ fnnzQ, FVEC fc,
      FMAT fA, FNNZ fnnzA, FMAT fB, FNNZ fnnzB,
      FVEC fb,
      FMAT fC, FNNZ fnnzC, FMAT fD, FNNZ fnnzD,
      FVEC fclow, FVEC ficlow, FVEC fcupp, FVEC ficupp,
      FVEC fxlow, FVEC fixlow, FVEC fxupp, FVEC fixupp,
      bool deleteUserData=false);

	// full callback constructor including linking constraints
	StochInputNode(void* user_data, int id,
	   FNNZ n, FNNZ my, FNNZ myl, FNNZ mz, FNNZ mzl,
	   FMAT fQ, FNNZ fnnzQ, FVEC fc,
	   FMAT fA, FNNZ fnnzA, FMAT fB, FNNZ fnnzB, FMAT fBl, FNNZ fnnzBl,
	   FVEC fb, FVEC fbl,
	   FMAT fC, FNNZ fnnzC, FMAT fD, FNNZ fnnzD, FMAT fDl, FNNZ fnnzDl,
	   FVEC fclow, FVEC ficlow, FVEC fcupp, FVEC ficupp,
	   FVEC fdllow, FVEC fidllow, FVEC fdlupp, FVEC fidlupp,
	   FVEC fxlow, FVEC fixlow, FVEC fxupp, FVEC fixupp,
	   bool deleteUserData=false);

  ~StochInputNode();

  protected:

    int id{-1};
    int n{-1};
    int my{-1};
    int myl{-1};
    int mz{-1};
    int mzl{-1};

    int nnzQ{-1};
    int nnzA{-1};
    int nnzB{-1};
    int nnzBl{-1};
    int nnzC{-1};
    int nnzD{-1};
    int nnzDl{-1};

    //callback functions nCall, myCall, mzCall, mylCall, mzlCall can be nullptr if data is provided through int n, my, mz, myl, mzl
    FNNZ nCall{};
    FNNZ myCall{};
    FNNZ mzCall{};
    FNNZ mylCall{};
    FNNZ mzlCall{};

    FNNZ fnnzQ{};
    FNNZ fnnzA{};
    FNNZ fnnzB{};
    FNNZ fnnzBl{};
    FNNZ fnnzC{};
    FNNZ fnnzD{};
    FNNZ fnnzDl{};

    FMAT fQ{};
    FMAT fA{};
    FMAT fB{};
    FMAT fBl{};
    FMAT fC{};
    FMAT fD{};
    FMAT fDl{};

    FVEC fc{};
    FVEC fb{};
    FVEC fbl{};

    FVEC fclow{};
    FVEC fcupp{};
    FVEC ficlow{};
    FVEC ficupp{};

    FVEC fdllow{};
    FVEC fdlupp{};
    FVEC fidllow{};
    FVEC fidlupp{};

    FVEC fxlow{};
    FVEC fxupp{};
    FVEC fixlow{};
    FVEC fixupp{};
    
    void* user_data{};

    bool deleteUserData{false};

  }; // end of inner class StochInputNode

  friend class StochTree;

 public:

  StochInputTree(const StochInputNode& root);
  StochInputTree(StochInputNode* root);
  virtual ~StochInputTree();

  void AddChild(const StochInputNode  &node);
  void AddChild(StochInputTree *subTree);

  std::vector<StochInputTree*> children;
 protected:
  StochInputTree();
  StochInputNode* nodeInput{};

};

#endif
