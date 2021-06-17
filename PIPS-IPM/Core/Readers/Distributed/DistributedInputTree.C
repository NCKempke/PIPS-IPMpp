#include "DistributedInputTree.h"
#include <cstdlib>

//***********************************************************
//************************** TREE ***************************
//***********************************************************
DistributedInputTree::DistributedInputTree(std::unique_ptr<DistributedInputNode> root) : nodeInput{std::move(root)} {
}

void DistributedInputTree::add_child(std::unique_ptr<DistributedInputNode> node) {
   children.push_back(std::make_unique<DistributedInputTree>(std::move(node)));
}

void DistributedInputTree::add_child(std::unique_ptr<DistributedInputTree> subTree) {
   children.push_back(std::move(subTree));
}

//***********************************************************
//************************** NODE ***************************
//***********************************************************
DistributedInputTree::DistributedInputNode::DistributedInputNode(int id_, int n_, int my_, int myl_, int mz_, int mzl_) : id{id_}, n{n_}, my{my_}, myl{myl_}, mz{mz_},
      mzl{mzl_} {}

DistributedInputTree::DistributedInputNode::DistributedInputNode(void* user_data_, int id_, int n_, int my_, int mz_, FMAT fQ_, FNNZ fnnzQ_, FVEC fc_, FMAT fA_,
      FNNZ fnnzA_, FMAT fB_, FNNZ fnnzB_, FVEC fb_, FMAT fC_, FNNZ fnnzC_, FMAT fD_, FNNZ fnnzD_, FVEC fclow_, FVEC ficlow_, FVEC fcupp_,
      FVEC ficupp_, FVEC fxlow_, FVEC fixlow_, FVEC fxupp_, FVEC fixupp_, bool deleteUserData_/*=false*/) : id(id_), n(n_), my(my_), mz(mz_),
      fnnzQ(fnnzQ_), fnnzA(fnnzA_), fnnzB(fnnzB_), fnnzC(fnnzC_), fnnzD(fnnzD_), fQ(fQ_), fA(fA_), fB(fB_), fC(fC_), fD(fD_), fc(fc_), fb(fb_),
      fclow(fclow_), fcupp(fcupp_), ficlow(ficlow_), ficupp(ficupp_), fxlow(fxlow_), fxupp(fxupp_), fixlow(fixlow_), fixupp(fixupp_),
      user_data(user_data_), deleteUserData(deleteUserData_) {}


// full callback constructor without constraints
DistributedInputTree::DistributedInputNode::DistributedInputNode(void* user_data_, int id_, FNNZ n_, FNNZ my_, FNNZ mz_, FMAT fQ_, FNNZ fnnzQ_, FVEC fc_, FMAT fA_,
      FNNZ fnnzA_, FMAT fB_, FNNZ fnnzB_, FVEC fb_, FMAT fC_, FNNZ fnnzC_, FMAT fD_, FNNZ fnnzD_, FVEC fclow_, FVEC ficlow_, FVEC fcupp_,
      FVEC ficupp_, FVEC fxlow_, FVEC fixlow_, FVEC fxupp_, FVEC fixupp_, bool deleteUserData_/*=false*/) : id(id_), nCall(n_), myCall(my_),
      mzCall(mz_), fnnzQ(fnnzQ_), fnnzA(fnnzA_), fnnzB(fnnzB_), fnnzC(fnnzC_), fnnzD(fnnzD_), fQ(fQ_), fA(fA_), fB(fB_), fC(fC_), fD(fD_), fc(fc_),
      fb(fb_), fclow(fclow_), fcupp(fcupp_), ficlow(ficlow_), ficupp(ficupp_), fxlow(fxlow_), fxupp(fxupp_), fixlow(fixlow_), fixupp(fixupp_),
      user_data(user_data_), deleteUserData(deleteUserData_) {}

// full callback constructor including linking constraints
DistributedInputTree::DistributedInputNode::DistributedInputNode(void* user_data_, int id_, FNNZ n_, FNNZ my_, FNNZ myl_, FNNZ mz_, FNNZ mzl_, FMAT fQ_, FNNZ fnnzQ_,
      FVEC fc_, FMAT fA_, FNNZ fnnzA_, FMAT fB_, FNNZ fnnzB_, FMAT fBl_, FNNZ fnnzBl_, FVEC fb_, FVEC fbl_, FMAT fC_, FNNZ fnnzC_, FMAT fD_,
      FNNZ fnnzD_, FMAT fDl_, FNNZ fnnzDl_, FVEC fclow_, FVEC ficlow_, FVEC fcupp_, FVEC ficupp_, FVEC fdllow_, FVEC fidllow_, FVEC fdlupp_,
      FVEC fidlupp_, FVEC fxlow_, FVEC fixlow_, FVEC fxupp_, FVEC fixupp_, bool deleteUserData_/*=false*/) : id(id_), nCall(n_), myCall(my_),
      mzCall(mz_), mylCall(myl_), mzlCall(mzl_), fnnzQ(fnnzQ_), fnnzA(fnnzA_), fnnzB(fnnzB_), fnnzBl(fnnzBl_), fnnzC(fnnzC_), fnnzD(fnnzD_),
      fnnzDl(fnnzDl_), fQ(fQ_), fA(fA_), fB(fB_), fBl(fBl_), fC(fC_), fD(fD_), fDl(fDl_), fc(fc_), fb(fb_), fbl(fbl_), fclow(fclow_), fcupp(fcupp_),
      ficlow(ficlow_), ficupp(ficupp_), fdllow(fdllow_), fdlupp(fdlupp_), fidllow(fidllow_), fidlupp(fidlupp_), fxlow(fxlow_), fxupp(fxupp_),
      fixlow(fixlow_), fixupp(fixupp_), user_data(user_data_), deleteUserData(deleteUserData_) {}

DistributedInputTree::DistributedInputNode::DistributedInputNode(int id_) : id(id_) {}


DistributedInputTree::DistributedInputNode::~DistributedInputNode() {
   if (deleteUserData)
      free(user_data);
}
