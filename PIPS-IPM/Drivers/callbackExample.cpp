#include <InteriorPointMethod.h>
#include "DistributedInputTree.h"
#include "PIPSIPMppInterface.hpp"
#include "MehrotraStrategy.hpp"
#include "mpi.h"

#define LINKING_CONS 1

extern "C" typedef int (* FNNZ)(void* user_data, int id, int* nnz);

/* Row-major format */
extern "C" typedef int (* FMAT)(void* user_data, int id, int* krowM, int* jcolM, double* M);

extern "C" typedef int (* FVEC)(void* user_data, int id, double* vec, int len);


/** Problem parameters and data */
class ProbData {
public: //data
   int nScenarios;

public: //methods
   ProbData(int nScenarios) {
      this->nScenarios = nScenarios;
   };

//  ~ProbData();
};

extern "C" {


int nSize(void*, int id, int* nnz) {
   if (id == 2)
      *nnz = 4;
   else
      *nnz = 2;

   return 0;
}

int mySize(void*, int, int* nnz) {
   *nnz = 2;

   return 0;
}

int mzSize(void*, int, int* nnz) {
   *nnz = 1;

   return 0;
}


int mylSize(void*, int, int* nnz) {

   *nnz = 2;

   return 0;
}

int mzlSize(void*, int, int* nnz) {
   *nnz = 1;

   return 0;
}


int nnzMatEqStage1(void*, int id, int* nnz) {
   if (id == 0)
      *nnz = 2;
   else
      *nnz = 2;

   return 0;
}


int nnzMatIneqStage1(void*, int id, int* nnz) {
   if (id == 0)
      *nnz = 1;
   else
      *nnz = 1;

   return 0;
}

int nnzMatEqStage2(void*, int id, int* nnz) {
   if (id == 0)
      *nnz = 0;
   else
      *nnz = 2;
   return 0;
}


int nnzMatIneqStage2(void*, int id, int* nnz) {
   if (id == 0)
      *nnz = 0;
   else
      *nnz = 1;
   return 0;
}

int nnzMatEqLink(void*, int id, int* nnz) {
   *nnz = 3;

   if (id == 2)
      *nnz = 4;

   return 0;
}

int nnzMatIneqLink(void*, int, int* nnz) {
   *nnz = 1;

   return 0;
}


int nnzAllZero(void*, int, int* nnz) {
   *nnz = 0;
   return 0;
}

int vecAllZero(void*, int, double* vec, int len) {
   for (int i = 0; i < len; i++)
      vec[i] = 0.0;

   return 0;
}

int vecEqRhs(void*, int id, double* vec, int) {
   if (id == 0) {
      vec[0] = 2.0;
      vec[1] = 7.0;
   }
   else if (id == 1) {
      vec[0] = 3.0;
      vec[1] = 7.0;
   }
   else if (id == 2) {
      vec[0] = 2.0;
      vec[1] = 7.0;
   }

   return 0;
}

int vecIneqRhs(void*, int, double* vec, int len) {
   for (int i = 0; i < len; i++)
      vec[i] = 5.0;

   return 0;
}

int vecIneqRhsLink(void*, int, double* vec, int) {
   vec[0] = 4.0;

   return 0;
}


int vecIneqRhsActive(void*, int, double* vec, int len) {
   for (int i = 0; i < len; i++)
      vec[i] = 1.0;

   return 0;
}

int vecIneqRhsActiveLink(void*, int, double* vec, int) {
   vec[0] = 1.0;
   return 0;
}

int vecObj(void*, int, double* vec, int len) {
   for (int i = 0; i < len; i++)
      vec[i] = 2.0;

   return 0;
}

int vecXlb(void*, int, double* vec, int len) {
   for (int i = 0; i < len; i++)
      vec[i] = 0.0;

   return 0;
}

int vecXlbActive(void*, int, double* vec, int len) {
   for (int i = 0; i < len; i++)
      vec[i] = 1.0;

   return 0;
}

int vecLinkRhs(void*, int, double* vec, int) {
   vec[0] = 6.0;
   vec[1] = 4.0;
   return 0;
}

int matAllZero(void*, int, int*, int*, double*) {
   return 0;
}

int matEqStage1(void*, int id, int* krowM, int* jcolM, double* M) {
   if (id == 0) {
      M[0] = 2.0;
      M[1] = 7.0;

      krowM[0] = 0;
      krowM[1] = 1;
      krowM[2] = 2;


      jcolM[0] = 0;
      jcolM[1] = 1;
   }
   else if (id == 1) {
      M[0] = 2.0;
      M[1] = 5.0;

      krowM[0] = 0;
      krowM[1] = 1;
      krowM[2] = 2;


      jcolM[0] = 0;
      jcolM[1] = 1;
   }
   else if (id == 2) {
      M[0] = 2.0;
      M[1] = 4.0;

      krowM[0] = 0;
      krowM[1] = 1;
      krowM[2] = 2;


      jcolM[0] = 0;
      jcolM[1] = 1;
   }

   return 0;
}

int matIneqLink(void*, int, int* krowM, int* jcolM, double* M) {
   M[0] = 1.0;

   krowM[0] = 0;
   jcolM[0] = 0;
   krowM[1] = 1;

   return 0;
}

int matIneqStage1(void*, int, int* krowM, int* jcolM, double* M) {
   M[0] = 2.0;

   krowM[0] = 0;
   krowM[1] = 1;

   jcolM[0] = 0;

   return 0;
}

int matEqStage2(void*, int id, int* krowM, int* jcolM, double* M) {
   if (id == 0) {
      krowM[0] = 0;
      krowM[1] = 0;
   }
   else if (id == 1) {
      M[0] = 1.0;
      M[1] = 2.0;

      krowM[0] = 0;
      krowM[1] = 1;
      krowM[2] = 2;

      jcolM[0] = 0;
      jcolM[1] = 1;
   }
   else {
      M[0] = 0.0;
      M[1] = 3.0;

      krowM[0] = 0;
      krowM[1] = 1;
      krowM[2] = 2;

      jcolM[0] = 0;
      jcolM[1] = 1;
   }

   return 0;
}

int matIneqStage2(void*, int id, int* krowM, int* jcolM, double* M) {
   if (id == 0) {
      krowM[0] = 0;
      krowM[1] = 0;
   }
   else {
      M[0] = 3.0;

      krowM[0] = 0;
      krowM[1] = 1;

      jcolM[0] = 0;
   }

   return 0;
}

int matEqLink(void*, int id, int* krowM, int* jcolM, double* M) {
   if (id == 2) {
      M[0] = 1.0;
      M[1] = 1.0;
      M[2] = 1.0;
      M[3] = 1.0;

      krowM[0] = 0;
      krowM[1] = 2;
      krowM[2] = 4;

      jcolM[0] = 1;
      jcolM[1] = 3;
      jcolM[2] = 2;
      jcolM[3] = 3;
   }
   else {
      M[0] = 1.0;
      M[1] = 1.0;
      M[2] = 1.0;

      krowM[0] = 0;
      krowM[1] = 2;
      krowM[2] = 3;

      jcolM[0] = 0;
      jcolM[1] = 1;
      jcolM[2] = 0;
   }

   return 0;
}

} /* extern C */


int main(int argc, char** argv) {

   MPI_Init(&argc, &argv);

   int nScenarios = 2;

   // set callbacks
   FNNZ nCall = &nSize;
   FNNZ myCall = &mySize;
   FNNZ mzCall = &mzSize;
   FNNZ mylCall = &mylSize;
   FNNZ mzlCall = &mzlSize;
   FNNZ fnnzQ = &nnzAllZero;
   FNNZ fnnzA = &nnzMatEqStage1;
   FNNZ fnnzB = &nnzMatEqStage2;

   FNNZ fnnzC = &nnzMatIneqStage1;//FNNZ fnnzC = &nnzAllZero;
   FNNZ fnnzD = &nnzMatIneqStage2;//FNNZ fnnzD = &nnzAllZero;

#if LINKING_CONS
   FNNZ fnnzBl = &nnzMatEqLink;
   FVEC fbl = &vecLinkRhs;
   FMAT fBl = &matEqLink;
   FMAT fDl = &matIneqLink;
   FNNZ fnnzDl = &nnzMatIneqLink;

   FVEC fdlupp = &vecIneqRhsLink;
   FVEC fdllow = &vecAllZero;
   FVEC fidlupp = &vecIneqRhsActiveLink;
   FVEC fidllow = &vecAllZero;

#else
   FNNZ fnnzBl = &nnzAllZero;
   FVEC fbl = &vecAllZero;
   FMAT fBl = &matAllZero;

   FMAT fDl = &matAllZero;
   FNNZ fnnzDl = &nnzAllZero;

   FVEC fdlupp = &vecAllZero;
   FVEC fdllow = &vecAllZero;
   FVEC fidlupp = &vecAllZero;
   FVEC fidllow = &vecAllZero;
#endif

   FVEC fc = &vecObj;
   FVEC fb = &vecEqRhs;

   FVEC fclow = &vecAllZero;
   FVEC fcupp = &vecIneqRhs;//FVEC fcupp = &vecAllZero;
   FVEC fxlow = &vecXlb;
   FVEC fxupp = &vecAllZero;
   FVEC ficlow = &vecAllZero;
   FVEC fixlow = &vecXlbActive;
   FVEC ficupp = &vecIneqRhsActive;//FVEC ficupp = vecAllZero;
   FVEC fixupp = &vecAllZero;


   FMAT fQ = &matAllZero;
   FMAT fA = &matEqStage1;
   FMAT fB = &matEqStage2;
   FMAT fC = &matIneqStage1;//FMAT fC = &matAllZero;
   FMAT fD = &matIneqStage2;//FMAT fD = &matAllZero;

   ProbData probData(nScenarios);

   int rank;
   int size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);


#if LINKING_CONS
   //build the problem tree
   DistributedInputTree::DistributedInputNode dataLinkCons(&probData, 0, nCall, myCall, mylCall, mzCall, mzlCall, fQ, fnnzQ, fc, fA, fnnzA, fB, fnnzB, fBl,
         fnnzBl, fb, fbl, fC, fnnzC, fD, fnnzD, fDl, fnnzDl, fclow, ficlow, fcupp, ficupp, fdllow, fidllow, fdlupp, fidlupp, fxlow, fixlow, fxupp,
         fixupp, false);

   DistributedInputTree* root = new DistributedInputTree(dataLinkCons);
#else

   int nx0 = 2;
   int my0 = 2;
   int mz0 = 1;

   //build the problem tree
   DistributedInputTree::DistributedInputNode data(&probData, 0,
                   nx0,my0,mz0, //myl0, mzl0
                   fQ, fnnzQ, fc,
                   fA, fnnzA,
                   fB, fnnzB,
                  //fBL, fnnzBL,
                   fb, //fbl
                   fC, fnnzC,
                   fD, fnnzD,
                  //fDL, fnnzDl,
                   fclow, ficlow, fcupp, ficupp,
                  //fdllow, fidllow, fdlupp, fidlupp,
                   fxlow, fixlow, fxupp, fixupp, false );

   DistributedInputTree* root = new DistributedInputTree(data);
#endif


   for (int id = 1; id <= nScenarios; id++) {
#if LINKING_CONS
      DistributedInputTree::DistributedInputNode dataLinkConsChild(&probData, id, nCall, myCall, mylCall, mzCall, mzlCall, fQ, fnnzQ, fc, fA, fnnzA, fB, fnnzB,
            fBl, fnnzBl, fb, fbl, fC, fnnzC, fD, fnnzD, fDl, fnnzDl, fclow, ficlow, fcupp, ficupp, fdllow, fidllow, fdlupp, fidlupp, fxlow, fixlow,
            fxupp, fixupp, false);

      root->AddChild(new DistributedInputTree(dataLinkConsChild));
#else
      int nx = 2;
      int my = 2;
      int mz = 1;

      DistributedInputTree::DistributedInputNode data(&probData, id,
                nx, my, mz, //myl, mzl
                fQ, fnnzQ, fc,
                fA, fnnzA,
                fB, fnnzB,
                //fBL, fnnzBL,
                fb, //fbl
                fC, fnnzC,
                fD, fnnzD,
                //fDL, fnnzDL,
                fclow, ficlow, fcupp, ficupp,
                //fdllow, fidllow, fdlupp, fidlupp,
                fxlow, fixlow, fxupp, fixupp, false);

      root->AddChild(new DistributedInputTree(data));
#endif

   }

   if (rank == 0)
      std::cout << "Using a total of " << size << " MPI processes." << std::endl;

   /* use BiCGStab for outer solve */
   pips_options::set_int_parameter("INNER_SC_SOLVE", 0);
   PIPSIPMppInterface<InteriorPointMethod> pipsIpm(root, PRIMAL, MPI_COMM_WORLD, SCALER_GEO_STOCH, PRESOLVER_NONE);

   if (rank == 0)
      std::cout << "PIPSIPMppInterface created" << std::endl;

   if (rank == 0)
      std::cout << "solving..." << std::endl;

   pipsIpm.run();

   const double objective = pipsIpm.getObjective();
   if (rank == 0)
      std::cout << "solving finished ... objective value: " << objective << std::endl;

   delete root;

   MPI_Finalize();

   return 0;
}
