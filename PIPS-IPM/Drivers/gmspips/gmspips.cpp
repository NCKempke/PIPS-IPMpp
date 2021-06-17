#if defined(GMS_PIPS)

#include "PIPSIPMppInterface.hpp"
#include "DistributedInputTree.h"

#include "PIPSIPMppOptions.h"
#include "PreprocessType.h"
#include "MehrotraStrategyType.h"

#endif
#if defined(GMS_MPI)

#include "mpi.h"

#endif

#include "gmspipsio.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>

#include <iostream>

extern "C" typedef int (* FNNZ)(void* user_data, int id, int* nnz);

/* Row-major format */
extern "C" typedef int (* FMAT)(void* user_data, int id, int* krowM, int* jcolM, double* M);

extern "C" typedef int (* FVEC)(void* user_data, int id, double* vec, int len);

extern "C" {

static int gmsRank = 0;
static int size = 1;
static bool allGDX = false;
static char fileName[256];
static char GDXDirectory[256];
static char* pGDXDirectory = NULL;
static int numBlocks = 0;
FILE* fLog;

#define checkAndAlloc(blk)                                                        \
if (!blocks[blk])                                                                 \
{                                                                                 \
   int rc;                                                                        \
   fprintf(fLog,"Block %d read on gmsRank %d\n", blk, gmsRank);                   \
   blocks[blk] = (GMSPIPSBlockData_t*) malloc(sizeof(GMSPIPSBlockData_t));        \
   if ( !allGDX )                                                                 \
   {                                                                              \
      char fname[256];                                                            \
      int r = snprintf(fname, 256, "%s%d.gdx", fileName, blk);                    \
      if( r < 0 )  abort();                                                       \
      rc = readBlock(numBlocks,blk,0,1,fname,pGDXDirectory,blocks[blk]);          \
   }                                                                              \
   else                                                                           \
      rc = readBlock(numBlocks,blk,0,1,fileName,pGDXDirectory,blocks[blk]);       \
   if (rc) {fprintf(fLog,"Block %d read on gmsRank %d failed rc=%d\n", blk, gmsRank, rc); return rc;} \
}

#define nCB(nType)                                                   \
int fsize##nType(void* user_data, int id, int* nnz)                  \
{                                                                    \
   GMSPIPSBlockData_t** blocks = (GMSPIPSBlockData_t**) user_data;   \
   checkAndAlloc(id);                                                \
   GMSPIPSBlockData_t* blk = blocks[id];                             \
   assert(blk);                                                      \
   *nnz = blk->nType;                                                \
   fprintf(fLog,"nCB blk=%d " #nType " %d\n",id,*nnz);               \
   return 0;                                                         \
}

nCB(ni)
nCB(mA)
nCB(mC)
nCB(mBL)
nCB(mDL)

#define nnzCB(nnzType)                                               \
int fnonzero##nnzType(void* user_data, int id, int* nnz)             \
{                                                                    \
   GMSPIPSBlockData_t** blocks = (GMSPIPSBlockData_t**) user_data;   \
   checkAndAlloc(id);                                                \
   GMSPIPSBlockData_t* blk = blocks[id];                             \
   assert(blk);                                                      \
   *nnz = blk->nnz##nnzType;                                         \
   fprintf(fLog,"nnzCB blk=%d " #nnzType " %d\n",id,*nnz);           \
   return 0;                                                         \
}

nnzCB(A)
nnzCB(B)
nnzCB(C)
nnzCB(D)
nnzCB(BL)
nnzCB(DL)

#define vecCB(vecType, size)                                         \
int fvec##vecType(void* user_data, int id, double* vec, int len)     \
{                                                                    \
   GMSPIPSBlockData_t** blocks = (GMSPIPSBlockData_t**) user_data;   \
   checkAndAlloc(id);                                                \
   GMSPIPSBlockData_t* blk = blocks[id];                             \
   assert(blk);                                                      \
   fprintf(fLog,"vecCB blk=%d " #vecType " len=%d allocLen %d\n",id,blk->size,len); \
   assert(len == blk->size);                                         \
   for( int i = 0; i < len; i++ ) {                                  \
      vec[i] = blk->vecType[i];                                      \
      fprintf(fLog,"  i=%d val=%g\n",i,vec[i]);                      \
   }                                                                 \
   return 0;                                                         \
}

vecCB(c, ni)
vecCB(xlow, ni)
vecCB(ixlow, ni)
vecCB(xupp, ni)
vecCB(ixupp, ni)
vecCB(b, mA)
vecCB(clow, mC)
vecCB(iclow, mC)
vecCB(cupp, mC)
vecCB(icupp, mC)
vecCB(bL, mBL)
vecCB(dlow, mDL)
vecCB(idlow, mDL)
vecCB(dupp, mDL)
vecCB(idupp, mDL)


#define matCB(mat, mmat)                                                   \
int fmat##mat(void* user_data, int id, int* krowM, int* jcolM, double* M) \
{                                                                         \
   GMSPIPSBlockData_t** blocks = (GMSPIPSBlockData_t**) user_data;        \
   checkAndAlloc(id);                                                     \
   GMSPIPSBlockData_t* blk = blocks[id];                                  \
   assert(blk);                                                           \
   fprintf(fLog,"matCB blk=%d " #mat " mLen %d nzLen %ld\n",id,blk->m##mmat,blk->nnz##mat); \
   if ( 0==blk->m##mmat )                                                 \
   {                                                                      \
     fprintf(fLog," empty\n"); \
     krowM[0] = 0;                                                        \
     return 0;                                                            \
   }                                                                      \
                                                                          \
   assert(blk->rm##mat);                                                  \
   for( int i = 0; i <= blk->m##mmat; i++ ) {                             \
      krowM[i] = blk->rm##mat[i];                                         \
      fprintf(fLog,"  i=%d krowM=%d\n",i,krowM[i]);                       \
   }                                                                      \
                                                                          \
   for( int k = 0; k < blk->nnz##mat; k++ ) {                             \
      jcolM[k] = blk->ci##mat[k];                                         \
      M[k] = blk->val##mat[k];                                            \
      fprintf(fLog,"  k=%d jcolM=%d M=%g\n",k,jcolM[k],M[k]);             \
   }                                                                      \
   return 0;                                                              \
}

matCB(A, A)
matCB(B, A)
matCB(C, C)
matCB(D, C)
matCB(BL, BL)
matCB(DL, DL)

int fnonzeroQ(void*, int, int* nnz) {
   *nnz = 0;
   return 0;
}

int fmatQ(void* user_data, int id, int* krowM, int*, double*) {
   GMSPIPSBlockData_t* blk = ((GMSPIPSBlockData_t**) user_data)[id];
   assert(blk);

   for (int i = 0; i <= blk->ni; i++)
      krowM[i] = 0;

   return 0;
}

}

#if defined(GMS_PIPS)
static void setParams(ScalerType& scaler_type, bool& stepDiffLp, bool& presolve, bool& printsol, bool& hierarchical, const char* paramname) {
   if (strcmp(paramname, "scale") == 0 || strcmp(paramname, "scaleEqui") == 0)
      scaler_type = ScalerType::SCALER_EQUI_STOCH;
   else if (strcmp(paramname, "scaleGeo") == 0)
      scaler_type = ScalerType::SCALER_GEO_STOCH;
   else if (strcmp(paramname, "scaleGeoEqui") == 0)
      scaler_type = ScalerType::SCALER_GEO_EQUI_STOCH;
   else if (strcmp(paramname, "stepLp") == 0)
      stepDiffLp = true;
   else if (strcmp(paramname, "presolve") == 0)
      presolve = true;
   else if (strcmp(paramname, "printsol") == 0)
      printsol = true;
   else if (strcmp(paramname, "hierarchical") == 0)
      hierarchical = true;
}

#endif
int main(int argc, char** argv) {

#if defined(GMS_MPI)
   MPI_Init(&argc, &argv);
   MPI_Barrier(MPI_COMM_WORLD);

   const double t0 = MPI_Wtime();
#endif

   initGMSPIPSIO();

   GMSPIPSBlockData_t** blocks;
#if defined(GMS_PIPS)
   ScalerType scaler_type = ScalerType::SCALER_NONE;
#endif

   bool primal_dual_step_length = false;
   bool presolve = false;
   bool printsol = false;
   bool hierarchical = false;

   if ((argc < 3) || (argc > 9)) {
      std::cout << "Usage: " << argv[0]
                << " numBlocks all.gdx|blockstem [GDXLibDir] [scale] [stepLp] [presolve] [printsol] [hierarchical_approach]\n";
      exit(1);
   }

   allGDX = strstr(argv[2], ".gdx") != nullptr;
   numBlocks = atoi(argv[1]);
   strcpy(fileName, argv[2]);
   if (argc >= 4) {
      strcpy(GDXDirectory, argv[3]);
      pGDXDirectory = &GDXDirectory[0];
   }

#if defined(GMS_PIPS)
   for (int i = 5; i <= argc; i++)
      setParams(scaler_type, primal_dual_step_length, presolve, printsol, hierarchical, argv[i - 1]);
#endif

   blocks = (GMSPIPSBlockData_t**) calloc(numBlocks, sizeof(GMSPIPSBlockData_t*));

   FNNZ fsni = &fsizeni;
   FNNZ fsmA = &fsizemA;
   FNNZ fsmC = &fsizemC;
   FNNZ fsmBL = &fsizemBL;
   FNNZ fsmDL = &fsizemDL;
   FNNZ fnnzQ = &fnonzeroQ;
   FNNZ fnnzA = &fnonzeroA;
   FNNZ fnnzB = &fnonzeroB;
   FNNZ fnnzC = &fnonzeroC;
   FNNZ fnnzD = &fnonzeroD;
   FNNZ fnnzBL = &fnonzeroBL;
   FNNZ fnnzDL = &fnonzeroDL;
   FVEC fc = &fvecc;
   FVEC fxlow = &fvecxlow;
   FVEC fixlow = &fvecixlow;
   FVEC fxupp = &fvecxupp;
   FVEC fixupp = &fvecixupp;
   FVEC fb = &fvecb;
   FVEC fclow = &fvecclow;
   FVEC ficlow = &fveciclow;
   FVEC fcupp = &fveccupp;
   FVEC ficupp = &fvecicupp;
   FVEC fbL = &fvecbL;
   FVEC fdlow = &fvecdlow;
   FVEC fidlow = &fvecidlow;
   FVEC fdupp = &fvecdupp;
   FVEC fidupp = &fvecidupp;

   FMAT fA = &fmatA;
   FMAT fB = &fmatB;
   FMAT fC = &fmatC;
   FMAT fD = &fmatD;
   FMAT fBL = &fmatBL;
   FMAT fDL = &fmatDL;
   FMAT fQ = &fmatQ;

#if defined(GMS_PIPS)
   //build the problem tree
   DistributedInputTree::DistributedInputNode data(blocks, 0,
         fsni, fsmA, fsmBL, fsmC, fsmDL,
         fQ, fnnzQ, fc, fA, fnnzA, fB, fnnzB,
         fBL, fnnzBL,
         fb,
         fbL,
         fC, fnnzC, fD, fnnzD,
         fDL, fnnzDL,
         fclow, ficlow, fcupp, ficupp,
         fdlow, fidlow, fdupp, fidupp,
         fxlow, fixlow, fxupp, fixupp, false);
   std::unique_ptr<DistributedInputTree> root = std::make_unique<DistributedInputTree>(data);
#endif
   for (int blk = 1; blk < numBlocks; blk++) {

#if defined(GMS_PIPS)
      DistributedInputTree::DistributedInputNode data(blocks, blk,
            fsni, fsmA, fsmBL, fsmC, fsmDL,
            fQ, fnnzQ, fc, fA, fnnzA, fB, fnnzB,
            fBL, fnnzBL,
            fb,
            fbL,
            fC, fnnzC, fD, fnnzD,
            fDL, fnnzDL,
            fclow, ficlow, fcupp, ficupp,
            fdlow, fidlow, fdupp, fidupp,
            fxlow, fixlow, fxupp, fixupp, false);

      root->AddChild(new DistributedInputTree(data));
#endif
   }

#if defined(GMS_MPI)
   MPI_Comm_rank(MPI_COMM_WORLD, &gmsRank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   char fbuf[256];
#if defined (GMS_LOG)
   sprintf(fbuf,"log%d.txt", gmsRank);
#else
   sprintf(fbuf, "/dev/null");
#endif
   fLog = fopen(fbuf, "w+");
   fprintf(fLog, "PIPS Log for gmsRank %d\n", gmsRank);
#endif

   if (gmsRank == 0)
      std::cout << "Using a total of " << size << " MPI processes.\n";


#if defined(GMS_PIPS)
   if (hierarchical) {
      if (gmsRank == 0)
         std::cout << "Using Hierarchical approach\n";
      pipsipmpp_options::activate_hierarchial_approach();
   }

   pipsipmpp_options::set_int_parameter("OUTER_SOLVE", 2);
   if (gmsRank == 0)
      std::cout << "Using outer BICGSTAB\n";

   if (gmsRank == 0 && pipsipmpp_options::get_int_parameter("INNER_SC_SOLVE") == 2)
      std::cout << "Using inner BICGSTAB\n";

   std::vector<double> primalSolVec;
   std::vector<double> dualSolEqVec;
   std::vector<double> dualSolIneqVec;
   std::vector<double> dualSolVarBounds;

   std::vector<double> eqValues;
   std::vector<double> ineqValues;

   pipsipmpp_options::set_bool_parameter("GONDZIO_ADAPTIVE_LINESEARCH", !primal_dual_step_length);
   if (primal_dual_step_length && gmsRank == 0) {
      std::cout << "Different steplengths in primal and dual direction are used.\n";
   }

   // create the PIPS-IPM++ interface
   PIPSIPMppInterface pipsIpm(root.get(), primal_dual_step_length ? MehrotraStrategyType::PRIMAL_DUAL : MehrotraStrategyType::PRIMAL, MPI_COMM_WORLD, scaler_type,
         presolve ? PresolverType::PRESOLVER_STOCH : PresolverType::PRESOLVER_NONE);

   if (gmsRank == 0) {
      std::cout << "PIPSIPMppInterface created\n";
      std::cout << "solving...\n";
   }

   // run PIPS-IPM++
   pipsIpm.run();
   double objective = pipsIpm.getObjective();

   if (presolve) {
      pipsIpm.postsolveComputedSolution();
   }
   if (printsol) {
      primalSolVec = pipsIpm.gatherPrimalSolution();
      dualSolEqVec = pipsIpm.gatherDualSolutionEq();
      dualSolIneqVec = pipsIpm.gatherDualSolutionIneq();
      dualSolVarBounds = pipsIpm.gatherDualSolutionVarBounds();
      eqValues = pipsIpm.gatherEqualityConsValues();
      ineqValues = pipsIpm.gatherInequalityConsValues();
   }


   if (gmsRank == 0)
      std::cout << "solving finished. \n ---Objective value: " << objective << "\n";

   if (printsol && gmsRank == 0) {
      int rc;

      rc = writeSolution(fileName, primalSolVec.size(), dualSolEqVec.size(), dualSolIneqVec.size(), objective, &primalSolVec[0], &dualSolVarBounds[0],
            &eqValues[0], &ineqValues[0], &dualSolEqVec[0], &dualSolIneqVec[0], pGDXDirectory);
      if (0 == rc)
         std::cout << "Solution written to " << fileName << "_sol.gdx\n";
      else if (-1 == rc)
         std::cout << "Could not access " << fileName << ".map\n";
      else
         std::cout << "Other error writing solution: rc=" << rc << "\n";
   }
#endif


   for (int blk = 0; blk < numBlocks; blk++) {
      freeBlock(blocks[blk]);
      free(blocks[blk]);
   }
   free(blocks);

#if defined(GMS_MPI)
   MPI_Barrier(MPI_COMM_WORLD);
   const double t1 = MPI_Wtime();

   if (gmsRank == 0)
      std::cout << "---total time (in sec.): " << t1 - t0 << "\n";

   MPI_Finalize();
#endif
   fclose(fLog);

   return 0;
}
