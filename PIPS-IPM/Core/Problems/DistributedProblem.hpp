#ifndef DISTRIBUTEDQP_H
#define DISTRIBUTEDQP_H

#include "Problem.hpp"
#include "DistributedResiduals.hpp"
#include "DistributedVariables.h"
#include "DistributedSymmetricMatrix.h"
#include "SparseSymmetricMatrix.h"
#include "DistributedMatrix.h"
#include "DistributedVector.h"
#include "PIPSIPMppOptions.h"
#include "DistributedTreeCallbacks.h"
#include "pipschecks.h"

#include <vector>
#include <memory>

class DistributedTree;

class DistributedProblem : public Problem {
protected:
   DistributedProblem() = default;

public:
   /** constructor that sets up pointers to the data objects that are passed as arguments */
   DistributedProblem(const DistributedTree* stochNode, std::shared_ptr<Vector<double>> c,
      std::shared_ptr<SymmetricMatrix> Q, std::shared_ptr<Vector<double>> xlow,
      std::shared_ptr<Vector<double>> ixlow, std::shared_ptr<Vector<double>> xupp,
      std::shared_ptr<Vector<double>> ixupp, std::shared_ptr<GeneralMatrix> A, std::shared_ptr<Vector<double>> bA,
      std::shared_ptr<GeneralMatrix> C, std::shared_ptr<Vector<double>> clow, std::shared_ptr<Vector<double>> iclow,
      std::shared_ptr<Vector<double>> cupp,
      std::shared_ptr<Vector<double>> icupp,
      std::shared_ptr<Vector<double>> integrality,
      bool add_children = true, bool is_hierarchy_root = false,
      bool is_hierarchy_inner_root = false,
      bool is_hierarchy_inner_leaf = false);

   std::vector<DistributedProblem*> children;

   void add_child(DistributedProblem* child);

   const DistributedTree* stochNode{};

   [[nodiscard]] Permutation getLinkVarsPermInv() const;

   [[nodiscard]] Permutation getLinkConsEqPermInv() const;

   [[nodiscard]] Permutation getLinkConsIneqPermInv() const;

public:
   [[nodiscard]] int getLocalnx() const;

   [[nodiscard]] int getLocalmy() const;

   [[nodiscard]] int getLocalmyl() const;

   [[nodiscard]] int getLocalmz() const;

   [[nodiscard]] int getLocalmzl() const;

   void getLocalSizes(int& nx, int& my, int& mz, int& myl, int& mzl) const;

   void getLocalNnz(int& nnzQ, int& nnzB, int& nnzD) const;

   [[nodiscard]] int getN0LinkVars() const { return n0LinkVars; }

   // returns upper bound on number of non-zeroes in Schur complement
   int getSchurCompMaxNnz() const;

   // distributed version
   int getSchurCompMaxNnzDist(int blocksStart, int blocksEnd) const;

   [[nodiscard]] bool exploitingLinkStructure() const { return useLinkStructure; };

   std::unique_ptr<SparseSymmetricMatrix> createSchurCompSymbSparseUpper() const;

   // distributed version
   std::unique_ptr<SparseSymmetricMatrix> createSchurCompSymbSparseUpperDist(int blocksStart, int blocksEnd) const;

   const SparseSymmetricMatrix& getLocalQ() const;

   const SparseMatrix& getLocalCrossHessian() const;

   const SparseMatrix& getLocalA() const;

   const SparseMatrix& getLocalB() const;

   const SparseMatrix& getLocalF() const;

   const SparseMatrix& getLocalC() const;

   const SparseMatrix& getLocalD() const;

   const SparseMatrix& getLocalG() const;

   const StripMatrix& getLocalGBorder() const;

   const StripMatrix& getLocalFBorder() const;


   void printLinkVarsStats();

   void printLinkConsStats();

   void activateLinkStructureExploitation();

   DistributedResiduals* getResidsUnperm(const Residuals& resids, const Problem& unpermData) const;

   DistributedVariables* getVarsUnperm(const Variables& vars, const Problem& unpermData) const;

   bool isRootNodeInSync() const;

protected:
   static void removeN0LinkVarsIn2Links(std::vector<int>& n_blocks_per_link_var, const DistributedMatrix& Astoch,
      const DistributedMatrix& Cstoch,
      const std::vector<int>& linkStartBlockIdA, const std::vector<int>& linkStartBlockIdC);

   static Permutation
   getChildLinkConsFirstOwnLinkConsLastPermutation(const std::vector<unsigned int>& map_block_subtree,
      const std::vector<int>& linkStartBlockId,
      int n_links_after_split);

   void addChildrenForSplit();

   void splitData();

   void reorderLinkingConstraintsAccordingToSplit();

   void splitDataAndAddAsChildLayer();

   DistributedProblem* shaveBorderFromDataAndCreateNewTop(const DistributedTree& tree);

public:
   DistributedProblem* shaveDenseBorder(const DistributedTree& tree);

   void splitDataAccordingToTree();

   void recomputeSize();

   void splitStringMatricesAccordingToSubtreeStructure();

   int getNGlobalVars() const { return n_global_linking_vars; };

   int getNGlobalEQConss() const { return n_global_eq_linking_conss; };

   int getNGlobalINEQConss() const { return n_global_ineq_linking_conss; };

   void write_to_streamDense(std::ostream& out) const override;

   std::unique_ptr<Problem> clone_full() const override {
      return clone_full(false);
   }
   std::unique_ptr<DistributedProblem> clone_full(bool switchToDynamicStorage = false) const;

   double evaluate_objective(const Variables& variables) const override;

   void
   cleanUpPresolvedData(const DistributedVector<int>& rowNnzVecA, const DistributedVector<int>& rowNnzVecC,
      const DistributedVector<int>& colNnzVec);

   // marker that indicates whether a Schur complement row is (2-link) local
   const std::vector<bool>& getSCrowMarkerLocal() const;

   // marker that indicates whether a Schur complement row is (2-link) local and owned by this MPI process
   const std::vector<bool>& getSCrowMarkerMyLocal() const;

   // number of sparse 2-link equality rows
   int n2linkRowsEq() const;

   // number of sparse 2-link inequality rows
   int n2linkRowsIneq() const;

   // start and end positions for local 2-links in Schur complement that are non-zero if only
   // blocks greater equal blocksStart and smaller blocksEnd are considered
   void getSCrangeMarkers(int blocksStart, int blocksEnd, int& local2linksStartEq, int& local2linksEndEq,
      int& local2linksStartIneq,
      int& local2linksEndIneq) const;

   // start and end positions for local 2-links in Schur complement that are owned by
   // blocks greater equal blocksStart and smaller blocksEnd are considered
   void getSCrangeMarkersMy(int blocksStart, int blocksEnd, int& local2linksStartEq, int& local2linksEndEq,
      int& local2linksStartIneq,
      int& local2linksEndIneq) const;

   bool isHierarchySparseTopLayerOnlyTwolinks() const {
      return (pipsipmpp_options::get_bool_parameter("HIERARCHICAL")) &&
         (pipsipmpp_options::get_int_parameter("HIERARCHICAL_APPROACH_N_LAYERS") > 1) &&
         threshold_global_cons <= 2 && threshold_global_vars <= 0;
   };

   bool isHierarchyRoot() const { return is_hierarchy_root; };

   bool isHierarchyInnerRoot() const { return is_hierarchy_inner_root; };

   bool isHierarchyInnerLeaf() const { return is_hierarchy_inner_leaf; };

   bool hasRAC() const { return has_RAC; };

   ~DistributedProblem() override;

protected:
   void createChildren();

   void destroyChildren();

private:
   int n0LinkVars{0};

   /* trimming everything that is not a 2 link and everything that is not a 0vec at the moment - makes border computations sparser */
   /* a constraint is global if it is in more than threshold_global_cons blocks */
   constexpr static int threshold_global_cons{2};
   constexpr static bool shave_nonlocal_2links{true};
   constexpr static int threshold_global_vars{1}; // TODO: adapt properly
   constexpr static int nLinkStats{6};
   constexpr static double minStructuredLinksRatio{0.5};

   static Permutation get0VarsLastGlobalsFirstPermutation(std::vector<int>& linkVarsNnzCount, int& n_globals);

   static Permutation
   getAscending2LinkFirstGlobalsLastPermutation(std::vector<int>& linkStartBlockId, std::vector<int>& n_blocks_per_row,
      size_t nBlocks,
      int& n_globals);

   // returns number of block rows
   static int getSCdiagBlocksNRows(const std::vector<int>& linkStartBlockLengths);

   // returns number of non-zero block rows within specified range
   static int getSCdiagBlocksNRows(const std::vector<int>& linkStartBlockLengths, int blocksStart, int blocksEnd);

   // returns number of owned block rows within specified range
   static int getSCdiagBlocksNRowsMy(const std::vector<int>& linkStartBlockLengths, int blocksStart, int blocksEnd);

   // max nnz in Schur complement diagonal block signified by given vector
   static int getSCdiagBlocksMaxNnz(size_t nRows, const std::vector<int>& linkStartBlockLengths);

   // distributed version
   static int getSCdiagBlocksMaxNnzDist(size_t nRows, const std::vector<int>& linkStartBlockLengths, int blocksStart,
      int blocksEnd);

   // max nnz in Schur complement mixed block signified by given vectors
   static int getSCmixedBlocksMaxNnz(size_t nRows, size_t nCols, const std::vector<int>& linkStartBlockLength_Left,
      const std::vector<int>& linkStartBlockLength_Right);

   // distributed version
   static int getSCmixedBlocksMaxNnzDist(size_t nRows, size_t nCols, const std::vector<int>& linkStartBlockLength_Left,
      const std::vector<int>& linkStartBlockLength_Right, int blocksStart, int blocksEnd);

   // number of sparse 2-link rows
   static int n2linksRows(const std::vector<int>& linkStartBlockLengths);

   static std::vector<int> get2LinkLengthsVec(const std::vector<int>& linkStartBlocks, const size_t nBlocks);

   /* a two link must be in two blocks directly after one another */
   const bool is_hierarchy_root{false};
   bool is_hierarchy_inner_root{false};
   bool is_hierarchy_inner_leaf{false};
   /* only for leafs: indicate whether A_mat belongs to the child (true) or is part of the border (false) */
   bool has_RAC{true};

   bool useLinkStructure{false};

   int n_global_linking_vars{-1};
   int n_global_eq_linking_conss{-1};
   int n_global_ineq_linking_conss{-1};

   /* number non-empty of blocks for each linking column */
   std::vector<int> n_blocks_per_link_var;

   /* which blocks do the individual two-links start in */
   std::vector<int> linkStartBlockIdA;
   std::vector<int> n_blocks_per_link_row_A;

   std::vector<int> linkStartBlockIdC;
   std::vector<int> n_blocks_per_link_row_C;

public:
   const std::vector<int>& getTwoLinksStartBlockA() { return linkStartBlockLengthsA; }

   const std::vector<int>& getTwoLinksStartBlockC() { return linkStartBlockLengthsC; }

private:
   /* how many two-links start in block i */
   std::vector<int> linkStartBlockLengthsA;
   std::vector<int> linkStartBlockLengthsC;

   Permutation linkVarsPermutation;
   Permutation linkConsPermutationA;
   Permutation linkConsPermutationC;

   /* get initialized lazily */
   mutable std::vector<bool> isSCrowLocal;
   mutable std::vector<bool> isSCrowMyLocal;

   void initDistMarker(int blocksStart, int blocksEnd) const;

   void permuteLinkStructureDetection(const Permutation& perm_A, const Permutation& perm_C);

   void permuteLinkingVars(const Permutation& perm);

   void permuteLinkingCons(const Permutation& permA, const Permutation& permC);
};


#endif
