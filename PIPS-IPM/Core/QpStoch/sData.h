#ifndef DATAQPSTOCH
#define DATAQPSTOCH

#include "QpGenData.h"
#include "sResiduals.h"
#include "sVars.h"
#include "StochSymMatrix.h"
#include "SparseSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "DoubleMatrixHandle.h"

#include "sTreeCallbacks.h"

#include "pipschecks.h"
#include "pipsport.h"

#include <vector>
#include <memory>

class sTree;
class LinearAlgebraPackage;

class sData : public QpGenData {
   protected:
      sData() = default;
   public:
  /** constructor that sets up pointers to the data objects that are
      passed as arguments */
  sData( const sTree* stochNode,
	 OoqpVector * c, SymMatrix * Q,
	 OoqpVector * xlow, OoqpVector * ixlow,
	 OoqpVector * xupp, OoqpVector * ixupp,
	 GenMatrix * A, OoqpVector * bA,
	 GenMatrix * C,
	 OoqpVector * clow, OoqpVector * iclow,
	 OoqpVector * cupp, OoqpVector * ciupp,
	 bool add_children = true, bool is_hierarchy_root = false, bool is_hierarchy_inner_root = false,
    bool is_hierarchy_inner_leaf = false );

  std::vector<sData*> children;
  void AddChild(sData* child);
  const sTree* stochNode{};

  PERMUTATION getLinkVarsPermInv() const;
  PERMUTATION getLinkConsEqPermInv() const;
  PERMUTATION getLinkConsIneqPermInv() const;

public:
  int getLocalnx() const;
  int getLocalmy() const;
  int getLocalmyl() const;
  int getLocalmz() const;
  int getLocalmzl() const;
  int getLocalSizes(int& nx, int& my, int& mz) const;
  int getLocalSizes(int& nx, int& my, int& mz, int& myl, int& mzl) const;

  int getLocalNnz(int& nnzQ, int& nnzB, int& nnzD);
  int getN0LinkVars() { return n0LinkVars; }

  // returns upper bound on number of non-zeroes in Schur complement
  int getSchurCompMaxNnz();

  // distributed version
  int getSchurCompMaxNnzDist(int blocksStart, int blocksEnd);
  bool exploitingLinkStructure() { return useLinkStructure; };

  SparseSymMatrix* createSchurCompSymbSparseUpper();

  // distributed version
  SparseSymMatrix* createSchurCompSymbSparseUpperDist(int blocksStart, int blocksEnd);

  SparseSymMatrix& getLocalQ();
  SparseGenMatrix& getLocalCrossHessian();
  SparseGenMatrix& getLocalA();
  SparseGenMatrix& getLocalB();
  SparseGenMatrix& getLocalF();
  SparseGenMatrix& getLocalC();
  SparseGenMatrix& getLocalD();
  SparseGenMatrix& getLocalG();
  StringGenMatrix& getLocalGBorder();
  StringGenMatrix& getLocalFBorder();


  void printLinkVarsStats();
  void printLinkConsStats();

  void activateLinkStructureExploitation();

  sResiduals* getResidsUnperm(const sResiduals& resids, const sData& unpermData) const;
  sVars* getVarsUnperm(const sVars& vars, const sData& unpermData) const;

  bool isRootNodeInSync() const;

protected:
  static void removeN0LinkVarsIn2Links( std::vector<int>& n_blocks_per_link_var, const StochGenMatrix& Astoch,
        const StochGenMatrix& Cstoch, const std::vector<int>& linkStartBlockIdA,
        const std::vector<int>& linkStartBlockIdC );

  static PERMUTATION getChildLinkConsFirstOwnLinkConsLastPermutation( const std::vector<unsigned int>& map_block_subtree,
        const std::vector<int>& linkStartBlockId, int n_links_after_split );

  void addChildrenForSplit();
  void splitData();
  void reorderLinkingConstraintsAccordingToSplit();
  void splitDataAndAddAsChildLayer();

  sData* shaveBorderFromDataAndCreateNewTop( const sTree* tree );

 public:
  sData* shaveDenseBorder( const sTree* tree );
  void splitDataAccordingToTree();
  void recomputeSize();
  void splitStringMatricesAccordingToSubtreeStructure();

  int getNGlobalVars() const { return n_global_linking_vars; };
  int getNGlobalEQConss() const { return n_global_eq_linking_conss; };
  int getNGlobalINEQConss() const { return n_global_ineq_linking_conss; };

  virtual void writeToStreamDense( std::ostream& out ) const;
  void writeMPSformat( std::ostream& out );
  void writeMPSColumns( std::ostream& out );

  virtual sData* cloneFull(bool switchToDynamicStorage = false) const;

  double objectiveValue( const QpGenVars * vars ) const override;
  void createScaleFromQ();

  void cleanUpPresolvedData(const StochVectorBase<int>& rowNnzVecA, const StochVectorBase<int>& rowNnzVecC, const StochVectorBase<int>& colNnzVec);

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
        int& local2linksStartIneq, int& local2linksEndIneq);

  // start and end positions for local 2-links in Schur complement that are owned by
  // blocks greater equal blocksStart and smaller blocksEnd are considered
  void getSCrangeMarkersMy(int blocksStart, int blocksEnd, int& local2linksStartEq, int& local2linksEndEq,
        int& local2linksStartIneq, int& local2linksEndIneq);

  bool isHierarchySparseTopLayerOnlyTwolinks() const { return pips_options::getBoolParameter("HIERARCHICAL") && threshold_global_cons <= 1 && threshold_global_vars == 0; };

  bool isHierarchyRoot() const { return is_hierarchy_root; };
  bool isHierarchyInnerRoot() const { return is_hierarchy_inner_root; };
  bool isHierarchyInnerLeaf() const { return is_hierarchy_inner_leaf; };
  bool hasRAC() const { return has_RAC; };

  ~sData() override;

 protected:
  void createChildren();
  void destroyChildren();

 private:
  int n0LinkVars{0};

  /* trimming everything that is not a 2 link and everything that is not a 0vec at the moment - makes border computations sparser */
  constexpr static int threshold_global_cons{1};
  constexpr static int threshold_global_vars{0}; // TODO: adapt properly
  constexpr static int nLinkStats{6};
  constexpr static double minStructuredLinksRatio{0.5};
  static PERMUTATION get0VarsLastGlobalsFirstPermutation(std::vector<int>& linkVarsNnzCount, int& n_globals);
  static PERMUTATION getAscending2LinkFirstGlobalsLastPermutation(std::vector<int>& linkStartBlockId,
        std::vector<int>& n_blocks_per_row, size_t nBlocks, int& n_globals);

  // returns number of block rows
  static int getSCdiagBlocksNRows(const std::vector<int>& linkStartBlockLengths);

  // returns number of non-zero block rows within specified range
  static int getSCdiagBlocksNRows(const std::vector<int>& linkStartBlockLengths,
        int blocksStart, int blocksEnd);

  // returns number of owned block rows within specified range
  static int getSCdiagBlocksNRowsMy(const std::vector<int>& linkStartBlockLengths,
        int blocksStart, int blocksEnd);

  // max nnz in Schur complement diagonal block signified by given vector
  static int getSCdiagBlocksMaxNnz(size_t nRows,
        const std::vector<int>& linkStartBlockLengths);

  // distributed version
  static int getSCdiagBlocksMaxNnzDist(size_t nRows,
        const std::vector<int>& linkStartBlockLengths, int blocksStart, int blocksEnd);

  // max nnz in Schur complement mixed block signified by given vectors
  static int  getSCmixedBlocksMaxNnz(size_t nRows, size_t nCols,
        const std::vector<int>& linkStartBlockLength_Left,
        const std::vector<int>& linkStartBlockLength_Right);

  // distributed version
  static int  getSCmixedBlocksMaxNnzDist(size_t nRows, size_t nCols,
        const std::vector<int>& linkStartBlockLength_Left,
        const std::vector<int>& linkStartBlockLength_Right, int blocksStart, int blocksEnd);

  // number of sparse 2-link rows
  static int n2linksRows(const std::vector<int>& linkStartBlockLengths);

  static std::vector<int> get2LinkLengthsVec(const std::vector<int>& linkStartBlocks, const size_t nBlocks);

  /* a two link must be in two blocks directly after one another */
  const bool is_hierarchy_root{false};
  bool is_hierarchy_inner_root{false};
  bool is_hierarchy_inner_leaf{false};
  /* only for leafs: indicate wether A_mat stroed in this child belongs to the child (true) or is part of the border (false) */
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
  const std::vector<int>& getTwoLinksStartBlockA()
  { return linkStartBlockLengthsA; }
  const std::vector<int>& getTwoLinksStartBlockC()
  { return linkStartBlockLengthsC; }

 private:
  /* how many two-links start in block i */
  std::vector<int> linkStartBlockLengthsA;
  std::vector<int> linkStartBlockLengthsC;

  PERMUTATION linkVarsPermutation;
  PERMUTATION linkConsPermutationA;
  PERMUTATION linkConsPermutationC;
  std::vector<bool> isSCrowLocal;
  std::vector<bool> isSCrowMyLocal;

  void initDistMarker(int blocksStart, int blocksEnd);

  void permuteLinkStructureDetection( const PERMUTATION& perm_A, const PERMUTATION& perm_C );
  void permuteLinkingVars( const PERMUTATION& perm );
  void permuteLinkingCons( const PERMUTATION& permA, const PERMUTATION& permC );
public:
  void printRanges() const;
};


#endif
