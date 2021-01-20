#ifndef STOCHGENMATRIX_H
#define STOCHGENMATRIX_H

#include "StochVector_fwd.h"
#include "OoqpVector_fwd.h"
#include "DoubleMatrix.h"
#include "SparseGenMatrix.h"
#include "BorderedGenMatrix.h"
#include "StringGenMatrix.h"
#include "sTree.h"

#include "pipsport.h"
#include "mpi.h"

#include <vector>

class StochGenMatrix : public GenMatrix {
protected:
  StochGenMatrix() = default;
public:

  StochGenMatrix( GenMatrix* Amat, GenMatrix* Bmat, GenMatrix* Blmat, MPI_Comm mpiComm_, bool inner_leaf = false, bool inner_root = false );

  /** Constructs a matrix having local A and B blocks having the sizes and number of nz specified by
   *  A_m, A_n, A_nnz and B_m, B_n, B_nnz.
   *  Also sets the global sizes to 'global_m' and 'global_n'.
   *  The matrix that will be created  has no children, just local data.
   */
  StochGenMatrix(long long global_m, long long global_n,
		 int A_m, int A_n, int A_nnz,
		 int B_m, int B_n, int B_nnz,
		 MPI_Comm mpiComm_ );

  /** Constructs a matrix with local A, B, and Bl (linking constraints) blocks having the sizes and number of nz specified by
      A_m, A_n, A_nnz, B_m, B_n, B_nnz, and Bl_m, Bl_n, Bl_nnz. Otherwise, identical to the above constructor */
  StochGenMatrix(long long global_m, long long global_n,
		 int A_m, int A_n, int A_nnz,
		 int B_m, int B_n, int B_nnz,
		 int Bl_m, int Bl_n, int Bl_nnz,
		 MPI_Comm mpiComm_ );

  /** Constructs a matrix with local A, B, and Bl (linking constraints) blocks set to nullptr */
  StochGenMatrix(long long global_m, long long global_n, MPI_Comm mpiComm_);

  // constructor for combining scenarios
  virtual ~StochGenMatrix();

  GenMatrix* cloneEmptyRows( bool switchToDynamicStorage = false ) const override;
  GenMatrix* cloneFull( bool switchToDynamicStorage = false ) const override;

  virtual void AddChild( StochGenMatrix* child );

  std::vector<StochGenMatrix*> children;
  GenMatrix* Amat{};
  GenMatrix* Bmat{};
  GenMatrix* Blmat{};

  long long m{-1};
  long long n{-1};
  MPI_Comm mpiComm{MPI_COMM_NULL};
  int iAmDistrib{false};

  /* is this matrix an inner matrix of the matrix hierarchy - if not, then its children hold the local Amat, Bmat and Blmat */
  const bool inner_leaf{false};
  const bool inner_root{false};
 private:
  bool hasSparseMatrices() const;

  /** trans mult method for children with linking constraints */
  virtual void transMult2( double beta, StochVector& y, double alpha, StochVector& x, const OoqpVector* xvecl ) const;

  virtual void mult2( double beta,  StochVector& y, double alpha, StochVector& x, OoqpVector* yparentl_ );

  /** column scale method for children */
  virtual void columnScale2( const OoqpVector& vec );

  /** row scale method for children */
  virtual void rowScale2( const OoqpVector& vec, const OoqpVector* linkingvec );

  virtual void getNnzPerRow(OoqpVectorBase<int>& nnzVec, OoqpVectorBase<int>* linkParent);

  virtual void getNnzPerCol(OoqpVectorBase<int>& nnzVec, OoqpVectorBase<int>* linkParent);

  virtual void addRowSums( OoqpVector& sumVec, OoqpVector* linkParent ) const;
  virtual void addColSums( OoqpVector& sumVec, OoqpVector* linkParent ) const;

  /** internal method needed for handling linking constraints */
  virtual void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, const OoqpVector* colScaleParent, OoqpVector& minmaxVec, OoqpVector* linkParent);

  virtual void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, const OoqpVector* rowScaleParent, OoqpVector& minmaxVec, OoqpVector* minmaxParent );


  virtual void initTransposedChild(bool dynamic);
  virtual void initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec,
    const OoqpVectorBase<int>* rowLinkVec, const OoqpVectorBase<int>* colParentVec);

  virtual void permuteLinkingVarsChild(const std::vector<unsigned int>& permvec);

  virtual void getLinkVarsNnzChild(std::vector<int>& vec) const;

 public:
  virtual void updateTransposed();

  void getSize( long long& m, long long& n ) const override;
  void getSize( int& m, int& n ) const override;

  /** The actual number of structural non-zero elements in this sparse
   *  matrix. This includes so-called "accidental" zeros, elements that
   *  are treated as non-zero even though their value happens to be zero.
   */  
  int numberOfNonZeros() const override;

  int isKindOf( int matType ) const override;

  void atPutDense( int, int, double*, int, int, int ) { assert( "Not implemented" && 0 ); };
  void fromGetDense( int, int, double*, int, int, int ) { assert( "Not implemented" && 0 ); };

  void columnScale( const OoqpVector& vec ) override;
  void rowScale( const OoqpVector& vec ) override;
  void symmetricScale( const OoqpVector& ) override { assert( "Not implemented" && 0 ); };
  void scalarMult( double num) override;

  void fromGetSpRow( int, int, double[], int, int[], int&, int, int& ) override { assert( "Not implemented" && 0 ); };
  void atPutSubmatrix( int, int, DoubleMatrix&, int, int, int, int ) override { assert( "Not implemented" && 0 ); };
  void atPutSpRow( int, double[], int, int[], int& ) override { assert( "Not implemented" && 0 ); };
  void putSparseTriple( int[], int, int[], double[], int& ) override { assert( "Not implemented" && 0 ); };

  void getDiagonal( OoqpVector& vec ) override;
  void setToDiagonal( const OoqpVector& vec ) override;

  /** y = beta * y + alpha * this * x */
  void mult ( double beta,  OoqpVector& y,
                      double alpha, const OoqpVector& x ) const override;

  void transMult ( double beta,   OoqpVector& y,
                           double alpha,  const OoqpVector& x ) const override;

  double abmaxnorm() const override;
  double abminnormNonZero( double tol = 1e-30 ) const override;

  virtual void getLinkVarsNnz(std::vector<int>& vec) const;

  void writeToStream( std::ostream& ) const override { assert( "Not implemented" && 0 ); };

  void writeToStreamDense( std::ostream& out ) const override
  { writeToStreamDense( out, 0 ); };
  virtual void writeToStreamDense( std::ostream& out, int offset ) const;
  virtual void writeToStreamDenseBordered( const StringGenMatrix& border, std::ostream& out, int offset = 0 )const;

  void writeMPSformatRows( std::ostream& out, int rowType, OoqpVector* irhs) const override;

  void randomize( double, double, double* ) override { assert( "Not implemented" && 0 ); };

  /** initialize (dynamic) transposed matrices for A, B, Bl */
  virtual void initTransposed(bool dynamic = false);
  virtual void deleteTransposed();

  void atPutDiagonal( int, OoqpVector& ) override { assert( "Not implemented" && 0 ); };
  void fromGetDiagonal( int, OoqpVector& ) override { assert( "Not implemented" && 0 ); };
  void matTransDMultMat( OoqpVector&, SymMatrix** ) override { assert( "Not implemented" && 0 ); };
  void matTransDinvMultMat( OoqpVector&, SymMatrix** ) override { assert( "Not implemented" && 0 ); };

  void getNnzPerRow(OoqpVectorBase<int>& nnzVec) override
  {
     getNnzPerRow(nnzVec, nullptr);
  };

  void getNnzPerCol(OoqpVectorBase<int>& nnzVec) override
  {
     getNnzPerCol(nnzVec, nullptr);
  };

  /** fill vector with absolute minimum/maximum value of each row */
  void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, OoqpVector& minmaxVec ) override
  {
     getRowMinMaxVec(getMin, initializeVec, colScaleVec, nullptr, minmaxVec, nullptr);
  };

  /** fill vector with absolute minimum/maximum value of each column */
  void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, OoqpVector& minmaxVec ) override
  {
     getColMinMaxVec(getMin, initializeVec, rowScaleVec, nullptr, minmaxVec, nullptr);
  };

  void addRowSums( OoqpVector& sumVec ) const override
     { addRowSums(sumVec, nullptr); };
  void addColSums( OoqpVector& sumVec ) const override
     { addColSums(sumVec, nullptr); };

  virtual void initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec)
  {
     initStaticStorageFromDynamic(rowNnzVec, colNnzVec, nullptr, nullptr);
  };
  virtual void freeDynamicStorage();
  virtual void recomputeSize( StochGenMatrix* parent = nullptr );

  /** returns Simple Vector indicating which linking rows have entries in exactly two blocks (indicated by 1.0 versus 0.0)*/
  virtual void get2LinkStartBlocksAndCountsNew( std::vector<int>& block_start, std::vector<int>& block_count) const;
  virtual std::vector<int> get2LinkStartBlocks() const;

  virtual void updateKLinkVarsCount(std::vector<int>& linkCount) const;
  virtual void updateKLinkConsCount(std::vector<int>& linkCount) const;

  virtual void permuteLinkingVars(const std::vector<unsigned int>& permvec);
  virtual void permuteLinkingCons(const std::vector<unsigned int>& permvec);

  virtual bool isRootNodeInSync() const;

  virtual int appendRow( const StochGenMatrix& matrix_row, int child, int row, bool linking );

  /** calculate vec^T * row where row is linking or not in child child and with row index row
   *  for a linking row only the available blocks will be multiplied - currently only possible for dynamic storage! (since 
   *  this was its foremost usecase)
   */
  virtual double localRowTimesVec( const StochVector& vec, int child, int row, bool linking ) const;

  /* y += alpha * RowAt(child, row, linking) */
  virtual void axpyWithRowAt( double alpha, StochVector* y, SimpleVector* y_linking, int child, int row, bool linking) const;
  virtual void axpyWithRowAtPosNeg( double alpha, StochVector* y_pos, SimpleVector* y_link_pos, StochVector* y_neg, SimpleVector* y_link_neg, int child, int row, bool linking ) const;

  virtual BorderedGenMatrix* raiseBorder( int m_conss, int n_vars );

  virtual StringGenMatrix* shaveLinkingConstraints( unsigned int n_conss );
  virtual void splitMatrix( const std::vector<int>& twolinks_start_in_block, const std::vector<unsigned int>& map_blocks_children, unsigned int n_links_in_root,
        const std::vector<MPI_Comm>& child_comms );

protected:

  virtual void writeToStreamDenseChild( std::ostream& out, int offset) const;
  virtual void writeToStreamDenseBorderedChild( const StringGenMatrix& border_left, std::ostream& out, int offset = 0 ) const;

  virtual void writeToStreamDenseRowLink( std::ostream& out, int rowidx) const;

  bool amatEmpty() const;
  virtual void shaveBorder(int m_conss, int n_vars, StringGenMatrix* border_left, StringGenMatrix* border_bottom);
  virtual StringGenMatrix* shaveLeftBorder( int n_vars );
  virtual StringGenMatrix* shaveLeftBorderChild( int n_vars );
};


/**
 * Dummy Class 
 */

class StochGenDummyMatrix : public StochGenMatrix {

protected:

public:

  StochGenDummyMatrix()
    : StochGenMatrix(0, 0, 0, 0, 0, 0, 0, 0, MPI_COMM_NULL) {};

  ~StochGenDummyMatrix() override = default;

  void AddChild( StochGenMatrix* ) override {};

 public:
  void updateTransposed() override {};

  void getSize( int& m, int& n ) const override { m = 0; n = 0; }
  void getSize( long long& m, long long& n ) const override { m = 0; n = 0; }

  GenMatrix* cloneEmptyRows( bool ) const override { return new StochGenDummyMatrix(); };
  GenMatrix* cloneFull( bool ) const  override { return new StochGenDummyMatrix(); };


  /** The actual number of structural non-zero elements in this sparse
   *  matrix. This includes so-called "accidental" zeros, elements that
   *  are treated as non-zero even though their value happens to be zero.
   */  
  int numberOfNonZeros() const override { return 0; };

  int isKindOf( int matType ) const override;

  void atPutDense( int, int, double*, int, int, int ) override {};
  void fromGetDense( int, int, double*, int, int, int ) override {};
  void columnScale( const OoqpVector& ) override {};
  void rowScale( const OoqpVector& ) override {};
  void symmetricScale( const OoqpVector& ) override {};
  void scalarMult( double ) override {};
  void fromGetSpRow( int, int, double[], int, int[], int&, int, int& ) override {};

  void atPutSubmatrix( int, int, DoubleMatrix&, int, int, int, int ) override {};

  void atPutSpRow( int, double[], int, int[], int& ) override {};

  void putSparseTriple( int[], int, int[], double[], int& ) override {};

  void getDiagonal( OoqpVector& ) override {};
  void setToDiagonal( const OoqpVector& ) override {};

  void mult ( double, OoqpVector&, double, const OoqpVector& ) const override {};
  void mult2 ( double, StochVector&, double, StochVector&, OoqpVector* ) override {};

  void transMult ( double, OoqpVector&, double, const OoqpVector& ) const override {};
  void transMult2 ( double, StochVector&, double, StochVector&, const OoqpVector* ) const override {};

  double abmaxnorm() const override { return 0.0; };
  double abminnormNonZero( double ) const override { return std::numeric_limits<double>::infinity(); };

  void permuteLinkingVarsChild( const std::vector<unsigned int>& )  override {};
  void getLinkVarsNnzChild( std::vector<int>& ) const override {};

  void getLinkVarsNnz( std::vector<int>& ) const override {};
  void writeToStream( std::ostream& ) const override {};

  void writeToStreamDense( std::ostream& ) const override {};
  void writeToStreamDense( std::ostream&, int ) const override{};
  void writeToStreamDenseBordered( const StringGenMatrix&, std::ostream&, int) const override{};

  void writeMPSformatRows( std::ostream&, int, OoqpVector* ) const override {};

 protected:
  void writeToStreamDenseChild( std::ostream&, int ) const override {};
  void writeToStreamDenseBorderedChild( const StringGenMatrix&, std::ostream&, int ) const override {};

  void writeToStreamDenseRowLink( std::ostream&, int ) const override {};

 public:
  void randomize( double, double, double* ) override {};

  void atPutDiagonal( int, OoqpVector& ) override {};
  void fromGetDiagonal( int, OoqpVector& ) override {};

  void initTransposedChild( bool ) override {};

  void columnScale2( const OoqpVector& ) override {};
  void rowScale2( const OoqpVector&, const OoqpVector* ) override {};

  void initTransposed( bool ) override {};
  void deleteTransposed() override {};

  void getNnzPerRow( OoqpVectorBase<int>&, OoqpVectorBase<int>* ) override {};
  void getNnzPerCol( OoqpVectorBase<int>&, OoqpVectorBase<int>* ) override {};
  void getNnzPerRow( OoqpVectorBase<int>& ) override {};
  void getNnzPerCol( OoqpVectorBase<int>& ) override {};

  void getRowMinMaxVec( bool, bool, const OoqpVector*, const OoqpVector*, OoqpVector&, OoqpVector* )override {};
  void getColMinMaxVec( bool, bool, const OoqpVector*, const OoqpVector*, OoqpVector&, OoqpVector* )override {};
  void getRowMinMaxVec( bool , bool, const OoqpVector*, OoqpVector& )override {};

  void getColMinMaxVec( bool, bool, const OoqpVector*, OoqpVector& )override {};

  void addRowSums( OoqpVector&, OoqpVector* ) const override {};
  void addColSums( OoqpVector&, OoqpVector* ) const override {};
  void addRowSums( OoqpVector& ) const override {};
  void addColSums( OoqpVector& ) const override {};

  void freeDynamicStorage() override {};
  void initStaticStorageFromDynamic( const OoqpVectorBase<int>&, const OoqpVectorBase<int>&) {};
  void initStaticStorageFromDynamic( const OoqpVectorBase<int>&, const OoqpVectorBase<int>&, const OoqpVectorBase<int>*, const OoqpVectorBase<int>* ) {};

  std::vector<int> get2LinkStartBlocks() const override { return std::vector<int>(); };

  void updateKLinkVarsCount( std::vector<int>& ) const override {};
  void updateKLinkConsCount( std::vector<int>& ) const override {};

  void permuteLinkingVars( const std::vector<unsigned int>& ) override {};
  void permuteLinkingCons( const std::vector<unsigned int>& ) override {};

  bool isRootNodeInSync() const override { return true; };

  int appendRow( const StochGenMatrix&, int, int, bool ) override { assert( 0 && "CANNOT APPEND ROW TO DUMMY MATRIX"); return -1; };
  double localRowTimesVec( const StochVector&, int, int, bool ) const override { assert( 0 && "CANNOT MULTIPLY ROW WITH DUMMY MATRIX"); return -1; };

  void axpyWithRowAt( double, StochVector*, SimpleVector*, int, int, bool ) const override {};
  void axpyWithRowAtPosNeg( double, StochVector*, SimpleVector*, StochVector*, SimpleVector*, int, int, bool ) const override {};

  BorderedGenMatrix* raiseBorder( int, int ) override { assert(0 && "CANNOT SHAVE BORDER OFF OF A DUMMY MATRIX"); return nullptr; };
  StringGenMatrix* shaveLinkingConstraints( unsigned int ) override { return new StringGenDummyMatrix(); };
  void splitMatrix( const std::vector<int>&, const std::vector<unsigned int>&, unsigned int, const std::vector<MPI_Comm>& ) override
     { assert(0 && "CANNOT SHAVE BORDER OFF OF A DUMMY MATRIX"); };

  void recomputeSize( StochGenMatrix* ) override {};
 protected:
  void shaveBorder( int, int, StringGenMatrix* border_left, StringGenMatrix* border_bottom) override
  { border_left->addChild(new StringGenDummyMatrix()); border_bottom->addChild(new StringGenDummyMatrix()); };
  StringGenMatrix* shaveLeftBorder( int ) override { return new StringGenDummyMatrix(); };
  StringGenMatrix* shaveLeftBorderChild( int ) override { return new StringGenDummyMatrix(); };

};


typedef SmartPointer<StochGenMatrix> StochGenMatrixHandle;

#endif
