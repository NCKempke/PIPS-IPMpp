/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SPARSESTORAGE_H
#define SPARSESTORAGE_H

#include "DoubleMatrix.h"
#include "SparseStorageHandle.h"
#include "OoqpVectorHandle.h"
#include "pipsport.h"

#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

/** A class for managing the matrix elements used by sparse matrices.
 *  @ingroup SparseLinearAlgebra
 */
class SparseStorage : public DoubleStorage {
private:
  /** store absolute non-zero minimum entry of row i and first[i] in first[i]; empty rows get value 0.0  */
  void getRowMinVec(const double* colScaleVec, double* vec) const;

  /** store absolute non-zero maximum entry of row i and first[i] in first[i]; empty rows get value 0.0  */
  void getRowMaxVec(const double* colScaleVec, double* vec) const;

  class index_sort
  {
     private:
       const int* indices;
       const int maxsize;

     public:
       index_sort(const int* indices, int maxsize) : indices(indices), maxsize(maxsize) {}

       bool operator()(int i, int j) const
       {
          assert(i < maxsize && j < maxsize);
          return (indices[i] < indices[j]);
       }
  };
protected:
  int neverDeleteElts;
  
public:
  static int instances;

  int m{};
  int n{};
  int len{};
  int* jcolM{};
  int* krowM{};
  double* M{};

  SparseStorage( int m_, int n_, int len_ );
  SparseStorage( int m_, int n_, int len_,
		 int * krowM_, int * jcolM_, double * M_,
		 int deleteElts=0);

  void copyFrom(int * krowM_, int * jcolM_, double * M_) const;

  void shiftRows( int row, int shift, int& info );
  void getSize( int& m, int& n ) const override;
  int rows() const { return m; }
  int cols() const { return n; }

  bool isValid() const;
  bool isSorted() const;


  int length() { return len; };
  int numberOfNonZeros() const {	return krowM[m]; };
  void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent ) override;
  void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent ) override;

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info );

  void getDiagonal( OoqpVector& vec ) override;
  void setToDiagonal( const OoqpVector& vec ) override;

  void columnScale( const OoqpVector& vec ) override;
  void rowScale( const OoqpVector& vec ) override;
  void symmetricScale( const OoqpVector& vec ) override;
  void scalarMult( double num ) override;

  void atPutSpRow( int col, double A[], int lenA, int irowA[],
			   int& info ) override;

  void fromGetSpRow( int row, int col,
			     double A[], int lenA, int irowA[], int& nnz,
			     int rowExtent, int& info ) override;

  virtual void randomize( double alpha, double beta, double * seed );

  virtual void clear();

  virtual void getTransposePat( int row, int col, int rowExtent, int colExtent,
				int kpat[], int krowM[], int jcolM[] );
  virtual void getFromPat( double data[], int n, int kpat[] );
  virtual void mult( double beta,  double y[], int incy,
		     double alpha, const double x[], int incx ) const;
  virtual void multSym( double beta,  double y[], int incy,
             double alpha, const double x[], int incx ) const;

  virtual void transMult ( double beta,  double y[], int incy,
			   double alpha, const double x[], int incx ) const;
  virtual void transMultD( double beta, double y[], int incy, double alpha, const double x[], const double d[], int incxd ) const;

  void atPutDiagonal( int idiag, OoqpVector& v ) override;
  void fromGetDiagonal( int idiag, OoqpVector& v ) override;

  virtual void atPutDiagonal( int idiag,
			      double x[], int incx, int extent );

  virtual void writeToStream(std::ostream& out) const;
  void writeNNZpatternToStreamDense( std::ostream& out ) const;
  virtual void writeToStreamDense(std::ostream& out) const;
  virtual void writeToStreamDenseRow( std::ostream& out, int rowidx) const;

  virtual void symmetrize( int& info);
  double abmaxnorm() const override;
  double abminnormNonZero( double tol = 1e-30 ) const override;

  /** Computes the sparsity pattern of MtM = M^T * D * M 
   *  where D=diag(d) is a diagonal matrix and M=this.
   *  
   *  Find the nonzero pattern of the product matrix, allocate it and returns it.
   *
   *  Also allocates, builds and returns this^T since it is needed later for
   *   numerical multiplication.
   */
  void matTransDSymbMultMat(double* /*d*/,
			    int* krowMt, int* jcolMt, double* /*dMt*/,
			    int** krowMtM, int** jcolMtM, double** dMtM); 
			    

  /** Numerical multiplication MtM = M^T * D * M  where 
   *  D=diag(d) is a diagonal matrix and M=this.
   *  M^T and MtM buffers should be allocated before calling this method by calling
   *  method matTransDSymbMultMat.
   */
  void matTransDMultMat(double* d, 
			int* krowMt, int* jcolMt, double* dMt,
			int* krowMtM, int* jcolMtM, double* dMtM);
  void matTransDinvMultMat(double* d, 
			int* krowMt, int* jcolMt, double* dMt,
			int* krowMtM, int* jcolMtM, double* dMtM);
  /** Builds the transpose: Mt = this^T */
  void transpose(int* krowMt, int* jcolMt, double* dMt) const;

  void reduceToLower();

  void multMatSymUpper( double beta, SparseStorage& y,
                double alpha, const double x[], int yrow, int ycolstart ) const;

  void transMultLower( double beta,  double y[],
			       double alpha, double x[], int firstrow );
  void transMultMat( double beta,  double* Y, int ny, int ldy,
						 double alpha, double *X, int ldx);
  void transMultMatLower( double* Y, int ny, int firstrow,
						 double alpha, double *X, int ldx);

  /** Y <- alpha* M^T X + beta*Y, where M is this 
   * Special update function, computes only the elements in Y that are lower
   * triangular elements in a larger matrix that contains Y (see impl file for 
   * more details)
   */
  void transMultMatLower( double beta,  double* Y, int ny, int ldy,
			  double alpha, double *X, int ldx, int colStart);

  void fromGetColBlock(int col, double *A, int lda, int colExtent, int* colSparsity, bool &allzero);

  void fromGetRowsBlock( double* rows_array_dense, size_t row_start, size_t n_rows, size_t array_line_size, size_t array_line_offest, int* row_sparsity ) const;
  void fromGetRowsBlock(const int* rowIndices, int nRows, int arrayLineSize, int arrayLineOffset, double* rowsArrayDense, int* rowSparsity = nullptr) const;

  /** add nnz per row to given array (of size nRows) */
  void addNnzPerRow(int* vec) const { addNnzPerRow(vec, 0, m); };
  void addNnzPerRow(int* vec, int begin_rows, int end_rows) const;

  void getLinkVarsNnz(std::vector<int>& vec) const;

  /** add abs. sum per row to given array (of size nRows) */
  void addRowSums( double* vec ) const;

  /** store absolute non-zero minimum/maximum entry of row i and first[i] in first[i];
   *  empty rows get value 0.0 for maximization and <double>::max() for minimization  */
  void getRowMinMaxVec(bool getMin, const double* colScaleVec, double* vec) const;

  void permuteRows(const std::vector<unsigned int>& permvec);
  void permuteCols(const std::vector<unsigned int>& permvec);

  void sortCols();

  void dump(const std::string& filename);

  void deleteEmptyRowsCols(const int* nnzRowVec, const int* nnzColVec);

  void getSparseTriplet_c2fortran(int*& irn, int*& jcn, double*& val) const;

  void getSparseTriplet_fortran2fortran(int*& irn, int*& jcn, double*& val) const;

  void deleteEmptyRows(int*& orgIndex);

  // should be used with care! other methods might nor work correctly todo: add flag to check in other methods
  void c2fortran();

  void fortran2c();

  bool fortranIndexed() const;

  void set2FortranIndexed();

  void deleteZeroRowsColsSym(int*& new2orgIdx);


  /*
   * computes the full sparse matrix representation from a upper triangular symmetric sparse representation
   *
   * Must be square, the storage for the full representation will be allocated within the matrix and must be released later
   */
  void fullMatrixFromUpperTriangular(int*& rowPtrFull, int*& colIdxFull, double*& valuesFull) const;

  virtual SparseStorage* shaveLeft( int n_cols );
  virtual SparseStorage* shaveSymLeftBottom( int n );
  virtual SparseStorage* shaveBottom( int n_rows );
  virtual void dropNEmptyRowsBottom( int n_rows );
  virtual void dropNEmptyRowsTop( int n_rows );

  virtual ~SparseStorage();

private:
  bool isFortranIndexed{false};
};

#endif
