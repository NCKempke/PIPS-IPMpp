#ifndef STOCHSYMMATRIX_H
#define STOCHSYMMATRIX_H

#include "AbstractMatrix.h"
#include "SparseSymmetricMatrix.h"
#include "SparseMatrix.h"
#include "StripMatrix.h"
#include "pipsport.h"

#include <vector>
#include <iostream>
#include <fstream>

#include "mpi.h"

class BorderedSymmetricMatrix;


/*
 * [ Q0  R1^T ... RN^T ]
 * [ R1  Q1            ]
 * [ R2     Q2         ]
 * [ .         .       ]
 * [ .           .     ]
 * [ RN            QN  ]
 */

class DistributedSymmetricMatrix : public SymmetricMatrix {

private:

   // note: also used for dummy class!
   virtual void deleteEmptyRowsCols(const Vector<int>& nnzVec, const Vector<int>* linkParent);
   virtual void writeToStreamDenseChild(std::stringstream& out, int offset) const;

public:
   DistributedSymmetricMatrix(SymmetricMatrix* diag, SparseMatrix* border, MPI_Comm mpiComm);

   /** Constructs a matrix with local size 'local_n' having 'local_nnz' local nonzeros
       and set the global size and the id to to 'global_n' and 'id', respectively.
       The parameter 'id' is used for output/debug purposes only.
       The created matrix will have no children.*/
   DistributedSymmetricMatrix(long long global_n, int local_n, int local_nnz, MPI_Comm mpiComm);

   ~DistributedSymmetricMatrix() override;

   std::vector<DistributedSymmetricMatrix*> children;
   SymmetricMatrix* diag{};
   SparseMatrix* border{};

   long long n{0};
   MPI_Comm mpiComm{MPI_COMM_NULL};
   int iAmDistrib{0};

   void AddChild(DistributedSymmetricMatrix* child);

   [[nodiscard]] SymmetricMatrix* clone() const override;

   [[nodiscard]] int is_a(int type) const override;

   void fromGetDense(int, int, double*, int, int, int) const override { assert(false && "Not implemented"); };
   void symAtPutSpRow(int, const double[], int, const int[], int&) override { assert(false && "Not implemented"); };

   [[nodiscard]] std::pair<long long, long long> n_rows_columns() const override;
   [[nodiscard]] long long n_rows() const override;
   [[nodiscard]] long long n_columns() const override;
   [[nodiscard]] long long size() const override;

   void symAtPutSubmatrix(int, int, const AbstractMatrix&, int, int, int, int) override { assert(false && "Not implemented"); };;
   void fromGetSpRow(int, int, double[], int, int[], int&, int, int&) const override { assert(false && "Not implemented"); };

   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;
   void transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   [[nodiscard]] double inf_norm() const override;
   using AbstractMatrix::abminnormNonZero;
   [[nodiscard]] double abminnormNonZero(double tol) const override;

   void writeToStream(std::ostream&) const override { assert(false && "Not implemented"); };

   void writeToStreamDense(std::ostream& out) const override;

   void getDiagonal(Vector<double>& vec) const override;
   void setToDiagonal(const Vector<double>& vec) override;
   void atPutDiagonal(int idiag, const Vector<double>& v) override;
   void atAddDiagonal(int idiag, const Vector<double>& v) override;
   void fromGetDiagonal(int, Vector<double>& x) const override;

   void putSparseTriple(const int[], int, const int[], const double[], int&) override { assert(false && "Not implemented"); };

   void symmetricScale(const Vector<double>& vec) override;
   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;

   void scalarMult(double num) override;

   // note: also used for dummy class!
   virtual void deleteEmptyRowsCols(const Vector<int>& nnzVec) {
      deleteEmptyRowsCols(nnzVec, nullptr);
   }

   // TODO specify border bottom and left..
   virtual BorderedSymmetricMatrix* raiseBorder(int n_vars);
   virtual void splitMatrix(const std::vector<unsigned int>& map_blocks_children, const std::vector<MPI_Comm>& child_comms);

   void recomputeSize();
protected:
   virtual StripMatrix* shaveBorder(int n_vars);
   virtual StripMatrix* shaveBorder2(int n_vars);

   DistributedSymmetricMatrix* parent{};
};

/** 
 * Dummy stochastic symmetric matrix
 */

class StochSymDummyMatrix : public DistributedSymmetricMatrix {

private:
   void writeToStreamDenseChild(std::stringstream&, int) const override {};

public:

   StochSymDummyMatrix() : DistributedSymmetricMatrix(0, 0, 0, MPI_COMM_NULL) {};

   ~StochSymDummyMatrix() override = default;

   [[nodiscard]] SymmetricMatrix* clone() const override { return new StochSymDummyMatrix(); };

   [[nodiscard]] int is_a(int type) const override;

   [[nodiscard]] std::pair<long long, long long> n_rows_columns() const override { return {0,0}; };

   [[nodiscard]] long long n_rows() const override { return size(); };
   [[nodiscard]] long long n_columns() const override { return size(); };
   [[nodiscard]] long long size() const override { return 0; };

   void mult(double, Vector<double>&, double, const Vector<double>&) const override {};
   void transMult(double, Vector<double>&, double, const Vector<double>&) const override {};

   [[nodiscard]] double inf_norm() const override { return 0.0; }
   [[nodiscard]] double abminnormNonZero(double) const override { return std::numeric_limits<double>::infinity(); }

   void writeToStreamDense(std::ostream&) const override {};

   void getDiagonal(Vector<double>&) const override {};
   void setToDiagonal(const Vector<double>&) override {};
   void atPutDiagonal(int, const Vector<double>&) override {};
   void atAddDiagonal(int, const Vector<double>&) override {};
   void fromGetDiagonal(int, Vector<double>&) const override {};

   void symmetricScale(const Vector<double>&) override {};
   void columnScale(const Vector<double>&) override {};
   void rowScale(const Vector<double>&) override {};
   void scalarMult(double) override {};

   BorderedSymmetricMatrix* raiseBorder(int) override {
      assert(0 && "CANNOT SHAVE BORDER OFF OF A DUMMY MATRIX");
      return nullptr;
   };
   void splitMatrix(const std::vector<unsigned int>&, const std::vector<MPI_Comm>&) override {};

protected:
   StripMatrix* shaveBorder(int) override { return new StringGenDummyMatrix(); };
   StripMatrix* shaveBorder2(int) override { return new StringGenDummyMatrix(); };

};
#endif
