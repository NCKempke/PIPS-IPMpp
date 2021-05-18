#ifndef DISTRIBUTEDVECTOR_H
#define DISTRIBUTEDVECTOR_H

#include "Vector.hpp"
#include "SmartPointer.h"
#include "SimpleVector.h"
#include "mpi.h"
#include <vector>
#include <memory>

class DistributedTree;

class DistributedQP;

enum class VectorType { PRIMAL, DUAL_Y, DUAL_Z };

template<typename T>
class DistributedVector : public Vector<T> {
public:
   DistributedVector(Vector<T>* first, Vector<T>* last, MPI_Comm mpi_comm);

   DistributedVector(int n, MPI_Comm mpiComm);
   DistributedVector(int n, int nl, MPI_Comm mpiComm);

   // this virtual AND override need to stay - intel has a bug in its compiler triggering warnings if either is deleted...
   virtual ~DistributedVector() override;

   virtual void AddChild(DistributedVector<T>* child);

   // TODO : use unique pointers
   /** The data for this node. */
   Vector<T>* first{};

   /** The linking constraint data for this node. */
   Vector<T>* last{};

   /** Children of this node */
   std::vector<DistributedVector<T>*> children;

   /** Links to this vectors parent.
    *  Needed when we multiply a matrix with this vector to get the appropriate linking first part.
    */
   DistributedVector<T>* parent{};

   /* MPI communicator */
   const MPI_Comm mpiComm{MPI_COMM_NULL};
   /* flag used to indicate if the children are distributed or not. */
   const bool iAmDistrib{false};
   const bool iAmSpecial{false};

   Vector<T>* clone() const override;
   /* copy vector entries as well */
   Vector<T>* cloneFull() const override;

   void jointCopyFrom(const Vector<T>& vx, const Vector<T>& vy, const Vector<T>& vz) override;
   void jointCopyTo(Vector<T>& vx, Vector<T>& vy, Vector<T>& vz) const override;

   [[nodiscard]] bool isKindOf(int kind) const override;
   void setToZero() override;
   void setToConstant(T c) override;
   [[nodiscard]] bool isZero() const override;

   void copyFrom(const Vector<T>& v) override;
   void copyFromAbs(const Vector<T>& v) override;
   [[nodiscard]] double two_norm() const override;
   T inf_norm() const override;
   T one_norm() const override;
   void min(T& m, int& index) const override;
   void max(T& m, int& index) const override;
   void absminVecUpdate(Vector<T>& absminvec) const override;
   void absmaxVecUpdate(Vector<T>& absmaxvec) const override;
   void absmin(T& m) const override;
   void absminNonZero(T& m, T zero_eps) const override;
   T stepbound(const Vector<T>& v, T maxStep) const override;
   T find_blocking(const Vector<T>& wstep_vec, const Vector<T>& u_vec, const Vector<T>& ustep_vec, T maxStep, T* w_elt, T* wstep_elt, T* u_elt,
         T* ustep_elt, int& first_or_second) const override;
   void find_blocking_pd(const Vector<T>& wstep_vec, const Vector<T>& u_vec, const Vector<T>& ustep_vec, T& maxStepPri, T& maxStepDual, T& w_elt_p,
         T& wstep_elt_p, T& u_elt_p, T& ustep_elt_p, T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d, bool& primalBlocking,
         bool& dualBlocking) const override;

   void componentMult(const Vector<T>& v) override;
   void componentDiv(const Vector<T>& v) override;
   bool componentEqual(const Vector<T>& v, T tol) const override;
   bool componentNotEqual(const T val, T const tol) const override;

   void setNotIndicatedEntriesToVal(const T val, const Vector<T>& ind) override;

   void scalarMult(T num) override;
   void writeToStream(std::ostream& out, int offset = 0) const override;
   void writefToStream(std::ostream& out, const char format[]) const override;

   void scale(T alpha) override;

   /** this += alpha * x */
   void axpy(T alpha, const Vector<T>& x) override;
   /** this += alpha * x * z */
   void axzpy(T alpha, const Vector<T>& x, const Vector<T>& z) override;
   /** this += alpha * x / z */
   void axdzpy(T alpha, const Vector<T>& x, const Vector<T>& z) override;

   void addConstant(T c) override;
   void gondzioProjection(T rmin, T rmax) override;
   T dotProductWith(const Vector<T>& v) const override;
   T dotProductSelf(T scaleFactor) const override;

   /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
    */
   T shiftedDotProductWith(T alpha, const Vector<T>& mystep, const Vector<T>& yvec, T beta, const Vector<T>& ystep) const override;
   void negate() override;
   void invert() override;
   void invertSave(T zeroReplacementVal = 0.0) override;
   void applySqrt() override;
   void roundToPow2() override;

   [[nodiscard]] bool allPositive() const override;
   bool allOf(const std::function<bool(const T&)>& pred) const override;

   bool matchesNonZeroPattern(const Vector<T>& select) const override;
   void selectNonZeros(const Vector<T>& select) override;
   void selectPositive() override;
   void selectNegative() override;
   [[nodiscard]] long long numberOfNonzeros() const override;
   void add_constant(T c, const Vector<T>& select) override;
   void writefSomeToStream(std::ostream&, const char[], const Vector<T>&) const override { assert(false && "Not yet implemented"); };
   void axdzpy(T alpha, const Vector<T>& x, const Vector<T>& z, const Vector<T>& select) override;

   bool somePositive(const Vector<T>& select) const override;
   void divideSome(const Vector<T>& div, const Vector<T>& select) override;
   void copyIntoArray(T[]) const override { assert("Not supported" && 0); };
   void copyFromArray(const T[]) override { assert("Not supported" && 0); };
   void copyFromArray(const char[]) override { assert("Not supported" && 0); };
   virtual void permuteVec0Entries(const std::vector<unsigned int>& permvec);
   virtual void permuteLinkingEntries(const std::vector<unsigned int>& permvec);
   virtual std::vector<T> gatherStochVector() const;

   /** remove entries i for which select[i] == 0 */
   void removeEntries(const Vector<int>& select) override;

   [[nodiscard]] virtual int getSize() const { return this->n; };
   [[nodiscard]] int getNnzs() const override;

   virtual void recomputeSize();

   [[nodiscard]] virtual bool isRootNodeInSync() const;

   virtual void split(const std::vector<unsigned int>& map_blocks_children, const std::vector<MPI_Comm>& child_comms,
         const std::vector<int>& twolinks_start_in_block = std::vector<int>(), int n_links_in_root = -1);
   virtual DistributedVector<T>* raiseBorder(int n_first_to_shave, int n_last_to_shave);
   virtual void collapseFromHierarchical(const DistributedQP& data_hier, const DistributedTree& tree_hier, VectorType type, bool empty_vec = false);
   virtual void appendHierarchicalToThis(SimpleVector<T>* new_vec, SimpleVector<T>* new_vecl, std::vector<DistributedVector<T>*>& new_children,
         const DistributedTree& tree_hier, const DistributedQP& data_hier, VectorType type, bool empty_vec);

   virtual Vector<T>* getLinkingVecNotHierarchicalTop() const;

   void pushAwayFromZero(double tol, double amount, const Vector<T>* select) override;
   void getSumCountIfSmall(double tol, double& sum_small, int& n_close, const Vector<T>* select) const override;
protected:
   DistributedVector() = default;

private:
   void
   pushSmallComplementarityPairs(Vector<T>& other_vec_in, const Vector<T>& select_in, double tol_this, double tol_other, double tol_pairs) override;
};

/** DUMMY VERSION
 *
 */
template<typename T>
class DistributedDummyVector : public DistributedVector<T> {
public:

   DistributedDummyVector() : DistributedVector<T>(0, MPI_COMM_NULL) {};
   ~DistributedDummyVector() override = default;

   void AddChild(DistributedVector<T>*) override {};

   DistributedVector<T>* clone() const override { return new DistributedDummyVector<T>(); }
   DistributedVector<T>* cloneFull() const override { return new DistributedDummyVector<T>(); }

   void jointCopyFrom(const Vector<T>&, const Vector<T>&, const Vector<T>&) override {};
   void jointCopyTo(Vector<T>&, Vector<T>&, Vector<T>&) const override {};

   bool isKindOf(int kind) const override { return kind == kStochDummy || kind == kStochVector; }
   bool isZero() const override { return true; };
   void setToZero() override {};
   void setToConstant(T) override {};
   void copyFrom(const Vector<T>&) override {};
   void copyFromAbs(const Vector<T>&) override {};
   double two_norm() const override { return 0.0; }
   T inf_norm() const override { return 0.0; }
   T one_norm() const override { return 0.0; }
   void min(T&, int&) const override {};
   void max(T&, int&) const override {};
   void absminVecUpdate(Vector<T>&) const override {};
   void absmaxVecUpdate(Vector<T>&) const override {};
   void absmin(T& m) const override { m = std::numeric_limits<T>::infinity(); };
   void absminNonZero(T& m, T) const override { m = std::numeric_limits<T>::infinity(); };
   T stepbound(const Vector<T>&, T maxStep) const override { return maxStep; }
   T find_blocking(const Vector<T>&, const Vector<T>&, const Vector<T>&, T maxStep, T*, T*, T*, T*, int&) const override { return maxStep; };
   void
   find_blocking_pd(const Vector<T>&, const Vector<T>&, const Vector<T>&, T&, T&, T&, T&, T&, T&, T&, T&, T&, T&, bool&, bool&) const override {};

   void componentMult(const Vector<T>&) override {};
   void componentDiv(const Vector<T>&) override {};
   bool componentEqual(const Vector<T>& v, T) const override {
      if (!v.isKindOf(kStochDummy))
         std::cout << "one should never end up here" << std::endl;
      return v.isKindOf(kStochDummy);
   };
   bool componentNotEqual(const T, const T) const override { return true; };
   void setNotIndicatedEntriesToVal(T, const Vector<T>&) override {};

   void scalarMult(T) override {};
   void writeToStream(std::ostream&, int) const override {};
   void writefToStream(std::ostream&, const char[]) const override {};
   void scale(T) override {};

   /** this += alpha * x */
   void axpy(T, const Vector<T>&) override {};
   /** this += alpha * x * z */
   void axzpy(T, const Vector<T>&, const Vector<T>&) override {};
   /** this += alpha * x / z */
   void axdzpy(T, const Vector<T>&, const Vector<T>&) override {};

   void addConstant(T) override {};
   void gondzioProjection(T, T) override {};
   T dotProductWith(const Vector<T>&) const override { return 0.0; }
   T dotProductSelf(T) const override { return 0.0; };

   /** Return the inner product <this + alpha * mystep, yvec + beta * ystep > */
   T shiftedDotProductWith(T, const Vector<T>&, const Vector<T>&, T, const Vector<T>&) const override { return 0.0; }
   void negate() override {};
   void invert() override {};
   void invertSave(T) override {};
   void applySqrt() override {};
   void roundToPow2() override {};
   bool allPositive() const override { return true; };
   bool allOf(const std::function<bool(const T&)>&) const override { return true; };

   bool matchesNonZeroPattern(const Vector<T>&) const override { return true; }
   void selectNonZeros(const Vector<T>&) override {};
   void selectPositive() override {};
   void selectNegative() override {};
   long long numberOfNonzeros() const override { return 0; }
   void add_constant(T, const Vector<T>&) override {};
   void writefSomeToStream(std::ostream&, const char[], const Vector<T>&) const override {};
   void axdzpy(T, const Vector<T>&, const Vector<T>&, const Vector<T>&) override {};

   bool somePositive(const Vector<T>&) const override { return 1; }
   void divideSome(const Vector<T>&, const Vector<T>&) override {};
   void copyIntoArray(T[]) const override {};
   void copyFromArray(const T[]) override {};
   void copyFromArray(const char[]) override {};

   void removeEntries(const Vector<int>&) override {};
   void permuteVec0Entries(const std::vector<unsigned int>&) override {};
   void permuteLinkingEntries(const std::vector<unsigned int>&) override {};
   std::vector<T> gatherStochVector() const override { return std::vector<T>(0); };

   int getSize() const override { return 0; };
   int getNnzs() const override { return 0; };

   void recomputeSize() override {};

   bool isRootNodeInSync() const override { return true; };

   void split(const std::vector<unsigned int>&, const std::vector<MPI_Comm>&, const std::vector<int>&, int) override {};
   DistributedVector<T>* raiseBorder(int, int) override {
      assert(0 && "This should never be attempted");
      return nullptr;
   };

   void appendHierarchicalToThis(SimpleVector<T>* new_vec, SimpleVector<T>* new_vecl, std::vector<DistributedVector<T>*>& new_children,
         const DistributedTree& tree_hier, const DistributedQP& data_hier, VectorType type, bool empty_vec) override;

   Vector<T>* getLinkingVecNotHierarchicalTop() const override {
      assert(false && "Should not end up here");
      return nullptr;
   };

   void pushAwayFromZero(double, double, const Vector<T>*) override {};
   void getSumCountIfSmall(double, double&, int&, const Vector<T>*) const override {};

   void pushSmallComplementarityPairs(Vector<T>&, const Vector<T>&, double, double, double) override {};
};

#endif
