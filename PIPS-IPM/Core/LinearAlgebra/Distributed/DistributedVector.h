#ifndef DISTRIBUTEDVECTOR_H
#define DISTRIBUTEDVECTOR_H

#include "Vector.hpp"
#include "DenseVector.hpp"
#include "mpi.h"
#include <vector>
#include <memory>

class DistributedTree;

class DistributedProblem;

enum class VectorType { PRIMAL, DUAL_Y, DUAL_Z };

template<typename T>
class DistributedVector : public Vector<T> {
public:
   DistributedVector(std::unique_ptr<Vector<T>> first, std::unique_ptr<Vector<T>> last, MPI_Comm mpi_comm);

   DistributedVector(int n, MPI_Comm mpiComm);
   DistributedVector(int n, int nl, MPI_Comm mpiComm);

   // this virtual AND override need to stay - Intel's icpc has a bug in its compiler triggering warnings if either is deleted...
   virtual ~DistributedVector() override = default;

   virtual void AddChild(std::shared_ptr<DistributedVector<T>> child);

   /** The data for this node. */
   std::shared_ptr<Vector<T>> first{};

   /** The linking constraint data for this node. */
   std::shared_ptr<Vector<T>> last{};

   /** Children of this node */
   std::vector<std::shared_ptr<DistributedVector<T>>> children;

   /** Links to this vectors parent.
    *  Needed when we multiply a matrix with this vector to get the appropriate linking first part.
    */
   DistributedVector<T>* parent{};

   /* MPI communicator */
   const MPI_Comm mpiComm{MPI_COMM_NULL};
   /* flag used to indicate if the children are distributed or not. */
   const bool iAmDistrib{false};
   const bool iAmSpecial{false};

   [[nodiscard]] Vector<T>* clone() const override;
   /* copy vector entries as well */
   [[nodiscard]] Vector<T>* clone_full() const override;

   void jointCopyFrom(const Vector<T>& vx, const Vector<T>& vy, const Vector<T>& vz) override;
   void jointCopyTo(Vector<T>& vx, Vector<T>& vy, Vector<T>& vz) const override;

   [[nodiscard]] bool isKindOf(int kind) const override;
   void setToZero() override;
   void setToConstant(T c) override;
   [[nodiscard]] bool isZero() const override;

   void copyFrom(const Vector<T>& v) override;
   void copyFromAbs(const Vector<T>& v) override;
   [[nodiscard]] double two_norm() const override;
   [[nodiscard]] T inf_norm() const override;
   [[nodiscard]] T one_norm() const override;
   void min(T& m, int& index) const override;
   void max(T& m, int& index) const override;
   void absminVecUpdate(Vector<T>& absminvec) const override;
   void absmaxVecUpdate(Vector<T>& absmaxvec) const override;
   void absmin(T& m) const override;
   void absminNonZero(T& m, T zero_eps) const override;
   [[nodiscard]] T fraction_to_boundary(const Vector<T>& v, T fraction) const override;
   [[nodiscard]] T find_blocking(const Vector<T>& wstep_vec, const Vector<T>& u_vec, const Vector<T>& ustep_vec, T maxStep, T* w_elt, T* wstep_elt, T* u_elt,
         T* ustep_elt, int& first_or_second) const override;
   void find_blocking_pd(const Vector<T>& wstep_vec, const Vector<T>& u_vec, const Vector<T>& ustep_vec, T& maxStepPri, T& maxStepDual, T& w_elt_p,
         T& wstep_elt_p, T& u_elt_p, T& ustep_elt_p, T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d, bool& primalBlocking,
         bool& dualBlocking) const override;

   void componentMult(const Vector<T>& v) override;
   void componentDiv(const Vector<T>& v) override;
   [[nodiscard]] bool componentEqual(const Vector<T>& v, T tol) const override;
   [[nodiscard]] bool componentNotEqual(const T val, T const tol) const override;

   void setNotIndicatedEntriesToVal(const T val, const Vector<T>& ind) override;

   void scalarMult(T num) override;
   void write_to_stream(std::ostream& out, int offset = 0) const override;

   void scale(T alpha) override;

   /** this += alpha * x */
   void add(T alpha, const Vector<T>& x) override;
   /** this += alpha * x * z */
   void add_product(T alpha, const Vector<T>& x, const Vector<T>& z) override;
   /** this += alpha * x / z */
   void add_quotient(T alpha, const Vector<T>& x, const Vector<T>& z) override;

   void add_constant(T c) override;
   void gondzioProjection(T rmin, T rmax) override;
   [[nodiscard]] T dotProductWith(const Vector<T>& v) const override;
   [[nodiscard]] T dotProductSelf(T scaleFactor) const override;
   [[nodiscard]] T scaled_dot_product_self(const Vector<T>& scale) const override;

   /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
    */
   [[nodiscard]] T shiftedDotProductWith(T alpha, const Vector<T>& mystep, const Vector<T>& yvec, T beta, const Vector<T>& ystep) const override;
   void negate() override;
   void invert() override;
   void safe_invert(T zero_replacement_value = 0.0) override;
   void sqrt() override;
   void roundToPow2() override;

   [[nodiscard]] bool all_positive() const override;

   void transform(const std::function<T(const T&)>& transformation) override;
    void transform_value(const std::function<T(const T&, const T&, const T&, const T&)>& transformation, const Vector<T>& lower_bounds_, const
    Vector<T>& upper_bounds_, const Vector<T>& integrality_) override;
    void fix_values(const Vector<T>& integrality_, double value) override;

   [[nodiscard]] T sum_reduce(const std::function<T(const T& a, const T& b)>& reduce) const override;
   [[nodiscard]] bool all_of(const std::function<bool(const T&)>& pred) const override;

   [[nodiscard]] bool matchesNonZeroPattern(const Vector<T>& select) const override;
   void selectNonZeros(const Vector<T>& select) override;
   void selectPositive() override;
   void selectNegative() override;
   [[nodiscard]] long long number_nonzeros() const override;
   void add_constant(T c, const Vector<T>& select) override;
   void add_quotient(T alpha, const Vector<T>& x, const Vector<T>& z, const Vector<T>& select) override;

   [[nodiscard]] bool are_positive(const Vector<T>& select) const override;
   void divideSome(const Vector<T>& div, const Vector<T>& select) override;
   void copyIntoArray(T[]) const override { assert("Not supported" && 0); };
   void copyFromArray(const T[]) override { assert("Not supported" && 0); };
   void copyFromArray(const char[]) override { assert("Not supported" && 0); };
   virtual void permuteVec0Entries(const std::vector<unsigned int>& permvec);
   virtual void permuteLinkingEntries(const std::vector<unsigned int>& permvec);
   [[nodiscard]] virtual std::vector<T> gatherStochVector() const;

   /** remove entries i for which select[i] == 0 */
   void removeEntries(const Vector<int>& select) override;

   [[nodiscard]] virtual int getSize() const { return this->n; };
   [[nodiscard]] int getNnzs() const override;

   virtual void recomputeSize();

   [[nodiscard]] virtual bool isRootNodeInSync() const;

   virtual void split(const std::vector<unsigned int>& map_blocks_children, const std::vector<MPI_Comm>& child_comms,
         const std::vector<int>& twolinks_start_in_block = std::vector<int>(), int n_links_in_root = -1);
   void move_first_to_parent();
   virtual DistributedVector<T>* raiseBorder(int n_first_to_shave, int n_last_to_shave);
   virtual void collapseFromHierarchical(const DistributedProblem& data_hier, const DistributedTree& tree_hier, VectorType type, bool empty_vec = false);
   virtual void appendHierarchicalToThis(DenseVector<T>* new_vec, DenseVector<T>* new_vecl, std::vector<std::shared_ptr<DistributedVector<T>>>& new_children,
         const DistributedTree& tree_hier, const DistributedProblem& data_hier, VectorType type, bool empty_vec);

   [[nodiscard]] virtual Vector<T>* getLinkingVecNotHierarchicalTop() const;

   void pushAwayFromZero(double tol, double amount, const Vector<T>* select) override;
   void getSumCountIfSmall(double tol, double& sum_small, int& n_close, const Vector<T>* select) const override;

   [[nodiscard]] double barrier_directional_derivative(const Vector<T>& x, const Vector<T>& bound, const Vector<T>& bound_indicator) const override;
   [[nodiscard]] double barrier_directional_derivative(const Vector<T>& x, double bound, const Vector<T>& bound_indicator) const override;

   [[nodiscard]] std::tuple<double, double, double, double> find_abs_nonzero_max_min_pair_a_by_b_plus_c_by_d(const Vector<T>& a,
      const Vector<T>& b, const Vector<T>& select_ab, bool use_ab, const Vector<T>& c, const Vector<T>& d, const Vector<T>& select_cd, bool use_cd, bool find_min) const override;

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

   void AddChild(std::shared_ptr<DistributedVector<T>>) override {};

   [[nodiscard]] DistributedVector<T>* clone() const override { return new DistributedDummyVector<T>(); }
   [[nodiscard]] DistributedVector<T>* clone_full() const override { return new DistributedDummyVector<T>(); }

   void jointCopyFrom(const Vector<T>&, const Vector<T>&, const Vector<T>&) override {};
   void jointCopyTo(Vector<T>&, Vector<T>&, Vector<T>&) const override {};

   [[nodiscard]] bool isKindOf(int kind) const override { return kind == kStochDummy || kind == kStochVector; }
   [[nodiscard]] bool isZero() const override { return true; };
   void setToZero() override {};
   void setToConstant(T) override {};
   void copyFrom(const Vector<T>&) override {};
   void copyFromAbs(const Vector<T>&) override {};
   [[nodiscard]] double two_norm() const override { return 0.0; }
   [[nodiscard]] T inf_norm() const override { return 0.0; }
   [[nodiscard]] T one_norm() const override { return 0.0; }
   void min(T&, int&) const override {};
   void max(T&, int&) const override {};
   void absminVecUpdate(Vector<T>&) const override {};
   void absmaxVecUpdate(Vector<T>&) const override {};
   void absmin(T& m) const override { m = std::numeric_limits<T>::infinity(); };
   void absminNonZero(T& m, T) const override { m = std::numeric_limits<T>::infinity(); };
   T fraction_to_boundary(const Vector<T>&, T fraction) const override { return T{1}; }
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
   void write_to_stream(std::ostream&, int) const override {};
   void scale(T) override {};

   /** this += alpha * x */
   void add(T, const Vector<T>&) override {};
   /** this += alpha * x * z */
   void add_product(T, const Vector<T>&, const Vector<T>&) override {};
   /** this += alpha * x / z */
   void add_quotient(T, const Vector<T>&, const Vector<T>&) override {};

   void add_constant(T) override {};
   void gondzioProjection(T, T) override {};
   T dotProductWith(const Vector<T>&) const override { return 0.0; }
   T dotProductSelf(T) const override { return 0.0; };
   [[nodiscard]] T scaled_dot_product_self(const Vector<T>&) const override { return 0.0; };

   /** Return the inner product <this + alpha * mystep, yvec + beta * ystep > */
   T shiftedDotProductWith(T, const Vector<T>&, const Vector<T>&, T, const Vector<T>&) const override { return 0.0; }
   void negate() override {};
   void invert() override {};
   void safe_invert(T) override {};
   void sqrt() override {};
   void roundToPow2() override {};
   [[nodiscard]] bool all_positive() const override { return true; };

   void transform(const std::function<T(const T&)>& transformation) override {};
   [[nodiscard]] T sum_reduce(const std::function<T(const T& a, const T& b)>& reduce) const override { return T{}; };
   [[nodiscard]] bool all_of(const std::function<bool(const T&)>&) const override { return true; };

   [[nodiscard]] bool matchesNonZeroPattern(const Vector<T>&) const override { return true; }
   void selectNonZeros(const Vector<T>&) override {};
   void selectPositive() override {};
   void selectNegative() override {};
   [[nodiscard]] long long number_nonzeros() const override { return 0; }
   void add_constant(T, const Vector<T>&) override {};
   void add_quotient(T, const Vector<T>&, const Vector<T>&, const Vector<T>&) override {};

   [[nodiscard]] bool are_positive(const Vector<T>&) const override { return 1; }
   void divideSome(const Vector<T>&, const Vector<T>&) override {};
   void copyIntoArray(T[]) const override {};
   void copyFromArray(const T[]) override {};
   void copyFromArray(const char[]) override {};

   void removeEntries(const Vector<int>&) override {};
   void permuteVec0Entries(const std::vector<unsigned int>&) override {};
   void permuteLinkingEntries(const std::vector<unsigned int>&) override {};
   [[nodiscard]] std::vector<T> gatherStochVector() const override { return std::vector<T>(0); };

   [[nodiscard]] int getSize() const override { return 0; };
   [[nodiscard]] int getNnzs() const override { return 0; };

   void recomputeSize() override {};

   [[nodiscard]] bool isRootNodeInSync() const override { return true; };

   void split(const std::vector<unsigned int>&, const std::vector<MPI_Comm>&, const std::vector<int>&, int) override {};
   DistributedVector<T>* raiseBorder(int, int) override {
      assert(0 && "This should never be attempted");
      return nullptr;
   };

   void appendHierarchicalToThis(DenseVector<T>* new_vec, DenseVector<T>* new_vecl, std::vector<std::shared_ptr<DistributedVector<T>>>& new_children,
         const DistributedTree& tree_hier, const DistributedProblem& data_hier, VectorType type, bool empty_vec) override;

   Vector<T>* getLinkingVecNotHierarchicalTop() const override {
      assert(false && "Should not end up here");
      return nullptr;
   };

   void pushAwayFromZero(double, double, const Vector<T>*) override {};
   void getSumCountIfSmall(double, double&, int&, const Vector<T>*) const override {};

   void pushSmallComplementarityPairs(Vector<T>&, const Vector<T>&, double, double, double) override {};

   [[nodiscard]] double barrier_directional_derivative(const Vector<T>& x, const Vector<T>& bound, const Vector<T>& bound_indicator) const override { return 0.;};
   [[nodiscard]] double barrier_directional_derivative(const Vector<T>& x, double bound, const Vector<T>& bound_indicator) const override { return 0.; };

   [[nodiscard]] std::tuple<double, double, double, double> find_abs_nonzero_max_min_pair_a_by_b_plus_c_by_d(const Vector<T>& a,
      const Vector<T>& b, const Vector<T>& select_ab, bool use_ab, const Vector<T>& c, const Vector<T>& d, const Vector<T>& select_cd, bool use_cd, bool find_min) const override {
      if (find_min) {
         return {std::numeric_limits<double>::infinity(), 1.0, std::numeric_limits<double>::infinity(), 1.0};
      } else {
         return {0.0, 1.0, 0.0, 1.0};
      }
   }
};

#endif
