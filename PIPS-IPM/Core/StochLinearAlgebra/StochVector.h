#ifndef STOCHVECTOR_H
#define STOCHVECTOR_H

#include "OoqpVector.h"
#include "SimpleVector.h"
#include "StochVector_fwd.h"
#include "StochVectorHandle.h"

#include "mpi.h"

#include <vector>

class sTree;
class sData;

enum class VectorType{ PRIMAL, DUAL_Y, DUAL_Z };

template <typename T>
class StochVectorBase : public OoqpVectorBase<T> {


public:
  StochVectorBase( OoqpVectorBase<T>* vec, OoqpVectorBase<T>* vecl, MPI_Comm mpi_comm);

  StochVectorBase( int n, MPI_Comm mpiComm );
  StochVectorBase( int n, int nl, MPI_Comm mpiComm );
  virtual ~StochVectorBase() override;

  virtual void AddChild(StochVectorBase<T>* child);

  // TODO : use unique pointers
  /** The data for this node. */
  OoqpVectorBase<T>* vec{};

  /** The linking constraint data for this node. */
  OoqpVectorBase<T>* vecl{};

  /** Children of this node */
  std::vector<StochVectorBase<T>*> children;

  /** Links to this vectors parent.
   *  Needed when we multiply a matrix with this vector to get the appropriate linking vec part.
   */
  StochVectorBase<T>* parent{};

public:
  /* MPI communicator */
  const MPI_Comm mpiComm{MPI_COMM_NULL};
  /* flag used to indicate if the children are distributed or not. */
  const bool iAmDistrib{false};
  const bool iAmSpecial{false};

  OoqpVectorBase<T>* clone() const override;
  /* copy vector entries as well */
  OoqpVectorBase<T>* cloneFull() const override;

  void jointCopyFrom(const OoqpVectorBase<T>& vx, const OoqpVectorBase<T>& vy, const OoqpVectorBase<T>& vz) override;
  void jointCopyTo(OoqpVectorBase<T>& vx, OoqpVectorBase<T>& vy, OoqpVectorBase<T>& vz) const override;

  bool isKindOf( int kind ) const override;
  void setToZero() override;
  void setToConstant( T c ) override;
  bool isZero() const override;

   void randomize( T, T, T* ) override { assert( "Not implemented" && 0 ); };
   void copyFrom( const OoqpVectorBase<T>& v ) override;
   void copyFromAbs(const OoqpVectorBase<T>& v) override;
   double twonorm() const override;
   T infnorm() const override;
   T onenorm() const override;
   void min( T& m, int& index ) const override;
   void max( T& m, int& index ) const override;
   void absminVecUpdate( OoqpVectorBase<T>& absminvec) const override;
   void absmaxVecUpdate( OoqpVectorBase<T>& absmaxvec) const override;
   void absmin( T& m ) const override;
   void absminNonZero(T& m, T zero_eps) const override;
   T stepbound(const OoqpVectorBase<T> & v, T maxStep ) const override;
   T findBlocking(const OoqpVectorBase<T> & wstep_vec,
			      const OoqpVectorBase<T> & u_vec,
			      const OoqpVectorBase<T> & ustep_vec,
			      T maxStep,
			      T *w_elt,
			      T *wstep_elt,
			      T *u_elt,
			      T *ustep_elt,
			      int& first_or_second) const override;
   void findBlocking_pd(const OoqpVectorBase<T>& wstep_vec,
               const OoqpVectorBase<T>& u_vec,
  			      const OoqpVectorBase<T>& ustep_vec,
  			      T& maxStepPri, T& maxStepDual,
  			      T& w_elt_p, T& wstep_elt_p, T& u_elt_p, T& ustep_elt_p,
  				   T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d,
				   bool& primalBlocking, bool& dualBlocking) const override;

   void componentMult( const OoqpVectorBase<T>& v ) override;
   void componentDiv ( const OoqpVectorBase<T>& v ) override;
   bool componentEqual( const OoqpVectorBase<T>& v , T tol) const override;
   bool componentNotEqual( const T val, T const tol ) const override;

   void setNotIndicatedEntriesToVal( const T val, const OoqpVectorBase<T>& ind ) override;

   void scalarMult( T num ) override;
   void writeToStream(std::ostream& out, int offset = 0 ) const override;
   void writefToStream( std::ostream& out, const char format[] ) const override;

   void writeMPSformatOnlyRhs( std::ostream&, const std::string, const OoqpVectorBase<T>*) const override {};
   void writeMPSformatBoundsWithVar( std::ostream&, const std::string, const OoqpVectorBase<T>*, bool ) const override {};
   void writeMPSformatRhs( std::ostream& out, int rowType, const OoqpVectorBase<T>* irhs) const override;
   void writeMPSformatBounds( std::ostream& out, const OoqpVectorBase<T>* ix, bool upperBound) const override;

   void scale( T alpha ) override;

  /** this += alpha * x */
   void axpy  ( T alpha, const OoqpVectorBase<T>& x ) override;
  /** this += alpha * x * z */
   void axzpy ( T alpha, const OoqpVectorBase<T>& x, const OoqpVectorBase<T>& z ) override;
  /** this += alpha * x / z */
   void axdzpy( T alpha, const OoqpVectorBase<T>& x, const OoqpVectorBase<T>& z ) override;

   void addConstant( T c ) override;
   void gondzioProjection( T rmin, T rmax ) override;
   T dotProductWith( const OoqpVectorBase<T>& v ) const override;
   T dotProductSelf( T scaleFactor) const override;

  /** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
   */
   T shiftedDotProductWith( T alpha, const OoqpVectorBase<T>& mystep,
					const OoqpVectorBase<T>& yvec,
					T beta, const OoqpVectorBase<T>& ystep ) const override;
   void negate() override;
   void invert() override;
   void invertSave( T zeroReplacementVal = 0.0 ) override;
   void applySqrt() override;
   void roundToPow2() override;

   bool allPositive() const override;
   bool allOf( const std::function<bool(const T&)>& pred ) const override;

   bool matchesNonZeroPattern( const OoqpVectorBase<T>& select ) const override;
   void selectNonZeros( const OoqpVectorBase<T>& select ) override;
   void selectPositive() override;
   void selectNegative() override;
   long long numberOfNonzeros() const override;
   void addSomeConstants( T c, const OoqpVectorBase<T>& select ) override;
   void writefSomeToStream( std::ostream&, const char[], const OoqpVectorBase<T>& ) const override { assert( false && "Not yet implemented" ); };
   void axdzpy( T alpha, const OoqpVectorBase<T>& x, const OoqpVectorBase<T>& z, const OoqpVectorBase<T>& select ) override;

   bool somePositive( const OoqpVectorBase<T>& select ) const override;
   void divideSome( const OoqpVectorBase<T>& div, const OoqpVectorBase<T>& select ) override;
   void copyIntoArray( T[] ) const override { assert( "Not supported" && 0 ); };
   void copyFromArray( const T[] ) override { assert( "Not supported" && 0 ); };
   void copyFromArray( const char[] ) override { assert( "Not supported" && 0 ); };
   virtual void permuteVec0Entries(const std::vector<unsigned int>& permvec);
   virtual void permuteLinkingEntries(const std::vector<unsigned int>& permvec);
   virtual std::vector<T> gatherStochVector() const;

   /** remove entries i for which select[i] == 0 */
   void removeEntries( const OoqpVectorBase<int>& select ) override;

   virtual int getSize() const { return this->n; };
   int getNnzs() const override;

   virtual void recomputeSize();

   virtual bool isRootNodeInSync() const;

   virtual void split( const std::vector<unsigned int>& map_blocks_children, const std::vector<MPI_Comm>& child_comms,
         const std::vector<int>& twolinks_start_in_block = std::vector<int>(), int n_links_in_root = -1);
   virtual StochVectorBase<T>* raiseBorder( int n_vars, bool linking_part, bool shave_top );
   virtual void collapseFromHierarchical( const sData& data_hier, const sTree& tree_hier, VectorType type, bool empty_vec = false );
   virtual void appendHierarchicalToThis( SimpleVectorBase<T>* new_vec, SimpleVectorBase<T>* new_vecl,
         std::vector<StochVectorBase<T>*>& new_children, const sTree& tree_hier, const sData& data_hier, VectorType type, bool empty_vec );

   virtual OoqpVectorBase<T>* getLinkingVecNotHierarchicalTop() const;

   void pushAwayFromZero( double tol, double amount, const OoqpVectorBase<T>* select ) override;
   void getSumCountIfSmall( double tol, double& sum_small, int& n_close, const OoqpVectorBase<T>* select ) const override;
protected:
   StochVectorBase() = default;

private:
   void appendChildrenToThis();
   void pushSmallComplementarityPairs( OoqpVectorBase<T>& other_vec_in, const OoqpVectorBase<T>& select_in, double tol_this, double tol_other, double tol_pairs ) override;
};

/** DUMMY VERSION
 *
 */
template <typename T>
class StochDummyVectorBase : public StochVectorBase<T> {
public:

   StochDummyVectorBase() : StochVectorBase<T>(0, MPI_COMM_NULL) {};
  ~StochDummyVectorBase() override = default;

  void AddChild(StochVectorBase<T>* ) override {};

   StochVectorBase<T>* clone() const override { return new StochDummyVectorBase<T>();}
   StochVectorBase<T>* cloneFull() const override { return new StochDummyVectorBase<T>();}

   void jointCopyFrom(const OoqpVectorBase<T>&, const OoqpVectorBase<T>&, const OoqpVectorBase<T>&) override {};
   void jointCopyTo(OoqpVectorBase<T>&, OoqpVectorBase<T>&, OoqpVectorBase<T>&) const override {};

   bool isKindOf( int kind ) const override {return kind == kStochDummy || kind == kStochVector;}
   bool isZero() const override { return true; };
   void setToZero() override {};
   void setToConstant( T ) override {};
   void randomize( T, T, T* ) override {};
   void copyFrom( const OoqpVectorBase<T>& ) override {};
   void copyFromAbs(const OoqpVectorBase<T>& ) override {};
   double twonorm() const override { return 0.0; }
   T infnorm() const override { return 0.0; }
   T onenorm() const override { return 0.0; }
   void min( T&, int& ) const override {};
   void max( T&, int& ) const override {};
   void absminVecUpdate( OoqpVectorBase<T>& ) const override {};
   void absmaxVecUpdate( OoqpVectorBase<T>& ) const override {};
   void absmin( T& m ) const override { m = std::numeric_limits<T>::infinity(); };
   void absminNonZero( T& m, T ) const override { m = std::numeric_limits<T>::infinity(); };
   T stepbound( const OoqpVectorBase<T>& , T maxStep ) const override { return maxStep; }
   T findBlocking( const OoqpVectorBase<T>&, const OoqpVectorBase<T>&, const OoqpVectorBase<T>&, T maxStep, T*, T*, T*, T*, int&) const override {return maxStep;}

   void findBlocking_pd(const OoqpVectorBase<T>&, const OoqpVectorBase<T>&, const OoqpVectorBase<T>&,
         T&, T&, T&, T&, T&, T&, T&, T&, T&, T&, bool&, bool&) const override {};

   void componentMult( const OoqpVectorBase<T>& ) override {};
   void componentDiv ( const OoqpVectorBase<T>& ) override {};
   bool componentEqual( const OoqpVectorBase<T>& v, T) const override
      { if(!v.isKindOf(kStochDummy)) std::cout << "one should never end up here" << std::endl; return v.isKindOf(kStochDummy); };
   bool componentNotEqual( const T, const T) const override { return true; };
   void setNotIndicatedEntriesToVal( T, const OoqpVectorBase<T>& ) override {};

   void scalarMult( T ) override {};
   void writeToStream( std::ostream&, int ) const override {};
   void writefToStream( std::ostream&, const char[] ) const override {};

   void writeMPSformatOnlyRhs(std::ostream&, const std::string, const OoqpVectorBase<T>* ) const override {};
   void writeMPSformatRhs(std::ostream&, int, const OoqpVectorBase<T>* ) const override {};
   void writeMPSformatBounds(std::ostream&, const OoqpVectorBase<T>* , bool) const override {};
   void writeMPSformatBoundsWithVar(std::ostream&, const std::string, const OoqpVectorBase<T>*, bool) const override {};

   void scale( T ) override {};

  /** this += alpha * x */
   void axpy  ( T, const OoqpVectorBase<T>& ) override {};
  /** this += alpha * x * z */
   void axzpy ( T, const OoqpVectorBase<T>&, const OoqpVectorBase<T>& ) override {};
  /** this += alpha * x / z */
   void axdzpy( T, const OoqpVectorBase<T>&, const OoqpVectorBase<T>& ) override {};

   void addConstant( T ) override {};
   void gondzioProjection( T, T ) override {};
   T dotProductWith( const OoqpVectorBase<T>& ) const override {return 0.0;}
   T dotProductSelf( T ) const override {return 0.0;};

  /** Return the inner product <this + alpha * mystep, yvec + beta * ystep > */
   T shiftedDotProductWith( T , const OoqpVectorBase<T>&, const OoqpVectorBase<T>&, T,  const OoqpVectorBase<T>& ) const override {return 0.0;}
   void negate()override {};
   void invert()override {};
   void invertSave( T ) override {};
   void applySqrt() override {};
   void roundToPow2() override {};
   bool allPositive() const override { return true; };
   bool allOf( const std::function<bool(const T&)>& ) const override { return true; };

   bool matchesNonZeroPattern( const OoqpVectorBase<T>& ) const override {return true;}
   void selectNonZeros( const OoqpVectorBase<T>& ) override {};
   void selectPositive() override {};
   void selectNegative() override {};
   long long numberOfNonzeros() const override { return 0; }
   void addSomeConstants( T, const OoqpVectorBase<T>& ) override {};
   void writefSomeToStream( std::ostream&, const char[], const OoqpVectorBase<T>& ) const override {};
   void axdzpy( T, const OoqpVectorBase<T>&, const OoqpVectorBase<T>&, const OoqpVectorBase<T>& ) override {};

   bool somePositive( const OoqpVectorBase<T>& ) const override { return 1; }
   void divideSome( const OoqpVectorBase<T>&, const OoqpVectorBase<T>& )override {};
   void copyIntoArray( T[] ) const override {};
   void copyFromArray( const T[] ) override {};
   void copyFromArray( const char[] ) override {};

   void removeEntries( const OoqpVectorBase<int>& ) override {};
   void permuteVec0Entries( const std::vector<unsigned int>& ) override {};
   void permuteLinkingEntries( const std::vector<unsigned int>& ) override {};
   std::vector<T> gatherStochVector() const override {return std::vector<T>(0);};

   int getSize() const override { return 0; };
   int getNnzs() const override { return 0; };

   void recomputeSize() override{};

   bool isRootNodeInSync() const override { return true; };

   void split( const std::vector<unsigned int>&, const std::vector<MPI_Comm>&,
         const std::vector<int>&, int ) override {};
   StochVectorBase<T>* raiseBorder( int, bool, bool ) override { assert( 0 && "This should never be attempted" ); return nullptr; };

   void appendHierarchicalToThis( SimpleVectorBase<T>* new_vec, SimpleVectorBase<T>* new_vecl, std::vector<StochVectorBase<T>*>& new_children,
         const sTree& tree_hier, const sData& data_hier, VectorType type, bool empty_vec ) override;

   OoqpVectorBase<T>* getLinkingVecNotHierarchicalTop() const override { assert( false && "Should not end up here"); return nullptr; };

   void pushAwayFromZero( double, double, const OoqpVectorBase<T>* ) override {};
   void getSumCountIfSmall( double, double&, int&, const OoqpVectorBase<T>* ) const override {};

   void pushSmallComplementarityPairs( OoqpVectorBase<T>&, const OoqpVectorBase<T>&, double, double, double) override {};
};

#endif
