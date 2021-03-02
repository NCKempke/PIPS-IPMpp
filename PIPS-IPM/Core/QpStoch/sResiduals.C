#include "sResiduals.h"
#include "sTree.h"
#include "StochVector.h"

sResiduals::sResiduals( OoqpVector * rQ_,    OoqpVector * rA_,
			OoqpVector * rC_,    OoqpVector * rz_, 
			OoqpVector * rt_,    OoqpVector * rlambda_, 
			OoqpVector * ru_,    OoqpVector * rpi_, 
			OoqpVector * rv_,    OoqpVector * rgamma_, 
			OoqpVector * rw_,    OoqpVector * rphi_, 
			OoqpVector * ixlow_, double nxlowGlobal,
			OoqpVector * ixupp_, double nxuppGlobal,
			OoqpVector * iclow_, double mclowGlobal, 
			OoqpVector * icupp_, double mcuppGlobal)
{
  SpReferTo( ixlow, ixlow_ );
  nxlow = nxlowGlobal;

  SpReferTo( ixupp, ixupp_ );
  nxupp = nxuppGlobal;

  SpReferTo( iclow, iclow_ );
  mclow = mclowGlobal;

  SpReferTo( icupp, icupp_ );
  mcupp = mcuppGlobal;

  SpReferTo( rQ, rQ_ );
  SpReferTo( rA, rA_ );
  SpReferTo( rC, rC_ );
  SpReferTo( rz, rz_ );
  SpReferTo( rt , rt_ );

  SpReferTo( rlambda, rlambda_ );
  SpReferTo( ru    , ru_ );
  SpReferTo( rpi   , rpi_ );
  SpReferTo( rv    , rv_ );
  SpReferTo( rgamma, rgamma_ );
  SpReferTo( rw  , rw_ );
  SpReferTo( rphi, rphi_ );

  createChildren();
}


sResiduals::sResiduals( const sTree* tree, OoqpVector * ixlow_, OoqpVector * ixupp_, OoqpVector * iclow_, OoqpVector * icupp_ )
{

  SpReferTo( ixlow, ixlow_ );
  nxlow = ixlow->numberOfNonzeros();

  SpReferTo( ixupp, ixupp_ );
  nxupp = ixupp->numberOfNonzeros();

  SpReferTo( iclow, iclow_ );
  mclow = iclow->numberOfNonzeros();

  SpReferTo( icupp, icupp_ );
  mcupp = icupp->numberOfNonzeros();

  const bool empty_vector = true;
  rQ = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector() );
  rA = OoqpVectorHandle( (OoqpVector*) tree->newDualYVector() );
  rC = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );

  rz = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );
  if ( mclow > 0 ) {
    rt      = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );
    rlambda = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );
  } else {
    rt      = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector(empty_vector) );
    rlambda = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector(empty_vector) );
  }

  if ( mcupp > 0 ) {
    ru     = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );
    rpi    = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );
  } else {
    ru     = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector(empty_vector) );
    rpi    = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector(empty_vector) );
  }

  if( nxlow > 0 ) {
    rv     = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector() );
    rgamma = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector() );
  } else {
    rv     = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector(empty_vector) );
    rgamma = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector(empty_vector) );
  }

  if( nxupp > 0 ) {
    rw   = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector() );
    rphi = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector() );
  } else {
    rw   = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector(empty_vector) );
    rphi = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector(empty_vector) );
  }
  
  createChildren();
}

sResiduals::sResiduals( const sResiduals& res ) : QpGenResiduals( res )
{
   for(unsigned int i = 0; i < res.children.size(); ++i)
   {
       children.push_back( new sResiduals( *res.children[i]) );
   }
}

sResiduals::~sResiduals()
{
   for( sResiduals* child : children )
      delete child;
}

void sResiduals::AddChild(sResiduals* child)
{
  children.push_back(child);
}

void sResiduals::createChildren()
{
  StochVector& rQSt = dynamic_cast<StochVector&>(*rQ);

  StochVector& rASt = dynamic_cast<StochVector&>(*rA); 
  StochVector& rCSt = dynamic_cast<StochVector&>(*rC);    
  StochVector& rzSt = dynamic_cast<StochVector&>(*rz); 
  StochVector& rtSt = dynamic_cast<StochVector&>(*rt);    
  StochVector& rlambdaSt = dynamic_cast<StochVector&>(*rlambda); 
  StochVector& ruSt = dynamic_cast<StochVector&>(*ru);    
  StochVector& rpiSt = dynamic_cast<StochVector&>(*rpi); 
  StochVector& rvSt = dynamic_cast<StochVector&>(*rv);    
  StochVector& rgammaSt = dynamic_cast<StochVector&>(*rgamma); 
  StochVector& rwSt = dynamic_cast<StochVector&>(*rw);    
  StochVector& rphiSt = dynamic_cast<StochVector&>(*rphi); 
  StochVector& ixlowSt = dynamic_cast<StochVector&>(*ixlow);
  StochVector& ixuppSt = dynamic_cast<StochVector&>(*ixupp);
  StochVector& iclowSt = dynamic_cast<StochVector&>(*iclow);
  StochVector& icuppSt = dynamic_cast<StochVector&>(*icupp);

  const size_t nChildren = rASt.children.size();
  for (size_t it = 0; it < nChildren; it++) {

      assert(nChildren == rASt.children.size());
      assert(nChildren == rzSt.children.size());
      assert(nChildren == rCSt.children.size());
      assert(nChildren == rtSt.children.size());
      assert(nChildren == rlambdaSt.children.size());
      assert(nChildren == ruSt.children.size());
      assert(nChildren == rpiSt.children.size());
      assert(nChildren == rvSt.children.size());
      assert(nChildren == rgammaSt.children.size());
      assert(nChildren == rwSt.children.size());
      assert(nChildren == rphiSt.children.size());
      assert(nChildren == ixlowSt.children.size());
      assert(nChildren == ixuppSt.children.size());
      assert(nChildren == iclowSt.children.size());
      assert(nChildren == icuppSt.children.size());
 
    AddChild(new sResiduals(rQSt.children[it],    rASt.children[it],
			    rCSt.children[it],    rzSt.children[it], 
			    rtSt.children[it],    rlambdaSt.children[it], 
			    ruSt.children[it],    rpiSt.children[it], 
			    rvSt.children[it],    rgammaSt.children[it], 
			    rwSt.children[it],    rphiSt.children[it], 
			    ixlowSt.children[it], nxlow,
			    ixuppSt.children[it], nxupp,
			    iclowSt.children[it], mclow, 
			    icuppSt.children[it], mcupp));
  }
}

void sResiduals::collapseHierarchicalStructure(const sData& data_hier, const sTree* tree_hier, OoqpVectorHandle ixlow_, OoqpVectorHandle ixupp_,
      OoqpVectorHandle iclow_, OoqpVectorHandle icupp_)
{
   dynamic_cast<StochVector&>(*rQ).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);

   const bool empty_vec = true;
   if( nxlow > 0 )
   {
      dynamic_cast<StochVector&>(*rv).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);
      dynamic_cast<StochVector&>(*rgamma).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);
   }
   else
   {
      dynamic_cast<StochVector&>(*rv).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL, empty_vec);
      dynamic_cast<StochVector&>(*rgamma).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL, empty_vec);
   }

   if( nxupp > 0 )
   {
      dynamic_cast<StochVector&>(*rw).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);
      dynamic_cast<StochVector&>(*rphi).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL);
   }
   else
   {
      dynamic_cast<StochVector&>(*rw).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL, empty_vec);
      dynamic_cast<StochVector&>(*rphi).collapseFromHierarchical(data_hier, *tree_hier, VectorType::PRIMAL, empty_vec);
   }

   dynamic_cast<StochVector&>(*rA).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Y);

   dynamic_cast<StochVector&>(*rC).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
   dynamic_cast<StochVector&>(*rz).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);

   if( mcupp > 0 )
   {
      dynamic_cast<StochVector&>(*ru).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
      dynamic_cast<StochVector&>(*rpi).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
   }
   else
   {
      dynamic_cast<StochVector&>(*ru).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z, empty_vec);
      dynamic_cast<StochVector&>(*rpi).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z, empty_vec);
   }

   if ( mclow > 0 )
   {
      dynamic_cast<StochVector&>(*rt).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
      dynamic_cast<StochVector&>(*rlambda).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z);
   }
   else
   {
      dynamic_cast<StochVector&>(*rt).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z, empty_vec);
      dynamic_cast<StochVector&>(*rlambda).collapseFromHierarchical(data_hier, *tree_hier, VectorType::DUAL_Z, empty_vec);
   }


   ixlow = ixlow_;
   ixupp = ixupp_;
   iclow = iclow_;
   icupp = icupp_;

   for (size_t c = 0; c < children.size(); c++)
     delete children[c];

   children.clear();
   createChildren();
}

void sResiduals::permuteVec0Entries( const std::vector<unsigned int>& perm, bool resids_only )
{
   if( !resids_only )
   {
      dynamic_cast<StochVector&>(*ixlow).permuteVec0Entries(perm);
      dynamic_cast<StochVector&>(*ixupp).permuteVec0Entries(perm);
   }

   dynamic_cast<StochVector&>(*rQ).permuteVec0Entries(perm);

   if( nxlow > 0 )
   {
      dynamic_cast<StochVector&>(*rv).permuteVec0Entries(perm);
      dynamic_cast<StochVector&>(*rgamma).permuteVec0Entries(perm);
   }

   if( nxupp > 0 )
   {
      dynamic_cast<StochVector&>(*rw).permuteVec0Entries(perm);
      dynamic_cast<StochVector&>(*rphi).permuteVec0Entries(perm);
   }
}

void sResiduals::permuteEqLinkingEntries( const std::vector<unsigned int>& perm )
{
   dynamic_cast<StochVector&>(*rA).permuteLinkingEntries(perm);
}

void sResiduals::permuteIneqLinkingEntries( const std::vector<unsigned int>& perm, bool resids_only )
{
   if( !resids_only )
   {
      dynamic_cast<StochVector&>(*iclow).permuteLinkingEntries(perm);
      dynamic_cast<StochVector&>(*icupp).permuteLinkingEntries(perm);
   }

   dynamic_cast<StochVector&>(*rC).permuteLinkingEntries(perm);
   dynamic_cast<StochVector&>(*rz).permuteLinkingEntries(perm);

   if( mcupp > 0 )
   {
      dynamic_cast<StochVector&>(*ru).permuteLinkingEntries(perm);
      dynamic_cast<StochVector&>(*rpi).permuteLinkingEntries(perm);
   }

   if ( mclow > 0 )
   {
      dynamic_cast<StochVector&>(*rt).permuteLinkingEntries(perm);
      dynamic_cast<StochVector&>(*rlambda).permuteLinkingEntries(perm);
   }
}

bool sResiduals::isRootNodeInSync() const
{
   bool in_sync = true;
   const int my_rank = PIPS_MPIgetRank( MPI_COMM_WORLD );
   if( !dynamic_cast<const StochVector&>(*rQ).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "rQ not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*rC).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "rC not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*rA).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "rA not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*rz).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "rz not in sync" << std::endl;
      in_sync = false;
   }

   if( !dynamic_cast<const StochVector&>(*rt).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "rt not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*rlambda).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "rlambda not in sync" << std::endl;
      in_sync = false;
   }

   if( !dynamic_cast<const StochVector&>(*ru).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "ru not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*rpi).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "rpi not in sync" << std::endl;
      in_sync = false;
   }

   if( !dynamic_cast<const StochVector&>(*rv).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "rv not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*rgamma).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "rgamma not in sync" << std::endl;
      in_sync = false;
   }

   if( !dynamic_cast<const StochVector&>(*rw).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "rw not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*rphi).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "rphi not in sync" << std::endl;
      in_sync = false;
   }

   MPI_Barrier(MPI_COMM_WORLD);
   return in_sync;

}
