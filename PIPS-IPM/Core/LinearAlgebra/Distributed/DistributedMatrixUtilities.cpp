//
// Created by nils-christian on 27.05.21.
//

#include "DistributedMatrixUtilities.h"

SparseMatrix* getSparseGenMatrixFromStochMat(const DistributedMatrix& sMat, int smat_node, BlockType block_type) {
   assert(-1 <= smat_node && smat_node < static_cast<int>(sMat.children.size()));

   if (smat_node == -1) {
      if (block_type == BL_MAT) {
         assert(sMat.Blmat->is_a(kSparseGenMatrix));
         return dynamic_cast<SparseMatrix*>(sMat.Blmat.get());
      } else {
         assert(block_type == B_MAT);
         assert(sMat.Bmat->is_a(kSparseGenMatrix));
         return dynamic_cast<SparseMatrix*>(sMat.Bmat.get());
      }
   } else {
      if (block_type == A_MAT) {
         assert(sMat.children[smat_node]->Amat->is_a(kSparseGenMatrix));
         return dynamic_cast<SparseMatrix*>(sMat.children[smat_node]->Amat.get());
      } else if (block_type == B_MAT) {
         assert(sMat.children[smat_node]->Bmat->is_a(kSparseGenMatrix));
         return dynamic_cast<SparseMatrix*>(sMat.children[smat_node]->Bmat.get());
      } else if (block_type == BL_MAT) {
         assert(sMat.children[smat_node]->Blmat->is_a(kSparseGenMatrix));
         return dynamic_cast<SparseMatrix*>(sMat.children[smat_node]->Blmat.get());
      }
   }
   return nullptr;
}

SparseSymmetricMatrix& getSparseSymmetricMatrixFromStochMat(const DistributedSymmetricMatrix& sMat, int smat_node)
{
   if (smat_node == -1) {
      assert(sMat.diag);
      return dynamic_cast<SparseSymmetricMatrix&>(*sMat.diag);
   } else {
      assert(0 <= smat_node && smat_node <= static_cast<int>(sMat.children.size()));
      assert(sMat.children[smat_node]);
      assert(sMat.children[smat_node]->diag);
      return dynamic_cast<SparseSymmetricMatrix&>(*sMat.children[smat_node]);
   }
}
