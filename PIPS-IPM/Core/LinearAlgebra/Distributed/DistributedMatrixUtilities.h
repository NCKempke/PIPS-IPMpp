/*
 * StochMatrixUtilities.h
 *
 *  Created on: 27.11.2019
 *      Author: bzfkempk
 */
#ifndef DISTRIBUTEDMATRIXUTILITIES_H
#define DISTRIBUTEDMATRIXUTILITIES_H

#include "SparseMatrix.h"
#include "DistributedMatrix.h"
#include "SystemType.h"

#include <vector>

inline SparseMatrix* getSparseGenMatrixFromStochMat(const DistributedMatrix& sMat, int smat_node, BlockType block_type) {
   assert(-1 <= smat_node && smat_node < static_cast<int>(sMat.children.size()));

   if (smat_node == -1) {
      if (block_type == BL_MAT) {
         assert(sMat.Blmat->is_a(kSparseGenMatrix));
         return dynamic_cast<SparseMatrix*>(sMat.Blmat);
      }
      else {
         assert(block_type == B_MAT);
         assert(sMat.Bmat->is_a(kSparseGenMatrix));
         return dynamic_cast<SparseMatrix*>(sMat.Bmat);
      }
   }
   else {
      if (block_type == A_MAT) {
         assert(sMat.children[smat_node]->Amat->is_a(kSparseGenMatrix));
         return dynamic_cast<SparseMatrix*>(sMat.children[smat_node]->Amat);
      }
      else if (block_type == B_MAT) {
         assert(sMat.children[smat_node]->Bmat->is_a(kSparseGenMatrix));
         return dynamic_cast<SparseMatrix*>(sMat.children[smat_node]->Bmat);
      }
      else if (block_type == BL_MAT) {
         assert(sMat.children[smat_node]->Blmat->is_a(kSparseGenMatrix));
         return dynamic_cast<SparseMatrix*>(sMat.children[smat_node]->Blmat);
      }
   }
   return nullptr;
}

#endif /* DISTRIBUTEDMATRIXUTILITIES_H */
