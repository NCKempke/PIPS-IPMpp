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

SparseMatrix* getSparseGenMatrixFromStochMat(const DistributedMatrix& sMat, int smat_node, BlockType block_type);
SparseSymmetricMatrix& getSparseSymmetricDiagFromStochMat(const DistributedSymmetricMatrix& sMat, int smat_node);
SparseMatrix& getSparseBorderFromStochMat(const DistributedSymmetricMatrix& sMat, int smat_node);

#endif /* DISTRIBUTEDMATRIXUTILITIES_H */
