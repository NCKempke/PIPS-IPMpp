#include "sData.h"
#include "sTree.h"
#include "sTreeCallbacks.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "SparseLinearAlgebraPackage.h"
#include "mpi.h"

#include "pipsport.h"
#include "StochOptions.h"
#include "BorderedSymMatrix.h"

#include <iomanip>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <functional>

static
std::vector<unsigned int> getInversePermutation(const std::vector<unsigned int>& perm)
{
   size_t size = perm.size();
   std::vector<unsigned int> perm_inv(size, 0);

   for( size_t i = 0; i < size; i++ )
      perm_inv[perm[i]] = i;

   return perm_inv;
}

static bool blockIsInRange(int block, int blocksStart, int blocksEnd)
{
   return ((block >= (blocksStart - 1) && block < blocksEnd) || block == -1);
}

static int nnzTriangular(int size)
{
   assert(size >= 0);
   return ((1 + size) * size) / 2;
}

static void appendRowDense(int start, int end, int& nnz, int* jcolM)
{
   assert(jcolM);
   assert(nnz >= 0);
   assert(start >= 0 && end >= start);

   for( int i = start; i < end; i++ )
      jcolM[nnz++] = i;
}

static void appendRowSparse(int startColIdx, int endColIdx, int colOffset, const int* jcolM_append, int& nnz, int* jcolM)
{
   assert(jcolM);
   assert(nnz >= 0);
   assert(startColIdx >= 0 && startColIdx <= endColIdx);

   for( int c = startColIdx; c < endColIdx; c++ )
      jcolM[nnz++] = colOffset + jcolM_append[c];
}


static int appendDiagBlocks(const std::vector<int>& linkStartBlockId, const std::vector<int>& linkStartBlockLengths, int borderstart, int bordersize, int rowSC,
                                          int rowBlock, int& blockStartrow, int& nnz, int* jcolM)
{
   assert(rowBlock >= blockStartrow && blockStartrow >= 0 && borderstart >= 0 && bordersize >= 0 && nnz >= 0);

   const int block = linkStartBlockId[rowBlock];
   const int currlength = (block >= 0) ? linkStartBlockLengths[block] : bordersize;

   assert(currlength >= 1);

   // add diagonal block (possibly up to the border)

   int rownnz = currlength - (rowBlock - blockStartrow);

   for( int i = 0; i < rownnz; ++i )
      jcolM[nnz++] = rowSC + i;

   // with offdiagonal blocks?
   if( block >= 0 )
   {
      // add right off-diagonal block and border part

      const int nextlength = linkStartBlockLengths[block + 1];

      assert(nextlength >= 0);
      assert(block != int(linkStartBlockLengths.size()) - 2 || nextlength == 0);

      for( int i = rownnz; i < rownnz + nextlength; ++i )
         jcolM[nnz++] = rowSC + i;

      rownnz += nextlength + bordersize;

      for( int i = borderstart; i < borderstart + bordersize; ++i )
         jcolM[nnz++] = i;

      // last row of current block?
      if( rowBlock + 1 == blockStartrow + currlength )
         blockStartrow = rowBlock + 1;
   }

   return rownnz;
}

static void appendDiagBlocksDist(const std::vector<int>& linkStartBlockId, const std::vector<int>& linkStartBlockLengths, int borderstart, int bordersize, int rowSC,
                                          int rowBlock, int blocksStart, int blocksEnd, int& blockStartrow, int& nnz, int* jcolM)
{
   assert(rowBlock >= blockStartrow && blockStartrow >= 0 && borderstart >= 0 && bordersize >= 0 && nnz >= 0);

   const int block = linkStartBlockId[rowBlock];
   const int lastBlock = blocksEnd - 1;

   if( blockIsInRange(block, blocksStart, blocksEnd) )
   {
      const int currlength = (block >= 0) ? linkStartBlockLengths[block] : bordersize;
      const int rownnz = currlength - (rowBlock - blockStartrow);

      assert(currlength >= 1);

      // add diagonal block (possibly up to the border)
      for( int i = 0; i < rownnz; ++i )
         jcolM[nnz++] = rowSC + i;

      // with off-diagonal blocks? (at sparse part)
      if( block >= 0 )
      {
         // add right off-diagonal block and border part

         const int nextlength = (block == lastBlock) ? 0 : linkStartBlockLengths[block + 1];

         assert(nextlength >= 0);
         assert(block != int(linkStartBlockLengths.size()) - 1 || nextlength == 0);

         for( int i = rownnz; i < rownnz + nextlength; ++i )
            jcolM[nnz++] = rowSC + i;

         for( int i = borderstart; i < borderstart + bordersize; ++i )
            jcolM[nnz++] = i;

         // last row of current block?
         if( rowBlock + 1 == blockStartrow + currlength )
            blockStartrow = rowBlock + 1;
      }
   }

   // at sparse part?
   if( block >= 0 )
   {
      const int currlength = linkStartBlockLengths[block];
      assert(currlength >= 1);

      // last row of current block?
      if( rowBlock + 1 == blockStartrow + currlength )
         blockStartrow = rowBlock + 1;
   }
}


static int appendMixedBlocks(const std::vector<int>& linkStartBlockId_Left,
      const std::vector<int>& linkStartBlockId_Right,
      const std::vector<int>& linkStartBlockLengths_Left,
      const std::vector<int>& linkStartBlockLengths_Right,
      int colStartIdxSC, int bordersize_cols,
      int rowIdx, int& colIdxOffset, int& rowBlockStartIdx, int& nnz, int* jcolM)
{
   assert(rowIdx >= rowBlockStartIdx && rowBlockStartIdx >= 0 && colStartIdxSC >= 0 && nnz >= 0);
   assert(bordersize_cols >= 0 && colIdxOffset >= 0);
   assert(linkStartBlockLengths_Left.size() == linkStartBlockLengths_Right.size());

   const int nCols = int(linkStartBlockId_Right.size());
   const int block = linkStartBlockId_Left[rowIdx];

   assert(nCols >= bordersize_cols);

   int rownnz;

   // sparse row?
   if( block >= 0 )
   {
      const int length_Right = linkStartBlockLengths_Right[block];
      const int colStartIdxBorderSC = colStartIdxSC + nCols - bordersize_cols;
      int colStartIdx = colStartIdxSC + colIdxOffset;
      rownnz = 0;

      assert(length_Right >= 0);

      // 1) left off-diagonal block (not for first block)
      if( block >= 1 )
      {
         const int prevlength_Right = linkStartBlockLengths_Right[block - 1];
         assert(prevlength_Right >= 0);

         for( int i = 0; i < prevlength_Right; ++i )
            jcolM[nnz++] = (colStartIdx + i);

         rownnz += prevlength_Right;
         colStartIdx += prevlength_Right;
      }

      // 2) diagonal block
      for( int i = 0; i < length_Right; ++i )
         jcolM[nnz++] = (colStartIdx + i);

      rownnz += length_Right;
      colStartIdx += length_Right;

      // 3) right off-diagonal block (not for last block)
      if( int(linkStartBlockLengths_Left.size()) != block + 1 )
      {
         const int nextlength_Right = linkStartBlockLengths_Right[block + 1];

         for( int i = 0; i < nextlength_Right; ++i )
            jcolM[nnz++] = (colStartIdx + i);

         rownnz += nextlength_Right;
         colStartIdx += nextlength_Right;
      }

      assert(colStartIdx <= colStartIdxBorderSC);

      // 4) right border
      for( int i = 0; i < bordersize_cols; ++i )
         jcolM[nnz++] = colStartIdxBorderSC + i;

      rownnz += bordersize_cols;

      // last row of current block?
      if( rowIdx + 1 == rowBlockStartIdx + linkStartBlockLengths_Left[block] )
      {
         rowBlockStartIdx = rowIdx + 1;

         if( block >= 1 )
            colIdxOffset += linkStartBlockLengths_Right[block - 1];
      }
      else
      {
         assert(block == linkStartBlockId_Left[rowIdx + 1]);
      }
   }
   else
   {
      // append fully dense row
      for( int i = 0; i < nCols; ++i )
         jcolM[nnz++] = colStartIdxSC + i;

      rownnz = nCols;
   }

   return rownnz;
}


static void appendMixedBlocksDist(const std::vector<int>& linkStartBlockId_Left,
      const std::vector<int>& linkStartBlockId_Right,
      const std::vector<int>& linkStartBlockLengths_Left,
      const std::vector<int>& linkStartBlockLengths_Right,
      int colStartIdxSC, int bordersize_cols,
      int rowIdx, int blocksStart, int blocksEnd,
      int& colIdxOffset, int& rowBlockStartIdx, int& nnz, int* jcolM)
{
   assert(rowIdx >= rowBlockStartIdx && rowBlockStartIdx >= 0 && colStartIdxSC >= 0 && nnz >= 0);
   assert(bordersize_cols >= 0 && colIdxOffset >= 0);
   assert(linkStartBlockLengths_Left.size() == linkStartBlockLengths_Right.size());

   const int block = linkStartBlockId_Left[rowIdx];
   const bool blockInRange = blockIsInRange(block, blocksStart, blocksEnd);
   const int nCols = int(linkStartBlockId_Right.size());

   assert(nCols >= bordersize_cols);

   // sparse row?
   if( block >= 0 )
   {
      if( blockInRange )
      {
         const int length_Right = linkStartBlockLengths_Right[block];
         const int colStartIdxBorderSC = colStartIdxSC + nCols - bordersize_cols;
         int colStartIdx = colStartIdxSC + colIdxOffset;

         assert(length_Right >= 0);

         // 1) left off-diagonal block (not for first block)
         if( block >= 1 && block != blocksStart - 1 )
         {
            const int prevlength_Right = linkStartBlockLengths_Right[block - 1];
            assert(prevlength_Right >= 0);

            for( int i = 0; i < prevlength_Right; ++i )
               jcolM[nnz++] = (colStartIdx + i);
         }

         if( block >= 1 )
         {
            const int prevlength_Right = linkStartBlockLengths_Right[block - 1];
            assert(prevlength_Right >= 0);

            colStartIdx += prevlength_Right;
         }

         // 2) diagonal block
         for( int i = 0; i < length_Right; ++i )
            jcolM[nnz++] = (colStartIdx + i);

         colStartIdx += length_Right;

         // 3) right off-diagonal block (not for last block)
         if( block != blocksEnd - 1 )
         {
            const int nextlength_Right = linkStartBlockLengths_Right[block + 1];

            for( int i = 0; i < nextlength_Right; ++i )
               jcolM[nnz++] = (colStartIdx + i);

            colStartIdx += nextlength_Right;
         }

         assert(colStartIdx <= colStartIdxBorderSC);

         // 4) right border
         for( int i = 0; i < bordersize_cols; ++i )
            jcolM[nnz++] = colStartIdxBorderSC + i;
      }

      // last row of current block?
      if( rowIdx + 1 == rowBlockStartIdx + linkStartBlockLengths_Left[block] )
      {
         rowBlockStartIdx = rowIdx + 1;

         if( block >= 1 )
            colIdxOffset += linkStartBlockLengths_Right[block - 1];
      }
      else
      {
         assert(block == linkStartBlockId_Left[rowIdx + 1]);
      }
   }
   else if( blockInRange )
   {
      assert(block == -1);

      // append dense row, but skip entries not in range
      for( int i = 0; i < nCols; ++i )
      {
         const int blockRight = linkStartBlockId_Right[i];
         const bool blockRightInRange = blockIsInRange(blockRight, blocksStart, blocksEnd);

         if( blockRightInRange )
            jcolM[nnz++] = colStartIdxSC + i;
      }
   }
}


void sData::getSCrangeMarkers(int blocksStart, int blocksEnd, int& local2linksStartEq, int& local2linksEndEq,
      int& local2linksStartIneq, int& local2linksEndIneq)
{
   const int blocksStartReal = (blocksStart > 0) ? (blocksStart - 1) : blocksStart;
   const int nx0 = getLocalnx();
   const int my0 = getLocalmy();
   const int myl = getLocalmyl();
   local2linksStartEq = nx0 + my0;
   local2linksStartIneq = nx0 + my0 + myl;

   for( int block = 0; block < blocksStartReal; ++block )
   {
      const int lengthEq = linkStartBlockLengthsA[block];
      const int lengthIneq = linkStartBlockLengthsC[block];

      assert(lengthEq >= 0 && lengthIneq >= 0);

      local2linksStartEq += lengthEq;
      local2linksStartIneq += lengthIneq;
   }

   local2linksEndEq = local2linksStartEq + getSCdiagBlocksNRows(linkStartBlockLengthsA, blocksStart, blocksEnd);
   local2linksEndIneq = local2linksStartIneq + getSCdiagBlocksNRows(linkStartBlockLengthsC, blocksStart, blocksEnd);
}

void sData::getSCrangeMarkersMy(int blocksStart, int blocksEnd, int& local2linksStartEq, int& local2linksEndEq,
      int& local2linksStartIneq, int& local2linksEndIneq)
{
   const int nx0 = getLocalnx();
   const int my0 = getLocalmy();
   const int myl = getLocalmyl();
   local2linksStartEq = nx0 + my0;
   local2linksStartIneq = nx0 + my0 + myl;

   for( int block = 0; block < blocksStart; ++block )
   {
      const int lengthEq = linkStartBlockLengthsA[block];
      const int lengthIneq = linkStartBlockLengthsC[block];

      assert(lengthEq >= 0 && lengthIneq >= 0);

      local2linksStartEq += lengthEq;
      local2linksStartIneq += lengthIneq;
   }

   local2linksEndEq = local2linksStartEq + getSCdiagBlocksNRowsMy(linkStartBlockLengthsA, blocksStart, blocksEnd);
   local2linksEndIneq = local2linksStartIneq + getSCdiagBlocksNRowsMy(linkStartBlockLengthsC, blocksStart, blocksEnd);
}


int sData::getSCdiagBlocksNRows(const std::vector<int>& linkStartBlockLengths,
      int blocksStart, int blocksEnd)
{
   assert(blocksStart >= 0 && blocksStart < blocksEnd);
   assert(blocksEnd <= int(linkStartBlockLengths.size()));

   int nRowsRange = 0;
   const int blocksStartReal = (blocksStart > 0) ? (blocksStart - 1) : blocksStart;

   // main loop, going over specified 2-link blocks
   for( int block = blocksStartReal; block < blocksEnd; ++block )
   {
      const int length = linkStartBlockLengths[block];
      assert(length >= 0);

      nRowsRange += length;
   }
   return nRowsRange;
}


int sData::getSCdiagBlocksNRowsMy(const std::vector<int>& linkStartBlockLengths,
      int blocksStart, int blocksEnd)
{
   assert(blocksStart >= 0 && blocksStart < blocksEnd);
   assert(blocksEnd <= int(linkStartBlockLengths.size()));

   int nRowsRange = 0;

   // main loop, going over specified 2-link blocks
   for( int block = blocksStart; block < blocksEnd; ++block )
   {
      const int length = linkStartBlockLengths[block];
      assert(length >= 0);

      nRowsRange += length;
   }
   return nRowsRange;
}

int sData::getSCdiagBlocksNRows(const std::vector<int>& linkStartBlockLengths)
{
   return (getSCdiagBlocksNRows(linkStartBlockLengths, 0, int(linkStartBlockLengths.size())));
}

int sData::getSCdiagBlocksMaxNnz(size_t nRows, const std::vector<int>& linkStartBlockLengths)
{
   const size_t nBlocks = linkStartBlockLengths.size();
   size_t nRowsSparse = 0;

   int nnz = 0;

   // main loop, going over all 2-link blocks
   for( size_t block = 0; block < nBlocks; ++block )
   {
      if( linkStartBlockLengths[block] == 0 )
         continue;

      const int length = linkStartBlockLengths[block];
      const int nextlength = linkStartBlockLengths[block + 1];

      assert(length > 0);
      assert(nextlength >= 0);
      assert(block != linkStartBlockLengths.size() - 2 || nextlength == 0);

      nRowsSparse += size_t(length);

      // diagonal block
      nnz += nnzTriangular(length);

      // (one) off-diagonal block
      nnz += length * nextlength;
   }

   // any rows left?
   if( nRowsSparse < nRows )
   {
      const size_t nRowsDense = nRows - nRowsSparse;
      nnz += nnzTriangular(nRowsDense) + nRowsDense * nRowsSparse;
   }

   return nnz;
}


int sData::getSCdiagBlocksMaxNnzDist(size_t nRows, const std::vector<int>& linkStartBlockLengths, int blocksStart, int blocksEnd)
{
#ifndef NDEBUG
   const int nblocks = int(linkStartBlockLengths.size());
#endif
   assert(blocksStart >= 0);
   assert(blocksStart < blocksEnd);
   assert(blocksEnd <= nblocks);
   assert(nblocks >= 2);

   const int nRowsSparse = getSCdiagBlocksNRows(linkStartBlockLengths);
   const int nRowsSparseRange = getSCdiagBlocksNRows(linkStartBlockLengths, blocksStart, blocksEnd);

   int nnz = 0;

   // main loop, going over specified 2-link blocks
   for( int block = blocksStart; block < blocksEnd; ++block )
   {
      const int length = linkStartBlockLengths[block];

      if( length == 0 )
         continue;

      const int prevlength = (block == 0) ? 0 : linkStartBlockLengths[block - 1];

      assert(length > 0);
      assert(prevlength >= 0);
      assert(block != nblocks - 1); // length should be 0 for last block

      // diagonal block
      nnz += nnzTriangular(length);

      // above off-diagonal block
      nnz += prevlength * length;
   }

   if( blocksStart > 0 )
   {
      const int prevlength = linkStartBlockLengths[blocksStart - 1];

      nnz += nnzTriangular(prevlength);
   }

   // any rows left?
   if( nRowsSparse < int(nRows) )
   {
      const int nRowsDense = int(nRows) - nRowsSparse;
      nnz += nnzTriangular(nRowsDense);
      nnz += nRowsDense * nRowsSparseRange;
   }

   return nnz;
}

int sData::getSCmixedBlocksMaxNnz(size_t nRows, size_t nCols,
      const std::vector<int>& linkStartBlockLength_Left,
      const std::vector<int>& linkStartBlockLength_Right)
{
   assert(linkStartBlockLength_Left.size() == linkStartBlockLength_Right.size());

   const size_t nBlocks = linkStartBlockLength_Left.size();
   size_t nRowsSparse = 0;
   size_t nColsSparse = 0;

   int nnz = 0;

   // main loop, going over all 2-link blocks
   for( size_t block = 0; block < nBlocks; ++block )
   {
      const int length_Left = linkStartBlockLength_Left[block];
      const int length_Right = linkStartBlockLength_Right[block];
      assert(length_Left >= 0 && length_Right >= 0);

      nRowsSparse += size_t(length_Left);
      nColsSparse += size_t(length_Right);

      // diagonal block
      nnz += length_Left * length_Right;

      if( block == 0 )
         continue;

      const int prevlength_Left = linkStartBlockLength_Left[block - 1];
      const int prevlength_Right = linkStartBlockLength_Right[block - 1];

      assert(prevlength_Left >= 0 && prevlength_Right >= 0);

      // left off-diagonal block
      nnz += length_Left * prevlength_Right;

      // upper off-diagonal block
      nnz += prevlength_Left * length_Right;
   }

   // dense border?
   if( nRowsSparse < nRows || nColsSparse < nCols )
   {
      assert(nRowsSparse <= nRows && nColsSparse <= nCols);

      const size_t nRowsDense = nRows - nRowsSparse;
      const size_t nColsDense = nCols - nColsSparse;

      nnz += nRowsDense * nColsSparse; // lower border part without right border
      nnz += nColsDense * nRows;       // complete right border
   }

   return nnz;
}


int sData::getSCmixedBlocksMaxNnzDist(size_t nRows, size_t nCols,
      const std::vector<int>& linkStartBlockLength_Left,
      const std::vector<int>& linkStartBlockLength_Right,
      int blocksStart, int blocksEnd)
{
   assert(linkStartBlockLength_Left.size() == linkStartBlockLength_Right.size());
   assert(blocksStart >= 0);
   assert(blocksStart < blocksEnd);

#ifndef NDEBUG
   const int nBlocks = int(linkStartBlockLength_Left.size());
#endif
   const int blockLast = blocksEnd - 1;
   const int blocksStartReal =  (blocksStart > 0) ? (blocksStart - 1) : 0;
   const int nRowsSparse = getSCdiagBlocksNRows(linkStartBlockLength_Left);
   const int nColsSparse = getSCdiagBlocksNRows(linkStartBlockLength_Right);
   const int nRowsSparseRange = getSCdiagBlocksNRows(linkStartBlockLength_Left, blocksStart, blocksEnd);
   const int nColsSparseRange = getSCdiagBlocksNRows(linkStartBlockLength_Right, blocksStart, blocksEnd);

   assert(nBlocks > 1 && blocksEnd <= nBlocks);
   assert(linkStartBlockLength_Left[nBlocks - 1] == 0 && linkStartBlockLength_Right[nBlocks - 1] == 0);

   int nnz = 0;

   // main loop, going over all 2-link blocks
   for( int block = blocksStartReal; block < blocksEnd; ++block )
   {
      const int length_Left = linkStartBlockLength_Left[block];
      const int length_Right = linkStartBlockLength_Right[block];

      // left off-diagonal block
      if( block != blocksStartReal )
      {
         assert(block >= 1);
         assert(linkStartBlockLength_Right[block - 1] >= 0);

         nnz += length_Left * length_Right;
      }

      // diagonal block
      nnz += length_Left * length_Right;

      // right off-diagonal block
      if( block != blockLast )
      {
         assert(block < nBlocks - 1);

         const int nextlength_Right = linkStartBlockLength_Right[block + 1];

         assert(nextlength_Right >= 0);

         nnz += length_Left * nextlength_Right;
      }
   }

   // dense (right or lower) border?
   if( nRowsSparse < int(nRows) || nColsSparse < int(nCols) )
   {
      assert(nRowsSparse <= int(nRows) && nColsSparse <= int(nCols));

      const int nRowsDense = int(nRows) - nRowsSparse;
      const int nColsDense = int(nCols) - nColsSparse;

      nnz += nRowsDense * nColsSparseRange;  // lower left border part (without right border)
      nnz += nRowsSparseRange * nColsDense;  // upper right border
      nnz += nRowsDense * nColsDense;        // lower right border

   }

   return nnz;
}


int sData::n2linksRows(const std::vector<int>& linkStartBlockLengths)
{
   int n = 0;

   for( size_t i = 0; i < linkStartBlockLengths.size(); ++i )
      n += linkStartBlockLengths[i];

   return n;
}

std::vector<int> sData::get2LinkLengthsVec(const std::vector<int>& linkStartBlockId, const size_t nBlocks)
{
   std::vector<int> linkStartBlockLengths(nBlocks, 0);

   const size_t nlinks = linkStartBlockId.size();

   for( size_t i = 0; i < nlinks; i++ )
   {
      const int block = linkStartBlockId[i];

      if( block >= 0 )
      {
         assert(size_t(block) < nBlocks);
         linkStartBlockLengths[block]++;
      }
   }
   assert(nBlocks == 0 || linkStartBlockLengths[nBlocks - 1] == 0);

   return linkStartBlockLengths;
}

SparseSymMatrix* sData::createSchurCompSymbSparseUpper()
{
   assert(children.size() > 0);

   const int nx0 = getLocalnx();
   const int my0 = getLocalmy();
   const int myl = getLocalmyl();
   const int mzl = getLocalmzl();
   const int sizeSC = nx0 + my0 + myl + mzl;
   const int nnz = getSchurCompMaxNnz();

   assert(nnz > 0);
   assert(myl >= 0 && mzl >= 0);

   int* krowM = new int[sizeSC + 1];
   int* jcolM = new int[nnz];
   double* M = new double[nnz];

   krowM[0] = 0;

   // get B_0^T (resp. A_0^T)
   SparseGenMatrix& Btrans = getLocalB().getTranspose();
   int* const startRowBtrans = Btrans.krowM();
   int* const colidxBtrans = Btrans.jcolM();

#ifndef NDEBUG
      int bm, bn;
      Btrans.getSize(bm, bn);
      assert(bm == nx0 && bn == my0);
#endif

   const int nx0NonZero = nx0 - n0LinkVars;
   int nnzcount = 0;

   assert(nx0NonZero >= 0);

   // dense square block, B_0^T, and dense border blocks todo: add space for CDCt
   for( int i = 0; i < nx0NonZero; ++i )
   {
      const int blength = startRowBtrans[i + 1] - startRowBtrans[i];
      assert(blength >= 0);

      krowM[i + 1] = krowM[i] + (nx0 - i) + blength + myl + mzl;

      appendRowDense(i, nx0, nnzcount, jcolM);

      appendRowSparse(startRowBtrans[i], startRowBtrans[i + 1], nx0, colidxBtrans, nnzcount, jcolM);

      appendRowDense(nx0 + my0, nx0 + my0 + myl + mzl, nnzcount, jcolM);

      assert(nnzcount == krowM[i + 1]);
   }

   // dense square block and rest of B_0, F_0^T, G_0^T
   for( int i = nx0NonZero; i < nx0; ++i )
   {
      appendRowDense(i, nx0, nnzcount, jcolM);

      appendRowSparse(startRowBtrans[i], startRowBtrans[i + 1], nx0, colidxBtrans, nnzcount, jcolM);

      if( myl > 0 )
      {
         SparseGenMatrix& Ft = getLocalF().getTranspose();
         const int* startRowFtrans = Ft.krowM();
         const int* colidxFtrans = Ft.jcolM();

         appendRowSparse(startRowFtrans[i], startRowFtrans[i + 1], nx0 + my0, colidxFtrans, nnzcount, jcolM);
      }

      if( mzl > 0 )
      {
         SparseGenMatrix& Gt = getLocalG().getTranspose();
         const int* startRowGtrans = Gt.krowM();
         const int* colidxGtrans = Gt.jcolM();

         appendRowSparse(startRowGtrans[i], startRowGtrans[i + 1], nx0 + my0 + myl, colidxGtrans, nnzcount, jcolM);
      }

      krowM[i + 1] = nnzcount;
   }

   // empty rows; put diagonal for PARDISO
   for( int i = nx0; i < nx0 + my0; ++i )
   {
      const int rowStartIdx = krowM[i];

      jcolM[rowStartIdx] = i;
      krowM[i + 1] = rowStartIdx + 1;
   }

   nnzcount += my0;

   // equality linking: sparse diagonal blocks, and mixed rows
   int blockStartrow = 0;
   const int n2linksRowsEq = n2linkRowsEq();
   const int bordersizeEq = linkStartBlockIdA.size() - n2linksRowsEq;
   const int borderstartEq = nx0 + my0 + n2linksRowsEq;
   const int n2linksRowsIneq = n2linkRowsIneq();
   const int bordersizeIneq = linkStartBlockIdC.size() - n2linksRowsIneq;
   const int borderstartIneq = nx0 + my0 + myl + n2linksRowsIneq;

   assert(bordersizeEq >= 0 && n2linksRowsEq <= myl);
   assert(bordersizeIneq >= 0 && n2linksRowsIneq <= mzl);

   for( int i = nx0 + my0, j = 0, colIdxOffset = 0, blockStartrowMix = 0; i < nx0 + my0 + myl; ++i, ++j )
   {
      int blockrownnz = appendDiagBlocks(linkStartBlockIdA, linkStartBlockLengthsA, borderstartEq, bordersizeEq, i, j,
            blockStartrow, nnzcount, jcolM);

      blockrownnz += appendMixedBlocks(linkStartBlockIdA, linkStartBlockIdC, linkStartBlockLengthsA, linkStartBlockLengthsC,
            (nx0 + my0 + myl), bordersizeIneq, j, colIdxOffset, blockStartrowMix, nnzcount, jcolM);

      assert(blockStartrowMix == blockStartrow);

      krowM[i + 1] = krowM[i] + blockrownnz;
   }

   // inequality linking: dense border block and sparse diagonal blocks
   blockStartrow = 0;

   for( int i = nx0 + my0 + myl, j = 0; i < nx0 + my0 + myl + mzl; ++i, ++j )
   {
       const int blockrownnz = appendDiagBlocks(linkStartBlockIdC, linkStartBlockLengthsC, borderstartIneq, bordersizeIneq, i, j, blockStartrow, nnzcount, jcolM);

       krowM[i + 1] = krowM[i] + blockrownnz;
   }

   assert(nnzcount == nnz);

   return (new SparseSymMatrix(sizeSC, nnz, krowM, jcolM, M, 1, false));
}


SparseSymMatrix* sData::createSchurCompSymbSparseUpperDist(int blocksStart, int blocksEnd)
{
   assert(children.size() > 0);

   const int nx0 = getLocalnx();
   const int my0 = getLocalmy();
   const int myl = getLocalmyl();
   const int mzl = getLocalmzl();
   const int mylLocal = myl - getSCdiagBlocksNRows(linkStartBlockLengthsA)
      + getSCdiagBlocksNRows(linkStartBlockLengthsA, blocksStart, blocksEnd);
   const int mzlLocal = mzl - getSCdiagBlocksNRows(linkStartBlockLengthsC)
      + getSCdiagBlocksNRows(linkStartBlockLengthsC, blocksStart, blocksEnd);
   const int sizeSC = nx0 + my0 + myl + mzl;
   const int nnz = getSchurCompMaxNnzDist(blocksStart, blocksEnd);

   assert(getSchurCompMaxNnzDist(0, linkStartBlockLengthsA.size()) == getSchurCompMaxNnz());
   assert(blocksStart >= 0 && blocksStart < blocksEnd);
   assert(nnz > 0);
   assert(myl >= 0 && mzl >= 0);

   int* const krowM = new int[sizeSC + 1];
   int* const jcolM = new int[nnz];
   double* const M = new double[nnz];

   krowM[0] = 0;

   // get B_0^T (resp. A_0^T)
   SparseGenMatrix& Btrans = getLocalB().getTranspose();
   int* const startRowBtrans = Btrans.krowM();
   int* const colidxBtrans = Btrans.jcolM();

#ifndef NDEBUG
      int bm, bn;
      Btrans.getSize(bm, bn);
      assert(bm == nx0 && bn == my0);
#endif

   const int nx0NonZero = nx0 - n0LinkVars;
   const int n2linksRowsEq = n2linkRowsEq();
   const int bordersizeEq = linkStartBlockIdA.size() - n2linksRowsEq;
   const int borderstartEq = nx0 + my0 + n2linksRowsEq;
   const int n2linksRowsIneq = n2linkRowsIneq();
   const int bordersizeIneq = linkStartBlockIdC.size() - n2linksRowsIneq;
   const int borderstartIneq = nx0 + my0 + myl + n2linksRowsIneq;
   int local2linksStartEq;
   int local2linksEndEq;
   int local2linksStartIneq;
   int local2linksEndIneq;

   this->getSCrangeMarkers(blocksStart, blocksEnd, local2linksStartEq, local2linksEndEq,
         local2linksStartIneq, local2linksEndIneq);

   assert(nx0NonZero >= 0);
   assert(bordersizeEq >= 0 && n2linksRowsEq <= myl);
   assert(bordersizeIneq >= 0 && n2linksRowsIneq <= mzl);
   assert(local2linksStartEq >= nx0 + my0 && local2linksEndEq <= borderstartEq);
   assert(local2linksStartIneq >= nx0 + my0 + myl && local2linksEndIneq <= borderstartIneq);

   int nnzcount = 0;

   // dense square block, B_0^T, and dense border blocks todo: add space for CDCt
   for( int i = 0; i < nx0NonZero; ++i )
   {
      const int blength = startRowBtrans[i + 1] - startRowBtrans[i];
      assert(blength >= 0);

      krowM[i + 1] = krowM[i] + (nx0 - i) + blength + mylLocal + mzlLocal;

      appendRowDense(i, nx0, nnzcount, jcolM);
      appendRowSparse(startRowBtrans[i], startRowBtrans[i + 1], nx0, colidxBtrans, nnzcount, jcolM);

      appendRowDense(local2linksStartEq, local2linksEndEq, nnzcount, jcolM);
      appendRowDense(borderstartEq, borderstartEq + bordersizeEq, nnzcount, jcolM);
      appendRowDense(local2linksStartIneq, local2linksEndIneq, nnzcount, jcolM);
      appendRowDense(borderstartIneq, borderstartIneq + bordersizeIneq, nnzcount, jcolM);

      assert(nnzcount == krowM[i + 1]);
   }

   // dense square block and rest of B_0, F_0^T, G_0^T
   for( int i = nx0NonZero; i < nx0; ++i )
   {
      appendRowDense(i, nx0, nnzcount, jcolM);

      appendRowSparse(startRowBtrans[i], startRowBtrans[i + 1], nx0, colidxBtrans, nnzcount, jcolM);

      if( myl > 0 )
      {
         SparseGenMatrix& Ft = getLocalF().getTranspose();
         const int* startRowFtrans = Ft.krowM();
         const int* colidxFtrans = Ft.jcolM();

         appendRowSparse(startRowFtrans[i], startRowFtrans[i + 1], nx0 + my0, colidxFtrans, nnzcount, jcolM);
      }

      if( mzl > 0 )
      {
         SparseGenMatrix& Gt = getLocalG().getTranspose();
         const int* startRowGtrans = Gt.krowM();
         const int* colidxGtrans = Gt.jcolM();

         appendRowSparse(startRowGtrans[i], startRowGtrans[i + 1], nx0 + my0 + myl, colidxGtrans, nnzcount, jcolM);
      }

      krowM[i + 1] = nnzcount;
   }

   // empty rows; put diagonal for PARDISO
   for( int i = nx0; i < nx0 + my0; ++i )
   {
      const int rowStartIdx = krowM[i];

      jcolM[rowStartIdx] = i;
      krowM[i + 1] = rowStartIdx + 1;
   }

   nnzcount += my0;

   // equality linking: sparse diagonal blocks, and mixed rows
   int blockStartrow = 0;

   for( int i = nx0 + my0, j = 0, colIdxOffset = 0, blockStartrowMix = 0; i < nx0 + my0 + myl; ++i, ++j )
   {
      appendDiagBlocksDist(linkStartBlockIdA, linkStartBlockLengthsA, borderstartEq, bordersizeEq, i, j,
            blocksStart, blocksEnd, blockStartrow, nnzcount, jcolM);

      appendMixedBlocksDist(linkStartBlockIdA, linkStartBlockIdC, linkStartBlockLengthsA, linkStartBlockLengthsC, (nx0 + my0 + myl),
            bordersizeIneq, j, blocksStart, blocksEnd, colIdxOffset, blockStartrowMix, nnzcount, jcolM);

      assert(blockStartrowMix == blockStartrow);

      krowM[i + 1] = nnzcount;
   }

   // inequality linking: dense border block and sparse diagonal blocks

   blockStartrow = 0;

   for( int i = nx0 + my0 + myl, j = 0; i < nx0 + my0 + myl + mzl; ++i, ++j )
   {
      appendDiagBlocksDist(linkStartBlockIdC, linkStartBlockLengthsC, borderstartIneq, bordersizeIneq, i, j,
            blocksStart, blocksEnd, blockStartrow, nnzcount, jcolM);

      krowM[i + 1] = nnzcount;
   }

   assert(nnzcount == nnz);

   this->initDistMarker(blocksStart, blocksEnd);

   return (new SparseSymMatrix(sizeSC, nnzcount, krowM, jcolM, M, 1, false));
}

std::vector<unsigned int> sData::get0VarsLastGlobalsFirstPermutation(std::vector<int>& link_vars_n_blocks, int& n_globals)
{
   const size_t n_link_vars = link_vars_n_blocks.size();
   n_globals = 0;

   if( n_link_vars == 0 )
      return std::vector<unsigned int>();

   std::vector<unsigned int> permvec(n_link_vars, 0);

   int count = 0;
   int back_count = n_link_vars - 1;

   for( size_t i = 0; i < n_link_vars; ++i )
   {
      assert( count <= back_count );
      assert( link_vars_n_blocks[i] >= 0 );

      if( link_vars_n_blocks[i] > threshold_global_vars )
      {
         ++n_globals;
         permvec[count++] = i;
      }
      else if( link_vars_n_blocks[i] == 0 )
         permvec[back_count--] = i;
   }

   for( size_t i = 0; i < n_link_vars; ++i )
   {
      if( link_vars_n_blocks[i] > 0 && link_vars_n_blocks[i] <= threshold_global_vars )
      {
         assert( count <= back_count );
         permvec[count++] = i;
      }
   }
   assert(count == back_count + 1);

   permuteVector(permvec, link_vars_n_blocks);

#ifndef NDEBUG
   int n_globals_copy = 0;
   int phase = 0;
   for( size_t i = 0; i < n_link_vars; ++i )
   {
      if( phase == 0 )
      {
         if( link_vars_n_blocks[i] <= threshold_global_vars )
         {
            ++phase;
            --i;
         }
         else
            ++n_globals_copy;
      }
      else if( phase == 1 )
      {
         if( link_vars_n_blocks[i] == 0 )
         {
            ++phase;
            --i;
         }
         else
         {
            assert( 0 < link_vars_n_blocks[i] );
            assert( link_vars_n_blocks[i] <= threshold_global_vars );
         }
      }
      else if( phase == 2 )
         assert( link_vars_n_blocks[i] == 0 );
   }
   assert( n_globals_copy == n_globals );
#endif

   return permvec;
}

std::vector<unsigned int> sData::getAscending2LinkFirstGlobalsLastPermutation(std::vector<int>& linkStartBlockId,
      std::vector<int>& n_blocks_per_row, size_t nBlocks, int& n_globals)
{
   assert( linkStartBlockId.size() == n_blocks_per_row.size() );
   const size_t n_links = linkStartBlockId.size();
   n_globals = 0;

   if( n_links == 0 )
      return std::vector<unsigned int>();

   std::vector<unsigned int> permvec(n_links, 0);
   std::vector<int> w(nBlocks + 1, 0);

   /* count the 2-links per block - the ones starting at block -1 are no 2-links and are counted in w[0] */
   for( size_t i = 0; i < n_links; ++i )
   {
      assert( -1 <= linkStartBlockId[i] );
      assert( linkStartBlockId[i] < int(nBlocks) );

      w[linkStartBlockId[i] + 1]++;
   }

   /* set w[i] to the amount of preceding 2-links, so the start of the 2-links starting in block i */
   size_t n_two_links = 0;
   for( size_t i = 1; i <= nBlocks; ++i )
   {
      n_two_links += w[i];
      w[i] = n_two_links;
   }

   assert( n_two_links + w[0] == n_links);
   w[0] = 0; // former n non 2-links

   /* sort 2-links ascending to front */
   for( size_t i = 0; i < n_links; ++i )
   {
      /* index of 2_link_start of nBlocks if not a 2 link */
      const int two_link_start = (linkStartBlockId[i] >= 0) ? linkStartBlockId[i] : int(nBlocks);

      assert(w[two_link_start] <= int(n_links));
      assert(permvec[w[two_link_start]] == 0);

      permvec[w[two_link_start]] = i;

      /* move start of i-2-link block */
      w[two_link_start]++;
   }

   /* permvec now moves 2-links ascending and the rest to the end */
   /* now permute global (long) linking constraints further to the end */
#ifndef NDEBUG
   for( size_t i = 1; i < n_two_links; i++ )
      assert( linkStartBlockId[permvec[i - 1]] <= linkStartBlockId[permvec[i]] );
   for( size_t i = n_two_links; i < n_links; ++i )
      assert( linkStartBlockId[permvec[i]] == -1 );
#endif

   /* got through non-2-links from front and back and swap all globals to the back */
   int front_pointer = n_two_links;
   assert( n_links > 0 );
   int end_pointer = n_links - 1;

   while( front_pointer <= end_pointer )
   {
      assert( linkStartBlockId[permvec[front_pointer]] == -1 );
      assert( linkStartBlockId[permvec[end_pointer]] == -1 );

      if( n_blocks_per_row[permvec[front_pointer]] <= threshold_global_cons )
         ++front_pointer;
      else if( n_blocks_per_row[permvec[end_pointer]] > threshold_global_cons )
         --end_pointer;
      else
      {
         assert( front_pointer < end_pointer );
         std::swap( permvec[end_pointer], permvec[front_pointer] );
         assert( n_blocks_per_row[permvec[front_pointer]] <= threshold_global_cons );
         assert( n_blocks_per_row[permvec[end_pointer]] > threshold_global_cons );

         ++front_pointer;
         --end_pointer;
      }
   }
   n_globals = n_links - front_pointer;

   assert( permutationIsValid(permvec) );

   permuteVector(permvec, n_blocks_per_row);
   permuteVector(permvec, linkStartBlockId);

#ifndef NDEBUG
   int phase = 0;
   int n_globals_copy = 0;
   for( size_t i = 0; i < linkStartBlockId.size(); ++i )
   {
      /* first ones are ascending 2-links */
      if( phase == 0 )
      {
         if( linkStartBlockId[i] == -1 )
         {
            ++phase;
            --i;
         }
         else
         {
            assert( n_blocks_per_row[i] == 2 );
            if( i > 1 )
               assert( linkStartBlockId[i - 1] <= linkStartBlockId[i] );
         }
      }
      /* n-links up to the threshold */
      else if( phase == 1 )
      {
         if( n_blocks_per_row[i] > threshold_global_cons )
         {
            ++phase;
            --i;
         }
         else
            assert( linkStartBlockId[i] == -1);
      }
      /* global linking constraints */
      else
      {
         ++n_globals_copy;
         assert( n_blocks_per_row[i] > threshold_global_cons );
      }
   }
   assert( n_globals == n_globals_copy );
#endif

   return permvec;
}

sData::sData(const sTree* tree_, OoqpVector * c_in, SymMatrix * Q_in,
        OoqpVector * xlow_in, OoqpVector * ixlow_in,
        OoqpVector * xupp_in, OoqpVector * ixupp_in,
        GenMatrix  * A_in, OoqpVector * bA_in,
        GenMatrix  * C_in,
        OoqpVector * clow_in, OoqpVector * iclow_in,
        OoqpVector * cupp_in, OoqpVector * icupp_in,
        bool add_children, bool is_hierarchy_root
        )
  : QpGenData(SparseLinearAlgebraPackage::soleInstance(),
         c_in, Q_in, xlow_in, ixlow_in, xupp_in, ixupp_in,
         A_in, bA_in, C_in, clow_in, iclow_in, cupp_in, icupp_in),
         is_hierarchy_root( is_hierarchy_root ),
         useLinkStructure( false )
{
  stochNode = tree_;

  if( add_children )
     createChildren();
  n0LinkVars = 0;
}

void sData::writeToStreamDense( std::ostream& out ) const
{
   const int myRank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   if( myRank == 0 ) out <<  "A: \n";
   (*A).writeToStreamDense(out);
   if( myRank == 0 ) out <<  "C: \n";
   (*C).writeToStreamDense(out);
   if( myRank == 0 ) out <<  "obj: \n";
   (*g).writeToStream(out);
   if( myRank == 0 ) out <<  "bA: \n";
   (*bA).writeToStream(out);
   if( myRank == 0 ) out <<  "xupp: \n";
   (*bux).writeToStream(out);
   if( myRank == 0 ) out <<  "ixupp: \n";
   (*ixupp).writeToStream(out);
   if( myRank == 0 ) out <<  "xlow: \n";
   (*blx).writeToStream(out);
   if( myRank == 0 ) out <<  "ixlow: \n";
   (*ixlow).writeToStream(out);
   if( myRank == 0 ) out <<  "cupp: \n";
   (*bu).writeToStream(out);
   if( myRank == 0 ) out <<  "icupp: \n";
   (*icupp).writeToStream(out);
   if( myRank == 0 ) out <<  "clow: \n";
   (*bl).writeToStream(out);
   if( myRank == 0 ) out <<  "iclow: \n";
   (*iclow).writeToStream(out);
}

/** Write the LP in MPS format. Only works if not distributed. */
void sData::writeMPSformat( std::ostream& out)
{
   // Note: only write the inequalities that have a finite rhs
   // (because no specified rhs of a row implies rhs=0).
   // Also, variable coefficients with indices in inequalitites with
   // inifnite rhs are not written because these rows do not appear in the MPs model.

   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   if( world_size > 1 )
   {
      std::cout << "MPS format writer only available using one Process!\n";
      return;
   }
   std::cout << "Writing MPS format...\n";

   out <<  "NAME PIPS_to_MPS\n";
   out << "ROWS\n";
   out << " N COST\n";

   // write all row names and if they are E, L or G
   (*A).writeMPSformatRows(out, 0, nullptr);
   (*C).writeMPSformatRows(out, 1, icupp);
   (*C).writeMPSformatRows(out, 2, iclow);

   // write all variable names
   out <<  "COLUMNS \n";
   writeMPSColumns(out);

   // write all rhs / lhs
   out <<  "RHS \n";

   (*bA).writeMPSformatRhs(out, 0, nullptr);
   (*bu).writeMPSformatRhs(out, 1, icupp);
   (*bl).writeMPSformatRhs(out, 2, iclow);

   // write all variable bounds
   out <<  "BOUNDS \n";
   (*bux).writeMPSformatBounds(out, ixupp, true);
   (*blx).writeMPSformatBounds(out, ixlow, false);

   out <<  "ENDATA \n";

   cout<<"Finished writing MPS format.\n";
}

void sData::writeMPSColumns(ostream& out)
{
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   assert( world_size == 1 );

   int n;
   string varName;
   string rowNameStub;
   string rowNameStubLT;
   string rowNameStubGT;
   StochVector& gStoch = dynamic_cast<StochVector&>(*g);
   StochVector& icuppStoch = dynamic_cast<StochVector&>(*icupp);
   StochVector& iclowStoch = dynamic_cast<StochVector&>(*iclow);
   StochGenMatrix& AStoch = dynamic_cast<StochGenMatrix&>(*A);
   SparseGenMatrix& ASparseTrans = dynamic_cast<SparseGenMatrix*>(AStoch.Bmat)->getTranspose();
   StochGenMatrix& CStoch = dynamic_cast<StochGenMatrix&>(*C);
   SparseGenMatrix& CSparseTrans = dynamic_cast<SparseGenMatrix*>(CStoch.Bmat)->getTranspose();

   SimpleVector* gSimple = dynamic_cast<SimpleVector*>(gStoch.vec);
   n = gSimple->length();

   std::stringstream sstmCol;
   std::stringstream sstmRow;


   // linking variables:
   for( int col = 0; col<n; col++ )
   {
      sstmCol.clear();
      sstmCol.str("");
      sstmCol << " var_L_" << col;
      varName = sstmCol.str();

      // cost coefficients:
      rowNameStub = "COST";
      if( gSimple->elements()[col] != 0 )
         out<<varName<< " " << rowNameStub << " " << gSimple->elements()[col] <<"\n";

      // coefficients in A_0:
      rowNameStub = "row_E_R_";
      for( int k = ASparseTrans.krowM()[col]; k<ASparseTrans.krowM()[col+1]; k++ )
         out<<varName<< " " << rowNameStub << ASparseTrans.jcolM()[k] << " " << ASparseTrans.M()[k] <<"\n";

      // coefficients in F_0:
      if( AStoch.Blmat )
      {
         SparseGenMatrix& ABlmatSparseTrans = dynamic_cast<SparseGenMatrix*>(AStoch.Blmat)->getTranspose();
         rowNameStub = "row_E_L_";
         for( int k = ABlmatSparseTrans.krowM()[col]; k<ABlmatSparseTrans.krowM()[col+1]; k++ )
            out<<varName<< " " << rowNameStub << ABlmatSparseTrans.jcolM()[k] << " " << ABlmatSparseTrans.M()[k] <<"\n";
         dynamic_cast<SparseGenMatrix*>(AStoch.Blmat)->deleteTransposed();
      }
      // coefficients in A_i:
      for( size_t it = 0; it < children.size(); it++ )
      {
         SparseGenMatrix& AChildSparseTrans = dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Amat)->getTranspose();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_E_"<<(int)it<<"_";
         rowNameStub = sstmRow.str();
         for( int k = AChildSparseTrans.krowM()[col]; k<AChildSparseTrans.krowM()[col+1]; k++ )
            out<<varName<< " " << rowNameStub << AChildSparseTrans.jcolM()[k] << " " << AChildSparseTrans.M()[k] <<"\n";
         dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Amat)->deleteTransposed();
      }

      // coefficients in C_0:
      rowNameStubLT = "row_L_R_";
      rowNameStubGT = "row_G_R_";
      for( int k = CSparseTrans.krowM()[col]; k<CSparseTrans.krowM()[col+1]; k++ )
      {
         int rowIdx = CSparseTrans.jcolM()[k];
         if( dynamic_cast<SimpleVector*>(icuppStoch.vec)->elements()[rowIdx] != 0.0)
            out<<varName<< " " << rowNameStubLT << rowIdx << " " << CSparseTrans.M()[k] <<"\n";
         if( dynamic_cast<SimpleVector*>(iclowStoch.vec)->elements()[rowIdx] != 0.0)
            out<<varName<< " " << rowNameStubGT << rowIdx << " " << CSparseTrans.M()[k] <<"\n";
      }
      // coefficients in G_0:
      if( CStoch.Blmat )
      {
         SparseGenMatrix& CBlmatSparseTrans = dynamic_cast<SparseGenMatrix*>(CStoch.Blmat)->getTranspose();
         rowNameStubLT = "row_L_L_";
         rowNameStubGT = "row_G_L_";
         for( int k = CBlmatSparseTrans.krowM()[col]; k<CBlmatSparseTrans.krowM()[col+1]; k++ )
         {
            int rowIdx = CBlmatSparseTrans.jcolM()[k];
            if( dynamic_cast<SimpleVector*>(icuppStoch.vecl)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubLT << rowIdx << " " << CBlmatSparseTrans.M()[k] <<"\n";
            if( dynamic_cast<SimpleVector*>(iclowStoch.vecl)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubGT << rowIdx << " " << CBlmatSparseTrans.M()[k] <<"\n";
         }
         dynamic_cast<SparseGenMatrix*>(CStoch.Blmat)->deleteTransposed();
      }
      // coefficients in C_i:
      for( size_t it = 0; it < children.size(); it++ )
      {
         SparseGenMatrix& CChildSparseTrans = dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Amat)->getTranspose();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_L_"<<(int)it<<"_";
         rowNameStubLT = sstmRow.str();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_G_"<<(int)it<<"_";
         rowNameStubGT = sstmRow.str();
         for( int k = CChildSparseTrans.krowM()[col]; k<CChildSparseTrans.krowM()[col+1]; k++ )
         {
            int rowIdx = CChildSparseTrans.jcolM()[k];
            if( dynamic_cast<SimpleVector*>(icuppStoch.children[it]->vec)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubLT << rowIdx << " " << CChildSparseTrans.M()[k] <<"\n";
            if( dynamic_cast<SimpleVector*>(iclowStoch.children[it]->vec)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubGT << rowIdx << " " << CChildSparseTrans.M()[k] <<"\n";
         }
         dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Amat)->deleteTransposed();
      }
   }

   // non-linking variables:
   for( size_t it = 0; it < children.size(); it++ )
   {
      SimpleVector* gSimple = dynamic_cast<SimpleVector*>(gStoch.children[it]->vec);
      n = gSimple->length();

      for( int col = 0; col<n; col++ )
      {
         sstmCol.clear();
         sstmCol.str("");
         sstmCol << " var_"<<(int)it <<"_" << col;
         varName = sstmCol.str();

         // coeffs in COST:
         rowNameStub = "COST";
         if( gSimple->elements()[col] != 0 )
            out<<varName<< " " << rowNameStub << " " << gSimple->elements()[col] <<"\n";

         // coeffs in A_i:
         SparseGenMatrix& AChildSparseTrans = dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Bmat)->getTranspose();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_E_"<<(int)it<<"_";
         rowNameStub = sstmRow.str();
         for( int k = AChildSparseTrans.krowM()[col]; k<AChildSparseTrans.krowM()[col+1]; k++ )
            out<<varName<< " " << rowNameStub << AChildSparseTrans.jcolM()[k] << " " << AChildSparseTrans.M()[k] <<"\n";
         dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Bmat)->deleteTransposed();

         // coefficients in D_i:
         SparseGenMatrix& CChildSparseTrans = dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Bmat)->getTranspose();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_L_"<<(int)it<<"_";
         rowNameStubLT = sstmRow.str();
         sstmRow.clear();
         sstmRow.str("");
         sstmRow << "row_G_"<<(int)it<<"_";
         rowNameStubGT = sstmRow.str();
         for( int k = CChildSparseTrans.krowM()[col]; k<CChildSparseTrans.krowM()[col+1]; k++ )
         {
            int rowIdx = CChildSparseTrans.jcolM()[k];
            if( dynamic_cast<SimpleVector*>(icuppStoch.children[it]->vec)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubLT << CChildSparseTrans.jcolM()[k] << " " << CChildSparseTrans.M()[k] <<"\n";
            if( dynamic_cast<SimpleVector*>(iclowStoch.children[it]->vec)->elements()[rowIdx] != 0.0)
               out<<varName<< " " << rowNameStubGT << CChildSparseTrans.jcolM()[k] << " " << CChildSparseTrans.M()[k] <<"\n";
         }
         dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Bmat)->deleteTransposed();

         // coefficients in F_i:
         if( dynamic_cast<StochGenMatrix*>(AStoch.children[it])->Blmat )
         {
            SparseGenMatrix& ABlmatSparseTrans = dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Blmat)->getTranspose();
            rowNameStub = "row_E_L_";
            for( int k = ABlmatSparseTrans.krowM()[col]; k<ABlmatSparseTrans.krowM()[col+1]; k++ )
               out<<varName<< " " << rowNameStub << ABlmatSparseTrans.jcolM()[k] << " " << ABlmatSparseTrans.M()[k] <<"\n";
            dynamic_cast<SparseGenMatrix*>(AStoch.children[it]->Blmat)->deleteTransposed();
         }

         // coefficients in G_i:
         if( dynamic_cast<StochGenMatrix*>(CStoch.children[it])->Blmat )
         {
            SparseGenMatrix& CBlmatSparseTrans = dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Blmat)->getTranspose();
            rowNameStubLT = "row_L_L_";
            rowNameStubGT = "row_G_L_";
            for( int k = CBlmatSparseTrans.krowM()[col]; k<CBlmatSparseTrans.krowM()[col+1]; k++ )
            {
               int rowIdx = CBlmatSparseTrans.jcolM()[k];
               if( dynamic_cast<SimpleVector*>(icuppStoch.vecl)->elements()[rowIdx] != 0.0)
                  out<<varName<< " " << rowNameStubLT << rowIdx << " " << CBlmatSparseTrans.M()[k] <<"\n";
               if( dynamic_cast<SimpleVector*>(iclowStoch.vecl)->elements()[rowIdx] != 0.0)
                  out<<varName<< " " << rowNameStubGT << rowIdx << " " << CBlmatSparseTrans.M()[k] <<"\n";
            }
            dynamic_cast<SparseGenMatrix*>(CStoch.children[it]->Blmat)->deleteTransposed();
         }
      }
   }

   // delete transposed matrices:
   dynamic_cast<SparseGenMatrix*>(AStoch.Bmat)->deleteTransposed();
   dynamic_cast<SparseGenMatrix*>(CStoch.Bmat)->deleteTransposed();

}

sData* sData::cloneFull(bool switchToDynamicStorage) const
{
   // todo Q is empty!
   StochSymMatrixHandle Q_clone(dynamic_cast<const StochSymMatrix&>(*Q).clone());
   StochGenMatrixHandle A_clone(dynamic_cast<const StochGenMatrix&>(*A).cloneFull(switchToDynamicStorage));
   StochGenMatrixHandle C_clone(dynamic_cast<const StochGenMatrix&>(*C).cloneFull(switchToDynamicStorage));

   StochVectorHandle c_clone (dynamic_cast<StochVector*>(g->cloneFull()));
   StochVectorHandle bA_clone ( dynamic_cast<StochVector*>(bA->cloneFull()));
   StochVectorHandle xupp_clone (dynamic_cast<StochVector*>(bux->cloneFull()));
   StochVectorHandle ixupp_clone (dynamic_cast<StochVector*>(ixupp->cloneFull()));
   StochVectorHandle xlow_clone (dynamic_cast<StochVector*>(blx->cloneFull()));
   StochVectorHandle ixlow_clone (dynamic_cast<StochVector*>(ixlow->cloneFull()));
   StochVectorHandle cupp_clone (dynamic_cast<StochVector*>(bu->cloneFull()));
   StochVectorHandle icupp_clone (dynamic_cast<StochVector*>(icupp->cloneFull()));
   StochVectorHandle clow_clone (dynamic_cast<StochVector*>(bl->cloneFull()));
   StochVectorHandle iclow_clone (dynamic_cast<StochVector*>(iclow->cloneFull()));

   const sTree* tree_clone = stochNode;

   // TODO : proper copy ctor..
   sData* clone = new sData(tree_clone, c_clone, Q_clone, xlow_clone,
         ixlow_clone, xupp_clone, ixupp_clone, A_clone, bA_clone,
         C_clone, clow_clone, iclow_clone, cupp_clone, icupp_clone );

   return clone;
}

void
sData::createChildren()
{
  //follow the structure of one of the tree objects and create the same
  //structure for this class, and link this object with the corresponding 
  //vectors and matrices
  StochVector& gSt     = dynamic_cast<StochVector&>(*g);
  StochSymMatrix& QSt  = dynamic_cast<StochSymMatrix&>(*Q);
  
  StochVector& xlowSt  = dynamic_cast<StochVector&>(*blx); 
  StochVector& ixlowSt = dynamic_cast<StochVector&>(*ixlow); 
  StochVector& xuppSt  = dynamic_cast<StochVector&>(*bux); 
  StochVector& ixuppSt = dynamic_cast<StochVector&>(*ixupp);
  StochGenMatrix& ASt  = dynamic_cast<StochGenMatrix&>(*A); 
  StochVector& bASt    = dynamic_cast<StochVector&>(*bA);
  StochGenMatrix& CSt  = dynamic_cast<StochGenMatrix&>(*C);
  StochVector& clowSt  = dynamic_cast<StochVector&>(*bl); 
  StochVector& iclowSt = dynamic_cast<StochVector&>(*iclow);
  StochVector& cuppSt  = dynamic_cast<StochVector&>(*bu); 
  StochVector& icuppSt = dynamic_cast<StochVector&>(*icupp); 
  
  for(size_t it=0; it<gSt.children.size(); it++) {
    AddChild(new sData(stochNode->getChildren()[it],
	       gSt.children[it], QSt.children[it],
	       xlowSt.children[it], ixlowSt.children[it],
	       xuppSt.children[it], ixuppSt.children[it],
	       ASt.children[it], bASt.children[it],
	       CSt.children[it],
	       clowSt.children[it], iclowSt.children[it],
	       cuppSt.children[it], icuppSt.children[it] )
    );
  }
}

void sData::destroyChildren()
{
   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->destroyChildren();
      delete children[it];
   }
   children.clear();
}

void sData::getLinkConsSplitPermutations(std::vector<unsigned int>& perm_A, std::vector<unsigned int>& perm_C)
{
   perm_A.resize( linkConsPermutationA.size() );
   perm_C.resize( linkConsPermutationC.size() );

   assert( linkStartBlockIdA.size() == linkConsPermutationA.size() );
   assert( linkStartBlockIdC.size() == linkConsPermutationC.size() );

   assert( stochNode->myl() >= 0 );
   assert( stochNode->mzl() >= 0 );
   const int n_global_links_eq_after_split = stochNode->myl();
   const int n_global_links_ineq_after_split = stochNode->mzl();

   assert( linkStartBlockIdA.size() >= static_cast<size_t>(n_global_links_eq_after_split + n_global_eq_linking_conss) );
   assert( linkStartBlockIdC.size() >= static_cast<size_t>(n_global_links_ineq_after_split + n_global_ineq_linking_conss) );

#ifndef NDEBUG
   int last_map_block = -1;
#endif
   const std::vector<unsigned int>& map_block_subcomm = dynamic_cast<const sTreeCallbacks*>(stochNode)->getMapBlockSubTrees();
   assert( false && "TODO: check..");
   {
      size_t pos_local_twolinks_A = 0;
      const size_t end_local_twolinks_A = static_cast<size_t>(linkStartBlockIdA.size() - n_global_links_eq_after_split - n_global_eq_linking_conss);
      size_t pos_global_linksA = end_local_twolinks_A;

      for( size_t i = 0; i < linkStartBlockIdA.size(); ++i )
      {
         assert( last_map_block <= linkStartBlockIdA[i] );

         if( linkStartBlockIdA[i] == - 1 )
         {
            assert( pos_local_twolinks_A == end_local_twolinks_A );
            perm_A[pos_global_linksA] = i;
            ++pos_global_linksA;
         }
         else
         {
            assert( 0 <= linkStartBlockIdA[i] );
            const size_t start_block_link_i = static_cast<size_t>(linkStartBlockIdA[i]);

            assert( start_block_link_i < map_block_subcomm.size() );

            if( start_block_link_i == map_block_subcomm.size() - 1 )
               perm_A[pos_local_twolinks_A] = i;
            else if( map_block_subcomm[start_block_link_i] != map_block_subcomm[start_block_link_i + 1] )
            {
               perm_A[pos_global_linksA] = i;
               ++pos_global_linksA;
            }
            else
            {
               assert( map_block_subcomm[start_block_link_i] == map_block_subcomm[start_block_link_i + 1] );
               perm_A[pos_local_twolinks_A] = i;
               ++pos_local_twolinks_A;
            }
         }
      }
      assert( pos_local_twolinks_A == end_local_twolinks_A );
      assert( pos_global_linksA == linkStartBlockIdA.size() );
   }

   {
      size_t pos_local_twolinks_C = 0;
      const size_t end_local_twolinks_C = static_cast<size_t>(linkStartBlockIdC.size() - n_global_links_ineq_after_split - n_global_ineq_linking_conss);
      size_t pos_global_linksC = end_local_twolinks_C;

      for( size_t i = 0; i < linkStartBlockIdC.size(); ++i )
      {
         assert( last_map_block <= linkStartBlockIdC[i] );

         if( linkStartBlockIdC[i] == - 1 )
         {
            assert( pos_local_twolinks_C == end_local_twolinks_C );
            perm_C[pos_global_linksC] = i;
            ++pos_global_linksC;
         }
         else
         {
            assert( 0 <= linkStartBlockIdC[i] );
            const size_t start_block_link_i = static_cast<size_t>(linkStartBlockIdC[i]);
            assert( start_block_link_i < map_block_subcomm.size() );

            if( start_block_link_i == map_block_subcomm.size() - 1 )
               perm_C[pos_local_twolinks_C] = i;
            else if( map_block_subcomm[start_block_link_i] != map_block_subcomm[start_block_link_i + 1] )
            {
               perm_C[pos_global_linksC] = i;
               ++pos_global_linksC;
            }
            else
            {
               assert( map_block_subcomm[start_block_link_i] == map_block_subcomm[start_block_link_i + 1] );
               perm_C[pos_local_twolinks_C] = i;
               ++pos_local_twolinks_C;
            }
         }
      }
      assert( pos_local_twolinks_C == end_local_twolinks_C );
      assert( pos_global_linksC == linkStartBlockIdC.size() );
   }
}

void sData::reorderLinkingConstraintsAccordingToSplit()
{
   /* assert that distributed Schur complement has not yet been initialized */
   assert( isSCrowLocal.size() == 0 );
   assert( isSCrowMyLocal.size() == 0 );

   // TODO assert tree compatibility (one split at exactly this node)

   std::vector<unsigned int> perm_A;
   std::vector<unsigned int> perm_C;
   getLinkConsSplitPermutations(perm_A, perm_C);

   assert( permutationIsValid(perm_A) );
   assert( permutationIsValid(perm_C) );

   /* which blocks do the individual two-links start in */
   permuteLinkingCons(perm_A, perm_C);
   permuteVector(perm_A, linkStartBlockIdA);
   permuteVector(perm_A, n_blocks_per_link_row_A);

   permuteVector(perm_C, linkStartBlockIdC);
   permuteVector(perm_C, n_blocks_per_link_row_C);

   permuteVector(perm_A, linkConsPermutationA);
   permuteVector(perm_C, linkConsPermutationC);
}

void sData::splitDataAccordingToTree()
{
   assert( stochNode->isHierarchicalInner() );
   reorderLinkingConstraintsAccordingToSplit();

   splitDataAccordingToTree();

   assert( false && "TODO: implement");
}

sData* sData::switchToHierarchicalData( const sTree* tree )
{
   assert( tree->isHierarchicalRoot() );
   assert( tree->nChildren() == 1 );

//   this->splitDataAccordingToTree( dynamic_cast<const sTreeCallbacks&>(*tree->children[0]) );

   // TODO : split lower levels according to tree structure
   // TODO : make method for top layer and split operation

   const int my_rank = PIPS_MPIgetRank();
   if( my_rank == 0 )
   {
      std::cout << "Switching to hierarchical data ..." << "\n";
      std::cout << "Trimming " << n_global_linking_vars << " vars, " << n_global_eq_linking_conss << " dense equalities, and " <<
            n_global_ineq_linking_conss << " inequalities for the border" << "\n";
   }

   SymMatrixHandle Q_hier( dynamic_cast<StochSymMatrix&>(*Q).raiseBorder(n_global_linking_vars) );

   GenMatrixHandle A_hier( dynamic_cast<StochGenMatrix&>(*A).raiseBorder(n_global_eq_linking_conss, n_global_linking_vars) );
   GenMatrixHandle C_hier( dynamic_cast<StochGenMatrix&>(*C).raiseBorder(n_global_ineq_linking_conss, n_global_linking_vars) );

   /* we ordered global linking vars first and global linking rows to the end */
   StochVectorHandle g_hier( dynamic_cast<StochVector&>(*g).raiseBorder(n_global_linking_vars, false, true) );
   StochVectorHandle bux_hier( dynamic_cast<StochVector&>(*bux).raiseBorder(n_global_linking_vars, false, true) );
   StochVectorHandle ixupp_hier( dynamic_cast<StochVector&>(*ixupp).raiseBorder(n_global_linking_vars, false, true) );
   StochVectorHandle blx_hier( dynamic_cast<StochVector&>(*blx).raiseBorder(n_global_linking_vars, false, true) );
   StochVectorHandle ixlow_hier( dynamic_cast<StochVector&>(*ixlow).raiseBorder(n_global_linking_vars, false, true) );

   StochVectorHandle bA_hier( dynamic_cast<StochVector&>(*bA).raiseBorder(n_global_eq_linking_conss, true, false) );

   StochVectorHandle bu_hier( dynamic_cast<StochVector&>(*bu).raiseBorder(n_global_ineq_linking_conss, true, false) );
   StochVectorHandle icupp_hier( dynamic_cast<StochVector&>(*icupp).raiseBorder(n_global_ineq_linking_conss, true, false) );
   StochVectorHandle bl_hier( dynamic_cast<StochVector&>(*bl).raiseBorder(n_global_ineq_linking_conss, true, false) );
   StochVectorHandle iclow_hier( dynamic_cast<StochVector&>(*iclow).raiseBorder(n_global_ineq_linking_conss, true, false) );

   // TODO what is this?
   //StochVector* sc_hier = dynamic_cast<StochVector&>(*sc).shaveBorder(-1);

   sData* hierarchical_top = new sData(tree, g_hier.ptr_unsave(), Q_hier.ptr_unsave(), blx_hier.ptr_unsave(),
         ixlow_hier.ptr_unsave(), bux_hier.ptr_unsave(), ixupp_hier.ptr_unsave(),
         A_hier.ptr_unsave(), bA_hier.ptr_unsave(), C_hier.ptr_unsave(), bl_hier.ptr_unsave(),
         iclow_hier.ptr_unsave(), bu_hier.ptr_unsave(), icupp_hier.ptr_unsave(),
         false, true);

   assert( ixlow_hier->vec );
   assert( ixupp_hier->vec );
   this->nxlow -= ixlow_hier->vec->numberOfNonzeros();
   this->nxupp -= ixupp_hier->vec->numberOfNonzeros();

   assert( iclow_hier->vecl );
   assert( icupp_hier->vecl );
   this->mclow -= iclow_hier->vecl->numberOfNonzeros();
   this->mcupp -= icupp_hier->vecl->numberOfNonzeros();

   long long dummy;
   nx = g->length();
   A->getSize( my, dummy );
   C->getSize( mz, dummy );

   /* adapt vectors and global link sizes - we pushed these up */

   /* global linking variables have been ordered to the front, global linking constraints to the end of the matrices */
   /* linking vars */
   hierarchical_top->n_blocks_per_link_var = this->n_blocks_per_link_var;
   this->n_blocks_per_link_var.erase(n_blocks_per_link_var.begin(), n_blocks_per_link_var.begin() + n_global_linking_vars);

   /* Amat linking cons */
   hierarchical_top->linkStartBlockIdA = this->linkStartBlockIdA;
   this->linkStartBlockIdA.erase(linkStartBlockIdA.end() - n_global_eq_linking_conss, linkStartBlockIdA.end() );
   hierarchical_top->n_blocks_per_link_row_A = this->n_blocks_per_link_row_A;
   this->n_blocks_per_link_row_A.erase(n_blocks_per_link_row_A.end() - n_global_eq_linking_conss, n_blocks_per_link_row_A.end());

   /* Cmat linking cons */
   hierarchical_top->linkStartBlockIdC = this->linkStartBlockIdC;
   this->linkStartBlockIdC.erase(linkStartBlockIdC.end() - n_global_ineq_linking_conss, linkStartBlockIdC.end() );
   hierarchical_top->n_blocks_per_link_row_C = this->n_blocks_per_link_row_C;
   this->n_blocks_per_link_row_C.erase(n_blocks_per_link_row_C.end() - n_global_ineq_linking_conss, n_blocks_per_link_row_C.end() );

   assert( isSCrowLocal.size() == 0 );
   assert( isSCrowMyLocal.size() == 0 );

   hierarchical_top->n_global_eq_linking_conss = n_global_eq_linking_conss;
   this->n_global_eq_linking_conss = 0;

   hierarchical_top->n_global_ineq_linking_conss = n_global_ineq_linking_conss;
   this->n_global_ineq_linking_conss = 0;

   hierarchical_top->n_global_linking_vars = n_global_linking_vars;
   this->n_global_linking_vars = 0;

   hierarchical_top->useLinkStructure = false;

   hierarchical_top->children.push_back(this);
   stochNode = tree->getChildren()[0];

   // TODO: implement recursive layering of linear system
   //   this->splitIntoMultiple();

   if( my_rank == 0 )
      std::cout << "Hierarchical data built" << "\n";

   return hierarchical_top;
}

void sData::permuteLinkingCons(const std::vector<unsigned int>& permA, const std::vector<unsigned int>& permC)
{
   assert( permutationIsValid(permA) );
   assert( permutationIsValid(permC) );
   assert( !is_hierarchy_root );

   dynamic_cast<StochGenMatrix&>(*A).permuteLinkingCons(permA);
   dynamic_cast<StochGenMatrix&>(*C).permuteLinkingCons(permC);
   dynamic_cast<StochVector&>(*bA).permuteLinkingEntries(permA);
   dynamic_cast<StochVector&>(*bl).permuteLinkingEntries(permC);
   dynamic_cast<StochVector&>(*bu).permuteLinkingEntries(permC);
   dynamic_cast<StochVector&>(*iclow).permuteLinkingEntries(permC);
   dynamic_cast<StochVector&>(*icupp).permuteLinkingEntries(permC);
}

void sData::permuteLinkingVars(const std::vector<unsigned int>& perm)
{
   assert( permutationIsValid(linkVarsPermutation) );
   assert( !is_hierarchy_root );

   dynamic_cast<StochGenMatrix&>(*A).permuteLinkingVars(perm);
   dynamic_cast<StochGenMatrix&>(*C).permuteLinkingVars(perm);
   dynamic_cast<StochVector&>(*g).permuteVec0Entries(perm);
   dynamic_cast<StochVector&>(*bux).permuteVec0Entries(perm);
   dynamic_cast<StochVector&>(*blx).permuteVec0Entries(perm);
   dynamic_cast<StochVector&>(*ixupp).permuteVec0Entries(perm);
   dynamic_cast<StochVector&>(*ixlow).permuteVec0Entries(perm);
}

sVars* sData::getVarsUnperm(const sVars& vars, const sData& unpermData) const
{
   sVars* unperm_vars = new sVars(vars);

   if( is_hierarchy_root )
      unperm_vars->collapseHierarchicalStructure( unpermData.stochNode, unpermData.ixlow, unpermData.ixupp, unpermData.iclow, unpermData.icupp );

   assert( unperm_vars->children.size() == unpermData.children.size() );

   const std::vector<unsigned int> perm_inv_link_vars = getLinkVarsPermInv();
   const std::vector<unsigned int> perm_inv_link_cons_eq = getLinkConsEqPermInv();
   const std::vector<unsigned int> perm_inv_link_cons_ineq = getLinkConsIneqPermInv();

   if( perm_inv_link_vars.size() != 0 )
      unperm_vars->permuteVec0Entries( perm_inv_link_vars, true );

   if( perm_inv_link_cons_eq.size() != 0 )
      unperm_vars->permuteEqLinkingEntries( perm_inv_link_cons_eq );

   if( perm_inv_link_cons_ineq.size() != 0 )
      unperm_vars->permuteIneqLinkingEntries( perm_inv_link_cons_ineq, true );

   return unperm_vars;
}

sResiduals* sData::getResidsUnperm(const sResiduals& resids, const sData& unpermData) const
{
   sResiduals* unperm_resids = new sResiduals(resids);

   if( is_hierarchy_root )
      unperm_resids->collapseHierarchicalStructure( unpermData.ixlow, unpermData.ixupp, unpermData.iclow, unpermData.icupp );

   assert( unperm_resids->children.size() == unpermData.children.size() );

   const std::vector<unsigned int> perm_inv_link_vars = this->getLinkVarsPermInv();
   const std::vector<unsigned int> perm_inv_link_cons_eq = this->getLinkConsEqPermInv();
   const std::vector<unsigned int> perm_inv_link_cons_ineq = this->getLinkConsIneqPermInv();

   /* when using the hierarchical approach the unpermute is done in collapsHierarchicalStructure already */
   const bool do_not_permut_bounds = is_hierarchy_root ? true : false;

   if( perm_inv_link_vars.size() != 0 )
      unperm_resids->permuteVec0Entries( perm_inv_link_vars, do_not_permut_bounds );

   if( perm_inv_link_cons_eq.size() != 0 )
      unperm_resids->permuteEqLinkingEntries( perm_inv_link_cons_eq );

   if( perm_inv_link_cons_ineq.size() != 0 )
      unperm_resids->permuteIneqLinkingEntries( perm_inv_link_cons_ineq, do_not_permut_bounds );

   return unperm_resids;
}

void sData::activateLinkStructureExploitation()
{
   assert( !stochNode->isHierarchicalRoot() );

   if( useLinkStructure )
      return;
   useLinkStructure = true;

   const int myrank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   /* don't attempt to use linking structure when there actually is no linking constraints */
   if( stochNode->myl() == 0 && stochNode->mzl() == 0 )
   {
      if( pips_options::getBoolParameter( "HIERARCHICAL" ) )
      {
         if( myrank == 0 )
            std::cout << "No linking constraints found - hierarchical approach cannot be used" << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      }

      useLinkStructure = false;
      if( myrank == 0 )
         std::cout << "no linking constraints so no linking structure found\n";
      return;
   }

   const int nx0 = getLocalnx();

   const StochGenMatrix& Astoch = dynamic_cast<const StochGenMatrix&>(*A);
   const StochGenMatrix& Cstoch = dynamic_cast<const StochGenMatrix&>(*C);

   n_blocks_per_link_var = std::vector<int>(nx0, 0);
   Astoch.updateKLinkVarsCount( n_blocks_per_link_var );

   std::vector<int> tmp = std::vector<int>(nx0, 0); // to avoid doubling through second Allreduce inside the count functions
   Cstoch.updateKLinkVarsCount( tmp );

   std::transform( n_blocks_per_link_var.begin(), n_blocks_per_link_var.end(), tmp.begin(), n_blocks_per_link_var.begin(), std::plus<int>() );

   Astoch.get2LinkStartBlocksAndCountsNew(linkStartBlockIdA, n_blocks_per_link_row_A);
   Cstoch.get2LinkStartBlocksAndCountsNew(linkStartBlockIdC, n_blocks_per_link_row_C);

#ifndef NDEBUG
   std::vector<int> linkStart_A2 = Astoch.get2LinkStartBlocks();
   assert( linkStartBlockIdA.size() == linkStart_A2.size() );
   for( size_t i = 0; i < linkStart_A2.size(); ++i )
   {
      assert( linkStart_A2[i] == linkStartBlockIdA[i] );
      if( linkStart_A2[i] != linkStartBlockIdA[i] && myrank == 0)
         std::cout << "New : " << linkStart_A2[i] << " != " << linkStartBlockIdA[i] << " old" << "\n";
   }

   std::vector<int> linkStart_C2 = Cstoch.get2LinkStartBlocks();
   assert( linkStartBlockIdC.size() == linkStart_C2.size() );
   for( size_t i = 0; i < linkStart_C2.size(); ++i )
   {
      assert( linkStart_C2[i] == linkStartBlockIdC[i] );
      if( linkStart_C2[i] != linkStartBlockIdC[i] && myrank == 0)
         std::cout << "New : " << linkStart_C2[i] << " != " << linkStartBlockIdC[i] << " old" << "\n";
   }
#endif

   linkStartBlockLengthsA = get2LinkLengthsVec(linkStartBlockIdA, stochNode->nChildren());
   linkStartBlockLengthsC = get2LinkLengthsVec(linkStartBlockIdC, stochNode->nChildren());

   printLinkConsStats();
   printLinkVarsStats();

   int n2LinksEq = 0;
   int n2LinksIneq = 0;

   n0LinkVars = std::count_if(n_blocks_per_link_var.begin(), n_blocks_per_link_var.end(), [](int blocks){
      return (blocks == 0);
   } );

   n2LinksEq = std::count_if(linkStartBlockIdA.begin(), linkStartBlockIdA.end(), [](int blocks){
      return (blocks >= 0);
   } );

   n2LinksIneq = std::count_if(linkStartBlockIdC.begin(), linkStartBlockIdC.end(), [](int blocks){
      return (blocks >= 0);
   } );

#ifndef NDEBUG
   int n0LinkVars_cpy = 0;
   for( size_t i = 0; i < n_blocks_per_link_var.size(); ++i )
      if( n_blocks_per_link_var[i] == 0 )
         n0LinkVars_cpy++;

   int n2LinksEq_cpy = 0;
   for( size_t i = 0; i < linkStartBlockIdA.size(); ++i )
      if( linkStartBlockIdA[i] >= 0 )
         n2LinksEq_cpy++;

   int n2LinksIneq_cpy = 0;
   for( size_t i = 0; i < linkStartBlockIdC.size(); ++i )
      if( linkStartBlockIdC[i] >= 0 )
         n2LinksIneq_cpy++;

   assert( n0LinkVars_cpy == n0LinkVars );
   assert( n2LinksEq_cpy == n2LinksEq );
   assert( n2LinksIneq_cpy == n2LinksIneq );
   assert(n2LinksEq == n2linkRowsEq());
   assert(n2LinksIneq == n2linkRowsIneq());
#endif

   const double ratio = (n2LinksEq + n2LinksIneq + n0LinkVars) / double(linkStartBlockIdA.size() + linkStartBlockIdC.size() + n_blocks_per_link_var.size());
   if( myrank == 0 )
   {
      std::cout << "number of 0-link variables: " << n0LinkVars << " (out of "
            << nx0 << " link variables)\n";
      std::cout << "number of equality 2-links: " << n2LinksEq << " (out of "
            << linkStartBlockIdA.size() << " equalities)\n";
      std::cout << "number of inequality 2-links: " << n2LinksIneq << " (out of "
            << linkStartBlockIdC.size() << " equalities)\n";

      std::cout << "ratio: " << ratio << "\n";
   }


   if( !pips_options::getBoolParameter( "HIERARCHICAL" ) )
   {
      if( ratio < minStructuredLinksRatio )
      {
         if( myrank == 0 )
            std::cout << "not enough linking structure found ( required ratio : " << minStructuredLinksRatio << "\n";
         useLinkStructure = false;
      }
   }

   if( useLinkStructure )
   {
      assert(linkStartBlockIdA.size() == unsigned(stochNode->myl()));
      assert(linkStartBlockIdC.size() == unsigned(stochNode->mzl()));

   #ifndef NDEBUG
      const int myl = stochNode->myl();
      const int mzl = stochNode->mzl();
      assert(myl >= 0 && mzl >= 0 && (mzl + myl > 0));
   #endif

      assert( linkConsPermutationA.size() == 0 );
      assert( linkConsPermutationC.size() == 0 );

      const size_t nBlocks = dynamic_cast<StochVector&>(*g).children.size();

      // compute permutation vectors
      linkConsPermutationA = getAscending2LinkFirstGlobalsLastPermutation(linkStartBlockIdA, n_blocks_per_link_row_A, nBlocks, n_global_eq_linking_conss);
      linkConsPermutationC = getAscending2LinkFirstGlobalsLastPermutation(linkStartBlockIdC, n_blocks_per_link_row_C, nBlocks, n_global_ineq_linking_conss);
      permuteLinkingCons(linkConsPermutationA, linkConsPermutationC);

      assert(linkVarsPermutation.size() == 0);
      linkVarsPermutation = get0VarsLastGlobalsFirstPermutation(n_blocks_per_link_var, n_global_linking_vars);
      permuteLinkingVars(linkVarsPermutation);
   }
}

void sData::AddChild(sData* child)
{
   children.push_back(child);
}

double sData::objectiveValue(const QpGenVars * vars) const
{
   const StochVector& x = dynamic_cast<const StochVector&>(*vars->x);
   OoqpVectorHandle temp(x.clone());

   this->getg(*temp);
   this->Qmult(1.0, *temp, 0.5, *vars->x);

   return temp->dotProductWith(*vars->x);
}

void
sData::createScaleFromQ()
{

   assert("Not implemented!" && 0);

   // Stuff the diagonal elements of Q into the vector "sc"
   this->getDiagonalOfQ(*sc);

   // Modifying scVector is equivalent to modifying sc
   /*SimpleVector & scVector = dynamic_cast<SimpleVector &>(*sc);

    int scLength = scVector.length();

    for( int i = 0; i < scLength; i++){
    if( scVector[i] > 1)
    scVector[i] = 1.0/sqrt( scVector[i]);
    else
    scVector[i] = 1.0;
    }
    */
}

void sData::printLinkVarsStats()
{
   int n = getLocalnx();

   std::vector<int> linkCountA(n, 0);
   std::vector<int> linkCountC(n, 0);
   std::vector<int> linkCount0(n, 0);
   std::vector<int> linkCountLC(n, 0);

   StochGenMatrix& Astoch = dynamic_cast<StochGenMatrix&>(*A);
   StochGenMatrix& Cstoch = dynamic_cast<StochGenMatrix&>(*C);

   Astoch.updateKLinkVarsCount(linkCountA);
   Cstoch.updateKLinkVarsCount(linkCountC);

   Astoch.Bmat->getTranspose().updateNonEmptyRowsCount(linkCount0);
   Astoch.Bmat->deleteTransposed();
   Cstoch.Bmat->getTranspose().updateNonEmptyRowsCount(linkCount0);
   Cstoch.Bmat->deleteTransposed();

   if( Astoch.Blmat )
   {
      Astoch.Blmat->getTranspose().updateNonEmptyRowsCount(linkCountLC);
      Astoch.Blmat->deleteTransposed();
   }

   if( Cstoch.Blmat )
   {
      Cstoch.Blmat->getTranspose().updateNonEmptyRowsCount(linkCountLC);
      Cstoch.Blmat->deleteTransposed();
   }

   const int rank = PIPS_MPIgetRank();

   if( rank == 0 )
   {
      std::vector<int> linkSizes(nLinkStats, 0);

      int count0 = 0;
      int countLC = 0;
      int count0LC = 0;

      for( int i = 0; i < n; i++ )
      {
         const int linkCountAB = linkCountA[i] + linkCountC[i];
         assert(linkCountAB >= 0 && linkCount0[i] >= 0 && linkCountLC[i] >= 0);
         assert(linkCount0[i] <= 2 && linkCountLC[i] <= 2);

         if( linkCountAB < nLinkStats )
            linkSizes[size_t(linkCountAB)]++;

         if( linkCountAB == 0 && linkCountLC[i] == 0 && linkCount0[i] != 0 )
            count0++;

         if( linkCountAB == 0 && linkCount0[i] == 0 && linkCountLC[i] != 0 )
            countLC++;

         if( linkCountAB == 0 && (linkCount0[i] != 0 || linkCountLC[i] != 0) )
            count0LC++;
      }


      int nlocal = 0;
      for( int i = 0; i < nLinkStats; i++ )
         if( linkSizes[i] != 0 )
         {
            nlocal += linkSizes[i];
            std::cout << i << "-link vars: " << linkSizes[i] << "\n";
         }

      assert(n - nlocal >= 0);
      std::cout << "---total linking variables: " << n << " (global: " << n - nlocal << ")" <<   "\n";

      std::cout << "   Block0 exclusive vars " << count0 << "\n";
      std::cout << "   LC exclusive vars " << countLC << "\n";
      std::cout << "   Block0 or LC vars " << count0LC  << "\n";
   }
}

void sData::printLinkConsStats()
{
   int myl = getLocalmyl();
   int mzl = getLocalmzl();

   int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if( myl > 0 )
   {
      std::vector<int> linkCount(myl, 0);

      dynamic_cast<StochGenMatrix&>(*A).updateKLinkConsCount(linkCount);

      if( rank == 0 )
      {
         std::vector<int> linkSizes(nLinkStats, 0);

         for( int i = 0; i < myl; i++ )
            if( linkCount[i] < nLinkStats )
            {
               assert(linkCount[i] >= 0);
               linkSizes[size_t(linkCount[i])]++;
            }

         int nlocal = 0;
         for( int i = 0; i < nLinkStats; i++ )
            if( linkSizes[i] != 0 )
            {
               nlocal += linkSizes[i];
               std::cout << "equality " <<  i << "-link cons: " << linkSizes[i] << "\n";
            }
         std::cout << "---total equality linking constraints: " << myl << " (global: " << myl - nlocal << ")" <<   "\n";

      }
   }

   if( mzl > 0 )
   {
      std::vector<int> linkCount(mzl, 0);

      dynamic_cast<StochGenMatrix&>(*C).updateKLinkConsCount(linkCount);

      if( rank == 0 )
      {
         std::vector<int> linkSizes(nLinkStats, 0);

         for( int i = 0; i < mzl; i++ )
            if( linkCount[i] < nLinkStats )
            {
               assert(linkCount[i] >= 0);
               linkSizes[size_t(linkCount[i])]++;
            }

         int nlocal = 0;
         for( int i = 0; i < nLinkStats; i++ )
            if( linkSizes[i] != 0 )
            {
               nlocal += linkSizes[i];
               std::cout << "inequality " <<  i << "-link cons: " << linkSizes[i] << "\n";
            }
         std::cout << "---total inequality linking constraints: " << mzl << " (global: " << mzl - nlocal << ")" <<   "\n";
      }
   }
}

sData::~sData()
{
   for( size_t it = 0; it < children.size(); it++ )
      delete children[it];
}

std::vector<unsigned int> sData::getLinkVarsPermInv() const
{
   if( is_hierarchy_root )
      return this->children[0]->getLinkVarsPermInv();
   else
      return getInversePermutation(linkVarsPermutation);
}
std::vector<unsigned int> sData::getLinkConsEqPermInv() const
{
   if( is_hierarchy_root )
      return this->children[0]->getLinkConsEqPermInv();
   else
      return getInversePermutation(linkConsPermutationA);
}
std::vector<unsigned int> sData::getLinkConsIneqPermInv() const
{
   if( is_hierarchy_root )
      return this->children[0]->getLinkConsIneqPermInv();
   else
      return getInversePermutation(linkConsPermutationC);
}

int sData::getLocalnx() const
{
   long long my, nx;
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
      const BorderedGenMatrix& Abd = dynamic_cast<const BorderedGenMatrix&>(*A);
      Abd.border_left->getSize(my, nx);
   }
   else
   {
      const StochGenMatrix& Ast = dynamic_cast<const StochGenMatrix&>(*A);
      Ast.Bmat->getSize(my, nx);
   }
   return nx;
}

int sData::getLocalmy() const
{
   long long my, nx;
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
      const BorderedGenMatrix& Cbd = dynamic_cast<const BorderedGenMatrix&>(*C);
      Cbd.getSize(my, nx);
   }
   else
   {
      const StochGenMatrix& Ast = dynamic_cast<const StochGenMatrix&>(*A);
      Ast.Bmat->getSize(my, nx);
   }
   return my;
}

int sData::getLocalmyl() const
{
   long long myl, nxl;
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
//      const BorderedGenMatrix& Abd = dynamic_cast<const BorderedGenMatrix&>(*A);
//      Abd.Blmat->getSize(myl, nxl);
   }
   else
   {
      const StochGenMatrix& Ast = dynamic_cast<const StochGenMatrix&>(*A);
      Ast.Blmat->getSize(myl, nxl);
   }
   return myl;
}

int sData::getLocalmz() const
{
   long long mz, nx;
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
   }
   else
   {
      const StochGenMatrix& Cst = dynamic_cast<const StochGenMatrix&>(*C);
      Cst.Bmat->getSize(mz, nx);
   }
   return mz;
}

int sData::getLocalmzl() const
{
   long long mzl, nxl;
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
   }
   else
   {
      const StochGenMatrix& Cst = dynamic_cast<const StochGenMatrix&>(*C);
      Cst.Blmat->getSize(mzl, nxl);
   }
   return mzl;
}

int sData::getLocalSizes(int& nx, int& my, int& mz, int& myl, int& mzl) const
{
   long long nx_loc, my_loc, mz_loc, myl_loc, mzl_loc;

   if( is_hierarchy_root )
   {
      const BorderedGenMatrix& Abd = dynamic_cast<const BorderedGenMatrix&>(*A);
      assert(Abd.border_left->mat);
      assert(Abd.border_left->mat_link);
      Abd.border_left->mat->getSize(my_loc, nx_loc);
      Abd.bottom_left_block->getSize(myl_loc, nx_loc);

      const BorderedGenMatrix& Cbd = dynamic_cast<const BorderedGenMatrix&>(*C);
      assert(Cbd.border_left->mat);
      assert(Cbd.border_left->mat_link);
      Cbd.border_left->mat->getSize(mz_loc, nx_loc);
      Cbd.bottom_left_block->getSize(mzl_loc, nx_loc);
   }
   else
   {
      const StochGenMatrix& Ast = dynamic_cast<const StochGenMatrix&>(*A);
      Ast.Blmat->getSize(myl_loc, nx_loc);
      Ast.Bmat->getSize(my_loc, nx_loc);

      const StochGenMatrix& Cst = dynamic_cast<const StochGenMatrix&>(*C);
      Cst.Blmat->getSize(mzl_loc, nx_loc);
      Cst.Bmat->getSize(mz_loc, nx_loc);
   }

   nx = nx_loc;
   my = my_loc;
   mz = mz_loc;
   myl = myl_loc;
   mzl = mzl_loc;
   return 0;
}

int sData::getLocalSizes(int& nx, int& my, int& mz) const
{
   long long nx_loc, my_loc, mz_loc;
   if( is_hierarchy_root )
   {
      const BorderedGenMatrix& Abd = dynamic_cast<const BorderedGenMatrix&>(*A);
      assert(Abd.border_left->mat);
      Abd.border_left->mat->getSize(my_loc, nx_loc);

      const BorderedGenMatrix& Cbd = dynamic_cast<const BorderedGenMatrix&>(*C);
      assert(Cbd.border_left->mat);
      Cbd.border_left->mat->getSize(mz_loc, nx_loc);
   }
   else
   {
      const StochGenMatrix& Ast = dynamic_cast<const StochGenMatrix&>(*A);
      Ast.Bmat->getSize(my_loc, nx_loc);

      const StochGenMatrix& Cst = dynamic_cast<const StochGenMatrix&>(*C);
      Cst.Bmat->getSize(mz_loc, nx_loc);
   }

   nx = nx_loc;
   my = my_loc;
   mz = mz_loc;
   return 0;
}

int sData::getLocalNnz(int& nnzQ, int& nnzB, int& nnzD)
{
   if( is_hierarchy_root )
      assert( 0 && "TODO : implement");
   const StochSymMatrix& Qst = dynamic_cast<const StochSymMatrix&>(*Q);
   const StochGenMatrix& Ast = dynamic_cast<const StochGenMatrix&>(*A);
   const StochGenMatrix& Cst = dynamic_cast<const StochGenMatrix&>(*C);

   nnzQ = Qst.diag->getStorageRef().len + Qst.border->getStorageRef().len;
   nnzB = Ast.Bmat->getStorageRef().len;
   nnzD = Cst.Bmat->getStorageRef().len;
   return 0;
}


/*
 * At this stage we expect the Schur Complement to be of the form
 *                                                                                             nx  my0  mz0  myl  mzl
 *         [  Xsymi   0     0    X1iT FiT      X1iT GiT  ]     [ Q0  A0T  C0T  F0T  G0T  ]   [  x   x    x    x    x ]
 *         [   0      0     0       0             0      ]     [ A0   0    0    0    0   ]   [  x   0    0    0    0 ]
 * - SUM_i [   0      0     0       0             0      ]  +  [ C0   0   Om0   0    0   ] = [  x   0    x    0    0 ] =: global Schur complement (symmetric)
 *         [ Fi X1i   0     0   Fi K11i FiT  Fi K11i GiT ]     [ F0   0    0    0    0   ]   [  x   0    0    x    x ]
 *         [ Gi X1i   0     0   Gi K11i FiT  Gi K11i GiT ]     [ G0   0    0    0  OmN+1 ]   [  x   0    0    x    x ]
 *
 * since we know, that x_3 = Om0^-1 ( b_3 - C0 x1 ) we can substitute x_3 in this system and arrive at the smaller
 *
 *         [  Xsymi   0      X1iT FiT      X1iT GiT  ]     [ Q0  A0T  F0T  G0T  ]   [  x   x    x    x ]
 *         [   0      0         0             0      ]  +  [ A0   0    0    0   ] = [  x   0    0    0 ] =: global Schur complement (symmetric)
 * - SUM_i [ Fi X1i   0     Fi K11i FiT  Fi K11i GiT ]     [ F0   0    0    0   ]   [  x   0    x    x ]
 *         [ Gi X1i   0     Gi K11i FiT  Gi K11i GiT ]     [ G0   0    0  OmN+1 ]   [  x   0    x    x ]
 *
 * additionally we detect n0linkVars - Varaibles only appearing in one linking constraint block (dual to the A0/B0 block).
 * These will not get dense since no one adds to them - we add them separately with their correct number of non-zeros:
 *
 *                                                                                                   nx  my0  myl-myl00 mzl-mzl00 myl00 mzl00
 *         [  Xsymi   0   X1iT FiT      X1iT GiT  0  0 ]     [ Q0  A0T  F0T  G0T F00T   G00T  ]   [   x   x      x         x        x     x  ]
 *         [   0      0      0             0      0  0 ]     [ A0   0    0    0    0     0    ]   [   x   0      0         0        0     0  ]
 * - SUM_i [ Fi X1i   0  Fi K11i FiT  Fi K11i GiT 0  0 ]  +  [ F0   0    0    0    0     0    ] = [   x   0      x         x        0     0  ]
 *         [ Gi X1i   0  Gi K11i FiT  Gi K11i GiT 0  0 ]     [ G0   0    0  OmN+1  0     0    ]   [   x   0      x         x        0     0  ]
 *         [    0     0      0             0      0  0 ]     [ F00  0    0    0    0     0    ]   [   x   0      0         0        0     0  ]
 *         [    0     0      0             0      0  0 ]     [ G00  0    0    0    0  OmN+100 ]   [   x   0      0         0        0     x  ]
 *
 *
 *                       [ Qi BiT DiT ]^-1     [ K11 K12 K13 ]
 * Where Ki = (Ki)_lk =  [ Bi  0   0  ]      = [ K21 K22 K23 ]
 *                       [ Di  0   0  ]        [ K31 K32 K33 ] symmetric, and Om0 OmN+1 are diagonal and Xsym symmetric too.
 *
 *
 * Structure:
 * SUM_i Fi K11i FiT = SUM_i Fi K11i GiT = SUM_i Gi K11i FiT = SUM_i Gi K11i GiT
 *
 *      [ x  x  x  x  .  .  x  x  x ]
 *      [ x  x  x                   ]
 *      [ x  x  x  x                ]
 *    = [ .     x  .  .             ]
 *      [ .        .  .  .          ]
 *      [ .           .  .  .       ]
 *      [ .              .  .  x    ]
 *      [ .                 x  x  x ]
 *      [ x                    x  x ]
 *
 *      sizes depending on the blocks involved.
 *
 * For the (distributed) computation of the Schur Complement we only need to communicate
 *
 *       [  Xsymi   X1iT FiT      X1iT GiT  ]     [ Q0  F0T  G0T  ]   [  x   x   x ]
 * SUM_i [ Fi X1i  Fi K11i FiT  Fi K11i GiT ]  +  [ F0   0    0   ] = [  x   x   x ]
 *       [ Gi X1i  Gi K11i FiT  Gi K11i GiT ]     [ G0   0  OmN+1 ]   [  x   x   x ]
 *
 *       and even more structure can be exploited : define n0LinkVars variables that do only appear in the A_0(B_0) / C_0(D_0) block and in non of
 *       the other A_i, C_i
 *
 *       X1i = K11iRi + K12i Ai + K13i Ci

 *       Since Fi X1i := Fi (K11i Ri + K12i Ai + K12i Ci) = Fi (K11i [0, Ri] + K12i [0, Ai] + K13i [0, Ci] )
 *           =  [     0        Fi (K11i Ri + K12i Ai + K13i Ci] myl
 *                 n0LinkVars             n - n0LinkVars
 *
 *       So we need not to store the n0LinkVars (same for sym?) part either (neither for F, not for G), thus:
 *
 *                                                                      nx  myl mzl
 *       [  Xsymi   X1iT FiT      X1iT GiT  ]     [ Q0  F0T  G0T  ]   [  x   x   x ]
 * SUM_i [ Fi X1i  Fi K11i FiT  Fi K11i GiT ]  +  [ F0   0    0   ] = [  x   x   x ]
 *       [ Gi X1i  Gi K11i FiT  Gi K11i GiT ]     [ G0   0  OmN+1 ]   [  x   x   x ]
 *
 *
 */
int sData::getSchurCompMaxNnz()
{
   if( is_hierarchy_root )
      assert( 0 && "not available in hierarchy root");
   assert(children.size() > 0);

   const int n0 = getLocalnx();
   const int my = getLocalmy();
   const int myl = getLocalmyl();
   const int mzl = getLocalmzl();

#ifndef NDEBUG
   {
      int mB, nB;
      getLocalB().getSize(mB, nB);
      assert(mB == my  && nB == n0);
   }
#endif

   int nnz = 0;

   assert(n0 >= n0LinkVars);

   /* Xsym */
   nnz += nnzTriangular(n0);

   // add B_0 (or A_0, depending on notation)
   nnz += getLocalB().numberOfNonZeros();

   // add borders
   /* X1iT FiT */
   nnz += myl * (n0 - n0LinkVars);
   /* X1iT GiT */
   nnz += mzl * (n0 - n0LinkVars);

   // (empty) diagonal
   nnz += my;

   // add linking equality parts
   nnz += getSCdiagBlocksMaxNnz(linkStartBlockIdA.size(), linkStartBlockLengthsA);

   // add linking inequality parts
   nnz += getSCdiagBlocksMaxNnz(linkStartBlockIdC.size(), linkStartBlockLengthsC);

   // add linking mixed parts
   nnz += getSCmixedBlocksMaxNnz(linkStartBlockIdA.size(), linkStartBlockIdC.size(),
                                 linkStartBlockLengthsA, linkStartBlockLengthsC);

   if( myl > 0 )
   {
      SparseGenMatrix& Ft = getLocalF().getTranspose();
      const int* startRowFtrans = Ft.krowM();
      nnz += startRowFtrans[n0] - startRowFtrans[n0 - n0LinkVars];
   }

   if( mzl > 0 )
   {
      SparseGenMatrix& Gt = getLocalG().getTranspose();
      const int* startRowGtrans = Gt.krowM();
      nnz += startRowGtrans[n0] - startRowGtrans[n0 - n0LinkVars];
   }

   return nnz;
}


int sData::getSchurCompMaxNnzDist(int blocksStart, int blocksEnd)
{
   assert(children.size() > 0);

   const int n0 = getLocalnx();
   const int my = getLocalmy();
   const int myl = getLocalmyl();
   const int mzl = getLocalmzl();
   const int mylLocal = myl - getSCdiagBlocksNRows(linkStartBlockLengthsA)
      + getSCdiagBlocksNRows(linkStartBlockLengthsA, blocksStart, blocksEnd);
   const int mzlLocal = mzl - getSCdiagBlocksNRows(linkStartBlockLengthsC)
      + getSCdiagBlocksNRows(linkStartBlockLengthsC, blocksStart, blocksEnd);

#ifndef NDEBUG
   {
      int mB, nB;
      getLocalB().getSize(mB, nB);
      assert(mB == my  && nB == n0);
   }
#endif

   int nnz = 0;

   assert(n0 >= n0LinkVars);

   // sum up half of dense square
   nnz += nnzTriangular(n0);

   // add B_0 (or A_0, depending on notation)
   nnz += getLocalB().numberOfNonZeros();

   // add borders
   nnz += mylLocal * (n0 - n0LinkVars);
   nnz += mzlLocal * (n0 - n0LinkVars);

   // (empty) diagonal
   nnz += my;

   // add linking equality parts
   nnz += getSCdiagBlocksMaxNnzDist(linkStartBlockIdA.size(), linkStartBlockLengthsA, blocksStart, blocksEnd);

   // add linking inequality parts
   nnz += getSCdiagBlocksMaxNnzDist(linkStartBlockIdC.size(), linkStartBlockLengthsC, blocksStart, blocksEnd);

   // add linking mixed parts
   nnz += getSCmixedBlocksMaxNnzDist(linkStartBlockIdA.size(), linkStartBlockIdC.size(),
                                 linkStartBlockLengthsA, linkStartBlockLengthsC, blocksStart, blocksEnd);

   if( myl > 0 )
   {
      SparseGenMatrix& Ft = getLocalF().getTranspose();
      const int* startRowFtrans = Ft.krowM();
      nnz += startRowFtrans[n0] - startRowFtrans[n0 - n0LinkVars];
   }

   if( mzl > 0 )
   {
      SparseGenMatrix& Gt = getLocalG().getTranspose();
      const int* startRowGtrans = Gt.krowM();
      nnz += startRowGtrans[n0] - startRowGtrans[n0 - n0LinkVars];
   }

   return nnz;
}

SparseSymMatrix& sData::getLocalQ()
{
   StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*Q);
   return *Qst.diag;
}

SparseGenMatrix&
sData::getLocalCrossHessian()
{
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
   }
   StochSymMatrix& Qst = dynamic_cast<StochSymMatrix&>(*Q);
   return *Qst.border;
}

// T_i x_0 + W_i x_i = b_i

// This is T_i
SparseGenMatrix&
sData::getLocalA()
{
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
   }
   StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
   return *Ast.Amat;
}

// This is W_i:
SparseGenMatrix&
sData::getLocalB()
{
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
   }
   StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
   return *Ast.Bmat;
}

// This is F_i (linking equality matrix):
SparseGenMatrix&
sData::getLocalF()
{
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
   }
   StochGenMatrix& Ast = dynamic_cast<StochGenMatrix&>(*A);
   return *Ast.Blmat;
}

// low_i <= C_i x_0 + D_i x_i <= upp_i

// This is C_i
SparseGenMatrix&
sData::getLocalC()
{
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
   }
   StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
   return *Cst.Amat;
}

// This is D_i
SparseGenMatrix&
sData::getLocalD()
{
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
   }
   StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
   return *Cst.Bmat;
}

// This is G_i (linking inequality matrix):
SparseGenMatrix&
sData::getLocalG()
{
   if( is_hierarchy_root )
   {
      assert( 0 && "TODO : implement");
   }
   StochGenMatrix& Cst = dynamic_cast<StochGenMatrix&>(*C);
   return *Cst.Blmat;
}

void sData::cleanUpPresolvedData(const StochVectorBase<int>& rowNnzVecA, const StochVectorBase<int>& rowNnzVecC,
      const StochVectorBase<int>& colNnzVec)
{
   StochSymMatrix& Q_stoch = dynamic_cast<StochSymMatrix&>(*Q);
   // todo only works if Q is empty - not existent
   Q_stoch.deleteEmptyRowsCols(colNnzVec);

   // clean up equality system
   StochGenMatrix& A_stoch = dynamic_cast<StochGenMatrix&>(*A);
   StochVector& b_Astoch = dynamic_cast<StochVector&>(*bA);

   A_stoch.initStaticStorageFromDynamic(rowNnzVecA, colNnzVec);
   A_stoch.freeDynamicStorage();
   A_stoch.recomputeSize();

   b_Astoch.removeEntries(rowNnzVecA);

   // clean up inequality system and x
   StochGenMatrix& C_stoch = dynamic_cast<StochGenMatrix&>(*C);
   StochVector& g_stoch = dynamic_cast<StochVector&>(*g);

   StochVector& blx_stoch = dynamic_cast<StochVector&>(*blx);
   StochVector& ixlow_stoch = dynamic_cast<StochVector&>(*ixlow);
   StochVector& bux_stoch = dynamic_cast<StochVector&>(*bux);
   StochVector& ixupp_stoch = dynamic_cast<StochVector&>(*ixupp);

   StochVector& bl_stoch = dynamic_cast<StochVector&>(*bl);
   StochVector& iclow_stoch = dynamic_cast<StochVector&>(*iclow);
   StochVector& bu_stoch = dynamic_cast<StochVector&>(*bu);
   StochVector& icupp_stoch = dynamic_cast<StochVector&>(*icupp);

   C_stoch.initStaticStorageFromDynamic(rowNnzVecC, colNnzVec);
   C_stoch.freeDynamicStorage();
   C_stoch.recomputeSize();

   g_stoch.removeEntries(colNnzVec);

   blx_stoch.removeEntries(colNnzVec);
   ixlow_stoch.removeEntries(colNnzVec);
   bux_stoch.removeEntries(colNnzVec);
   ixupp_stoch.removeEntries(colNnzVec);

   bl_stoch.removeEntries(rowNnzVecC);
   iclow_stoch.removeEntries(rowNnzVecC);
   bu_stoch.removeEntries(rowNnzVecC);
   icupp_stoch.removeEntries(rowNnzVecC);

   long long dummy;
   nx = g_stoch.length();
   A_stoch.getSize( my, dummy );
   C_stoch.getSize( mz, dummy );

   nxlow = ixlow_stoch.numberOfNonzeros();
   nxupp = ixupp_stoch.numberOfNonzeros();
   mclow = iclow_stoch.numberOfNonzeros();
   mcupp = icupp_stoch.numberOfNonzeros();
}

void sData::initDistMarker(int blocksStart, int blocksEnd)
{
   assert(isSCrowLocal.size() == 0);
   assert(isSCrowMyLocal.size() == 0);

   assert(linkStartBlockIdA.size() > 0 || linkStartBlockIdC.size() > 0);
   assert(blocksStart >= 0 && blocksStart < blocksEnd && blocksEnd <= int(linkStartBlockLengthsA.size()));

   const int nx0 = getLocalnx();
   const int my0 = getLocalmy();
   const int myl = getLocalmyl();
   const int mzl = getLocalmzl();
   const int sizeSC = nx0 + my0 + myl + mzl;

   assert(sizeSC > 0);

   isSCrowLocal.resize(sizeSC);
   isSCrowMyLocal.resize(sizeSC);

   for( int i = 0; i < nx0; i++ )
   {
      isSCrowLocal[i] = false;
      isSCrowMyLocal[i] = false;
   }

   for( int i = nx0; i < nx0 + my0; i++ )
   {
      isSCrowLocal[i] = false;
      isSCrowMyLocal[i] = false;
   }

   // equality linking
   for( int i = nx0 + my0, j = 0; i < nx0 + my0 + myl; i++, j++ )
   {
      assert( unsigned(j) < linkStartBlockIdA.size() );
      const int block = linkStartBlockIdA[j];
      isSCrowLocal[i] = (block != -1);
      isSCrowMyLocal[i] = (block >= blocksStart && block < blocksEnd);
   }

   // inequality linking
   for( int i = nx0 + my0 + myl, j = 0; i < nx0 + my0 + myl + mzl; i++, j++ )
   {
      assert( unsigned(j) < linkStartBlockIdC.size() );
      const int block = linkStartBlockIdC[j];
      isSCrowLocal[i] = (block != -1);
      isSCrowMyLocal[i] = (block >= blocksStart && block < blocksEnd);
   }
}

const std::vector<bool>& sData::getSCrowMarkerLocal() const
{
   assert(isSCrowLocal.size() != 0);

   return isSCrowLocal;
}

const std::vector<bool>& sData::getSCrowMarkerMyLocal() const
{
   assert(isSCrowMyLocal.size() != 0);

   return isSCrowMyLocal;
}

int sData::n2linkRowsEq() const
{
   return n2linksRows(linkStartBlockLengthsA);
}

int sData::n2linkRowsIneq() const
{
   return n2linksRows(linkStartBlockLengthsC);
}

// is root node data of sData object same on all procs?
bool sData::isRootNodeInSync() const
{
   bool in_sync = true;

   /* matrix Q */
   // todo

   /* matrix A */
   if(!dynamic_cast<const StochGenMatrix&>(*A).isRootNodeInSync())
   {
      std::cout << "ERROR: matrix A corrupted!\n";
      in_sync = false;
   }

   /* matrix C */
   if( !dynamic_cast<const StochGenMatrix&>(*C).isRootNodeInSync() )
   {
      std::cout << "ERROR: matrix C corrupted!\n";
      in_sync = false;
   }

   /* objective g */
   if( !dynamic_cast<const StochVector&>(*g).isRootNodeInSync() )
   {
      std::cout << "ERROR: objective vector corrupted!\n";
      in_sync = false;
   }

   /* rhs equality bA */
   if( !dynamic_cast<const StochVector&>(*bA).isRootNodeInSync() )
   {
      std::cout << "ERROR: rhs of A corrupted!\n";
      in_sync = false;
   }

   /* upper bounds x bux */
   if( !dynamic_cast<const StochVector&>(*bux).isRootNodeInSync() )
   {
      std::cout << "ERROR: upper bounds x corrupted!\n";
      in_sync = false;
   }

   /* index for upper bounds x ixupp */
   if( !dynamic_cast<const StochVector&>(*ixupp).isRootNodeInSync() )
   {
      std::cout << "ERROR: index upper bounds x corrupted!\n";
      in_sync = false;
   }

   /* lower bounds x blx */
   if( !dynamic_cast<const StochVector&>(*blx).isRootNodeInSync() )
   {
      std::cout << "ERROR: lower bounds x corrupted!\n";
      in_sync = false;
   }

   /* index for lower bounds x ixlow */
   if( !dynamic_cast<const StochVector&>(*ixlow).isRootNodeInSync() )
   {
      std::cout << "ERROR: index lower bounds x corrupted!\n";
      in_sync = false;
   }

   /* upper bounds C bu */
   if( !dynamic_cast<const StochVector&>(*bu).isRootNodeInSync() )
   {
      std::cout << "ERROR: rhs C corrupted!\n";
      in_sync = false;
   }

   /* index upper bounds C icupp */
   if( !dynamic_cast<const StochVector&>(*icupp).isRootNodeInSync() )
   {
      std::cout << "ERROR: index rhs C corrupted!\n";
      in_sync = false;
   }

   /* lower bounds C bl */
   if( !dynamic_cast<const StochVector&>(*bl).isRootNodeInSync() )
   {
      std::cout << "ERROR: lower bounds C corrupted!\n";
      in_sync = false;
   }

   /* index for lower bounds C iclow */
   if( !dynamic_cast<const StochVector&>(*iclow).isRootNodeInSync() )
   {
      std::cout << "ERROR: index lower bounds C corrupted!\n";
      in_sync = false;
   }

   /* scale sc */
   // todo

   return in_sync;
}

void sData::printRanges() const
{
   /* objective */
   double absmin_objective; g->absminNonZero( absmin_objective, 0.0 );
   assert( absmin_objective >= 0 );
   const double absmax_objective = g->infnorm();
   assert( absmax_objective >= 0 );

   /* matrix range */
   const double absmax_A = A->abmaxnorm();
   const double absmax_C = C->abmaxnorm();

   const double absmin_A = A->abminnormNonZero();
   const double absmin_C = C->abminnormNonZero();

   const double mat_min = std::min( absmin_A, absmin_C );
   const double mat_max = std::max( absmax_A, absmax_C );

   /* rhs range */
   double absmin_bA; bA->absminNonZero( absmin_bA, 0.0 );
   double absmin_bl; bl->absminNonZero( absmin_bl, 0.0 );
   double absmin_bu; bu->absminNonZero( absmin_bu, 0.0 );

   const double absmax_bA = bA->infnorm();
   const double absmax_bl = bl->infnorm();
   const double absmax_bu = bu->infnorm();

   const double rhs_min = std::min( absmin_bA, std::min( absmin_bl, absmin_bu ) );
   const double rhs_max = std::max( absmax_bA, std::max( absmax_bl, absmax_bu ) );

   /* bounds range */
   double absmin_blx; blx->absminNonZero( absmin_blx, 0.0 );
   double absmin_bux; bux->absminNonZero( absmin_bux, 0.0 );

   const double absmax_blx = blx->infnorm();
   const double absmax_bux = bux->infnorm();

   const double bounds_min = std::min( absmin_blx, absmin_bux );
   const double bounds_max = std::max( absmax_blx, absmax_bux );

   if( PIPS_MPIgetRank() == 0 )
   {
      const double inf = std::numeric_limits<double>::infinity();

      const std::streamsize pre_old = std::cout.precision();
      std::cout << std::setprecision(0) << std::scientific;
      std::cout << "Matrix range    [" << mat_min << ", " << mat_max << "]\n";
      std::cout << "Objective range [" << ( absmin_objective == inf ? 0.0 : absmin_objective ) << ", " << absmax_objective << "]\n";
      std::cout << "Bounds range    [" << ( bounds_min == inf ? 0.0 : bounds_min ) << ", " << bounds_max << "]\n";
      std::cout << "RhsLhs range    [" << ( rhs_min == inf ? 0.0 : rhs_min ) << ", " << rhs_max << "]\n";
      std::cout << std::setprecision(pre_old) << std::defaultfloat;
   }
}
