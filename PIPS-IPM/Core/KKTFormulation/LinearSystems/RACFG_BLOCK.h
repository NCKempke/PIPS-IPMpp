//
// Created by nils-christian on 15.06.21.
//

#ifndef PIPSIPMPP_RACFG_BLOCK_H
#define PIPSIPMPP_RACFG_BLOCK_H

#include "SparseMatrix.h"
#include "StripMatrix.h"

#include <memory>

template<typename T>
struct RACFG_BLOCK {
private:
   static const std::unique_ptr<const T> dummy;

public:
   const bool use_local_RAC{};
   const bool has_RAC{};
   /* represents a block like
    * [ R_i 0 F_i^T G_i^T ]             [ R_i^T A_i^T C_i^T ]
    * [ A_i 0   0     0   ] or possibly [   0     0     0   ]
    * [ C_i 0   0     0   ]             [  F_i    0     0   ]
    *                                   [  G_i    0     0   ]
    */
   const T& R;
   const T& A;
   const T& C;
   const T& F;
   const T& G;

   /* n_empty_rows gives the distance between RAC and F,G blocks */
   const int n_empty_rows;

   [[nodiscard]] bool isEmpty() const;

   RACFG_BLOCK(const T& R, const T& A, const T& C, int n_empty_rows, const T& F, const T& G) : has_RAC{true}, R{R}, A{A}, C{C}, F{F}, G{G}, n_empty_rows{n_empty_rows} {
      assert(n_empty_rows >= 0);
   };

   RACFG_BLOCK(int n_empty_rows, const T& F, const T& G, bool use_local_RAC) : use_local_RAC{use_local_RAC}, has_RAC{false},
   R{*RACFG_BLOCK<T>::dummy}, A{*RACFG_BLOCK<T>::dummy}, C{*RACFG_BLOCK<T>::dummy},
      F{F}, G{G}, n_empty_rows{n_empty_rows} { assert(n_empty_rows >= 0); };

   RACFG_BLOCK(const RACFG_BLOCK<T>& block) : use_local_RAC{block.use_local_RAC}, has_RAC{block.has_RAC}, R{block.R}, A{block.A}, C{block.C},
      F{block.F}, G{block.G}, n_empty_rows{block.n_empty_rows} { assert(n_empty_rows >= 0); };
};

using BorderLinsys = RACFG_BLOCK<StripMatrix>;
using BorderBiBlock = RACFG_BLOCK<SparseMatrix>;

BorderLinsys getChild(BorderLinsys& border, unsigned int i);

#endif //PIPSIPMPP_RACFG_BLOCK_H
