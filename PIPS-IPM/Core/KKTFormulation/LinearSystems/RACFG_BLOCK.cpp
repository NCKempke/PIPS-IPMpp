//
// Created by nils-christian on 15.06.21.
//

#include "RACFG_BLOCK.h"

template<>
const std::unique_ptr<const SparseMatrix> RACFG_BLOCK<SparseMatrix>::dummy = std::make_unique<SparseMatrix>();

template<>
const std::unique_ptr<const StripMatrix> RACFG_BLOCK<StripMatrix>::dummy = std::make_unique<StripMatrix>();

BorderLinsys getChild(BorderLinsys& border, unsigned int i) {
   const bool is_dummy = border.F.children[i]->is_a(kStringGenDummyMatrix);
   assert(i < border.F.children.size());
   if (border.has_RAC) {
      if (!is_dummy && border.F.children[i]->first->is_a(kStripMatrix))
         return BorderLinsys(dynamic_cast<StripMatrix&>(*border.R.children[i]->first),
            dynamic_cast<StripMatrix&>(*border.A.children[i]->first), dynamic_cast<StripMatrix&>(*border.C.children[i]->first),
            border.n_empty_rows, dynamic_cast<StripMatrix&>(*border.F.children[i]->first),
            dynamic_cast<StripMatrix&>(*border.G.children[i]->first));
      else
         return BorderLinsys(*border.R.children[i], *border.A.children[i], *border.C.children[i], border.n_empty_rows, *border.F.children[i],
            *border.G.children[i]);
   }
   else {
      if (!is_dummy && border.F.children[i]->first->is_a(kStripMatrix))
         return BorderLinsys(border.n_empty_rows, dynamic_cast<StripMatrix&>(*border.F.children[i]->first),
            dynamic_cast<StripMatrix&>(*border.G.children[i]->first), border.use_local_RAC);
      else
         return BorderLinsys(border.n_empty_rows, *border.F.children[i], *border.G.children[i], border.use_local_RAC);
   }
}

template<>
bool BorderLinsys::isEmpty() const {
   if (use_local_RAC)
      return false;
   else {
      if (F.numberOfNonZeros() == 0 && G.numberOfNonZeros() == 0) {
         if (has_RAC)
            return R.numberOfNonZeros() == 0 && A.numberOfNonZeros() == 0 && C.numberOfNonZeros() == 0;
         else
            return true;
      } else
         return false;
   }
}

template<>
bool BorderBiBlock::isEmpty() const {
   if (use_local_RAC)
      return false;
   else {
      if (F.numberOfNonZeros() == 0 && G.numberOfNonZeros() == 0) {
         if (has_RAC)
            return R.numberOfNonZeros() == 0 && A.numberOfNonZeros() == 0 && C.numberOfNonZeros() == 0;
         else
            return true;
      } else
         return false;
   }
}
