#include "p3io.h"
#include "exceptions.h"
#include "math_p3.h"

SYSTEM_double MATH_P3_minsingle = 1.5e-45;
SYSTEM_double MATH_P3_maxsingle = 3.4e38;
SYSTEM_double MATH_P3_mindouble = 5.0e-324;
SYSTEM_double MATH_P3_maxdouble = 1.7e308;
typedef struct MATH_P3_ti64rec_S {
  union{
    struct{
      SYSTEM_double x;
    } _c1;
    struct{
      SYSTEM_int64 i64;
    } _c2;
  } _u;
} MATH_P3_ti64rec;

static SYSTEM_int64 MATH_P3_signmask;
static SYSTEM_int64 MATH_P3_expomask;
static SYSTEM_int64 MATH_P3_mantmask;
/**** C code included from math_p3.pas(105:1): 33 lines ****/
#if   defined(AIX)
# include <float.h>
# include <fpxcp.h>
# include <fptrap.h>

#elif defined(__APPLE__)
# include <fenv.h>
/* SPD, Oct 2007: N.B.: to get more exceptions raised, set more bits in the input
 * to get fewer exceptions and more non-stop arithmetic, set fewer bits
 * or maybe it is the other way around on DII, still testing
 */

#elif defined(BGP) || defined(__linux__)
# include <fpu_control.h>
# include <fenv.h>
#elif defined(SOL)
# include <ieeefp.h>
#endif
#if defined(__linux__) && defined(__GNUC__) /* Linux GCC compilers */
# if ! defined(FE_DENORMAL)
#  define FE_DENORMAL __FE_DENORM
# endif
#endif

#define ADD2MASK(X) _P3SET_p(result,_len_ret,result,X)
#define ISINMASK(X) (_P3SET_ic(5,X,mask))
#define EX_INVALIDOP   _P3set1("\001")
#define EX_DENORMAL    _P3set1("\002")
#define EX_ZERODIVIDE  _P3set1("\004")
#define EX_OVERFLOW    _P3set1("\010")
#define EX_UNDERFLOW   _P3set1("\020")
#define EX_PRECISION   _P3set1("\040")


Function(SYSTEM_double ) MATH_P3_arccos(
  SYSTEM_double x)
{
  SYSTEM_double result;

  /**** C code included from math_p3.pas(147:1): 1 lines ****/
  result = acos(x);
  return result;
}  /* arccos */

Function(SYSTEM_double ) MATH_P3_arcsin(
  SYSTEM_double x)
{
  SYSTEM_double result;

  /**** C code included from math_p3.pas(158:1): 1 lines ****/
  result = asin(x);
  return result;
}  /* arcsin */

Function(SYSTEM_double ) MATH_P3_arctan2(
  SYSTEM_double y,
  SYSTEM_double x)
{
  SYSTEM_double result;

  /**** C code included from math_p3.pas(181:1): 1 lines ****/
  result = atan2(y,x);
  return result;
}  /* arctan2 */

Function(SYSTEM_double ) MATH_P3_tan(
  SYSTEM_double x)
{
  SYSTEM_double result;

  /**** C code included from math_p3.pas(191:1): 1 lines ****/
  result = tan(x);
  return result;
}  /* tan */

Function(SYSTEM_double ) MATH_P3_intpower(
  SYSTEM_double x,
  SYSTEM_integer i)
{
  SYSTEM_double result;
  SYSTEM_integer y;

  y = SYSTEM_abs_i(i);
  result = 1.0;
  while (y > 0) {
    while (!SYSTEM_odd(y)) {
      y = ValueCast(SYSTEM_uint32,y) >> 1;
      x = x * x;
    }
    _P3dec0(y);
    result = result * x;
  }
  if (i < 0) 
    result = 1.0 /  result;
  return result;
}  /* intpower */

Function(SYSTEM_double ) MATH_P3_power(
  SYSTEM_double x,
  SYSTEM_double y)
{
  SYSTEM_double result;

  if (y == 0.0) { 
    result = 1.0;
  } else 
    if (x == 0.0 && y > 0.0) { 
      result = 0.0;
    } else 
      if (SYSTEM_frac(y) == 0.0 && SYSTEM_abs_r(y) <= SYSTEM_maxint) { 
        result = MATH_P3_intpower(x,SYSTEM_trunc(y));
      } else 
        result = SYSTEM_exp(y * SYSTEM_ln(x));
  return result;
}  /* power */
typedef SYSTEM_int8 MATH_P3_tfoorange;

Function(SYSTEM_double ) MATH_P3_roundto(
  SYSTEM_double x,
  SYSTEM_integer i)
{
  SYSTEM_double result;
  typedef SYSTEM_uint8 _sub_0ROUNDTO;
  typedef SYSTEM_double tfactors[2];
  typedef tfactors _arr_1ROUNDTO[41];
  static _arr_1ROUNDTO lfactorarray = {{1e-20, 1e20}, {1e-19, 1e19}, {1e-18, 1e18}, 
    {1e-17, 1e17}, {1e-16, 1e16}, {1e-15, 1e15}, {1e-14, 1e14}, {1e-13, 1e13}, 
    {1e-12, 1e12}, {1e-11, 1e11}, {1e-10, 1e10}, {1e-09, 1e09}, {1e-08, 1e08}, 
    {1e-07, 1e07}, {1e-06, 1e06}, {1e-05, 1e05}, {1e-04, 1e04}, {1e-03, 1e03}, 
    {1e-02, 1e02}, {1e-01, 1e01}, {1, 1}, {1e01, 1e-01}, {1e02, 1e-02}, 
    {1e03, 1e-03}, {1e04, 1e-04}, {1e05, 1e-05}, {1e06, 1e-06}, {1e07, 1e-07}, 
    {1e08, 1e-08}, {1e09, 1e-09}, {1e10, 1e-10}, {1e11, 1e-11}, {1e12, 1e-12}, 
    {1e13, 1e-13}, {1e14, 1e-14}, {1e15, 1e-15}, {1e16, 1e-16}, {1e17, 1e-17}, 
    {1e18, 1e-18}, {1e19, 1e-19}, {1e20, 1e-20}};
  static SYSTEM_double roundmax = 4503599627370496.;
  SYSTEM_double t, lo, up, dlo, dup;
  SYSTEM_int64 i64;

  if (i <  -20 || i > 20) 
    _P3_RAISE(ValueCast(EXCEPTIONS_einvalidop,
      SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,
      _P3alloc_object(&EXCEPTIONS_einvalidop_CD)),_P3str1("\035RoundTo arg 2 is out of range"))));
  t = SYSTEM_abs_r(x) * lfactorarray[i -  -20][1];
  if (t < roundmax) {
    i64 = SYSTEM_round(t);
    if (i64 > t) 
      _P3dec0(i64);
    lo = i64;
    up = lo + 1;
    dlo = t - lo;
    dup = up - t;
    if (dlo < dup) { 
      t = lo;
    } else 
      if (dlo > dup) { 
        t = up;
      } else 
        if (SYSTEM_odd(i64)) { 
          t = up;
        } else 
          t = lo;
  } 
  result = t * lfactorarray[i -  -20][0];
  if (x < 0) 
    result = -result;
  return result;
}  /* roundto */

Function(SYSTEM_boolean ) MATH_P3_isnan(
  SYSTEM_double avalue)
{
  SYSTEM_boolean result;
  MATH_P3_ti64rec i64rec;
  SYSTEM_int64 mantissa;

  result = SYSTEM_false;
  i64rec._u._c1.x = avalue;
  if (ValueCast(SYSTEM_uint64,i64rec._u._c2.i64 & MATH_P3_expomask) >> 52 == 2047) {
    mantissa = i64rec._u._c2.i64 & MATH_P3_mantmask;
    if (mantissa != 0) 
      result = SYSTEM_true;
  } 
  return result;
}  /* isnan */

Function(SYSTEM_boolean ) MATH_P3_isinfinite(
  SYSTEM_double avalue)
{
  SYSTEM_boolean result;
  MATH_P3_ti64rec i64rec;
  SYSTEM_int64 mantissa;

  result = SYSTEM_false;
  i64rec._u._c1.x = avalue;
  if (ValueCast(SYSTEM_uint64,i64rec._u._c2.i64 & MATH_P3_expomask) >> 52 == 2047) {
    mantissa = i64rec._u._c2.i64 & MATH_P3_mantmask;
    if (mantissa == 0) 
      result = SYSTEM_true;
  } 
  return result;
}  /* isinfinite */

Function(SYSTEM_double ) MATH_P3_lnxp1(
  SYSTEM_double x)
{
  SYSTEM_double result;

  /**** C code included from math_p3.pas(334:1): 1 lines ****/
  result = log1p(x);
  return result;
}  /* lnxp1 */

Function(_P3set_elem *) MATH_P3_getexceptionmask(
  _P3set_elem *result,
  SYSTEM_uint8 _len_ret)
{
  _P3SET_copy(result,_len_ret,_P3empty_set);
  /**** C code included from math_p3.pas(355:1): 85 lines ****/
#if defined(_WIN32)
{
  unsigned int cw;

  cw = _control87(0,0);
  if (cw & _EM_INVALID   ) ADD2MASK(EX_INVALIDOP );
  if (cw & _EM_DENORMAL  ) ADD2MASK(EX_DENORMAL  );
  if (cw & _EM_ZERODIVIDE) ADD2MASK(EX_ZERODIVIDE);
  if (cw & _EM_OVERFLOW  ) ADD2MASK(EX_OVERFLOW  );
  if (cw & _EM_UNDERFLOW ) ADD2MASK(EX_UNDERFLOW );
  if (cw & _EM_INEXACT   ) ADD2MASK(EX_PRECISION );
}
#elif defined(AIX)
  if (!fp_is_enabled(TRP_INVALID    )) ADD2MASK(EX_INVALIDOP );
                                       ADD2MASK(EX_DENORMAL  );
  if (!fp_is_enabled(TRP_DIV_BY_ZERO)) ADD2MASK(EX_ZERODIVIDE);
  if (!fp_is_enabled(TRP_OVERFLOW   )) ADD2MASK(EX_OVERFLOW  );
  if (!fp_is_enabled(TRP_UNDERFLOW  )) ADD2MASK(EX_UNDERFLOW );
  if (!fp_is_enabled(TRP_INEXACT    )) ADD2MASK(EX_PRECISION );
#elif defined (BGP) 
{
  fpu_control_t cw;

  _FPU_GETCW(cw);
  if (!(cw & _FPU_MASK_IM)) ADD2MASK(EX_INVALIDOP );
                            ADD2MASK(EX_DENORMAL  );
  if (!(cw & _FPU_MASK_ZM)) ADD2MASK(EX_ZERODIVIDE);
  if (!(cw & _FPU_MASK_OM)) ADD2MASK(EX_OVERFLOW  );
  if (!(cw & _FPU_MASK_UM)) ADD2MASK(EX_UNDERFLOW );
  if (!(cw & _FPU_MASK_XM)) ADD2MASK(EX_PRECISION );
}
#elif defined(__APPLE__) && defined(__GNUC__) /* Mac, GCC compilers */
{
  fenv_t fenv;
  unsigned int ex;

  (void) fegetenv (&fenv);
  ex = fenv.__control & FE_ALL_EXCEPT;
  if (ex & FE_INVALID   ) ADD2MASK(EX_INVALIDOP );
                          ADD2MASK(EX_DENORMAL  );
  if (ex & FE_DIVBYZERO ) ADD2MASK(EX_ZERODIVIDE);
  if (ex & FE_OVERFLOW  ) ADD2MASK(EX_OVERFLOW  );
  if (ex & FE_UNDERFLOW ) ADD2MASK(EX_UNDERFLOW );
  if (ex & FE_INEXACT   ) ADD2MASK(EX_PRECISION );
}
#elif defined(__APPLE__) /* Mac, presumably Intel compilers */
{
  fenv_t fenv;
  unsigned short cw;

  (void) fegetenv (&fenv);
  cw = fenv.__control;
  if (cw & FE_INVALID   ) ADD2MASK(EX_INVALIDOP );
                          ADD2MASK(EX_DENORMAL  );
  if (cw & FE_DIVBYZERO ) ADD2MASK(EX_ZERODIVIDE);
  if (cw & FE_OVERFLOW  ) ADD2MASK(EX_OVERFLOW  );
  if (cw & FE_UNDERFLOW ) ADD2MASK(EX_UNDERFLOW );
  if (cw & FE_INEXACT   ) ADD2MASK(EX_PRECISION );
}
#elif defined(__linux__)
{
  fpu_control_t cw;

  _FPU_GETCW(cw);
  /* this code depends on the bits in the set element result being
   * ordered in a certain way (like Delphi does)
   */
  *result = cw & 0x3f;
}
#elif defined(SOL)
{
  fp_except flags;

  flags = fpgetmask();
  if (!(flags & FP_X_INV)) ADD2MASK(EX_INVALIDOP );
                           ADD2MASK(EX_DENORMAL  );
  if (!(flags & FP_X_DZ )) ADD2MASK(EX_ZERODIVIDE);
  if (!(flags & FP_X_OFL)) ADD2MASK(EX_OVERFLOW  );
  if (!(flags & FP_X_UFL)) ADD2MASK(EX_UNDERFLOW );
  if (!(flags & FP_X_IMP)) ADD2MASK(EX_PRECISION );
}
#else
# error "This OS not yet implemented"
  is_not_implemented;
#endif
  return result;
}  /* getexceptionmask */

Function(_P3set_elem *) MATH_P3_setexceptionmask(
  _P3set_elem *result,
  SYSTEM_uint8 _len_ret,
  const _P3set_elem *mask)
{
  _P3SET_copy(result,_len_ret,_P3empty_set);
  /**** C code included from math_p3.pas(470:1): 177 lines ****/

#if defined(_WIN32)
{
  unsigned int cw, tcw;

  cw = _control87(0,0);
  if (cw & _EM_INVALID   ) ADD2MASK(EX_INVALIDOP );
  if (cw & _EM_DENORMAL  ) ADD2MASK(EX_DENORMAL  );
  if (cw & _EM_ZERODIVIDE) ADD2MASK(EX_ZERODIVIDE);
  if (cw & _EM_OVERFLOW  ) ADD2MASK(EX_OVERFLOW  );
  if (cw & _EM_UNDERFLOW ) ADD2MASK(EX_UNDERFLOW );
  if (cw & _EM_INEXACT   ) ADD2MASK(EX_PRECISION );
  tcw = 0;
  if (ISINMASK(MATH_P3_exinvalidop   )) tcw |= _EM_INVALID   ;
  if (ISINMASK(MATH_P3_exdenormalized)) tcw |= _EM_DENORMAL  ;
  if (ISINMASK(MATH_P3_exzerodivide  )) tcw |= _EM_ZERODIVIDE;
  if (ISINMASK(MATH_P3_exoverflow    )) tcw |= _EM_OVERFLOW  ;
  if (ISINMASK(MATH_P3_exunderflow   )) tcw |= _EM_UNDERFLOW ;
  if (ISINMASK(MATH_P3_exprecision   )) tcw |= _EM_INEXACT   ;
  (void) _control87 (tcw,_MCW_EM);
}
#elif defined(AIX)
{
  fptrap_t traps;

  if (!fp_is_enabled(TRP_INVALID    )) ADD2MASK(EX_INVALIDOP );
                                       ADD2MASK(EX_DENORMAL  );
  if (!fp_is_enabled(TRP_DIV_BY_ZERO)) ADD2MASK(EX_ZERODIVIDE);
  if (!fp_is_enabled(TRP_OVERFLOW   )) ADD2MASK(EX_OVERFLOW  );
  if (!fp_is_enabled(TRP_UNDERFLOW  )) ADD2MASK(EX_UNDERFLOW );
  if (!fp_is_enabled(TRP_INEXACT    )) ADD2MASK(EX_PRECISION );
  traps = 0;
  if (!ISINMASK(MATH_P3_exinvalidop )) traps |= TRP_INVALID  ;
  if (!ISINMASK(MATH_P3_exzerodivide)) traps |= TRP_DIV_BY_ZERO;
  if (!ISINMASK(MATH_P3_exoverflow  )) traps |= TRP_OVERFLOW ;
  if (!ISINMASK(MATH_P3_exunderflow )) traps |= TRP_UNDERFLOW;
  if (!ISINMASK(MATH_P3_exprecision )) traps |= TRP_INEXACT  ;
  fp_enable(traps);
}

#elif defined (BGP)
{
  fpu_control_t cw;

  /* Need to get control word to preserve rounding mode and other bits */
  _FPU_GETCW(cw);

  if (!(cw & _FPU_MASK_IM)) ADD2MASK(EX_INVALIDOP );
                            ADD2MASK(EX_DENORMAL  );
  if (!(cw & _FPU_MASK_ZM)) ADD2MASK(EX_ZERODIVIDE);
  if (!(cw & _FPU_MASK_OM)) ADD2MASK(EX_OVERFLOW  );
  if (!(cw & _FPU_MASK_UM)) ADD2MASK(EX_UNDERFLOW );
  if (!(cw & _FPU_MASK_XM)) ADD2MASK(EX_PRECISION );

  if (!ISINMASK(MATH_P3_exinvalidop )) cw |=  _FPU_MASK_IM;
  else                                 cw &= ~_FPU_MASK_IM;
  if (!ISINMASK(MATH_P3_exzerodivide)) cw |=  _FPU_MASK_ZM;
  else                                 cw &= ~_FPU_MASK_ZM;
  if (!ISINMASK(MATH_P3_exoverflow  )) cw |=  _FPU_MASK_OM;
  else                                 cw &= ~_FPU_MASK_OM;
  if (!ISINMASK(MATH_P3_exunderflow )) cw |=  _FPU_MASK_UM;
  else                                 cw &= ~_FPU_MASK_UM;
  if (!ISINMASK(MATH_P3_exprecision )) cw |=  _FPU_MASK_XM;
  else                                 cw &= ~_FPU_MASK_XM;
  _FPU_SETCW(cw);
}
#elif defined(__APPLE__) && defined(__GNUC__) /* Mac, GCC compilers */
/* try something different for GNU compilers */
{
  fenv_t fenv;
  unsigned int oldEx, newEx;

  (void) fegetenv (&fenv);
  oldEx = fenv.__control & FE_ALL_EXCEPT;
  if (oldEx & FE_INVALID   ) ADD2MASK(EX_INVALIDOP );
                             ADD2MASK(EX_DENORMAL  );
  if (oldEx & FE_DIVBYZERO ) ADD2MASK(EX_ZERODIVIDE);
  if (oldEx & FE_OVERFLOW  ) ADD2MASK(EX_OVERFLOW  );
  if (oldEx & FE_UNDERFLOW ) ADD2MASK(EX_UNDERFLOW );
  if (oldEx & FE_INEXACT   ) ADD2MASK(EX_PRECISION );

  newEx = 0;
  if (ISINMASK(MATH_P3_exinvalidop )) newEx |= FE_INVALID  ;
  if (ISINMASK(MATH_P3_exzerodivide)) newEx |= FE_DIVBYZERO;
  if (ISINMASK(MATH_P3_exoverflow  )) newEx |= FE_OVERFLOW ;
  if (ISINMASK(MATH_P3_exunderflow )) newEx |= FE_UNDERFLOW;
  if (ISINMASK(MATH_P3_exprecision )) newEx |= FE_INEXACT  ;
  fenv.__control &= ~FE_ALL_EXCEPT;
  fenv.__control |= newEx;
# if __WORDSIZE == 64
#  if 0
  fenv.__mxcsr   &= ~(FE_ALL_EXCEPT << 7);
  fenv.__mxcsr   |= (newEx << 7);
#  endif
# endif
  (void) fesetenv (&fenv);
} /* DEG,DIG */
#elif defined(__APPLE__) /* Mac, presumably Intel compilers */
{
  fenv_t fenv;
  unsigned short cw;

  (void) fegetenv (&fenv);
  cw = fenv.__control;
  if (cw & FE_INVALID   ) ADD2MASK(EX_INVALIDOP );
                          ADD2MASK(EX_DENORMAL  );
  if (cw & FE_DIVBYZERO ) ADD2MASK(EX_ZERODIVIDE);
  if (cw & FE_OVERFLOW  ) ADD2MASK(EX_OVERFLOW  );
  if (cw & FE_UNDERFLOW ) ADD2MASK(EX_UNDERFLOW );
  if (cw & FE_INEXACT   ) ADD2MASK(EX_PRECISION );

  if (!ISINMASK(MATH_P3_exinvalidop )) cw |= FE_INVALID  ;
  if (!ISINMASK(MATH_P3_exzerodivide)) cw |= FE_DIVBYZERO;
  if (!ISINMASK(MATH_P3_exoverflow  )) cw |= FE_OVERFLOW ;
  if (!ISINMASK(MATH_P3_exunderflow )) cw |= FE_UNDERFLOW;
  if (!ISINMASK(MATH_P3_exprecision )) cw |= FE_INEXACT  ;
  fenv.__control = cw;
  (void) fesetenv (&fenv);
} /* DEI,DII */
#elif defined(__linux__)
{
#if 1
  fenv_t fenv;
  unsigned short cw,newcw;

  (void) fegetenv (&fenv);
  cw = fenv.__control_word;
  if (cw & FE_INVALID   ) ADD2MASK(EX_INVALIDOP );
  if (cw & FE_DENORMAL  ) ADD2MASK(EX_DENORMAL  );
  if (cw & FE_DIVBYZERO ) ADD2MASK(EX_ZERODIVIDE);
  if (cw & FE_OVERFLOW  ) ADD2MASK(EX_OVERFLOW  );
  if (cw & FE_UNDERFLOW ) ADD2MASK(EX_UNDERFLOW );
  if (cw & FE_INEXACT   ) ADD2MASK(EX_PRECISION );
  newcw = 0;
  if (ISINMASK(MATH_P3_exinvalidop   )) newcw |= FE_INVALID  ;
  if (ISINMASK(MATH_P3_exdenormalized)) newcw |= FE_DENORMAL ;
  if (ISINMASK(MATH_P3_exzerodivide  )) newcw |= FE_DIVBYZERO;
  if (ISINMASK(MATH_P3_exoverflow    )) newcw |= FE_OVERFLOW ;
  if (ISINMASK(MATH_P3_exunderflow   )) newcw |= FE_UNDERFLOW;
  if (ISINMASK(MATH_P3_exprecision   )) newcw |= FE_INEXACT  ;
  fenv.__control_word = newcw;
  (void) fesetenv (&fenv);
#else
  fpu_control_t cw1, tcw;

  /* this code depends on the bits in the set element result being
   * ordered in a certain way (like Delphi does)
   */
  _FPU_GETCW(cw1);
  tcw = (cw1 & ~0x3f) | ((*mask) & 0x3f);
  _FPU_SETCW(tcw);
  *result = cw1 & 0x3f;
#endif
}
#elif defined(SOL)
{
  fp_except flgs;

  flgs = fpgetmask();
  if (!(flgs & FP_X_INV)) ADD2MASK(EX_INVALIDOP );
                          ADD2MASK(EX_DENORMAL  );
  if (!(flgs & FP_X_DZ )) ADD2MASK(EX_ZERODIVIDE);
  if (!(flgs & FP_X_OFL)) ADD2MASK(EX_OVERFLOW  );
  if (!(flgs & FP_X_UFL)) ADD2MASK(EX_UNDERFLOW );
  if (!(flgs & FP_X_IMP)) ADD2MASK(EX_PRECISION );
  flgs = 0;
  if (!(ISINMASK(MATH_P3_exinvalidop ))) flgs |= FP_X_INV;
  if (!(ISINMASK(MATH_P3_exzerodivide))) flgs |= FP_X_DZ ;
  if (!(ISINMASK(MATH_P3_exoverflow  ))) flgs |= FP_X_OFL;
  if (!(ISINMASK(MATH_P3_exunderflow ))) flgs |= FP_X_UFL;
  if (!(ISINMASK(MATH_P3_exprecision ))) flgs |= FP_X_IMP;
  (void) fpsetmask(flgs);
}
#else
# error "This OS not yet implemented"
  is_not_implemented;
#endif
  return result;
}  /* setexceptionmask */

Procedure MATH_P3_setexceptionmask2p3(void)
{
}  /* setexceptionmask2p3 */

Procedure MATH_P3_clearexceptions(void)
{
  /**** C code included from math_p3.pas(681:1): 16 lines ****/
#if defined(_WIN32)
  (void) _clearfp();
#elif defined(AIX)
  fp_clr_flag (FP_INVALID | FP_OVERFLOW | FP_UNDERFLOW |
               FP_DIV_BY_ZERO | FP_INEXACT);
#elif defined(__APPLE__)
{
assert(0);
}
#elif defined(BGP) || defined(__linux__)
  (void) feclearexcept (FE_ALL_EXCEPT);
#elif defined(SOL)
  (void) fpsetsticky(0);
#else
# error "This OS not yet implemented"
#endif
}  /* clearexceptions */

/* unit math_p3 */
void _Init_Module_math_p3(void)
{
  MATH_P3_signmask =  SYSTEM_minint;
  MATH_P3_signmask = ValueCast(SYSTEM_uint64,MATH_P3_signmask) << 32;
  MATH_P3_expomask = 2146435072;
  MATH_P3_expomask = ValueCast(SYSTEM_uint64,MATH_P3_expomask) << 32;
  MATH_P3_mantmask = ~(MATH_P3_signmask | MATH_P3_expomask);
} /* _Init_Module_math_p3 */

void _Final_Module_math_p3(void)
{
} /* _Final_Module_math_p3 */

