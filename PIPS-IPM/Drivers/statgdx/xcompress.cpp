#include "p3io.h"
#include "p3platform.h"
#include "system_p3.h"
#include "p3utils.h"
#include "p3process.h"
#include "p3library.h"
#include "exceptions.h"
#include "math_p3.h"
#include "p3ieeefp.h"
#include "sysutils_p3.h"
#include "p3threads.h"
#include "idglobal_p3.h"
#include "clibtypes.h"
#include "gmslibname.h"
#include "xcompress.h"


Prototype Function(SYSTEM_integer) ( CDECL* XCOMPRESS_tcompress)(SYSTEM_pointer pdest, CLIBTYPES_clib_ulong* ldest, SYSTEM_pointer psrc,
      CLIBTYPES_clib_ulong lsrc);


Prototype Function(SYSTEM_integer) ( CDECL* XCOMPRESS_tuncompress)(SYSTEM_pointer pdest, CLIBTYPES_clib_ulong* ldest, SYSTEM_pointer psrc,
      CLIBTYPES_clib_ulong lsrc);


Prototype Function(XCOMPRESS_pgzfile) ( CDECL* XCOMPRESS_tgzreadopen)(SYSTEM_P3_pchar fn, SYSTEM_P3_pchar mode);


Prototype Function(SYSTEM_integer) ( CDECL* XCOMPRESS_tgzread)(XCOMPRESS_pgzfile pgz, SYSTEM_pointer pdest, SYSTEM_cardinal ldest);


Prototype Function(SYSTEM_integer) ( CDECL* XCOMPRESS_tgzreadclose)(XCOMPRESS_pgzfile pgz);

static P3LIBRARY_tlibhandle XCOMPRESS_zlibhandle;
static XCOMPRESS_tcompress XCOMPRESS_pcompress;
static XCOMPRESS_tuncompress XCOMPRESS_puncompress;
static XCOMPRESS_tgzreadopen XCOMPRESS_pgzreadopen;
static XCOMPRESS_tgzread XCOMPRESS_pgzread;
static XCOMPRESS_tgzreadclose XCOMPRESS_pgzreadclose;

Function(SYSTEM_boolean) XCOMPRESS_zlibdllloaded(void) {
   SYSTEM_boolean result;

   result = XCOMPRESS_zlibhandle != ValueCast(SYSTEM_pointer, 0);
   return result;
}  /* zlibdllloaded */

static Function(SYSTEM_pointer) loadentry(const SYSTEM_ansichar* n, SYSTEM_ansichar* _2wfn, SYSTEM_ansichar* _2loadmsg) {
   SYSTEM_pointer result;

   if (_P3strcmpN(_2loadmsg, _P3str1("\000"))) {
      result = NULL;
   }
   else {
      {
         SYSTEM_shortstring _t1;

         result = P3LIBRARY_p3getprocaddress(XCOMPRESS_zlibhandle, SYSUTILS_P3_lowercase(_t1, 255, n));
      }
      if (result == NULL) {
         _P3STR_255 _t1;
         _P3STR_255 _t2;

         _P3strcat(_2loadmsg, 255, _P3strcat(_t2, 255, _P3strcat(_t1, 255, _P3str1("\021Entry not found: "), n), _P3str1("\004 in ")), _2wfn);
      }
   }
   return result;
}  /* loadentry */

Function(SYSTEM_boolean) XCOMPRESS_loadzliblibrary(const SYSTEM_ansichar* fn, SYSTEM_ansichar* loadmsg) {
   SYSTEM_boolean result;
   SYSTEM_shortstring wfn, basename;
   SYSTEM_shortstring path;

   _P3strclr(loadmsg);
   if (XCOMPRESS_zlibhandle == ValueCast(SYSTEM_pointer, 0)) {
      {
         SYSTEM_shortstring _t1;

         _P3strcpy(path, 255, SYSUTILS_P3_extractfilepath(_t1, 255, fn));
      }
      {
         SYSTEM_shortstring _t1;

         _P3strcpy(basename, 255, SYSUTILS_P3_extractfilename(_t1, 255, fn));
      }
      if (_P3strcmpE(basename, _P3str1("\000")))
         _P3strcpy(basename, 255, _P3str1("\010gmszlib1"));
      {
         SYSTEM_shortstring _t1;

         _P3strcpy(wfn, 255, GMSLIBNAME_gamslibnamep3(_t1, 255, basename));
      }
      _P3strcat(wfn, 255, path, wfn);
      XCOMPRESS_zlibhandle = P3LIBRARY_p3loadlibrary(wfn, loadmsg);
      if (XCOMPRESS_zlibhandle != ValueCast(SYSTEM_pointer, 0) && _P3strcmpE(loadmsg, _P3str1("\000"))) {
         XCOMPRESS_pcompress = ValueCast(XCOMPRESS_tcompress, loadentry(_P3str1("\010compress"), wfn, loadmsg));
         XCOMPRESS_puncompress = ValueCast(XCOMPRESS_tuncompress, loadentry(_P3str1("\012uncompress"), wfn, loadmsg));
         XCOMPRESS_pgzreadopen = ValueCast(XCOMPRESS_tgzreadopen, loadentry(_P3str1("\006gzopen"), wfn, loadmsg));
         XCOMPRESS_pgzread = ValueCast(XCOMPRESS_tgzread, loadentry(_P3str1("\006gzread"), wfn, loadmsg));
         XCOMPRESS_pgzreadclose = ValueCast(XCOMPRESS_tgzreadclose, loadentry(_P3str1("\007gzclose"), wfn, loadmsg));
      }
   }
   if (_P3strcmpN(loadmsg, _P3str1("\000"))) {
      XCOMPRESS_pcompress = NULL;
      XCOMPRESS_puncompress = NULL;
      XCOMPRESS_pgzreadopen = NULL;
      XCOMPRESS_pgzread = NULL;
      XCOMPRESS_pgzreadclose = NULL;
   }
   result = ValueCast(SYSTEM_pointer, XCOMPRESS_pcompress) != NULL;
   return result;
}  /* loadzliblibrary */

Procedure XCOMPRESS_unloadzliblibrary(void) {
   if (XCOMPRESS_zlibhandle != ValueCast(SYSTEM_pointer, 0)) {
      P3LIBRARY_p3freelibrary(XCOMPRESS_zlibhandle);
      XCOMPRESS_zlibhandle = ValueCast(SYSTEM_pointer, 0);
   }
   XCOMPRESS_pcompress = NULL;
   XCOMPRESS_puncompress = NULL;
   XCOMPRESS_pgzreadopen = NULL;
   XCOMPRESS_pgzread = NULL;
   XCOMPRESS_pgzreadclose = NULL;
}  /* unloadzliblibrary */

Function(SYSTEM_integer) XCOMPRESS_compress(SYSTEM_pointer pdest, CLIBTYPES_clib_ulong* ldest, SYSTEM_pointer psrc, CLIBTYPES_clib_ulong lsrc) {
   SYSTEM_integer result;

   result = (*XCOMPRESS_pcompress)(pdest, ldest, psrc, lsrc);
   return result;
}  /* compress */

Function(SYSTEM_integer) XCOMPRESS_uncompress(SYSTEM_pointer pdest, CLIBTYPES_clib_ulong* ldest, SYSTEM_pointer psrc, CLIBTYPES_clib_ulong lsrc) {
   SYSTEM_integer result;

   result = (*XCOMPRESS_puncompress)(pdest, ldest, psrc, lsrc);
   return result;
}  /* uncompress */

Function(XCOMPRESS_pgzfile) XCOMPRESS_gzreadopen(const SYSTEM_ansichar* fn) {
   XCOMPRESS_pgzfile result;
   SYSTEM_shortstring sfn;
   SYSTEM_shortstring smode;

   _P3strcat(sfn, 255, fn, _P3str1("\001\000"));
   _P3strcpy(smode, 255, _P3str1("\003rb\000"));
   result = (*XCOMPRESS_pgzreadopen)(ValueCast(SYSTEM_P3_pansichar, &sfn[1]), ValueCast(SYSTEM_P3_pansichar, &smode[1]));
   return result;
}  /* gzreadopen */

Function(SYSTEM_integer) XCOMPRESS_gzread(XCOMPRESS_pgzfile pgz, SYSTEM_untyped* buf, SYSTEM_longword ldest) {
   SYSTEM_integer result;

   result = (*XCOMPRESS_pgzread)(pgz, ValueCast(SYSTEM_pointer, buf), ldest);
   return result;
}  /* gzread */

Function(SYSTEM_integer) XCOMPRESS_gzreadclose(XCOMPRESS_pgzfile* pgz) {
   SYSTEM_integer result;

   result = (*XCOMPRESS_pgzreadclose)(*pgz);
   *pgz = NULL;
   return result;
}  /* gzreadclose */

/* unit xcompress */
void _Init_Module_xcompress(void) {
   XCOMPRESS_zlibhandle = ValueCast(SYSTEM_pointer, 0);
   XCOMPRESS_unloadzliblibrary();
} /* _Init_Module_xcompress */

void _Final_Module_xcompress(void) {
} /* _Final_Module_xcompress */

