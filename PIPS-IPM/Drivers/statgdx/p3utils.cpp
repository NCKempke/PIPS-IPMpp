#include "p3io.h"
#include "system_p3.h"
#include "p3private.h"
#include "p3platform.h"
#include "exceptions.h"
#include "sysutils_p3.h"
#include "math_p3.h"
#include "p3library.h"
#include "p3utils.h"

/**** C code included from p3utils.pas(239:1): 59 lines ****/
#if defined(_WIN32)
# include <winsock2.h>
# include <io.h>
# include <psapi.h>  /* enough if we run on Windows 7 or later */
# include <iphlpapi.h>
# include <shlobj.h>

typedef BOOL (WINAPI * GetFileSizeEx_t) (HANDLE h, PLARGE_INTEGER fileSize);
GetFileSizeEx_t pGetFileSizeEx = NULL;
int triedGetFileSizeEx = 0;

typedef BOOL (WINAPI * SetFilePointerEx_t)
    (HANDLE h, LARGE_INTEGER distance,
     PLARGE_INTEGER newPointer, DWORD whence);
SetFilePointerEx_t pSetFilePointerEx = NULL;
int triedSetFilePointerEx = 0;

# if ! defined(_WIN64)

WINBASEAPI BOOL WINAPI
GetFileSizeEx (HANDLE h, PLARGE_INTEGER fileSize);
WINBASEAPI BOOL WINAPI
SetFilePointerEx (HANDLE h, LARGE_INTEGER distance,
                  PLARGE_INTEGER newPointer, DWORD whence);
# endif /* ! defined(_WIN64) */
#else

# include <fcntl.h>
# include <sys/types.h>
# include <sys/stat.h>

/* these next for socket commo */
# include <sys/socket.h>

# if (defined(__linux__) || defined(__APPLE__) || defined(__HOS_AIX__) || defined(__sparc) || defined(__sun__)) /* at least, maybe for others too */

#  include <arpa/inet.h>

#  if defined(__linux__)

#   include <sys/ioctl.h>
#   include <net/if.h>

#  elif defined(__APPLE__)
#   include <sys/ioctl.h>
#   include <sys/sysctl.h>
#   include <net/if.h>
#   include <net/if_dl.h>
#  endif
# endif

# include <netinet/in.h>
# include <netdb.h>

# include <sys/utsname.h>
# include <pwd.h>

# if defined(__APPLE__)
#  include <libproc.h>
# elif defined(__sparc)
#  include <procfs.h>
# endif

#endif

#include <locale.h>

Function(SYSTEM_integer) P3UTILS_p3chmod(const SYSTEM_ansichar* path, SYSTEM_integer mode) {
   SYSTEM_integer result;

   /**** C code included from p3utils.pas(310:1): 13 lines ****/
#if defined(_WIN32)
   result = 0;
#else
   {
      char filename[256];
      int len;
      /* */
      len = path[0];
      memcpy(filename, path + 1, len);
      filename[len] = '\0';
      result = chmod(filename, mode);
   }
#endif
   return result;
}  /* p3chmod */

Function(SYSTEM_double) P3UTILS_realtrunc(SYSTEM_double x) {
   SYSTEM_double result;

   result = SYSTEM_int(x);
   return result;
}  /* realtrunc */

Function(SYSTEM_double) P3UTILS_realround(SYSTEM_double x) {
   SYSTEM_double result;

   if (x >= 0) {
      result = SYSTEM_int(x + 0.5);
   }
   else
      result = SYSTEM_int(x - 0.5);
   return result;
}  /* realround */

static Function(SYSTEM_double) myroundto(SYSTEM_double x, SYSTEM_integer i) {
   SYSTEM_double result;
   SYSTEM_double z;

   if (i == 0) {
      if (x > 0) {
         result = SYSTEM_int(x + 0.5);
      }
      else
         result = SYSTEM_int(x - 0.5);
   }
   else if (i > 0) {
      z = MATH_P3_intpower(10, i);
      if (x > 0) {
         result = SYSTEM_int(x * z + 0.5) / z;
      }
      else
         result = SYSTEM_int(x * z - 0.5) / z;
   }
   else {
      z = MATH_P3_intpower(10, -i);
      if (x > 0) {
         result = SYSTEM_int(x / z + 0.5) * z;
      }
      else
         result = SYSTEM_int(x / z - 0.5) * z;
   }
   return result;
}  /* myroundto */

Function(SYSTEM_ansichar *) P3UTILS_floattoe(SYSTEM_ansichar* result, SYSTEM_uint8 _len_ret, SYSTEM_double y, SYSTEM_integer decimals) {
   SYSTEM_integer e, i, j, k, n;
   SYSTEM_double x;
   SYSTEM_shortstring s;

   x = SYSTEM_abs_r(y);
   n = 0;
   if (x != 0) {
      while (x >= 10.0) {
         n = n + 1;
         x = x / 10.0;
      }
      while (x < 1.0) {
         n = n - 1;
         x = x * 10.0;
      }
      x = myroundto(x, decimals);
      x = x * MATH_P3_intpower(10.0, n);
   }
   _P3str_d0(x, s, 255);
   k = SYSUTILS_P3_lastdelimiter(_P3str1("\002+-"), s);
   j = SYSTEM_pos(_P3str1("\001."), s);
   if (k - j - 2 < decimals)
      decimals = k - j - 2;
   _P3strcpy(result, _len_ret, _P3str1("\002  "));
   if (y < 0)
      result[2] = _P3char('-');
   {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;
      _P3STR_255 _t3;
      _P3STR_3 _t4;

      _P3strcat(result, _len_ret, _P3strcat(_t3, 255, _P3strcat(_t2, 255, result, SYSTEM_copy(_t1, 255, s, j - 1, decimals + 2)), _P3str1("\001E")),
            _P3ch2str(_t4, 1, s[k]));
   }
   {
      SYSTEM_shortstring _t1;

      _P3val_i(SYSTEM_copy(_t1, 255, s, k, 5), e, &i);
   }
   if (e < 0)
      e = -e;
   if (e > 99) {
      {
         SYSTEM_shortstring _t1;

         _P3strcat(result, _len_ret, result, VariableCast(SYSTEM_shortstring, SYSUTILS_P3_inttostr(_t1, 255, e), SYSTEM_shortstring));
      }
   }
   else {
      SYSTEM_shortstring _t1;

      _P3strcat(result, _len_ret, result, SYSTEM_copy(_t1, 255, s, ValueCast(SYSTEM_int32, SYSTEM_length(s)) - 1, 2));
   }
   return result;
}  /* floattoe */

Function(SYSTEM_ansichar *) P3UTILS_replacefileext(SYSTEM_ansichar* result, SYSTEM_uint8 _len_ret, const SYSTEM_ansichar* filename,
      const SYSTEM_ansichar* extension) {
   {
      SYSTEM_shortstring _t1;

      _P3strcpy(result, _len_ret, VariableCast(SYSTEM_shortstring, SYSUTILS_P3_changefileext(_t1, 255, filename, extension), SYSTEM_shortstring));
   }
   return result;
}  /* replacefileext */

Function(SYSTEM_ansichar *) P3UTILS_completefileext(SYSTEM_ansichar* result, SYSTEM_uint8 _len_ret, const SYSTEM_ansichar* filename,
      const SYSTEM_ansichar* extension) {
   {
      SYSTEM_shortstring _t1;

      if (_P3strcmpE(SYSUTILS_P3_extractfileext(_t1, 255, filename), _P3str1("\000"))) {
         {
            SYSTEM_shortstring _t1;

            _P3strcpy(result, _len_ret,
                  VariableCast(SYSTEM_shortstring, SYSUTILS_P3_changefileext(_t1, 255, filename, extension), SYSTEM_shortstring));
         }
      }
      else
         _P3strcpy(result, _len_ret, filename);
   }
   return result;
}  /* completefileext */

Function(SYSTEM_ansichar *) P3UTILS_paramstrzero(SYSTEM_ansichar* result, SYSTEM_uint8 _len_ret) {
   SYSTEM_P3_paramstr(result, _len_ret, 0);
   return result;
}  /* paramstrzero */
/**** C code included from p3utils.pas(451:1): 26 lines ****/

#define _EINVALIDCAST_RAISE_E(msg) _P3_RAISE(ValueCast(EXCEPTIONS_einvalidcast,SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_einvalidcast_CD)),msg)));

#if !defined(_WIN32)

static SYSTEM_boolean setEnvironmentVariable(const SYSTEM_char* name, const SYSTEM_char* value) {
   int rc;

   if (NULL == name || '\0' == *name) {
      return SYSTEM_false;
   }
   if (NULL == value) {          /* delete name from the environment */
      /* unsetenv will do the job here */
      unsetenv((const char*) name);
      rc = 0;
      return rc ? SYSTEM_false : SYSTEM_true;
   }

   rc = setenv((const char*) name, (const char*) value, 1);

   return rc ? SYSTEM_false : SYSTEM_true;
}                                       /* setEnvironmentVariable */
#endif /* ! defined(_WIN32) */

Function(SYSTEM_boolean) P3UTILS_prefixpath(const SYSTEM_ansichar* s) {
   SYSTEM_boolean result;
   static _P3STR_7 cpath = {5, 'P', 'A', 'T', 'H', '\000'};
   SYSTEM_integer slen, plen;
   SYSTEM_P3_pansichar tptr;

   slen = SYSTEM_length(s);
   if (slen == 0) {
      result = SYSTEM_true;
   }
   else {
      tptr = NULL;
      /**** C code included from p3utils.pas(526:1): 42 lines ****/
      {
         char* p;

#if defined(_WIN32)
         plen = GetEnvironmentVariable((char *)cpath+1,NULL,0);
         p = (char *) malloc (slen + 1 + plen);
         if (NULL == p) {
           return SYSTEM_false;
         }
         memcpy(p, (char *)s+1, slen);
         if (plen > 0) {
           int tlen;
           p[slen] = SYSUTILS_P3_pathsep;
           tlen = GetEnvironmentVariable((char *)cpath+1,p+slen+1,plen);
           assert(tlen == plen -1);
         }
         else {
           p[slen] = '\0';
         }
         result = SetEnvironmentVariable((char *)cpath+1,p);
         free(p);
#else
         tptr = (SYSTEM_char*) getenv((char*) cpath + 1);
         plen = 0;
         if (NULL != tptr)
            plen = strlen((char*) tptr);
         p = (char*) malloc(slen + 1 + plen + 1);
         if (NULL == p) {
            return SYSTEM_false;
         }
         memcpy(p, (char*) s + 1, slen);
         if (plen > 0) {
            p[slen] = SYSUTILS_P3_pathsep;
            memcpy(p + slen + 1, (char*) tptr, plen);
            p[slen + 1 + plen] = '\0';
         }
         else
            p[slen] = '\0';
         result = setEnvironmentVariable(cpath + 1, (SYSTEM_char*) p);
         free(p);
#endif /* #if defined(_WIN32) .. else .. */
      }
   }
   return result;
}  /* prefixpath */

Function(SYSTEM_ansichar *) P3UTILS_loadpathvarname(SYSTEM_ansichar* result, SYSTEM_uint8 _len_ret) {
   _P3strclr(result);
   switch (P3PLATFORM_osplatform()) {
      case P3PLATFORM_osaix:
         _P3strcpy(result, _len_ret, _P3str1("\007LIBPATH"));
         break;
      case P3PLATFORM_osdarwin_x64:
         _P3strcpy(result, _len_ret, _P3str1("\021DYLD_LIBRARY_PATH"));
         break;
      case P3PLATFORM_osbluegene:
      case P3PLATFORM_oslinux:
      case P3PLATFORM_oslinux86_64:
      case P3PLATFORM_ossunos_sparc32:
      case P3PLATFORM_ossunos_sparc64:
      case P3PLATFORM_ossunos_i86pc:
         _P3strcpy(result, _len_ret, _P3str1("\017LD_LIBRARY_PATH"));
         break;
      default:
         break;
   }
   return result;
}  /* loadpathvarname */

Function(SYSTEM_boolean) P3UTILS_prefixloadpath(const SYSTEM_ansichar* dir) {
   SYSTEM_boolean result;
   SYSTEM_integer slen, plen;
   SYSTEM_P3_pansichar tptr;
   SYSTEM_shortstring s;
   SYSTEM_shortstring ldpath;

   result = SYSTEM_false;
   slen = SYSTEM_length(dir);
   if (slen == 0) {
      {
         SYSTEM_shortstring _t1;
         SYSTEM_shortstring _t2;
         SYSTEM_shortstring _t3;

         _P3strcpy(s, 255, SYSUTILS_P3_excludetrailingpathdelimiter(_t1, 255, SYSUTILS_P3_extractfilepath(_t2, 255, P3UTILS_paramstrzero(_t3, 255))));
      }
      slen = SYSTEM_length(s);
   }
   else
      _P3strcpy(s, 255, dir);
   P3UTILS_loadpathvarname(ldpath, 255);
   if (0 == SYSTEM_length(ldpath)) {
      result = SYSTEM_true;
      return result;
   }
   _P3strcat(ldpath, 255, ldpath, _P3str1("\001\000"));
   tptr = NULL;
   /**** C code included from p3utils.pas(629:1): 24 lines ****/
#if !defined (_WIN32)
   {
      char* p;

      tptr = (SYSTEM_char*) getenv((char*) ldpath + 1);
      plen = 0;
      if (NULL != tptr)
         plen = strlen((char*) tptr);
      p = (char*) malloc(slen + 1 + plen + 1);
      if (NULL == p) {
         return SYSTEM_false;
      }
      memcpy(p, (char*) s + 1, slen);
      if (plen > 0) {
         p[slen] = SYSUTILS_P3_pathsep;
         memcpy(p + slen + 1, (char*) tptr, plen);
         p[slen + 1 + plen] = '\0';
      }
      else
         p[slen] = '\0';
      result = setEnvironmentVariable(ldpath + 1, (SYSTEM_char*) p);
      free(p);
   }
#endif /* if ! defined (_WIN32) */
   return result;
}  /* prefixloadpath */

Function(SYSTEM_boolean) P3UTILS_prefixenv(const SYSTEM_ansichar* dir, const SYSTEM_ansichar* evname) {
   SYSTEM_boolean result;
   SYSTEM_shortstring trimdir;
   P3PRIVATE_shortstrbuf dirbuf, evnamebuf;
   SYSTEM_P3_pansichar dirptr, evnameptr;
   SYSTEM_cardinal dlen;
   SYSTEM_P3_pansichar tptr;

   {
      SYSTEM_shortstring _t1;

      _P3strcpy(trimdir, 255, SYSUTILS_P3_trim(_t1, 255, dir));
   }
   dlen = SYSTEM_length(trimdir);
   if (0 == dlen) {
      result = SYSTEM_true;
      return result;
   }
   dirptr = P3PRIVATE_strtostrbuf(trimdir, dirbuf);
   evnameptr = P3PRIVATE_strtostrbuf(evname, evnamebuf);
   tptr = NULL;
   /**** C code included from p3utils.pas(751:1): 78 lines ****/
   {
      char* p;
      size_t evLen;

#if defined(_WIN32)
      char *oldPtr, *op;
      size_t evSiz;
      evSiz = GetEnvironmentVariable((const char *)evnameptr,NULL,0);
      if (0 == evSiz) {  /* not set: just set it */
        result = SetEnvironmentVariable((const char *)evnameptr,(const char *)dirptr);
        return result;
      }
      oldPtr = (char *) malloc (evSiz);
      if (NULL == oldPtr)
        return SYSTEM_false;
      evLen = GetEnvironmentVariable((const char *)evnameptr,oldPtr,evSiz);
      assert((evSiz-1)==evLen);
      op = oldPtr;
      while (*op) {
        if (! isspace(*op))
          break;
        op++, evLen--;
      }

      if ('\0' == *op) {  /* if empty, just set it */
        result = SetEnvironmentVariable((const char *)evnameptr,(const char *)dirptr);
        free(oldPtr);
        return result;
      }
      /* check if dir is the first element */
      if ((evLen >= (size_t)dlen) &&
          (0==strncmp(op,(char *)dirptr,dlen)) &&
          ( ('\0' == op[dlen]) || (SYSUTILS_P3_pathsep == op[dlen]))
         ) { /* evName already starts with dir */
        free(oldPtr);
        return SYSTEM_true;
      }
      p = (char *) malloc (dlen + 1 + evLen + 1);
      if (NULL == p) {
        free(oldPtr);
        return SYSTEM_false;
      }
      memcpy(p,(char *)dirptr,dlen);
      p[dlen] = SYSUTILS_P3_pathsep;
      memcpy(p+dlen+1,op,evLen);
      p[dlen+1+evLen] = '\0';
      result = SetEnvironmentVariable((const char *)evnameptr,p);
      free(oldPtr);
      free(p);

#else

      tptr = (SYSTEM_char*) getenv((char*) evnameptr);
      evLen = 0;
      if (NULL != tptr)
         evLen = strlen((char*) tptr);
      if (0 == evLen) { /* env var is not set or empty */
         result = setEnvironmentVariable(evnameptr, dirptr);
      }
      else {
         if (((ssize_t) evLen >= dlen) && (0 == strncmp((char*) tptr, (char*) dirptr, dlen)) &&
             (('\0' == tptr[dlen]) || (SYSUTILS_P3_pathsep == tptr[dlen]))) /* evName already starts with dir */
            return SYSTEM_true;
         p = (char*) malloc(dlen + 1 + evLen + 1);
         if (NULL == p) {
            return SYSTEM_false;
         }
         memcpy(p, (char*) dirptr, dlen);
         p[dlen] = SYSUTILS_P3_pathsep;
         memcpy(p + dlen + 1, (char*) tptr, evLen);
         p[dlen + 1 + evLen] = '\0';
         result = setEnvironmentVariable(evnameptr, (const SYSTEM_char*) p);
         free(p);
      }
#endif /* #if defined(_WIN32) .. else .. */
   }
   return result;
}  /* prefixenv */

Function(SYSTEM_boolean) P3UTILS_p3setenv(const SYSTEM_ansichar* name, const SYSTEM_ansichar* val) {
   SYSTEM_boolean result;
   P3PRIVATE_shortstrbuf namebuf, valbuf;
   SYSTEM_P3_pansichar nameptr, valptr;

   nameptr = P3PRIVATE_strtostrbuf(name, namebuf);
   valptr = P3PRIVATE_strtostrbuf(val, valbuf);
   /**** C code included from p3utils.pas(850:1): 5 lines ****/
#if defined(_WIN32)
   result = SetEnvironmentVariable((char *)nameptr,(char *)valptr);
#else
   result = setEnvironmentVariable(nameptr, valptr);
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* p3setenv */

Procedure P3UTILS_p3unsetenv(const SYSTEM_ansichar* name) {
   P3PRIVATE_shortstrbuf namebuf;
   SYSTEM_P3_pansichar nameptr;

   nameptr = P3PRIVATE_strtostrbuf(name, namebuf);
   /**** C code included from p3utils.pas(876:1): 5 lines ****/
#if defined(_WIN32)
   (void) SetEnvironmentVariable((char *)nameptr,NULL);
#else
   (void) setEnvironmentVariable(nameptr, NULL);
#endif /* if defined(_WIN32) .. else .. */
}  /* p3unsetenv */

Function(SYSTEM_boolean) P3UTILS_p3issetenv(const SYSTEM_ansichar* name) {
   SYSTEM_boolean result;
   P3PRIVATE_shortstrbuf namebuf;
   SYSTEM_P3_pansichar nameptr;

   nameptr = P3PRIVATE_strtostrbuf(name, namebuf);
   /**** C code included from p3utils.pas(902:1): 5 lines ****/
#if defined(_WIN32)
   result = (GetEnvironmentVariable((char *)nameptr,NULL,0) != 0);
#else
   result = (getenv((char*) nameptr) != NULL);
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* p3issetenv */

Function(SYSTEM_boolean) P3UTILS_p3setenvpc(const SYSTEM_ansichar* name, SYSTEM_P3_pansichar val) {
   SYSTEM_boolean result;
   P3PRIVATE_shortstrbuf namebuf;
   SYSTEM_P3_pansichar nameptr;

   nameptr = P3PRIVATE_strtostrbuf(name, namebuf);
   /**** C code included from p3utils.pas(927:1): 5 lines ****/
#if defined(_WIN32)
   result = SetEnvironmentVariable((char *)nameptr,(char *)val);
#else
   result = setEnvironmentVariable(nameptr, val);
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* p3setenvpc */

Function(SYSTEM_cardinal) P3UTILS_p3getenvpc(const SYSTEM_ansichar* name, SYSTEM_P3_pansichar buf, SYSTEM_cardinal bufsiz) {
   SYSTEM_cardinal result;
   P3PRIVATE_shortstrbuf namebuf;
   SYSTEM_P3_pansichar nameptr;

   nameptr = P3PRIVATE_strtostrbuf(name, namebuf);
   /**** C code included from p3utils.pas(952:1): 21 lines ****/
#if defined(_WIN32)
   result = GetEnvironmentVariableA((char *)nameptr,(char *)buf, bufsiz);
#else
   {
      const char* p;

      p = getenv((char*) nameptr);
      if (NULL == p) { /* no match in the environment */
         result = 0;
      }
      else {
         size_t psiz = strlen(p) + 1;
         if (psiz <= bufsiz) {  /* it fits: copy it over */
            (void) memmove((char*) buf, p, psiz);
            result = psiz - 1;
         }
         else
            result = psiz;
      }
   }
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* p3getenvpc */

Procedure P3UTILS_p3setconsoletitle(const SYSTEM_ansichar* s) {
   P3PRIVATE_shortstrbuf namebuf;
   SYSTEM_P3_pansichar nameptr;

   nameptr = P3PRIVATE_strtostrbuf(s, namebuf);
   /**** C code included from p3utils.pas(992:1): 5 lines ****/
#if defined(_WIN32)
   SetConsoleTitle((char *)nameptr);
#else
   /* do nothing for now */
#endif /* if defined(_WIN32) .. else .. */
}  /* p3setconsoletitle */

Procedure P3UTILS_p3nopopups(void) {
   /**** C code included from p3utils.pas(1013:1): 5 lines ****/
#if defined(_WIN32)
   SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
#else
   /* do nothing */
#endif /* if defined(_WIN32) .. else .. */
}  /* p3nopopups */
/**** C code included from p3utils.pas(1024:1): 62 lines ****/
#if defined(_WIN32)
/* map the Windows codes returned by GetLastError to libc codes
 * use when making Windows API calls with P3, since P3 ioresult-ish codes
 * are expected to be libc codes on all platforms
 */
static int
win2c (int rc)
{
  int result;

  switch (rc) {
    case ERROR_FILE_NOT_FOUND:
    case ERROR_PATH_NOT_FOUND:
      result = ENOENT;
      break;
    case ERROR_TOO_MANY_OPEN_FILES:
      result = EMFILE;
      break;
    case ERROR_ACCESS_DENIED:
      result = EACCES;
      break;
    case ERROR_INVALID_HANDLE:
      result = EBADF;
      break;
    case ERROR_NOT_ENOUGH_MEMORY:
      result = ENOMEM;
      break;
    case ERROR_INVALID_ACCESS:
      result = EACCES;
      break;
    case ERROR_NO_MORE_FILES:
      result = ENFILE;
      break;
    case ERROR_SEEK_ON_DEVICE:
      result = ESPIPE;
      break;
    case ERROR_INVALID_PARAMETER:
    case ERROR_NEGATIVE_SEEK:
      result = EINVAL;
      break;
    default:
      result = 0; /* no guessing */
  }              /* case */
  return result;
} /* win2c */

static const DWORD accessMode[3] = {
  GENERIC_READ,
  GENERIC_WRITE,
  GENERIC_READ | GENERIC_WRITE };
/* this works for GDX so we do it: it is kind of silly to use the
 * shareMode var then but why not? */
static const DWORD shareMode[3] = {
  FILE_SHARE_READ | FILE_SHARE_WRITE,
  FILE_SHARE_READ | FILE_SHARE_WRITE,
  FILE_SHARE_READ | FILE_SHARE_WRITE };
static const DWORD createHow[3] = {
  OPEN_EXISTING,
  CREATE_ALWAYS,
  OPEN_ALWAYS };

#endif /* if defined(_WIN32) */

static Function(SYSTEM_boolean) P3UTILS_isvalidhandle(P3UTILS_tp3filehandle h) {
   SYSTEM_boolean result;

   /**** C code included from p3utils.pas(1096:1): 5 lines ****/
#if defined(_WIN32)
   result = h && (INVALID_HANDLE_VALUE != (HANDLE) h);
#else
   result = (SYSTEM_int64) h > 0;
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* isvalidhandle */

Function(SYSTEM_integer) P3UTILS_p3fileopen(const SYSTEM_ansichar* fname, P3UTILS_tp3fileopenaction mode, P3UTILS_tp3filehandle* h) {
   SYSTEM_integer result;
   P3PRIVATE_shortstrbuf namebuf;
   SYSTEM_P3_pansichar nameptr;

   nameptr = P3PRIVATE_strtostrbuf(fname, namebuf);
   /**** C code included from p3utils.pas(1178:1): 81 lines ****/
#if defined(_WIN32)
   {
     DWORD lowMode;
     HANDLE hFile;
     int rc;

     lowMode = mode & 3;
     if (3 == lowMode) {
       *h = (P3UTILS_tp3filehandle) INVALID_HANDLE_VALUE;
       return ERROR_INVALID_PARAMETER;
     }

     if (nameptr[0] == 0) {
        if (P3UTILS_p3openread == mode)
           hFile = GetStdHandle(STD_INPUT_HANDLE);
        else if (P3UTILS_p3openwrite == mode)
           hFile = GetStdHandle(STD_OUTPUT_HANDLE);
        else {
           *h = (P3UTILS_tp3filehandle) INVALID_HANDLE_VALUE;
           return ERROR_INVALID_PARAMETER;
        }
     }
     else
       hFile = CreateFile ((char *)nameptr, accessMode[lowMode], shareMode[lowMode], NULL,
                         createHow[lowMode], FILE_ATTRIBUTE_NORMAL, NULL);
     if (INVALID_HANDLE_VALUE == hFile) {
       *h = (P3UTILS_tp3filehandle) INVALID_HANDLE_VALUE;
       rc = GetLastError();
       result = win2c(rc);
       if (0 == result) { /* ouch: just pick a likely but non-specific code */
         result = EACCES;
       }
     }
     else {
       *h = (P3UTILS_tp3filehandle) hFile;
       result = 0;
     }
   }
#else
   {
      struct stat statBuf;
      int fd, flags, rc;

      if (nameptr[0] == 0) {
         if (P3UTILS_p3openread == mode)
            *h = (P3UTILS_tp3filehandle) STDIN_FILENO;
         else if (P3UTILS_p3openwrite == mode)
            *h = (P3UTILS_tp3filehandle) STDOUT_FILENO;
         else {
            *h = (P3UTILS_tp3filehandle) 0;
            return -1;
         }
         return 0;
      }
      flags = mode & 3;
      if (flags > 0)            /* write-only or read-write */
         flags |= O_CREAT;
      if (flags & 1)
         flags |= O_TRUNC;
      fd = open((char*) nameptr, flags, 0666);
      if (-1 == fd) {
         *h = (P3UTILS_tp3filehandle) 0;
         return errno;
      }
      result = 0;
      /* before calling this a success, check for directory on read-only */
      if (P3UTILS_p3openread == mode) {
         rc = fstat(fd, &statBuf);
         if (rc)
            result = errno;
         else if (S_ISDIR(statBuf.st_mode))
            result = EISDIR;
      }
      if (result) {
         close(fd);
         return result;
      }

      *h = (P3UTILS_tp3filehandle) (SYSTEM_nativeuint) fd;
   }
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* p3fileopen */

Function(SYSTEM_integer) P3UTILS_p3fileclose(P3UTILS_tp3filehandle* h) {
   SYSTEM_integer result;

   /**** C code included from p3utils.pas(1280:1): 28 lines ****/
   if (!P3UTILS_isvalidhandle(h)) {
      return EBADF;
   }
#if defined(_WIN32)
   {
     if (CloseHandle((HANDLE)*h)) {
       result = 0; /* success */
     }
     else {
       int rc = GetLastError();
       result = win2c(rc);
       if (0 == result) { /* ouch: just pick a likely but non-specific code */
         result = EIO;
       }
     }
     *h = (P3UTILS_tp3filehandle) INVALID_HANDLE_VALUE;
   }
#else
   {
      SYSTEM_int64 h64;

      h64 = (SYSTEM_int64) *h;
      result = 0;
      if (close((int) h64))   /* error!! */
         result = errno;
      *h = (P3UTILS_tp3filehandle) 0;
   }
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* p3fileclose */

Function(SYSTEM_integer) P3UTILS_p3fileread(P3UTILS_tp3filehandle h, SYSTEM_untyped* buffer, SYSTEM_longword buflen, SYSTEM_longword* numread) {
   SYSTEM_integer result;

   /**** C code included from p3utils.pas(1325:1): 31 lines ****/
#if defined(_WIN32)
   {
     if (ReadFile((HANDLE)h, buffer, buflen, (LPDWORD)numread, NULL)) {
       result = 0; /* success */
     }
     else {
       int rc = GetLastError();
       result = win2c(rc);
       if (0 == result) { /* ouch: just pick a likely but non-specific code */
         result = EIO;
       }
     }
   }
#else
   {
      int rc;
      SYSTEM_int64 h64;

      h64 = (SYSTEM_int64) h;
      result = 0;
      rc = read((int) h64, buffer, buflen);
      if (rc < 0) {
         result = errno;
         *numread = 0;
      }
      else {
         result = 0;
         *numread = (SYSTEM_longword) rc;
      }
   }
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* p3fileread */

Function(SYSTEM_integer)
P3UTILS_p3filewrite(P3UTILS_tp3filehandle h, const SYSTEM_untyped* buffer, SYSTEM_longword buflen, SYSTEM_longword* numwritten) {
   SYSTEM_integer result;

   /**** C code included from p3utils.pas(1373:1): 31 lines ****/
#if defined(_WIN32)
   {
     if (WriteFile((HANDLE)h, buffer, buflen, (LPDWORD)numwritten, NULL)) {
       result = 0; /* success */
     }
     else {
       int rc = GetLastError();
       result = win2c(rc);
       if (0 == result) { /* ouch: just pick a likely but non-specific code */
         result = EIO;
       }
     }
   }
#else
   {
      int rc;
      SYSTEM_int64 h64;

      h64 = (SYSTEM_int64) h;
      result = 0;
      rc = write((int) h64, buffer, buflen);
      if (rc < 0) {
         result = errno;
         *numwritten = 0;
      }
      else {
         result = 0;
         *numwritten = (SYSTEM_longword) rc;
      }
   }
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* p3filewrite */

Function(SYSTEM_integer) P3UTILS_p3filegetsize(P3UTILS_tp3filehandle h, SYSTEM_int64* filesize) {
   SYSTEM_integer result;

   *filesize = -1;
   /**** C code included from p3utils.pas(1427:1): 49 lines ****/
   if (!P3UTILS_isvalidhandle(h)) {
      return EBADF;
   }
#if defined(_WIN32)
   {
     BOOL frc;

     if (! triedGetFileSizeEx) {
       pGetFileSizeEx = (GetFileSizeEx_t) GetProcAddress(
         GetModuleHandle("kernel32"),"GetFileSizeEx");
       triedGetFileSizeEx = 1;
     }
     if (pGetFileSizeEx) {
       frc = pGetFileSizeEx((HANDLE)h, (PLARGE_INTEGER) filesize);
     }
     else {
       DWORD tt;

       frc = GetFileSize((HANDLE)h, &tt);
       *filesize = tt;
     }

     if (frc)
       result = 0;
     else {
       int rc = GetLastError();
       result = win2c(rc);
       if (0 == result) { /* ouch: just pick a likely but non-specific code */
         result = EACCES;
       }
     }
   }
#else
   {
      struct stat statBuf;
      int rc;
      SYSTEM_int64 h64;

      h64 = (SYSTEM_int64) h;
      rc = fstat((int) h64, &statBuf);
      if (rc) {
         result = errno;
      }
      else {
         result = 0;
         *filesize = statBuf.st_size;
      }
   }
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* p3filegetsize */

Function(SYSTEM_integer) P3UTILS_p3filesetpointer(P3UTILS_tp3filehandle h, SYSTEM_int64 distance, SYSTEM_int64* newpointer, SYSTEM_longword whence) {
   SYSTEM_integer result;

   /**** C code included from p3utils.pas(1499:1): 75 lines ****/
   if (!P3UTILS_isvalidhandle(h)) {
      return EBADF;
   }
#if defined(_WIN32)
   {
     int rc;
     BOOL frc;
     LARGE_INTEGER d; /* declared as a union - compiler rejects a cast */

     if (! triedSetFilePointerEx) {
       pSetFilePointerEx = (SetFilePointerEx_t) GetProcAddress(
         GetModuleHandle("kernel32"),"SetFilePointerEx");
       triedSetFilePointerEx = 1;
     }
     d.QuadPart = distance;
     if (pSetFilePointerEx) {
       frc = pSetFilePointerEx((HANDLE)h, d, (PLARGE_INTEGER) newpointer, whence);
     }
     else {
       DWORD tt, trc;;

       trc = SetFilePointer((HANDLE)h, (int)distance, NULL, whence);
       if (INVALID_SET_FILE_POINTER == trc) {
         frc = 0;
         *newpointer = 0;
       }
       else {
         *newpointer = trc;
         frc = 1;
       }
     }
     if (frc)
       result = 0;
     else {
       rc = GetLastError();
       result = win2c(rc);
       if (0 == result) {
         result = EINVAL;
       }
     }
   }
#else
   {
      int w = -1;
      off_t newPos, offset;
      SYSTEM_int64 h64;

      switch (whence) {
         case P3UTILS_p3_file_begin:
            w = SEEK_SET;
            break;
         case P3UTILS_p3_file_current:
            w = SEEK_CUR;
            break;
         case P3UTILS_p3_file_end:
            w = SEEK_END;
            break;
         default:
            return EINVAL;
      } /* switch */

      /* check if conversion to off_t loses info */
      offset = (off_t) distance;
#if !defined(__alpha) /* no overflow possible on OSF, so not defined */
      if (offset != distance)
         return EOVERFLOW;
#endif
      h64 = (SYSTEM_int64) h;
      newPos = lseek((int) h64, offset, w);
      if ((off_t) -1 == newPos)
         return errno;
      *newpointer = newPos;
      result = 0;
   }
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* p3filesetpointer */

Function(SYSTEM_integer) P3UTILS_p3filegetpointer(P3UTILS_tp3filehandle h, SYSTEM_int64* filepointer) {
   SYSTEM_integer result;

   /**** C code included from p3utils.pas(1594:1): 56 lines ****/
   if (!P3UTILS_isvalidhandle(h)) {
      return EBADF;
   }
#if defined(_WIN32)
   {
     int rc;
     BOOL frc;
     LARGE_INTEGER d; /* declared as a union - compiler rejects a cast */

     if (! triedSetFilePointerEx) {
       pSetFilePointerEx = (SetFilePointerEx_t) GetProcAddress(
         GetModuleHandle("kernel32"),"SetFilePointerEx");
       triedSetFilePointerEx = 1;
     }
     d.QuadPart = 0;
     if (pSetFilePointerEx) {
       frc = pSetFilePointerEx((HANDLE)h, d, (PLARGE_INTEGER) filepointer,
                               P3UTILS_p3_file_current);
     }
     else {
       DWORD tt, trc;;

       trc = SetFilePointer((HANDLE)h, 0, NULL, P3UTILS_p3_file_current);
       if (INVALID_SET_FILE_POINTER == trc) {
         frc = 0;
         *filepointer = 0;
       }
       else {
         *filepointer = trc;
         frc = 1;
       }
     }

     if (frc)
       result = 0;
     else {
       rc = GetLastError();
       result = win2c(rc);
       if (0 == result) {
         result = EINVAL;
       }
     }
   }
#else
   {
      off_t newPos;
      SYSTEM_int64 h64;

      h64 = (SYSTEM_int64) h;
      newPos = lseek((int) h64, 0, SEEK_CUR);
      if ((off_t) -1 == newPos)
         return errno;
      *filepointer = newPos;
      result = 0;
   }
#endif /* if defined(_WIN32) .. else .. */
   return result;
}  /* p3filegetpointer */

Function(SYSTEM_int64) P3UTILS_p3allocmemsize64(void) {
   SYSTEM_int64 result;

   /**** C code included from p3utils.pas(1662:1): 1 lines ****/
   result = SYSTEM_allocmemsize64;
   return result;
}  /* p3allocmemsize64 */

Function(SYSTEM_pointer) P3UTILS_p3allocmem64(SYSTEM_int64 size) {
   SYSTEM_pointer result;

   result = NULL;
   if (size <= 0)
      return result;
   /**** C code included from p3utils.pas(1683:1): 13 lines ****/
   if (8 == sizeof(result)) {
      _P3_new64(&result, size);
      (void) memset(result, 0, (size_t) size);
   }
   else {
      if (size > SYSTEM_maxint) {
         _P3_Exception(_P3_EXC_CODE_OUTOFMEMORY, "");
      }
      else {
         _P3_new(&result, (SYSTEM_longint) size);
         (void) memset(result, 0, (size_t) size);
      }
   }
   return result;
}  /* p3allocmem64 */

Procedure P3UTILS_p3fillchar64(SYSTEM_untyped* p, SYSTEM_int64 size, SYSTEM_byte fillvalue) {
   /**** C code included from p3utils.pas(1728:1): 11 lines ****/
   if (size <= 0)
      return;
   if (8 == sizeof(p)) {
      (void) memset(p, (int) fillvalue, (size_t) size);
   }
   else {
      size_t sz = (size_t) size;
      if ((SYSTEM_int64) sz != size)
         _EINVALIDCAST_RAISE_E(_P3str1("\034p3FillChar64: size too large"));
      (void) memset(p, (int) fillvalue, sz);
   }
}  /* p3fillchar64 */

Procedure P3UTILS_p3freemem64(SYSTEM_pointer* p, SYSTEM_int64 size) {
   /**** C code included from p3utils.pas(1753:1): 11 lines ****/
   if (8 == sizeof(p)) {
      _P3_free64(*p, size);
   }
   else {
      if (size < 0 || size > SYSTEM_maxint) {
         _P3_free(*p, 0);
      }
      else {
         _P3_free(*p, (SYSTEM_longint) size);
      }
   }
}  /* p3freemem64 */

Procedure P3UTILS_p3getmem64(SYSTEM_pointer* p, SYSTEM_int64 size) {
   if (size <= 0) {
      *p = NULL;
      return;
   }
   /**** C code included from p3utils.pas(1782:1): 11 lines ****/
   if (8 == sizeof(*p)) {
      _P3_new64(p, size);
   }
   else {
      if (size > SYSTEM_maxint) {
         _P3_Exception(_P3_EXC_CODE_OUTOFMEMORY, "");
      }
      else {
         _P3_new(p, (SYSTEM_longint) size);
      }
   }
}  /* p3getmem64 */

Procedure P3UTILS_p3move64(const SYSTEM_untyped* source, SYSTEM_untyped* dest, SYSTEM_int64 sz) {
   /**** C code included from p3utils.pas(1812:1): 9 lines ****/
   if (8 == sizeof(source)) {
      (void) memmove(dest, source, (size_t) sz);
   }
   else {
      size_t ssz = (size_t) sz;
      if ((SYSTEM_int64) ssz != sz)
         _EINVALIDCAST_RAISE_E(_P3str1("\030p3Move64: size too large"));
      (void) memmove(dest, source, ssz);
   }
}  /* p3move64 */

Procedure P3UTILS_p3reallocmem64(SYSTEM_pointer* p, SYSTEM_int64 size) {
   if (size <= 0)
      if (*p == NULL) {
         return;
      }
      else {
         _P3freemem0(*p);
         *p = NULL;
      }
   /**** C code included from p3utils.pas(1848:1): 11 lines ****/
   if (8 == sizeof(p)) {
      SYSTEM_reallocmem64(p, size);
   }
   else {
      if (size > SYSTEM_maxint) {
         _P3_Exception(_P3_EXC_CODE_OUTOFMEMORY, "");
      }
      else {
         SYSTEM_reallocmem(p, (SYSTEM_longint) size);
      }
   }
}  /* p3reallocmem64 */

Procedure P3UTILS_p3getfromurl(const SYSTEM_ansichar* servername, const SYSTEM_ansichar* filename, SYSTEM_word port, P3UTILS_thavedatacb havedata,
      SYSTEM_pointer usermem, SYSTEM_ansichar* msg) {
   cnstdef { maxbuf = 4096 };
   typedef SYSTEM_uint16 _sub_1P3GETFROMURL;
   typedef SYSTEM_ansichar _arr_0P3GETFROMURL[4096];
   _arr_0P3GETFROMURL buffer;
   SYSTEM_integer len;

   _P3strclr(msg);
   _P3strcpy(msg, 255, _P3str1("\041Not implemented for this platform"));
   /**** C code included from p3utils.pas(1964:1): 115 lines ****/
#if defined(_WIN32)
#define WSACLEANUP   WSACancelBlockingCall(); WSACleanup(); return
   {
     WSADATA wsaData;
     SOCKET conn;
     struct hostent *host;
     char sbuf[300];
     unsigned int addr;
     struct sockaddr_in server;
     int rc;

     rc = WSAStartup (0x101, &wsaData);
     if (rc) {
       strcpy ((char *)msg, "\026Winsock startup failed");
       return;
     }
     conn = socket (AF_INET, SOCK_STREAM, IPPROTO_TCP);
     if (conn < 0) {
       strcpy ((char *)msg, "\017No socket found");
       WSACLEANUP;
     }
     len = servername[0];
     strncpy (sbuf, (char *)servername+1, len);
     sbuf[len] = '\0';
     addr = inet_addr(sbuf);
     if (INADDR_NONE == addr)
       host = gethostbyname (sbuf);
     else
       host = gethostbyaddr ((char *) &addr, sizeof(addr), PF_INET);
     if (NULL == host) {
       closesocket (conn);
       strcpy ((char *)msg, "\017Unknown address");
       WSACLEANUP;
     }
     memcpy (&server.sin_addr.s_addr, host->h_addr, host->h_length);
     server.sin_family=AF_INET;
     server.sin_port=htons(port);
     if (connect (conn, (struct sockaddr *) &server, sizeof(server))) {
       closesocket(conn);
       strcpy ((char *)msg, "\021Could not connect");
       WSACLEANUP;
     }
     sprintf(sbuf,"GET /%.*s HTTP/1.0\r\n\r\n", *filename, filename+1);
     len = send (conn, sbuf, sizeof(sbuf), 0);
     for ( ; ; ) {
       len = recv (conn, (char *)buffer, sizeof(buffer), 0);
       if (len < 0) {
         closesocket(conn);
         strcpy ((char *)msg, "\024Error receiving data");
         WSACLEANUP;
       }
       if (0 == len)
         break;
       if (! havedata(buffer, len, usermem)) {
         strcpy ((char *)msg, "\023Stopped by callback");
         WSACLEANUP;
       }
     }  /* receive loop */
     closesocket(conn);
     *msg = '\0';
   }
#else
#if !defined(INADDR_NONE)
# define INADDR_NONE 0xffffffff
#endif
   {
      struct hostent* host;
      in_addr_t addr;
      struct sockaddr_in server;
      int conn;
      char sbuf[300];

      conn = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
      if (conn < 0) {
         strcpy((char*) msg, "\017No socket found");
         return;
      }
      len = servername[0];
      strncpy(sbuf, (char*) servername + 1, len);
      sbuf[len] = '\0';
      addr = inet_addr(sbuf);
      if (INADDR_NONE == addr)
         host = gethostbyname(sbuf);
      else
         host = gethostbyaddr((void*) &addr, sizeof(addr), PF_INET);
      if (NULL == host) {
         strcpy((char*) msg, "\017Unknown address");
         return;
      }
      memcpy(&server.sin_addr.s_addr, host->h_addr_list[0], host->h_length);
      server.sin_family = AF_INET;
      server.sin_port = htons(port);
      if (connect(conn, (struct sockaddr*) &server, sizeof(server))) {
         strcpy((char*) msg, "\021Could not connect");
         return;
      }
      sprintf(sbuf, "GET /%.*s HTTP/1.0\r\n\r\n", *filename, filename + 1);
      len = send(conn, sbuf, sizeof(sbuf), 0);
      for (;;) {
         len = recv(conn, buffer, sizeof(buffer), 0);
         if (len < 0) {
            strcpy((char*) msg, "\024Error receiving data");
            return;
         }
         if (0 == len)
            break;
         if (!havedata(buffer, len, usermem)) {
            strcpy((char*) msg, "\023Stopped by callback");
            return;
         }
      }  /* receive loop */
      (void) shutdown(conn, SHUT_RDWR);
      *msg = '\0';
   }
#endif /* if defined(_WIN32) .. else .. */
}  /* p3getfromurl */

Function(SYSTEM_boolean) P3UTILS_p3getfirstmacaddress(SYSTEM_ansichar* mac) {
   SYSTEM_boolean result;

   result = SYSTEM_false;
   /**** C code included from p3utils.pas(2093:1): 167 lines ****/
#if defined(__linux__)
   {
      struct ifreq ifr;
      struct ifconf ifc;
      char buf[1024];
      int success = 0, sock;

      sock = socket(AF_INET, SOCK_DGRAM, IPPROTO_IP);
      if (sock == -1)
         return result;

      ifc.ifc_len = sizeof(buf);
      ifc.ifc_buf = buf;
      if (ioctl(sock, SIOCGIFCONF, &ifc) == -1)
         return result;

      {
         struct ifreq* it = ifc.ifc_req;
         const struct ifreq* const end = it + (ifc.ifc_len / sizeof(struct ifreq));

         for (; it != end; ++it) {
            strcpy(ifr.ifr_name, it->ifr_name);
            if (ioctl(sock, SIOCGIFFLAGS, &ifr) == 0) {
               if (ifr.ifr_flags & IFF_LOOPBACK) /* don't count loopback */
                  continue;
               if (ioctl(sock, SIOCGIFHWADDR, &ifr) == 0) {
                  success = 1;
                  break;
               }
            }
         } /* loop over interfaces */
      }

      if (success) {
         unsigned char mb[6];
         memcpy(mb, ifr.ifr_hwaddr.sa_data, 6);
         sprintf((char*) mac + 1, "%02x:%02x:%02x:%02x:%02x:%02x", mb[0], mb[1], mb[2], mb[3], mb[4], mb[5]);
         mac[0] = 17;
         result = SYSTEM_true;
      }
   } /* if __linux__ */

#elif defined(__APPLE__)
   {
     char prevName[IF_NAMESIZE];
     int mib[6];
     int sock;
     int halfDone = 0;  /* true if we have a MAC number for an interface that is down */
     struct ifconf ifc;
     char buf[1024];
     char buf2[1024];
     unsigned char *mp;
     struct ifreq ifr;
     struct ifreq *it, *end;
     size_t recLen, sz;
     struct if_msghdr *ifm;
     struct sockaddr_dl *sdl;

     sock = socket(PF_INET, SOCK_DGRAM, IPPROTO_IP);
     if (sock < 0)
       return SYSTEM_false;
     ifc.ifc_len = sizeof(buf);
     ifc.ifc_buf = buf;
     if (ioctl(sock, SIOCGIFCONF, &ifc))
       return SYSTEM_false;
#if 0
     if ((ifc.ifc_len + sizeof(struct ifreq)) > sizeof(buf))
       return SYSTEM_false;  /* we should have used a bigger buffer */
#endif
     it = ifc.ifc_req;
     end = (struct ifreq *) ((unsigned char *) ifc.ifc_buf + ifc.ifc_len);

     mib[0] = CTL_NET;
     mib[1] = AF_ROUTE;
     mib[2] = 0;
     mib[3] = AF_LINK;
     mib[3] = 0;
     mib[4] = NET_RT_IFLIST;
     prevName[0] = '\0';
     for ( ;  it < end;  it = (struct ifreq *) ((unsigned char *) it + recLen)) {
       recLen = _SIZEOF_ADDR_IFREQ(*it);
       if (0==strcmp(it->ifr_name, prevName)) /* just checked it already */
         continue;
       (void) strcpy (prevName, it->ifr_name);
       (void) strcpy (ifr.ifr_name, it->ifr_name);
       if (ioctl(sock, SIOCGIFFLAGS, &ifr))  /* we should always get flags but if not skip ahead */
         continue;
       if (ifr.ifr_flags & IFF_LOOPBACK) /* always skip loopback interfaces */
         continue;
       if (halfDone && (0 == (ifr.ifr_flags & IFF_UP) ) )
         continue;  /* we already have a MAC address for a down interface */
       mib[5] = if_nametoindex(it->ifr_name);
       if (0 == mib[5])
         continue;      /* no valid index found */
       sz = sizeof(buf2);
       if (sysctl(mib, 6, buf2, &sz, NULL, 0))
         continue;     /* sysctl call failed */
       ifm = (struct if_msghdr *) buf2;
       /* printf ("msglen 0 = %d\n", ifm->ifm_msglen); */
       sdl = (struct sockaddr_dl *) (ifm +1);
       if (RTM_IFINFO != ifm->ifm_type)
         continue;     /* WTF */
       mp = (unsigned char *) LLADDR(sdl);
       sprintf((char*) mac+1,"%02x:%02x:%02x:%02x:%02x:%02x",
               mp[0], mp[1], mp[2], mp[3], mp[4], mp[5]);
       mac[0] = 17;
       if (0 != (ifr.ifr_flags & IFF_UP) )
         return SYSTEM_true;
       else
         halfDone = 1;
     } /* loop over interfaces */

   } /* if __APPLE__ */

#elif defined(_WIN32)
   {
     ULONG bufSiz, prevBufSiz;
     ULONG flags = GAA_FLAG_INCLUDE_PREFIX;
     prevBufSiz = bufSiz = 4096;
     int nTries = 0, maxTries = 3;
     PIP_ADAPTER_ADDRESSES addrBuf = NULL;
     PIP_ADAPTER_ADDRESSES currAddr;
     DWORD dwrc, ifType;
     unsigned char *mp;
     int halfDone = 0; /* if we have a MAC number for an interface that is down */
     do {
       addrBuf = (IP_ADAPTER_ADDRESSES *) malloc(bufSiz);
       if (NULL == addrBuf)
         return SYSTEM_false;
       dwrc = GetAdaptersAddresses(AF_INET, flags, NULL, addrBuf, &bufSiz);
       if (ERROR_BUFFER_OVERFLOW == dwrc) {
         prevBufSiz = bufSiz;
         free(addrBuf);
         addrBuf = NULL;
       }
       nTries++;
     } while ((ERROR_BUFFER_OVERFLOW == dwrc) && (nTries < maxTries));
     if (NO_ERROR != dwrc) {
       if (addrBuf)
         free(addrBuf);
       return SYSTEM_false;
     }
     for (currAddr = addrBuf;  currAddr;  currAddr = currAddr->Next) {
       ifType = currAddr->IfType;
       if ( (IF_TYPE_ETHERNET_CSMACD != ifType) &&
            (IF_TYPE_IEEE80211 != ifType) )
         continue;
       if (halfDone && (IfOperStatusUp != currAddr->OperStatus))
         continue;  /* we already have a MAC address for a down interface */
       if (6 != currAddr->PhysicalAddressLength)
         continue;
       mp = (unsigned char *) currAddr->PhysicalAddress;
       sprintf((char*) mac+1,"%02x:%02x:%02x:%02x:%02x:%02x",
               mp[0], mp[1], mp[2], mp[3], mp[4], mp[5]);
       mac[0] = 17;
       if (IfOperStatusUp == currAddr->OperStatus) {
         free(addrBuf);
         return SYSTEM_true;
       }
       else
         halfDone = 1;
     }
     free(addrBuf);
   } /* if _WIN32 */

#endif
   return result;
}  /* p3getfirstmacaddress */

Function(SYSTEM_integer) P3UTILS_p3getwindowsversion(void) {
   SYSTEM_integer result;
   SYSTEM_pointer libhandle;
   SYSTEM_shortstring loadmsg;
   SYSTEM_pointer wine_get_version;

   result = P3UTILS_os_unknown;
   libhandle = P3LIBRARY_p3loadlibrary(_P3str1("\011ntdll.dll"), loadmsg);
   if (libhandle != NULL) {
      result = P3UTILS_os_winnt;
      wine_get_version = P3LIBRARY_p3getprocaddress(libhandle, _P3str1("\020wine_get_version"));
      if (wine_get_version != NULL)
         result = P3UTILS_os_wine;
      P3LIBRARY_p3freelibrary(libhandle);
   }
   else {
      libhandle = P3LIBRARY_p3loadlibrary(_P3str1("\014kernel32.dll"), loadmsg);
      if (libhandle != NULL) {
         result = P3UTILS_os_win9x;
         wine_get_version = P3LIBRARY_p3getprocaddress(libhandle, _P3str1("\021GetDKrnl32Version"));
         if (wine_get_version != NULL)
            result = P3UTILS_os_hx;
         P3LIBRARY_p3freelibrary(libhandle);
      }
   }
   return result;
}  /* p3getwindowsversion */

Function(SYSTEM_ansichar *) P3UTILS_p3getcomputername(SYSTEM_ansichar* result, SYSTEM_uint8 _len_ret) {
   _P3strcpy(result, _len_ret, _P3str1("\007unknown"));
   /**** C code included from p3utils.pas(2316:1): 27 lines ****/
#if defined(_WIN32)
   {
     char compName[256];
     DWORD n;

     n = (DWORD) sizeof(compName);
     if (GetComputerName(compName, &n)) {
       /* success, copy compName to result */
       *result = (SYSTEM_byte) n;
       (void) memcpy ((char *)(result+1), compName, n);
     }
   }
#else
   {
      int rc, n;
      struct utsname uts;

      rc = uname(&uts);
      if (rc >= 0) {
         n = strlen(uts.nodename);
         if (n > 255)
            n = 255;
         *result = (SYSTEM_byte) n;
         (void) memcpy((char*) (result + 1), uts.nodename, n);
      }
   }
#endif /* #if defined(_WIN32) .. else .. */
   return result;
}  /* p3getcomputername */

Function(SYSTEM_ansichar *) P3UTILS_p3getusername(SYSTEM_ansichar* result, SYSTEM_uint8 _len_ret) {
   _P3strcpy(result, _len_ret, _P3str1("\007unknown"));
   /**** C code included from p3utils.pas(2366:1): 47 lines ****/
#if defined(_WIN32)
   {
     char userName[256];
     DWORD n;

     n = (DWORD) sizeof(userName);
     if (GetUserName(userName, &n)) {
       /* success, copy userName to result */
       n--;  /* output n includes space for the null byte */
       *result = (SYSTEM_byte) n;
       (void) memcpy ((char *)(result+1), userName, n);
     }
   }
#else
   {
      int n;
      char loginName[256];
      char* p = NULL;

# if (defined(__USE_XOPEN) || defined(_XOPEN_SOURCE) || defined(_POSIX_PTHREAD_SEMANTICS) || defined(__EXTENSIONS__)) && !defined(__APPLE__)
      /* the cuserid routine is preferred */
      p = cuserid(loginName);
      if (p) {
         loginName[sizeof(loginName) - 1] = '\0';
         n = strlen(loginName);
         *result = (SYSTEM_byte) n;
         (void) memcpy((char*) (result + 1), loginName, n);
      }
# else
      {
       int rc;
#  if (defined(SOL) && ! defined(_POSIX_PTHREAD_SEMANTICS))
       p = getlogin_r (loginName, sizeof(loginName));
       rc = (NULL == p);
#  else
       rc = getlogin_r (loginName, sizeof(loginName));
#  endif
       if (0 == rc) { /* success */
         loginName[sizeof(loginName)-1] = '\0';
         n = strlen(loginName);
         *result = (SYSTEM_byte) n;
         (void) memcpy ((char *)(result+1), loginName, n);
       }
      }
# endif
   }
#endif /* #if defined(_WIN32) .. else .. */
   return result;
}  /* p3getusername */

Function(SYSTEM_boolean) P3UTILS_p3senddatamessage(SYSTEM_boolean broadcast, const SYSTEM_ansichar* _ftmp1, const SYSTEM_ansichar* _ftmp2) {
   SYSTEM_shortstring wintitle;
   SYSTEM_shortstring data;
   SYSTEM_boolean result;

   _P3strcpy(wintitle, 255, _ftmp1);
   _P3strcpy(data, 255, _ftmp2);
   result = SYSTEM_false;
   /**** C code included from p3utils.pas(2455:1): 38 lines ****/
#if defined(_WIN32)
   {
      HWND receiver;
      COPYDATASTRUCT cds;
      LRESULT r;
      char dBuf[256], tBuf[256];
      unsigned int n;

      n = data[0];
      if (n > 255)
         n = 255;
      (void) memcpy (dBuf, (char *)data+1, n);
      dBuf[n] = '\0';
      cds.dwData = 0;
      cds.lpData = dBuf;
      cds.cbData = n+1; /* include the null byte in this count */
      r = 0;
      if (broadcast) {
         r = SendMessage (HWND_BROADCAST, (UINT)WM_COPYDATA, (WPARAM)0,
                          (LPARAM)&cds);
      }
      else {
         n = wintitle[0];
         if (n > 255)
            n = 255;
         if (n > 0) {
            (void) memcpy (tBuf, (char *)wintitle+1, n);
            tBuf[n] = '\0';
            receiver = FindWindow (NULL, tBuf);
            if (receiver)
               r = SendMessage (receiver, (UINT)WM_COPYDATA, (WPARAM)0,
                                (LPARAM)&cds);
         }
      }
      result = (SYSTEM_boolean) r;
   }
#else
#endif /* #if defined(_WIN32) .. else .. */
   return result;
}  /* p3senddatamessage */

Function(SYSTEM_ansichar *) P3UTILS_p3pushdeflocale(SYSTEM_ansichar* result, SYSTEM_uint8 _len_ret) {
   _P3strcpy(result, _len_ret, _P3str1("\001C"));
   /**** C code included from p3utils.pas(2502:1): 10 lines ****/
   {
      char* s;
      s = setlocale(LC_NUMERIC, NULL);
      if (('C' != s[0]) || ('\0' != s[1])) { /* not what we want */
         /* so save the current locale and reset */
         (void) strcpy((char*) result + 1, s);
         *result = strlen(s);
         (void) setlocale(LC_NUMERIC, "C");
      }
   }
   return result;
}  /* p3pushdeflocale */

Procedure P3UTILS_p3popdeflocale(const SYSTEM_ansichar* s) {
   if (_P3char('\001') == s[0] && _P3char('C') == s[1])
      return;
   /**** C code included from p3utils.pas(2525:1): 10 lines ****/
   {
      char buf[32];
      int n = sizeof(buf) - 1;

      if (*s < n)
         n = *s;
      (void) memcpy(buf, s + 1, n);
      buf[n] = '\0';
      (void) setlocale(LC_NUMERIC, buf);
   }
}  /* p3popdeflocale */
/**** C code included from p3utils.pas(2543:1): 1 lines ****/
#include "p3Custom2.h"

Function(SYSTEM_integer) P3UTILS_p3getexecname(SYSTEM_ansichar* execname, SYSTEM_ansichar* msg) {
   SYSTEM_integer result;

   result = 9;
   _P3strclr(execname);
   _P3strcpy(msg, 255, _P3str1("\027P3: not yet implemented"));
   /**** C code included from p3utils.pas(2556:1): 1 lines ****/
   result = xGetExecName(execname, msg);
   return result;
}  /* p3getexecname */

Function(SYSTEM_integer) P3UTILS_p3getlibname(SYSTEM_ansichar* libname, SYSTEM_ansichar* msg) {
   SYSTEM_integer result;

   if (!SYSTEM_islibrary()) {
      result = 2;
      _P3strclr(libname);
      _P3strcpy(msg, 255, _P3str1("\031Not called from a library"));
      return result;
   }
   result = 9;
   _P3strclr(libname);
   _P3strcpy(msg, 255, _P3str1("\027P3: not yet implemented"));
   /**** C code included from p3utils.pas(2598:1): 1 lines ****/
   result = xGetLibName(libname, msg);
   return result;
}  /* p3getlibname */

Function(SYSTEM_integer) P3UTILS_p3someioresult(void) {
   SYSTEM_integer result;

   /**** C code included from p3utils.pas(2636:1): 1 lines ****/
   result = EIO;
   return result;
}  /* p3someioresult */

Function(SYSTEM_boolean) P3UTILS_p3getmemoryinfo(SYSTEM_int64* rss, SYSTEM_int64* vss) {
   SYSTEM_boolean result;

   *rss = 0;
   *vss = 0;
   /**** C code included from p3utils.pas(2711:1): 68 lines ****/
#if defined(_WIN32)
   {
     PROCESS_MEMORY_COUNTERS info;
     int ok;
     ok = GetProcessMemoryInfo (GetCurrentProcess( ), &info, sizeof(info));
     if (!ok)
       return SYSTEM_false;  /* failure */
     *rss = (SYSTEM_int64) info.PagefileUsage;
     *vss = (SYSTEM_int64) info.WorkingSetSize;
     result = SYSTEM_true; /* success */
   }

#elif defined(__linux)
   {
      FILE* fp;
      int n;
      unsigned long urss, uvss;
      size_t sz;

      fp = fopen("/proc/self/statm", "r");
      if (!fp)
         return SYSTEM_false;  /* failure */
      /* first two are VmSize, VmRSS */
      n = fscanf(fp, "%lu %lu", &uvss, &urss);
      fclose(fp);
      if (2 != n)
         return SYSTEM_false;  /* failure */
      sz = sysconf(_SC_PAGESIZE);
      *rss = sz * urss;
      *vss = sz * uvss;
      result = SYSTEM_true; /* success */
   }

#elif defined(__APPLE__)
   {
     int ret;
     struct proc_taskinfo procTaskInfo;

     ret = proc_pidinfo ((int) getpid(), PROC_PIDTASKINFO, 0,
                          (void *) &procTaskInfo, sizeof(procTaskInfo));
     if (ret < (int)sizeof(procTaskInfo))
       return SYSTEM_false;  /* failure */
     *rss = (SYSTEM_int64) procTaskInfo.pti_resident_size;
     *vss = (SYSTEM_int64) procTaskInfo.pti_virtual_size;
     result = SYSTEM_true; /* success */
   }

#elif defined(__sparc)
   {
     struct psinfo psinfo;
     int fd = -1;
     ssize_t n;

     if ( (fd = open( "/proc/self/psinfo", O_RDONLY)) == -1)
       return SYSTEM_false;  /* failure */
     n = read (fd, &psinfo, sizeof(psinfo));
     close (fd);
     if (sizeof(psinfo) != n)
       return SYSTEM_false;  /* failure */
     *rss = (SYSTEM_int64) psinfo.pr_rssize * 1024;
     *vss = (SYSTEM_int64) psinfo.pr_size   * 1024;
     result = SYSTEM_true; /* success */
   }

#else
     result = SYSTEM_false; /* fail */

#endif
   return result;
}  /* p3getmemoryinfo */

Function(SYSTEM_boolean) P3UTILS_p3getmemoryinfoex(SYSTEM_cardinal pid, SYSTEM_int64* rss, SYSTEM_int64* vss) {
   SYSTEM_boolean result;

   *rss = 0;
   *vss = 0;
   /**** C code included from p3utils.pas(2821:1): 87 lines ****/
#if defined(_WIN32)
   {
     PROCESS_MEMORY_COUNTERS info;
     int ok;
     HANDLE h;
     if ((SYSTEM_cardinal) P3UTILS_pid_self == pid)
       h = GetCurrentProcess();
     else {
       h = OpenProcess (PROCESS_ALL_ACCESS, FALSE, pid);
       if (NULL == h)
         return SYSTEM_false;  /* failure */
     }
     ok = GetProcessMemoryInfo (h, &info, sizeof(info));
     if (!ok)
       return SYSTEM_false;  /* failure */
     *rss = (SYSTEM_int64) info.PagefileUsage;
     *vss = (SYSTEM_int64) info.WorkingSetSize;
     result = SYSTEM_true; /* success */
   }

#elif defined(__linux)
   {
      FILE* fp;
      char buf[32];
      unsigned long urss, uvss;
      size_t sz;
      int n;

      if ((SYSTEM_cardinal) P3UTILS_pid_self == pid)
         strcpy(buf, "/proc/self/statm");
      else
         snprintf(buf, sizeof(buf), "/proc/%d/statm", (int) pid);
      fp = fopen(buf, "r");
      if (!fp)
         return SYSTEM_false;  /* failure */
      /* first two are VmSize, VmRSS */
      n = fscanf(fp, "%lu %lu", &uvss, &urss);
      fclose(fp);
      if (2 != n)
         return SYSTEM_false;  /* failure */
      sz = sysconf(_SC_PAGESIZE);
      *rss = sz * urss;
      *vss = sz * uvss;
      result = SYSTEM_true; /* success */
   }

#elif defined(__APPLE__)
   {
     int ret, iPid;
     struct proc_taskinfo procTaskInfo;

     if ((SYSTEM_cardinal) P3UTILS_pid_self == pid)
       iPid = (int) getpid();
     else {
       iPid = (int) pid;
     }
     ret = proc_pidinfo (iPid, PROC_PIDTASKINFO, 0,
                         (void *) &procTaskInfo, sizeof(procTaskInfo));
     if (ret < (int)sizeof(procTaskInfo))
       return SYSTEM_false;  /* failure */
     *rss = (SYSTEM_int64) procTaskInfo.pti_resident_size;
     *vss = (SYSTEM_int64) procTaskInfo.pti_virtual_size;
     result = SYSTEM_true; /* success */
   }

#elif defined(__sparc)
   {
     struct psinfo psinfo;
     int fd = -1;
     ssize_t n;

     return SYSTEM_false;  /* failure */
     if ( (fd = open( "/proc/self/psinfo", O_RDONLY)) == -1)
       return SYSTEM_false;  /* failure */
     n = read (fd, &psinfo, sizeof(psinfo));
     close (fd);
     if (sizeof(psinfo) != n)
       return SYSTEM_false;  /* failure */
     *rss = (SYSTEM_int64) psinfo.pr_rssize * 1024;
     *vss = (SYSTEM_int64) psinfo.pr_size   * 1024;
     result = SYSTEM_true; /* success */
   }

#else
     result = SYSTEM_false; /* fail */

#endif
   return result;
}  /* p3getmemoryinfoex */

static Function(SYSTEM_boolean) P3UTILS_homeplus(const SYSTEM_ansichar* dd1, const SYSTEM_ansichar* dd2, SYSTEM_ansichar* s) {
   SYSTEM_boolean result;
   SYSTEM_cardinal len, bufsiz;
   typedef SYSTEM_uint8 _sub_1HOMEPLUS;
   typedef SYSTEM_ansichar _arr_0HOMEPLUS[256];
   _arr_0HOMEPLUS buf;

   bufsiz = sizeof(_arr_0HOMEPLUS);
   result = SYSTEM_false;
   len = P3UTILS_p3getenvpc(_P3str1("\004HOME"), ValueCast(SYSTEM_P3_pansichar, &buf[0]), bufsiz);
   if (0 == len)
      return result;
   if (len >= sizeof(SYSTEM_shortstring))
      return result;
   SYSTEM_move(&buf[0], &s[1], len);
   _P3setlength(s, len, 255);
   if (SYSTEM_length(dd1) > 0) {
      _P3inc1(len, SYSTEM_length(dd1));
      if (len >= sizeof(SYSTEM_shortstring))
         return result;
      _P3strcat(s, 255, s, dd1);
   }
   if (SYSTEM_length(dd2) > 0) {
      _P3inc1(len, SYSTEM_length(dd2));
      if (len >= sizeof(SYSTEM_shortstring))
         return result;
      _P3strcat(s, 255, s, dd2);
   }
   result = SYSTEM_true;
   return result;
}  /* homeplus */

Function(SYSTEM_boolean) P3UTILS_p3writablelocation(P3UTILS_tp3location loctype, const SYSTEM_ansichar* appname, SYSTEM_ansichar* locname) {
   SYSTEM_boolean result;
   SYSTEM_shortstring dd;
   SYSTEM_cardinal len, bufsiz;
   typedef SYSTEM_uint8 _sub_1P3WRITABLELOCATION;
   typedef SYSTEM_ansichar _arr_0P3WRITABLELOCATION[256];
   _arr_0P3WRITABLELOCATION buf;

   _P3strclr(locname);
   result = SYSTEM_false;
   bufsiz = sizeof(_arr_0P3WRITABLELOCATION);
   if (P3PLATFORM_osfiletype() == P3PLATFORM_osfilewin) { ;
   }
   else if (P3PLATFORM_osplatform() == P3PLATFORM_osdarwin_x64) {
      if (P3UTILS_p3config == loctype) {
         if (!P3UTILS_homeplus(_P3str1("\024/Library/Preferences"), _P3str1("\000"), locname))
            return result;
      }
      else if (P3UTILS_p3appconfig == loctype) {
         _P3strclr(dd);
         if (SYSTEM_length(appname) > 0) {
            _P3STR_3 _t1;

            _P3strcat(dd, 255, _P3ch2str(_t1, 1, SYSUTILS_P3_pathdelim), appname);
         }
         if (!P3UTILS_homeplus(_P3str1("\024/Library/Preferences"), dd, locname))
            return result;
      }
      else if (_P3SET_in_1(loctype, P3UTILS_p3data, _P3SET_in_1(loctype, P3UTILS_p3appdata, _P3SET_equal(loctype, P3UTILS_p3applocaldata)))) {
         _P3strclr(dd);
         if (SYSTEM_length(appname) > 0) {
            _P3STR_3 _t1;

            _P3strcat(dd, 255, _P3ch2str(_t1, 1, SYSUTILS_P3_pathdelim), appname);
         }
         if (!P3UTILS_homeplus(_P3str1("\034/Library/Application Support"), dd, locname))
            return result;
      }
      else if (P3UTILS_p3documents == loctype) {
         if (!P3UTILS_homeplus(_P3str1("\012/Documents"), _P3str1("\000"), locname))
            return result;
      }
      else
         return result;
      result = SYSTEM_true;
      return result;
   }
   else {
      if (P3UTILS_p3config == loctype) {
         len = P3UTILS_p3getenvpc(_P3str1("\017XDG_CONFIG_HOME"), ValueCast(SYSTEM_P3_pansichar, &buf[0]), bufsiz);
         if (len >= sizeof(SYSTEM_shortstring))
            return result;
         if (len > 0) {
            SYSTEM_move(&buf[0], &locname[1], len);
            _P3setlength(locname, len, 255);
         }
         else if (!P3UTILS_homeplus(_P3str1("\010/.config"), _P3str1("\000"), locname))
            return result;
      }
      else if (P3UTILS_p3appconfig == loctype) {
         len = P3UTILS_p3getenvpc(_P3str1("\017XDG_CONFIG_HOME"), ValueCast(SYSTEM_P3_pansichar, &buf[0]), bufsiz);
         if (len >= sizeof(SYSTEM_shortstring))
            return result;
         if (len > 0) {
            SYSTEM_move(&buf[0], &locname[1], len);
            _P3setlength(locname, len, 255);
         }
         else {
            _P3strclr(dd);
            if (SYSTEM_length(appname) > 0) {
               _P3STR_3 _t1;

               _P3strcat(dd, 255, _P3ch2str(_t1, 1, SYSUTILS_P3_pathdelim), appname);
            }
            if (!P3UTILS_homeplus(_P3str1("\010/.config"), dd, locname))
               return result;
         }
      }
      else if (_P3SET_in_1(loctype, P3UTILS_p3data, _P3SET_in_1(loctype, P3UTILS_p3appdata, _P3SET_equal(loctype, P3UTILS_p3applocaldata)))) {
         len = P3UTILS_p3getenvpc(_P3str1("\015XDG_DATA_HOME"), ValueCast(SYSTEM_P3_pansichar, &buf[0]), bufsiz);
         if (len >= sizeof(SYSTEM_shortstring))
            return result;
         if (len > 0) {
            SYSTEM_move(&buf[0], &locname[1], len);
            _P3setlength(locname, len, 255);
         }
         else {
            _P3strclr(dd);
            if (SYSTEM_length(appname) > 0) {
               _P3STR_3 _t1;

               _P3strcat(dd, 255, _P3ch2str(_t1, 1, SYSUTILS_P3_pathdelim), appname);
            }
            if (!P3UTILS_homeplus(_P3str1("\015/.local/share"), dd, locname))
               return result;
         }
      }
      else if (P3UTILS_p3documents == loctype) {
         if (!P3UTILS_homeplus(_P3str1("\012/Documents"), _P3str1("\000"), locname))
            return result;
      }
      else
         return result;
      result = SYSTEM_true;
      return result;
   }
   /**** C code included from p3utils.pas(3117:1): 38 lines ****/
#if defined(_WIN32)
   {
     size_t bufLen;
     DWORD r;
     HRESULT h;
     TCHAR folderName[MAX_PATH];

     if ((P3UTILS_p3config == loctype) || (P3UTILS_p3appconfig == loctype)
          || (P3UTILS_p3data == loctype) || (P3UTILS_p3applocaldata == loctype)) {
       r = GetEnvironmentVariableA("LOCALAPPDATA", (char *)folderName, MAX_PATH);
       if ((r > 0) && (r < sizeof(SYSTEM_shortstring))) { /* found it */
         result = SYSTEM_true;
         *locname = (SYSTEM_byte) r;
         (void) memcpy ((char *)(locname+1), folderName, r);
       }
     }
     else if (P3UTILS_p3appdata == loctype) {
       r = GetEnvironmentVariableA("APPDATA", (char *)folderName, MAX_PATH);
       if ((r > 0) && (r < sizeof(SYSTEM_shortstring))) { /* found it */
         result = SYSTEM_true;
         *locname = (SYSTEM_byte) r;
         (void) memcpy ((char *)(locname+1), folderName, r);
       }
     }
     else if (P3UTILS_p3documents == loctype) {
       h = SHGetFolderPathA (NULL, CSIDL_PERSONAL, NULL, 0, folderName);
       if (S_OK == h) {
         bufLen = strlen(folderName);
         if (bufLen >= sizeof(SYSTEM_shortstring))
           return result;
         result = SYSTEM_true;
         *locname = (SYSTEM_byte) bufLen;
         (void) memcpy ((char *)(locname+1), folderName, bufLen);
       }
     }
   } /* if _WIN32 */

#endif
   if (SYSTEM_length(appname) > 0 && _P3SET_in_1(loctype, P3UTILS_p3config, _P3SET_in_1(loctype, P3UTILS_p3appconfig,
         _P3SET_in_1(loctype, P3UTILS_p3data, _P3SET_in_1(loctype, P3UTILS_p3appdata, _P3SET_equal(loctype, P3UTILS_p3applocaldata)))))) {
      _P3STR_3 _t1;
      _P3STR_255 _t2;

      _P3strcat(locname, 255, _P3strcat(_t2, 255, locname, _P3ch2str(_t1, 1, SYSUTILS_P3_pathdelim)), appname);
   }
   return result;
}  /* p3writablelocation */

Function(SYSTEM_boolean)
P3UTILS_p3standardlocations(P3UTILS_tp3location loctype, const SYSTEM_ansichar* appname, SYSTEM_integer* loccount, SYSTEM_shortstring* locnames,
      SYSTEM_integer* ecount) {
   SYSTEM_boolean result;
   SYSTEM_integer rc, n;
   SYSTEM_cardinal k, dpos;
   SYSTEM_cardinal buflen, bufsiz;
   SYSTEM_shortstring execname, execpath, msg;
   typedef SYSTEM_uint16 _sub_1P3STANDARDLOCATIONS;
   typedef SYSTEM_ansichar _arr_0P3STANDARDLOCATIONS[1024];
   _arr_0P3STANDARDLOCATIONS buf;

   *loccount = 0;
   *ecount = 0;
   result = P3UTILS_p3writablelocation(loctype, appname, locnames[0]);
   if (result)
      _P3inc0(*loccount);
   if (P3UTILS_p3documents == loctype)
      return result;
   if (P3PLATFORM_osfiletype() == P3PLATFORM_osfilewin) {
      if (_P3SET_in_1(loctype, P3UTILS_p3config, _P3SET_equal(loctype, P3UTILS_p3appconfig))) {
         _P3inc0(*loccount);
         if (SYSTEM_length(appname) > 0) {
            {
               _P3STR_3 _t1;
               _P3STR_15 _t2;

               _P3strcat(locnames[*loccount - 1], 255, _P3strcat(_t2, 15, _P3str1("\016C:\\ProgramData"), _P3ch2str(_t1, 1, SYSUTILS_P3_pathdelim)),
                     appname);
            }
         }
         else
            _P3strcpy(locnames[*loccount - 1], 255, _P3str1("\016C:\\ProgramData"));
      }
      else if (_P3SET_in_1(loctype, P3UTILS_p3data, _P3SET_in_1(loctype, P3UTILS_p3appdata, _P3SET_equal(loctype, P3UTILS_p3applocaldata)))) {
         _P3inc0(*loccount);
         if (SYSTEM_length(appname) > 0) {
            {
               _P3STR_3 _t1;
               _P3STR_15 _t2;

               _P3strcat(locnames[*loccount - 1], 255, _P3strcat(_t2, 15, _P3str1("\016C:\\ProgramData"), _P3ch2str(_t1, 1, SYSUTILS_P3_pathdelim)),
                     appname);
            }
         }
         else
            _P3strcpy(locnames[*loccount - 1], 255, _P3str1("\016C:\\ProgramData"));
         rc = P3UTILS_p3getexecname(execname, msg);
         if (0 != rc) {
            _P3inc0(*ecount);
            return result;
         }
         {
            SYSTEM_shortstring _t1;

            _P3strcpy(execpath, 255, SYSUTILS_P3_extractfilepath(_t1, 255, execname));
         }
         _P3inc0(*loccount);
         {
            SYSTEM_shortstring _t1;

            _P3strcpy(locnames[*loccount - 1], 255, SYSUTILS_P3_excludetrailingpathdelimiter(_t1, 255, execpath));
         }
         _P3inc0(*loccount);
         _P3strcat(locnames[*loccount - 1], 255, execpath, _P3str1("\004data"));
         if (SYSTEM_length(appname) > 0) {
            _P3inc0(*loccount);
            {
               _P3STR_3 _t1;
               _P3STR_255 _t2;

               _P3strcat(locnames[*loccount - 1], 255, _P3strcat(_t2, 255, locnames[*loccount - 1 - 1], _P3ch2str(_t1, 1, SYSUTILS_P3_pathdelim)),
                     appname);
            }
         }
      }
   }
   else if (P3PLATFORM_osplatform() == P3PLATFORM_osdarwin_x64) {
      if (_P3SET_in_1(loctype, P3UTILS_p3data, _P3SET_in_1(loctype, P3UTILS_p3appdata, _P3SET_equal(loctype, P3UTILS_p3applocaldata)))) {
         _P3strcpy(msg, 255, _P3str1("\034/Library/Application Support"));
         if (SYSTEM_length(appname) > 0) {
            _P3STR_255 _t1;

            _P3strcat(msg, 255, _P3strcat(_t1, 255, msg, _P3str1("\001/")), appname);
         }
         _P3inc0(*loccount);
         _P3strcpy(locnames[*loccount - 1], 255, msg);
         rc = P3UTILS_p3getexecname(execname, msg);
         if (0 != rc) {
            _P3inc0(*ecount);
            return result;
         }
         {
            SYSTEM_shortstring _t1;
            SYSTEM_shortstring _t2;

            _P3strcpy(execpath, 255, SYSUTILS_P3_excludetrailingpathdelimiter(_t1, 255, SYSUTILS_P3_extractfilepath(_t2, 255, execname)));
         }
         if (SYSUTILS_P3_lastdelimiter(_P3str1("\001/"), execpath) >= 2) {
            _P3inc0(*loccount);
            {
               SYSTEM_shortstring _t1;

               _P3strcat(locnames[*loccount - 1], 255, SYSUTILS_P3_extractfilepath(_t1, 255, execpath), _P3str1("\011Resources"));
            }
         }
         else
            _P3inc0(*ecount);
      }
   }
   else {
      bufsiz = sizeof(_arr_0P3STANDARDLOCATIONS);
      if (P3UTILS_p3config == loctype) {
         buflen = P3UTILS_p3getenvpc(_P3str1("\017XDG_CONFIG_DIRS"), ValueCast(SYSTEM_P3_pansichar, &buf[0]), bufsiz);
         if (buflen >= bufsiz) {
            _P3inc0(*ecount);
            return result;
         }
         else if (buflen > 0) {
            k = 0;
            dpos = 0;
            do {
               while (buf[dpos] != _P3char('\000') && buf[dpos] != _P3char(':')) {

                  _P3inc0(dpos);
               }
               n = dpos - k;
               if (n >= sizeof(SYSTEM_shortstring)) {
                  _P3inc0(*ecount);
               }
               else if (n > 0)
                  if (*loccount >= 8) {
                     _P3inc0(*ecount);
                  }
                  else {
                     _P3inc0(*loccount);
                     SYSTEM_move(&buf[k], &locnames[*loccount - 1][1], n);
                     _P3setlength(locnames[*loccount - 1], n, 255);
                  }
               _P3inc0(dpos);
               k = dpos;
            } while (!(k > buflen));
         }
         else {
            _P3inc0(*loccount);
            _P3strcpy(locnames[*loccount - 1], 255, _P3str1("\010/etc/xdg"));
         }
      }
      else if (P3UTILS_p3appconfig == loctype) {
         buflen = P3UTILS_p3getenvpc(_P3str1("\017XDG_CONFIG_DIRS"), ValueCast(SYSTEM_P3_pansichar, &buf[0]), bufsiz);
         if (buflen >= bufsiz) {
            _P3inc0(*ecount);
            return result;
         }
         else if (buflen > 0) {
            _P3strclr(msg);
            if (SYSTEM_length(appname) > 0)
               _P3strcat(msg, 255, _P3str1("\001/"), appname);
            k = 0;
            dpos = 0;
            do {
               while (buf[dpos] != _P3char('\000') && buf[dpos] != _P3char(':')) {

                  _P3inc0(dpos);
               }
               n = dpos - k;
               if (n + SYSTEM_length(msg) >= sizeof(SYSTEM_shortstring)) {
                  _P3inc0(*ecount);
               }
               else if (n > 0)
                  if (*loccount >= 8) {
                     _P3inc0(*ecount);
                  }
                  else {
                     _P3inc0(*loccount);
                     SYSTEM_move(&buf[k], &locnames[*loccount - 1][1], n);
                     _P3setlength(locnames[*loccount - 1], n, 255);
                     _P3strcat(locnames[*loccount - 1], 255, locnames[*loccount - 1], msg);
                  }
               _P3inc0(dpos);
               k = dpos;
            } while (!(k > buflen));
         }
         else {
            _P3strcpy(msg, 255, _P3str1("\010/etc/xdg"));
            if (SYSTEM_length(appname) > 0) {
               _P3STR_255 _t1;

               _P3strcat(msg, 255, _P3strcat(_t1, 255, msg, _P3str1("\001/")), appname);
            }
            _P3inc0(*loccount);
            _P3strcpy(locnames[*loccount - 1], 255, msg);
         }
      }
      else if (_P3SET_in_1(loctype, P3UTILS_p3data, _P3SET_in_1(loctype, P3UTILS_p3appdata, _P3SET_equal(loctype, P3UTILS_p3applocaldata)))) {
         buflen = P3UTILS_p3getenvpc(_P3str1("\015XDG_DATA_DIRS"), ValueCast(SYSTEM_P3_pansichar, &buf[0]), bufsiz);
         if (buflen >= bufsiz) {
            _P3inc0(*ecount);
            return result;
         }
         else if (buflen > 0) {
            _P3strclr(msg);
            if (SYSTEM_length(appname) > 0)
               _P3strcat(msg, 255, _P3str1("\001/"), appname);
            k = 0;
            dpos = 0;
            do {
               while (buf[dpos] != _P3char('\000') && buf[dpos] != _P3char(':')) {

                  _P3inc0(dpos);
               }
               n = dpos - k;
               if (n > 1 && _P3char('/') == buf[dpos - 1])
                  _P3dec0(n);
               if (n + SYSTEM_length(msg) >= sizeof(SYSTEM_shortstring)) {
                  _P3inc0(*ecount);
               }
               else if (n > 0)
                  if (*loccount >= 8) {
                     _P3inc0(*ecount);
                  }
                  else {
                     _P3inc0(*loccount);
                     SYSTEM_move(&buf[k], &locnames[*loccount - 1][1], n);
                     _P3setlength(locnames[*loccount - 1], n, 255);
                     _P3strcat(locnames[*loccount - 1], 255, locnames[*loccount - 1], msg);
                  }
               _P3inc0(dpos);
               k = dpos;
            } while (!(k > buflen));
         }
         else {
            _P3strcpy(msg, 255, _P3str1("\020/usr/local/share"));
            if (SYSTEM_length(appname) > 0) {
               _P3STR_255 _t1;

               _P3strcat(msg, 255, _P3strcat(_t1, 255, msg, _P3str1("\001/")), appname);
            }
            _P3inc0(*loccount);
            _P3strcpy(locnames[*loccount - 1], 255, msg);
            _P3strcpy(msg, 255, _P3str1("\012/usr/share"));
            if (SYSTEM_length(appname) > 0) {
               _P3STR_255 _t1;

               _P3strcat(msg, 255, _P3strcat(_t1, 255, msg, _P3str1("\001/")), appname);
            }
            _P3inc0(*loccount);
            _P3strcpy(locnames[*loccount - 1], 255, msg);
         }
      }
   }
   return result;
}  /* p3standardlocations */

/* unit p3utils */
void _Init_Module_p3utils(void) {
} /* _Init_Module_p3utils */

void _Final_Module_p3utils(void) {
} /* _Final_Module_p3utils */

