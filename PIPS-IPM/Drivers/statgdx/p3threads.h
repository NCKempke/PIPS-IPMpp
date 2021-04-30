#ifndef _P3___p3threads___H
#define _P3___p3threads___H

cnstdef {P3THREADS_p3_semsize = 64};
cnstdef {P3THREADS_p3_muxsize = 64};
cnstdef {P3THREADS_p3_condvarsize = 64};
/**** C code included from p3threads.pas(56:1): 55 lines ****/
#define FHANDLE(who)             (who)->P3THREADS_tp3thread_DOT_fhandle
#define FTHREADID(who)           (who)->P3THREADS_tp3thread_DOT_fthreadid
#define FCREATESUSPENDEDSEM(who) (who)->P3THREADS_tp3thread_DOT_fcreatesuspendedsem
#define FINITIALSUSPENDDONE(who) (who)->P3THREADS_tp3thread_DOT_finitialsuspenddone
#define FMUX(who)                (who)->P3THREADS_tp3mutex_DOT_fmuxbuf
#define FCV(who)                 (who)->P3THREADS_tp3condvar_DOT_fcvbuf
#define EMSG_HELPER_CREATE(msg,tmp,i)   _P3strcat(msg,255,_P3str1("\027Thread creation error: "),SYSUTILS_P3_syserrormessage(tmp,255,i))
#define _P3THREAD_RAISE_E(msg) _P3_RAISE(ValueCast(P3THREADS_ep3thread,SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,_P3alloc_object(&P3THREADS_ep3thread_CD)),msg)));
#define _P3MUTEX_RAISE_E(msg) _P3_RAISE(ValueCast(P3THREADS_ep3mutex,SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,_P3alloc_object(&P3THREADS_ep3mutex_CD)),msg)));
#define _P3CONDVAR_RAISE_E(msg) _P3_RAISE(ValueCast(P3THREADS_ep3condvar,SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,_P3alloc_object(&P3THREADS_ep3condvar_CD)),msg)));


#if defined(_WIN32)
#  define snprintf _snprintf
#else
#  if (! defined(_POSIX_THREADS))
#    error "_POSIX_THREADS must be defined on this platform"
#  endif
#  if (_POSIX_THREADS < 1)
     /* will this be enough? SOL is only 1, not higher */
#    error "_POSIX_THREADS should be >= 1 on this platform"
#  endif
#  if 0
#  if (_POSIX_THREADS < 200112L)
     /* will this be enough? */
#    error "_POSIX_THREADS should be >= 200112L on this platform"
#  endif
#  if (_POSIX_THREADS < 200809L)
     /* this is a little greedy: maybe we will get this though */
#    error "_POSIX_THREADS should be >= 200809L on this platform"
#  endif
#  endif
#  if defined(__WORDSIZE)
#    if 64 == __WORDSIZE
#      define __P3_SIZEOF_POINTER__ 8
#    else
#      define __P3_SIZEOF_POINTER__ 4
#    endif
#  elif defined(__SIZEOF_POINTER__)
#    define __P3_SIZEOF_POINTER__ __SIZEOF_POINTER__
#  elif defined(__sparcv9)
#    define __P3_SIZEOF_POINTER__ 8
#  elif defined(__sparc)
/*  check __sparc after __sparcv9, both are defined for 64-bit */
#    define __P3_SIZEOF_POINTER__ 4
#  endif
#  if (! defined(__P3_HAVE_SEM_INIT__))
     /* define it to be 1 (default) or 0 (special case for MacOSX) */
#    if defined(__APPLE__)
#      define __P3_HAVE_SEM_INIT__ 0
#    else
#      define __P3_HAVE_SEM_INIT__ 1
#    endif
#  endif
#endif
typedef struct P3THREADS_ep3thread_OD_S* P3THREADS_ep3thread; /* sy_class */
typedef struct P3THREADS_ep3thread_OD_S {  /* Objects of 'ep3thread' */
  SYSTEM_classreference_t CD;  /* = &P3THREADS_ep3thread_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} P3THREADS_ep3thread_OD;

extern void * const P3THREADS_ep3thread_VT[];
extern const SYSTEM_classdescriptor_t P3THREADS_ep3thread_CD;


typedef struct P3THREADS_ep3mutex_OD_S* P3THREADS_ep3mutex; /* sy_class */
typedef struct P3THREADS_ep3mutex_OD_S {  /* Objects of 'ep3mutex' */
  SYSTEM_classreference_t CD;  /* = &P3THREADS_ep3mutex_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} P3THREADS_ep3mutex_OD;

extern void * const P3THREADS_ep3mutex_VT[];
extern const SYSTEM_classdescriptor_t P3THREADS_ep3mutex_CD;


typedef struct P3THREADS_ep3condvar_OD_S* P3THREADS_ep3condvar; /* sy_class */
typedef struct P3THREADS_ep3condvar_OD_S {  /* Objects of 'ep3condvar' */
  SYSTEM_classreference_t CD;  /* = &P3THREADS_ep3condvar_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} P3THREADS_ep3condvar_OD;

extern void * const P3THREADS_ep3condvar_VT[];
extern const SYSTEM_classdescriptor_t P3THREADS_ep3condvar_CD;


typedef struct P3THREADS_ep3sharedmutex_OD_S* P3THREADS_ep3sharedmutex; /* sy_class */
typedef struct P3THREADS_ep3sharedmutex_OD_S {  /* Objects of 'ep3sharedmutex' */
  SYSTEM_classreference_t CD;  /* = &P3THREADS_ep3sharedmutex_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} P3THREADS_ep3sharedmutex_OD;

extern void * const P3THREADS_ep3sharedmutex_VT[];
extern const SYSTEM_classdescriptor_t P3THREADS_ep3sharedmutex_CD;



Prototype Function(SYSTEM_integer ) (*P3THREADS_tp3threadfunc)(
SYSTEM_pointer parameter);

typedef struct P3THREADS_tp3threadrec_S *P3THREADS_pp3threadrec;
typedef struct P3THREADS_tp3threadrec_S {
  P3THREADS_tp3threadfunc func;
  SYSTEM_pointer parameter;
} P3THREADS_tp3threadrec;

typedef struct P3THREADS_tp3handlerec_S {
  union{
    struct{
      SYSTEM_longword h32,dummy;
    } _c1;
    struct{
      SYSTEM_int64 h64;
    } _c2;
    struct{
      SYSTEM_pointer ptr;
    } _c3;
    struct{
      SYSTEM_cardinal dw32;
    } _c4;
    struct{
      SYSTEM_int64 dw64;
    } _c5;
  } _u;
} P3THREADS_tp3handlerec;

typedef SYSTEM_uint8 _sub_1P3THREADS;
typedef SYSTEM_byte _arr_0P3THREADS[64];
typedef struct P3THREADS_tp3thread_OD_S* P3THREADS_tp3thread; /* sy_class */
typedef struct P3THREADS_tp3thread_OD_S {  /* Objects of 'tp3thread' */
  SYSTEM_classreference_t CD;  /* = &P3THREADS_tp3thread_CD */
  P3THREADS_tp3handlerec P3THREADS_tp3thread_DOT_fhandle;
  P3THREADS_tp3handlerec P3THREADS_tp3thread_DOT_fthreadid;
  _arr_0P3THREADS P3THREADS_tp3thread_DOT_fcreatesuspendedsem;
  SYSTEM_boolean P3THREADS_tp3thread_DOT_finitialsuspenddone;
  SYSTEM_boolean P3THREADS_tp3thread_DOT_fterminated;
  SYSTEM_boolean P3THREADS_tp3thread_DOT_ffreeonterminate;
  SYSTEM_boolean P3THREADS_tp3thread_DOT_ffinished;
  SYSTEM_integer P3THREADS_tp3thread_DOT_freturnvalue;
} P3THREADS_tp3thread_OD;


Function(SYSTEM_int64 ) P3THREADS_tp3thread_DOT_getthreadid(
  P3THREADS_tp3thread self);

Function(SYSTEM_boolean ) P3THREADS_tp3thread_DOT_semisnull(
  P3THREADS_tp3thread self);

Function(SYSTEM_boolean ) P3THREADS_tp3thread_DOT_issuspended(
  P3THREADS_tp3thread self);

Prototype Procedure (*P3THREADS_tp3thread_DOT_execute_T)(
  P3THREADS_tp3thread self);

Constructor(P3THREADS_tp3thread ) P3THREADS_tp3thread_DOT_create(
  P3THREADS_tp3thread self,
  SYSTEM_boolean createsuspended);

Destructor(P3THREADS_tp3thread ) P3THREADS_tp3thread_DOT_destroy(
  P3THREADS_tp3thread self);

Procedure P3THREADS_tp3thread_DOT_resume(
  P3THREADS_tp3thread self);

Procedure P3THREADS_tp3thread_DOT_terminate(
  P3THREADS_tp3thread self);

Function(SYSTEM_longword ) P3THREADS_tp3thread_DOT_waitfor(
  P3THREADS_tp3thread self);
extern void * const P3THREADS_tp3thread_VT[];
extern const SYSTEM_classdescriptor_t P3THREADS_tp3thread_CD;


typedef SYSTEM_uint8 _sub_2P3THREADS;
typedef SYSTEM_byte P3THREADS_tmuxbuf[64];
typedef struct P3THREADS_tp3mutex_OD_S* P3THREADS_tp3mutex; /* sy_class */
typedef struct P3THREADS_tp3mutex_OD_S {  /* Objects of 'tp3mutex' */
  SYSTEM_classreference_t CD;  /* = &P3THREADS_tp3mutex_CD */
  P3THREADS_tmuxbuf P3THREADS_tp3mutex_DOT_fmuxbuf;
} P3THREADS_tp3mutex_OD;


Constructor(P3THREADS_tp3mutex ) P3THREADS_tp3mutex_DOT_create(
  P3THREADS_tp3mutex self);

Destructor(P3THREADS_tp3mutex ) P3THREADS_tp3mutex_DOT_destroy(
  P3THREADS_tp3mutex self);

Procedure P3THREADS_tp3mutex_DOT_lock(
  P3THREADS_tp3mutex self);

Function(SYSTEM_boolean ) P3THREADS_tp3mutex_DOT_trylock(
  P3THREADS_tp3mutex self);

Procedure P3THREADS_tp3mutex_DOT_unlock(
  P3THREADS_tp3mutex self);
extern void * const P3THREADS_tp3mutex_VT[];
extern const SYSTEM_classdescriptor_t P3THREADS_tp3mutex_CD;


typedef SYSTEM_int64 P3THREADS_tcvtime;
typedef SYSTEM_uint8 _sub_3P3THREADS;
typedef SYSTEM_byte P3THREADS_tcondvarbuf[64];
typedef struct P3THREADS_tp3condvar_OD_S* P3THREADS_tp3condvar; /* sy_class */
typedef struct P3THREADS_tp3condvar_OD_S {  /* Objects of 'tp3condvar' */
  SYSTEM_classreference_t CD;  /* = &P3THREADS_tp3condvar_CD */
  P3THREADS_tcondvarbuf P3THREADS_tp3condvar_DOT_fcvbuf;
} P3THREADS_tp3condvar_OD;


Constructor(P3THREADS_tp3condvar ) P3THREADS_tp3condvar_DOT_create(
  P3THREADS_tp3condvar self);

Destructor(P3THREADS_tp3condvar ) P3THREADS_tp3condvar_DOT_destroy(
  P3THREADS_tp3condvar self);

Procedure P3THREADS_tp3condvar_DOT_notifyone(
  P3THREADS_tp3condvar self);

Procedure P3THREADS_tp3condvar_DOT_notifyall(
  P3THREADS_tp3condvar self);

Procedure P3THREADS_tp3condvar_DOT_wait(
  P3THREADS_tp3condvar self,
  P3THREADS_tp3mutex mx);

Function(SYSTEM_boolean ) P3THREADS_tp3condvar_DOT_timedwaitabs(
  P3THREADS_tp3condvar self,
  P3THREADS_tp3mutex mx,
  P3THREADS_tcvtime abstime);
extern void * const P3THREADS_tp3condvar_VT[];
extern const SYSTEM_classdescriptor_t P3THREADS_tp3condvar_CD;


typedef struct P3THREADS_tsmstate_S {
  SYSTEM_cardinal nshared;
  SYSTEM_boolean exclusive;
  SYSTEM_boolean upgrade;
  SYSTEM_boolean exclusivewaitingblocked;
  SYSTEM_boolean upgradewaitingblocked;
} P3THREADS_tsmstate;

typedef struct P3THREADS_tp3sharedmutex_OD_S* P3THREADS_tp3sharedmutex; /* sy_class */
typedef struct P3THREADS_tp3sharedmutex_OD_S {  /* Objects of 'tp3sharedmutex' */
  SYSTEM_classreference_t CD;  /* = &P3THREADS_tp3sharedmutex_CD */
  P3THREADS_tsmstate P3THREADS_tp3sharedmutex_DOT_fstate;
  P3THREADS_tp3mutex P3THREADS_tp3sharedmutex_DOT_statechange;
  P3THREADS_tp3condvar P3THREADS_tp3sharedmutex_DOT_shared;
  P3THREADS_tp3condvar P3THREADS_tp3sharedmutex_DOT_upgrade;
  P3THREADS_tp3condvar P3THREADS_tp3sharedmutex_DOT_exclusive;
} P3THREADS_tp3sharedmutex_OD;


Procedure P3THREADS_tp3sharedmutex_DOT_releasewaiters(
  P3THREADS_tp3sharedmutex self);

Constructor(P3THREADS_tp3sharedmutex ) 
  P3THREADS_tp3sharedmutex_DOT_create(
  P3THREADS_tp3sharedmutex self);

Destructor(P3THREADS_tp3sharedmutex ) 
  P3THREADS_tp3sharedmutex_DOT_destroy(
  P3THREADS_tp3sharedmutex self);

Procedure P3THREADS_tp3sharedmutex_DOT_lock(
  P3THREADS_tp3sharedmutex self);

Function(SYSTEM_boolean ) P3THREADS_tp3sharedmutex_DOT_trylock(
  P3THREADS_tp3sharedmutex self);

Procedure P3THREADS_tp3sharedmutex_DOT_unlock(
  P3THREADS_tp3sharedmutex self);

Procedure P3THREADS_tp3sharedmutex_DOT_lockshared(
  P3THREADS_tp3sharedmutex self);

Function(SYSTEM_boolean ) P3THREADS_tp3sharedmutex_DOT_trylockshared(
  P3THREADS_tp3sharedmutex self);

Procedure P3THREADS_tp3sharedmutex_DOT_unlockshared(
  P3THREADS_tp3sharedmutex self);

Procedure P3THREADS_tp3sharedmutex_DOT_lockupgrade(
  P3THREADS_tp3sharedmutex self);

Function(SYSTEM_boolean ) P3THREADS_tp3sharedmutex_DOT_trylockupgrade(
  P3THREADS_tp3sharedmutex self);

Procedure P3THREADS_tp3sharedmutex_DOT_unlockupgrade(
  P3THREADS_tp3sharedmutex self);

Procedure P3THREADS_tp3sharedmutex_DOT_unlockupgrade_and_lock(
  P3THREADS_tp3sharedmutex self);

Procedure P3THREADS_tp3sharedmutex_DOT_unlock_and_lockupgrade(
  P3THREADS_tp3sharedmutex self);

Procedure P3THREADS_tp3sharedmutex_DOT_unlock_and_lockshared(
  P3THREADS_tp3sharedmutex self);

Procedure P3THREADS_tp3sharedmutex_DOT_unlockupgrade_and_lockshared(
  P3THREADS_tp3sharedmutex self);
extern void * const P3THREADS_tp3sharedmutex_VT[];
extern const SYSTEM_classdescriptor_t P3THREADS_tp3sharedmutex_CD;



Procedure P3THREADS_incdt(
  SYSTEM_P3_tdatetime *dt,
  SYSTEM_cardinal ticks);

Procedure P3THREADS_decdt(
  SYSTEM_P3_tdatetime *dt,
  SYSTEM_cardinal ticks);

Procedure P3THREADS_gettimespec(
  SYSTEM_P3_tdatetime dt,
  SYSTEM_int64 *tv_sec,
  SYSTEM_int64 *tv_nsec);

Function(P3THREADS_tcvtime ) P3THREADS_encodedatecv(
  SYSTEM_word year,
  SYSTEM_word month,
  SYSTEM_word day);

Function(P3THREADS_tcvtime ) P3THREADS_encodetimecv(
  SYSTEM_word hour,
  SYSTEM_word _min,
  SYSTEM_word sec,
  SYSTEM_word msec);

Function(SYSTEM_boolean ) P3THREADS_decodedatefullycv(
  P3THREADS_tcvtime cvt,
  SYSTEM_word *year,
  SYSTEM_word *month,
  SYSTEM_word *day,
  SYSTEM_word *dow);

Procedure P3THREADS_decodedatecv(
  P3THREADS_tcvtime cvt,
  SYSTEM_word *year,
  SYSTEM_word *month,
  SYSTEM_word *day);

Procedure P3THREADS_decodetimecv(
  P3THREADS_tcvtime cvt,
  SYSTEM_word *hour,
  SYSTEM_word *_min,
  SYSTEM_word *sec,
  SYSTEM_word *msec);

Function(P3THREADS_tcvtime ) P3THREADS_nowcv(void);

Function(P3THREADS_tcvtime ) P3THREADS_datecv(void);

Function(P3THREADS_tcvtime ) P3THREADS_timecv(void);

Procedure P3THREADS_inccvtimemillis(
  P3THREADS_tcvtime *cvt,
  SYSTEM_cardinal ticks);

Procedure P3THREADS_deccvtimemillis(
  P3THREADS_tcvtime *cvt,
  SYSTEM_cardinal ticks);

Function(SYSTEM_cardinal ) P3THREADS_tickdeltacv(
  P3THREADS_tcvtime cvt);

extern void _Init_Module_p3threads(void);
extern void _Final_Module_p3threads(void);

#endif /* ! defined _P3___p3threads___H */
