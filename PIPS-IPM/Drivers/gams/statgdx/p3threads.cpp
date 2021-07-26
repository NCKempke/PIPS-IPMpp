#include "p3io.h"
#include "system_p3.h"
#include "exceptions.h"
#include "sysutils_p3.h"
#include "p3threads.h"


void * const P3THREADS_ep3thread_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'ep3thread' */
const SYSTEM_classdescriptor_t P3THREADS_ep3thread_CD = {
  _P3str1("\011ep3thread"), 
  &SYSTEM_exception_CD, NULL, 0, 
  sizeof(P3THREADS_ep3thread_OD), P3THREADS_ep3thread_VT, NULL};


void * const P3THREADS_ep3mutex_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'ep3mutex' */
const SYSTEM_classdescriptor_t P3THREADS_ep3mutex_CD = {
  _P3str1("\010ep3mutex"), 
  &SYSTEM_exception_CD, NULL, 0, 
  sizeof(P3THREADS_ep3mutex_OD), P3THREADS_ep3mutex_VT, NULL};


void * const P3THREADS_ep3condvar_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'ep3condvar' */
const SYSTEM_classdescriptor_t P3THREADS_ep3condvar_CD = {
  _P3str1("\012ep3condvar"), 
  &SYSTEM_exception_CD, NULL, 0, 
  sizeof(P3THREADS_ep3condvar_OD), P3THREADS_ep3condvar_VT, NULL};


void * const P3THREADS_ep3sharedmutex_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'ep3sharedmutex' */
const SYSTEM_classdescriptor_t P3THREADS_ep3sharedmutex_CD = {
  _P3str1("\016ep3sharedmutex"), 
  &SYSTEM_exception_CD, NULL, 0, 
  sizeof(P3THREADS_ep3sharedmutex_OD), P3THREADS_ep3sharedmutex_VT, NULL};


void * const P3THREADS_tp3thread_VT[] = {(void*)&
  P3THREADS_tp3thread_DOT_destroy, (void*)&_P3_abstract_call1};

/* Class descriptor for 'tp3thread' */
const SYSTEM_classdescriptor_t P3THREADS_tp3thread_CD = {
  _P3str1("\011tp3thread"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(P3THREADS_tp3thread_OD), P3THREADS_tp3thread_VT, NULL};


void * const P3THREADS_tp3mutex_VT[] = {(void*)&
  P3THREADS_tp3mutex_DOT_destroy};

/* Class descriptor for 'tp3mutex' */
const SYSTEM_classdescriptor_t P3THREADS_tp3mutex_CD = {
  _P3str1("\010tp3mutex"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(P3THREADS_tp3mutex_OD), P3THREADS_tp3mutex_VT, NULL};


void * const P3THREADS_tp3condvar_VT[] = {(void*)&
  P3THREADS_tp3condvar_DOT_destroy};

/* Class descriptor for 'tp3condvar' */
const SYSTEM_classdescriptor_t P3THREADS_tp3condvar_CD = {
  _P3str1("\012tp3condvar"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(P3THREADS_tp3condvar_OD), P3THREADS_tp3condvar_VT, NULL};


void * const P3THREADS_tp3sharedmutex_VT[] = {(void*)&
  P3THREADS_tp3sharedmutex_DOT_destroy};

/* Class descriptor for 'tp3sharedmutex' */
const SYSTEM_classdescriptor_t P3THREADS_tp3sharedmutex_CD = {
  _P3str1("\016tp3sharedmutex"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(P3THREADS_tp3sharedmutex_OD), P3THREADS_tp3sharedmutex_VT, NULL};

cnstdef {P3THREADS_date_delta_cv_dt = 109205};
cnstdef {P3THREADS_date_delta_cv_unix = 134774};
cnstdef {P3THREADS_secs_per_day = 86400};
cnstdef {P3THREADS_secs_per_year = 31536000};
cnstdef {P3THREADS_nicks_per_usec = 10};
cnstdef {P3THREADS_nicks_per_msec = 10000};
cnstdef {P3THREADS_nicks_per_sec = 10000000};
/**** C code included from p3threads.pas(336:1): 9 lines ****/
#if defined(_WIN32)
typedef union cvTime {
  __int64 i64;
  FILETIME ft;
} cvTime_t;
#else
# include <sys/time.h>
# include <time.h>
#endif
static SYSTEM_int64 P3THREADS_nickspermin;
static SYSTEM_int64 P3THREADS_nicksperhour;
static SYSTEM_int64 P3THREADS_nicksperday;
/**** C code included from p3threads.pas(357:1): 77 lines ****/
#define nicksperday          P3THREADS_nicksperday
#define date_delta_cv_unix   P3THREADS_date_delta_cv_unix
#define nicks_per_sec        P3THREADS_nicks_per_sec
#define nicks_per_usec       P3THREADS_nicks_per_usec
static SYSTEM_char * mkExceptMsg (SYSTEM_integer errnum, const char *txt, SYSTEM_shortstring msg)
{
  SYSTEM_shortstring tmp;
  int n;

  (void) SYSUTILS_P3_syserrormessage (tmp,255,errnum);
  n = snprintf ((char *)msg+1, 255, "%s: %.*s (%d)", txt, (int)(*tmp), tmp+1, errnum);
  /* stupid MS _snprintf */
  if ((n > 255) || (n < 0))
    n = 255;
  *msg = (SYSTEM_byte) n;
  return msg;
} /* mkExceptMsg */

#if ! defined(_WIN32)

typedef union pthreadPtrRec {
  void *vv;
  pthread_t pp;
} pthreadPtr_t;

static pthread_t
voidPtr2threadHandle (void *v)
{
  pthreadPtr_t r;

  r.vv = v;
  return r.pp;
} /* voidPtr2threadHandle */

static void *
threadHandle2voidPtr (pthread_t p)
{
  pthreadPtr_t r;

  r.vv = NULL;
  r.pp = p;
  return r.vv;
} /* threadHandle2voidPtr */

#define DELTA_DATETIME2TIMESPEC 25569; /* days(now) - DELTA = days(gettimeofday) */

static void datetime2timespec (double dt, struct timespec *ts)
{
  SYSTEM_int64 i64;
  SYSTEM_double r;

  ts->tv_sec = 0;
  ts->tv_nsec = 0;
  i64 = (SYSTEM_int64) dt;
  r = dt - i64;
  i64 -= DELTA_DATETIME2TIMESPEC;
  if (i64 < 0)
    return;
  i64 *= (60 * 60 * 24);
  r *= (60 * 60 * 24);
  ts->tv_sec = (time_t) r;
  r = r - ts->tv_sec;
  ts->tv_sec += i64;
  ts->tv_nsec = (long) (r * 1e9);
} /* datetime2timespec */

static void cvtime2timespec (SYSTEM_int64 ctv, struct timespec *ts)
{
  SYSTEM_int64 t;

  t = ctv - nicksperday * date_delta_cv_unix;
  ts->tv_sec = t / nicks_per_sec;
  t = t % nicks_per_sec;
  ts->tv_nsec = t * 100;
} /* cvtime2timespec */

#endif

Procedure P3THREADS_gettimespec(
  SYSTEM_P3_tdatetime dt,
  SYSTEM_int64 *tv_sec,
  SYSTEM_int64 *tv_nsec)
{
  *tv_sec =  -1;
  *tv_nsec =  -1;
  /**** C code included from p3threads.pas(443:1): 9 lines ****/
#if ! defined(_WIN32)
{
  struct timespec ts;

  datetime2timespec (dt, &ts);
  *tv_sec  = ts.tv_sec;
  *tv_nsec = ts.tv_nsec;
}
#endif
}  /* gettimespec */

static Function(SYSTEM_cardinal )  STDCALL P3THREADS_threadwrapper(
  SYSTEM_pointer arg)
{
  SYSTEM_cardinal result;
  P3THREADS_pp3threadrec p;
  P3THREADS_tp3threadfunc f;
  SYSTEM_pointer ptr;

  p = ValueCast(P3THREADS_pp3threadrec,arg);
  f = p->func;
  ptr = p->parameter;
  _P3freemem(p);
  result = (*f)(ptr);
  return result;
}  /* threadwrapper */
/**** C code included from p3threads.pas(480:1): 34 lines ****/
#if defined(_WIN32)
static DWORD STDCALL
threadWrapper(SYSTEM_pointer arg)
{
  SYSTEM_integer result;
  P3THREADS_pp3threadrec p;
  P3THREADS_tp3threadfunc f;
  SYSTEM_pointer ptr;

  p = ValueCast(P3THREADS_pp3threadrec,arg);
  f = p->func;
  ptr = p->parameter;
  _P3freemem(p);
  result = (*f)(ptr);
  return result;
} /* threadwrapper */
#else

C_LINKAGE(static void * threadWrapper(SYSTEM_pointer arg);)
static void *
threadWrapper(SYSTEM_pointer arg)
{
  P3THREADS_pp3threadrec p;
  P3THREADS_tp3threadfunc f;
  SYSTEM_pointer ptr;

  p = ValueCast(P3THREADS_pp3threadrec,arg);
  f = p->func;
  ptr = p->parameter;
  _P3freemem(p);
  (void) (*f)(ptr);
  return NULL;
} /* threadWrapper */
#endif

static Function(SYSTEM_integer ) P3THREADS_threadproc(
  P3THREADS_tp3thread thread)
{
  SYSTEM_integer result;
  SYSTEM_boolean dofreeonterm;

  /**** C code included from p3threads.pas(525:1): 18 lines ****/
#if ! defined(_WIN32)
  /* do this unconditionally: we always start suspended */
# if 1 == __P3_HAVE_SEM_INIT__
{
  sem_t *createSuspendedSem;

  createSuspendedSem = (sem_t *) &(FCREATESUSPENDEDSEM(thread));
  (void) sem_wait (createSuspendedSem);
}
# else
{
  pthread_mutex_t *createSuspendedMux;

  createSuspendedMux = (pthread_mutex_t *) &(FCREATESUSPENDEDSEM(thread));
  (void) pthread_mutex_lock (createSuspendedMux);  // use mutex like a binary semaphore
}
# endif
#endif
  if (!thread->P3THREADS_tp3thread_DOT_fterminated) 
    VirtMethodCall(thread, P3THREADS_tp3thread_DOT_execute_T, 1, (
      thread));
  dofreeonterm = thread->P3THREADS_tp3thread_DOT_ffreeonterminate;
  result = thread->P3THREADS_tp3thread_DOT_freturnvalue;
  thread->P3THREADS_tp3thread_DOT_ffinished = SYSTEM_true;
  if (dofreeonterm) 
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,thread));
  /**** C code included from p3threads.pas(556:1): 5 lines ****/
#if defined(_WIN32)
  ExitThread(result);
#else
  pthread_exit((void *)(size_t)result);
#endif
  return result;
}  /* threadproc */

Function(SYSTEM_int64 ) P3THREADS_tp3thread_DOT_getthreadid(
  P3THREADS_tp3thread self)
{
  SYSTEM_int64 result;

  /**** C code included from p3threads.pas(572:1): 20 lines ****/
#if defined(_WIN32)
#  if defined(_WIN64)
     result = FTHREADID(self)._u._c5.dw64;
#  else
     result = FTHREADID(self)._u._c4.dw32;
#  endif
#else
#  if defined(__P3_SIZEOF_POINTER__)
#    if 8 == __P3_SIZEOF_POINTER__
       result = FTHREADID(self)._u._c5.dw64;
#    else
       result = FTHREADID(self)._u._c4.dw32;
#    endif
#  else
     if (8 == sizeof(void *))
       result = FTHREADID(self)._u._c5.dw64;
     else
       result = FTHREADID(self)._u._c4.dw32;
#  endif
#endif
  return result;
}  /* getthreadid */

Function(SYSTEM_boolean ) P3THREADS_tp3thread_DOT_semisnull(
  P3THREADS_tp3thread self)
{
  SYSTEM_boolean result;
  SYSTEM_integer i;

  result = SYSTEM_false;
  for (i = 1;i <= 63;++i) {
    if (0 != self->P3THREADS_tp3thread_DOT_fcreatesuspendedsem[i]) 
      return result;
  }
  result = SYSTEM_true;
  return result;
}  /* semisnull */

Function(SYSTEM_boolean ) P3THREADS_tp3thread_DOT_issuspended(
  P3THREADS_tp3thread self)
{
  SYSTEM_boolean result;

  result = !self->P3THREADS_tp3thread_DOT_finitialsuspenddone;
  return result;
}  /* issuspended */

Constructor(P3THREADS_tp3thread ) P3THREADS_tp3thread_DOT_create(
  P3THREADS_tp3thread self,
  SYSTEM_boolean createsuspended)
{
  P3THREADS_pp3threadrec p;
  SYSTEM_shortstring msg, tmp;
  SYSTEM_integer nbytes;

  ValueCast(P3THREADS_tp3thread,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  if (!createsuspended) 
    _P3_RAISE(ValueCast(P3THREADS_ep3thread,
      SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,
      _P3alloc_object(&P3THREADS_ep3thread_CD)),_P3str1("\063Thread creation error: createSuspended must be true"))));
  nbytes = 0;
  /**** C code included from p3threads.pas(626:1): 7 lines ****/
#if ! defined(_WIN32)
# if 1 == __P3_HAVE_SEM_INIT__
  nbytes = sizeof(sem_t);
# else
  nbytes = sizeof(pthread_mutex_t);
# endif
#endif
  if (nbytes > 64) {
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;
      _P3STR_255 _t3;
      SYSTEM_shortstring _t4;
      _P3STR_255 _t5;

      _P3strcat(msg,255,_P3strcat(_t5,255,_P3strcat(_t3,255,
        _P3strcat(_t2,255,_P3str1("\064internal Tp3Thread error: sizeof(semaphore/mutex) = "),
        SYSUTILS_P3_inttostr(_t1,255,nbytes)),_P3str1("\023 but only reserved ")),
        SYSUTILS_P3_inttostr(_t4,255,P3THREADS_p3_semsize)),_P3str1("\006 bytes"));
    }
    _P3_RAISE(ValueCast(P3THREADS_ep3thread,
      SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,
      _P3alloc_object(&P3THREADS_ep3thread_CD)),msg)));
  } 
  nbytes = 0;
  /**** C code included from p3threads.pas(646:1): 3 lines ****/
#if ! defined(_WIN32)
  nbytes = sizeof(pthread_t);
#endif
  if (nbytes > ValueCast(SYSTEM_int32,sizeof(SYSTEM_pointer))) {
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;
      _P3STR_255 _t3;
      SYSTEM_shortstring _t4;

      _P3strcat(msg,255,_P3strcat(_t3,255,_P3strcat(_t2,255,_P3str1("\077internal Tp3Thread error: sizeof(pthread_t) > sizeof(pointer): "),
        SYSUTILS_P3_inttostr(_t1,255,nbytes)),_P3str1("\003 > ")),
        SYSUTILS_P3_inttostr(_t4,255,sizeof(SYSTEM_pointer)));
    }
    _P3_RAISE(ValueCast(P3THREADS_ep3thread,
      SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,
      _P3alloc_object(&P3THREADS_ep3thread_CD)),msg)));
  } 
  _P3new(p);
  p->func = ValueCast(P3THREADS_tp3threadfunc,&P3THREADS_threadproc);
  p->parameter = ValueCast(SYSTEM_pointer,self);
  /**** C code included from p3threads.pas(689:1): 73 lines ****/
#if defined(_WIN32)
{
  HANDLE h;
  DWORD tid;
  int rc;

  h = CreateThread (NULL, 0, &threadWrapper, p, CREATE_SUSPENDED, &tid);
  if (0 == h) {
    rc = GetLastError();
    EMSG_HELPER_CREATE(msg,tmp,rc);
    _P3THREAD_RAISE_E(msg);
  }
#  if defined(_WIN64)
  FHANDLE(self)._u._c2.h64 = (__int64) h;
  FTHREADID(self)._u._c5.dw64 = (__int64) tid;
#  else
  FHANDLE(self)._u._c1.h32 = (SYSTEM_longword) h;
  FTHREADID(self)._u._c4.dw32 = tid;

#  endif
}
#else
{
  int rc;
  pthread_t pth;

# if 1 == __P3_HAVE_SEM_INIT__
  sem_t *createSuspendedSem = (sem_t *) &(FCREATESUSPENDEDSEM(self));
  rc = sem_init (createSuspendedSem, 0, 0);
  if (rc) {
    rc = errno;
    EMSG_HELPER_CREATE(msg,tmp,rc);
    _P3THREAD_RAISE_E(msg);
  }
# else
  pthread_mutex_t *createSuspendedMux = (pthread_mutex_t *) &(FCREATESUSPENDEDSEM(self));
  rc = pthread_mutex_init (createSuspendedMux, NULL);
  if (rc) {
    rc = errno;
    EMSG_HELPER_CREATE(msg,tmp,rc);
    _P3THREAD_RAISE_E(msg);
  }
  rc = pthread_mutex_lock (createSuspendedMux);
  if (rc) {
    rc = errno;
    EMSG_HELPER_CREATE(msg,tmp,rc);
    _P3THREAD_RAISE_E(msg);
  }
#endif

  rc = pthread_create (&pth, NULL, &threadWrapper, p);
  if (rc) {
# if 1 == __P3_HAVE_SEM_INIT__
    (void) sem_destroy (createSuspendedSem);
# else
    (void) pthread_mutex_unlock (createSuspendedMux);
    (void) pthread_mutex_destroy (createSuspendedMux);
# endif
    EMSG_HELPER_CREATE(msg,tmp,rc);
    _P3THREAD_RAISE_E(msg);
  }
#if 0
  FTHREADID(self)._u._c3.ptr = (void *) pth;
#else
  FTHREADID(self)._u._c3.ptr = threadHandle2voidPtr (pth);
#endif
  /*
  ErrCode := BeginThread(nil, @ThreadProc, Pointer(Self), FThreadID);
  if ErrCode <> 0 then
    raise EThread.CreateResFmt(@SThreadCreateError, [SysErrorMessage(ErrCode)]);
   */
}
#endif
  return self;
}  /* create */

Destructor(P3THREADS_tp3thread ) P3THREADS_tp3thread_DOT_destroy(
  P3THREADS_tp3thread self)
{
  if (self->P3THREADS_tp3thread_DOT_fthreadid._u._c2.h64 != 0 && !
    self->P3THREADS_tp3thread_DOT_ffinished) {
    P3THREADS_tp3thread_DOT_terminate(self);
    P3THREADS_tp3thread_DOT_resume(self);
    P3THREADS_tp3thread_DOT_waitfor(self);
  } 
  /**** C code included from p3threads.pas(791:1): 31 lines ****/
#if defined(_WIN32)
#  if defined(_WIN64)
     if (FHANDLE(self)._u._c2.h64)
       CloseHandle ((HANDLE) FHANDLE(self)._u._c2.h64);
#  else
     if (FHANDLE(self)._u._c1.h32)
       CloseHandle ((HANDLE) FHANDLE(self)._u._c1.h32);
#  endif
#else
  {
    pthread_t pth;

#if 0
    pth = (pthread_t) FTHREADID(self)._u._c3.ptr;
#else
    pth = voidPtr2threadHandle (FTHREADID(self)._u._c3.ptr);
#endif
    if (pth) {
      pthread_detach(pth);
    }
    if (! P3THREADS_tp3thread_DOT_semisnull(self)) {
# if 1 == __P3_HAVE_SEM_INIT__
      sem_t *createSuspendedSem = (sem_t *) &(FCREATESUSPENDEDSEM(self));
      (void) sem_destroy (createSuspendedSem);
#else
      pthread_mutex_t *createSuspendedMux = (pthread_mutex_t *) &(FCREATESUSPENDEDSEM(self));
      (void) pthread_mutex_destroy (createSuspendedMux);
#endif
    }
  }
#endif
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure P3THREADS_tp3thread_DOT_resume(
  P3THREADS_tp3thread self)
{
  /**** C code included from p3threads.pas(862:1): 39 lines ****/
#if defined(_WIN32)
{
  int suspendCount;
  HANDLE h;

  if (! FINITIALSUSPENDDONE(self)) {
#  if defined(_WIN64)
    h = (HANDLE) FHANDLE(self)._u._c2.h64;
#  else
    h = (HANDLE) FHANDLE(self)._u._c1.h32;
#  endif
    FINITIALSUSPENDDONE(self) = SYSTEM_true;
    suspendCount = ResumeThread (h);
    /* assert(1 == suspendCount); */
  }
}
#else
# if 1 == __P3_HAVE_SEM_INIT__
{
  sem_t *createSuspendedSem;

  if (! FINITIALSUSPENDDONE(self)) {
    FINITIALSUSPENDDONE(self) = SYSTEM_true;
    createSuspendedSem = (sem_t *) &(FCREATESUSPENDEDSEM(self));
    sem_post(createSuspendedSem);
  }
}
# else
{
  pthread_mutex_t *createSuspendedMux;

  if (! FINITIALSUSPENDDONE(self)) {
    FINITIALSUSPENDDONE(self) = SYSTEM_true;
    createSuspendedMux = (pthread_mutex_t *) &(FCREATESUSPENDEDSEM(self));
    (void) pthread_mutex_unlock (createSuspendedMux);
  }
}
# endif
#endif
}  /* resume */

Procedure P3THREADS_tp3thread_DOT_terminate(
  P3THREADS_tp3thread self)
{
  self->P3THREADS_tp3thread_DOT_fterminated = SYSTEM_true;
}  /* terminate */

Function(SYSTEM_longword ) P3THREADS_tp3thread_DOT_waitfor(
  P3THREADS_tp3thread self)
{
  SYSTEM_longword result;
  SYSTEM_integer errnum;
  SYSTEM_shortstring msg;

  /**** C code included from p3threads.pas(939:1): 45 lines ****/
#if defined(_WIN32)
{
  int waitResult;
  HANDLE h;
  BOOL brc;

#if defined(_WIN64)
  h = (HANDLE) FHANDLE(self)._u._c2.h64;
#else
  h = (HANDLE) FHANDLE(self)._u._c1.h32;
#endif
  waitResult = WaitForSingleObject (h, INFINITE);
  brc = GetExitCodeThread(h, (LPDWORD)&result);
  if (! brc) {
    errnum = GetLastError();
    _P3THREAD_RAISE_E(mkExceptMsg(errnum,"Tp3Thread.waitFor",msg));
  }
}
#else
{
  pthread_t pth;
  void *vp = NULL;

#if 0
  pth = (pthread_t) FTHREADID(self)._u._c3.ptr;
#else
  pth = voidPtr2threadHandle (FTHREADID(self)._u._c3.ptr);
#endif
  FTHREADID(self)._u._c2.h64 = 0;
  errnum = pthread_join (pth, &vp);

  if (errnum) {
    _P3THREAD_RAISE_E(mkExceptMsg(errnum,"Tp3Thread.waitFor",msg));
  }
#if 0  /* this is rejected by C++ */
  result = (SYSTEM_longword) vp;
#else
  {
  SYSTEM_int64 r64;
  r64 = (SYSTEM_int64) vp;
  result = (SYSTEM_longword) r64;
  }
#endif
}
#endif
  return result;
}  /* waitfor */

Constructor(P3THREADS_tp3mutex ) P3THREADS_tp3mutex_DOT_create(
  P3THREADS_tp3mutex self)
{
  SYSTEM_shortstring msg;
  SYSTEM_integer nbytes;

  ValueCast(P3THREADS_tp3mutex,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  nbytes = 0;
  /**** C code included from p3threads.pas(999:1): 5 lines ****/
#if defined(_WIN32)
  nbytes = sizeof(CRITICAL_SECTION);
#else
  nbytes = sizeof(pthread_mutex_t);
#endif
  if (nbytes > 64) {
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;
      _P3STR_255 _t3;
      SYSTEM_shortstring _t4;
      _P3STR_255 _t5;

      _P3strcat(msg,255,_P3strcat(_t5,255,_P3strcat(_t3,255,
        _P3strcat(_t2,255,_P3str1("\065internal Tp3Mutex error: sizeof(os_specific_mutex) = "),
        SYSUTILS_P3_inttostr(_t1,255,nbytes)),_P3str1("\023 but only reserved ")),
        SYSUTILS_P3_inttostr(_t4,255,P3THREADS_p3_muxsize)),_P3str1("\006 bytes"));
    }
    _P3_RAISE(ValueCast(P3THREADS_ep3mutex,SYSTEM_exception_DOT_create(ValueCast(
      SYSTEM_exception,_P3alloc_object(&P3THREADS_ep3mutex_CD)),msg)));
  } 
  /**** C code included from p3threads.pas(1018:1): 15 lines ****/
#if defined(_WIN32)
{
  CRITICAL_SECTION *pCS = (CRITICAL_SECTION *) &(FMUX(self));
  InitializeCriticalSection (pCS);
}
#else
{
  pthread_mutex_t *pMux = (pthread_mutex_t *) &(FMUX(self));
  int rc;
  rc = pthread_mutex_init (pMux, NULL);
  if (rc) {
    _P3MUTEX_RAISE_E(mkExceptMsg(rc,"Tp3Mutex.Create",msg));
  }
}
#endif
  return self;
}  /* create */

Destructor(P3THREADS_tp3mutex ) P3THREADS_tp3mutex_DOT_destroy(
  P3THREADS_tp3mutex self)
{
  /**** C code included from p3threads.pas(1043:1): 22 lines ****/
#if defined(_WIN32)
{
  CRITICAL_SECTION *pCS = (CRITICAL_SECTION *) &(FMUX(self));
  DeleteCriticalSection (pCS);
}
#else
{
  pthread_mutex_t *pMux = (pthread_mutex_t *) &(FMUX(self));

# if 1
  (void) pthread_mutex_destroy (pMux);
# else
  {
    int rc = pthread_mutex_destroy (pMux);
    if (rc) {
      SYSTEM_shortstring msg;
      _P3MUTEX_RAISE_E(mkExceptMsg(rc,"Tp3Mutex.Destroy",msg));
    }
  }
# endif
}
#endif
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure P3THREADS_tp3mutex_DOT_lock(
  P3THREADS_tp3mutex self)
{
  /**** C code included from p3threads.pas(1076:1): 17 lines ****/
#if defined(_WIN32)
{
  CRITICAL_SECTION *pCS = (CRITICAL_SECTION *) &(FMUX(self));
  EnterCriticalSection (pCS);
}
#else
{
  pthread_mutex_t *pMux = (pthread_mutex_t *) &(FMUX(self));
  int rc;

  rc = pthread_mutex_lock (pMux);
  if (rc) {
    SYSTEM_shortstring msg;
    _P3MUTEX_RAISE_E(mkExceptMsg(rc,"Tp3Mutex.lock",msg));
  }
}
#endif
}  /* lock */

Function(SYSTEM_boolean ) P3THREADS_tp3mutex_DOT_trylock(
  P3THREADS_tp3mutex self)
{
  SYSTEM_boolean result;

  /**** C code included from p3threads.pas(1109:1): 22 lines ****/
#if defined(_WIN32)
{
  CRITICAL_SECTION *pCS = (CRITICAL_SECTION *) &(FMUX(self));
  BOOL gotIt;
  gotIt = TryEnterCriticalSection (pCS);
  result = gotIt;
}
#else
{
  pthread_mutex_t *pMux = (pthread_mutex_t *) &(FMUX(self));
  int rc;

  result = SYSTEM_true;
  rc = pthread_mutex_trylock (pMux);
  if (EBUSY == rc)
    result = SYSTEM_false;
  else if (rc) {
    SYSTEM_shortstring msg;
    _P3MUTEX_RAISE_E(mkExceptMsg(rc,"Tp3Mutex.trylock",msg));
  }
}
#endif
  return result;
}  /* trylock */

Procedure P3THREADS_tp3mutex_DOT_unlock(
  P3THREADS_tp3mutex self)
{
  SYSTEM_shortstring msg;

  /**** C code included from p3threads.pas(1145:1): 16 lines ****/
#if defined(_WIN32)
{
  CRITICAL_SECTION *pCS = (CRITICAL_SECTION *) &(FMUX(self));
  LeaveCriticalSection (pCS);
}
#else
{
  pthread_mutex_t *pMux = (pthread_mutex_t *) &(FMUX(self));
  int rc;

  rc = pthread_mutex_unlock (pMux);
  if (rc) {
    _P3MUTEX_RAISE_E(mkExceptMsg(rc,"Tp3Mutex.unlock",msg));
  }
}
#endif
}  /* unlock */

Constructor(P3THREADS_tp3condvar ) P3THREADS_tp3condvar_DOT_create(
  P3THREADS_tp3condvar self)
{
  SYSTEM_shortstring msg;
  SYSTEM_integer nbytes;

  ValueCast(P3THREADS_tp3condvar,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  nbytes = 0;
  /**** C code included from p3threads.pas(1218:1): 5 lines ****/
#if defined(_WIN32)
  nbytes = sizeof(CONDITION_VARIABLE);
#else
  nbytes = sizeof(pthread_cond_t);
#endif
  if (nbytes > 64) {
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;
      _P3STR_255 _t3;
      SYSTEM_shortstring _t4;
      _P3STR_255 _t5;

      _P3strcat(msg,255,_P3strcat(_t5,255,_P3strcat(_t3,255,
        _P3strcat(_t2,255,_P3str1("\071internal Tp3CondVar error: sizeof(os_specific_condVar) = "),
        SYSUTILS_P3_inttostr(_t1,255,nbytes)),_P3str1("\023 but only reserved ")),
        SYSUTILS_P3_inttostr(_t4,255,P3THREADS_p3_condvarsize)),_P3str1("\006 bytes"));
    }
    _P3_RAISE(ValueCast(P3THREADS_ep3condvar,
      SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,
      _P3alloc_object(&P3THREADS_ep3condvar_CD)),msg)));
  } 
  /**** C code included from p3threads.pas(1233:1): 15 lines ****/
#if defined(_WIN32)
{
  CONDITION_VARIABLE *pCV = (CONDITION_VARIABLE *) &(FCV(self));
  InitializeConditionVariable (pCV);
}
#else
{
  pthread_cond_t *pCV = (pthread_cond_t *) &(FCV(self));
  int rc;
  rc = pthread_cond_init (pCV, NULL);
  if (rc) {
    _P3CONDVAR_RAISE_E(mkExceptMsg(rc,"Tp3CondVar.Create",msg));
  }
}
#endif
  return self;
}  /* create */

Destructor(P3THREADS_tp3condvar ) P3THREADS_tp3condvar_DOT_destroy(
  P3THREADS_tp3condvar self)
{
  /**** C code included from p3threads.pas(1261:1): 16 lines ****/
#if defined(_WIN32)
  /* nothing to do: no DeleteConditionVariable or similar */
{
}
#else
{
  pthread_cond_t *pCV = (pthread_cond_t *) &(FCV(self));
# if 0
  int rc;
  rc = pthread_cond_destroy (pCV);
  assert(0==rc); /* but why handle the error on destroy? */
# else
  (void) pthread_cond_destroy (pCV);
# endif
}
#endif
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure P3THREADS_tp3condvar_DOT_notifyone(
  P3THREADS_tp3condvar self)
{
  SYSTEM_shortstring msg;

  /**** C code included from p3threads.pas(1293:1): 16 lines ****/

#if defined(_WIN32)
{
  CONDITION_VARIABLE *pCV = (CONDITION_VARIABLE *) &(FCV(self));
  WakeConditionVariable (pCV);
}
#else
{
  pthread_cond_t *pCV = (pthread_cond_t *) &(FCV(self));
  int rc;
  rc = pthread_cond_signal (pCV);
  if (rc) {
    _P3CONDVAR_RAISE_E(mkExceptMsg(rc,"Tp3CondVar.notifyOne",msg));
  }
}
#endif
}  /* notifyone */

Procedure P3THREADS_tp3condvar_DOT_notifyall(
  P3THREADS_tp3condvar self)
{
  SYSTEM_shortstring msg;

  /**** C code included from p3threads.pas(1324:1): 15 lines ****/
#if defined(_WIN32)
{
  CONDITION_VARIABLE *pCV = (CONDITION_VARIABLE *) &(FCV(self));
  WakeAllConditionVariable (pCV);
}
#else
{
  pthread_cond_t *pCV = (pthread_cond_t *) &(FCV(self));
  int rc;
  rc = pthread_cond_broadcast (pCV);
  if (rc) {
    _P3CONDVAR_RAISE_E(mkExceptMsg(rc,"Tp3CondVar.notifyAll",msg));
  }
}
#endif
}  /* notifyall */

Procedure P3THREADS_tp3condvar_DOT_wait(
  P3THREADS_tp3condvar self,
  P3THREADS_tp3mutex mx)
{
  SYSTEM_shortstring msg;

  /**** C code included from p3threads.pas(1355:1): 17 lines ****/
#if defined(_WIN32)
{
  CONDITION_VARIABLE *pCV = (CONDITION_VARIABLE *) &(FCV(self));
  CRITICAL_SECTION *pCS = (CRITICAL_SECTION *) &(FMUX(mx));
  SleepConditionVariableCS (pCV, pCS, INFINITE);
}
#else
{
  pthread_cond_t *pCV = (pthread_cond_t *) &(FCV(self));
  pthread_mutex_t *pMux = (pthread_mutex_t *) &(FMUX(mx));
  int rc;
  rc = pthread_cond_wait (pCV, pMux);
  if (rc) {
    _P3CONDVAR_RAISE_E(mkExceptMsg(rc,"Tp3CondVar.wait",msg));
  }
}
#endif
}  /* wait */

Function(SYSTEM_boolean ) P3THREADS_tp3condvar_DOT_timedwaitabs(
  P3THREADS_tp3condvar self,
  P3THREADS_tp3mutex mx,
  P3THREADS_tcvtime abstime)
{
  SYSTEM_boolean result;
  SYSTEM_shortstring msg;

  result = SYSTEM_true;
  /**** C code included from p3threads.pas(1403:1): 40 lines ****/
#if defined(_WIN32)
{
  SYSTEM_cardinal ticks;
  DWORD dw;
  CONDITION_VARIABLE *pCV = (CONDITION_VARIABLE *) &(FCV(self));
  CRITICAL_SECTION *pCS = (CRITICAL_SECTION *) &(FMUX(mx));

  ticks = P3THREADS_tickdeltacv(abstime);
  if (! SleepConditionVariableCS (pCV, pCS, ticks)) {
    dw = GetLastError();
    if ((ERROR_TIMEOUT == dw) || (WAIT_TIMEOUT == dw))
      result = SYSTEM_false;
    else {
      _P3CONDVAR_RAISE_E(mkExceptMsg(dw,"Tp3CondVar.timedWaitAbs",msg));
    }
  }
}
#else
{
  pthread_cond_t *pCV = (pthread_cond_t *) &(FCV(self));
  pthread_mutex_t *pMux = (pthread_mutex_t *) &(FMUX(mx));
  int rc;
  struct timespec ts;

  cvtime2timespec (abstime, &ts);
  rc = pthread_cond_timedwait (pCV, pMux, &ts);
  /*
  fprintf (stderr, "tv_sec  = %ld\n", (long) ts.tv_sec);
  fprintf (stderr, "tv_nsec = %ld\n", (long) ts.tv_nsec);
  */
  if (rc == ETIMEDOUT) {
    result = SYSTEM_false;
    return result;
  }
  if (rc) {
    _P3CONDVAR_RAISE_E(mkExceptMsg(rc,"Tp3CondVar.timedWaitAbs",msg));
  }
  return result;
}
#endif
  return result;
}  /* timedwaitabs */

Constructor(P3THREADS_tp3sharedmutex ) 
  P3THREADS_tp3sharedmutex_DOT_create(
  P3THREADS_tp3sharedmutex self)
{
  self->P3THREADS_tp3sharedmutex_DOT_statechange = ValueCast(
    P3THREADS_tp3mutex,P3THREADS_tp3mutex_DOT_create(ValueCast(
    P3THREADS_tp3mutex,_P3alloc_object(&P3THREADS_tp3mutex_CD))));
  self->P3THREADS_tp3sharedmutex_DOT_shared = ValueCast(
    P3THREADS_tp3condvar,P3THREADS_tp3condvar_DOT_create(ValueCast(
    P3THREADS_tp3condvar,_P3alloc_object(&P3THREADS_tp3condvar_CD))));
  self->P3THREADS_tp3sharedmutex_DOT_upgrade = ValueCast(
    P3THREADS_tp3condvar,P3THREADS_tp3condvar_DOT_create(ValueCast(
    P3THREADS_tp3condvar,_P3alloc_object(&P3THREADS_tp3condvar_CD))));
  self->P3THREADS_tp3sharedmutex_DOT_exclusive = ValueCast(
    P3THREADS_tp3condvar,P3THREADS_tp3condvar_DOT_create(ValueCast(
    P3THREADS_tp3condvar,_P3alloc_object(&P3THREADS_tp3condvar_CD))));
  return self;
}  /* create */

Destructor(P3THREADS_tp3sharedmutex ) 
  P3THREADS_tp3sharedmutex_DOT_destroy(
  P3THREADS_tp3sharedmutex self)
{
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    P3THREADS_tp3sharedmutex_DOT_statechange));
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    P3THREADS_tp3sharedmutex_DOT_shared));
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    P3THREADS_tp3sharedmutex_DOT_upgrade));
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    P3THREADS_tp3sharedmutex_DOT_exclusive));
  return self;
}  /* destroy */

Procedure P3THREADS_tp3sharedmutex_DOT_releasewaiters(
  P3THREADS_tp3sharedmutex self)
{
  self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusivewaitingblocked = 
    self->P3THREADS_tp3sharedmutex_DOT_fstate.upgradewaitingblocked;
  P3THREADS_tp3condvar_DOT_notifyone(self->
    P3THREADS_tp3sharedmutex_DOT_exclusive);
  P3THREADS_tp3condvar_DOT_notifyall(self->
    P3THREADS_tp3sharedmutex_DOT_shared);
}  /* releasewaiters */

Procedure P3THREADS_tp3sharedmutex_DOT_lock(
  P3THREADS_tp3sharedmutex self)
{
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  while (self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared > 0 || 
    self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive) {
    self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusivewaitingblocked = 
      SYSTEM_true;
    P3THREADS_tp3condvar_DOT_wait(self->
      P3THREADS_tp3sharedmutex_DOT_exclusive,self->
      P3THREADS_tp3sharedmutex_DOT_statechange);
  }
  self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive = SYSTEM_true;
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
}  /* lock */

Function(SYSTEM_boolean ) P3THREADS_tp3sharedmutex_DOT_trylock(
  P3THREADS_tp3sharedmutex self)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  while (self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared > 0 || 
    self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive) {
    P3THREADS_tp3mutex_DOT_unlock(self->
      P3THREADS_tp3sharedmutex_DOT_statechange);
    return result;
  }
  self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive = SYSTEM_true;
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  result = SYSTEM_true;
  return result;
}  /* trylock */

Procedure P3THREADS_tp3sharedmutex_DOT_unlock(
  P3THREADS_tp3sharedmutex self)
{
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive = SYSTEM_false;
  P3THREADS_tp3sharedmutex_DOT_releasewaiters(self);
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
}  /* unlock */

Procedure P3THREADS_tp3sharedmutex_DOT_lockshared(
  P3THREADS_tp3sharedmutex self)
{
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  while (self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive || self->
    P3THREADS_tp3sharedmutex_DOT_fstate.exclusivewaitingblocked) {

    P3THREADS_tp3condvar_DOT_wait(self->
      P3THREADS_tp3sharedmutex_DOT_shared,self->
      P3THREADS_tp3sharedmutex_DOT_statechange);
}
  _P3inc0(self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared);
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
}  /* lockshared */

Function(SYSTEM_boolean ) P3THREADS_tp3sharedmutex_DOT_trylockshared(
  P3THREADS_tp3sharedmutex self)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  while (self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive || self->
    P3THREADS_tp3sharedmutex_DOT_fstate.exclusivewaitingblocked) {
    P3THREADS_tp3mutex_DOT_unlock(self->
      P3THREADS_tp3sharedmutex_DOT_statechange);
    return result;
  }
  _P3inc0(self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared);
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  result = SYSTEM_true;
  return result;
}  /* trylockshared */

Procedure P3THREADS_tp3sharedmutex_DOT_unlockshared(
  P3THREADS_tp3sharedmutex self)
{
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  _P3dec0(self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared);
  if (0 == self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared) 
    if (self->P3THREADS_tp3sharedmutex_DOT_fstate.
      upgradewaitingblocked) {
      SYSTEM_assert(self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade,_P3str1("\000"));
      SYSTEM_assert(!self->P3THREADS_tp3sharedmutex_DOT_fstate.
        exclusive,_P3str1("\000"));
      SYSTEM_assert(self->P3THREADS_tp3sharedmutex_DOT_fstate.
        upgradewaitingblocked,_P3str1("\000"));
      self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade = SYSTEM_false;
      self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive = 
        SYSTEM_true;
      self->P3THREADS_tp3sharedmutex_DOT_fstate.upgradewaitingblocked = 
        SYSTEM_false;
      P3THREADS_tp3condvar_DOT_notifyone(self->
        P3THREADS_tp3sharedmutex_DOT_upgrade);
    } 
  P3THREADS_tp3sharedmutex_DOT_releasewaiters(self);
  if (0 == self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared) 
    SYSTEM_assert(!self->P3THREADS_tp3sharedmutex_DOT_fstate.
      exclusivewaitingblocked,_P3str1("\000"));
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
}  /* unlockshared */

Procedure P3THREADS_tp3sharedmutex_DOT_lockupgrade(
  P3THREADS_tp3sharedmutex self)
{
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  while (self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive || self->
    P3THREADS_tp3sharedmutex_DOT_fstate.exclusivewaitingblocked || 
    self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade) {

    P3THREADS_tp3condvar_DOT_wait(self->
      P3THREADS_tp3sharedmutex_DOT_shared,self->
      P3THREADS_tp3sharedmutex_DOT_statechange);
}
  _P3inc0(self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared);
  self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade = SYSTEM_true;
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
}  /* lockupgrade */

Function(SYSTEM_boolean ) P3THREADS_tp3sharedmutex_DOT_trylockupgrade(
  P3THREADS_tp3sharedmutex self)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  while (self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive || self->
    P3THREADS_tp3sharedmutex_DOT_fstate.exclusivewaitingblocked || 
    self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade) {
    P3THREADS_tp3mutex_DOT_unlock(self->
      P3THREADS_tp3sharedmutex_DOT_statechange);
    return result;
  }
  _P3inc0(self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared);
  self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade = SYSTEM_true;
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  result = SYSTEM_true;
  return result;
}  /* trylockupgrade */

Procedure P3THREADS_tp3sharedmutex_DOT_unlockupgrade(
  P3THREADS_tp3sharedmutex self)
{
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade = SYSTEM_false;
  _P3dec0(self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared);
  if (0 == self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared) 
    P3THREADS_tp3sharedmutex_DOT_releasewaiters(self);
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
}  /* unlockupgrade */

Procedure P3THREADS_tp3sharedmutex_DOT_unlockupgrade_and_lock(
  P3THREADS_tp3sharedmutex self)
{
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  _P3dec0(self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared);
  SYSTEM_assert(!self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive,_P3str1("\000"));
  SYSTEM_assert(self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade,_P3str1("\000"));
  while (self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared > 0) {
    self->P3THREADS_tp3sharedmutex_DOT_fstate.upgradewaitingblocked = 
      SYSTEM_true;
    self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusivewaitingblocked = 
      SYSTEM_true;
    SYSTEM_assert(self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade,_P3str1("\000"));
    P3THREADS_tp3condvar_DOT_wait(self->
      P3THREADS_tp3sharedmutex_DOT_upgrade,self->
      P3THREADS_tp3sharedmutex_DOT_statechange);
    SYSTEM_assert(!self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade,_P3str1("\000"));
    SYSTEM_assert(self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive,_P3str1("\000"));
  }
  self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade = SYSTEM_false;
  self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive = SYSTEM_true;
  SYSTEM_assert(self->P3THREADS_tp3sharedmutex_DOT_fstate.
    upgradewaitingblocked == SYSTEM_false,_P3str1("\000"));
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
}  /* unlockupgrade_and_lock */

Procedure P3THREADS_tp3sharedmutex_DOT_unlock_and_lockupgrade(
  P3THREADS_tp3sharedmutex self)
{
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive = SYSTEM_false;
  self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade = SYSTEM_true;
  _P3inc0(self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared);
  P3THREADS_tp3sharedmutex_DOT_releasewaiters(self);
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
}  /* unlock_and_lockupgrade */

Procedure P3THREADS_tp3sharedmutex_DOT_unlock_and_lockshared(
  P3THREADS_tp3sharedmutex self)
{
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  self->P3THREADS_tp3sharedmutex_DOT_fstate.exclusive = SYSTEM_false;
  _P3inc0(self->P3THREADS_tp3sharedmutex_DOT_fstate.nshared);
  P3THREADS_tp3sharedmutex_DOT_releasewaiters(self);
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
}  /* unlock_and_lockshared */

Procedure P3THREADS_tp3sharedmutex_DOT_unlockupgrade_and_lockshared(
  P3THREADS_tp3sharedmutex self)
{
  P3THREADS_tp3mutex_DOT_lock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
  self->P3THREADS_tp3sharedmutex_DOT_fstate.upgrade = SYSTEM_false;
  P3THREADS_tp3sharedmutex_DOT_releasewaiters(self);
  P3THREADS_tp3mutex_DOT_unlock(self->
    P3THREADS_tp3sharedmutex_DOT_statechange);
}  /* unlockupgrade_and_lockshared */
static SYSTEM_double P3THREADS_ticks_per_day = 24.0 * 60.0 * 60.0 * 1000.0;

Procedure P3THREADS_incdt(
  SYSTEM_P3_tdatetime *dt,
  SYSTEM_cardinal ticks)
{
  if (*dt < 0) 
    return;
  *dt = *dt + ticks /  P3THREADS_ticks_per_day;
}  /* incdt */

Procedure P3THREADS_decdt(
  SYSTEM_P3_tdatetime *dt,
  SYSTEM_cardinal ticks)
{
  if (*dt < 0) 
    return;
  *dt = *dt - ticks /  P3THREADS_ticks_per_day;
  if (*dt < 0) 
    *dt = 0;
}  /* decdt */

Function(P3THREADS_tcvtime ) P3THREADS_encodedatecv(
  SYSTEM_word year,
  SYSTEM_word month,
  SYSTEM_word day)
{
  P3THREADS_tcvtime result;
  SYSTEM_integer i;
  SYSTEM_integer ydays;
  SYSUTILS_P3_pdaytable daytable;

  result = 1;
  if (year < 1601 || year > 9999) 
    return result;
  if (month < 1 || month > 12) 
    return result;
  daytable = ValueCast(SYSUTILS_P3_pdaytable,SYSUTILS_P3_monthdays[
    SYSUTILS_P3_isleapyear(year)]);
  if (day < 1 || day > (*daytable)[month - 1]) 
    return result;
  result = ValueCast(SYSTEM_int32,day) - 1;
  { register SYSTEM_int32 _stop = ValueCast(SYSTEM_int32,month) - 1;
    if ((i = 1) <=  _stop) do {
      _P3inc1(result,(*daytable)[i - 1]);
    } while (i++ !=  _stop);

  }
  _P3dec1(year,1601);
  ydays = ValueCast(SYSTEM_int32,year) * 365;
  ydays = ydays + ValueCast(SYSTEM_int32,year) /  4 - ValueCast(
    SYSTEM_int32,year) /  100 + ValueCast(SYSTEM_int32,year) /  400;
  _P3inc1(result,ydays);
  result = result * P3THREADS_nicksperday;
  return result;
}  /* encodedatecv */

Function(P3THREADS_tcvtime ) P3THREADS_encodetimecv(
  SYSTEM_word hour,
  SYSTEM_word _min,
  SYSTEM_word sec,
  SYSTEM_word msec)
{
  P3THREADS_tcvtime result;

  result = 1;
  if (hour < 1 || hour > 24) 
    return result;
  if (_min < 1 || _min > 60) 
    return result;
  if (sec < 1 || sec > 60) 
    return result;
  if (msec < 1 || msec > 1000) 
    return result;
  result = (ValueCast(SYSTEM_int32,hour) - 1) * 60 + (ValueCast(
    SYSTEM_int32,_min) - 1);
  result = result * 60 + (ValueCast(SYSTEM_int32,sec) - 1);
  result = result * 1000 + (ValueCast(SYSTEM_int32,msec) - 1);
  result = result * P3THREADS_nicks_per_msec;
  return result;
}  /* encodetimecv */

Function(SYSTEM_boolean ) P3THREADS_decodedatefullycv(
  P3THREADS_tcvtime cvt,
  SYSTEM_word *year,
  SYSTEM_word *month,
  SYSTEM_word *day,
  SYSTEM_word *dow)
{
  SYSTEM_boolean result;
  SYSTEM_int64 d;
  SYSTEM_P3_tdatetime dt;

  result = SYSTEM_false;
  d = cvt /  P3THREADS_nicksperday;
  if (d < P3THREADS_date_delta_cv_dt) {
    *year = 1899;
    *month = 12;
    *day = 30;
    *dow = 0;
  } else {
    dt = d - P3THREADS_date_delta_cv_dt;
    result = SYSUTILS_P3_decodedatefully(dt,year,month,day,dow);
  }
  return result;
}  /* decodedatefullycv */

Procedure P3THREADS_decodedatecv(
  P3THREADS_tcvtime cvt,
  SYSTEM_word *year,
  SYSTEM_word *month,
  SYSTEM_word *day)
{
  SYSTEM_word dow;

  P3THREADS_decodedatefullycv(cvt,year,month,day,&dow);
}  /* decodedatecv */

Procedure P3THREADS_decodetimecv(
  P3THREADS_tcvtime cvt,
  SYSTEM_word *hour,
  SYSTEM_word *_min,
  SYSTEM_word *sec,
  SYSTEM_word *msec)
{
  SYSTEM_int64 nk;
  SYSTEM_int64 t;

  nk = cvt % P3THREADS_nicksperday;
  t = nk /  P3THREADS_nicksperhour;
  *hour = t;
  nk = nk % P3THREADS_nicksperhour;
  t = nk /  P3THREADS_nickspermin;
  *_min = t;
  nk = nk % P3THREADS_nickspermin;
  t = nk /  P3THREADS_nicks_per_sec;
  *sec = t;
  nk = nk % P3THREADS_nicks_per_sec;
  t = nk /  P3THREADS_nicks_per_msec;
  *msec = t;
}  /* decodetimecv */

Function(P3THREADS_tcvtime ) P3THREADS_nowcv(void)
{
  P3THREADS_tcvtime result;

  result = 0;
  /**** C code included from p3threads.pas(1775:1): 19 lines ****/
#if defined(_WIN32)
{
  cvTime_t cvt;

  GetSystemTimeAsFileTime (&(cvt.ft));
  result = cvt.i64;
}
#else
{
  struct timeval tv;
  int rc;

  rc = gettimeofday (&tv, NULL);
  if (rc)
    return result;
  result = (SYSTEM_int64)tv.tv_sec * nicks_per_sec + (SYSTEM_int64)tv.tv_usec * nicks_per_usec
         + nicksperday * date_delta_cv_unix;
}
#endif
  return result;
}  /* nowcv */

Function(P3THREADS_tcvtime ) P3THREADS_datecv(void)
{
  P3THREADS_tcvtime result;

  result = P3THREADS_nowcv();
  result = result /  P3THREADS_nicksperday;
  return result;
}  /* datecv */

Function(P3THREADS_tcvtime ) P3THREADS_timecv(void)
{
  P3THREADS_tcvtime result;

  result = P3THREADS_nowcv();
  result = result % P3THREADS_nicksperday;
  return result;
}  /* timecv */

Procedure P3THREADS_inccvtimemillis(
  P3THREADS_tcvtime *cvt,
  SYSTEM_cardinal ticks)
{
  _P3inc1(*cvt,ticks * P3THREADS_nicks_per_msec);
}  /* inccvtimemillis */

Procedure P3THREADS_deccvtimemillis(
  P3THREADS_tcvtime *cvt,
  SYSTEM_cardinal ticks)
{
  _P3dec1(*cvt,ticks * P3THREADS_nicks_per_msec);
  if (*cvt < 0) 
    *cvt = 0;
}  /* deccvtimemillis */

Function(SYSTEM_cardinal ) P3THREADS_tickdeltacv(
  P3THREADS_tcvtime cvt)
{
  SYSTEM_cardinal result;
  P3THREADS_tcvtime d, r;

  result = 0;
  d = cvt - P3THREADS_nowcv();
  if (d < 0) 
    return result;
  result = d /  P3THREADS_nicks_per_msec;
  r = d % P3THREADS_nicks_per_msec;
  if (r > 0) 
    _P3inc0(result);
  return result;
}  /* tickdeltacv */

/* unit p3threads */
void _Init_Module_p3threads(void)
{
  P3THREADS_nickspermin = 60;
  P3THREADS_nickspermin = P3THREADS_nickspermin * 
    P3THREADS_nicks_per_sec;
  P3THREADS_nicksperhour = P3THREADS_nickspermin * 60;
  P3THREADS_nicksperday = P3THREADS_nicksperhour * 24;
} /* _Init_Module_p3threads */

void _Final_Module_p3threads(void)
{
} /* _Final_Module_p3threads */

