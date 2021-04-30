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
#include "gmsspecs.h"
#include "gmsgen.h"
#include "strutilx.h"
#include "gmsobj.h"
#include "gmsheapnew.h"


void* const GMSHEAPNEW_tbigblockmgr_VT[] = {(void*) &GMSHEAPNEW_tbigblockmgr_DOT_destroy};

/* Class descriptor for 'tbigblockmgr' */
const SYSTEM_classdescriptor_t GMSHEAPNEW_tbigblockmgr_CD = {_P3str1("\014tbigblockmgr"), &SYSTEM_tobject_CD, NULL, 0,
                                                             sizeof(GMSHEAPNEW_tbigblockmgr_OD), GMSHEAPNEW_tbigblockmgr_VT, NULL};


void* const GMSHEAPNEW_theapmgr_VT[] = {(void*) &GMSHEAPNEW_theapmgr_DOT_destroy};

/* Class descriptor for 'theapmgr' */
const SYSTEM_classdescriptor_t GMSHEAPNEW_theapmgr_CD = {_P3str1("\010theapmgr"), &SYSTEM_tobject_CD, NULL, 0, sizeof(GMSHEAPNEW_theapmgr_OD),
                                                         GMSHEAPNEW_theapmgr_VT, NULL};


Constructor(GMSHEAPNEW_tbigblockmgr) GMSHEAPNEW_tbigblockmgr_DOT_create(GMSHEAPNEW_tbigblockmgr self, const SYSTEM_ansichar* name) {
   ValueCast(GMSHEAPNEW_tbigblockmgr, SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject, self)));
   self->GMSHEAPNEW_tbigblockmgr_DOT_memoryreportproc = NULL;
   _P3getmem(self->GMSHEAPNEW_tbigblockmgr_DOT_spname, ValueCast(SYSTEM_int32, SYSTEM_length(name)) + 1);
   _P3strcpy(*self->GMSHEAPNEW_tbigblockmgr_DOT_spname, 255, name);
   self->GMSHEAPNEW_tbigblockmgr_DOT_freelist = ValueCast(GMSOBJ_txlist,
         SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject, _P3alloc_object(&GMSOBJ_txlist_CD))));
   self->GMSHEAPNEW_tbigblockmgr_DOT_memorylimit = 1e200;
   self->GMSHEAPNEW_tbigblockmgr_DOT_totalhighmark = 0.0;
   self->GMSHEAPNEW_tbigblockmgr_DOT_totalmemory = 0.0;
   self->GMSHEAPNEW_tbigblockmgr_DOT_showosmem = 0;
   return self;
}  /* create */

Destructor(GMSHEAPNEW_tbigblockmgr) GMSHEAPNEW_tbigblockmgr_DOT_destroy(GMSHEAPNEW_tbigblockmgr self) {
   GMSHEAPNEW_tbigblockmgr_DOT_xclear(self);
   _P3freemem(self->GMSHEAPNEW_tbigblockmgr_DOT_spname);
   SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject, self->GMSHEAPNEW_tbigblockmgr_DOT_freelist));
   SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject, self));
   return self;
}  /* destroy */

Procedure GMSHEAPNEW_tbigblockmgr_DOT_xclear(GMSHEAPNEW_tbigblockmgr self) {
   SYSTEM_integer n;

   {
      register SYSTEM_int32 _stop = self->GMSHEAPNEW_tbigblockmgr_DOT_freelist->GMSOBJ_txlist_DOT_fcount - 1;
      if ((n = 0) <= _stop)
         do {
            _P3freemem2(GMSOBJ_txlist_DOT_get(self->GMSHEAPNEW_tbigblockmgr_DOT_freelist, n), GMSHEAPNEW_bigblocksize);
         } while (n++ != _stop);

   }
   GMSHEAPNEW_tbigblockmgr_DOT_reducememorysize(self, self->GMSHEAPNEW_tbigblockmgr_DOT_freelist->GMSOBJ_txlist_DOT_fcount * GMSHEAPNEW_bigblocksize);
   GMSOBJ_txlist_DOT_clear(self->GMSHEAPNEW_tbigblockmgr_DOT_freelist);
}  /* xclear */

Procedure GMSHEAPNEW_tbigblockmgr_DOT_reducememorysize(GMSHEAPNEW_tbigblockmgr self, SYSTEM_int64 delta) {
   self->GMSHEAPNEW_tbigblockmgr_DOT_othermemory = self->GMSHEAPNEW_tbigblockmgr_DOT_othermemory - delta;
   self->GMSHEAPNEW_tbigblockmgr_DOT_totalmemory = self->GMSHEAPNEW_tbigblockmgr_DOT_totalmemory - delta;
   if (_P3assigned(self->GMSHEAPNEW_tbigblockmgr_DOT_memoryreportproc))
      (*self->GMSHEAPNEW_tbigblockmgr_DOT_memoryreportproc)(GMSHEAPNEW_tbigblockmgr_DOT_memoryusedmb(self));
}  /* reducememorysize */

Procedure GMSHEAPNEW_tbigblockmgr_DOT_increasememorysize(GMSHEAPNEW_tbigblockmgr self, SYSTEM_int64 delta) {
   if (self->GMSHEAPNEW_tbigblockmgr_DOT_totalmemory + delta > self->GMSHEAPNEW_tbigblockmgr_DOT_memorylimit)
      _P3_RAISE(ValueCast(EXCEPTIONS_eoutofmemory,
            SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception, _P3alloc_object(&EXCEPTIONS_eoutofmemory_CD)),
                  _P3str1("\053Requested memory exceeds assigned HeapLimit"))));
   self->GMSHEAPNEW_tbigblockmgr_DOT_othermemory = self->GMSHEAPNEW_tbigblockmgr_DOT_othermemory + delta;
   if (self->GMSHEAPNEW_tbigblockmgr_DOT_othermemory > self->GMSHEAPNEW_tbigblockmgr_DOT_highmark)
      self->GMSHEAPNEW_tbigblockmgr_DOT_highmark = self->GMSHEAPNEW_tbigblockmgr_DOT_othermemory;
   self->GMSHEAPNEW_tbigblockmgr_DOT_totalmemory = self->GMSHEAPNEW_tbigblockmgr_DOT_totalmemory + delta;
   if (self->GMSHEAPNEW_tbigblockmgr_DOT_totalmemory > self->GMSHEAPNEW_tbigblockmgr_DOT_totalhighmark)
      self->GMSHEAPNEW_tbigblockmgr_DOT_totalhighmark = self->GMSHEAPNEW_tbigblockmgr_DOT_totalmemory;
   if (_P3assigned(self->GMSHEAPNEW_tbigblockmgr_DOT_memoryreportproc))
      (*self->GMSHEAPNEW_tbigblockmgr_DOT_memoryreportproc)(GMSHEAPNEW_tbigblockmgr_DOT_memoryusedmb(self));
}  /* increasememorysize */

Function(SYSTEM_double) GMSHEAPNEW_tbigblockmgr_DOT_memorylimitmb(GMSHEAPNEW_tbigblockmgr self) {
   SYSTEM_double result;

   result = self->GMSHEAPNEW_tbigblockmgr_DOT_memorylimit / 1e6;
   return result;
}  /* memorylimitmb */

Function(SYSTEM_double) GMSHEAPNEW_tbigblockmgr_DOT_memoryusedmb(GMSHEAPNEW_tbigblockmgr self) {
   SYSTEM_double result;
   SYSTEM_int64 rss, vss;

   result = self->GMSHEAPNEW_tbigblockmgr_DOT_totalmemory / 1e6;
   if (self->GMSHEAPNEW_tbigblockmgr_DOT_showosmem == 1 && P3UTILS_p3getmemoryinfo(&rss, &vss)) {
      result = rss / 1e6;
   }
   else if (self->GMSHEAPNEW_tbigblockmgr_DOT_showosmem == 2 && P3UTILS_p3getmemoryinfo(&rss, &vss))
      result = vss / 1e6;
   return result;
}  /* memoryusedmb */

Function(SYSTEM_pointer) GMSHEAPNEW_tbigblockmgr_DOT_getbigblock(GMSHEAPNEW_tbigblockmgr self) {
   SYSTEM_pointer result;

   result = GMSOBJ_txlist_DOT_getlast(self->GMSHEAPNEW_tbigblockmgr_DOT_freelist);
   if (result != NULL) {
      GMSOBJ_txlist_DOT_delete(self->GMSHEAPNEW_tbigblockmgr_DOT_freelist, self->GMSHEAPNEW_tbigblockmgr_DOT_freelist->GMSOBJ_txlist_DOT_fcount - 1);
   }
   else {
      GMSHEAPNEW_tbigblockmgr_DOT_increasememorysize(self, GMSHEAPNEW_bigblocksize);
      _P3getmem(result, GMSHEAPNEW_bigblocksize);
   }
   return result;
}  /* getbigblock */

Procedure GMSHEAPNEW_tbigblockmgr_DOT_releasebigblock(GMSHEAPNEW_tbigblockmgr self, SYSTEM_pointer p) {
   GMSOBJ_txlist_DOT_add(self->GMSHEAPNEW_tbigblockmgr_DOT_freelist, p);
}  /* releasebigblock */

Function(SYSTEM_ansichar *) GMSHEAPNEW_tbigblockmgr_DOT_getname(SYSTEM_ansichar* result, SYSTEM_uint8 _len_ret, GMSHEAPNEW_tbigblockmgr self) {
   _P3strcpy(result, _len_ret, *self->GMSHEAPNEW_tbigblockmgr_DOT_spname);
   return result;
}  /* getname */

Procedure GMSHEAPNEW_tbigblockmgr_DOT_getbigstats(GMSHEAPNEW_tbigblockmgr self, SYSTEM_int64* sizeothermemory, SYSTEM_int64* sizehighmark,
      SYSTEM_int64* cntfree) {
   *sizeothermemory = self->GMSHEAPNEW_tbigblockmgr_DOT_othermemory;
   *sizehighmark = self->GMSHEAPNEW_tbigblockmgr_DOT_highmark;
   *cntfree = self->GMSHEAPNEW_tbigblockmgr_DOT_freelist->GMSOBJ_txlist_DOT_fcount;
}  /* getbigstats */

Procedure GMSHEAPNEW_theapmgr_DOT_prvclear(GMSHEAPNEW_theapmgr self) {
   SYSTEM_integer n;
   GMSHEAPNEW_theapslotnr slot;

   while (self->GMSHEAPNEW_theapmgr_DOT_wrkbuffs->GMSOBJ_txlist_DOT_fcount > 0) {

      GMSHEAPNEW_theapmgr_DOT_releaseworkbuffer(self, ValueCast(GMSHEAPNEW_plargeblock,
            GMSOBJ_txlist_DOT_get(self->GMSHEAPNEW_theapmgr_DOT_wrkbuffs, self->GMSHEAPNEW_theapmgr_DOT_wrkbuffs->GMSOBJ_txlist_DOT_fcount - 1)));
   }
   GMSOBJ_txlist_DOT_clear(self->GMSHEAPNEW_theapmgr_DOT_wrkbuffs);
   self->GMSHEAPNEW_theapmgr_DOT_workbuffer = NULL;
   {
      register SYSTEM_int32 _stop = self->GMSHEAPNEW_theapmgr_DOT_active->GMSOBJ_txlist_DOT_fcount - 1;
      if ((n = 0) <= _stop)
         do {
            _P3freemem0(GMSOBJ_txlist_DOT_get(self->GMSHEAPNEW_theapmgr_DOT_active, n));
         } while (n++ != _stop);

   }
   GMSOBJ_txlist_DOT_clear(self->GMSHEAPNEW_theapmgr_DOT_active);
   for (slot = 1; slot <= (SYSTEM_uint8) GMSHEAPNEW_lastslot; ++slot) {
      {
         register GMSHEAPNEW_tslotrecord* _W2 = &self->GMSHEAPNEW_theapmgr_DOT_slots[slot - 1];
         _W2->getcount = 0;
         _W2->freecount = 0;
         _W2->listcount = 0;
         _W2->firstfree = NULL;

      }
   }
   GMSHEAPNEW_theapmgr_DOT_reducememorysize(self, self->GMSHEAPNEW_theapmgr_DOT_othermemory);
   self->GMSHEAPNEW_theapmgr_DOT_otherget = 0;
   self->GMSHEAPNEW_theapmgr_DOT_otherfree = 0;
   self->GMSHEAPNEW_theapmgr_DOT_otherget64 = 0;
   self->GMSHEAPNEW_theapmgr_DOT_otherfree64 = 0;
   self->GMSHEAPNEW_theapmgr_DOT_realloccnt = 0;
   self->GMSHEAPNEW_theapmgr_DOT_reallocused = 0;
   self->GMSHEAPNEW_theapmgr_DOT_realloccnt64 = 0;
   self->GMSHEAPNEW_theapmgr_DOT_reallocused64 = 0;
}  /* prvclear */

Procedure GMSHEAPNEW_theapmgr_DOT_clear(GMSHEAPNEW_theapmgr self) {
   GMSHEAPNEW_theapmgr_DOT_prvclear(self);
}  /* clear */

Constructor(GMSHEAPNEW_theapmgr) GMSHEAPNEW_theapmgr_DOT_create(GMSHEAPNEW_theapmgr self, const SYSTEM_ansichar* name) {
   ValueCast(GMSHEAPNEW_theapmgr, SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject, self)));
   {
      _P3STR_255 _t1;

      self->GMSHEAPNEW_theapmgr_DOT_blockmgr = ValueCast(GMSHEAPNEW_tbigblockmgr,
            GMSHEAPNEW_tbigblockmgr_DOT_create(ValueCast(GMSHEAPNEW_tbigblockmgr, _P3alloc_object(&GMSHEAPNEW_tbigblockmgr_CD)),
                  _P3strcat(_t1, 255, _P3str1("\006BBMgr_"), name)));
   }
   _P3getmem(self->GMSHEAPNEW_theapmgr_DOT_spname, ValueCast(SYSTEM_int32, SYSTEM_length(name)) + 1);
   _P3strcpy(*self->GMSHEAPNEW_theapmgr_DOT_spname, 255, name);
   self->GMSHEAPNEW_theapmgr_DOT_wrkbuffs = ValueCast(GMSOBJ_txlist,
         SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject, _P3alloc_object(&GMSOBJ_txlist_CD))));
   self->GMSHEAPNEW_theapmgr_DOT_active = ValueCast(GMSOBJ_txlist,
         SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject, _P3alloc_object(&GMSOBJ_txlist_CD))));
   GMSHEAPNEW_theapmgr_DOT_prvclear(self);
   return self;
}  /* create */

Destructor(GMSHEAPNEW_theapmgr) GMSHEAPNEW_theapmgr_DOT_destroy(GMSHEAPNEW_theapmgr self) {
   GMSHEAPNEW_theapmgr_DOT_prvclear(self);
   SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject, self->GMSHEAPNEW_theapmgr_DOT_wrkbuffs));
   SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject, self->GMSHEAPNEW_theapmgr_DOT_active));
   _P3freemem(self->GMSHEAPNEW_theapmgr_DOT_spname);
   SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject, self->GMSHEAPNEW_theapmgr_DOT_blockmgr));
   SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject, self));
   return self;
}  /* destroy */

Function(GMSHEAPNEW_plargeblock) GMSHEAPNEW_theapmgr_DOT_getworkbuffer(GMSHEAPNEW_theapmgr self) {
   GMSHEAPNEW_plargeblock result;

   _P3getmem(result, sizeof(GMSHEAPNEW_tlargeblock));
   {
      register GMSHEAPNEW_tlargeblock* _W2 = result;
      _W2->freeslots = 65536;
      _W2->initialptr = ValueCast(SYSTEM_P3_pbyte, GMSHEAPNEW_tbigblockmgr_DOT_getbigblock(self->GMSHEAPNEW_theapmgr_DOT_blockmgr));
      _W2->currptr = _W2->initialptr;

   }
   GMSOBJ_txlist_DOT_add(self->GMSHEAPNEW_theapmgr_DOT_wrkbuffs, result);
   return result;
}  /* getworkbuffer */

Procedure GMSHEAPNEW_theapmgr_DOT_releaseworkbuffer(GMSHEAPNEW_theapmgr self, GMSHEAPNEW_plargeblock p) {
   GMSHEAPNEW_tbigblockmgr_DOT_releasebigblock(self->GMSHEAPNEW_theapmgr_DOT_blockmgr, p->initialptr);
   GMSOBJ_txlist_DOT_remove(self->GMSHEAPNEW_theapmgr_DOT_wrkbuffs, p);
   _P3freemem(p);
}  /* releaseworkbuffer */

Procedure GMSHEAPNEW_theapmgr_DOT_prvgmsfreemem(GMSHEAPNEW_theapmgr self, SYSTEM_pointer p, SYSTEM_word slot) {
   {
      register GMSHEAPNEW_tslotrecord* _W2 = &self->GMSHEAPNEW_theapmgr_DOT_slots[slot - 1];
      _P3inc0(_W2->freecount);
      _P3inc0(_W2->listcount);
      (ValueCast(GMSHEAPNEW_psmallblock, p))->nextsmallblock = _W2->firstfree;
      _W2->firstfree = ValueCast(GMSHEAPNEW_psmallblock, p);

   }
}  /* prvgmsfreemem */

Procedure GMSHEAPNEW_theapmgr_DOT_gmsfreemem(GMSHEAPNEW_theapmgr self, SYSTEM_pointer p, SYSTEM_word slot) {
   GMSHEAPNEW_theapmgr_DOT_prvgmsfreemem(self, p, slot);
}  /* gmsfreemem */

Function(SYSTEM_pointer) GMSHEAPNEW_theapmgr_DOT_prvgmsgetmem(GMSHEAPNEW_theapmgr self, SYSTEM_word slot) {
   SYSTEM_pointer result;

   {
      register GMSHEAPNEW_tslotrecord* _W2 = &self->GMSHEAPNEW_theapmgr_DOT_slots[slot - 1];
      _P3inc0(_W2->getcount);
      result = _W2->firstfree;
      if (result != NULL) {
         _W2->firstfree = (ValueCast(GMSHEAPNEW_psmallblock, result))->nextsmallblock;
         _P3dec0(_W2->listcount);
         return result;
      }

   }
   if (self->GMSHEAPNEW_theapmgr_DOT_workbuffer == NULL)
      self->GMSHEAPNEW_theapmgr_DOT_workbuffer = GMSHEAPNEW_theapmgr_DOT_getworkbuffer(self);
   {
      register GMSHEAPNEW_tlargeblock* _W2 = self->GMSHEAPNEW_theapmgr_DOT_workbuffer;
      if (_W2->freeslots >= ValueCast(SYSTEM_int32, slot)) {
         result = _W2->currptr;
         _P3inc1(_W2->currptr, ValueCast(SYSTEM_int32, slot) * GMSHEAPNEW_heapgranularity);
         _P3dec1(_W2->freeslots, slot);
         return result;
      }
      if (_W2->freeslots > 0) {
         register GMSHEAPNEW_tslotrecord* _W3 = &self->GMSHEAPNEW_theapmgr_DOT_slots[_W2->freeslots - 1];
         _P3inc0(_W3->listcount);
         (ValueCast(GMSHEAPNEW_psmallblock, _W2->currptr))->nextsmallblock = _W3->firstfree;
         _W3->firstfree = ValueCast(GMSHEAPNEW_psmallblock, _W2->currptr);

      }

   }
   self->GMSHEAPNEW_theapmgr_DOT_workbuffer = GMSHEAPNEW_theapmgr_DOT_getworkbuffer(self);
   {
      register GMSHEAPNEW_tlargeblock* _W2 = self->GMSHEAPNEW_theapmgr_DOT_workbuffer;
      result = _W2->currptr;
      _P3inc1(_W2->currptr, ValueCast(SYSTEM_int32, slot) * GMSHEAPNEW_heapgranularity);
      _P3dec1(_W2->freeslots, slot);

   }
   return result;
}  /* prvgmsgetmem */

Function(SYSTEM_pointer) GMSHEAPNEW_theapmgr_DOT_gmsgetmem(GMSHEAPNEW_theapmgr self, SYSTEM_word slot) {
   SYSTEM_pointer result;

   result = GMSHEAPNEW_theapmgr_DOT_prvgmsgetmem(self, slot);
   return result;
}  /* gmsgetmem */

Procedure GMSHEAPNEW_theapmgr_DOT_prvxfreemem(GMSHEAPNEW_theapmgr self, SYSTEM_pointer p, SYSTEM_integer size) {
   if (size > 0)
      if (size <= 256) {
         GMSHEAPNEW_theapmgr_DOT_prvgmsfreemem(self, p, (size - 1) / GMSHEAPNEW_heapgranularity + 1);
      }
      else {
         _P3inc0(self->GMSHEAPNEW_theapmgr_DOT_otherfree);
         GMSOBJ_txlist_DOT_remove(self->GMSHEAPNEW_theapmgr_DOT_active, p);
         GMSHEAPNEW_theapmgr_DOT_reducememorysize(self, size);
         _P3freemem0(p);
      }
}  /* prvxfreemem */

Procedure GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_theapmgr self, SYSTEM_pointer p, SYSTEM_integer size) {
   GMSHEAPNEW_theapmgr_DOT_prvxfreemem(self, p, size);
}  /* xfreemem */

Procedure GMSHEAPNEW_theapmgr_DOT_xfreememnc(GMSHEAPNEW_theapmgr self, SYSTEM_pointer p, SYSTEM_integer size) {
   if (size > 0) {
      _P3inc0(self->GMSHEAPNEW_theapmgr_DOT_otherfree);
      GMSHEAPNEW_theapmgr_DOT_reducememorysize(self, size);
      _P3freemem0(p);
   }
}  /* xfreememnc */

Procedure GMSHEAPNEW_theapmgr_DOT_xfreememandnil(GMSHEAPNEW_theapmgr self, SYSTEM_pointer* p, SYSTEM_integer size) {
   GMSHEAPNEW_theapmgr_DOT_prvxfreemem(self, *p, size);
   *p = NULL;
}  /* xfreememandnil */

Procedure GMSHEAPNEW_theapmgr_DOT_xfreemem64andnil(GMSHEAPNEW_theapmgr self, SYSTEM_pointer* p, SYSTEM_int64 size) {
   GMSHEAPNEW_theapmgr_DOT_prvxfreemem64(self, *p, size);
   *p = NULL;
}  /* xfreemem64andnil */

Procedure GMSHEAPNEW_theapmgr_DOT_prvxfreemem64(GMSHEAPNEW_theapmgr self, SYSTEM_pointer p, SYSTEM_int64 size) {
   if (size > 0)
      if (size <= 256) {
         GMSHEAPNEW_theapmgr_DOT_prvgmsfreemem(self, p, (size - 1) / GMSHEAPNEW_heapgranularity + 1);
      }
      else {
         _P3inc0(self->GMSHEAPNEW_theapmgr_DOT_otherfree64);
         GMSOBJ_txlist_DOT_remove(self->GMSHEAPNEW_theapmgr_DOT_active, p);
         GMSHEAPNEW_theapmgr_DOT_reducememorysize(self, size);
         P3UTILS_p3freemem64(&p, size);
      }
}  /* prvxfreemem64 */

Procedure GMSHEAPNEW_theapmgr_DOT_xfreemem64(GMSHEAPNEW_theapmgr self, SYSTEM_pointer p, SYSTEM_int64 size) {
   GMSHEAPNEW_theapmgr_DOT_prvxfreemem64(self, p, size);
}  /* xfreemem64 */

Function(SYSTEM_pointer) GMSHEAPNEW_theapmgr_DOT_prvxgetmem(GMSHEAPNEW_theapmgr self, SYSTEM_integer size) {
   SYSTEM_pointer result;

   if (size <= 0) {
      result = NULL;
   }
   else if (size <= 256) {
      result = GMSHEAPNEW_theapmgr_DOT_prvgmsgetmem(self, (size - 1) / GMSHEAPNEW_heapgranularity + 1);
   }
   else {
      _P3inc0(self->GMSHEAPNEW_theapmgr_DOT_otherget);
      GMSHEAPNEW_theapmgr_DOT_increasememorysize(self, size);
      _P3getmem(result, size);
      GMSOBJ_txlist_DOT_add(self->GMSHEAPNEW_theapmgr_DOT_active, result);
   }
   return result;
}  /* prvxgetmem */

Function(SYSTEM_pointer) GMSHEAPNEW_theapmgr_DOT_xgetmem(GMSHEAPNEW_theapmgr self, SYSTEM_integer size) {
   SYSTEM_pointer result;

   result = GMSHEAPNEW_theapmgr_DOT_prvxgetmem(self, size);
   return result;
}  /* xgetmem */

Function(SYSTEM_pointer) GMSHEAPNEW_theapmgr_DOT_prvxgetmemnc(GMSHEAPNEW_theapmgr self, SYSTEM_integer size) {
   SYSTEM_pointer result;

   if (size <= 0) {
      result = NULL;
   }
   else {
      _P3inc0(self->GMSHEAPNEW_theapmgr_DOT_otherget);
      GMSHEAPNEW_theapmgr_DOT_increasememorysize(self, size);
      _P3getmem(result, size);
   }
   return result;
}  /* prvxgetmemnc */

Function(SYSTEM_pointer) GMSHEAPNEW_theapmgr_DOT_xgetmemnc(GMSHEAPNEW_theapmgr self, SYSTEM_integer size) {
   SYSTEM_pointer result;

   result = GMSHEAPNEW_theapmgr_DOT_prvxgetmemnc(self, size);
   return result;
}  /* xgetmemnc */

Function(SYSTEM_pointer) GMSHEAPNEW_theapmgr_DOT_xallocmem(GMSHEAPNEW_theapmgr self, SYSTEM_integer size) {
   SYSTEM_pointer result;

   result = GMSHEAPNEW_theapmgr_DOT_prvxgetmem(self, size);
   if (result != NULL)
      SYSTEM_P3_fillchar(ValueCast(SYSTEM_P3_pbyte, result), size, 0);
   return result;
}  /* xallocmem */

Function(SYSTEM_pointer) GMSHEAPNEW_theapmgr_DOT_xallocmemnc(GMSHEAPNEW_theapmgr self, SYSTEM_integer size) {
   SYSTEM_pointer result;

   result = GMSHEAPNEW_theapmgr_DOT_prvxgetmemnc(self, size);
   if (result != NULL)
      SYSTEM_P3_fillchar(ValueCast(SYSTEM_P3_pbyte, result), size, 0);
   return result;
}  /* xallocmemnc */

Function(SYSTEM_pointer) GMSHEAPNEW_theapmgr_DOT_prvxgetmem64(GMSHEAPNEW_theapmgr self, SYSTEM_int64 size) {
   SYSTEM_pointer result;

   if (size <= 0) {
      result = NULL;
   }
   else if (size <= 256) {
      result = GMSHEAPNEW_theapmgr_DOT_prvgmsgetmem(self, (size - 1) / GMSHEAPNEW_heapgranularity + 1);
   }
   else {
      _P3inc0(self->GMSHEAPNEW_theapmgr_DOT_otherget64);
      GMSHEAPNEW_theapmgr_DOT_increasememorysize(self, size);
      P3UTILS_p3getmem64(&result, size);
      GMSOBJ_txlist_DOT_add(self->GMSHEAPNEW_theapmgr_DOT_active, result);
   }
   return result;
}  /* prvxgetmem64 */

Function(SYSTEM_pointer) GMSHEAPNEW_theapmgr_DOT_xgetmem64(GMSHEAPNEW_theapmgr self, SYSTEM_int64 size) {
   SYSTEM_pointer result;

   result = GMSHEAPNEW_theapmgr_DOT_prvxgetmem64(self, size);
   return result;
}  /* xgetmem64 */

Function(SYSTEM_pointer) GMSHEAPNEW_theapmgr_DOT_xallocmem64(GMSHEAPNEW_theapmgr self, SYSTEM_int64 size) {
   SYSTEM_pointer result;

   result = GMSHEAPNEW_theapmgr_DOT_prvxgetmem64(self, size);
   if (result != NULL)
      P3UTILS_p3fillchar64(ValueCast(SYSTEM_P3_pbyte, result), size, 0);
   return result;
}  /* xallocmem64 */

Procedure GMSHEAPNEW_theapmgr_DOT_xreallocmem(GMSHEAPNEW_theapmgr self, SYSTEM_pointer* p, SYSTEM_integer oldsize, SYSTEM_integer newsize) {
   SYSTEM_pointer pnew;

   _P3inc0(self->GMSHEAPNEW_theapmgr_DOT_realloccnt);
   _P3dec1(self->GMSHEAPNEW_theapmgr_DOT_reallocused, oldsize);
   _P3inc1(self->GMSHEAPNEW_theapmgr_DOT_reallocused, newsize);
   if (newsize <= 0) {
      if (oldsize > 0 && *p != NULL)
         GMSHEAPNEW_theapmgr_DOT_prvxfreemem(self, *p, oldsize);
      pnew = NULL;
   }
   else if (*p == NULL || oldsize <= 0) {
      pnew = GMSHEAPNEW_theapmgr_DOT_prvxgetmem(self, newsize);
   }
   else if (oldsize == newsize) {
      pnew = *p;
   }
   else if (oldsize > 256 && newsize > 256) {
      pnew = *p;
      GMSOBJ_txlist_DOT_remove(self->GMSHEAPNEW_theapmgr_DOT_active, *p);
      SYSTEM_reallocmem(&pnew, newsize);
      GMSOBJ_txlist_DOT_add(self->GMSHEAPNEW_theapmgr_DOT_active, pnew);
      if (newsize > oldsize) {
         GMSHEAPNEW_theapmgr_DOT_increasememorysize(self, newsize - oldsize);
      }
      else
         GMSHEAPNEW_theapmgr_DOT_reducememorysize(self, oldsize - newsize);
   }
   else {
      pnew = GMSHEAPNEW_theapmgr_DOT_prvxgetmem(self, newsize);
      if (oldsize <= newsize) {
         SYSTEM_move(ValueCast(SYSTEM_P3_pbyte, *p), ValueCast(SYSTEM_P3_pbyte, pnew), oldsize);
      }
      else
         SYSTEM_move(ValueCast(SYSTEM_P3_pbyte, *p), ValueCast(SYSTEM_P3_pbyte, pnew), newsize);
      GMSHEAPNEW_theapmgr_DOT_prvxfreemem(self, *p, oldsize);
   }
   *p = pnew;
}  /* xreallocmem */

Procedure GMSHEAPNEW_theapmgr_DOT_xreallocmemnc(GMSHEAPNEW_theapmgr self, SYSTEM_pointer* p, SYSTEM_integer oldsize, SYSTEM_integer newsize) {
   _P3inc0(self->GMSHEAPNEW_theapmgr_DOT_realloccnt);
   _P3dec1(self->GMSHEAPNEW_theapmgr_DOT_reallocused, oldsize);
   _P3inc1(self->GMSHEAPNEW_theapmgr_DOT_reallocused, newsize);
   SYSTEM_reallocmem(p, newsize);
}  /* xreallocmemnc */

Procedure GMSHEAPNEW_theapmgr_DOT_xreallocmem64(GMSHEAPNEW_theapmgr self, SYSTEM_pointer* p, SYSTEM_int64 oldsize, SYSTEM_int64 newsize) {
   SYSTEM_pointer pnew;

   _P3inc0(self->GMSHEAPNEW_theapmgr_DOT_realloccnt64);
   _P3dec1(self->GMSHEAPNEW_theapmgr_DOT_reallocused64, oldsize);
   _P3inc1(self->GMSHEAPNEW_theapmgr_DOT_reallocused64, newsize);
   if (newsize <= 0) {
      if (oldsize > 0 && *p != NULL)
         GMSHEAPNEW_theapmgr_DOT_prvxfreemem64(self, *p, oldsize);
      pnew = NULL;
   }
   else if (*p == NULL || oldsize <= 0) {
      pnew = GMSHEAPNEW_theapmgr_DOT_prvxgetmem64(self, newsize);
   }
   else if (oldsize == newsize) {
      pnew = *p;
   }
   else if (oldsize > 256 && newsize > 256) {
      pnew = *p;
      GMSOBJ_txlist_DOT_remove(self->GMSHEAPNEW_theapmgr_DOT_active, *p);
      P3UTILS_p3reallocmem64(&pnew, newsize);
      GMSOBJ_txlist_DOT_add(self->GMSHEAPNEW_theapmgr_DOT_active, pnew);
      if (newsize > oldsize) {
         GMSHEAPNEW_theapmgr_DOT_increasememorysize(self, newsize - oldsize);
      }
      else
         GMSHEAPNEW_theapmgr_DOT_reducememorysize(self, oldsize - newsize);
   }
   else {
      pnew = GMSHEAPNEW_theapmgr_DOT_prvxgetmem64(self, newsize);
      if (oldsize <= newsize) {
         SYSTEM_move(ValueCast(SYSTEM_P3_pbyte, *p), ValueCast(SYSTEM_P3_pbyte, pnew), oldsize);
      }
      else
         SYSTEM_move(ValueCast(SYSTEM_P3_pbyte, *p), ValueCast(SYSTEM_P3_pbyte, pnew), newsize);
      GMSHEAPNEW_theapmgr_DOT_prvxfreemem64(self, *p, oldsize);
   }
   *p = pnew;
}  /* xreallocmem64 */

Procedure GMSHEAPNEW_theapmgr_DOT_increasememorysize(GMSHEAPNEW_theapmgr self, SYSTEM_int64 delta) {
   GMSHEAPNEW_tbigblockmgr_DOT_increasememorysize(self->GMSHEAPNEW_theapmgr_DOT_blockmgr, delta);
   self->GMSHEAPNEW_theapmgr_DOT_othermemory = self->GMSHEAPNEW_theapmgr_DOT_othermemory + delta;
   if (self->GMSHEAPNEW_theapmgr_DOT_othermemory > self->GMSHEAPNEW_theapmgr_DOT_highmark)
      self->GMSHEAPNEW_theapmgr_DOT_highmark = self->GMSHEAPNEW_theapmgr_DOT_othermemory;
}  /* increasememorysize */

Procedure GMSHEAPNEW_theapmgr_DOT_reducememorysize(GMSHEAPNEW_theapmgr self, SYSTEM_int64 delta) {
   GMSHEAPNEW_tbigblockmgr_DOT_reducememorysize(self->GMSHEAPNEW_theapmgr_DOT_blockmgr, delta);
   self->GMSHEAPNEW_theapmgr_DOT_othermemory = self->GMSHEAPNEW_theapmgr_DOT_othermemory - delta;
}  /* reducememorysize */

Procedure GMSHEAPNEW_theapmgr_DOT_prvgetslotcnts(GMSHEAPNEW_theapmgr self, GMSHEAPNEW_theapslotnr slot, SYSTEM_int64* cntget, SYSTEM_int64* cntfree,
      SYSTEM_int64* cntavail) {
   {
      register GMSHEAPNEW_tslotrecord* _W2 = &self->GMSHEAPNEW_theapmgr_DOT_slots[slot - 1];
      *cntget = _W2->getcount;
      *cntfree = _W2->freecount;
      *cntavail = _W2->listcount;

   }
}  /* prvgetslotcnts */

Procedure GMSHEAPNEW_theapmgr_DOT_getslotcnts(GMSHEAPNEW_theapmgr self, GMSHEAPNEW_theapslotnr slot, SYSTEM_int64* cntget, SYSTEM_int64* cntfree,
      SYSTEM_int64* cntavail) {
   GMSHEAPNEW_theapmgr_DOT_prvgetslotcnts(self, slot, cntget, cntfree, cntavail);
}  /* getslotcnts */

Procedure
GMSHEAPNEW_theapmgr_DOT_getblockstats(GMSHEAPNEW_theapmgr self, SYSTEM_int64* cntwrkbuffs, SYSTEM_int64* cntactive, SYSTEM_int64* sizeothermemory,
      SYSTEM_int64* sizehighmark) {
   *cntwrkbuffs = self->GMSHEAPNEW_theapmgr_DOT_wrkbuffs->GMSOBJ_txlist_DOT_fcount;
   *cntactive = self->GMSHEAPNEW_theapmgr_DOT_active->GMSOBJ_txlist_DOT_fcount;
   *sizeothermemory = self->GMSHEAPNEW_theapmgr_DOT_othermemory;
   *sizehighmark = self->GMSHEAPNEW_theapmgr_DOT_highmark;
}  /* getblockstats */

Procedure GMSHEAPNEW_theapmgr_DOT_getotherstats(GMSHEAPNEW_theapmgr self, SYSTEM_boolean do64, SYSTEM_int64* cntget, SYSTEM_int64* cntfree,
      SYSTEM_int64* cntrealloc, SYSTEM_int64* sizerused) {
   if (do64) {
      *cntget = self->GMSHEAPNEW_theapmgr_DOT_otherget64;
      *cntfree = self->GMSHEAPNEW_theapmgr_DOT_otherfree64;
      *cntrealloc = self->GMSHEAPNEW_theapmgr_DOT_realloccnt64;
      *sizerused = self->GMSHEAPNEW_theapmgr_DOT_reallocused64;
   }
   else {
      *cntget = self->GMSHEAPNEW_theapmgr_DOT_otherget;
      *cntfree = self->GMSHEAPNEW_theapmgr_DOT_otherfree;
      *cntrealloc = self->GMSHEAPNEW_theapmgr_DOT_realloccnt;
      *sizerused = self->GMSHEAPNEW_theapmgr_DOT_reallocused;
   }
}  /* getotherstats */

Function(SYSTEM_int64) GMSHEAPNEW_theapmgr_DOT_getfreeslotspace(GMSHEAPNEW_theapmgr self) {
   SYSTEM_int64 result;
   GMSHEAPNEW_theapslotnr slot;
   SYSTEM_int64 cntget, cntfree, cntavail;
   SYSTEM_int64 sizefree;

   result = 0;
   for (slot = 1; slot <= 32; ++slot) {
      GMSHEAPNEW_theapmgr_DOT_prvgetslotcnts(self, slot, &cntget, &cntfree, &cntavail);
      sizefree = cntavail * slot * GMSHEAPNEW_heapgranularity;
      _P3inc1(result, sizefree);
   }
   return result;
}  /* getfreeslotspace */

Function(SYSTEM_ansichar *) GMSHEAPNEW_theapmgr_DOT_getname(SYSTEM_ansichar* result, SYSTEM_uint8 _len_ret, GMSHEAPNEW_theapmgr self) {
   _P3strcpy(result, _len_ret, *self->GMSHEAPNEW_theapmgr_DOT_spname);
   return result;
}  /* getname */

Function(SYSTEM_boolean) GMSHEAPNEW_theapmgr_DOT_setmemorylimit(GMSHEAPNEW_theapmgr self, SYSTEM_double limit) {
   SYSTEM_boolean result;

   self->GMSHEAPNEW_theapmgr_DOT_blockmgr->GMSHEAPNEW_tbigblockmgr_DOT_memorylimit = limit;
   result = limit >= self->GMSHEAPNEW_theapmgr_DOT_blockmgr->GMSHEAPNEW_tbigblockmgr_DOT_totalmemory;
   return result;
}  /* setmemorylimit */

Procedure GMSHEAPNEW_theapmgr_DOT_setmemoryreportproc(GMSHEAPNEW_theapmgr self, GMSHEAPNEW_tmemoryreportproc f) {
   self->GMSHEAPNEW_theapmgr_DOT_blockmgr->GMSHEAPNEW_tbigblockmgr_DOT_memoryreportproc = f;
}  /* setmemoryreportproc */

/* unit gmsheapnew */
void _Init_Module_gmsheapnew(void) {
} /* _Init_Module_gmsheapnew */

void _Final_Module_gmsheapnew(void) {
} /* _Final_Module_gmsheapnew */

