#if ! defined(_GDXSTATIC_H_)
#     define  _GDXSTATIC_H_

#include "gclgms.h"

#define gdxAcronymAdd          cgdxacronymadd
#define gdxAcronymCount         gdxacronymcount
#define gdxAcronymGetInfo      cgdxacronymgetinfo
#define gdxAcronymGetMapping    gdxacronymgetmapping
#define gdxAcronymIndex         gdxacronymindex
#define gdxAcronymName         cgdxacronymname
#define gdxAcronymNextNr        gdxacronymnextnr
#define gdxAcronymSetInfo      cgdxacronymsetinfo
#define gdxAcronymValue         gdxacronymvalue
#define gdxAddAlias            cgdxaddalias
#define gdxAddSetText          cgdxaddsettext
#define gdxAutoConvert          gdxautoconvert
#define gdxClose                gdxclose
#define gdxDataErrorCount       gdxdataerrorcount
#define gdxDataErrorRecord      gdxdataerrorrecord
#define gdxDataErrorRecordX     gdxdataerrorrecordx
#define gdxDataReadDone         gdxdatareaddone
#define gdxDataReadFilteredStart  gdxdatareadfilteredstart
#define gdxDataReadMap            gdxdatareadmap
#define gdxDataReadMapStart       gdxdatareadmapstart
#define gdxDataReadRaw            gdxdatareadraw
#define gdxDataReadRawFast        gdxdatareadrawfast
#define gdxDataReadRawFastFilt   cgdxdatareadrawfastfilt
#define gdxDataReadRawStart       gdxdatareadrawstart
#define gdxDataReadSlice         cgdxdatareadslice
#define gdxDataReadSliceStart     gdxdatareadslicestart
#define gdxDataReadStr           cgdxdatareadstr
#define gdxDataReadStrStart       gdxdatareadstrstart
#define gdxDataSliceUELS         cgdxdatasliceuels
#define gdxDataWriteDone          gdxdatawritedone
#define gdxDataWriteMap           gdxdatawritemap
#define gdxDataWriteMapStart     cgdxdatawritemapstart
#define gdxDataWriteRaw           gdxdatawriteraw
#define gdxDataWriteRawStart     cgdxdatawriterawstart
#define gdxDataWriteStr          cgdxdatawritestr
#define gdxDataWriteStrStart     cgdxdatawritestrstart
#define gdxGetDLLVersion         cgdxgetdllversion
#define gdxErrorCount             gdxerrorcount
#define gdxErrorStr              cgdxerrorstr
#define gdxFileInfo               gdxfileinfo
#define gdxFileVersion           cgdxfileversion
#define gdxFilterExists           gdxfilterexists
#define gdxFilterRegister         gdxfilterregister
#define gdxFilterRegisterDone     gdxfilterregisterdone
#define gdxFilterRegisterStart    gdxfilterregisterstart
#define gdxFindSymbol            cgdxfindsymbol
#define gdxGetElemText           cgdxgetelemtext
#define gdxGetLastError           gdxgetlasterror
#define gdxGetMemoryUsed          gdxgetmemoryused
#define gdxGetSpecialValues       gdxgetspecialvalues
#define gdxGetUEL                cgdxgetuel
#define gdxMapValue               gdxmapvalue
#define gdxOpenAppend            cgdxopenappend
#define gdxOpenRead              cgdxopenread
#define gdxOpenReadEx            cgdxopenreadex
#define gdxOpenWrite             cgdxopenwrite
#define gdxOpenWriteEx           cgdxopenwriteex
#define gdxResetSpecialValues     gdxresetspecialvalues
#define gdxSetHasText             gdxsethastext
#define gdxSetReadSpecialValues   gdxsetreadspecialvalues
#define gdxSetSpecialValues       gdxsetspecialvalues
#define gdxSetTextNodeNr          gdxsettextnodenr
#define gdxSetTraceLevel         cgdxsettracelevel
#define gdxSymbIndxMaxLength      gdxsymbindxmaxlength
#define gdxSymbMaxLength          gdxsymbmaxlength
#define gdxSymbolAddComment      cgdxsymboladdcomment
#define gdxSymbolGetComment      cgdxsymbolgetcomment
#define gdxSymbolGetDomain        gdxsymbolgetdomain
#define gdxSymbolGetDomainX      cgdxsymbolgetdomainx
#define gdxSymbolDim              gdxsymboldim
#define gdxSymbolInfo            cgdxsymbolinfo
#define gdxSymbolInfoX           cgdxsymbolinfox
#define gdxSymbolSetDomain       cgdxsymbolsetdomain
#define gdxSymbolSetDomainX      cgdxsymbolsetdomainx
#define gdxSystemInfo             gdxsysteminfo
#define gdxUELMaxLength           gdxuelmaxlength
#define gdxUELRegisterDone        gdxuelregisterdone
#define gdxUELRegisterMap        cgdxuelregistermap
#define gdxUELRegisterMapStart    gdxuelregistermapstart
#define gdxUELRegisterRaw        cgdxuelregisterraw
#define gdxUELRegisterRawStart    gdxuelregisterrawstart
#define gdxUELRegisterStr        cgdxuelregisterstr
#define gdxUELRegisterStrStart    gdxuelregisterstrstart
#define gdxUMFindUEL             cgdxumfinduel
#define gdxUMUelGet              cgdxumuelget
#define gdxUMUelInfo              gdxumuelinfo
#define gdxGetDomainElements      gdxgetdomainelements
#define gdxCurrentDim             gdxcurrentdim
#define gdxRenameUEL             cgdxrenameuel
#define gdxStoreDomainSets       gdxstoredomainsets
#define gdxStoreDomainSetsSet    gdxstoredomainsetsset

#if defined(_WIN32)
typedef __int64 INT64;
#elif defined(__LP64__) || defined(__axu__) || defined(_FCGLU_LP64_)
typedef signed long int INT64;
#else
typedef signed long long int INT64;
#endif

struct gdxRec;
typedef struct gdxRec *gdxHandle_t;

#if defined(__cplusplus)
extern "C" {
#endif

int gdxCreate    (gdxHandle_t *pgdx, char *msgBuf, int msgBufLen);
int gdxFree      (gdxHandle_t *pgdx);



/* function typedefs and pointer definitions */

typedef void (*TDataStoreProc_t) (const int Indx[], const double Vals[]);
typedef int (*TDataStoreFiltProc_t) (const int Indx[], const double Vals[], void *Uptr);
typedef void (*TDomainIndexProc_t) (int RawIndex, int MappedIndex, void *Uptr);

int gdxAcronymAdd (gdxHandle_t pgdx, const char *AName, const char *Txt, int AIndx);
int gdxAcronymCount (gdxHandle_t pgdx);
int gdxAcronymGetInfo (gdxHandle_t pgdx, int N, char *AName, char *Txt, int *AIndx);
int gdxAcronymGetMapping (gdxHandle_t pgdx, int N, int *orgIndx, int *newIndx, int *autoIndex);
int gdxAcronymIndex (gdxHandle_t pgdx, double V);
int gdxAcronymName (gdxHandle_t pgdx, double V, char *AName);
int gdxAcronymNextNr (gdxHandle_t pgdx, int NV);
int gdxAcronymSetInfo (gdxHandle_t pgdx, int N, const char *AName, const char *Txt, int AIndx);
double gdxAcronymValue (gdxHandle_t pgdx, int AIndx);
int gdxAddAlias (gdxHandle_t pgdx, const char *Id1, const char *Id2);
int gdxAddSetText (gdxHandle_t pgdx, const char *Txt, int *TxtNr);
int gdxAutoConvert (gdxHandle_t pgdx, int NV);
int gdxClose (gdxHandle_t pgdx);
int gdxDataErrorCount (gdxHandle_t pgdx);
int gdxDataErrorRecord (gdxHandle_t pgdx, int RecNr, int KeyInt[], double Values[]);
int gdxDataErrorRecordX (gdxHandle_t pgdx, int RecNr, int KeyInt[], double Values[]);
int gdxDataReadDone (gdxHandle_t pgdx);
int gdxDataReadFilteredStart (gdxHandle_t pgdx, int SyNr, const int FilterAction[], int *NrRecs);
int gdxDataReadMap (gdxHandle_t pgdx, int RecNr, int KeyInt[], double Values[], int *DimFrst);
int gdxDataReadMapStart (gdxHandle_t pgdx, int SyNr, int *NrRecs);
int gdxDataReadRaw (gdxHandle_t pgdx, int KeyInt[], double Values[], int *DimFrst);
int gdxDataReadRawFast (gdxHandle_t pgdx, int SyNr, TDataStoreProc_t DP, int *NrRecs);
int gdxDataReadRawFastFilt (gdxHandle_t pgdx, int SyNr, const char *UelFilterStr[], TDataStoreFiltProc_t DP);
int gdxDataReadRawStart (gdxHandle_t pgdx, int SyNr, int *NrRecs);
int gdxDataReadSlice (gdxHandle_t pgdx, const char *UelFilterStr[], int *Dimen, TDataStoreProc_t DP);
int gdxDataReadSliceStart (gdxHandle_t pgdx, int SyNr, int ElemCounts[]);
int gdxDataReadStr (gdxHandle_t pgdx, char *KeyStr[], double Values[], int *DimFrst);
int gdxDataReadStrStart (gdxHandle_t pgdx, int SyNr, int *NrRecs);
int gdxDataSliceUELS (gdxHandle_t pgdx, const int SliceKeyInt[], char *KeyStr[]);
int gdxDataWriteDone (gdxHandle_t pgdx);
int gdxDataWriteMap (gdxHandle_t pgdx, const int KeyInt[], const double Values[]);
int gdxDataWriteMapStart (gdxHandle_t pgdx, const char *SyId, const char *ExplTxt, int Dimen, int Typ, int UserInfo);
int gdxDataWriteRaw (gdxHandle_t pgdx, const int KeyInt[], const double Values[]);
int gdxDataWriteRawStart (gdxHandle_t pgdx, const char *SyId, const char *ExplTxt, int Dimen, int Typ, int UserInfo);
int gdxDataWriteStr (gdxHandle_t pgdx, const char *KeyStr[], const double Values[]);
int gdxDataWriteStrStart (gdxHandle_t pgdx, const char *SyId, const char *ExplTxt, int Dimen, int Typ, int UserInfo);
int gdxGetDLLVersion (gdxHandle_t pgdx, char *V);
int gdxErrorCount (gdxHandle_t pgdx);
int gdxErrorStr (gdxHandle_t pgdx, int ErrNr, char *ErrMsg);
int gdxFileInfo (gdxHandle_t pgdx, int *FileVer, int *ComprLev);
int gdxFileVersion (gdxHandle_t pgdx, char *FileStr, char *ProduceStr);
int gdxFilterExists (gdxHandle_t pgdx, int FilterNr);
int gdxFilterRegister (gdxHandle_t pgdx, int UelMap);
int gdxFilterRegisterDone (gdxHandle_t pgdx);
int gdxFilterRegisterStart (gdxHandle_t pgdx, int FilterNr);
int gdxFindSymbol (gdxHandle_t pgdx, const char *SyId, int *SyNr);
int gdxGetElemText (gdxHandle_t pgdx, int TxtNr, char *Txt, int *Node);
int gdxGetLastError (gdxHandle_t pgdx);
INT64 gdxGetMemoryUsed (gdxHandle_t pgdx);
int gdxGetSpecialValues (gdxHandle_t pgdx, double AVals[]);
int gdxGetUEL (gdxHandle_t pgdx, int UelNr, char *Uel);
int gdxMapValue (gdxHandle_t pgdx, double D, int *sv);
int gdxOpenAppend (gdxHandle_t pgdx, const char *FileName, const char *Producer, int *ErrNr);
int gdxOpenRead (gdxHandle_t pgdx, const char *FileName, int *ErrNr);
int gdxOpenReadEx (gdxHandle_t pgdx, const char *FileName, int ReadMode, int *ErrNr);
int gdxOpenWrite (gdxHandle_t pgdx, const char *FileName, const char *Producer, int *ErrNr);
int gdxOpenWriteEx (gdxHandle_t pgdx, const char *FileName, const char *Producer, int Compr, int *ErrNr);
int gdxResetSpecialValues (gdxHandle_t pgdx);
int gdxSetHasText (gdxHandle_t pgdx, int SyNr);
int gdxSetReadSpecialValues (gdxHandle_t pgdx, const double AVals[]);
int gdxSetSpecialValues (gdxHandle_t pgdx, const double AVals[]);
int gdxSetTextNodeNr (gdxHandle_t pgdx, int TxtNr, int Node);
int gdxSetTraceLevel (gdxHandle_t pgdx, int N, const char *s);
int gdxSymbIndxMaxLength (gdxHandle_t pgdx, int SyNr, int LengthInfo[]);
int gdxSymbMaxLength (gdxHandle_t pgdx);
int gdxSymbolAddComment (gdxHandle_t pgdx, int SyNr, const char *Txt);
int gdxSymbolGetComment (gdxHandle_t pgdx, int SyNr, int N, char *Txt);
int gdxSymbolGetDomain (gdxHandle_t pgdx, int SyNr, int DomainSyNrs[]);
int gdxSymbolGetDomainX (gdxHandle_t pgdx, int SyNr, char *DomainIDs[]);
int gdxSymbolDim (gdxHandle_t pgdx, int SyNr);
int gdxSymbolInfo (gdxHandle_t pgdx, int SyNr, char *SyId, int *Dimen, int *Typ);
int gdxSymbolInfoX (gdxHandle_t pgdx, int SyNr, int *RecCnt, int *UserInfo, char *ExplTxt);
int gdxSymbolSetDomain (gdxHandle_t pgdx, const char *DomainIDs[]);
int gdxSymbolSetDomainX (gdxHandle_t pgdx, int SyNr, const char *DomainIDs[]);
int gdxSystemInfo (gdxHandle_t pgdx, int *SyCnt, int *UelCnt);
int gdxUELMaxLength (gdxHandle_t pgdx);
int gdxUELRegisterDone (gdxHandle_t pgdx);
int gdxUELRegisterMap (gdxHandle_t pgdx, int UMap, const char *Uel);
int gdxUELRegisterMapStart (gdxHandle_t pgdx);
int gdxUELRegisterRaw (gdxHandle_t pgdx, const char *Uel);
int gdxUELRegisterRawStart (gdxHandle_t pgdx);
int gdxUELRegisterStr (gdxHandle_t pgdx, const char *Uel, int *UelNr);
int gdxUELRegisterStrStart (gdxHandle_t pgdx);
int gdxUMFindUEL (gdxHandle_t pgdx, const char *Uel, int *UelNr, int *UelMap);
int gdxUMUelGet (gdxHandle_t pgdx, int UelNr, char *Uel, int *UelMap);
int gdxUMUelInfo (gdxHandle_t pgdx, int *UelCnt, int *HighMap);
int gdxGetDomainElements (gdxHandle_t pgdx, int SyNr, int DimPos, int FilterNr, TDomainIndexProc_t DP, int *NrElem, void *Uptr);
int gdxCurrentDim (gdxHandle_t pgdx);
int gdxRenameUEL (gdxHandle_t pgdx, const char *OldName, const char *NewName);
int gdxStoreDomainSets (gdxHandle_t pgdx);
void gdxStoreDomainSetsSet (gdxHandle_t pgdx, const int x);


#if defined(__cplusplus)
}
#endif

#endif /* #if ! defined(_GDXSTATIC_H_) */
