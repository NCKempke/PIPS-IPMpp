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
#include "gmsglobx.h"

typedef SYSTEM_uint8 _sub_1GMSGLOBX;
typedef _P3STR_31 _arr_0GMSGLOBX[23];
static _arr_0GMSGLOBX GMSGLOBX_setconstantskeys = {{10,'M','o','d','e','l','T','y','p','e','s'}, {14,'G','a','m','s','P','a','r','a','m','e','t','e','r','s'}, {21,'G','a','m','s','P','a','r','a','m','e','t','e','r','S','y','n','o','n','y','m','s'}, {23,'G','a','m','s','P','a','r','a','m','e','t','e','r','S','y','n','o','n','y','m','M','a','p'}, {13,'D','o','l','l','a','r','O','p','t','i','o','n','s'}, {13,'G','a','m','s','F','u','n','c','t','i','o','n','s'}, {14,'S','y','s','t','e','m','S','u','f','f','i','x','e','s'}, {5,'E','m','p','t','y'}, {17,'P','r','e','d','e','f','i','n','e','d','S','y','m','b','o','l','s'}, {19,'G','U','S','S','M','o','d','e','l','A','t','t','r','i','b','u','t','e','s'}, {12,'S','e','t','C','o','n','s','t','a','n','t','s'}, {11,'S','o','l','v','e','r','N','a','m','e','s'}, {9,'P','l','a','t','f','o','r','m','s'}, {7,'V','e','n','d','o','r','s'}, {9,'T','L','L','i','c','e','n','s','e'}, {10,'C','o','m','p','o','n','e','n','t','s'}, {9,'C','l','i','p','C','o','d','e','s'}, {12,'G','a','m','s','L','i','c','e','n','s','e','s'}, {16,'G','a','m','s','L','i','c','e','n','s','e','T','y','p','e','s'}, {18,'C','o','m','p','o','n','e','n','t','S','o','l','v','e','r','M','a','p'}, {16,'C','l','i','p','C','o','m','p','o','n','e','n','t','M','a','p'}, {17,'S','o','l','v','e','r','P','l','a','t','f','o','r','m','M','a','p'}, {21,'S','o','l','v','e','r','T','y','p','e','P','l','a','t','f','o','r','m','M','a','p'}};

Function(SYSTEM_ansichar *) GMSGLOBX_setconstantstext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  switch (n) {
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\012ModelTypes"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\016GamsParameters"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\025GamsParameterSynonyms"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\027GamsParameterSynonymMap"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\015DollarOptions"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\015GamsFunctions"));
      break;
    case 7: 
      _P3strcpy(result,_len_ret,_P3str1("\016SystemSuffixes"));
      break;
    case 8: 
      _P3strcpy(result,_len_ret,_P3str1("\005Empty"));
      break;
    case 9: 
      _P3strcpy(result,_len_ret,_P3str1("\021PredefinedSymbols"));
      break;
    case 10: 
      _P3strcpy(result,_len_ret,_P3str1("\023GUSSModelAttributes"));
      break;
    case 11: 
      _P3strcpy(result,_len_ret,_P3str1("\014SetConstants"));
      break;
    case 12: 
      _P3strcpy(result,_len_ret,_P3str1("\013SolverNames"));
      break;
    case 13: 
      _P3strcpy(result,_len_ret,_P3str1("\011Platforms"));
      break;
    case 14: 
      _P3strcpy(result,_len_ret,_P3str1("\007Vendors"));
      break;
    case 15: 
      _P3strcpy(result,_len_ret,_P3str1("\011TLLicense"));
      break;
    case 16: 
      _P3strcpy(result,_len_ret,_P3str1("\012Components"));
      break;
    case 17: 
      _P3strcpy(result,_len_ret,_P3str1("\011ClipCodes"));
      break;
    case 18: 
      _P3strcpy(result,_len_ret,_P3str1("\014GamsLicenses"));
      break;
    case 19: 
      _P3strcpy(result,_len_ret,_P3str1("\020GamsLicenseTypes"));
      break;
    case 20: 
      _P3strcpy(result,_len_ret,_P3str1("\022ComponentSolverMap"));
      break;
    case 21: 
      _P3strcpy(result,_len_ret,_P3str1("\020ClipComponentMap"));
      break;
    case 22: 
      _P3strcpy(result,_len_ret,_P3str1("\021SolverPlatformMap"));
      break;
    case 23: 
      _P3strcpy(result,_len_ret,_P3str1("\025SolverTypePlatformMap"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030**** should never happen"));
  }
  return result;
}  /* setconstantstext */

Function(SYSTEM_integer ) GMSGLOBX_setconstantslookup(
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s;
  SYSTEM_integer result;

  _P3strcpy(s,255,_ftmp1);
  for (result = 1;result <= (SYSTEM_int32)
    GMSGLOBX_maxsetconstants;++result) {
    if (SYSUTILS_P3_sametext(GMSGLOBX_setconstantskeys[result - 1],
      s)) 
      return result;
  }
  result = 0;
  return result;
}  /* setconstantslookup */

Function(SYSTEM_ansichar *) GMSGLOBX_setconstantskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  if (n >= 1 && n <= 23) { 
    _P3strcpy(result,_len_ret,GMSGLOBX_setconstantskeys[n - 1]);
  } else 
    _P3strclr(result);
  return result;
}  /* setconstantskey */
typedef SYSTEM_uint8 _sub_3GMSGLOBX;
typedef _P3STR_15 _arr_2GMSGLOBX[153];
static _arr_2GMSGLOBX GMSGLOBX_solvernameskeys = {{2,'D','E'}, {2,'L','S'}, {2,'X','A'}, {8,'A','N','T','I','G','O','N','E'}, {5,'B','A','R','O','N'}, {6,'C','O','N','O','P','T'}, {5,'C','P','L','E','X'}, {11,'C','P','O','P','T','I','M','I','Z','E','R'}, {5,'D','E','C','I','S'}, {6,'D','I','C','O','P','T'}, {8,'A','L','P','H','A','E','C','P'}, {5,'M','P','S','G','E'}, {7,'G','L','O','M','I','Q','O'}, {6,'G','U','R','O','B','I'}, {5,'I','P','O','P','T'}, {6,'K','N','I','T','R','O'}, {3,'L','G','O'}, {5,'L','I','N','D','O'}, {11,'L','O','C','A','L','S','O','L','V','E','R'}, {5,'M','I','N','O','S'}, {5,'M','S','N','L','P'}, {8,'O','D','H','C','P','L','E','X'}, {4,'P','A','T','H'}, {3,'S','B','B'}, {4,'S','C','I','P'}, {5,'S','N','O','P','T'}, {6,'X','P','R','E','S','S'}, {4,'M','I','R','O'}, {6,'S','e','c','u','r','e'}, {5,'B','D','M','L','P'}, {5,'M','I','L','E','S'}, {5,'M','O','S','E','K'}, {4,'A','M','P','L'}, {6,'B','A','R','O','N','2'}, {6,'B','D','M','L','P','D'}, {5,'B','E','N','C','H'}, {6,'B','O','N','M','I','N'}, {3,'C','B','C'}, {7,'C','O','N','O','P','T','D'}, {7,'C','O','N','O','P','T','4'}, {7,'C','O','N','V','E','R','T'}, {8,'C','O','N','V','E','R','T','D'}, {12,'C','P','L','E','X','C','L','A','S','S','I','C'}, {6,'D','E','C','I','S','C'}, {6,'D','E','C','I','S','M'}, {7,'D','I','C','O','P','T','D'}, {8,'E','X','A','M','I','N','E','R'}, {9,'E','X','A','M','I','N','E','R','2'}, {7,'G','A','M','S','C','H','K'}, {10,'S','C','E','N','S','O','L','V','E','R'}, {4,'J','A','M','S'}, {7,'K','E','S','T','R','E','L'}, {9,'K','N','I','T','R','O','N','E','W'}, {4,'L','G','O','D'}, {11,'L','I','N','D','O','G','L','O','B','A','L'}, {5,'L','I','N','G','O'}, {13,'L','O','C','A','L','S','O','L','V','E','R','7','0'}, {7,'M','I','N','O','S','5','5'}, {9,'Q','U','A','D','M','I','N','O','S'}, {8,'M','P','E','C','D','U','M','P'}, {5,'N','L','P','E','C'}, {8,'O','s','i','C','p','l','e','x'}, {9,'O','s','i','G','u','r','o','b','i'}, {8,'O','s','i','M','o','s','e','k'}, {9,'O','s','i','X','p','r','e','s','s'}, {7,'P','A','T','H','N','L','P'}, {5,'P','A','T','H','C'}, {6,'P','A','T','H','V','I'}, {5,'P','Y','O','M','O'}, {6,'S','E','L','K','I','E'}, {4,'S','H','O','T'}, {6,'S','o','p','l','e','x'}, {3,'A','S','K'}, {7,'B','I','B','2','G','M','S'}, {7,'C','H','K','4','U','P','D'}, {8,'C','H','O','L','E','S','K','Y'}, {7,'C','S','V','2','G','D','X'}, {10,'E','I','G','E','N','V','A','L','U','E'}, {11,'E','I','G','E','N','V','E','C','T','O','R'}, {9,'E','N','D','E','C','R','Y','P','T'}, {12,'F','I','N','D','T','H','I','S','G','A','M','S'}, {7,'G','A','M','S','I','D','E'}, {10,'G','D','X','2','A','C','C','E','S','S'}, {10,'G','D','X','2','S','Q','L','I','T','E'}, {8,'G','D','X','2','V','E','D','A'}, {7,'G','D','X','2','X','L','S'}, {7,'G','D','X','C','O','P','Y'}, {7,'G','D','X','D','I','F','F'}, {7,'G','D','X','D','U','M','P'}, {8,'G','D','X','M','E','R','G','E'}, {6,'G','D','X','M','R','W'}, {7,'G','D','X','R','A','N','K'}, {9,'G','D','X','R','E','N','A','M','E'}, {6,'G','D','X','R','R','W'}, {8,'G','D','X','T','R','O','L','L'}, {9,'G','D','X','V','I','E','W','E','R'}, {6,'G','D','X','X','R','W'}, {9,'G','M','S','P','Y','T','H','O','N'}, {6,'G','M','S','Z','I','P'}, {7,'G','M','S','Z','L','I','B'}, {12,'H','A','R','U','T','I','L','I','T','I','E','S'}, {7,'H','E','X','D','U','M','P'}, {6,'I','N','V','E','R','T'}, {8,'M','C','F','I','L','T','E','R'}, {7,'M','D','B','2','G','M','S'}, {7,'M','S','G','R','W','I','N'}, {9,'M','O','D','E','L','2','T','E','X'}, {7,'M','P','S','2','G','M','S'}, {10,'M','S','A','P','P','A','V','A','I','L'}, {7,'S','C','E','N','R','E','D'}, {8,'S','C','E','N','R','E','D','2'}, {12,'S','H','E','L','L','E','X','E','C','U','T','E'}, {7,'S','Q','L','2','G','M','S'}, {6,'S','T','U','D','I','O'}, {5,'T','o','o','l','s'}, {7,'X','L','S','2','G','M','S'}, {7,'X','L','S','D','U','M','P'}, {7,'X','L','S','T','A','L','K'}, {4,'G','A','M','S'}, {9,'D','o','c','u','m','e','n','t','s'}, {4,'G','R','I','D'}, {6,'L','S','T','L','I','B'}, {12,'M','o','d','e','l','L','i','b','r','a','r','y'}, {11,'T','e','s','t','L','i','b','r','a','r','y'}, {11,'D','a','t','a','L','i','b','r','a','r','y'}, {10,'F','i','n','a','n','c','e','L','i','b'}, {10,'E','M','P','L','i','b','r','a','r','y'}, {10,'A','P','I','L','i','b','r','a','r','y'}, {12,'N','o','n','l','i','n','e','a','r','L','i','b'}, {12,'P','S','O','p','t','L','i','b','r','a','r','y'}, {6,'g','d','x','A','P','I'}, {6,'g','m','d','A','P','I'}, {6,'g','u','c','A','P','I'}, {7,'j','o','a','t','A','P','I'}, {9,'o','p','t','i','o','n','A','P','I'}, {7,'g','a','m','s','A','P','I'}, {9,'i','d','x','g','d','x','A','P','I'}, {7,'g','a','m','s','C','P','P'}, {10,'g','a','m','s','D','o','t','N','e','t'}, {10,'g','a','m','s','P','y','t','h','o','n'}, {8,'g','a','m','s','J','a','v','a'}, {6,'I','P','O','P','T','H'}, {7,'B','O','N','M','I','N','H'}, {7,'C','O','N','O','P','T','3'}, {10,'C','O','I','N','B','O','N','M','I','N'}, {7,'C','O','I','N','C','B','C'}, {9,'C','O','I','N','I','P','O','P','T'}, {8,'C','O','I','N','S','C','I','P'}, {6,'C','P','L','E','X','D'}, {6,'L','O','G','M','I','P'}, {6,'M','I','L','E','S','E'}, {6,'M','I','N','O','S','5'}, {9,'O','S','I','S','O','P','L','E','X'}};

Function(SYSTEM_ansichar *) GMSGLOBX_solvernamestext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  switch (n) {
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\002DE"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\002LS"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\002XA"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\010ANTIGONE"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\005BARON"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\006CONOPT"));
      break;
    case 7: 
      _P3strcpy(result,_len_ret,_P3str1("\005CPLEX"));
      break;
    case 8: 
      _P3strcpy(result,_len_ret,_P3str1("\013CPOPTIMIZER"));
      break;
    case 9: 
      _P3strcpy(result,_len_ret,_P3str1("\005DECIS"));
      break;
    case 10: 
      _P3strcpy(result,_len_ret,_P3str1("\006DICOPT"));
      break;
    case 11: 
      _P3strcpy(result,_len_ret,_P3str1("\010ALPHAECP"));
      break;
    case 12: 
      _P3strcpy(result,_len_ret,_P3str1("\005MPSGE"));
      break;
    case 13: 
      _P3strcpy(result,_len_ret,_P3str1("\007GLOMIQO"));
      break;
    case 14: 
      _P3strcpy(result,_len_ret,_P3str1("\006GUROBI"));
      break;
    case 15: 
      _P3strcpy(result,_len_ret,_P3str1("\005IPOPT"));
      break;
    case 16: 
      _P3strcpy(result,_len_ret,_P3str1("\006KNITRO"));
      break;
    case 17: 
      _P3strcpy(result,_len_ret,_P3str1("\003LGO"));
      break;
    case 18: 
      _P3strcpy(result,_len_ret,_P3str1("\005LINDO"));
      break;
    case 19: 
      _P3strcpy(result,_len_ret,_P3str1("\013LOCALSOLVER"));
      break;
    case 20: 
      _P3strcpy(result,_len_ret,_P3str1("\005MINOS"));
      break;
    case 21: 
      _P3strcpy(result,_len_ret,_P3str1("\005MSNLP"));
      break;
    case 22: 
      _P3strcpy(result,_len_ret,_P3str1("\010ODHCPLEX"));
      break;
    case 23: 
      _P3strcpy(result,_len_ret,_P3str1("\004PATH"));
      break;
    case 24: 
      _P3strcpy(result,_len_ret,_P3str1("\003SBB"));
      break;
    case 25: 
      _P3strcpy(result,_len_ret,_P3str1("\004SCIP"));
      break;
    case 26: 
      _P3strcpy(result,_len_ret,_P3str1("\005SNOPT"));
      break;
    case 27: 
      _P3strcpy(result,_len_ret,_P3str1("\006XPRESS"));
      break;
    case 28: 
      _P3strcpy(result,_len_ret,_P3str1("\004MIRO"));
      break;
    case 29: 
      _P3strcpy(result,_len_ret,_P3str1("\006Secure"));
      break;
    case 30: 
      _P3strcpy(result,_len_ret,_P3str1("\005BDMLP"));
      break;
    case 31: 
      _P3strcpy(result,_len_ret,_P3str1("\005MILES"));
      break;
    case 32: 
      _P3strcpy(result,_len_ret,_P3str1("\005MOSEK"));
      break;
    case 33: 
      _P3strcpy(result,_len_ret,_P3str1("\004AMPL"));
      break;
    case 34: 
      _P3strcpy(result,_len_ret,_P3str1("\006BARON2"));
      break;
    case 35: 
      _P3strcpy(result,_len_ret,_P3str1("\006BDMLPD"));
      break;
    case 36: 
      _P3strcpy(result,_len_ret,_P3str1("\005BENCH"));
      break;
    case 37: 
      _P3strcpy(result,_len_ret,_P3str1("\006BONMIN"));
      break;
    case 38: 
      _P3strcpy(result,_len_ret,_P3str1("\003CBC"));
      break;
    case 39: 
      _P3strcpy(result,_len_ret,_P3str1("\007CONOPTD"));
      break;
    case 40: 
      _P3strcpy(result,_len_ret,_P3str1("\007CONOPT4"));
      break;
    case 41: 
      _P3strcpy(result,_len_ret,_P3str1("\007CONVERT"));
      break;
    case 42: 
      _P3strcpy(result,_len_ret,_P3str1("\010CONVERTD"));
      break;
    case 43: 
      _P3strcpy(result,_len_ret,_P3str1("\014CPLEXCLASSIC"));
      break;
    case 44: 
      _P3strcpy(result,_len_ret,_P3str1("\006DECISC"));
      break;
    case 45: 
      _P3strcpy(result,_len_ret,_P3str1("\006DECISM"));
      break;
    case 46: 
      _P3strcpy(result,_len_ret,_P3str1("\007DICOPTD"));
      break;
    case 47: 
      _P3strcpy(result,_len_ret,_P3str1("\010EXAMINER"));
      break;
    case 48: 
      _P3strcpy(result,_len_ret,_P3str1("\011EXAMINER2"));
      break;
    case 49: 
      _P3strcpy(result,_len_ret,_P3str1("\007GAMSCHK"));
      break;
    case 50: 
      _P3strcpy(result,_len_ret,_P3str1("\012SCENSOLVER"));
      break;
    case 51: 
      _P3strcpy(result,_len_ret,_P3str1("\004JAMS"));
      break;
    case 52: 
      _P3strcpy(result,_len_ret,_P3str1("\007KESTREL"));
      break;
    case 53: 
      _P3strcpy(result,_len_ret,_P3str1("\011KNITRONEW"));
      break;
    case 54: 
      _P3strcpy(result,_len_ret,_P3str1("\004LGOD"));
      break;
    case 55: 
      _P3strcpy(result,_len_ret,_P3str1("\013LINDOGLOBAL"));
      break;
    case 56: 
      _P3strcpy(result,_len_ret,_P3str1("\005LINGO"));
      break;
    case 57: 
      _P3strcpy(result,_len_ret,_P3str1("\015LOCALSOLVER70"));
      break;
    case 58: 
      _P3strcpy(result,_len_ret,_P3str1("\007MINOS55"));
      break;
    case 59: 
      _P3strcpy(result,_len_ret,_P3str1("\011QUADMINOS"));
      break;
    case 60: 
      _P3strcpy(result,_len_ret,_P3str1("\010MPECDUMP"));
      break;
    case 61: 
      _P3strcpy(result,_len_ret,_P3str1("\005NLPEC"));
      break;
    case 62: 
      _P3strcpy(result,_len_ret,_P3str1("\010OsiCplex"));
      break;
    case 63: 
      _P3strcpy(result,_len_ret,_P3str1("\011OsiGurobi"));
      break;
    case 64: 
      _P3strcpy(result,_len_ret,_P3str1("\010OsiMosek"));
      break;
    case 65: 
      _P3strcpy(result,_len_ret,_P3str1("\011OsiXpress"));
      break;
    case 66: 
      _P3strcpy(result,_len_ret,_P3str1("\007PATHNLP"));
      break;
    case 67: 
      _P3strcpy(result,_len_ret,_P3str1("\005PATHC"));
      break;
    case 68: 
      _P3strcpy(result,_len_ret,_P3str1("\006PATHVI"));
      break;
    case 69: 
      _P3strcpy(result,_len_ret,_P3str1("\005PYOMO"));
      break;
    case 70: 
      _P3strcpy(result,_len_ret,_P3str1("\006SELKIE"));
      break;
    case 71: 
      _P3strcpy(result,_len_ret,_P3str1("\004SHOT"));
      break;
    case 72: 
      _P3strcpy(result,_len_ret,_P3str1("\006Soplex"));
      break;
    case 73: 
      _P3strcpy(result,_len_ret,_P3str1("\003ASK"));
      break;
    case 74: 
      _P3strcpy(result,_len_ret,_P3str1("\007BIB2GMS"));
      break;
    case 75: 
      _P3strcpy(result,_len_ret,_P3str1("\007CHK4UPD"));
      break;
    case 76: 
      _P3strcpy(result,_len_ret,_P3str1("\010CHOLESKY"));
      break;
    case 77: 
      _P3strcpy(result,_len_ret,_P3str1("\007CSV2GDX"));
      break;
    case 78: 
      _P3strcpy(result,_len_ret,_P3str1("\012EIGENVALUE"));
      break;
    case 79: 
      _P3strcpy(result,_len_ret,_P3str1("\013EIGENVECTOR"));
      break;
    case 80: 
      _P3strcpy(result,_len_ret,_P3str1("\011ENDECRYPT"));
      break;
    case 81: 
      _P3strcpy(result,_len_ret,_P3str1("\014FINDTHISGAMS"));
      break;
    case 82: 
      _P3strcpy(result,_len_ret,_P3str1("\007GAMSIDE"));
      break;
    case 83: 
      _P3strcpy(result,_len_ret,_P3str1("\012GDX2ACCESS"));
      break;
    case 84: 
      _P3strcpy(result,_len_ret,_P3str1("\012GDX2SQLITE"));
      break;
    case 85: 
      _P3strcpy(result,_len_ret,_P3str1("\010GDX2VEDA"));
      break;
    case 86: 
      _P3strcpy(result,_len_ret,_P3str1("\007GDX2XLS"));
      break;
    case 87: 
      _P3strcpy(result,_len_ret,_P3str1("\007GDXCOPY"));
      break;
    case 88: 
      _P3strcpy(result,_len_ret,_P3str1("\007GDXDIFF"));
      break;
    case 89: 
      _P3strcpy(result,_len_ret,_P3str1("\007GDXDUMP"));
      break;
    case 90: 
      _P3strcpy(result,_len_ret,_P3str1("\010GDXMERGE"));
      break;
    case 91: 
      _P3strcpy(result,_len_ret,_P3str1("\006GDXMRW"));
      break;
    case 92: 
      _P3strcpy(result,_len_ret,_P3str1("\007GDXRANK"));
      break;
    case 93: 
      _P3strcpy(result,_len_ret,_P3str1("\011GDXRENAME"));
      break;
    case 94: 
      _P3strcpy(result,_len_ret,_P3str1("\006GDXRRW"));
      break;
    case 95: 
      _P3strcpy(result,_len_ret,_P3str1("\010GDXTROLL"));
      break;
    case 96: 
      _P3strcpy(result,_len_ret,_P3str1("\011GDXVIEWER"));
      break;
    case 97: 
      _P3strcpy(result,_len_ret,_P3str1("\006GDXXRW"));
      break;
    case 98: 
      _P3strcpy(result,_len_ret,_P3str1("\011GMSPYTHON"));
      break;
    case 99: 
      _P3strcpy(result,_len_ret,_P3str1("\006GMSZIP"));
      break;
    case 100: 
      _P3strcpy(result,_len_ret,_P3str1("\007GMSZLIB"));
      break;
    case 101: 
      _P3strcpy(result,_len_ret,_P3str1("\014HARUTILITIES"));
      break;
    case 102: 
      _P3strcpy(result,_len_ret,_P3str1("\007HEXDUMP"));
      break;
    case 103: 
      _P3strcpy(result,_len_ret,_P3str1("\006INVERT"));
      break;
    case 104: 
      _P3strcpy(result,_len_ret,_P3str1("\010MCFILTER"));
      break;
    case 105: 
      _P3strcpy(result,_len_ret,_P3str1("\007MDB2GMS"));
      break;
    case 106: 
      _P3strcpy(result,_len_ret,_P3str1("\007MSGRWIN"));
      break;
    case 107: 
      _P3strcpy(result,_len_ret,_P3str1("\011MODEL2TEX"));
      break;
    case 108: 
      _P3strcpy(result,_len_ret,_P3str1("\007MPS2GMS"));
      break;
    case 109: 
      _P3strcpy(result,_len_ret,_P3str1("\012MSAPPAVAIL"));
      break;
    case 110: 
      _P3strcpy(result,_len_ret,_P3str1("\007SCENRED"));
      break;
    case 111: 
      _P3strcpy(result,_len_ret,_P3str1("\010SCENRED2"));
      break;
    case 112: 
      _P3strcpy(result,_len_ret,_P3str1("\014SHELLEXECUTE"));
      break;
    case 113: 
      _P3strcpy(result,_len_ret,_P3str1("\007SQL2GMS"));
      break;
    case 114: 
      _P3strcpy(result,_len_ret,_P3str1("\006STUDIO"));
      break;
    case 115: 
      _P3strcpy(result,_len_ret,_P3str1("\005Tools"));
      break;
    case 116: 
      _P3strcpy(result,_len_ret,_P3str1("\007XLS2GMS"));
      break;
    case 117: 
      _P3strcpy(result,_len_ret,_P3str1("\007XLSDUMP"));
      break;
    case 118: 
      _P3strcpy(result,_len_ret,_P3str1("\007XLSTALK"));
      break;
    case 119: 
      _P3strcpy(result,_len_ret,_P3str1("\004GAMS"));
      break;
    case 120: 
      _P3strcpy(result,_len_ret,_P3str1("\011Documents"));
      break;
    case 121: 
      _P3strcpy(result,_len_ret,_P3str1("\004GRID"));
      break;
    case 122: 
      _P3strcpy(result,_len_ret,_P3str1("\006LSTLIB"));
      break;
    case 123: 
      _P3strcpy(result,_len_ret,_P3str1("\014ModelLibrary"));
      break;
    case 124: 
      _P3strcpy(result,_len_ret,_P3str1("\013TestLibrary"));
      break;
    case 125: 
      _P3strcpy(result,_len_ret,_P3str1("\013DataLibrary"));
      break;
    case 126: 
      _P3strcpy(result,_len_ret,_P3str1("\012FinanceLib"));
      break;
    case 127: 
      _P3strcpy(result,_len_ret,_P3str1("\012EMPLibrary"));
      break;
    case 128: 
      _P3strcpy(result,_len_ret,_P3str1("\012APILibrary"));
      break;
    case 129: 
      _P3strcpy(result,_len_ret,_P3str1("\014NonlinearLib"));
      break;
    case 130: 
      _P3strcpy(result,_len_ret,_P3str1("\014PSOptLibrary"));
      break;
    case 131: 
      _P3strcpy(result,_len_ret,_P3str1("\006gdxAPI"));
      break;
    case 132: 
      _P3strcpy(result,_len_ret,_P3str1("\006gmdAPI"));
      break;
    case 133: 
      _P3strcpy(result,_len_ret,_P3str1("\006gucAPI"));
      break;
    case 134: 
      _P3strcpy(result,_len_ret,_P3str1("\007joatAPI"));
      break;
    case 135: 
      _P3strcpy(result,_len_ret,_P3str1("\011optionAPI"));
      break;
    case 136: 
      _P3strcpy(result,_len_ret,_P3str1("\007gamsAPI"));
      break;
    case 137: 
      _P3strcpy(result,_len_ret,_P3str1("\011idxgdxAPI"));
      break;
    case 138: 
      _P3strcpy(result,_len_ret,_P3str1("\007gamsCPP"));
      break;
    case 139: 
      _P3strcpy(result,_len_ret,_P3str1("\012gamsDotNet"));
      break;
    case 140: 
      _P3strcpy(result,_len_ret,_P3str1("\012gamsPython"));
      break;
    case 141: 
      _P3strcpy(result,_len_ret,_P3str1("\010gamsJava"));
      break;
    case 142: 
      _P3strcpy(result,_len_ret,_P3str1("\006IPOPTH"));
      break;
    case 143: 
      _P3strcpy(result,_len_ret,_P3str1("\007BONMINH"));
      break;
    case 144: 
      _P3strcpy(result,_len_ret,_P3str1("\007CONOPT3"));
      break;
    case 145: 
      _P3strcpy(result,_len_ret,_P3str1("\012COINBONMIN"));
      break;
    case 146: 
      _P3strcpy(result,_len_ret,_P3str1("\007COINCBC"));
      break;
    case 147: 
      _P3strcpy(result,_len_ret,_P3str1("\011COINIPOPT"));
      break;
    case 148: 
      _P3strcpy(result,_len_ret,_P3str1("\010COINSCIP"));
      break;
    case 149: 
      _P3strcpy(result,_len_ret,_P3str1("\006CPLEXD"));
      break;
    case 150: 
      _P3strcpy(result,_len_ret,_P3str1("\006LOGMIP"));
      break;
    case 151: 
      _P3strcpy(result,_len_ret,_P3str1("\006MILESE"));
      break;
    case 152: 
      _P3strcpy(result,_len_ret,_P3str1("\006MINOS5"));
      break;
    case 153: 
      _P3strcpy(result,_len_ret,_P3str1("\011OSISOPLEX"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030**** should never happen"));
  }
  return result;
}  /* solvernamestext */

Function(SYSTEM_integer ) GMSGLOBX_solvernameslookup(
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s;
  SYSTEM_integer result;

  _P3strcpy(s,255,_ftmp1);
  for (result = 1;result <= (SYSTEM_int32)GMSGLOBX_maxsolvernames;++
    result) {
    if (SYSUTILS_P3_sametext(GMSGLOBX_solvernameskeys[result - 1],
      s)) 
      return result;
  }
  result = 0;
  return result;
}  /* solvernameslookup */

Function(SYSTEM_ansichar *) GMSGLOBX_solvernameskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  if (n >= 1 && n <= 153) { 
    _P3strcpy(result,_len_ret,GMSGLOBX_solvernameskeys[n - 1]);
  } else 
    _P3strclr(result);
  return result;
}  /* solvernameskey */
typedef SYSTEM_uint8 _sub_5GMSGLOBX;
typedef _P3STR_3 _arr_4GMSGLOBX[14];
static _arr_4GMSGLOBX GMSGLOBX_platformskeys = {{3,'W','E','X'}, {3,'L','E','X'}, {3,'D','E','X'}, {3,'S','O','X'}, {3,'A','I','X'}, {3,'W','I','N'}, {3,'S','I','S'}, {3,'B','G','P'}, {3,'L','N','X'}, {3,'S','O','L'}, {3,'D','A','R'}, {3,'D','I','I'}, {3,'G','E','N'}, {3,'A','L','L'}};

Function(SYSTEM_ansichar *) GMSGLOBX_platformstext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  switch (n) {
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\024x86 64bit MS Windows"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\017x86 64bit Linux"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\022x86 64bit Mac OS X"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\023Sparc 64bit SOLARIS"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\023IBM Power 64bit AIX"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\041x86 32bit MS Windows 32 and 64bit"));
      break;
    case 7: 
      _P3strcpy(result,_len_ret,_P3str1("\021x86 64bit SOLARIS"));
      break;
    case 8: 
      _P3strcpy(result,_len_ret,_P3str1("\025IBM Power 64bit Linux"));
      break;
    case 9: 
      _P3strcpy(result,_len_ret,_P3str1("\026x86 Linux 32 and 64bit"));
      break;
    case 10: 
      _P3strcpy(result,_len_ret,_P3str1("\021Sun Sparc SOLARIS"));
      break;
    case 11: 
      _P3strcpy(result,_len_ret,_P3str1("\022Mac PowerPC Darwin"));
      break;
    case 12: 
      _P3strcpy(result,_len_ret,_P3str1("\022Mac Intel32 Darwin"));
      break;
    case 13: 
      _P3strcpy(result,_len_ret,_P3str1("\031Generic popular platforms"));
      break;
    case 14: 
      _P3strcpy(result,_len_ret,_P3str1("\015All platforms"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030**** should never happen"));
  }
  return result;
}  /* platformstext */

Function(SYSTEM_ansichar *) GMSGLOBX_platformstext2(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  switch (n) {
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\024x86 64bit MS Windows"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\017x86 64bit Linux"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\022x86 64bit Mac OS X"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\023Sparc 64bit SOLARIS"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\023IBM Power 64bit AIX"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\024x86 32bit MS Windows"));
      break;
    case 7: 
      _P3strcpy(result,_len_ret,_P3str1("\021x86 64bit SOLARIS"));
      break;
    case 8: 
      _P3strcpy(result,_len_ret,_P3str1("\025IBM Power 64bit Linux"));
      break;
    case 9: 
      _P3strcpy(result,_len_ret,_P3str1("\011x86 Linux"));
      break;
    case 10: 
      _P3strcpy(result,_len_ret,_P3str1("\021Sun Sparc SOLARIS"));
      break;
    case 11: 
      _P3strcpy(result,_len_ret,_P3str1("\022Mac PowerPC Darwin"));
      break;
    case 12: 
      _P3strcpy(result,_len_ret,_P3str1("\022Mac Intel32 Darwin"));
      break;
    case 13: 
      _P3strcpy(result,_len_ret,_P3str1("\031Generic popular platforms"));
      break;
    case 14: 
      _P3strcpy(result,_len_ret,_P3str1("\015All platforms"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030**** should never happen"));
  }
  return result;
}  /* platformstext2 */

Function(SYSTEM_integer ) GMSGLOBX_platformslookup(
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s;
  SYSTEM_integer result;

  _P3strcpy(s,255,_ftmp1);
  for (result = 1;result <= (SYSTEM_int32)GMSGLOBX_maxplatforms;++
    result) {
    if (SYSUTILS_P3_sametext(GMSGLOBX_platformskeys[result - 1],s)) 
      return result;
  }
  result = 0;
  return result;
}  /* platformslookup */

Function(SYSTEM_ansichar *) GMSGLOBX_platformskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  if (n >= 1 && n <= 14) { 
    _P3strcpy(result,_len_ret,GMSGLOBX_platformskeys[n - 1]);
  } else 
    _P3strclr(result);
  return result;
}  /* platformskey */
typedef SYSTEM_uint8 _sub_7GMSGLOBX;
typedef _P3STR_3 _arr_6GMSGLOBX[21];
static _arr_6GMSGLOBX GMSGLOBX_vendorskeys = {{1,'G'}, {1,'B'}, {1,'A'}, {1,'C'}, {1,'D'}, {1,'H'}, {1,'K'}, {1,'W'}, {1,'L'}, {1,'S'}, {1,'E'}, {1,'F'}, {1,'I'}, {1,'J'}, {1,'N'}, {1,'O'}, {1,'P'}, {1,'Q'}, {1,'R'}, {1,'#'}, {1,'x'}};

Function(SYSTEM_ansichar *) GMSGLOBX_vendorstext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  switch (n) {
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\034GAMS Development Corporation"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\023BOYOUN TMS CO., LTD"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\041ARKI Consulting & Development A/S"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\030Scientific Formosa China"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\031Delta Computer Consulting"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\027Scientific Formosa Inc."));
      break;
    case 7: 
      _P3strcpy(result,_len_ret,_P3str1("\042kfs, scientific management systems"));
      break;
    case 8: 
      _P3strcpy(result,_len_ret,_P3str1("\022Decision Ware Inc."));
      break;
    case 9: 
      _P3strcpy(result,_len_ret,_P3str1("\033AddLink Software Cientifico"));
      break;
    case 10: 
      _P3strcpy(result,_len_ret,_P3str1("\022GAMS Software GmbH"));
      break;
    case 11: 
      _P3strcpy(result,_len_ret,_P3str1("\045Beijing Tianyan Rongzhi Ltd./Turntech"));
      break;
    case 12: 
      _P3strcpy(result,_len_ret,_P3str1("\012Focus Corp"));
      break;
    case 13: 
      _P3strcpy(result,_len_ret,_P3str1("\014Hulinks Inc."));
      break;
    case 14: 
      _P3strcpy(result,_len_ret,_P3str1("\011JUCA Inc."));
      break;
    case 15: 
      _P3strcpy(result,_len_ret,_P3str1("\047MultiON Consulting, S.A. de C.V.MultiOn"));
      break;
    case 16: 
      _P3strcpy(result,_len_ret,_P3str1("\020Pitotech Co.,Ltd"));
      break;
    case 17: 
      _P3strcpy(result,_len_ret,_P3str1("\022SOFTWARE Shop Inc."));
      break;
    case 18: 
      _P3strcpy(result,_len_ret,_P3str1("\022Tegara Corporation"));
      break;
    case 19: 
      _P3strcpy(result,_len_ret,_P3str1("\045Amsterdam Optimization Modeling Group"));
      break;
    case 20: 
      _P3strcpy(result,_len_ret,_P3str1("\014Other Vendor"));
      break;
    case 21: 
      _P3strcpy(result,_len_ret,_P3str1("\016Unknown Vendor"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030**** should never happen"));
  }
  return result;
}  /* vendorstext */

Function(SYSTEM_integer ) GMSGLOBX_vendorslookup(
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s;
  SYSTEM_integer result;

  _P3strcpy(s,255,_ftmp1);
  for (result = 1;result <= (SYSTEM_int32)GMSGLOBX_maxvendors;++
    result) {
    if (SYSUTILS_P3_sametext(GMSGLOBX_vendorskeys[result - 1],s)) 
      return result;
  }
  result = 0;
  return result;
}  /* vendorslookup */

Function(SYSTEM_ansichar *) GMSGLOBX_vendorskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  if (n >= 1 && n <= 21) { 
    _P3strcpy(result,_len_ret,GMSGLOBX_vendorskeys[n - 1]);
  } else 
    _P3strclr(result);
  return result;
}  /* vendorskey */
typedef SYSTEM_uint8 _sub_9GMSGLOBX;
typedef _P3STR_15 _arr_8GMSGLOBX[6];
static _arr_8GMSGLOBX GMSGLOBX_tllicensekeys = {{10,'E','V','A','L','U','A','T','I','O','N'}, {6,'C','O','U','R','S','E'}, {5,'L','E','A','S','E'}, {6,'A','N','N','U','A','L'}, {15,'T','A','K','E','G','A','M','S','W','I','T','H','Y','O','U'}, {11,'A','P','P','L','I','C','A','T','I','O','N'}};

Function(SYSTEM_ansichar *) GMSGLOBX_tllicensetext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  switch (n) {
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\070Evaluation license: Not for commercial or production use"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\100Course license for use within the course and related course work"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\041GAMS month to month lease license"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\031GAMS annual lease license"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\100Take-GAMS-With-You license: Not for commercial or production use"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\036GAMS Application lease license"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030**** should never happen"));
  }
  return result;
}  /* tllicensetext */

Function(SYSTEM_integer ) GMSGLOBX_tllicenselookup(
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s;
  SYSTEM_integer result;

  _P3strcpy(s,255,_ftmp1);
  for (result = 1;result <= (SYSTEM_int32)GMSGLOBX_maxtllicense;++
    result) {
    if (SYSUTILS_P3_sametext(GMSGLOBX_tllicensekeys[result - 1],s)) 
      return result;
  }
  result = 0;
  return result;
}  /* tllicenselookup */

Function(SYSTEM_ansichar *) GMSGLOBX_tllicensekey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  if (n >= 1 && n <= 6) { 
    _P3strcpy(result,_len_ret,GMSGLOBX_tllicensekeys[n - 1]);
  } else 
    _P3strclr(result);
  return result;
}  /* tllicensekey */
typedef SYSTEM_uint8 _sub_11GMSGLOBX;
typedef _P3STR_15 _arr_10GMSGLOBX[46];
static _arr_10GMSGLOBX GMSGLOBX_componentskeys = {{2,'X','A'}, {13,'M','I','R','O','C','o','n','n','e','c','t','o','r'}, {14,'S','e','c','u','r','e','W','o','r','k','f','i','l','e'}, {8,'A','N','T','I','G','O','N','E'}, {5,'B','A','R','O','N'}, {9,'C','P','L','E','X','L','i','n','k'}, {6,'C','O','N','O','P','T'}, {5,'C','P','L','E','X'}, {9,'C','P','L','E','X','D','I','S','T'}, {11,'C','P','O','P','T','I','M','I','Z','E','R'}, {5,'D','E','C','I','S'}, {6,'D','I','C','O','P','T'}, {8,'A','L','P','H','A','E','C','P'}, {5,'M','P','S','G','E'}, {6,'0','0','B','a','s','e'}, {7,'G','L','O','M','I','Q','O'}, {10,'G','U','R','O','B','I','L','i','n','k'}, {6,'G','U','R','O','B','I'}, {10,'G','U','R','O','B','I','D','i','s','t'}, {5,'I','P','O','P','T'}, {6,'K','N','I','T','R','O'}, {8,'L','I','N','D','O','A','P','I'}, {3,'L','G','O'}, {5,'L','I','N','D','O'}, {11,'L','O','C','A','L','S','O','L','V','E','R'}, {14,'L','O','C','A','L','S','O','L','V','E','R','L','N','K'}, {5,'M','I','N','O','S'}, {9,'M','O','S','E','K','B','a','s','e'}, {9,'M','O','S','E','K','L','i','n','k'}, {8,'M','O','S','E','K','M','I','P'}, {5,'M','S','N','L','P'}, {11,'O','S','L','V','e','r','s','i','o','n','2'}, {8,'C','P','L','E','X','O','S','I'}, {8,'O','D','H','C','P','L','E','X'}, {7,'O','S','L','L','i','n','k'}, {5,'O','Q','N','L','P'}, {11,'O','S','L','V','e','r','s','i','o','n','1'}, {4,'P','A','T','H'}, {3,'S','B','B'}, {4,'S','C','I','P'}, {5,'O','S','L','S','E'}, {5,'S','N','O','P','T'}, {10,'X','P','R','E','S','S','L','i','n','k'}, {6,'X','P','R','E','S','S'}, {9,'X','P','R','E','S','S','S','L','P'}, {9,'X','P','R','E','S','S','A','L','L'}};

Function(SYSTEM_ansichar *) GMSGLOBX_componentstext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  switch (n) {
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\041GAMS/XA           Sunset Software"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\045GAMS/MIRO         GAMS MIRO Connector"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\045GAMS/Secure       GAMS Secure Version"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\046GAMS/ANTIGONE     Princeton University"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\050GAMS/BARON        University of Illinois"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\040GAMS/CPLEX Link   IBM-ILOG/CPLEX"));
      break;
    case 7: 
      _P3strcpy(result,_len_ret,_P3str1("\041GAMS/CONOPT       ARKI Consulting"));
      break;
    case 8: 
      _P3strcpy(result,_len_ret,_P3str1("\040GAMS/CPLEX        IBM-ILOG/CPLEX"));
      break;
    case 9: 
      _P3strcpy(result,_len_ret,_P3str1("\040GAMS/CPLEXDIST    IBM-ILOG/CPLEX"));
      break;
    case 10: 
      _P3strcpy(result,_len_ret,_P3str1("\040GAMS/CPOPTIMIZER  IBM-ILOG/CPLEX"));
      break;
    case 11: 
      _P3strcpy(result,_len_ret,_P3str1("\041GAMS/DECIS        G Infanger, Inc"));
      break;
    case 12: 
      _P3strcpy(result,_len_ret,_P3str1("\033GAMS/DICOPT       CAPD, CMU"));
      break;
    case 13: 
      _P3strcpy(result,_len_ret,_P3str1("\051GAMS/ALPHAECP     Abo University, Finland"));
      break;
    case 14: 
      _P3strcpy(result,_len_ret,_P3str1("\043GAMS/MPSGE        Thomas Rutherford"));
      break;
    case 15: 
      _P3strcpy(result,_len_ret,_P3str1("\047GAMS/Base Module  GAMS Development Corp"));
      break;
    case 16: 
      _P3strcpy(result,_len_ret,_P3str1("\046GAMS/GLOMIQO      Princeton University"));
      break;
    case 17: 
      _P3strcpy(result,_len_ret,_P3str1("\045GAMS/GUROBI Link  Gurobi Optimization"));
      break;
    case 18: 
      _P3strcpy(result,_len_ret,_P3str1("\045GAMS/GUROBI       Gurobi Optimization"));
      break;
    case 19: 
      _P3strcpy(result,_len_ret,_P3str1("\045GAMS/GUROBI Dist  Gurobi Optimization"));
      break;
    case 20: 
      _P3strcpy(result,_len_ret,_P3str1("\044GAMS/IPOPT        COIN-OR Foundation"));
      break;
    case 21: 
      _P3strcpy(result,_len_ret,_P3str1("\031GAMS/KNITRO       Artelys"));
      break;
    case 22: 
      _P3strcpy(result,_len_ret,_P3str1("\045GAMS/LINDOAPI     Lindo Systems, Inc."));
      break;
    case 23: 
      _P3strcpy(result,_len_ret,_P3str1("\054GAMS/LGO          Pinter Consulting Services"));
      break;
    case 24: 
      _P3strcpy(result,_len_ret,_P3str1("\045GAMS/LINDO        Lindo Systems, Inc."));
      break;
    case 25: 
      _P3strcpy(result,_len_ret,_P3str1("\035GAMS/LOCALSOLVER  LocalSolver"));
      break;
    case 26: 
      _P3strcpy(result,_len_ret,_P3str1("\035GAMS/LOCALSLNK    LocalSolver"));
      break;
    case 27: 
      _P3strcpy(result,_len_ret,_P3str1("\045GAMS/MINOS        Stanford University"));
      break;
    case 28: 
      _P3strcpy(result,_len_ret,_P3str1("\034GAMS/MOSEK Base   MOSEK Base"));
      break;
    case 29: 
      _P3strcpy(result,_len_ret,_P3str1("\034GAMS/MOSEK Link   MOSEK Link"));
      break;
    case 30: 
      _P3strcpy(result,_len_ret,_P3str1("\027GAMS/MOSEK        MOSEK"));
      break;
    case 31: 
      _P3strcpy(result,_len_ret,_P3str1("\047GAMS/MSNLP        Optimal Methods, Inc."));
      break;
    case 32: 
      _P3strcpy(result,_len_ret,_P3str1("\025GAMS/OSL V2       IBM"));
      break;
    case 33: 
      _P3strcpy(result,_len_ret,_P3str1("\044GAMS/CPLEXOSI     COIN-OR Foundation"));
      break;
    case 34: 
      _P3strcpy(result,_len_ret,_P3str1("\052GAMS/ODHCPLEX     Optimization Direct Inc."));
      break;
    case 35: 
      _P3strcpy(result,_len_ret,_P3str1("\025GAMS/OSL Link     IBM"));
      break;
    case 36: 
      _P3strcpy(result,_len_ret,_P3str1("\047GAMS/OQNLP        Optimal Methods, Inc."));
      break;
    case 37: 
      _P3strcpy(result,_len_ret,_P3str1("\025GAMS/OSL          IBM"));
      break;
    case 38: 
      _P3strcpy(result,_len_ret,_P3str1("\056GAMS/PATH         University Wisconsin-Madison"));
      break;
    case 39: 
      _P3strcpy(result,_len_ret,_P3str1("\041GAMS/SBB          ARKI Consulting"));
      break;
    case 40: 
      _P3strcpy(result,_len_ret,_P3str1("\034GAMS/SCIP         ZIB Berlin"));
      break;
    case 41: 
      _P3strcpy(result,_len_ret,_P3str1("\057GAMS/OSL/SE       IBM/OSL Stochastic Extensions"));
      break;
    case 42: 
      _P3strcpy(result,_len_ret,_P3str1("\050GAMS/SNOPT        Stanford & UC-SanDiego"));
      break;
    case 43: 
      _P3strcpy(result,_len_ret,_P3str1("\026GAMS/XPRESS Link  FICO"));
      break;
    case 44: 
      _P3strcpy(result,_len_ret,_P3str1("\026GAMS/XPRESS       FICO"));
      break;
    case 45: 
      _P3strcpy(result,_len_ret,_P3str1("\026GAMS/XPRESSSLP    FICO"));
      break;
    case 46: 
      _P3strcpy(result,_len_ret,_P3str1("\026GAMS/XPRESSALL    FICO"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030**** should never happen"));
  }
  return result;
}  /* componentstext */

Function(SYSTEM_integer ) GMSGLOBX_componentslookup(
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s;
  SYSTEM_integer result;

  _P3strcpy(s,255,_ftmp1);
  for (result = 1;result <= (SYSTEM_int32)GMSGLOBX_maxcomponents;++
    result) {
    if (SYSUTILS_P3_sametext(GMSGLOBX_componentskeys[result - 1],
      s)) 
      return result;
  }
  result = 0;
  return result;
}  /* componentslookup */

Function(SYSTEM_ansichar *) GMSGLOBX_componentskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  if (n >= 1 && n <= 46) { 
    _P3strcpy(result,_len_ret,GMSGLOBX_componentskeys[n - 1]);
  } else 
    _P3strclr(result);
  return result;
}  /* componentskey */
typedef SYSTEM_uint8 _sub_13GMSGLOBX;
typedef _P3STR_3 _arr_12GMSGLOBX[62];
static _arr_12GMSGLOBX GMSGLOBX_clipcodeskeys = {{2,'0','0'}, {2,'0','1'}, {2,'0','2'}, {2,'0','3'}, {2,'0','4'}, {2,'0','5'}, {2,'0','M'}, {2,'0','S'}, {2,'A','T'}, {2,'B','A'}, {2,'B','D'}, {2,'C','L'}, {2,'C','O'}, {2,'C','P'}, {2,'C','D'}, {2,'C','M'}, {2,'D','E'}, {2,'D','I'}, {2,'E','C'}, {2,'F','B'}, {2,'F','R'}, {2,'G','E'}, {2,'G','S'}, {2,'G','Q'}, {2,'G','L'}, {2,'G','U'}, {2,'G','D'}, {2,'I','P'}, {2,'K','N'}, {2,'L','A'}, {2,'L','D'}, {2,'L','G'}, {2,'L','I'}, {2,'L','S'}, {2,'L','L'}, {2,'L','O'}, {2,'M','5'}, {2,'M','B'}, {2,'M','C'}, {2,'M','L'}, {2,'M','K'}, {2,'M','N'}, {2,'M','S'}, {2,'N','A'}, {2,'O','2'}, {2,'O','C'}, {2,'O','D'}, {2,'O','L'}, {2,'O','Q'}, {2,'O','S'}, {2,'P','T'}, {2,'S','B'}, {2,'S','C'}, {2,'S','E'}, {2,'S','N'}, {2,'X','A'}, {2,'X','L'}, {2,'X','P'}, {2,'X','S'}, {2,'X','X'}, {2,'W','Z'}, {2,'Z','O'}};

Function(SYSTEM_ansichar *) GMSGLOBX_clipcodestext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  switch (n) {
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\037GAMS/Demo GAMS Development Corp"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\037GAMS      GAMS Development Corp"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\037GAMS/Run  GAMS Development Corp"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\037GAMS/App  GAMS Development Corp"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\037GAMS/Node GAMS Development Corp"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\037GAMS/Comm GAMS Development Corp"));
      break;
    case 7: 
      _P3strcpy(result,_len_ret,_P3str1("\037GAMS/MIRO GAMS Development Corp"));
      break;
    case 8: 
      _P3strcpy(result,_len_ret,_P3str1("\037GAMS/Sec  GAMS Development Corp"));
      break;
    case 9: 
      _P3strcpy(result,_len_ret,_P3str1("\036ANTIGONE  Princeton Univeristy"));
      break;
    case 10: 
      _P3strcpy(result,_len_ret,_P3str1("\064BARON     University of Illinois at Urbana-Champaign"));
      break;
    case 11: 
      _P3strcpy(result,_len_ret,_P3str1("\037BDMLP     GAMS Development Corp"));
      break;
    case 12: 
      _P3strcpy(result,_len_ret,_P3str1("\050CPLEX/L   Cplex Optimization (Link only)"));
      break;
    case 13: 
      _P3strcpy(result,_len_ret,_P3str1("\031CONOPT    ARKI Consulting"));
      break;
    case 14: 
      _P3strcpy(result,_len_ret,_P3str1("\076CPLEX     Cplex Optimization (license options in comment line)"));
      break;
    case 15: 
      _P3strcpy(result,_len_ret,_P3str1("\054CPLEXDIST Cplex Optimization Distributed MIP"));
      break;
    case 16: 
      _P3strcpy(result,_len_ret,_P3str1("\034CPOPTIM   Cplex Optimization"));
      break;
    case 17: 
      _P3strcpy(result,_len_ret,_P3str1("\033DECIS     G. Infanger, Inc."));
      break;
    case 18: 
      _P3strcpy(result,_len_ret,_P3str1("\052DICOPT    CAPD, Carnegie Mellon University"));
      break;
    case 19: 
      _P3strcpy(result,_len_ret,_P3str1("\041ALPHAECP  Abo University, Finland"));
      break;
    case 20: 
      _P3strcpy(result,_len_ret,_P3str1("\042FILTERBB  The University of Dundee"));
      break;
    case 21: 
      _P3strcpy(result,_len_ret,_P3str1("\035FREE      Maintained Freeware"));
      break;
    case 22: 
      _P3strcpy(result,_len_ret,_P3str1("\033MPS/GE    Thomas Rutherford"));
      break;
    case 23: 
      _P3strcpy(result,_len_ret,_P3str1("\037GAMS/Sys  GAMS Development Corp"));
      break;
    case 24: 
      _P3strcpy(result,_len_ret,_P3str1("\036GloMIQO   Princeton Univeristy"));
      break;
    case 25: 
      _P3strcpy(result,_len_ret,_P3str1("\051GUROBI/L  Gurobi Optimization (Link only)"));
      break;
    case 26: 
      _P3strcpy(result,_len_ret,_P3str1("\035GUROBI    Gurobi Optimization"));
      break;
    case 27: 
      _P3strcpy(result,_len_ret,_P3str1("\035GUROBI/D  Gurobi Optimization"));
      break;
    case 28: 
      _P3strcpy(result,_len_ret,_P3str1("\034IPOPT     COIN-OR Foundation"));
      break;
    case 29: 
      _P3strcpy(result,_len_ret,_P3str1("\034KNITRO    Artelys           "));
      break;
    case 30: 
      _P3strcpy(result,_len_ret,_P3str1("\056LAMPS     Advanced Mathematical Software, Inc."));
      break;
    case 31: 
      _P3strcpy(result,_len_ret,_P3str1("\035LINDO-API Lindo Systems, Inc."));
      break;
    case 32: 
      _P3strcpy(result,_len_ret,_P3str1("\044LGO       Pinter Consulting Services"));
      break;
    case 33: 
      _P3strcpy(result,_len_ret,_P3str1("\035LINDOGlob Lindo Systems, Inc."));
      break;
    case 34: 
      _P3strcpy(result,_len_ret,_P3str1("\025LOCALSol  LocalSolver"));
      break;
    case 35: 
      _P3strcpy(result,_len_ret,_P3str1("\041LOCALS/L  LocalSolver (Link only)"));
      break;
    case 36: 
      _P3strcpy(result,_len_ret,_P3str1("\036LOQO      Princeton University"));
      break;
    case 37: 
      _P3strcpy(result,_len_ret,_P3str1("\035MINOS     Stanford University"));
      break;
    case 38: 
      _P3strcpy(result,_len_ret,_P3str1("\034MOSEK     EKA Consulting ApS"));
      break;
    case 39: 
      _P3strcpy(result,_len_ret,_P3str1("\037MILES     GAMS Development Corp"));
      break;
    case 40: 
      _P3strcpy(result,_len_ret,_P3str1("\034MOSEK/L   EKA Consulting ApS"));
      break;
    case 41: 
      _P3strcpy(result,_len_ret,_P3str1("\034MOSEK/MIP EKA Consulting ApS"));
      break;
    case 42: 
      _P3strcpy(result,_len_ret,_P3str1("\031MSNLP     Optimal Methods"));
      break;
    case 43: 
      _P3strcpy(result,_len_ret,_P3str1("\041MOPS      Freie University Berlin"));
      break;
    case 44: 
      _P3strcpy(result,_len_ret,_P3str1("\025NONO      Not Allowed"));
      break;
    case 45: 
      _P3strcpy(result,_len_ret,_P3str1("\033OSL2      IBM/OSL Version 2"));
      break;
    case 46: 
      _P3strcpy(result,_len_ret,_P3str1("\034CPLEXOSI  COIN-OR Foundation"));
      break;
    case 47: 
      _P3strcpy(result,_len_ret,_P3str1("\042ODHCPLEX  Optimization Direct Inc."));
      break;
    case 48: 
      _P3strcpy(result,_len_ret,_P3str1("\026OSL/L     IBM/OSL Link"));
      break;
    case 49: 
      _P3strcpy(result,_len_ret,_P3str1("\054OQNLP/GRG OptTek Systems and Optimal Methods"));
      break;
    case 50: 
      _P3strcpy(result,_len_ret,_P3str1("\021OSL       IBM/OSL"));
      break;
    case 51: 
      _P3strcpy(result,_len_ret,_P3str1("\053PATH      University of Wisconsin - Madison"));
      break;
    case 52: 
      _P3strcpy(result,_len_ret,_P3str1("\031SBB       ARKI Consulting"));
      break;
    case 53: 
      _P3strcpy(result,_len_ret,_P3str1("\024SCIP      ZIB Berlin"));
      break;
    case 54: 
      _P3strcpy(result,_len_ret,_P3str1("\046OSLSE     IBM/OSL Stochastic Extension"));
      break;
    case 55: 
      _P3strcpy(result,_len_ret,_P3str1("\035SNOPT     Stanford University"));
      break;
    case 56: 
      _P3strcpy(result,_len_ret,_P3str1("\031XA        Sunset Software"));
      break;
    case 57: 
      _P3strcpy(result,_len_ret,_P3str1("\016XPRESS/L  FICO"));
      break;
    case 58: 
      _P3strcpy(result,_len_ret,_P3str1("\016XPRESS    FICO"));
      break;
    case 59: 
      _P3strcpy(result,_len_ret,_P3str1("\016XPRESSSLP FICO"));
      break;
    case 60: 
      _P3strcpy(result,_len_ret,_P3str1("\016XPRESSALL FICO"));
      break;
    case 61: 
      _P3strcpy(result,_len_ret,_P3str1("\043WHIZZARD  Ketron Management Science"));
      break;
    case 62: 
      _P3strcpy(result,_len_ret,_P3str1("\043ZOOM      XMP Optimization Software"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030**** should never happen"));
  }
  return result;
}  /* clipcodestext */

Function(SYSTEM_integer ) GMSGLOBX_clipcodeslookup(
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s;
  SYSTEM_integer result;

  _P3strcpy(s,255,_ftmp1);
  for (result = 1;result <= (SYSTEM_int32)GMSGLOBX_maxclipcodes;++
    result) {
    if (SYSUTILS_P3_sametext(GMSGLOBX_clipcodeskeys[result - 1],s)) 
      return result;
  }
  result = 0;
  return result;
}  /* clipcodeslookup */

Function(SYSTEM_ansichar *) GMSGLOBX_clipcodeskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  if (n >= 1 && n <= 62) { 
    _P3strcpy(result,_len_ret,GMSGLOBX_clipcodeskeys[n - 1]);
  } else 
    _P3strclr(result);
  return result;
}  /* clipcodeskey */
typedef SYSTEM_uint8 _sub_15GMSGLOBX;
typedef _P3STR_3 _arr_14GMSGLOBX[6];
static _arr_14GMSGLOBX GMSGLOBX_gamslicenseskeys = {{2,'0','0'}, {2,'0','1'}, {2,'0','2'}, {2,'0','3'}, {2,'0','4'}, {2,'0','5'}};

Function(SYSTEM_ansichar *) GMSGLOBX_gamslicensestext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  switch (n) {
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\011GAMS/Demo"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\021GAMS/Professional"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\015GAMS/Run Time"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\025GAMS/Secured Run Time"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\015GAMS/Nodelock"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\016GAMS/Community"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030**** should never happen"));
  }
  return result;
}  /* gamslicensestext */

Function(SYSTEM_integer ) GMSGLOBX_gamslicenseslookup(
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s;
  SYSTEM_integer result;

  _P3strcpy(s,255,_ftmp1);
  for (result = 1;result <= (SYSTEM_int32)
    GMSGLOBX_maxgamslicenses;++result) {
    if (SYSUTILS_P3_sametext(GMSGLOBX_gamslicenseskeys[result - 1],
      s)) 
      return result;
  }
  result = 0;
  return result;
}  /* gamslicenseslookup */

Function(SYSTEM_ansichar *) GMSGLOBX_gamslicenseskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  if (n >= 1 && n <= 6) { 
    _P3strcpy(result,_len_ret,GMSGLOBX_gamslicenseskeys[n - 1]);
  } else 
    _P3strclr(result);
  return result;
}  /* gamslicenseskey */
typedef SYSTEM_uint8 _sub_17GMSGLOBX;
typedef _P3STR_3 _arr_16GMSGLOBX[20];
static _arr_16GMSGLOBX GMSGLOBX_gamslicensetypeskeys = {{1,'1'}, {1,'2'}, {1,'3'}, {1,'4'}, {1,'5'}, {1,'6'}, {1,'7'}, {1,'G'}, {1,'B'}, {1,'A'}, {1,'C'}, {1,'D'}, {1,'H'}, {1,'K'}, {1,'L'}, {1,'E'}, {1,'F'}, {1,'I'}, {1,'J'}, {1,'M'}};

Function(SYSTEM_ansichar *) GMSGLOBX_gamslicensetypestext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  switch (n) {
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\013Single User"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\022Small MUD - 5 User"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\024Medium MUD - 10 User"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\023Large MUD - 20 User"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\012Multi-User"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\011Mainframe"));
      break;
    case 7: 
      _P3strcpy(result,_len_ret,_P3str1("\005Other"));
      break;
    case 8: 
      _P3strcpy(result,_len_ret,_P3str1("\015MUD - 90 User"));
      break;
    case 9: 
      _P3strcpy(result,_len_ret,_P3str1("\015MUD - 40 User"));
      break;
    case 10: 
      _P3strcpy(result,_len_ret,_P3str1("\015MUD - 30 User"));
      break;
    case 11: 
      _P3strcpy(result,_len_ret,_P3str1("\015MUD - 50 User"));
      break;
    case 12: 
      _P3strcpy(result,_len_ret,_P3str1("\015MUD - 60 User"));
      break;
    case 13: 
      _P3strcpy(result,_len_ret,_P3str1("\016MUD - 100 User"));
      break;
    case 14: 
      _P3strcpy(result,_len_ret,_P3str1("\032Single Machine - 2 sockets"));
      break;
    case 15: 
      _P3strcpy(result,_len_ret,_P3str1("\032Single Machine - 3 sockets"));
      break;
    case 16: 
      _P3strcpy(result,_len_ret,_P3str1("\015MUD - 70 User"));
      break;
    case 17: 
      _P3strcpy(result,_len_ret,_P3str1("\015MUD - 80 User"));
      break;
    case 18: 
      _P3strcpy(result,_len_ret,_P3str1("\043Single Machine License - deprecated"));
      break;
    case 19: 
      _P3strcpy(result,_len_ret,_P3str1("\031Single Machine - 1 socket"));
      break;
    case 20: 
      _P3strcpy(result,_len_ret,_P3str1("\032Single Machine - 4 sockets"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030**** should never happen"));
  }
  return result;
}  /* gamslicensetypestext */

Function(SYSTEM_integer ) GMSGLOBX_gamslicensetypeslookup(
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s;
  SYSTEM_integer result;

  _P3strcpy(s,255,_ftmp1);
  for (result = 1;result <= (SYSTEM_int32)
    GMSGLOBX_maxgamslicensetypes;++result) {
    if (SYSUTILS_P3_sametext(GMSGLOBX_gamslicensetypeskeys[result - 1],
      s)) 
      return result;
  }
  result = 0;
  return result;
}  /* gamslicensetypeslookup */

Function(SYSTEM_ansichar *) GMSGLOBX_gamslicensetypeskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  if (n >= 1 && n <= 20) { 
    _P3strcpy(result,_len_ret,GMSGLOBX_gamslicensetypeskeys[n - 1]);
  } else 
    _P3strclr(result);
  return result;
}  /* gamslicensetypeskey */
typedef SYSTEM_uint8 _sub_19GMSGLOBX;
typedef _P3STR_7 _arr_18GMSGLOBX[16];
static _arr_18GMSGLOBX GMSGLOBX_modeltypesxkeys = {{4,'N','O','N','E'}, {2,'L','P'}, {3,'M','I','P'}, {4,'R','M','I','P'}, {3,'N','L','P'}, {3,'M','C','P'}, {4,'M','P','E','C'}, {5,'R','M','P','E','C'}, {3,'C','N','S'}, {4,'D','N','L','P'}, {6,'R','M','I','N','L','P'}, {5,'M','I','N','L','P'}, {3,'Q','C','P'}, {5,'M','I','Q','C','P'}, {6,'R','M','I','Q','C','P'}, {3,'E','M','P'}};

Function(SYSTEM_ansichar *) GMSGLOBX_modeltypesxtext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  switch (n) {
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\004NONE"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\022Linear Programming"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\031Mixed-Integer Programming"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\041Relaxed Mixed-Integer Programming"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\026Non-Linear Programming"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\036Mixed Complementarity Problems"));
      break;
    case 7: 
      _P3strcpy(result,_len_ret,_P3str1("\062Mathematical Programs with Equilibrium Constraints"));
      break;
    case 8: 
      _P3strcpy(result,_len_ret,_P3str1("\072Relaxed Mathematical Programs with Equilibrium Constraints"));
      break;
    case 9: 
      _P3strcpy(result,_len_ret,_P3str1("\035Constrained Nonlinear Systems"));
      break;
    case 10: 
      _P3strcpy(result,_len_ret,_P3str1("\065Non-Linear Programming with Discontinuous Derivatives"));
      break;
    case 11: 
      _P3strcpy(result,_len_ret,_P3str1("\054Relaxed Mixed-Integer Non-Linear Programming"));
      break;
    case 12: 
      _P3strcpy(result,_len_ret,_P3str1("\044Mixed-Integer Non-Linear Programming"));
      break;
    case 13: 
      _P3strcpy(result,_len_ret,_P3str1("\042Quadratically Constrained Programs"));
      break;
    case 14: 
      _P3strcpy(result,_len_ret,_P3str1("\060Mixed Integer Quadratically Constrained Programs"));
      break;
    case 15: 
      _P3strcpy(result,_len_ret,_P3str1("\070Relaxed Mixed Integer Quadratically Constrained Programs"));
      break;
    case 16: 
      _P3strcpy(result,_len_ret,_P3str1("\036Extended Mathematical Programs"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030**** should never happen"));
  }
  return result;
}  /* modeltypesxtext */

Function(SYSTEM_integer ) GMSGLOBX_modeltypesxlookup(
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s;
  SYSTEM_integer result;

  _P3strcpy(s,255,_ftmp1);
  for (result = 1;result <= (SYSTEM_int32)GMSGLOBX_maxmodeltypesx;++
    result) {
    if (SYSUTILS_P3_sametext(GMSGLOBX_modeltypesxkeys[result - 1],
      s)) 
      return result;
  }
  result = 0;
  return result;
}  /* modeltypesxlookup */

Function(SYSTEM_ansichar *) GMSGLOBX_modeltypesxkey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  if (n >= 1 && n <= 16) { 
    _P3strcpy(result,_len_ret,GMSGLOBX_modeltypesxkeys[n - 1]);
  } else 
    _P3strclr(result);
  return result;
}  /* modeltypesxkey */
typedef SYSTEM_uint8 _sub_21GMSGLOBX;
typedef SYSTEM_uint8 _sub_23GMSGLOBX;
typedef SYSTEM_integer _arr_22GMSGLOBX[2];
typedef _arr_22GMSGLOBX _arr_20GMSGLOBX[170];
static _arr_20GMSGLOBX GMSGLOBX_componentsolvermaptuple = {{1, 3}, 
  {2, 28}, {3, 29}, {4, 4}, {4, 13}, 
  {5, 5}, {5, 34}, {6, 7}, {6, 43}, 
  {6, 62}, {6, 149}, {7, 6}, {7, 39}, 
  {7, 40}, {7, 144}, {8, 7}, {8, 43}, 
  {8, 62}, {8, 149}, {10, 8}, {11, 9}, 
  {11, 44}, {11, 45}, {12, 10}, {12, 46}, 
  {13, 11}, {14, 12}, {15, 1}, {15, 2}, 
  {15, 15}, {15, 30}, {15, 31}, {15, 33}, 
  {15, 35}, {15, 36}, {15, 37}, {15, 38}, 
  {15, 41}, {15, 42}, {15, 47}, {15, 48}, 
  {15, 49}, {15, 50}, {15, 51}, {15, 52}, 
  {15, 56}, {15, 60}, {15, 61}, {15, 63}, 
  {15, 64}, {15, 65}, {15, 69}, {15, 70}, 
  {15, 71}, {15, 73}, {15, 74}, {15, 75}, 
  {15, 76}, {15, 77}, {15, 78}, {15, 79}, 
  {15, 80}, {15, 81}, {15, 82}, {15, 83}, 
  {15, 84}, {15, 85}, {15, 86}, {15, 87}, 
  {15, 88}, {15, 89}, {15, 90}, {15, 91}, 
  {15, 92}, {15, 93}, {15, 94}, {15, 95}, 
  {15, 96}, {15, 97}, {15, 98}, {15, 99}, 
  {15, 100}, {15, 101}, {15, 102}, {15, 103}, 
  {15, 104}, {15, 105}, {15, 106}, {15, 107}, 
  {15, 108}, {15, 109}, {15, 110}, {15, 111}, 
  {15, 112}, {15, 113}, {15, 114}, {15, 115}, 
  {15, 116}, {15, 117}, {15, 118}, {15, 119}, 
  {15, 120}, {15, 121}, {15, 122}, {15, 123}, 
  {15, 124}, {15, 125}, {15, 126}, {15, 127}, 
  {15, 128}, {15, 129}, {15, 130}, {15, 131}, 
  {15, 132}, {15, 133}, {15, 134}, {15, 135}, 
  {15, 136}, {15, 137}, {15, 138}, {15, 139}, 
  {15, 140}, {15, 141}, {15, 142}, {15, 143}, 
  {15, 145}, {15, 146}, {15, 147}, {15, 150}, 
  {15, 151}, {16, 13}, {17, 14}, {18, 14}, 
  {19, 14}, {21, 16}, {21, 53}, {22, 18}, 
  {22, 55}, {23, 17}, {23, 54}, {24, 55}, 
  {25, 19}, {25, 57}, {26, 19}, {26, 57}, 
  {27, 20}, {27, 58}, {27, 59}, {27, 152}, 
  {28, 32}, {29, 32}, {30, 32}, {31, 21}, 
  {33, 62}, {34, 22}, {36, 21}, {38, 23}, 
  {38, 66}, {38, 67}, {38, 68}, {39, 24}, 
  {40, 25}, {40, 72}, {40, 148}, {40, 153}, 
  {42, 26}, {43, 27}, {44, 27}, {45, 27}, 
  {46, 27}};

Function(SYSTEM_integer ) GMSGLOBX_componentsolvermapmap(
  SYSTEM_integer i,
  SYSTEM_integer j)
{
  SYSTEM_integer result;

  result = 0;
  if (i >= 1 && i <= 170) 
    if (j >= 1 && j <= 2) 
      result = GMSGLOBX_componentsolvermaptuple[i - 1][j - 1];
  return result;
}  /* componentsolvermapmap */
typedef SYSTEM_uint8 _sub_25GMSGLOBX;
typedef SYSTEM_uint8 _sub_27GMSGLOBX;
typedef SYSTEM_integer _arr_26GMSGLOBX[2];
typedef _arr_26GMSGLOBX _arr_24GMSGLOBX[51];
static _arr_24GMSGLOBX GMSGLOBX_clipcomponentmaptuple = {{2, 15}, {3, 15}, 
  {4, 15}, {5, 15}, {6, 15}, {7, 2}, 
  {8, 3}, {9, 4}, {10, 5}, {12, 6}, 
  {13, 7}, {14, 8}, {15, 9}, {16, 10}, 
  {17, 11}, {18, 12}, {19, 13}, {22, 14}, 
  {23, 15}, {24, 16}, {25, 17}, {26, 18}, 
  {27, 19}, {28, 20}, {29, 21}, {31, 22}, 
  {32, 23}, {33, 24}, {34, 25}, {35, 26}, 
  {37, 27}, {38, 28}, {40, 29}, {41, 30}, 
  {42, 31}, {45, 32}, {46, 33}, {47, 34}, 
  {48, 35}, {49, 36}, {50, 37}, {51, 38}, 
  {52, 39}, {53, 40}, {54, 41}, {55, 42}, 
  {56, 1}, {57, 43}, {58, 44}, {59, 45}, 
  {60, 46}};

Function(SYSTEM_integer ) GMSGLOBX_clipcomponentmapmap(
  SYSTEM_integer i,
  SYSTEM_integer j)
{
  SYSTEM_integer result;

  result = 0;
  if (i >= 1 && i <= 51) 
    if (j >= 1 && j <= 2) 
      result = GMSGLOBX_clipcomponentmaptuple[i - 1][j - 1];
  return result;
}  /* clipcomponentmapmap */
typedef SYSTEM_uint16 _sub_29GMSGLOBX;
typedef SYSTEM_uint8 _sub_31GMSGLOBX;
typedef SYSTEM_integer _arr_30GMSGLOBX[2];
typedef _arr_30GMSGLOBX _arr_28GMSGLOBX[413];
static _arr_28GMSGLOBX GMSGLOBX_solverplatformmaptuple = {{1, 1}, {1, 2}, 
  {1, 3}, {2, 1}, {2, 2}, {2, 3}, 
  {3, 1}, {3, 2}, {4, 1}, {4, 2}, 
  {4, 3}, {5, 1}, {5, 2}, {5, 3}, 
  {6, 1}, {6, 2}, {6, 3}, {7, 1}, 
  {7, 2}, {7, 3}, {8, 1}, {8, 2}, 
  {8, 3}, {9, 1}, {9, 2}, {9, 3}, 
  {10, 1}, {10, 2}, {10, 3}, {11, 1}, 
  {11, 2}, {11, 3}, {12, 1}, {12, 2}, 
  {12, 3}, {13, 1}, {13, 2}, {13, 3}, 
  {14, 1}, {14, 2}, {14, 3}, {15, 1}, 
  {15, 2}, {15, 3}, {16, 1}, {16, 2}, 
  {16, 3}, {17, 1}, {17, 2}, {17, 3}, 
  {18, 1}, {18, 2}, {18, 3}, {19, 1}, 
  {19, 2}, {19, 3}, {20, 1}, {20, 2}, 
  {20, 3}, {21, 1}, {21, 2}, {21, 3}, 
  {22, 1}, {22, 2}, {23, 1}, {23, 2}, 
  {23, 3}, {24, 1}, {24, 2}, {24, 3}, 
  {25, 1}, {25, 2}, {25, 3}, {26, 1}, 
  {26, 2}, {26, 3}, {27, 1}, {27, 2}, 
  {27, 3}, {28, 1}, {28, 2}, {28, 3}, 
  {29, 1}, {29, 2}, {29, 3}, {30, 1}, 
  {30, 2}, {30, 3}, {31, 1}, {31, 2}, 
  {31, 3}, {32, 1}, {32, 2}, {32, 3}, 
  {33, 1}, {33, 2}, {33, 3}, {35, 1}, 
  {35, 2}, {35, 3}, {36, 1}, {36, 2}, 
  {36, 3}, {37, 1}, {37, 2}, {37, 3}, 
  {38, 1}, {38, 2}, {38, 3}, {39, 1}, 
  {39, 2}, {39, 3}, {40, 1}, {40, 2}, 
  {40, 3}, {41, 1}, {41, 2}, {41, 3}, 
  {42, 1}, {42, 2}, {42, 3}, {43, 1}, 
  {43, 2}, {43, 3}, {44, 1}, {44, 2}, 
  {44, 3}, {45, 1}, {45, 2}, {45, 3}, 
  {47, 1}, {47, 2}, {47, 3}, {48, 1}, 
  {48, 2}, {48, 3}, {49, 1}, {49, 2}, 
  {49, 3}, {50, 1}, {50, 2}, {50, 3}, 
  {51, 1}, {51, 2}, {51, 3}, {52, 1}, 
  {52, 2}, {52, 3}, {54, 1}, {54, 2}, 
  {54, 3}, {55, 1}, {55, 2}, {55, 3}, 
  {56, 1}, {56, 2}, {56, 3}, {57, 1}, 
  {57, 2}, {57, 3}, {58, 1}, {58, 2}, 
  {58, 3}, {59, 1}, {59, 2}, {59, 3}, 
  {60, 1}, {60, 2}, {60, 3}, {61, 1}, 
  {61, 2}, {61, 3}, {62, 1}, {62, 2}, 
  {62, 3}, {63, 1}, {63, 2}, {63, 3}, 
  {64, 1}, {64, 2}, {64, 3}, {65, 1}, 
  {65, 2}, {65, 3}, {66, 1}, {66, 2}, 
  {66, 3}, {67, 1}, {67, 2}, {67, 3}, 
  {69, 1}, {69, 2}, {69, 3}, {70, 1}, 
  {70, 2}, {70, 3}, {71, 1}, {71, 2}, 
  {71, 3}, {72, 1}, {72, 2}, {72, 3}, 
  {73, 1}, {74, 1}, {74, 2}, {74, 3}, 
  {75, 1}, {75, 2}, {75, 3}, {76, 1}, 
  {76, 2}, {76, 3}, {77, 1}, {77, 2}, 
  {77, 3}, {78, 1}, {78, 2}, {78, 3}, 
  {79, 1}, {79, 2}, {79, 3}, {80, 1}, 
  {80, 2}, {80, 3}, {81, 1}, {82, 1}, 
  {83, 1}, {84, 1}, {84, 2}, {84, 3}, 
  {85, 1}, {85, 2}, {85, 3}, {86, 1}, 
  {87, 1}, {87, 2}, {87, 3}, {88, 1}, 
  {88, 2}, {88, 3}, {89, 1}, {89, 2}, 
  {89, 3}, {90, 1}, {90, 2}, {90, 3}, 
  {91, 1}, {91, 2}, {91, 3}, {92, 1}, 
  {92, 2}, {92, 3}, {93, 1}, {93, 2}, 
  {93, 3}, {94, 1}, {94, 2}, {94, 3}, 
  {95, 1}, {95, 2}, {95, 3}, {96, 1}, 
  {97, 1}, {98, 1}, {98, 2}, {98, 3}, 
  {99, 1}, {99, 2}, {99, 3}, {100, 1}, 
  {100, 2}, {100, 3}, {101, 1}, {102, 1}, 
  {102, 2}, {102, 3}, {103, 1}, {103, 2}, 
  {103, 3}, {104, 1}, {104, 2}, {104, 3}, 
  {105, 1}, {106, 1}, {107, 1}, {107, 2}, 
  {107, 3}, {108, 1}, {108, 2}, {108, 3}, 
  {109, 1}, {110, 1}, {110, 2}, {110, 3}, 
  {111, 1}, {111, 2}, {111, 3}, {112, 1}, 
  {113, 1}, {114, 1}, {114, 2}, {114, 3}, 
  {115, 1}, {115, 2}, {115, 3}, {116, 1}, 
  {117, 1}, {118, 1}, {119, 1}, {119, 2}, 
  {119, 3}, {120, 1}, {120, 2}, {120, 3}, 
  {121, 1}, {121, 2}, {121, 3}, {122, 1}, 
  {122, 2}, {122, 3}, {123, 1}, {123, 2}, 
  {123, 3}, {124, 1}, {124, 2}, {124, 3}, 
  {125, 1}, {125, 2}, {125, 3}, {126, 1}, 
  {126, 2}, {126, 3}, {127, 1}, {127, 2}, 
  {127, 3}, {128, 1}, {128, 2}, {128, 3}, 
  {129, 1}, {129, 2}, {129, 3}, {130, 1}, 
  {130, 2}, {130, 3}, {131, 1}, {131, 2}, 
  {131, 3}, {132, 1}, {132, 2}, {132, 3}, 
  {133, 1}, {133, 2}, {133, 3}, {134, 1}, 
  {134, 2}, {134, 3}, {135, 1}, {135, 2}, 
  {135, 3}, {136, 1}, {136, 2}, {136, 3}, 
  {137, 1}, {137, 2}, {137, 3}, {138, 1}, 
  {138, 2}, {138, 3}, {139, 1}, {139, 2}, 
  {139, 3}, {140, 1}, {140, 2}, {140, 3}, 
  {141, 1}, {141, 2}, {141, 3}, {142, 1}, 
  {142, 2}, {142, 3}, {143, 1}, {143, 2}, 
  {143, 3}, {144, 1}, {144, 2}, {144, 3}, 
  {145, 1}, {145, 2}, {145, 3}, {146, 1}, 
  {146, 2}, {146, 3}, {147, 1}, {147, 2}, 
  {147, 3}, {148, 1}, {148, 2}, {148, 3}, 
  {149, 1}, {149, 2}, {149, 3}, {150, 1}, 
  {150, 2}, {150, 3}, {151, 1}, {151, 2}, 
  {151, 3}, {152, 1}, {152, 2}, {152, 3}, 
  {153, 1}, {153, 2}, {153, 3}};

Function(SYSTEM_integer ) GMSGLOBX_solverplatformmapmap(
  SYSTEM_integer i,
  SYSTEM_integer j)
{
  SYSTEM_integer result;

  result = 0;
  if (i >= 1 && i <= 413) 
    if (j >= 1 && j <= 2) 
      result = GMSGLOBX_solverplatformmaptuple[i - 1][j - 1];
  return result;
}  /* solverplatformmapmap */
typedef SYSTEM_uint16 _sub_33GMSGLOBX;
typedef SYSTEM_uint8 _sub_35GMSGLOBX;
typedef SYSTEM_integer _arr_34GMSGLOBX[3];
typedef _arr_34GMSGLOBX _arr_32GMSGLOBX[1387];
static _arr_32GMSGLOBX GMSGLOBX_solvertypeplatformmaptuple = {{1, 16, 1}, 
  {1, 16, 2}, {1, 16, 3}, {2, 2, 1}, 
  {2, 2, 2}, {2, 2, 3}, {2, 4, 1}, 
  {2, 4, 2}, {2, 4, 3}, {3, 2, 1}, 
  {3, 2, 2}, {3, 3, 1}, {3, 3, 2}, 
  {3, 4, 1}, {3, 4, 2}, {4, 5, 1}, 
  {4, 5, 2}, {4, 5, 3}, {4, 9, 1}, 
  {4, 9, 2}, {4, 9, 3}, {4, 10, 1}, 
  {4, 10, 2}, {4, 10, 3}, {4, 11, 1}, 
  {4, 11, 2}, {4, 11, 3}, {4, 12, 1}, 
  {4, 12, 2}, {4, 12, 3}, {4, 13, 1}, 
  {4, 13, 2}, {4, 13, 3}, {4, 14, 1}, 
  {4, 14, 2}, {4, 14, 3}, {4, 15, 1}, 
  {4, 15, 2}, {4, 15, 3}, {5, 2, 1}, 
  {5, 2, 2}, {5, 2, 3}, {5, 3, 1}, 
  {5, 3, 2}, {5, 3, 3}, {5, 4, 1}, 
  {5, 4, 2}, {5, 4, 3}, {5, 5, 1}, 
  {5, 5, 2}, {5, 5, 3}, {5, 9, 1}, 
  {5, 9, 2}, {5, 9, 3}, {5, 10, 1}, 
  {5, 10, 2}, {5, 10, 3}, {5, 11, 1}, 
  {5, 11, 2}, {5, 11, 3}, {5, 12, 1}, 
  {5, 12, 2}, {5, 12, 3}, {5, 13, 1}, 
  {5, 13, 2}, {5, 13, 3}, {5, 14, 1}, 
  {5, 14, 2}, {5, 14, 3}, {5, 15, 1}, 
  {5, 15, 2}, {5, 15, 3}, {6, 2, 1}, 
  {6, 2, 2}, {6, 2, 3}, {6, 4, 1}, 
  {6, 4, 2}, {6, 4, 3}, {6, 5, 1}, 
  {6, 5, 2}, {6, 5, 3}, {6, 9, 1}, 
  {6, 9, 2}, {6, 9, 3}, {6, 10, 1}, 
  {6, 10, 2}, {6, 10, 3}, {6, 11, 1}, 
  {6, 11, 2}, {6, 11, 3}, {6, 13, 1}, 
  {6, 13, 2}, {6, 13, 3}, {6, 15, 1}, 
  {6, 15, 2}, {6, 15, 3}, {7, 2, 1}, 
  {7, 2, 2}, {7, 2, 3}, {7, 3, 1}, 
  {7, 3, 2}, {7, 3, 3}, {7, 4, 1}, 
  {7, 4, 2}, {7, 4, 3}, {7, 13, 1}, 
  {7, 13, 2}, {7, 13, 3}, {7, 14, 1}, 
  {7, 14, 2}, {7, 14, 3}, {7, 15, 1}, 
  {7, 15, 2}, {7, 15, 3}, {8, 3, 1}, 
  {8, 3, 2}, {8, 3, 3}, {8, 12, 1}, 
  {8, 12, 2}, {8, 12, 3}, {8, 14, 1}, 
  {8, 14, 2}, {8, 14, 3}, {9, 16, 1}, 
  {9, 16, 2}, {9, 16, 3}, {10, 12, 1}, 
  {10, 12, 2}, {10, 12, 3}, {10, 14, 1}, 
  {10, 14, 2}, {10, 14, 3}, {11, 12, 1}, 
  {11, 12, 2}, {11, 12, 3}, {11, 14, 1}, 
  {11, 14, 2}, {11, 14, 3}, {13, 13, 1}, 
  {13, 13, 2}, {13, 13, 3}, {13, 14, 1}, 
  {13, 14, 2}, {13, 14, 3}, {13, 15, 1}, 
  {13, 15, 2}, {13, 15, 3}, {14, 2, 1}, 
  {14, 2, 2}, {14, 2, 3}, {14, 3, 1}, 
  {14, 3, 2}, {14, 3, 3}, {14, 4, 1}, 
  {14, 4, 2}, {14, 4, 3}, {14, 13, 1}, 
  {14, 13, 2}, {14, 13, 3}, {14, 14, 1}, 
  {14, 14, 2}, {14, 14, 3}, {14, 15, 1}, 
  {14, 15, 2}, {14, 15, 3}, {15, 2, 1}, 
  {15, 2, 2}, {15, 2, 3}, {15, 4, 1}, 
  {15, 4, 2}, {15, 4, 3}, {15, 5, 1}, 
  {15, 5, 2}, {15, 5, 3}, {15, 9, 1}, 
  {15, 9, 2}, {15, 9, 3}, {15, 10, 1}, 
  {15, 10, 2}, {15, 10, 3}, {15, 11, 1}, 
  {15, 11, 2}, {15, 11, 3}, {15, 13, 1}, 
  {15, 13, 2}, {15, 13, 3}, {15, 15, 1}, 
  {15, 15, 2}, {15, 15, 3}, {16, 2, 1}, 
  {16, 2, 2}, {16, 2, 3}, {16, 4, 1}, 
  {16, 4, 2}, {16, 4, 3}, {16, 5, 1}, 
  {16, 5, 2}, {16, 5, 3}, {16, 7, 1}, 
  {16, 7, 2}, {16, 7, 3}, {16, 8, 1}, 
  {16, 8, 2}, {16, 8, 3}, {16, 9, 1}, 
  {16, 9, 2}, {16, 9, 3}, {16, 10, 1}, 
  {16, 10, 2}, {16, 10, 3}, {16, 11, 1}, 
  {16, 11, 2}, {16, 11, 3}, {16, 12, 1}, 
  {16, 12, 2}, {16, 12, 3}, {16, 13, 1}, 
  {16, 13, 2}, {16, 13, 3}, {16, 14, 1}, 
  {16, 14, 2}, {16, 14, 3}, {16, 15, 1}, 
  {16, 15, 2}, {16, 15, 3}, {17, 2, 1}, 
  {17, 2, 2}, {17, 2, 3}, {17, 4, 1}, 
  {17, 4, 2}, {17, 4, 3}, {17, 5, 1}, 
  {17, 5, 2}, {17, 5, 3}, {17, 10, 1}, 
  {17, 10, 2}, {17, 10, 3}, {17, 11, 1}, 
  {17, 11, 2}, {17, 11, 3}, {17, 13, 1}, 
  {17, 13, 2}, {17, 13, 3}, {17, 15, 1}, 
  {17, 15, 2}, {17, 15, 3}, {18, 2, 1}, 
  {18, 2, 2}, {18, 2, 3}, {18, 3, 1}, 
  {18, 3, 2}, {18, 3, 3}, {18, 4, 1}, 
  {18, 4, 2}, {18, 4, 3}, {18, 5, 1}, 
  {18, 5, 2}, {18, 5, 3}, {18, 10, 1}, 
  {18, 10, 2}, {18, 10, 3}, {18, 11, 1}, 
  {18, 11, 2}, {18, 11, 3}, {18, 12, 1}, 
  {18, 12, 2}, {18, 12, 3}, {18, 13, 1}, 
  {18, 13, 2}, {18, 13, 3}, {18, 14, 1}, 
  {18, 14, 2}, {18, 14, 3}, {18, 15, 1}, 
  {18, 15, 2}, {18, 15, 3}, {18, 16, 1}, 
  {18, 16, 2}, {18, 16, 3}, {19, 3, 1}, 
  {19, 3, 2}, {19, 3, 3}, {19, 5, 1}, 
  {19, 5, 2}, {19, 5, 3}, {19, 9, 1}, 
  {19, 9, 2}, {19, 9, 3}, {19, 10, 1}, 
  {19, 10, 2}, {19, 10, 3}, {19, 11, 1}, 
  {19, 11, 2}, {19, 11, 3}, {19, 12, 1}, 
  {19, 12, 2}, {19, 12, 3}, {19, 13, 1}, 
  {19, 13, 2}, {19, 13, 3}, {19, 14, 1}, 
  {19, 14, 2}, {19, 14, 3}, {19, 15, 1}, 
  {19, 15, 2}, {19, 15, 3}, {20, 2, 1}, 
  {20, 2, 2}, {20, 2, 3}, {20, 4, 1}, 
  {20, 4, 2}, {20, 4, 3}, {20, 5, 1}, 
  {20, 5, 2}, {20, 5, 3}, {20, 9, 1}, 
  {20, 9, 2}, {20, 9, 3}, {20, 10, 1}, 
  {20, 10, 2}, {20, 10, 3}, {20, 11, 1}, 
  {20, 11, 2}, {20, 11, 3}, {20, 13, 1}, 
  {20, 13, 2}, {20, 13, 3}, {20, 15, 1}, 
  {20, 15, 2}, {20, 15, 3}, {21, 5, 1}, 
  {21, 5, 2}, {21, 5, 3}, {21, 10, 1}, 
  {21, 10, 2}, {21, 10, 3}, {21, 11, 1}, 
  {21, 11, 2}, {21, 11, 3}, {21, 13, 1}, 
  {21, 13, 2}, {21, 13, 3}, {21, 15, 1}, 
  {21, 15, 2}, {21, 15, 3}, {22, 3, 1}, 
  {22, 3, 2}, {22, 14, 1}, {22, 14, 2}, 
  {23, 6, 1}, {23, 6, 2}, {23, 6, 3}, 
  {23, 9, 1}, {23, 9, 2}, {23, 9, 3}, 
  {24, 12, 1}, {24, 12, 2}, {24, 12, 3}, 
  {24, 14, 1}, {24, 14, 2}, {24, 14, 3}, 
  {25, 3, 1}, {25, 3, 2}, {25, 3, 3}, 
  {25, 5, 1}, {25, 5, 2}, {25, 5, 3}, 
  {25, 9, 1}, {25, 9, 2}, {25, 9, 3}, 
  {25, 10, 1}, {25, 10, 2}, {25, 10, 3}, 
  {25, 11, 1}, {25, 11, 2}, {25, 11, 3}, 
  {25, 12, 1}, {25, 12, 2}, {25, 12, 3}, 
  {25, 13, 1}, {25, 13, 2}, {25, 13, 3}, 
  {25, 14, 1}, {25, 14, 2}, {25, 14, 3}, 
  {25, 15, 1}, {25, 15, 2}, {25, 15, 3}, 
  {26, 2, 1}, {26, 2, 2}, {26, 2, 3}, 
  {26, 4, 1}, {26, 4, 2}, {26, 4, 3}, 
  {26, 5, 1}, {26, 5, 2}, {26, 5, 3}, 
  {26, 9, 1}, {26, 9, 2}, {26, 9, 3}, 
  {26, 10, 1}, {26, 10, 2}, {26, 10, 3}, 
  {26, 11, 1}, {26, 11, 2}, {26, 11, 3}, 
  {26, 13, 1}, {26, 13, 2}, {26, 13, 3}, 
  {26, 15, 1}, {26, 15, 2}, {26, 15, 3}, 
  {27, 2, 1}, {27, 2, 2}, {27, 2, 3}, 
  {27, 3, 1}, {27, 3, 2}, {27, 3, 3}, 
  {27, 4, 1}, {27, 4, 2}, {27, 4, 3}, 
  {27, 5, 1}, {27, 5, 2}, {27, 5, 3}, 
  {27, 9, 1}, {27, 9, 2}, {27, 9, 3}, 
  {27, 10, 1}, {27, 10, 2}, {27, 10, 3}, 
  {27, 11, 1}, {27, 11, 2}, {27, 11, 3}, 
  {27, 12, 1}, {27, 12, 2}, {27, 12, 3}, 
  {27, 13, 1}, {27, 13, 2}, {27, 13, 3}, 
  {27, 14, 1}, {27, 14, 2}, {27, 14, 3}, 
  {27, 15, 1}, {27, 15, 2}, {27, 15, 3}, 
  {30, 2, 1}, {30, 2, 2}, {30, 2, 3}, 
  {30, 3, 1}, {30, 3, 2}, {30, 3, 3}, 
  {30, 4, 1}, {30, 4, 2}, {30, 4, 3}, 
  {31, 6, 1}, {31, 6, 2}, {31, 6, 3}, 
  {32, 2, 1}, {32, 2, 2}, {32, 2, 3}, 
  {32, 3, 1}, {32, 3, 2}, {32, 3, 3}, 
  {32, 4, 1}, {32, 4, 2}, {32, 4, 3}, 
  {32, 5, 1}, {32, 5, 2}, {32, 5, 3}, 
  {32, 10, 1}, {32, 10, 2}, {32, 10, 3}, 
  {32, 11, 1}, {32, 11, 2}, {32, 11, 3}, 
  {32, 12, 1}, {32, 12, 2}, {32, 12, 3}, 
  {32, 13, 1}, {32, 13, 2}, {32, 13, 3}, 
  {32, 14, 1}, {32, 14, 2}, {32, 14, 3}, 
  {32, 15, 1}, {32, 15, 2}, {32, 15, 3}, 
  {33, 2, 1}, {33, 2, 2}, {33, 2, 3}, 
  {33, 3, 1}, {33, 3, 2}, {33, 3, 3}, 
  {33, 4, 1}, {33, 4, 2}, {33, 4, 3}, 
  {33, 5, 1}, {33, 5, 2}, {33, 5, 3}, 
  {33, 6, 1}, {33, 6, 2}, {33, 6, 3}, 
  {33, 7, 1}, {33, 7, 2}, {33, 7, 3}, 
  {33, 8, 1}, {33, 8, 2}, {33, 8, 3}, 
  {33, 9, 1}, {33, 9, 2}, {33, 9, 3}, 
  {33, 10, 1}, {33, 10, 2}, {33, 10, 3}, 
  {33, 11, 1}, {33, 11, 2}, {33, 11, 3}, 
  {33, 12, 1}, {33, 12, 2}, {33, 12, 3}, 
  {35, 2, 1}, {35, 2, 2}, {35, 2, 3}, 
  {35, 4, 1}, {35, 4, 2}, {35, 4, 3}, 
  {36, 2, 1}, {36, 2, 2}, {36, 2, 3}, 
  {36, 3, 1}, {36, 3, 2}, {36, 3, 3}, 
  {36, 4, 1}, {36, 4, 2}, {36, 4, 3}, 
  {36, 5, 1}, {36, 5, 2}, {36, 5, 3}, 
  {36, 6, 1}, {36, 6, 2}, {36, 6, 3}, 
  {36, 7, 1}, {36, 7, 2}, {36, 7, 3}, 
  {36, 8, 1}, {36, 8, 2}, {36, 8, 3}, 
  {36, 9, 1}, {36, 9, 2}, {36, 9, 3}, 
  {36, 10, 1}, {36, 10, 2}, {36, 10, 3}, 
  {36, 11, 1}, {36, 11, 2}, {36, 11, 3}, 
  {36, 12, 1}, {36, 12, 2}, {36, 12, 3}, 
  {36, 13, 1}, {36, 13, 2}, {36, 13, 3}, 
  {36, 14, 1}, {36, 14, 2}, {36, 14, 3}, 
  {36, 15, 1}, {36, 15, 2}, {36, 15, 3}, 
  {37, 12, 1}, {37, 12, 2}, {37, 12, 3}, 
  {37, 14, 1}, {37, 14, 2}, {37, 14, 3}, 
  {38, 2, 1}, {38, 2, 2}, {38, 2, 3}, 
  {38, 3, 1}, {38, 3, 2}, {38, 3, 3}, 
  {38, 4, 1}, {38, 4, 2}, {38, 4, 3}, 
  {39, 2, 1}, {39, 2, 2}, {39, 2, 3}, 
  {39, 4, 1}, {39, 4, 2}, {39, 4, 3}, 
  {39, 5, 1}, {39, 5, 2}, {39, 5, 3}, 
  {39, 9, 1}, {39, 9, 2}, {39, 9, 3}, 
  {39, 10, 1}, {39, 10, 2}, {39, 10, 3}, 
  {39, 11, 1}, {39, 11, 2}, {39, 11, 3}, 
  {39, 13, 1}, {39, 13, 2}, {39, 13, 3}, 
  {39, 15, 1}, {39, 15, 2}, {39, 15, 3}, 
  {40, 2, 1}, {40, 2, 2}, {40, 2, 3}, 
  {40, 4, 1}, {40, 4, 2}, {40, 4, 3}, 
  {40, 5, 1}, {40, 5, 2}, {40, 5, 3}, 
  {40, 9, 1}, {40, 9, 2}, {40, 9, 3}, 
  {40, 10, 1}, {40, 10, 2}, {40, 10, 3}, 
  {40, 11, 1}, {40, 11, 2}, {40, 11, 3}, 
  {40, 13, 1}, {40, 13, 2}, {40, 13, 3}, 
  {40, 15, 1}, {40, 15, 2}, {40, 15, 3}, 
  {41, 2, 1}, {41, 2, 2}, {41, 2, 3}, 
  {41, 3, 1}, {41, 3, 2}, {41, 3, 3}, 
  {41, 4, 1}, {41, 4, 2}, {41, 4, 3}, 
  {41, 5, 1}, {41, 5, 2}, {41, 5, 3}, 
  {41, 6, 1}, {41, 6, 2}, {41, 6, 3}, 
  {41, 7, 1}, {41, 7, 2}, {41, 7, 3}, 
  {41, 8, 1}, {41, 8, 2}, {41, 8, 3}, 
  {41, 9, 1}, {41, 9, 2}, {41, 9, 3}, 
  {41, 10, 1}, {41, 10, 2}, {41, 10, 3}, 
  {41, 11, 1}, {41, 11, 2}, {41, 11, 3}, 
  {41, 12, 1}, {41, 12, 2}, {41, 12, 3}, 
  {41, 13, 1}, {41, 13, 2}, {41, 13, 3}, 
  {41, 14, 1}, {41, 14, 2}, {41, 14, 3}, 
  {41, 15, 1}, {41, 15, 2}, {41, 15, 3}, 
  {42, 2, 1}, {42, 2, 2}, {42, 2, 3}, 
  {42, 3, 1}, {42, 3, 2}, {42, 3, 3}, 
  {42, 4, 1}, {42, 4, 2}, {42, 4, 3}, 
  {42, 5, 1}, {42, 5, 2}, {42, 5, 3}, 
  {42, 6, 1}, {42, 6, 2}, {42, 6, 3}, 
  {42, 7, 1}, {42, 7, 2}, {42, 7, 3}, 
  {42, 8, 1}, {42, 8, 2}, {42, 8, 3}, 
  {42, 9, 1}, {42, 9, 2}, {42, 9, 3}, 
  {42, 10, 1}, {42, 10, 2}, {42, 10, 3}, 
  {42, 11, 1}, {42, 11, 2}, {42, 11, 3}, 
  {42, 12, 1}, {42, 12, 2}, {42, 12, 3}, 
  {42, 13, 1}, {42, 13, 2}, {42, 13, 3}, 
  {42, 14, 1}, {42, 14, 2}, {42, 14, 3}, 
  {42, 15, 1}, {42, 15, 2}, {42, 15, 3}, 
  {42, 16, 1}, {42, 16, 2}, {42, 16, 3}, 
  {43, 2, 1}, {43, 2, 2}, {43, 2, 3}, 
  {43, 3, 1}, {43, 3, 2}, {43, 3, 3}, 
  {43, 4, 1}, {43, 4, 2}, {43, 4, 3}, 
  {43, 13, 1}, {43, 13, 2}, {43, 13, 3}, 
  {43, 14, 1}, {43, 14, 2}, {43, 14, 3}, 
  {43, 15, 1}, {43, 15, 2}, {43, 15, 3}, 
  {44, 2, 1}, {44, 2, 2}, {44, 2, 3}, 
  {45, 2, 1}, {45, 2, 2}, {45, 2, 3}, 
  {47, 2, 1}, {47, 2, 2}, {47, 2, 3}, 
  {47, 3, 1}, {47, 3, 2}, {47, 3, 3}, 
  {47, 4, 1}, {47, 4, 2}, {47, 4, 3}, 
  {47, 5, 1}, {47, 5, 2}, {47, 5, 3}, 
  {47, 6, 1}, {47, 6, 2}, {47, 6, 3}, 
  {47, 7, 1}, {47, 7, 2}, {47, 7, 3}, 
  {47, 8, 1}, {47, 8, 2}, {47, 8, 3}, 
  {47, 10, 1}, {47, 10, 2}, {47, 10, 3}, 
  {47, 11, 1}, {47, 11, 2}, {47, 11, 3}, 
  {47, 12, 1}, {47, 12, 2}, {47, 12, 3}, 
  {47, 13, 1}, {47, 13, 2}, {47, 13, 3}, 
  {47, 14, 1}, {47, 14, 2}, {47, 14, 3}, 
  {47, 15, 1}, {47, 15, 2}, {47, 15, 3}, 
  {48, 2, 1}, {48, 2, 2}, {48, 2, 3}, 
  {48, 3, 1}, {48, 3, 2}, {48, 3, 3}, 
  {48, 4, 1}, {48, 4, 2}, {48, 4, 3}, 
  {48, 5, 1}, {48, 5, 2}, {48, 5, 3}, 
  {48, 6, 1}, {48, 6, 2}, {48, 6, 3}, 
  {48, 10, 1}, {48, 10, 2}, {48, 10, 3}, 
  {48, 11, 1}, {48, 11, 2}, {48, 11, 3}, 
  {48, 12, 1}, {48, 12, 2}, {48, 12, 3}, 
  {48, 13, 1}, {48, 13, 2}, {48, 13, 3}, 
  {48, 14, 1}, {48, 14, 2}, {48, 14, 3}, 
  {48, 15, 1}, {48, 15, 2}, {48, 15, 3}, 
  {49, 2, 1}, {49, 2, 2}, {49, 2, 3}, 
  {49, 3, 1}, {49, 3, 2}, {49, 3, 3}, 
  {49, 4, 1}, {49, 4, 2}, {49, 4, 3}, 
  {49, 5, 1}, {49, 5, 2}, {49, 5, 3}, 
  {49, 6, 1}, {49, 6, 2}, {49, 6, 3}, 
  {49, 10, 1}, {49, 10, 2}, {49, 10, 3}, 
  {49, 11, 1}, {49, 11, 2}, {49, 11, 3}, 
  {49, 12, 1}, {49, 12, 2}, {49, 12, 3}, 
  {49, 13, 1}, {49, 13, 2}, {49, 13, 3}, 
  {49, 14, 1}, {49, 14, 2}, {49, 14, 3}, 
  {49, 15, 1}, {49, 15, 2}, {49, 15, 3}, 
  {50, 2, 1}, {50, 2, 2}, {50, 2, 3}, 
  {50, 3, 1}, {50, 3, 2}, {50, 3, 3}, 
  {50, 4, 1}, {50, 4, 2}, {50, 4, 3}, 
  {50, 5, 1}, {50, 5, 2}, {50, 5, 3}, 
  {50, 6, 1}, {50, 6, 2}, {50, 6, 3}, 
  {50, 9, 1}, {50, 9, 2}, {50, 9, 3}, 
  {50, 10, 1}, {50, 10, 2}, {50, 10, 3}, 
  {50, 11, 1}, {50, 11, 2}, {50, 11, 3}, 
  {50, 12, 1}, {50, 12, 2}, {50, 12, 3}, 
  {50, 13, 1}, {50, 13, 2}, {50, 13, 3}, 
  {50, 14, 1}, {50, 14, 2}, {50, 14, 3}, 
  {50, 15, 1}, {50, 15, 2}, {50, 15, 3}, 
  {51, 16, 1}, {51, 16, 2}, {51, 16, 3}, 
  {52, 2, 1}, {52, 2, 2}, {52, 2, 3}, 
  {52, 3, 1}, {52, 3, 2}, {52, 3, 3}, 
  {52, 4, 1}, {52, 4, 2}, {52, 4, 3}, 
  {52, 5, 1}, {52, 5, 2}, {52, 5, 3}, 
  {52, 6, 1}, {52, 6, 2}, {52, 6, 3}, 
  {52, 7, 1}, {52, 7, 2}, {52, 7, 3}, 
  {52, 8, 1}, {52, 8, 2}, {52, 8, 3}, 
  {52, 9, 1}, {52, 9, 2}, {52, 9, 3}, 
  {52, 10, 1}, {52, 10, 2}, {52, 10, 3}, 
  {52, 11, 1}, {52, 11, 2}, {52, 11, 3}, 
  {52, 12, 1}, {52, 12, 2}, {52, 12, 3}, 
  {52, 13, 1}, {52, 13, 2}, {52, 13, 3}, 
  {52, 14, 1}, {52, 14, 2}, {52, 14, 3}, 
  {52, 15, 1}, {52, 15, 2}, {52, 15, 3}, 
  {52, 16, 1}, {52, 16, 2}, {52, 16, 3}, 
  {54, 2, 1}, {54, 2, 2}, {54, 2, 3}, 
  {54, 4, 1}, {54, 4, 2}, {54, 4, 3}, 
  {54, 5, 1}, {54, 5, 2}, {54, 5, 3}, 
  {54, 10, 1}, {54, 10, 2}, {54, 10, 3}, 
  {54, 11, 1}, {54, 11, 2}, {54, 11, 3}, 
  {54, 13, 1}, {54, 13, 2}, {54, 13, 3}, 
  {54, 15, 1}, {54, 15, 2}, {54, 15, 3}, 
  {55, 2, 1}, {55, 2, 2}, {55, 2, 3}, 
  {55, 3, 1}, {55, 3, 2}, {55, 3, 3}, 
  {55, 4, 1}, {55, 4, 2}, {55, 4, 3}, 
  {55, 5, 1}, {55, 5, 2}, {55, 5, 3}, 
  {55, 10, 1}, {55, 10, 2}, {55, 10, 3}, 
  {55, 11, 1}, {55, 11, 2}, {55, 11, 3}, 
  {55, 12, 1}, {55, 12, 2}, {55, 12, 3}, 
  {55, 13, 1}, {55, 13, 2}, {55, 13, 3}, 
  {55, 14, 1}, {55, 14, 2}, {55, 14, 3}, 
  {55, 15, 1}, {55, 15, 2}, {55, 15, 3}, 
  {56, 2, 1}, {56, 2, 2}, {56, 2, 3}, 
  {56, 3, 1}, {56, 3, 2}, {56, 3, 3}, 
  {56, 4, 1}, {56, 4, 2}, {56, 4, 3}, 
  {56, 5, 1}, {56, 5, 2}, {56, 5, 3}, 
  {56, 10, 1}, {56, 10, 2}, {56, 10, 3}, 
  {56, 11, 1}, {56, 11, 2}, {56, 11, 3}, 
  {56, 12, 1}, {56, 12, 2}, {56, 12, 3}, 
  {57, 3, 1}, {57, 3, 2}, {57, 3, 3}, 
  {57, 5, 1}, {57, 5, 2}, {57, 5, 3}, 
  {57, 9, 1}, {57, 9, 2}, {57, 9, 3}, 
  {57, 10, 1}, {57, 10, 2}, {57, 10, 3}, 
  {57, 11, 1}, {57, 11, 2}, {57, 11, 3}, 
  {57, 12, 1}, {57, 12, 2}, {57, 12, 3}, 
  {57, 13, 1}, {57, 13, 2}, {57, 13, 3}, 
  {57, 14, 1}, {57, 14, 2}, {57, 14, 3}, 
  {57, 15, 1}, {57, 15, 2}, {57, 15, 3}, 
  {58, 2, 1}, {58, 2, 2}, {58, 2, 3}, 
  {58, 4, 1}, {58, 4, 2}, {58, 4, 3}, 
  {58, 5, 1}, {58, 5, 2}, {58, 5, 3}, 
  {58, 9, 1}, {58, 9, 2}, {58, 9, 3}, 
  {58, 10, 1}, {58, 10, 2}, {58, 10, 3}, 
  {58, 11, 1}, {58, 11, 2}, {58, 11, 3}, 
  {58, 13, 1}, {58, 13, 2}, {58, 13, 3}, 
  {58, 15, 1}, {58, 15, 2}, {58, 15, 3}, 
  {59, 2, 1}, {59, 2, 2}, {59, 2, 3}, 
  {59, 4, 1}, {59, 4, 2}, {59, 4, 3}, 
  {60, 2, 1}, {60, 2, 2}, {60, 2, 3}, 
  {60, 3, 1}, {60, 3, 2}, {60, 3, 3}, 
  {60, 4, 1}, {60, 4, 2}, {60, 4, 3}, 
  {60, 5, 1}, {60, 5, 2}, {60, 5, 3}, 
  {60, 6, 1}, {60, 6, 2}, {60, 6, 3}, 
  {60, 7, 1}, {60, 7, 2}, {60, 7, 3}, 
  {60, 8, 1}, {60, 8, 2}, {60, 8, 3}, 
  {60, 9, 1}, {60, 9, 2}, {60, 9, 3}, 
  {60, 10, 1}, {60, 10, 2}, {60, 10, 3}, 
  {60, 11, 1}, {60, 11, 2}, {60, 11, 3}, 
  {60, 12, 1}, {60, 12, 2}, {60, 12, 3}, 
  {61, 6, 1}, {61, 6, 2}, {61, 6, 3}, 
  {61, 7, 1}, {61, 7, 2}, {61, 7, 3}, 
  {61, 8, 1}, {61, 8, 2}, {61, 8, 3}, 
  {62, 2, 1}, {62, 2, 2}, {62, 2, 3}, 
  {62, 3, 1}, {62, 3, 2}, {62, 3, 3}, 
  {62, 4, 1}, {62, 4, 2}, {62, 4, 3}, 
  {63, 2, 1}, {63, 2, 2}, {63, 2, 3}, 
  {63, 3, 1}, {63, 3, 2}, {63, 3, 3}, 
  {63, 4, 1}, {63, 4, 2}, {63, 4, 3}, 
  {64, 2, 1}, {64, 2, 2}, {64, 2, 3}, 
  {64, 3, 1}, {64, 3, 2}, {64, 3, 3}, 
  {64, 4, 1}, {64, 4, 2}, {64, 4, 3}, 
  {65, 2, 1}, {65, 2, 2}, {65, 2, 3}, 
  {65, 3, 1}, {65, 3, 2}, {65, 3, 3}, 
  {65, 4, 1}, {65, 4, 2}, {65, 4, 3}, 
  {66, 2, 1}, {66, 2, 2}, {66, 2, 3}, 
  {66, 4, 1}, {66, 4, 2}, {66, 4, 3}, 
  {66, 5, 1}, {66, 5, 2}, {66, 5, 3}, 
  {66, 10, 1}, {66, 10, 2}, {66, 10, 3}, 
  {66, 11, 1}, {66, 11, 2}, {66, 11, 3}, 
  {66, 13, 1}, {66, 13, 2}, {66, 13, 3}, 
  {66, 15, 1}, {66, 15, 2}, {66, 15, 3}, 
  {67, 6, 1}, {67, 6, 2}, {67, 6, 3}, 
  {67, 9, 1}, {67, 9, 2}, {67, 9, 3}, 
  {69, 2, 1}, {69, 2, 2}, {69, 2, 3}, 
  {69, 3, 1}, {69, 3, 2}, {69, 3, 3}, 
  {69, 4, 1}, {69, 4, 2}, {69, 4, 3}, 
  {69, 5, 1}, {69, 5, 2}, {69, 5, 3}, 
  {69, 6, 1}, {69, 6, 2}, {69, 6, 3}, 
  {69, 7, 1}, {69, 7, 2}, {69, 7, 3}, 
  {69, 8, 1}, {69, 8, 2}, {69, 8, 3}, 
  {69, 9, 1}, {69, 9, 2}, {69, 9, 3}, 
  {69, 10, 1}, {69, 10, 2}, {69, 10, 3}, 
  {69, 11, 1}, {69, 11, 2}, {69, 11, 3}, 
  {69, 12, 1}, {69, 12, 2}, {69, 12, 3}, 
  {70, 16, 1}, {70, 16, 2}, {70, 16, 3}, 
  {71, 12, 1}, {71, 12, 2}, {71, 12, 3}, 
  {71, 14, 1}, {71, 14, 2}, {71, 14, 3}, 
  {72, 2, 1}, {72, 2, 2}, {72, 2, 3}, 
  {72, 4, 1}, {72, 4, 2}, {72, 4, 3}, 
  {142, 2, 1}, {142, 2, 2}, {142, 2, 3}, 
  {142, 4, 1}, {142, 4, 2}, {142, 4, 3}, 
  {142, 5, 1}, {142, 5, 2}, {142, 5, 3}, 
  {142, 9, 1}, {142, 9, 2}, {142, 9, 3}, 
  {142, 10, 1}, {142, 10, 2}, {142, 10, 3}, 
  {142, 11, 1}, {142, 11, 2}, {142, 11, 3}, 
  {142, 13, 1}, {142, 13, 2}, {142, 13, 3}, 
  {142, 15, 1}, {142, 15, 2}, {142, 15, 3}, 
  {143, 12, 1}, {143, 12, 2}, {143, 12, 3}, 
  {143, 14, 1}, {143, 14, 2}, {143, 14, 3}, 
  {144, 2, 1}, {144, 2, 2}, {144, 2, 3}, 
  {144, 4, 1}, {144, 4, 2}, {144, 4, 3}, 
  {144, 5, 1}, {144, 5, 2}, {144, 5, 3}, 
  {144, 9, 1}, {144, 9, 2}, {144, 9, 3}, 
  {144, 10, 1}, {144, 10, 2}, {144, 10, 3}, 
  {144, 11, 1}, {144, 11, 2}, {144, 11, 3}, 
  {144, 13, 1}, {144, 13, 2}, {144, 13, 3}, 
  {144, 15, 1}, {144, 15, 2}, {144, 15, 3}, 
  {145, 12, 1}, {145, 12, 2}, {145, 12, 3}, 
  {145, 14, 1}, {145, 14, 2}, {145, 14, 3}, 
  {146, 2, 1}, {146, 2, 2}, {146, 2, 3}, 
  {146, 3, 1}, {146, 3, 2}, {146, 3, 3}, 
  {146, 4, 1}, {146, 4, 2}, {146, 4, 3}, 
  {147, 2, 1}, {147, 2, 2}, {147, 2, 3}, 
  {147, 4, 1}, {147, 4, 2}, {147, 4, 3}, 
  {147, 5, 1}, {147, 5, 2}, {147, 5, 3}, 
  {147, 9, 1}, {147, 9, 2}, {147, 9, 3}, 
  {147, 10, 1}, {147, 10, 2}, {147, 10, 3}, 
  {147, 11, 1}, {147, 11, 2}, {147, 11, 3}, 
  {147, 13, 1}, {147, 13, 2}, {147, 13, 3}, 
  {147, 15, 1}, {147, 15, 2}, {147, 15, 3}, 
  {148, 3, 1}, {148, 3, 2}, {148, 3, 3}, 
  {148, 5, 1}, {148, 5, 2}, {148, 5, 3}, 
  {148, 9, 1}, {148, 9, 2}, {148, 9, 3}, 
  {148, 10, 1}, {148, 10, 2}, {148, 10, 3}, 
  {148, 11, 1}, {148, 11, 2}, {148, 11, 3}, 
  {148, 12, 1}, {148, 12, 2}, {148, 12, 3}, 
  {148, 13, 1}, {148, 13, 2}, {148, 13, 3}, 
  {148, 14, 1}, {148, 14, 2}, {148, 14, 3}, 
  {148, 15, 1}, {148, 15, 2}, {148, 15, 3}, 
  {149, 2, 1}, {149, 2, 2}, {149, 2, 3}, 
  {149, 3, 1}, {149, 3, 2}, {149, 3, 3}, 
  {149, 4, 1}, {149, 4, 2}, {149, 4, 3}, 
  {149, 13, 1}, {149, 13, 2}, {149, 13, 3}, 
  {149, 14, 1}, {149, 14, 2}, {149, 14, 3}, 
  {149, 15, 1}, {149, 15, 2}, {149, 15, 3}, 
  {150, 16, 1}, {150, 16, 2}, {150, 16, 3}, 
  {151, 6, 1}, {151, 6, 2}, {151, 6, 3}, 
  {152, 2, 1}, {152, 2, 2}, {152, 2, 3}, 
  {152, 4, 1}, {152, 4, 2}, {152, 4, 3}, 
  {152, 5, 1}, {152, 5, 2}, {152, 5, 3}, 
  {152, 9, 1}, {152, 9, 2}, {152, 9, 3}, 
  {152, 10, 1}, {152, 10, 2}, {152, 10, 3}, 
  {152, 11, 1}, {152, 11, 2}, {152, 11, 3}, 
  {152, 13, 1}, {152, 13, 2}, {152, 13, 3}, 
  {152, 15, 1}, {152, 15, 2}, {152, 15, 3}, 
  {153, 2, 1}, {153, 2, 2}, {153, 2, 3}, 
  {153, 4, 1}, {153, 4, 2}, {153, 4, 3}};

Function(SYSTEM_integer ) GMSGLOBX_solvertypeplatformmapmap(
  SYSTEM_integer i,
  SYSTEM_integer j)
{
  SYSTEM_integer result;

  result = 0;
  if (i >= 1 && i <= 1387) 
    if (j >= 1 && j <= 3) 
      result = GMSGLOBX_solvertypeplatformmaptuple[i - 1][j - 1];
  return result;
}  /* solvertypeplatformmapmap */

Function(SYSTEM_ansichar *) GMSGLOBX_hostplatform(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  switch (P3PLATFORM_osplatform()) {
    case P3PLATFORM_oswindows64emt: 
      _P3strcpy(result,_len_ret,_P3str1("\003WEX"));
      break;
    case P3PLATFORM_oslinux86_64: 
      _P3strcpy(result,_len_ret,_P3str1("\003LEX"));
      break;
    case P3PLATFORM_osdarwin_x64: 
      _P3strcpy(result,_len_ret,_P3str1("\003DEX"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\003XXX"));
  }
  return result;
}  /* hostplatform */

/* unit gmsglobx */
void _Init_Module_gmsglobx(void)
{
} /* _Init_Module_gmsglobx */

void _Final_Module_gmsglobx(void)
{
} /* _Final_Module_gmsglobx */

