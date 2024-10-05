// Model.cpp
// Generated with TDL5 23.10.1.1

#include <stdio.h> // for sprintf and possibly fprintf
static char _msgbuf[1024];

#include <algorithm>
#include <numeric>

// following line stops math.h from defining old bessel functions like y1, etc.
#define _NO_OLDNAMES
#define  NO_OLDNAMES
#include <cmath>
#include <float.h>
#include <string.h>
#include <memory.h>
#include "MUtil.h"
#include "ModelAPI.h"
#include "Model.h"
#undef  _NO_OLDNAMES
#undef   NO_OLDNAMES

using std::min;
using std::max;
using std::isfinite;

extern "C"
{

ModelCallbackTable* _pcb = 0;

int GetNFixef()
{
    return NFIXEF;
}
int GetNSecondary()
{
    return NSECONDARY;
}
int GetNError()
{
    return NERROR;
}
int GetNRanef1()
{
    return NRANEF1;
}
int GetNRanef2()
{
    return NRANEF2;
}
int GetNRanef3()
{
    return NRANEF3;
}
int GetNRanef4()
{
    return NRANEF4;
}
int GetNRanef5()
{
    return NRANEF5;
}
int MultiLevelRanef()
{
    return MULTI_LEVEL_RANEF;
}
int AllowGaussianFit()
{
    return ALLOWGAUSSIANFIT;
}
int AllowLogTransform()
{
    return ALLOWLOGTRANSFORM;
}
int TrulyTimeBased()
{
    return TRULYTIMEBASED;
}
int ModelUniqueId()
{
    return MODELUNIQUEID;
}
int GetNderiv()
{
    return NDERIV;
}
int GetNurine()
{
    return NURINE;
}
int GetNtransit2()
{
    return NTRANSIT2;
}
int GetNtransit()
{
    return NTRANSIT;
}
int GetNevent()
{
    return NEVENT;
}
int GetNintegAll()
{
    return NINTEG_ALL;
}
int GetNcfmacro()
{
    return NCFMACRO;
}

// this allows sigma of current epsilon to be referred to in the model
#define sigma() _std

static double t;
static double WT;
static double _initial_WT;
static double AaDose;
static double AaInfDose;
static double AaInfRate;

bool HasDelay() {
  return false;
}

	// compartment integrators
#define     Aa       zzY[0]
#define zzR_Aa       zzR[0]
#define    _Aa_irate zzIRate[0]

#define     A1       zzY[1]
#define zzR_A1       zzR[1]
#define    _A1_irate zzIRate[1]

 // zzY[2] is extra integrator for infusions

int _i_CMultStdev = -1;	// index of CMultStdev in _fixef
static double ___CMultStdev = 0;		// value of CMultStdev if it is disabled
#define CMultStdev (*(_i_CMultStdev >= 0 ? &_fixef[_i_CMultStdev] : &___CMultStdev))
int _i_tvVmax = -1;	// index of tvVmax in _fixef
static double ___tvVmax = 0;		// value of tvVmax if it is disabled
#define tvVmax (*(_i_tvVmax >= 0 ? &_fixef[_i_tvVmax] : &___tvVmax))
int _i_tvKm = -1;	// index of tvKm in _fixef
static double ___tvKm = 0;		// value of tvKm if it is disabled
#define tvKm (*(_i_tvKm >= 0 ? &_fixef[_i_tvKm] : &___tvKm))
int _i_tvV = -1;	// index of tvV in _fixef
static double ___tvV = 0;		// value of tvV if it is disabled
#define tvV (*(_i_tvV >= 0 ? &_fixef[_i_tvV] : &___tvV))
int _i_tvKa = -1;	// index of tvKa in _fixef
static double ___tvKa = 0;		// value of tvKa if it is disabled
#define tvKa (*(_i_tvKa >= 0 ? &_fixef[_i_tvKa] : &___tvKa))
void InitializeIndirectionTable(){
	_i_CMultStdev = 0;
	_i_tvVmax = 1;
	_i_tvKm = 2;
	_i_tvV = 3;
	_i_tvKa = 4;
}
#define nVmax _ranef5[0]
#define nKm _ranef5[1]
#define nV _ranef5[2]
static double Vmax;
static double _initial_Vmax;
static double Km;
static double _initial_Km;
static double V;
static double _initial_V;
static double Ka;
static double _initial_Ka;
static double C;

void GetFixefName(int *_pi, char* _nm){
    _nm[0] = 0;

    const int fixefCount = 6;

    if (*_pi < 0 || *_pi >= fixefCount)
    {
        return;
    }

    static const char * names[fixefCount] = {
        "CMultStdev",
        "tvVmax",
        "tvKm",
        "tvV",
        "tvKa",
        "CEps"
    };

	strcpy(_nm, names[*_pi]);
}

void GetFixefUnits(int *_pi, char* _nm){
    _nm[0] = 0;

    const int fixefCount = 5;

    if (*_pi < 0 || *_pi >= fixefCount)
    {
        return;
    }

    static const char * names[fixefCount] = {
        "",
        "",
        "",
        "",
        ""
    };

	strcpy(_nm, names[*_pi]);
}

void GetSecondaryName(int *_pi, char* _nm){
	_nm[0] = 0;
	if (0){
	}
}

void InitEnable(){	// Initialize Enablement
	InitializeIndirectionTable();
}

// if gaussian alg is allowed, but not being used, set this to 0
// so that the main standard deviation is treated as a theta
int _bUsingGaussianAlgorithm = 1;
void SetUsingGaussianAlgorithm(int* _pBool){
	_bUsingGaussianAlgorithm = *_pBool;
}

int GetNThetaOld(){	// get number of theta minus possibly 1 sigma to estimate
	int _n = 5 + 1; // number of fixefs + number of epsilons
	// don't count any frozen or disabled thetas
	// if using a gaussian alg, don't treat the last stdev as a theta
	if (_bUsingGaussianAlgorithm)
		_n--;
	return _n;
}

int GetNTheta(){	// get number of theta to estimate
	int _n = 5 + 1; // number of fixefs + number of epsilons
	// don't count any frozen or disabled fixefs
	// don't count any frozen errors
	return _n;
}

int GetNSigma(){	// get number of sigma to estimate, 0 or 1
	int _n = 1; // number of free epsilons
	// if not using a gaussian alg, treat the last stdev as a theta
	if (! _bUsingGaussianAlgorithm)
		_n--;
	return (_n > 0 ? 1 : 0);
}

void GetRanef1Name(int *_pi, char* _nm){
	_nm[0] = 0;
	if (0){
	}
}

void GetRanef2Name(int *_pi, char* _nm){
	_nm[0] = 0;
	if (0);
}

void GetRanef3Name(int *_pi, char* _nm){
	_nm[0] = 0;
	if (0);
}

void GetRanef4Name(int *_pi, char* _nm){
	_nm[0] = 0;
	if (0);
}

void GetRanef5Name(int *_pi, char* _nm){
	_nm[0] = 0;
	if (0);
	else if (*_pi==0) strcpy(_nm, "nVmax");
	else if (*_pi==1) strcpy(_nm, "nKm");
	else if (*_pi==2) strcpy(_nm, "nV");
}

void GetErrorName(int *_pi, char* _nm){
	_nm[0] = 0;
	if (0);
	else if (*_pi==0) strcpy(_nm, "CEps");
}


static double** _args;
static int* _pbOK;

void SetArgsOK(double** __args, int *__pbOK, int *__piState){
	_args = __args;
	_pbOK = __pbOK;
	_piState = __piState;
}

bool IsOK()
{
	return !_pbOK || *_pbOK;
}

void NotOK(){
	if (_pbOK != 0L)(*_pbOK) = 0;
}

void SetCallbacks(ModelCallbackTable* _p){
	_pcb = _p;
}

double unif(){
	double result = 0.5;
	if (_pcb) result = _pcb->GetUniform();
	return result;
}

double norm(){
	double result = 0.0;
	if (_pcb) result = _pcb->GetNormal();
	return result;
}

double valid_strip(double res, const char * message)
{
    if (res == 0)
    {
        SetErrorMsg(message);
        NotOK();
    }

    return res;
}

const char* _t_Units;
const char* _CObs_Units;
const char* _Aa_Units;
const char* _WT_Units;
const char* _AaDose_Units;
const char* _AaInfDose_Units;
const char* _AaInfRate_Units;
void ClearUnits(){
	_t_Units = "";
	_CObs_Units = "";
	_Aa_Units = "";
	_WT_Units = "";
	_AaDose_Units = "";
	_AaInfDose_Units = "";
	_AaInfRate_Units = "";
}

void SetVarUnits(const char* _nm, const char* _units){
    if (_nm[0]=='t' && strcmp(_nm, "t")==0) {
        _t_Units = _units;
        return;
    }
    if (_nm[0]=='C' && strcmp(_nm, "CObs")==0) {
        _CObs_Units = _units;
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "Aa")==0) {
        _Aa_Units = _units;
        return;
    }
    if (_nm[0]=='W' && strcmp(_nm, "WT")==0) {
        _WT_Units = _units;
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaDose")==0) {
        _AaDose_Units = _units;
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaInfDose")==0) {
        _AaInfDose_Units = _units;
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaInfRate")==0) {
        _AaInfRate_Units = _units;
        return;
    }
}

void GetPlotVarInfo(int *_pi, int* _pflags, char* _nmObs, char* _nmPred){
	*_pflags = 0;
	*_nmObs = '\0';
	*_nmPred = '\0';
	if (0){
	} else if (*_pi == 0){
		*_pflags = 0x011f;
		strcpy(_nmObs, "CObs");
		strcpy(_nmPred, "C");
	} else if (*_pi == 1){
		*_pflags = 0x0019;
		strcpy(_nmObs, "");
		strcpy(_nmPred, "Vmax");
	} else if (*_pi == 2){
		*_pflags = 0x0019;
		strcpy(_nmObs, "");
		strcpy(_nmPred, "Km");
	} else if (*_pi == 3){
		*_pflags = 0x0019;
		strcpy(_nmObs, "");
		strcpy(_nmPred, "V");
	} else if (*_pi == 4){
		*_pflags = 0x0019;
		strcpy(_nmObs, "");
		strcpy(_nmPred, "Ka");
	} else if (*_pi == 5){
		*_pflags = 0x0019;
		strcpy(_nmObs, "");
		strcpy(_nmPred, "t");
	} else if (*_pi == 6){
		*_pflags = 0x011b;
		strcpy(_nmObs, "");
		strcpy(_nmPred, "WT");
	}
}

int GetNPlotVar(){return 7;}

void AddToCpt(double zzY[], int _iCpt, double _dv);	// forward declaration

int GetNumCovariates(){
	return 1;
}

int GetCovariateNumber(const char* _nm){
	int _i = -1;
	if (0);
	else if (_nm[0]=='W' && strcmp(_nm, "WT")==0) _i = 0;
	return _i;
}

void GetCovariateName(int* _pi, char* _nm){
	if (0);
	else if (*_pi == 0) strcpy(_nm, "WT");
}

void GetCovariateCategorical(int* _pi, int* _cat){
	if (0);
	else if (*_pi == 0) *_cat = 0;
}

void GetCovariateDirection(int* _pi, int* _dir){
	// 0=backward, 1=forward, 2=interpolate
	if (0);
	else if (*_pi == 0) *_dir = 1;
}

void SetCovariate(double _t, const char* _nm, double _v){
	if (0){
	} else if (_nm[0]=='W' && strcmp(_nm, "WT")==0){
		if (WT == _NA){
			_initial_WT = _v;
		}
		WT = _v;
	}
}

void SetCovariateInterp(double _t, const char* _nm, double _v1, double _dt, double _v2, int _bForward){
	if (0){
	} else if (_nm[0]=='W' && strcmp(_nm, "WT")==0){
		WT = (_bForward ? _v1 : _v2);
	}
}

int GetStripDoseName(char* _nm){
	int _rvalue = 0;
	if (0){
	}
	return _rvalue;
}

int GetNumStructural(){
	return 4;
}

int GetStructuralNumber(const char* _nm){
	int _i = -1;
	if (0);
	else if (_nm[0]=='V' && strcmp(_nm, "Vmax")==0) _i = 0;
	else if (_nm[0]=='K' && strcmp(_nm, "Km")==0) _i = 1;
	else if (_nm[0]=='V' && strcmp(_nm, "V")==0) _i = 2;
	else if (_nm[0]=='K' && strcmp(_nm, "Ka")==0) _i = 3;
	return _i;
}

void GetStructuralName(int* _pi, char* _nm){
	if (0);
	else if (*_pi == 0) strcpy(_nm, "Vmax");
	else if (*_pi == 1) strcpy(_nm, "Km");
	else if (*_pi == 2) strcpy(_nm, "V");
	else if (*_pi == 3) strcpy(_nm, "Ka");
}

void GetVarUnits(const char* _nm, char* _units){
	const char* _p = nullptr;
	_units[0] = 0;
    if (_nm[0]=='A' && strcmp(_nm, "Aa")==0){
        _p = _Aa_Units;
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaDose")==0){
        _p = _Aa_Units;
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaInfDose")==0){
        _p = _Aa_Units;
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaInfRate")==0){
        _p = _units_Div(_Aa_Units,_t_Units);
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='C' && strcmp(_nm, "CObs")==0){
        _p = _CObs_Units;
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='W' && strcmp(_nm, "WT")==0){
        _p = _WT_Units;
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='C' && strcmp(_nm, "C")==0){
        _p = _CObs_Units;
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='K' && strcmp(_nm, "Ka")==0){
        _p = _units_Div(__NumberString(1),_t_Units);
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='K' && strcmp(_nm, "Km")==0){
        _p = _CObs_Units;
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvKm")==0){
        _p = _CObs_Units;
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvKa")==0){
        _p = _units_Div(__NumberString(1),_t_Units);
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "A1")==0){
        _p = _Aa_Units;
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='V' && strcmp(_nm, "V")==0){
        _p = _units_Div(_Aa_Units,_CObs_Units);
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvV")==0){
        _p = _units_Div(_Aa_Units,_CObs_Units);
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='V' && strcmp(_nm, "Vmax")==0){
        _p = _units_Div(_Aa_Units,_t_Units);
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvVmax")==0){
        _p = _units_Div(_Aa_Units,_t_Units);
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='n' && strcmp(_nm, "nV")==0){
        _p = __NumberString(0);
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='n' && strcmp(_nm, "nKm")==0){
        _p = __NumberString(0);
        strcpy(_units, _p);
        return;
    }
    if (_nm[0]=='n' && strcmp(_nm, "nVmax")==0){
        _p = __NumberString(0);
        strcpy(_units, _p);
        return;
    }

    if (_nm[0]=='t' && strcmp(_nm, "t")==0){
		_p = _t_Units;
        strcpy(_units, _p);
        return;
    }
}


// routine to read out all the CF state vector
void GetAllCFState(double* _y){
}
void SetAllCFState(double* _y){
}
int GetNInteg(){
	return NINTEG_ALL;
}
void ClearSSIrrelevantState(double* _y){
}
void GetIntegrators(double* _y, double* zzY){
	int _i;
	_y[0] = Aa;
	_y[1] = A1;
}
void SetIntegrators(double* _y, double* zzY){
	int _i;
	Aa = _y[0];
	A1 = _y[1];
}
void SetTimeOfStateCF(double _t){
}
void SetTimeOfGetValCF(double _t){
}
void ResetSubject(){
	ResetTime();
	ClearActions();
	ResetInteg();
	ResetCF();
	ResetUrine();
	TDL4ResetEvent();
	ResetGradient();
}
void InitDelays(double _t){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
}
void InitCF(double _t){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	#if defined(USE_SETTIMECF)
		SetTimeOfStateCF(_t);
		SetTimeOfGetValCF(_t);
	#endif
}
// routine to capture CF state, change covariates, and re-initialize CF machines
// make sure all covariate-derived values are up to date before calling
void ReInitCF(double _t){
	double _state[NCFMACRO*6+1];
	GetAllCFState(_state); // get the compartment values
	InitCF(_t); // regenerate all the As
	SetAllCFState(_state); // put the values back in all the compartments
}
void GetNumCVMean(int* _pn){*_pn = 0;}
void GetCVMeanName(int* _pi, char* _nm){
	if (0);
}
void GetCVMean(int* _pi, double* _pv){
	if (0);
}
void SetCVMean(int* _pi, double* _pv){
	if (0);
}
void GetNumCVMedian(int* _pn){*_pn = 0;}
void GetCVMedianName(int* _pi, char* _nm){
	if (0);
}
void GetCVMedian(int* _pi, double* _pv){
	if (0);
}
void SetCVMedian(int* _pi, double* _pv){
	if (0);
}
void InitSubject(){
	ResetSubject();
	WT = _NA;
	_initial_WT = _NA;
	AaDose = 0;
	AaInfDose = 0;
	AaInfRate = 0;
	Vmax = _NA;
	_initial_Vmax = _NA;
	Km = _NA;
	_initial_Km = _NA;
	V = _NA;
	_initial_V = _NA;
	Ka = _NA;
	_initial_Ka = _NA;
	ClearErrorMsg();

	*_piState = ODESTATE_START;
}

void RestartSubjectSequences(double* zzY){
}

void RestartSubjectEventSims(double* zzY){
}
void EndSubject(int bSimulating){
	if (bSimulating){
	}
	// set covariates to their first settings, in case needed for secondary parameters
	SetCovariate(curTime, "WT", _initial_WT);
	Vmax = _initial_Vmax;
	Km = _initial_Km;
	V = _initial_V;
	Ka = _initial_Ka;
}
static double __CMultStdev_low  = 0.0;
static double __CMultStdev_init = 1.0;
static double __CMultStdev_high = 3.3e-151;
static double __tvVmax_low  = 100.0;
static double __tvVmax_init = 4000.0;
static double __tvVmax_high = 3.3e-151;
static double __tvKm_low  = 100.0;
static double __tvKm_init = 4000.0;
static double __tvKm_high = 3.3e-151;
static double __tvV_low  = 1.0;
static double __tvV_init = 50.0;
static double __tvV_high = 3.3e-151;
static double __tvKa_low  = 0.01000000000000000021;
static double __tvKa_init = 1.199999999999999956;
static double __tvKa_high = 3.0;
void RestoreFixefLowInitHighFromModel(){
	__CMultStdev_low  = 0.0;
	__CMultStdev_init = 1.0;
	__CMultStdev_high = 3.3e-151;
	__tvVmax_low  = 100.0;
	__tvVmax_init = 4000.0;
	__tvVmax_high = 3.3e-151;
	__tvKm_low  = 100.0;
	__tvKm_init = 4000.0;
	__tvKm_high = 3.3e-151;
	__tvV_low  = 1.0;
	__tvV_init = 50.0;
	__tvV_high = 3.3e-151;
	__tvKa_low  = 0.01000000000000000021;
	__tvKa_init = 1.199999999999999956;
	__tvKa_high = 3.0;
}
void SetFixefLowInitHigh(const char* _nm, double* _fixefLow, double* _fixefInit, double* _fixefHigh){
    if (strcmp(_nm, "CMultStdev")==0){
        __CMultStdev_low  = *_fixefLow;
        __CMultStdev_init = *_fixefInit;
        __CMultStdev_high = *_fixefHigh;
        return;
	}
    if (strcmp(_nm, "tvVmax")==0){
        __tvVmax_low  = *_fixefLow;
        __tvVmax_init = *_fixefInit;
        __tvVmax_high = *_fixefHigh;
        return;
	}
    if (strcmp(_nm, "tvKm")==0){
        __tvKm_low  = *_fixefLow;
        __tvKm_init = *_fixefInit;
        __tvKm_high = *_fixefHigh;
        return;
	}
    if (strcmp(_nm, "tvV")==0){
        __tvV_low  = *_fixefLow;
        __tvV_init = *_fixefInit;
        __tvV_high = *_fixefHigh;
        return;
	}
    if (strcmp(_nm, "tvKa")==0){
        __tvKa_low  = *_fixefLow;
        __tvKa_init = *_fixefInit;
        __tvKa_high = *_fixefHigh;
        return;
	}
}
void GetFixefLowInitHigh(const char* _nm, double* _fixefLow, double* _fixefInit, double* _fixefHigh){
	if (0){
    }
    if (strcmp(_nm, "CMultStdev")==0){
        *_fixefLow  = __CMultStdev_low;
        *_fixefInit = __CMultStdev_init;
        *_fixefHigh = __CMultStdev_high;
        return;
    }
    if (strcmp(_nm, "tvVmax")==0){
        *_fixefLow  = __tvVmax_low;
        *_fixefInit = __tvVmax_init;
        *_fixefHigh = __tvVmax_high;
        return;
    }
    if (strcmp(_nm, "tvKm")==0){
        *_fixefLow  = __tvKm_low;
        *_fixefInit = __tvKm_init;
        *_fixefHigh = __tvKm_high;
        return;
    }
    if (strcmp(_nm, "tvV")==0){
        *_fixefLow  = __tvV_low;
        *_fixefInit = __tvV_init;
        *_fixefHigh = __tvV_high;
        return;
    }
    if (strcmp(_nm, "tvKa")==0){
        *_fixefLow  = __tvKa_low;
        *_fixefInit = __tvKa_init;
        *_fixefHigh = __tvKa_high;
        return;
	}
}
void GetThetaLimits(double* _fixefLow, double* _fixefHigh){
	int _n = 0;
	_fixefLow[_n] = __CMultStdev_low; // set CMultStdev
	_fixefHigh[_n] = __CMultStdev_high;
	_n++;
	_fixefLow[_n] = __tvVmax_low; // set tvVmax
	_fixefHigh[_n] = __tvVmax_high;
	_n++;
	_fixefLow[_n] = __tvKm_low; // set tvKm
	_fixefHigh[_n] = __tvKm_high;
	_n++;
	_fixefLow[_n] = __tvV_low; // set tvV
	_fixefHigh[_n] = __tvV_high;
	_n++;
	_fixefLow[_n] = __tvKa_low; // set tvKa
	_fixefHigh[_n] = __tvKa_high;
	_n++;
	// if not using a gaussian alg, last theta is the main stdev
	if (! _bUsingGaussianAlgorithm){
		_fixefLow[_n] = 0; // stdev 0
		_fixefHigh[_n] = 1e100;
		_n++;
	}
}

void InitFixefErrorEstimate(){
	double* _fixef = _args[0];
	_fixef[0] = __CMultStdev_init; // set CMultStdev
	_fixef[1] = __tvVmax_init; // set tvVmax
	_fixef[2] = __tvKm_init; // set tvKm
	_fixef[3] = __tvV_init; // set tvV
	_fixef[4] = __tvKa_init; // set tvKa
	_fixef[5] = 0.5; // stdev 0
}

void GetFixefErrorInitial(double* _fixef){
	_fixef[0] = __CMultStdev_init; // set CMultStdev
	_fixef[1] = __tvVmax_init; // set tvVmax
	_fixef[2] = __tvKm_init; // set tvKm
	_fixef[3] = __tvV_init; // set tvV
	_fixef[4] = __tvKa_init; // set tvKa
	_fixef[5] = 0.5; // stdev 0
}

void GetFixefErrorLimits(double* _low, double* _high){
	_low[0] = __CMultStdev_low; // set CMultStdev
	_high[0] = __CMultStdev_high;
	_low[1] = __tvVmax_low; // set tvVmax
	_high[1] = __tvVmax_high;
	_low[2] = __tvKm_low; // set tvKm
	_high[2] = __tvKm_high;
	_low[3] = __tvV_low; // set tvV
	_high[3] = __tvV_high;
	_low[4] = __tvKa_low; // set tvKa
	_high[4] = __tvKa_high;
	_low[5] = 0; // stdev 0
	_high[5] = 1e100; // stdev 0
}

void GetFixefValue(int* _iFixef, double* _dValue){
	*_dValue = _args[0][*_iFixef];
	return;
}

double GetObsSigma(int* _piWhichObs)
{
	double* _fixef = _args[0];
	double _std = 0;
	if (0){
	} else if (*_piWhichObs == 0){
		_std = _fixef[5];
	}
	return _std;
}

const char * GetObsSigmaName(int* _piWhichObs)
{
	if (0){
	} else if (*_piWhichObs == 0){
		return "CEps";
	}
	return "";
}

void GetInitialStdev(double* _svec){
	int _i = 0;
	double _temp = 0;
	_svec[0] = 0.5;
}

void GetInitialOmega5(double* _omat){
	int _i = 0;
	double _temp = 0;
	for (_i=0; _i<9; _i++) _omat[_i] = 0;
	_omat[0+0+0*3] = 0.1000000000000000056;
	_omat[0+1+1*3] = 0.1000000000000000056;
	_omat[0+2+2*3] = 0.1000000000000000056;
}

double GetEtaVariance5(int _i){
	double _temp = 1;
	if (0){
	} else if (_i == 0){
		_temp = 0.1000000000000000056;
	} else if (_i == 1){
		_temp = 0.1000000000000000056;
	} else if (_i == 2){
		_temp = 0.1000000000000000056;
	}
	return _temp;
}

void TruncateOmega5(double* _omat){
	int _i = 0;
	double _temp = 0;
	double _mtemp[NRANEF5*NRANEF5+1];
	memcpy(_mtemp, _omat, sizeof(double)*NRANEF5*NRANEF5);
	memset(_omat, 0, sizeof(double)*NRANEF5*NRANEF5);
	_temp = _mtemp[0+0+0*3];
	_omat[0+0+0*3] = _temp;
	_temp = _mtemp[0+1+1*3];
	_omat[0+1+1*3] = _temp;
	_temp = _mtemp[0+2+2*3];
	_omat[0+2+2*3] = _temp;
}

// Caution: omega matrix is assumed to be in lower cholesky form
void UnpackParameters(double* _thvec, double* _omat, double* _pack, int* _ldOmat, int* _nPack){
	int _i = 0, _j, _n = 0, _i1 = 0, _j1 = 0, _ldo = *_ldOmat;
	double _temp = 0;
	// assert NRANEF5 <= ldo
	_i1 = GetNTheta();
	for (_i = 0; _i < _i1; _i++){
		_thvec[_n++] = _pack[_i];
	}
	_i = 0 + 0;
	_j = 0 + 0;
	_omat[_j * _ldo + _i] = _pack[_n++];
	_i = 0 + 1;
	_j = 0 + 1;
	_omat[_j * _ldo + _i] = _pack[_n++];
	_i = 0 + 2;
	_j = 0 + 2;
	_omat[_j * _ldo + _i] = _pack[_n++];
	// assert *_nPack == _n;
}

// Caution: omega matrix is assumed to be in lower cholesky form
void UnpackParameters1(double* _thvec, double* _omat, double* _svec, double* _pack, int* _ldOmat, int* _nPack){
	int _i = 0, _j, _n = 0, _i1 = 0, _j1 = 0, _ldo = *_ldOmat;
	double _temp = 0;
	// assert NRANEF5 <= ldo
	_i1 = GetNTheta();
	for (_i = 0; _i < _i1; _i++){
		_thvec[_n++] = _pack[_i];
	}
	if (_i1 > 0){
		_svec[0] = _pack[_i1-1];
	}
	_i = 0 + 0;
	_j = 0 + 0;
	_omat[_j * _ldo + _i] = _pack[_n++];
	_i = 0 + 1;
	_j = 0 + 1;
	_omat[_j * _ldo + _i] = _pack[_n++];
	_i = 0 + 2;
	_j = 0 + 2;
	_omat[_j * _ldo + _i] = _pack[_n++];
	// assert *_nPack == _n;
}

// Caution: omega matrix is assumed to be in lower cholesky form
void   PackParameters(double* _thvec, double* _omat, double* _pack, int* _ldOmat, int* _nPack){
	int _i = 0, _j, _n = 0, _i1 = 0, _ldo = *_ldOmat;
	// assert NTHETA == *nth
	_i1 = GetNTheta();
	for (_i = 0; _i < _i1; _i++){
		_pack[_i] = _thvec[_n++];
	}
	{
		_i = 0 + 0;
		_j = 0 + 0;
		_pack[_n++] = _omat[_j * _ldo + _i];
	}
	{
		_i = 0 + 1;
		_j = 0 + 1;
		_pack[_n++] = _omat[_j * _ldo + _i];
	}
	{
		_i = 0 + 2;
		_j = 0 + 2;
		_pack[_n++] = _omat[_j * _ldo + _i];
	}
	*_nPack = _n;
}

// Caution: omega matrix is assumed to be in lower cholesky form
void   PackParameters1(double* _thvec, double* _omat, double* _svec, double* _pack, int* _ldOmat, int* _nPack){
	int _i = 0, _j, _n = 0, _i1 = 0, _ldo = *_ldOmat;
	// assert NTHETA == *nth
	_i1 = GetNTheta();
	for (_i = 0; _i < _i1; _i++){
		_pack[_i] = _thvec[_n++];
	}
	if (_i1 > 0){
		_pack[_i1-1] = _svec[0];
	}
	{
		_i = 0 + 0;
		_j = 0 + 0;
		_pack[_n++] = _omat[_j * _ldo + _i];
	}
	{
		_i = 0 + 1;
		_j = 0 + 1;
		_pack[_n++] = _omat[_j * _ldo + _i];
	}
	{
		_i = 0 + 2;
		_j = 0 + 2;
		_pack[_n++] = _omat[_j * _ldo + _i];
	}
	*_nPack = _n;
}

// Caution: omega matrix is assumed to be in lower cholesky form
void   TypXParameters(double* _thvec, double* _omat, double* _pack, int* _ldOmat, int* _nPack){
	int _i = 0, _j, _n = 0, _i1 = 0, _ldo = *_ldOmat;
	// assert NTHETA == *nth
	_i1 = GetNTheta();
	for (_i = 0; _i < _i1; _i++){
		_pack[_i] = fabs(_thvec[_n++]);
	}
	{
		_i = 0 + 0;
		_j = 0 + 0;
		_pack[_n++] = fabs(_omat[_j * _ldo + _i]);
	}
	{
		_i = 0 + 1;
		_j = 0 + 1;
		_pack[_n++] = fabs(_omat[_j * _ldo + _i]);
	}
	{
		_i = 0 + 2;
		_j = 0 + 2;
		_pack[_n++] = fabs(_omat[_j * _ldo + _i]);
	}
	*_nPack = _n;
}

// Caution: omega matrix is assumed to be in lower cholesky form
void   TypXParameters1(double* _thvec, double* _omat, double* _svec, double* _pack, int* _ldOmat, int* _nPack){
	int _i = 0, _j, _n = 0, _i1 = 0, _ldo = *_ldOmat;
	// assert NTHETA == *nth
	_i1 = GetNTheta();
	for (_i = 0; _i < _i1; _i++){
		_pack[_i] = fabs(_thvec[_n++]);
	}
	if (_i1 > 0){
		_pack[_i1-1] = fabs(_svec[0]);
	}
	{
		_i = 0 + 0;
		_j = 0 + 0;
		_pack[_n++] = fabs(_omat[_j * _ldo + _i]);
	}
	{
		_i = 0 + 1;
		_j = 0 + 1;
		_pack[_n++] = fabs(_omat[_j * _ldo + _i]);
	}
	{
		_i = 0 + 2;
		_j = 0 + 2;
		_pack[_n++] = fabs(_omat[_j * _ldo + _i]);
	}
	*_nPack = _n;
}

// Caution: omega matrix is assumed to be in lower cholesky form
void UnpackOmega(double* _omat, double* _pack, int* _ldOmat, int* _nPack){
	int _i = 0, _j, _n = 0, _i1 = 0, _j1 = 0, _ldo = *_ldOmat;
	double _temp = 0;
	// assert NRANEF5 <= ldo
	_i = 0 + 0;
	_j = 0 + 0;
	_omat[_j * _ldo + _i] = _pack[_n++];
	_i = 0 + 1;
	_j = 0 + 1;
	_omat[_j * _ldo + _i] = _pack[_n++];
	_i = 0 + 2;
	_j = 0 + 2;
	_omat[_j * _ldo + _i] = _pack[_n++];
	// assert *_nPack == _n;
}

// Caution: omega matrix is assumed to be in lower cholesky form
void   PackOmega(double* _omat, double* _pack, int* _ldOmat, int* _nPack){
	int _i = 0, _j, _n = 0, _ldo = *_ldOmat;
	// assert NTHETA == *nth
	{
		_i = 0 + 0;
		_j = 0 + 0;
		_pack[_n++] = _omat[_j * _ldo + _i];
	}
	{
		_i = 0 + 1;
		_j = 0 + 1;
		_pack[_n++] = _omat[_j * _ldo + _i];
	}
	{
		_i = 0 + 2;
		_j = 0 + 2;
		_pack[_n++] = _omat[_j * _ldo + _i];
	}
	*_nPack = _n;
}

// Caution: omega matrix is assumed to be in lower cholesky form
void   TypXOmega(double* _omat, double* _pack, int* _ldOmat, int* _nPack){
	int _i = 0, _j, _n = 0, _ldo = *_ldOmat;
	// assert NTHETA == *nth
	{
		_i = 0 + 0;
		_j = 0 + 0;
		_pack[_n++] = fabs(_omat[_j * _ldo + _i]);
	}
	{
		_i = 0 + 1;
		_j = 0 + 1;
		_pack[_n++] = fabs(_omat[_j * _ldo + _i]);
	}
	{
		_i = 0 + 2;
		_j = 0 + 2;
		_pack[_n++] = fabs(_omat[_j * _ldo + _i]);
	}
	*_nPack = _n;
}

void GetFixefErrToThetaMap(int* _iaTheta, int* _nth){
	int _n = 0;
	_iaTheta[0] = _n++; // CMultStdev
	_iaTheta[1] = _n++; // tvVmax
	_iaTheta[2] = _n++; // tvKm
	_iaTheta[3] = _n++; // tvV
	_iaTheta[4] = _n++; // tvKa
	// unfrozen error sigmas 
	_iaTheta[5] = _n++; // stdev 0
	*_nth = _n;
}

void GetRanefToEtaMap(int* _iaEta, int* _nth){
	int _n = 0;
	_iaEta[0] = _n++; // ranef 0
	_iaEta[1] = _n++; // ranef 1
	_iaEta[2] = _n++; // ranef 2
	*_nth = _n;
}

void UnpackTheta(double* _fixef, double* _thvec){
	int _n = 0;
	// copy theta to fixef+error
	// unfrozen, enabled fixed effects 
	_fixef[0] = _thvec[_n++];
	_fixef[1] = _thvec[_n++];
	_fixef[2] = _thvec[_n++];
	_fixef[3] = _thvec[_n++];
	_fixef[4] = _thvec[_n++];
	// unfrozen error sigmas 
	_fixef[5] = _thvec[_n++];
}

void   PackTheta(double* _fixef, double* _thvec, int* _nth){
	int _n = 0;
	// copy fixef+error to theta
	// unfrozen, enabled fixed effects 
	_thvec[_n++] = _fixef[0];
	_thvec[_n++] = _fixef[1];
	_thvec[_n++] = _fixef[2];
	_thvec[_n++] = _fixef[3];
	_thvec[_n++] = _fixef[4];
	// unfrozen error sigmas 
	_thvec[_n++] = _fixef[5];
	*_nth = _n;
}

void   TypXTheta(double* _fixef, double* _thvec, int* _nth){
	int _n = 0;
	// copy fabs(fixef+error) to theta
	// unfrozen, enabled fixed effects 
	_thvec[_n++] = fabs(_fixef[0]);
	_thvec[_n++] = fabs(_fixef[1]);
	_thvec[_n++] = fabs(_fixef[2]);
	_thvec[_n++] = fabs(_fixef[3]);
	_thvec[_n++] = fabs(_fixef[4]);
	// unfrozen error sigmas 
	_thvec[_n++] = fabs(_fixef[5]);
	*_nth = _n;
}

// Caution: live portion of omega matrix is assumed to be packed in lower cholesky form
void   GetNPack(int* _nPack){
	int _n = 0;
	_n++; // CMultStdev
	_n++; // tvVmax
	_n++; // tvKm
	_n++; // tvV
	_n++; // tvKa
	_n++; // CEps
	{
		_n++;
	}
	{
		_n++;
	}
	{
		_n++;
	}
	*_nPack = _n;
}

// Caution: live portion of omega matrix is assumed to be packed in lower cholesky form
void   GetNPackOmega(int* _nPack){
	int _n = 0;
	{
		_n++;
	}
	{
		_n++;
	}
	{
		_n++;
	}
	*_nPack = _n;
}

void   GetNSeqOmega(int* _nSeq){
	int _n = 0;
	_n++; // block 0
	*_nSeq = _n;
}

void   GetOmegaSequenceInfo(int* _iSeq1, int* _nEta, int* _nSame){
	int _iSeqCount = 0;
	*_nEta = *_nSame = 0;
	if ((*_iSeq1-1) == _iSeqCount){
		*_nEta = 3;
		*_nSame = 0;
		return;
	}
	_iSeqCount++;
}

// Caution: live portion of omega matrix is assumed to be packed in lower cholesky form
void   GetNPackRanCov(int* _nPack){
	int _n = 0;
	{
		_n++;
	}
	{
		_n++;
	}
	{
		_n++;
	}
	*_nPack = _n;
}

// Caution: live portion of omega matrix is assumed to be packed in lower cholesky form
void   GetPackedIndexRanCov(int* _iRan, int* _jRan, int* _index){
	int _n = 0;
	*_index = -1;
	_n++; // CMultStdev
	_n++; // tvVmax
	_n++; // tvKm
	_n++; // tvV
	_n++; // tvKa
	_n++; // CEps
	{
		if (*_iRan==0 && *_jRan==0){*_index = _n; return;}
		_n++;
	}
	{
		if (*_iRan==1 && *_jRan==1){*_index = _n; return;}
		_n++;
	}
	{
		if (*_iRan==2 && *_jRan==2){*_index = _n; return;}
		_n++;
	}
}

static double _dMainSigma;
void SetMainSigma(double* _sig){
	_dMainSigma = *_sig;
}

static double _tempResult[3];
static double _tempArg[3];
int _bTempValid = 0;
void InitTemp(){ int _i;
	for (_i = 0; _i < 3; _i++) _tempResult[_i] = -999999999.0;
	for (_i = 0; _i < 3; _i++) _tempArg[_i]    = -999999999.0;
	_bTempValid = 1;
}
void EvalSParms(double* zzY){
	if (zzY == 0)
    {
        return;
    }
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	int _ix0 = 0, _ix1 = 0;	// used for table statement
	double _xFrac = 0;
	if (!_bTempValid) InitTemp();
	EvalGroup(zzY);
	if (WT == _NA){
		SetErrorMsg("Covariate 'WT' not set");
		NotOK();
	}
	Vmax = (tvVmax * _exp_cached(nVmax, &_tempArg[0], &_tempResult[0]));
	if (_initial_Vmax == _NA) _initial_Vmax = Vmax;
	Km = (tvKm * _exp_cached(nV, &_tempArg[1], &_tempResult[1]));
	if (_initial_Km == _NA) _initial_Km = Km;
	V = (tvV * _exp_cached(nKm, &_tempArg[2], &_tempResult[2]));
	if (_initial_V == _NA) _initial_V = V;
	Ka = tvKa;
	if (_initial_Ka == _NA) _initial_Ka = Ka;
	EvalAssign(zzY);
}
void EvalGroup(double* zzY){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	int _ix0 = 0, _ix1 = 0;	// used for table statement
	double _xFrac = 0;
	if (!_bTempValid) InitTemp();
}

static double calc_sum_x(int n_, double shape_, double x_[])
{
    double BinomialCoeff_ = 1;

    double res_ = x_[0];

    for (int i_ = 1; i_ < n_; ++i_)
    {
        BinomialCoeff_ = -BinomialCoeff_ * (shape_ - i_) / i_;
        res_ += BinomialCoeff_ * x_[i_];
    }

    return res_;
}

static double _lgamm_cached(double x)
{
    static double _cachedX = 1;
    static double _cachedGamma = lgamm(_cachedX);

    if (x != _cachedX)
    {
        _cachedX = x;
        _cachedGamma = lgamm(x);
    }

    return _cachedGamma;
}

void EvalAssign(double* zzY){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	int _ix0 = 0, _ix1 = 0;	// used for table statement
	double _xFrac = 0;
	if (!_bTempValid) InitTemp();
	C = (A1 / V);
}
void CalcSecondary(){
	double* _fixef = _args[0];
};

void GetSecondary(int* _pn, double* _sec){
	CalcSecondary();
	*_pn = 0;
};

void SetModelTime(double _t){
	t = _t;
}

static int GetGammaDelayDECnt(int i_)
{
    return 0;
}

int GetGammaDelayDEStart(int i_)
{
    return NINTEG_ALL + 1 + (i_ > 0 ? GetGammaDelayDECnt(i_ - 1) : 0);
}

int GetGammaDelayDECount()
{
    return 0;
}


static void gamma_deriv(int n_, double s_, double eps_, double ktr_, double x_[], double R_[], double irate_[])
{
    for (int i_ = 0; i_ < n_; ++i_)
    {
        R_[i_] = (s_ - (ktr_ + (i_ / (t + eps_))) * x_[i_]) + irate_[i_];
    }
}

void Deriv(long* _pN, double* _pt, double zzY[], double zzR[], double zzIRate[]){
	int _i = 0;
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	t = *_pt;
	#if defined(USE_SETTIMECF)
		SetTimeOfGetValCF(*_pt);
	#endif
	#if defined(THERE_ARE_INTERPOLATED_COVARIATES)
		EvalSParms(zzY);
	#else
		EvalAssign(zzY);
	#endif
	zzR_Aa = ( -(Ka) * Aa);
	zzR_Aa += _Aa_irate;
	zzR_A1 = ((Ka * Aa) - ((Vmax * C) / (Km + C)));
	zzR_A1 += _A1_irate;
	zzR[2] = 0; // extra null equation
}
	//  -(tvKa)
	// tvKa
	// 0.0
	//  -((((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) - (((tvVmax * exp(nVmax)) * (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))) / (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))))))  is not constant
	// 0.0
	// 0.0
void Jacobian(int* pN, double* pt, double zzY[], double zzIRate[], int* pML, int* pMU, double* _jacobian, int* pNRJ){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	int _nrj = *pNRJ;
	double __temp[3];
	__temp[0] = exp(nV);
	__temp[1] = exp(nKm);
	__temp[2] = exp(nVmax);
	memset(_jacobian, 0, sizeof(double)*_nrj*_nrj);
	// rate terms for Aa
		// w.r.t. Aa
	_jacobian[0 + _nrj*0] =  -(tvKa);
		// forcing rate for Aa
	_jacobian[0 + _nrj*2] = 0.0 + zzIRate[0];
	// rate terms for A1
		// w.r.t. Aa
	_jacobian[1 + _nrj*0] = tvKa;
		// w.r.t. A1
	_jacobian[1 + _nrj*1] =  -((((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) - (((tvVmax * __temp[2]) * (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))) / (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))))));
		// forcing rate for A1
	_jacobian[1 + _nrj*2] = 0.0 + zzIRate[1];
}
void ResetTime(){
	curTime = 0;
	curTimeLastDose = curTime;
	iCurWhichDose = 0;
	t = curTime;
}
void ResetCF(){
}
void Advance(int bGenGradients, int iODELevel, double *pTimeOuter, double _dt, int bTimeAdvances, double zzY[]){
	#if defined(USE_SETTIMECF)
		SetTimeOfStateCF(*pTimeOuter);
	#endif
	// truncate the level
	if (iODELevel > ODELEVEL_DOPRI) iODELevel = ODELEVEL_MATEXP;
	if (iODELevel < 1) iODELevel = 3;
	// don't have constant analytic jacobian: matexp defaults to RK
	if (iODELevel == ODELEVEL_MATEXP) iODELevel = ODELEVEL_NONSTIFF;
	ODESolve(bGenGradients, iODELevel, pTimeOuter, _dt, bTimeAdvances, zzY);
	EvalSParms(zzY);
}
int GetNumCompartments(){
	return 2;
}
int GetCompartmentNumber(const char* _nm){
	int _i = -1;
	if (0);
	else if (_nm[0]=='A' && strcmp(_nm, "Aa")==0) _i = 0;
	else if (_nm[0]=='A' && strcmp(_nm, "A1")==0) _i = 1;
	return _i;
}
void GetCompartmentName(int* _pi, char* _nm){
	_nm[0] = 0;
	if (0);
	else if (*_pi == 0) strcpy(_nm, "Aa");
	else if (*_pi == 1) strcpy(_nm, "A1");
}
void ScheduleBolus1(int _mark, double _tRel, int _iCpt, double _amt, int _addl){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	double _tlag = 0;
	double _rate = 1, _dt = 1, _bio = 1;
	if (0){
	} else if (_iCpt==0){
		if (AaDose == 0) AaDose = _amt;

		           InsertAction(_tRel        , 0, ZZACTION_DOSETIME, 0, 0, 0, (void*)0L);		// mark dose requested time
		           InsertAction(_tRel + _tlag, 0, ZZACTION_BOLUS1, _iCpt, _amt * _bio, _NA, (void*)0L);		// bolus
		if (_mark) InsertAction(_tRel        , 0,   ZZMARK_BOLUS1, _iCpt, _amt, 0, (void*)0L);		// mark bolus requested time
	}
}
void ScheduleIV1(int _mark, double _tRel, int _iCpt, double _amt, double _rate, int _addl){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	double _tlag = 0, _bio = 1, _dt;
	if (_rate == 0) return;
	_dt = _amt / _rate;
	if (0){
	} else if (_iCpt==0){
		if (AaInfDose == 0 && AaInfRate == 0){AaInfDose = _amt; AaInfRate = _rate;}

		           InsertAction(_tRel              , 0, ZZACTION_DOSETIME, 0, 0, 0, (void*)0L);		// mark dose time
		           InsertAction(_tRel + _tlag      , 0, ZZACTION_INFSTRT1, _iCpt, _rate * _bio, 0, (void*)0L);		// start/end inf
		if (_mark) InsertAction(_tRel              , 0,   ZZMARK_INFSTRT1, _iCpt, _amt, _rate, (void*)0L);		// start/end inf
		           InsertAction(_tRel + _tlag + _dt, 0, ZZACTION_INFEND1, _iCpt, _rate * _bio, 0, (void*)0L);
	}
}
void PerformBolus1(double* _pt, double zzY[], double zzR[], int _iCpt, double _dv, double _rate){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	} else if (_iCpt==0){
		zzY[_iCpt] += _dv;
	}
}
void PerformInfStrt1(double* _pt, double zzY[], double zzR[], double zzIRate[], int _iCpt, double _dv){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	} else if (_iCpt==0){
		zzIRate[_iCpt] += _dv;
	}
}
void PerformInfEnd1(double* _pt, double zzY[], double zzR[], double zzIRate[], int _iCpt, double _dv){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	} else if (_iCpt==0){
		zzIRate[_iCpt] += - _dv;
	}
}
void ScheduleBolus2(int _mark, double _tRel, int _iCpt, double _amt, int _addl){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	double _tlag = 0;
	double _rate = 1, _dt = 1, _bio = 1;
	if (0){
	}
}
void ScheduleIV2(int _mark, double _tRel, int _iCpt, double _amt, double _rate, int _addl){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	double _tlag = 0, _bio = 1, _dt;
	if (_rate == 0) return;
	_dt = _amt / _rate;
	if (0){
	}
}
void PerformBolus2(double* _pt, double zzY[], double zzR[], int _iCpt, double _dv, double _rate){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	}
}
void PerformInfStrt2(double* _pt, double zzY[], double zzR[], double zzIRate[], int _iCpt, double _dv){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	}
}
void PerformInfEnd2(double* _pt, double zzY[], double zzR[], double zzIRate[], int _iCpt, double _dv){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	}
}
void AddToCpt(double zzY[], int _iCpt, double _dv){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	} else if (_iCpt==0){
		zzY[_iCpt] += _dv;
	} else if (_iCpt==1){
		zzY[_iCpt] += _dv;
	}
}
int GetVarValue(const char* _nm, double* _pv, double zzY[]){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	EvalSParms(zzY);

    if (_nm[0]=='t' && strcmp(_nm, "t")==0) {
        *_pv = t;
        return 1;
    }
    if (_nm[0]=='W' && strcmp(_nm, "WT")==0) {
        *_pv = WT;
        return 1;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaDose")==0) {
        *_pv = AaDose;
        return 1;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaInfDose")==0) {
        *_pv = AaInfDose;
        return 1;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaInfRate")==0) {
        *_pv = AaInfRate;
        return 1;
    }
    if (_nm[0]=='A' && strcmp(_nm, "Aa")==0) {
        *_pv = Aa;
        return 1;
    }
    if (_nm[0]=='A' && strcmp(_nm, "A1")==0) {
        *_pv = A1;
        return 1;
    }
    if (_nm[0]=='C' && strcmp(_nm, "CMultStdev")==0) {
        *_pv = CMultStdev;
        return 1;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvVmax")==0) {
        *_pv = tvVmax;
        return 1;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvKm")==0) {
        *_pv = tvKm;
        return 1;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvV")==0) {
        *_pv = tvV;
        return 1;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvKa")==0) {
        *_pv = tvKa;
        return 1;
    }
    if (_nm[0]=='n' && strcmp(_nm, "nVmax")==0) {
        *_pv = nVmax;
        return 1;
    }
    if (_nm[0]=='n' && strcmp(_nm, "nKm")==0) {
        *_pv = nKm;
        return 1;
    }
    if (_nm[0]=='n' && strcmp(_nm, "nV")==0) {
        *_pv = nV;
        return 1;
    }
    if (_nm[0]=='V' && strcmp(_nm, "Vmax")==0) {
        *_pv = Vmax;
        return 1;
    }
    if (_nm[0]=='K' && strcmp(_nm, "Km")==0) {
        *_pv = Km;
        return 1;
    }
    if (_nm[0]=='V' && strcmp(_nm, "V")==0) {
        *_pv = V;
        return 1;
    }
    if (_nm[0]=='K' && strcmp(_nm, "Ka")==0) {
        *_pv = Ka;
        return 1;
    }
    if (_nm[0]=='C' && strcmp(_nm, "C")==0) {
        *_pv = C;
        return 1;
    }
    if (_nm[0]=='t' && strcmp(_nm, "t")==0) {
        *_pv = t;
        return 1;
    }
	return 0;
}

void _set_var_value(const char* _nm, double* _pv, double zzY[]){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];

    if (_nm[0]=='W' && strcmp(_nm, "WT")==0) {
        WT = *_pv;
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaDose")==0) {
        AaDose = *_pv;
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaInfDose")==0) {
        AaInfDose = *_pv;
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaInfRate")==0) {
        AaInfRate = *_pv;
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "Aa")==0) {
        Aa = *_pv;
        return;
    }
    if (_nm[0]=='A' && strcmp(_nm, "A1")==0) {
        A1 = *_pv;
        return;
    }
    if (_nm[0]=='C' && strcmp(_nm, "CMultStdev")==0) {
        CMultStdev = *_pv;
        return;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvVmax")==0) {
        tvVmax = *_pv;
        return;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvKm")==0) {
        tvKm = *_pv;
        return;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvV")==0) {
        tvV = *_pv;
        return;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvKa")==0) {
        tvKa = *_pv;
        return;
    }
    if (_nm[0]=='n' && strcmp(_nm, "nVmax")==0) {
        nVmax = *_pv;
        return;
    }
    if (_nm[0]=='n' && strcmp(_nm, "nKm")==0) {
        nKm = *_pv;
        return;
    }
    if (_nm[0]=='n' && strcmp(_nm, "nV")==0) {
        nV = *_pv;
        return;
    }
    if (_nm[0]=='V' && strcmp(_nm, "Vmax")==0) {
        Vmax = *_pv;
        return;
    }
    if (_nm[0]=='K' && strcmp(_nm, "Km")==0) {
        Km = *_pv;
        return;
    }
    if (_nm[0]=='V' && strcmp(_nm, "V")==0) {
        V = *_pv;
        return;
    }
    if (_nm[0]=='K' && strcmp(_nm, "Ka")==0) {
        Ka = *_pv;
        return;
    }
    if (_nm[0]=='C' && strcmp(_nm, "C")==0) {
        C = *_pv;
        return;
    }
    if (_nm[0]=='t' && strcmp(_nm, "t")==0) {
        t = *_pv;
        return;
    }
}

void SetVarValue(const char* _nm, double* _pv, double zzY[])
{
    EvalSParms(zzY);

    _set_var_value(_nm, _pv, zzY);

    EvalSParms(zzY);
}

int IsSaemRanef(int * idx)
{
	static const int is_saem[] = {1, 1, 1};
    return is_saem[* idx];
}

void GetSaemCMatrixSize(int* nBaseEta, int* nExtraEta, int* nTheta){
	*nBaseEta = *nTheta = *nExtraEta = 0;

	(*nTheta)++; // tvVmax
	(*nBaseEta)++; // nVmax

	(*nTheta)++; // tvV
	(*nBaseEta)++; // nKm

	(*nTheta)++; // tvKm
	(*nBaseEta)++; // nV

	if ((-1 == -1)||_ISENABLED(-1)){(*nExtraEta)++; (*nTheta)++;} // bare theta CMultStdev

	if ((-1 == -1)||_ISENABLED(-1)){(*nExtraEta)++; (*nTheta)++;} // bare theta tvKa

}
void GetSaemCMatrix(double* daCMat, int* theta2eta){
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	int _nBaseEta, _nExtraEta, _nTheta, _nEta;
	int _iEta = 0, _jTheta = 0;
	GetSaemCMatrixSize(&_nBaseEta, &_nExtraEta, &_nTheta);
	_nEta = _nBaseEta + _nExtraEta;
	// clear the C matrix
	memset(daCMat, 0, sizeof(double)*_nEta*_nTheta);

	// for structural parameter log_Vmax
	daCMat[_iEta + _jTheta * _nEta] = 1; // eta nVmax, theta log_tvVmax
	theta2eta[_jTheta] = _iEta;
	_jTheta++;
	_iEta++;

	// for structural parameter log_V
	daCMat[_iEta + _jTheta * _nEta] = 1; // eta nKm, theta log_tvV
	theta2eta[_jTheta] = _iEta;
	_jTheta++;
	_iEta++;

	// for structural parameter log_Km
	daCMat[_iEta + _jTheta * _nEta] = 1; // eta nV, theta log_tvKm
	theta2eta[_jTheta] = _iEta;
	_jTheta++;
	_iEta++;

	// for additional bare thetas
	{
		daCMat[_iEta + _jTheta * _nEta] = 1; // extra eta for bare theta CMultStdev
		theta2eta[_jTheta] = -1;
		_iEta++;
		_jTheta++;
	}

	{
		daCMat[_iEta + _jTheta * _nEta] = 1; // extra eta for bare theta tvKa
		theta2eta[_jTheta] = -1;
		_iEta++;
		_jTheta++;
	}

	// assert _iEta == _nEta
	// assert _jTheta == _nTheta
}
int GetCovariateIsInCMatrix(const char* name){
	int rval = 0;
	if (0);
	return rval;
}
void GetSaemThetaInfo(int* iaCColToITheta, int* baCColIsLog){
	int _jCTheta = 0;
	int _nTheta = 0;
	int _iaIFixefErrToITheta[5+1];
	int _i;
	for (_i = 0; _i < 5+1; _i++) _iaIFixefErrToITheta[_i] = -1;
	GetFixefErrToThetaMap(_iaIFixefErrToITheta, &_nTheta);

	iaCColToITheta[_jCTheta] = _iaIFixefErrToITheta[1]; // tvVmax
	baCColIsLog[_jCTheta] = 1;
	_jCTheta++;

	iaCColToITheta[_jCTheta] = _iaIFixefErrToITheta[3]; // tvV
	baCColIsLog[_jCTheta] = 1;
	_jCTheta++;

	iaCColToITheta[_jCTheta] = _iaIFixefErrToITheta[2]; // tvKm
	baCColIsLog[_jCTheta] = 1;
	_jCTheta++;

	// for additional bare thetas
	{
		iaCColToITheta[_jCTheta] = _iaIFixefErrToITheta[0]; // CMultStdev
		baCColIsLog[_jCTheta] = 0;
		_jCTheta++;
	}

	{
		iaCColToITheta[_jCTheta] = _iaIFixefErrToITheta[4]; // tvKa
		baCColIsLog[_jCTheta] = 0;
		_jCTheta++;
	}

}
void GetSaemNaivePoolThetas(int* nNPTheta, int* iaTheta){
	int _jTheta = 0;

	*nNPTheta = _jTheta;
}
int CanUseSAEM(){
	return 1;
}
void DoEndObsActions(double zzY[]){
}
#define GRADIENT_DERIV_START (GetGammaDelayDEStart(0) + GetGammaDelayDECount())

#define     _A1_nVmax       zzY1[0]
#define zzR__A1_nVmax       zzR1[0]

#define     _A1_nKm       zzY1[1]
#define zzR__A1_nKm       zzR1[1]

#define     _A1_nV       zzY1[2]
#define zzR__A1_nV       zzR1[2]

#define NINTEG_GRADIENT_DERIV 3

#define     _Aa_nVmax       0.0
#define     _Aa_nKm       0.0
#define     _Aa_nV       0.0

int GetNIntegGradientDeriv(){return GRADIENT_DERIV_START + NINTEG_GRADIENT_DERIV;}



void GradientDeriv(long* _pN, double* _pt, double zzY[], double zzR[], double zzIRate[])
{
    int _i = 0;
    double* _fixef = _args[0];
    double* _ranef1 = _args[1];
    double* _ranef2 = _args[2];
    double* _ranef3 = _args[3];
    double* _ranef4 = _args[4];
    double* _ranef5 = _args[5];

    double * zzY1 = zzY + GRADIENT_DERIV_START;
    double * zzR1 = zzR + GRADIENT_DERIV_START;

	t = *_pt;
	#if defined(USE_SETTIMECF)
		SetTimeOfGetValCF(*_pt);
	#endif
	#if defined(THERE_ARE_INTERPOLATED_COVARIATES)
		EvalSParms(zzY);
	#else
		EvalAssign(zzY);
	#endif
	zzR_Aa = ( -(Ka) * Aa);
	zzR_Aa += _Aa_irate;
	zzR_A1 = ((Ka * Aa) - ((Vmax * C) / (Km + C)));
	zzR_A1 += _A1_irate;
	zzR[2] = 0; // extra null equation
	zzR__A1_nVmax = (((Ka * _Aa_nVmax) + ( -(((((Km + (A1 / V)) * (Vmax / V)) - ((Vmax * (A1 / V)) / V)) / ((Km + (A1 / V)) * (Km + (A1 / V))))) * _A1_nVmax)) +  -((((Km + (A1 / V)) * ((A1 / V) * Vmax)) / ((Km + (A1 / V)) * (Km + (A1 / V))))));
	zzR__A1_nKm = (((Ka * _Aa_nKm) + ( -(((((Km + (A1 / V)) * (Vmax / V)) - ((Vmax * (A1 / V)) / V)) / ((Km + (A1 / V)) * (Km + (A1 / V))))) * _A1_nKm)) +  -(((((Km + (A1 / V)) * (Vmax * ( -((A1 * V)) / (V * V)))) - ((Vmax * (A1 / V)) * ( -((A1 * V)) / (V * V)))) / ((Km + (A1 / V)) * (Km + (A1 / V))))));
	zzR__A1_nV = (((Ka * _Aa_nV) + ( -(((((Km + (A1 / V)) * (Vmax / V)) - ((Vmax * (A1 / V)) / V)) / ((Km + (A1 / V)) * (Km + (A1 / V))))) * _A1_nV)) +  -(( -(((Vmax * (A1 / V)) * Km)) / ((Km + (A1 / V)) * (Km + (A1 / V))))));
}
	//[0,0]  -(tvKa)
	//[1,0] tvKa
	//[3,0] 0.0
	//[4,0] 0.0
	//[5,0] 0.0

	//[0,1] 0.0
	//[1,1]  -((((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) - (((tvVmax * exp(nVmax)) * (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))) / (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))))))  is not constant
	//[3,1] ((_A1_nVmax *  -(((((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) * ((((tvVmax * exp(nVmax)) / (tvV * exp(nKm))) / (tvV * exp(nKm))) - (((tvV * exp(nKm)) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) / ((tvV * exp(nKm)) * (tvV * exp(nKm)))))) - (((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) - (((tvVmax * exp(nVmax)) * (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))) * ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm))) + (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))))) / ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) * (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))))))) +  -(((((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) * ((((A1 / (tvV * exp(nKm))) * (tvVmax * exp(nVmax))) / (tvV * exp(nKm))) + (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))))) - ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((A1 / (tvV * exp(nKm))) * (tvVmax * exp(nVmax)))) * ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm))) + (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))))) / ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) * (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))))))))  is not constant
	//[4,1] ((_A1_nKm *  -(((((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) * ((((tvVmax * exp(nVmax)) / (tvV * exp(nKm))) / (tvV * exp(nKm))) - (((tvV * exp(nKm)) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) / ((tvV * exp(nKm)) * (tvV * exp(nKm)))))) - (((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) - (((tvVmax * exp(nVmax)) * (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))) * ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm))) + (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))))) / ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) * (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))))))) +  -(((((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) * (((((tvVmax * exp(nVmax)) * ( -((A1 * (tvV * exp(nKm)))) / ((tvV * exp(nKm)) * (tvV * exp(nKm))))) / (tvV * exp(nKm))) + (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) * ((((tvV * exp(nKm)) * (tvV * exp(nKm))) *  -((tvV * exp(nKm)))) / (((tvV * exp(nKm)) * (tvV * exp(nKm))) * ((tvV * exp(nKm)) * (tvV * exp(nKm)))))))) - ((( -((A1 * (tvV * exp(nKm)))) / ((tvV * exp(nKm)) * (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) + (((tvVmax * exp(nVmax)) * (A1 / (tvV * exp(nKm)))) * ((((tvV * exp(nKm)) * (tvV * exp(nKm))) *  -((tvV * exp(nKm)))) / (((tvV * exp(nKm)) * (tvV * exp(nKm))) * ((tvV * exp(nKm)) * (tvV * exp(nKm))))))))) - (((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) * ( -((A1 * (tvV * exp(nKm)))) / ((tvV * exp(nKm)) * (tvV * exp(nKm)))))) - (((tvVmax * exp(nVmax)) * (A1 / (tvV * exp(nKm)))) * ( -((A1 * (tvV * exp(nKm)))) / ((tvV * exp(nKm)) * (tvV * exp(nKm)))))) * ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm))) + (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))))) / ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) * (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))))))))  is not constant
	//[5,1] ((_A1_nV *  -(((((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) * ((((tvVmax * exp(nVmax)) / (tvV * exp(nKm))) / (tvV * exp(nKm))) - (((tvV * exp(nKm)) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) / ((tvV * exp(nKm)) * (tvV * exp(nKm)))))) - (((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) - (((tvVmax * exp(nVmax)) * (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))) * ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm))) + (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))))) / ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) * (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))))))) +  -(((((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) *  -(((tvKm * exp(nV)) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))))) - ( -((((tvVmax * exp(nVmax)) * (A1 / (tvV * exp(nKm)))) * (tvKm * exp(nV)))) * ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm))) + (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))))) / ((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm))))) * (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))))))))  is not constant

	//[0,3] 0.0
	//[1,3] 0.0
	//[3,3]  -((((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) - (((tvVmax * exp(nVmax)) * (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))) / (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))))))  is not constant
	//[4,3] 0.0
	//[5,3] 0.0

	//[0,4] 0.0
	//[1,4] 0.0
	//[3,4] 0.0
	//[4,4]  -((((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) - (((tvVmax * exp(nVmax)) * (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))) / (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))))))  is not constant
	//[5,4] 0.0

	//[0,5] 0.0
	//[1,5] 0.0
	//[3,5] 0.0
	//[4,5] 0.0
	//[5,5]  -((((((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvVmax * exp(nVmax)) / (tvV * exp(nKm)))) - (((tvVmax * exp(nVmax)) * (A1 / (tvV * exp(nKm)))) / (tvV * exp(nKm)))) / (((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))) * ((tvKm * exp(nV)) + (A1 / (tvV * exp(nKm)))))))  is not constant

	//[0,2] 0.0
	//[1,2] 0.0
void GradientJacobian(int* pN, double* pt, double zzY[], double zzIRate[], int* pML, int* pMU, double* _gjacobian, int* pNRJ){
	double * zzY1 = zzY + GRADIENT_DERIV_START;
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	int _nrj = *pNRJ;
	double __temp[3];
	__temp[0] = exp(nV);
	__temp[1] = exp(nKm);
	__temp[2] = exp(nVmax);
	memset(_gjacobian, 0, sizeof(double)*_nrj*_nrj);
	// rate terms for Aa
		// w.r.t. Aa
	_gjacobian[0 + _nrj*0] =  -(tvKa);
		// forcing rate for Aa
	_gjacobian[0 + _nrj*2] = 0.0 + zzIRate[0];
	// rate terms for A1
		// w.r.t. Aa
	_gjacobian[1 + _nrj*0] = tvKa;
		// w.r.t. A1
	_gjacobian[1 + _nrj*1] =  -((((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) - (((tvVmax * __temp[2]) * (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))) / (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))))));
		// forcing rate for A1
	_gjacobian[1 + _nrj*2] = 0.0 + zzIRate[1];
	// rate terms for zzR[2]
	// rate terms for _A1_nVmax
		// w.r.t. A1
	_gjacobian[3 + _nrj*1] = ((_A1_nVmax *  -(((((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) * ((((tvVmax * __temp[2]) / (tvV * __temp[1])) / (tvV * __temp[1])) - (((tvV * __temp[1]) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) / ((tvV * __temp[1]) * (tvV * __temp[1]))))) - (((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) - (((tvVmax * __temp[2]) * (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))) * ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1])) + (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))))) / ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) * (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))))))) +  -(((((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) * ((((A1 / (tvV * __temp[1])) * (tvVmax * __temp[2])) / (tvV * __temp[1])) + (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))))) - ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((A1 / (tvV * __temp[1])) * (tvVmax * __temp[2]))) * ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1])) + (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))))) / ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) * (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))))))));
		// w.r.t. _A1_nVmax
	_gjacobian[3 + _nrj*3] =  -((((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) - (((tvVmax * __temp[2]) * (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))) / (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))))));
	// rate terms for _A1_nKm
		// w.r.t. A1
	_gjacobian[4 + _nrj*1] = ((_A1_nKm *  -(((((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) * ((((tvVmax * __temp[2]) / (tvV * __temp[1])) / (tvV * __temp[1])) - (((tvV * __temp[1]) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) / ((tvV * __temp[1]) * (tvV * __temp[1]))))) - (((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) - (((tvVmax * __temp[2]) * (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))) * ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1])) + (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))))) / ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) * (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))))))) +  -(((((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) * (((((tvVmax * __temp[2]) * ( -((A1 * (tvV * __temp[1]))) / ((tvV * __temp[1]) * (tvV * __temp[1])))) / (tvV * __temp[1])) + (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvVmax * __temp[2]) * ((((tvV * __temp[1]) * (tvV * __temp[1])) *  -((tvV * __temp[1]))) / (((tvV * __temp[1]) * (tvV * __temp[1])) * ((tvV * __temp[1]) * (tvV * __temp[1]))))))) - ((( -((A1 * (tvV * __temp[1]))) / ((tvV * __temp[1]) * (tvV * __temp[1]))) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) + (((tvVmax * __temp[2]) * (A1 / (tvV * __temp[1]))) * ((((tvV * __temp[1]) * (tvV * __temp[1])) *  -((tvV * __temp[1]))) / (((tvV * __temp[1]) * (tvV * __temp[1])) * ((tvV * __temp[1]) * (tvV * __temp[1])))))))) - (((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvVmax * __temp[2]) * ( -((A1 * (tvV * __temp[1]))) / ((tvV * __temp[1]) * (tvV * __temp[1]))))) - (((tvVmax * __temp[2]) * (A1 / (tvV * __temp[1]))) * ( -((A1 * (tvV * __temp[1]))) / ((tvV * __temp[1]) * (tvV * __temp[1]))))) * ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1])) + (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))))) / ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) * (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))))))));
		// w.r.t. _A1_nKm
	_gjacobian[4 + _nrj*4] =  -((((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) - (((tvVmax * __temp[2]) * (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))) / (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))))));
	// rate terms for _A1_nV
		// w.r.t. A1
	_gjacobian[5 + _nrj*1] = ((_A1_nV *  -(((((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) * ((((tvVmax * __temp[2]) / (tvV * __temp[1])) / (tvV * __temp[1])) - (((tvV * __temp[1]) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) / ((tvV * __temp[1]) * (tvV * __temp[1]))))) - (((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) - (((tvVmax * __temp[2]) * (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))) * ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1])) + (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))))) / ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) * (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))))))) +  -(((((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) *  -(((tvKm * __temp[0]) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))))) - ( -((((tvVmax * __temp[2]) * (A1 / (tvV * __temp[1]))) * (tvKm * __temp[0]))) * ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1])) + (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))))) / ((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1])))) * (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))))))));
		// w.r.t. _A1_nV
	_gjacobian[5 + _nrj*5] =  -((((((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvVmax * __temp[2]) / (tvV * __temp[1]))) - (((tvVmax * __temp[2]) * (A1 / (tvV * __temp[1]))) / (tvV * __temp[1]))) / (((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))) * ((tvKm * __temp[0]) + (A1 / (tvV * __temp[1]))))));
}

void NewGetLL(int bGenGradients, double _dLL[], double _dLLdEta[], double zzY[], double zzIRate[], int* pnObs, double _t, const char* _nm, double _v1, int bBql)
{
    double * zzY1 = zzY + GRADIENT_DERIV_START;

    double _yval[12];
    double _probge[12];
    double _prob[12];
    double _probge1[12];
    double _prob1[12];
    double _dyde[12];
    double _dldy[12];

    const double t = _t;

    double* _fixef = _args[0];
    double* _ranef1 = _args[1];
    double* _ranef2 = _args[2];
    double* _ranef3 = _args[3];
    double* _ranef4 = _args[4];
    double* _ranef5 = _args[5];

    double _beta, _power, _zprob;

    double _ypred = 0, _var = 1, _sig = 1, _sqrtvf = 1, _std = 1;
    double _dLL1 = 0, _dLL2;
    double _dYdEta = 0, _dSdEta = 0, _dLLdY = 0, _temp = 0, _haz0 = 0, _haz1 = 0, _del = 1e-6;
    int _i = 0;

	double _dLLdEta1[3];
	memset(_dLLdEta1, 0, sizeof(_dLLdEta1));
	memset(_yval, 0, sizeof(_yval));
	memset(_probge, 0, sizeof(_probge));
	memset(_prob, 0, sizeof(_prob));
	memset(_probge1, 0, sizeof(_probge1));
	memset(_prob1, 0, sizeof(_prob1));
	memset(_dyde, 0, sizeof(_dyde));
	memset(_dldy, 0, sizeof(_dldy));

    if (_nm[0]=='C' && strcmp(_nm, "CObs")==0){
		_std = _fixef[5];
		_ypred = C;
		if (!isfinite(_ypred)){
			SetErrorMsg("Error (Model.cpp): non-finite _ypred");
		}
		_sqrtvf = sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0))));
		/* handle 0 observations and proportional models */
		if (_sqrtvf == 0 || (_sqrtvf < 1e-12 && _v1 == 0)){_ypred = 0; _sqrtvf = 1;}
		_sig = _sqrtvf;
		if (!isfinite(_sig)){
			SetErrorMsg("Error (Model.cpp): non-finite _sig");
		}
		if (_sig == 0){
			SetErrorMsg("Error (Model.cpp): zero _sig");
		}
		// get base LL
		if (ALLOWLOGTRANSFORM && bEnableLogTransform){
			if (bBql){
				_dLL1 = lphi(logfl(_v1) - logfl(_ypred), 1 * _std); // add current LL to log(P)
			} else {
				_dLL1 = lnorm(logfl(_v1) - logfl(_ypred), 1 * _std); // add current LL to log(P)
			}
		} else {
			if (bBql){
				_dLL1 = lphi((_v1 - _ypred), _sig * _std); // add current LL to log(P)
			} else {
				_dLL1 = lnorm((_v1 - _ypred), _sig * _std); // add current LL to log(P)
			}
		}

		if (bGenGradients){

			if (ALLOWLOGTRANSFORM && bEnableLogTransform){
				_dYdEta = ((V * _A1_nVmax) / (V * V)) / _ypred; // Get dY / dnVmax
				_dSdEta = 0; // Get dS / dnVmax
			} else {
				_dYdEta = ((V * _A1_nVmax) / (V * V)); // Get dY / dnVmax
				_dSdEta = 0.0; // Get dS / dnVmax
			}
			if (ALLOWLOGTRANSFORM && bEnableLogTransform){
				if (bBql){
					_dLL2 = lphi(logfl(_v1) - logfl(_ypred + _dYdEta*_del), (1+_dSdEta*_del)*_std);
					_dLLdEta1[0] = (_dLL2 - _dLL1)/_del;
				} else {
					_dLLdEta1[0]  = dLNormdEta(logfl(_v1) - logfl(_ypred), 1*_std, -_dYdEta, _dSdEta*_std);
				}
			} else {
				if (bBql){
					_dLL2 = lphi((_v1 - (_ypred + _dYdEta*_del)), (_sig+_dSdEta*_del)*_std);
					_dLLdEta1[0] = (_dLL2 - _dLL1)/_del;
				} else {
					_dLLdEta1[0]  = dLNormdEta((_v1 - _ypred), _sig*_std, -_dYdEta, _dSdEta*_std);
				}
			}

			if (ALLOWLOGTRANSFORM && bEnableLogTransform){
				_dYdEta = (((V * _A1_nKm) - (A1 * V)) / (V * V)) / _ypred; // Get dY / dnKm
				_dSdEta = 0; // Get dS / dnKm
			} else {
				_dYdEta = (((V * _A1_nKm) - (A1 * V)) / (V * V)); // Get dY / dnKm
				_dSdEta = ((0.5 * pow((1.0 + (pow((A1 / V), 2.0) * pow((CMultStdev / sigma()), 2.0))), -0.5)) * (pow((CMultStdev / sigma()), 2.0) * ((2.0 * (A1 / V)) * (((V * _A1_nKm) - (A1 * V)) / (V * V))))); // Get dS / dnKm
			}
			if (ALLOWLOGTRANSFORM && bEnableLogTransform){
				if (bBql){
					_dLL2 = lphi(logfl(_v1) - logfl(_ypred + _dYdEta*_del), (1+_dSdEta*_del)*_std);
					_dLLdEta1[1] = (_dLL2 - _dLL1)/_del;
				} else {
					_dLLdEta1[1]  = dLNormdEta(logfl(_v1) - logfl(_ypred), 1*_std, -_dYdEta, _dSdEta*_std);
				}
			} else {
				if (bBql){
					_dLL2 = lphi((_v1 - (_ypred + _dYdEta*_del)), (_sig+_dSdEta*_del)*_std);
					_dLLdEta1[1] = (_dLL2 - _dLL1)/_del;
				} else {
					_dLLdEta1[1]  = dLNormdEta((_v1 - _ypred), _sig*_std, -_dYdEta, _dSdEta*_std);
				}
			}

			if (ALLOWLOGTRANSFORM && bEnableLogTransform){
				_dYdEta = ((V * _A1_nV) / (V * V)) / _ypred; // Get dY / dnV
				_dSdEta = 0; // Get dS / dnV
			} else {
				_dYdEta = ((V * _A1_nV) / (V * V)); // Get dY / dnV
				_dSdEta = 0.0; // Get dS / dnV
			}
			if (ALLOWLOGTRANSFORM && bEnableLogTransform){
				if (bBql){
					_dLL2 = lphi(logfl(_v1) - logfl(_ypred + _dYdEta*_del), (1+_dSdEta*_del)*_std);
					_dLLdEta1[2] = (_dLL2 - _dLL1)/_del;
				} else {
					_dLLdEta1[2]  = dLNormdEta(logfl(_v1) - logfl(_ypred), 1*_std, -_dYdEta, _dSdEta*_std);
				}
			} else {
				if (bBql){
					_dLL2 = lphi((_v1 - (_ypred + _dYdEta*_del)), (_sig+_dSdEta*_del)*_std);
					_dLLdEta1[2] = (_dLL2 - _dLL1)/_del;
				} else {
					_dLLdEta1[2]  = dLNormdEta((_v1 - _ypred), _sig*_std, -_dYdEta, _dSdEta*_std);
				}
			}

			_VecAddTo(3, _dLLdEta1, _dLLdEta);
		}
		_dLL[0] += _dLL1;

		(*pnObs)++; // count all observations

        return;
    }
}

void NewGetPred(int bGenGradients, int* pnObs, double zzY[], double zzIRate[], double _t, const char* _nm, double _v1, bool doActions /*= true*/)
{
	double * zzY1 = zzY + GRADIENT_DERIV_START;
	double _dVFCorrectionFactor = 1;
	const double t = _t;
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	double _ypred = 0; // predicted value
	double _var = 1;
	double _sig = 1; // current standard deviation of Epsilon
	double _std = 1; // current standard deviation of Epsilon
	double _sqrtvf = 1; // dObs / dEpsilon (additive = 1, proportional = _ypred)
	double _vf = 1; // _sqrtvf*_sqrtvf (additive = 1, proportional = _ypred^2)
	double _dYdEta[3];
	double _dSdEta[3];
	_dYdEta[0] = _dYdEta[1] = _dYdEta[2] = 0;
	_dSdEta[0] = _dSdEta[1] = _dSdEta[2] = 0;

    if (_nm[0]=='C' && strcmp(_nm, "CObs")==0){
		_std = _fixef[5];
		_ypred = C;
		if (!isfinite(_ypred)){
			SetErrorMsg("Error (Model.cpp): non-finite _ypred");
		}
		_sqrtvf = sqrt((sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))) * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0))))));
		/* handle 0 observations and proportional models */
		if (_sqrtvf == 0 || (_sqrtvf < 1e-12 && _v1 == 0)){_ypred = 0; _sqrtvf = 1;}
		_sig = _sqrtvf;
		_dVFCorrectionFactor = GetPredVFCorrectionFactor(0);
		if (_dVFCorrectionFactor != 0) _sig /= _dVFCorrectionFactor;
		if (!isfinite(_sig)){
			SetErrorMsg("Error (Model.cpp): non-finite _sig");
		}
		if (_sig == 0){
			SetErrorMsg("Error (Model.cpp): zero _sig");
		}
		_var = _sig * _sig;
		if (ALLOWLOGTRANSFORM && bEnableLogTransform){
			RecordObservation(t, curTimeLastDose, iCurWhichDose, iCurWhichReset, logfl(_ypred), logfl(_v1), 1, 0);
		} else {
			RecordObservation(t, curTimeLastDose, iCurWhichDose, iCurWhichReset, _ypred, _v1, _var, 0);
		}

		if (bGenGradients){
			_dYdEta[0] = ((V * _A1_nVmax) / (V * V)); // Get dY / dnVmax
			_dSdEta[0] = 0.0; // Get dS / dnVmax
			_dYdEta[1] = (((V * _A1_nKm) - (A1 * V)) / (V * V)); // Get dY / dnKm
			_dSdEta[1] = ((0.5 * pow((1.0 + (pow((A1 / V), 2.0) * pow((CMultStdev / sigma()), 2.0))), -0.5)) * (pow((CMultStdev / sigma()), 2.0) * ((2.0 * (A1 / V)) * (((V * _A1_nKm) - (A1 * V)) / (V * V))))); // Get dS / dnKm
			_dYdEta[2] = ((V * _A1_nV) / (V * V)); // Get dY / dnV
			_dSdEta[2] = 0.0; // Get dS / dnV
			if (ALLOWLOGTRANSFORM && bEnableLogTransform){
				RecordObservationGradient(_dYdEta[0] / (_ypred + 1e-30), 0, 0);
				RecordObservationGradient(_dYdEta[1] / (_ypred + 1e-30), 0, 1);
				RecordObservationGradient(_dYdEta[2] / (_ypred + 1e-30), 0, 2);
			} else {
				RecordObservationGradient(_dYdEta[0], _dSdEta[0], 0);
				RecordObservationGradient(_dYdEta[1], _dSdEta[1], 1);
				RecordObservationGradient(_dYdEta[2], _dSdEta[2], 2);
			}
		}
		(*pnObs)++; // count all observations

        return;
    }
}

void GetLLGradientEta(
	  int* _pnObs
	, const unsigned char* *ppb
	, const SubjectId & subjId
	, double _dLL[]
	, double _dGrad[]
	)
{
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	int _nEta = 3;
	int _nObs = 0;
	double _dLLLocal[1];
	double _dGradLocal[3];
	memset(_dGrad, 0, sizeof(double)*3);
	_dLLLocal[0] = 0;
	memset(_dGradLocal, 0, sizeof(_dGradLocal));
	const unsigned char* pb0 = *ppb;
	*ppb = pb0;
	_nObs = 0;
	LLOneSubject1(TRUE, _nEta, &_nObs, ppb, subjId, _dLL, _dGrad, _pbOK);
}
void NewPredOneSubject(
	  int* _pnObs
	, const unsigned char* *ppb
	, const SubjectId & subjId
	, int* _pbOK
	, int* _nPrediction
	, double* __daXPrediction
	, double* __taLastDose
	, int  * __iaWhichDose
	, int  * __iaWhichReset
	, double* __daYPrediction
	, double* __daYObservation
	, double* __daVarianceFunction
	, double* _dYdK // comes out as _daYdK(iEta, iObs)
	, double* _dSdK // comes out as _daSdK(iEta, iObs)
	, int* _iaWhichObs
	)
{
	int _nEta = 3;
	int _nObs = 0;
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	NewPredOneSubject1(TRUE, _nEta, &_nObs, ppb, subjId, _pbOK, _nPrediction, __daXPrediction, __taLastDose, __iaWhichDose, __iaWhichReset, __daYPrediction, __daYObservation, __daVarianceFunction, _dYdK, _dSdK, _iaWhichObs);
}
const char * GetGradientErrorMessage()
{
	return "";
}

// jacobian info: 2=can matexp, 1=can stiffjac, 0=no jac, -1=no odes
void GetJacobianInfo(int* bUsingSyntheticGradients, int* _piJacobianInfo){
	if (*bUsingSyntheticGradients){
		* _piJacobianInfo = 1;
	} else {
		* _piJacobianInfo = 1;
	}
}
void GetMaxODELevel(int* bUsingSyntheticGradients, int* _piMaxODELevel){
	if (*bUsingSyntheticGradients){
		* _piMaxODELevel = ODELEVEL_STIFFJAC;
	} else {
		* _piMaxODELevel = ODELEVEL_STIFFJAC;
	}
}
int GetNumObservations(){
	return 1;
}

int GetObservationNumber(const char* _nm){
	int _i = -1;
	if (0);
	else if (strcmp(_nm, "CObs")==0) _i = 0;
	return _i;
}

const char * GetObservationName(int* _pi)
{
	if (0);
	else if (*_pi == 0) return "CObs";
    return "";
}

ObsType GetObservationType(int * _pi){
	ObsType _type = ObsType::Unknown;
	if (0);
	else if (*_pi == 0) _type = ObsType::Gaussian;
	return _type;
}

void GetObservationBqlLimit(int* _pi, double* _bqllimit){
	*_bqllimit = 0;
	if (0);
}

void GetLL(int bGenGradients, double _dLL[], double _dGrad[], double zzY[], double zzIRate[], int* pnObs, double _t, const char* _nm, double _v1, int bBql){
	double _dLL1 = 0;
	double _probge[12];
	const double t = _t;
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	double _ypred = 0, _sig = 1, _std = 1, _var = 1, _sqrtvf = 1;
	int _i = 0;

    if (_nm[0]=='C' && strcmp(_nm, "CObs")==0){
		_std = _fixef[5];
		_ypred = C;
		if (!isfinite(_ypred)){
			SetErrorMsg("Error (Model.cpp): non-finite _ypred");
		}
		if (ALLOWLOGTRANSFORM && bEnableLogTransform){
			_sqrtvf = 1;
		} else {
			_sqrtvf = sqrt((sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))) * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0))))));
		}
		/* handle 0 observations and proportional models */
		if (_sqrtvf == 0 || (_sqrtvf < 1e-12 && _v1 == 0)){_ypred = 0; _sqrtvf = 1;}
		_sig = _sqrtvf;
		_sig *= _fixef[5];
		if (!isfinite(_sig)){
			SetErrorMsg("Error (Model.cpp): non-finite _sig");
		}
		if (_sig == 0){
			SetErrorMsg("Error (Model.cpp): zero _sig");
		}
		_var = _sig * _sig;
		if (ALLOWLOGTRANSFORM && bEnableLogTransform){
			if (bBql){
				_dLL1 += lphi(logfl(_v1) - logfl(_ypred), _sig); // add current LL to log(P)
			}
			else {
				_dLL1 += lnorm(logfl(_v1) - logfl(_ypred), _sig); // add current LL to log(P)
			}
		} else {
			if (bBql){
				_dLL1 += lphi((_v1 - _ypred), _sig); // add current LL to log(P)
			} else {
				_dLL1 += lnorm((_v1 - _ypred), _sig); // add current LL to log(P)
			}
		}
		_dLL[0] += _dLL1;
		(*pnObs)++;
		if (_pcb){
			if (ALLOWLOGTRANSFORM && bEnableLogTransform){
				_pcb->OnObserve(curTime, curTimeLastDose, _nm, -1, logfl(_ypred), 0, logfl(_v1), _sig, _sqrtvf, bBql);
			} else {
				_pcb->OnObserve(curTime, curTimeLastDose, _nm, -1, _ypred, 0, _v1, _sig, _sqrtvf, bBql);
			}
		}

        return;
    }
}

void GetPred(int* pnObs, double zzY[], double zzIRate[], double _t, const char* _nm, double _v1){
	double _dVFCorrectionFactor = 1;
	double _probge[12];
	const double t = _t;
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	double _ypred = 0, _sig = 1, _std = 1, _var = 1, _sqrtvf = 1;
	int _i = 0;

    if (_nm[0]=='C' && strcmp(_nm, "CObs")==0){
		_std = _fixef[5];
		_ypred = C;
		if (!isfinite(_ypred)){
			SetErrorMsg("Error (Model.cpp): non-finite _ypred");
		}
		_sqrtvf = sqrt((sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))) * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0))))));
		/* handle 0 observations and proportional models */
		if (_sqrtvf == 0 || (_sqrtvf < 1e-12 && _v1 == 0)){_ypred = 0; _sqrtvf = 1;}
		_sig = _sqrtvf;
		_dVFCorrectionFactor = GetPredVFCorrectionFactor(0);
		if (_dVFCorrectionFactor != 0) _sig /= _dVFCorrectionFactor;
		if (!isfinite(_sig)){
			SetErrorMsg("Error (Model.cpp): non-finite _sig");
		}
		if (_sig == 0){
			SetErrorMsg("Error (Model.cpp): zero _sig");
		}
		_var = _sig * _sig;
		if (ALLOWLOGTRANSFORM && bEnableLogTransform){
			RecordObservation(t, curTimeLastDose, iCurWhichDose, iCurWhichReset, logfl(_ypred), logfl(_v1), 1, 0);
		} else {
			RecordObservation(t, curTimeLastDose, iCurWhichDose, iCurWhichReset, _ypred, _v1, _var, 0);
		}
		(*pnObs)++;

        return;
    }
}

int GetPredObsError(int _iObs){
	int _iError = -1;
	if (0){
	} else if (_iObs == 0){
		_iError = 0;
	}
	return _iError;
}
double GetPredVFCorrectionFactor(int _iError){
	double* _fixef = _args[0];
	int _nFixef = 5;
	int _iLastFreeError = 0;
	int _bErrorFrozen = 0;
	double _vfCorrectionFactor = 1;
	if (0){} // code to detect if _iError is frozen
	if (!_bErrorFrozen && _iError != _iLastFreeError){
		double _dThisSigma = _fixef[_nFixef + _iError];
		if (_dThisSigma != 0){
			_vfCorrectionFactor = _dMainSigma / _dThisSigma;
		}
	}
	return _vfCorrectionFactor;
}
void GetPredName(int* _iWhichObs, char* _nm){
	_nm[0] = '\0';
	if (0);
	else if (*_iWhichObs == 0) strcpy(_nm, "C");
}
void GetObsName(int* _iWhichObs, char* _nm){
	_nm[0] = '\0';
	if (0);
	else if (*_iWhichObs == 0) strcpy(_nm, "CObs");
}
void GetIVarName(int* _iWhichObs, char* _nm){
	_nm[0] = '\0';
	if (0);
	else if (*_iWhichObs == 0) strcpy(_nm, "t");
}
void GetObsVF(int* _iWhichObs, double* _vf, double zzY[]){
	#define t (_t + 0)
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	*_vf = 1;
	if (0);
	else if (*_iWhichObs == 0){
		double _std = _fixef[5];
		*_vf = (sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))) * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))));
	}
	#undef t
}
void GetNumGaussianObs(int* _nGaussianObs){
	(*_nGaussianObs) = 0;
	(*_nGaussianObs)++; /* CObs */
}
double _unifControlled(BOOL bSample){
	return ((bSample && _pcb) ? _pcb->GetUniform() : 0.5);
}
double _normControlled(BOOL bSample){
	return ((bSample && _pcb) ? _pcb->GetNormal() : 0.0);
}
void Simulate(double _t, const char* _nm, double _v1, int bBql, int bSampleEpsilon, int bForVPC, double zzY[], double zzIRate[]){
	double _probge[12];
	const double t = _t;
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	double _ypred = 0, _sig = 1, _std = 1, _var = 1, _sqrtvf = 1;
	int _i = 0;
	EvalSParms(zzY);

    if (_nm[0]=='C' && strcmp(_nm, "CObs")==0){
		_std = _fixef[5];
		_ypred = C;
		if (!isfinite(_ypred)){
			SetErrorMsg("Error (Model.cpp): non-finite _ypred");
		}
		if (ALLOWLOGTRANSFORM && bEnableLogTransform){
			_sqrtvf = 1;
		} else {
			_sqrtvf = sqrt((sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))) * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0))))));
		}
		/* handle 0 observations and proportional models */
		if (_sqrtvf == 0 || (_sqrtvf < 1e-12 && _v1 == 0)){_ypred = 0; _sqrtvf = 1;}
		_sig = _sqrtvf;
		if (!isfinite(_sig)){
			SetErrorMsg("Error (Model.cpp): non-finite _sig");
		}
		if (_sig == 0){
			SetErrorMsg("Error (Model.cpp): zero _sig");
		}
		_var = _sig * _sig;
		if (_pcb != 0){
			double _eps = (bSampleEpsilon ? _pcb->GetNormal() : 0); // get a N(0,1)
			double _ypred1 = 0;
			_eps *= _fixef[5];
			_ypred1 = (C + (_eps * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0))))));
			if (ALLOWLOGTRANSFORM && bEnableLogTransform){
				_pcb->OnObserve(curTime, curTimeLastDose, _nm, -1, logfl(_ypred), logfl(_ypred1), logfl(_v1), _sig, _sqrtvf, bBql);
			} else {
				_pcb->OnObserve(curTime, curTimeLastDose, _nm, -1, _ypred, _ypred1, _v1, _sig, _sqrtvf, bBql);
			}
		}

        return;
    }
}

void GetNumCategories(const char* _nm, int* pNCat){
	*pNCat = 0;
}

void SimulateForSimTbl(double _t, const char* _nm, double* _pObs, double* _pPred, double zzY[]){
	const double t = _t;
	double* _fixef = _args[0];
	double* _ranef1 = _args[1];
	double* _ranef2 = _args[2];
	double* _ranef3 = _args[3];
	double* _ranef4 = _args[4];
	double* _ranef5 = _args[5];
	double _eps = 0;
	EvalSParms(zzY);
	if (_pcb == 0) return;

    if (_nm[0]=='C' && strcmp(_nm, "CObs")==0){
        double _std = 0;
        _std = _fixef[5];
        // MD 151026: Tentative change: Turn off log transformation for sim tables
        if (FALSE && ALLOWLOGTRANSFORM && bEnableLogTransform){
            *_pPred = logfl((C + (_eps * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))))));
            _eps = _pcb->GetNormal() * _fixef[5]; // get a N(0,stdev)
            *_pObs = logfl((C + (_eps * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))))));
        } else {
            *_pPred = (C + (_eps * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0))))));
            _eps = _pcb->GetNormal() * _fixef[5]; // get a N(0,stdev)
            *_pObs = (C + (_eps * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0))))));
        }

        return;
    }
}


} // extern "C"
