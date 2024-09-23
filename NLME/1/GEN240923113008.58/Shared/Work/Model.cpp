// Model.cpp
// Generated with TDL5 23.10.1.0000

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
#include "DualNumber.h"
#include "dverk.h"
#include "mutil/closed_form.h"

#undef  _NO_OLDNAMES
#undef   NO_OLDNAMES

// this allows sigma of current epsilon to be referred to in the model
#define sigma() _std

using std::min;
using std::max;
using std::isfinite;

using DualNumber = DN<NRANEF5>;

static constexpr int GetGammaDelayDEStart(int i_);

static int GetNDeriv(){return GetGammaDelayDEStart(0) + GetGammaDelayDECount();}

std::vector<double> gIRate(GetNDeriv(), 0);

double * zzIRate = &gIRate[0];

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

inline bool isfinite(const DualNumber & x)
{
	return isfinite(x.a);
}

static DualNumber ilogit(const DualNumber & x)
{
    DualNumber e(x);

    if (e < -20.0)
    {
        e = DualNumber(-20.0);
    }

    if (e > 20.0)
    {
        e = DualNumber(20.0);
    }

    e = exp(x);
    return e / (1 + e);
}

inline void copy_grad(const DualNumber & ll, double grad[])
{
    if (!grad)
    {
        return;
    }

    for (int i = 0; i < NRANEF5; ++i)
    {
        grad[i] += ll.b[i];
    }
}

inline void copy_grad(double ll, double grad[])
{
}

struct temp_arg
{
	double dbl;
	DualNumber dual;

	double & operator = (double d)
	{
		return dbl = d;
	}

	DualNumber & operator = (const DualNumber & d)
    {
        return dual = d;
    }

	operator DualNumber & ()
	{
		return dual;
	}

    operator double & ()
    {
        return dbl;
    }

	void reset()
	{
		dbl = -999999999.0;
		dual = DualNumber(-999999999.0);
	}
};

struct cf_machines
{
	cf_m<DualNumber>	cfm_du;
    cf_m<double>		cfm_do;

	operator cf_m<double> * ()
	{
		return &cfm_do;
	}

    operator cf_m<DualNumber> * ()
    {
        return &cfm_du;
    }

    void reset()
    {
        memset(this, 0, sizeof(cf_machines));
    }
};

struct base_data
{
    virtual void init_ranef(const double * ran) = 0;

    virtual void init_subject() = 0;
    virtual void reset_subject() = 0;
    virtual void end_subject(int bSimulating) = 0;

    virtual void advance(int iODELevel, double *pTimeOuter, double _dt, int bTimeAdvances) = 0;

    virtual void init_delays(double _t) = 0;
    virtual void init_cf(double _t) = 0;
    virtual void reinit_cf(double _t) = 0;

    virtual void eval_sparams() = 0;

    virtual void get_secondary(int* _pn, double* _sec) = 0;

    virtual int  get_var_value(const char * _nm, double * _pv) = 0;
    virtual void set_var_value(const char * _nm, double * _pv) = 0;

    virtual void schedule_bolus1(int _mark, double _tRel, int _iCpt, double _amt, int _addl) = 0;
    virtual void perform_bolus1(double * _pt, int _iCpt, double _dv, double _rate) = 0;
    virtual void schedule_iv1(int _mark, double _tRel, int _iCpt, double _amt, double _rate, int _addl) = 0;
    virtual void perform_inf_start1(double * _pt, int _iCpt, double _dv) = 0;
    virtual void perform_inf_end1(double * _pt, int _iCpt, double _dv) = 0;
    virtual void schedule_bolus2(int _mark, double _tRel, int _iCpt, double _amt, int _addl) = 0;
    virtual void schedule_iv2(int _mark, double _tRel, int _iCpt, double _amt, double _rate, int _addl) = 0;
    virtual void perform_bolus2(double * _pt, int _iCpt, double _dv, double _rate) = 0;
    virtual void perform_inf_start2(double * _pt, int _iCpt, double _dv) = 0;
    virtual void perform_inf_end2(double * _pt, int _iCpt, double _dv) = 0;
    virtual void add_to_cpt(int _iCpt, double _dv) = 0;

    virtual void get_ll(double _dLL[], double _dGrad[], int * pnObs, double _t, const char * _nm, double _v1, int bBql) = 0;
    virtual void new_get_pred(int * pnObs, double _t, const char * _nm, double _v1, bool doActions) = 0;

    virtual void get_obs_vf(int* _iWhichObs, double* _vf) = 0;
    virtual void end_obs_actions() = 0;

    virtual void simulate(double _t, const char* _nm, double _v1, int bBql, int bSampleEpsilon, int bForVPC) = 0;
    virtual void simulate_for_tables(double _t, const char* _nm, double* _pObs, double* _pPred) = 0;

    virtual void restart_subject_sequences() = 0;

    virtual void set_time_of_state(double _t) = 0;
    virtual void set_time_of_getval(double _t) = 0;

    virtual void get_ss(const unsigned char * pc, int * pbOK) = 0;
};

ModelCallbackTable* _pcb = 0;

static double t;
static double WT;
static double _initial_WT;
static double AaDose;
static double AaInfDose;
static double AaInfRate;

static bool reset_thetas_to_init = true;

bool HasDelay()
{
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


template <typename NUM>
struct model_data :
	base_data
{
	typedef NUM type;
	static const bool is_dual;

	NUM * globalY;
	NUM * globalR;
    NUM * ranef5;

    NUM Vmax;
    NUM _initial_Vmax;
    NUM Km;
    NUM _initial_Km;
    NUM V;
    NUM _initial_V;
    NUM Ka;
    NUM _initial_Ka;
    NUM C;

    using CalcDerivPtr = void(long *, double *, NUM *, NUM *);

	model_data()
	{
        static std::vector<NUM> gY(GetNDeriv(), NUM(0));
		static std::vector<NUM> gR(GetNDeriv(), NUM(0));
        static std::vector<NUM> gRanef(GetNRanef5(), NUM(0));

        globalY = &gY[0];
        globalR = &gR[0];
        ranef5 = &gRanef[0];
	}

    void get_y(NUM y[]) const
    {
        for (int i = 0; i < NINTEG_ALL; i++)
        {
            y[i] = globalY[i];
        }
    }

    void set_y(const NUM y[])
    {
        for (int i = 0; i < NINTEG_ALL; i++)
        {
		    globalY[i] = y[i];
        }
    }

private:
	void reset_integ()
	{
		std::fill_n(globalY, GetNDeriv(), NUM(0));
        std::fill_n(zzIRate, GetGammaDelayDEStart(0) + GetGammaDelayDECount(), .0);
        globalY[GetNintegAll()] = NUM(1); // init dummy integrator for infusions
	}

    void reset_urine()
    {
        for (int i = NINTEG_ALL - NURINE - NEVENT; i < NINTEG_ALL - NEVENT; i++)
        {
            globalY[i] = NUM(0);
        }
    }

    void reset_event()
    {
        for (int i = NINTEG_ALL - NURINE; i < NINTEG_ALL; i++)
        {
            globalY[i] = NUM(0);
        }
    }

    void eval_assign(NUM * zzY);
    void eval_sparams(NUM * zzY);
    void eval_group(NUM * zzY);
    void _set_var_value(const char * _nm, double * _pv);
    void calc_secondary();

public:
    void deriv(long* _pN, double* _pt, NUM zzY[], NUM zzR[]);

    void init_ranef(const double * ran) override;
	void init_subject() override;
	void reset_subject() override;
    void end_subject(int sim) override;

    void advance(int iODELevel, double *pTimeOuter, double _dt, int bTimeAdvances) override;

    void init_delays(double _t) override;
    void init_cf(double _t) override;
    void reinit_cf(double _t) override;

	void eval_sparams() override
	{
		eval_sparams(globalY);
	}

    void get_secondary(int* _pn, double* _sec) override;

    int  get_var_value(const char * _nm, double * _pv) override;
    void set_var_value(const char * _nm, double * _pv) override;

    void schedule_bolus1(int _mark, double _tRel, int _iCpt, double _amt, int _addl) override;
    void perform_bolus1(double * _pt, int _iCpt, double _dv, double _rate) override;
    void schedule_iv1(int _mark, double _tRel, int _iCpt, double _amt, double _rate, int _addl) override;
    void perform_inf_start1(double * _pt, int _iCpt, double _dv) override;
    void perform_inf_end1(double * _pt, int _iCpt, double _dv) override;
    void schedule_bolus2(int _mark, double _tRel, int _iCpt, double _amt, int _addl) override;
    void schedule_iv2(int _mark, double _tRel, int _iCpt, double _amt, double _rate, int _addl) override;
    void perform_bolus2(double * _pt, int _iCpt, double _dv, double _rate) override;
    void perform_inf_start2(double * _pt, int _iCpt, double _dv) override;
    void perform_inf_end2(double * _pt, int _iCpt, double _dv) override;
    void add_to_cpt(int _iCpt, double _dv) override;

    void get_ll(double _dLL[], double _dGrad[], int * pnObs, double _t, const char * _nm, double _v1, int bBql) override;
    void new_get_pred(int * pnObs, double _t, const char * _nm, double _v1, bool doActions) override;

    void get_obs_vf(int* _iWhichObs, double* _vf) override;
    void end_obs_actions() override;

    void simulate(double _t, const char* _nm, double _v1, int bBql, int bSampleEpsilon, int bForVPC) override;
    void simulate_for_tables(double _t, const char* _nm, double* _pObs, double* _pPred) override;

    void restart_subject_sequences() override;

    void set_time_of_state(double _t) override;
    void set_time_of_getval(double _t) override;

    void get_all_cf_state(NUM * _y);
    void set_all_cf_state(NUM * _y);
	void clear_ss_irrelevant_state(NUM * _y);

    void get_ss(const unsigned char * pc, int * pbOK) override;
};


template<>
const bool model_data<DualNumber>::is_dual = true;
template<>
const bool model_data<double>::is_dual = false;

model_data<DualNumber> md;
model_data<double> mdo;

base_data * data = &mdo;

void switch_to_dual()
{
    data = &md;
}

void switch_to_real()
{
    data = &mdo;
}

void InitSubject()
{
    data->init_subject();
}

void ResetSubject()
{
    data->reset_subject();
}

void EndSubject(int sim)
{
    data->end_subject(sim);
}

void InitDelays(double _t)
{
    data->init_delays(_t);
}

void InitCF(double _t)
{
    data->init_cf(_t);
}

void ReInitCF(double _t)
{
    data->reinit_cf(_t);
}

void EvalSParms()
{
    data->eval_sparams();
}

void GetSecondary(int* _pn, double* _sec)
{
    data->get_secondary(_pn, _sec);
}

void ScheduleBolus1(int _mark, double _tRel, int _iCpt, double _amt, int _addl)
{
    data->schedule_bolus1(_mark, _tRel, _iCpt, _amt, _addl);
}

void ScheduleBolus2(int _mark, double _tRel, int _iCpt, double _amt, int _addl)
{
    data->schedule_bolus2(_mark, _tRel, _iCpt, _amt, _addl);
}

void ScheduleIV1(int _mark, double _tRel, int _iCpt, double _amt, double _rate, int _addl)
{
    data->schedule_iv1(_mark, _tRel, _iCpt, _amt, _rate, _addl);
}

void ScheduleIV2(int _mark, double _tRel, int _iCpt, double _amt, double _rate, int _addl)
{
    data->schedule_iv2(_mark, _tRel, _iCpt, _amt, _rate, _addl);
}

void PerformBolus1(double * _pt, int _iCpt, double _dv, double _rate)
{
    data->perform_bolus1(_pt, _iCpt, _dv, _rate);
}

void PerformBolus2(double * _pt, int _iCpt, double _dv, double _rate)
{
    data->perform_bolus2(_pt, _iCpt, _dv, _rate);
}

void PerformInfStrt1(double * _pt, int _iCpt, double _dv)
{
    data->perform_inf_start1(_pt, _iCpt, _dv);
}

void PerformInfEnd1(double * _pt, int _iCpt, double _dv)
{
    data->perform_inf_end1(_pt, _iCpt, _dv);
}

void PerformInfStrt2(double * _pt, int _iCpt, double _dv)
{
    data->perform_inf_start2(_pt, _iCpt, _dv);
}

void PerformInfEnd2(double * _pt, int _iCpt, double _dv)
{
    data->perform_inf_end2(_pt, _iCpt, _dv);
}

void AddToCpt(int _iCpt, double _dv)
{
    data->add_to_cpt(_iCpt, _dv);
}

int GetVarValue(const char * _nm, double * _pv)
{
    return data->get_var_value(_nm, _pv);
}

void SetVarValue(const char * _nm, double * _pv)
{
    data->set_var_value(_nm, _pv);
}

void GetLL(double _dLL[], double _dGrad[], int * pnObs, double _t, const char * _nm, double _v1, int bBql)
{
    data->get_ll(_dLL, _dGrad, pnObs, _t, _nm, _v1, bBql);
}

void NewGetPred(int * pnObs, double _t, const char * _nm, double _v1, bool doActions /*= true*/)
{
    data->new_get_pred(pnObs, _t, _nm, _v1, doActions);
}

void GetObsVF(int* _iWhichObs, double* _vf)
{
    data->get_obs_vf(_iWhichObs, _vf);
}

void DoEndObsActions()
{
    data->end_obs_actions();
}

void Simulate(double _t, const char* _nm, double _v1, int bBql, int bSampleEpsilon, int bForVPC)
{
    data->simulate(_t, _nm, _v1, bBql, bSampleEpsilon, bForVPC);
}

void SimulateForSimTbl(double _t, const char* _nm, double* _pObs, double* _pPred)
{
    data->simulate_for_tables(_t, _nm, _pObs, _pPred);
}

void RestartSubjectSequences()
{
    data->restart_subject_sequences();
}

void Advance(int iODELevel, double *pTimeOuter, double _dt, int bTimeAdvances)
{
    data->advance(iODELevel, pTimeOuter, _dt, bTimeAdvances);
}

void init_ranef(const double * ran)
{
    data->init_ranef(ran);
}

void SetTimeOfStateCF(double _t)
{
    data->set_time_of_state(_t);
}

void SetTimeOfGetValCF(double _t)
{
    data->set_time_of_getval(_t);
}

template <>
void model_data<double>::init_ranef(const double * ran)
{
    for (int i = 0; i < NRANEF5; ++i)
    {
        ranef5[i] = ran[i];
    }
}

template <>
void model_data<DualNumber>::init_ranef(const double * ran)
{
    for (int i = 0; i < DualNumber::ngrad; ++i)
    {
        ranef5[i] = DualNumber(ran[i]);
        ranef5[i].b[i] = 1;
    }
}

void DerivDual(long * _pN, double * _pt, DualNumber zzY[], DualNumber zzR[]);

static void cdverk(double * ptOuter, double dt, DualNumber zzy[])
{
    if (dt <= 0)
    {
        return;
    }

    const int neq = GetNintegAll() + 1 + GetGammaDelayDECount();

    static int old_m = -1;
    static double old_dt = -1;

    static std::vector<DualNumber> rwork(9 * neq, DualNumber());

    double t = *ptOuter;
    double tout = t + dt;

	extern int nMXSTEP;
	extern double dOdeRTol;

    int res = cdverk(neq, DerivDual, t, zzy, tout, nMXSTEP, dOdeRTol, neq, &rwork[0]);
}

static void ode_solve(int iODELevel, double * pTimeOuter, double _dt, int bTimeAdvances, double zzY[])
{
    ODESolve(iODELevel, pTimeOuter, _dt, bTimeAdvances, zzY);
}

static void ode_solve(int iODELevel, double * pTimeOuter, double _dt, int bTimeAdvances, DualNumber zzY[])
{
    cdverk(pTimeOuter, _dt, zzY);
}

template <typename NUM>
static NUM vdis(int n, const NUM a[], const NUM b[])
{
	NUM sum(0.0);

    for (int i = 0; i < n; i++)
    {
		NUM del = a[i] - b[i];
        sum += del * del;
    }

    return sqrt(sum);
}

bool _extraSS = false;

template <typename NUM>
void model_data<NUM>::get_ss(const unsigned char * pc, int * pbOK)
{
    static const size_t nIntegAll = NINTEG_ALL;
    static const size_t nCFMacro = NCFMACRO;

    static const size_t ARR_SIZE = nIntegAll + nCFMacro * 6 + 1;

    static NUM zvec[ARR_SIZE];
    static NUM y0[ARR_SIZE];
    static NUM y1[ARR_SIZE];
    static NUM y2[ARR_SIZE];
    static NUM y3[ARR_SIZE];
    static NUM vec[ARR_SIZE];

	NUM dy01, dy02, dy1z;
	NUM _a, _absa;

    double _tol = 1e-5;
    bool changeTime = HasDelay();
    int i;
    int j;
    int bTrySimpleWay = changeTime;
    int bSSFailed = 0;
    SaveEventQueue();
    double tempTime = curTime;

    // first run two cycles to get in the ballpark
    if (pbOK && *pbOK == 0)
    {
        return;
    }

    curTime = 0;

    new_ss();

    InterpSS(pc, pbOK, changeTime, TRUE);

    if (pbOK && *pbOK == 0)
    {
        return;
    }

    InterpSS(pc, pbOK, changeTime, FALSE);

    if (pbOK && *pbOK == 0)
    {
        return;
    }

    for (j = 0; j < 20 && !bTrySimpleWay; j++)
    {
		get_y(y0);

		get_all_cf_state(y0 + nIntegAll);
		clear_ss_irrelevant_state(y0);

        InterpSS(pc, pbOK, changeTime, FALSE);

        if (pbOK && *pbOK == 0)
        {
            return;
        }

		get_y(y1);

		get_all_cf_state(y1 + nIntegAll);
		clear_ss_irrelevant_state(y1);

        dy01 = vdis(nIntegAll + nCFMacro * 6, y0, y1);
        dy1z = vdis(nIntegAll + nCFMacro * 6, y1, zvec);

        // check for convergence
        if ((NUM)dy01 < _tol * dy1z)
        {
            break;
        }

        InterpSS(pc, pbOK, changeTime, FALSE);

        if (pbOK && *pbOK == 0)
        {
            return;
        }

		get_y(y2);

		get_all_cf_state(y2 + nIntegAll);
		clear_ss_irrelevant_state(y2);

        dy02 = vdis(nIntegAll + nCFMacro * 6, y0, y2);

        // guard against zero-divide
        if (2 * dy01 - dy02 == 0)
        {
            bTrySimpleWay = 1;
            break;
        }

        // get vector y2 - y0
        for (i = 0; i < nIntegAll + nCFMacro * 6; i++)
        {
            vec[i] = (y2[i] - y0[i]);
        }

        // extrapolate
        _a = (dy01 * dy01) / (2 * dy01 - dy02);
        // if it's suspiciously large, give up extrapolation
        _absa = fabs(_a);

        if (_absa > 5 * dy02)
        {
            bTrySimpleWay = 1;
            break;
        }

        // do the extrapolation
        for (i = 0; i < nIntegAll + nCFMacro * 6; i++)
        {
            y3[i] = y0[i] + vec[i] / dy02 * _a;
        }

        // if any closed-form component is negative, give up extrapolation
        for (i = nIntegAll; i < nIntegAll + nCFMacro * 6; i++)
        {
            if (y3[i] < 0)
            {
                bTrySimpleWay = 1;
                break;
            }
        }

		set_y(y3);

		set_all_cf_state(y3 + nIntegAll);
    }

    if (j >= 20)
    {
        bTrySimpleWay = 1;
    }

    if (bTrySimpleWay)
    {
        bSSFailed = 0;

		get_y(y0);

		get_all_cf_state(y0 + nIntegAll);
		clear_ss_irrelevant_state(y0);

        // allow it to cycle up to this number of times before giving up
        // in case it's a really bad but still possible case
        for (j = 1000; --j >= 0;)
        {
            if (pbOK && *pbOK == 0)
            {
                return;
            }

            InterpSS(pc, pbOK, changeTime, FALSE);

            if (pbOK && *pbOK == 0)
            {
                return;
            }

			get_y(y1);

			get_all_cf_state(y1 + nIntegAll);
			clear_ss_irrelevant_state(y1);

            dy01 = vdis(nIntegAll + nCFMacro * 6, y0, y1);
            dy1z = vdis(nIntegAll + nCFMacro * 6, y1, zvec);

            // check for convergence
            if ((NUM)dy01 < _tol * dy1z)
            {
                break;
            }

            for (i = 0; i < nIntegAll + nCFMacro * 6; i++)
            {
                y0[i] = y1[i];
            }
        }

        if (j < 0)
        {
            bSSFailed = 1;
        }
    }

    if (HasDelay())
    {
        _extraSS = true;

        InterpSS(pc, pbOK, changeTime, FALSE);

        _extraSS = false;
    }

    // it's a choice between failing and stopping with a result
    // that is pretty good but doesn't meet the tolerance
    // here I'm choosing to accept the pretty-good result
    if (bSSFailed) {
        SetErrorMessageWithParams("Error: unable to achieve steady-state.");

        if (pbOK != 0L)
        {
            (*pbOK) = 0;
        }
    }

    curTime = tempTime;

    RestoreEventQueue();
}

void GetSS(const unsigned char * pc, int * pbOK)
{
	data->get_ss(pc, pbOK);
}


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
    if (*_pi == 0){
        *_pflags = _TYP_CONTIN | _SWEEP_AVAIL | _CAN_OVERLAY | _CAN_TRELLIS | _CAN_LOG | _POINT_AVAIL;
        strcpy(_nmObs, "CObs");
        strcpy(_nmPred, "C");
        return;
    }
    if (*_pi == 1){
        *_pflags = _TYP_CONTIN | _SWEEP_AVAIL | _CAN_OVERLAY | _CAN_TRELLIS | _CAN_LOG;
        strcpy(_nmObs, "CObs");
        strcpy(_nmPred, "CMultStdev");
        return;
    }
    if (*_pi == 2){
        *_pflags = _TYP_CONTIN | _SWEEP_AVAIL | _POINT_AVAIL | _CAN_OVERLAY | _CAN_TRELLIS;
        strcpy(_nmObs, "");
        strcpy(_nmPred, "WT");
        return;
    }
    if (*_pi == 3){
        *_pflags = _SWEEP_AVAIL | _CAN_OVERLAY | _CAN_TRELLIS;
        strcpy(_nmObs, "");
        strcpy(_nmPred, "Vmax");
        return;
    }
    if (*_pi == 4){
        *_pflags = _SWEEP_AVAIL | _CAN_OVERLAY | _CAN_TRELLIS;
        strcpy(_nmObs, "");
        strcpy(_nmPred, "Km");
        return;
    }
    if (*_pi == 5){
        *_pflags = _SWEEP_AVAIL | _CAN_OVERLAY | _CAN_TRELLIS;
        strcpy(_nmObs, "");
        strcpy(_nmPred, "V");
        return;
    }
    if (*_pi == 6){
        *_pflags = _SWEEP_AVAIL | _CAN_OVERLAY | _CAN_TRELLIS;
        strcpy(_nmObs, "");
        strcpy(_nmPred, "Ka");
        return;
    }
    if (*_pi == 7){
        *_pflags = _SWEEP_AVAIL | _CAN_OVERLAY | _CAN_TRELLIS;
        strcpy(_nmObs, "");
        strcpy(_nmPred, "t");
        return;
    }
}

int GetNPlotVar(){return 8;}

void AddToCpt(int _iCpt, double _dv);	// forward declaration

int GetNumCovariates(){
	return 1;
}

int GetCovariateNumber(const char* _nm){
    if (_nm[0]=='W' && strcmp(_nm, "WT")==0) return 0;

    return -1;
}

const char * GetCovariateName(int i)
{
    if (i >= 1) return "";

    static const char * names[1] = {
        "WT",
    };

    return names[i];
}

int GetCovariateDirection(int i){
    // 0=backward, 1=forward, 2=interpolate
    if (i >= 1) return -1;

    static int dirs[1] = {
        1,
    };

    return dirs[i];
}

void SetCovariate(double _t, const char* _nm, double _v){
    if (_nm[0]=='W' && strcmp(_nm, "WT")==0){
        if (WT == _NA){
            _initial_WT = _v;
        }
        WT = _v;
        return;
    }
}

void SetCovariateInterp(double _t, const char* _nm, double _v1, double _dt, double _v2, int _bForward){
    if (_nm[0]=='W' && strcmp(_nm, "WT")==0){
        WT = (_bForward ? _v1 : _v2);
        return;
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

template <typename NUM>
void model_data<NUM>::get_all_cf_state(NUM * _y)
{
}

template <typename NUM>
void model_data<NUM>::set_all_cf_state(NUM * _y)
{
}

template <typename NUM>
void model_data<NUM>::clear_ss_irrelevant_state(NUM * _y)
{
}

template <typename NUM>
void model_data<NUM>::set_time_of_state(double _t)
{
}

template <typename NUM>
void model_data<NUM>::set_time_of_getval(double _t)
{
}

template <typename NUM>
void model_data<NUM>::reset_subject()
{
	ResetTime();
	ClearActions();
	reset_integ();
	ResetCF();
	reset_urine();
	reset_event();
}

template <typename NUM>
void model_data<NUM>::init_delays(double _t)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
}

template <typename NUM>
void model_data<NUM>::init_cf(double _t)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	#if defined(USE_SETTIMECF)
		SetTimeOfStateCF(_t);
		SetTimeOfGetValCF(_t);
	#endif
}
// routine to capture CF state, change covariates, and re-initialize CF machines
// make sure all covariate-derived values are up to date before calling
template <typename NUM>
void model_data<NUM>::reinit_cf(double _t)
{
	NUM _state[NCFMACRO*6+1];
	get_all_cf_state(_state); // get the compartment values
	init_cf(_t); // regenerate all the As
	set_all_cf_state(_state); // put the values back in all the compartments
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

template <typename NUM>
void model_data<NUM>::init_subject()
{
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

template <typename NUM>
void model_data<NUM>::restart_subject_sequences()
{
}

void RestartSubjectEventSims(){
}

template <typename NUM>
void model_data<NUM>::end_subject(int bSimulating)
{
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
	for (_i=0; _i<1; _i++) _omat[_i] = 0;
	_omat[0+0+0*1] = 0.5470000000000000417;
}

double GetEtaVariance5(int _i){
	double _temp = 1;
	if (0){
	} else if (_i == 0){
		_temp = 0.5470000000000000417;
	}
	return _temp;
}

void TruncateOmega5(double* _omat){
	int _i = 0;
	double _temp = 0;
	double _mtemp[NRANEF5*NRANEF5+1];
	memcpy(_mtemp, _omat, sizeof(double)*NRANEF5*NRANEF5);
	memset(_omat, 0, sizeof(double)*NRANEF5*NRANEF5);
	_temp = _mtemp[0+0+0*1];
	_omat[0+0+0*1] = _temp;
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
	*_nPack = _n;
}

// Caution: live portion of omega matrix is assumed to be packed in lower cholesky form
void   GetNPackOmega(int* _nPack){
	int _n = 0;
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
		*_nEta = 1;
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
}

static double _dMainSigma;
void SetMainSigma(double* _sig){
	_dMainSigma = *_sig;
}

static temp_arg _tempResult[1];
static temp_arg _tempArg[1];
int _bTempValid = 0;
void InitTemp()
{
	for (int _i = 0; _i < 1; _i++) _tempResult[_i].reset();
	for (int _i = 0; _i < 1; _i++) _tempArg[_i].reset();
	_bTempValid = 1;
}
template <typename NUM>
void model_data<NUM>::eval_group(NUM * zzY)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
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

template <typename NUM>
void model_data<NUM>::eval_assign(NUM * zzY)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	int _ix0 = 0, _ix1 = 0;	// used for table statement
	double _xFrac = 0;
	if (!_bTempValid) InitTemp();
	C = (A1 / V);
}
template <typename NUM>
void model_data<NUM>::eval_sparams(NUM * zzY)
{
	if (zzY == 0)
    {
        return;
    }
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	int _ix0 = 0, _ix1 = 0;	// used for table statement
	double _xFrac = 0;
	if (!_bTempValid) InitTemp();
	eval_group(zzY);
	if (WT == _NA){
		SetErrorMsg("Covariate 'WT' not set");
		NotOK();
	}
	Vmax = (tvVmax * _exp_cached(nVmax, _tempArg[0], _tempResult[0]));
	if (_initial_Vmax == _NA) _initial_Vmax = Vmax;
	Km = tvKm;
	if (_initial_Km == _NA) _initial_Km = Km;
	V = tvV;
	if (_initial_V == _NA) _initial_V = V;
	Ka = tvKa;
	if (_initial_Ka == _NA) _initial_Ka = Ka;
	eval_assign(zzY);
}
template <typename NUM>
void model_data<NUM>::calc_secondary()
{
	double* _fixef = _args[0];
	NUM * zzY = globalY;

}

template <typename NUM>
void model_data<NUM>::get_secondary(int* _pn, double* _sec)
{
	*_pn = 0;
}

void SetModelTime(double _t){
	t = _t;
}

static int GetGammaDelayDECnt(int i_)
{
    return 0;
}

constexpr int GetGammaDelayDEStart(int i_)
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

template <typename NUM>
void model_data<NUM>::deriv(long* _pN, double* _pt, NUM zzY[], NUM zzR[])
{
	int _i = 0;
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	t = *_pt;
	#if defined(USE_SETTIMECF)
		SetTimeOfGetValCF(*_pt);
	#endif
	#if defined(THERE_ARE_INTERPOLATED_COVARIATES)
		eval_sparams(zzY);
	#else
		eval_assign(zzY);
	#endif
	zzR_Aa = ( -(Ka) * Aa);
	zzR_Aa += _Aa_irate;
	zzR_A1 = ((Ka * Aa) - ((Vmax * C) / (Km + C)));
	zzR_A1 += _A1_irate;
	zzR[2] = 0; // extra null equation
}

void Deriv(long * _pN, double * _pt, double zzY[], double zzR[])
{
    mdo.deriv(_pN, _pt, zzY, zzR);
}

void DerivDual(long * _pN, double * _pt, DualNumber zzY[], DualNumber zzR[])
{
    md.deriv(_pN, _pt, zzY, zzR);
}

	//  -(tvKa)
	// tvKa
	// 0.0
	//  -(((((tvKm + (A1 / tvV)) * ((tvVmax * exp(nVmax)) / tvV)) - (((tvVmax * exp(nVmax)) * (A1 / tvV)) / tvV)) / ((tvKm + (A1 / tvV)) * (tvKm + (A1 / tvV)))))  is not constant
	// 0.0
	// 0.0
void Jacobian(int* pN, double* pt, double zzY[], int* pML, int* pMU, double* _jacobian, int* pNRJ){
	typedef double NUM;

	double* _fixef = _args[0];
	const double * _ranef5 = _args[5];
	int _nrj = *pNRJ;
	double __temp[1];
	__temp[0] = exp(nVmax);
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
	_jacobian[1 + _nrj*1] =  -(((((tvKm + (A1 / tvV)) * ((tvVmax * __temp[0]) / tvV)) - (((tvVmax * __temp[0]) * (A1 / tvV)) / tvV)) / ((tvKm + (A1 / tvV)) * (tvKm + (A1 / tvV)))));
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

template <typename NUM>
void model_data<NUM>::advance(int iODELevel, double *pTimeOuter, double _dt, int bTimeAdvances)
{
	NUM * zzY = globalY;

	#if defined(USE_SETTIMECF)
		SetTimeOfStateCF(*pTimeOuter);
	#endif
	// truncate the level
	if (iODELevel > ODELEVEL_MAX)
    {
        iODELevel = ODELEVEL_MATEXP;
    }

	if (iODELevel < 1)
    {
        iODELevel = ODELEVEL_NONSTIFF;
    }

	// don't have constant analytic jacobian: matexp defaults to RK
	if (iODELevel == ODELEVEL_MATEXP || iODELevel == ODELEVEL_HIGHAM)
    {
        iODELevel = ODELEVEL_NONSTIFF;
    }

	ode_solve(iODELevel, pTimeOuter, _dt, bTimeAdvances, zzY);

	eval_sparams(zzY);
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
template <typename NUM>
void model_data<NUM>::schedule_bolus1(int _mark, double _tRel, int _iCpt, double _amt, int _addl)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
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

template <typename NUM>
void model_data<NUM>::schedule_iv1(int _mark, double _tRel, int _iCpt, double _amt, double _rate, int _addl)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
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

template <typename NUM>
void model_data<NUM>::perform_bolus1(double* _pt, int _iCpt, double _dv, double _rate)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	NUM * zzR = globalR;
	NUM * zzY = globalY;

	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	} else if (_iCpt==0){
		zzY[_iCpt] += _dv;
	}
}

template <typename NUM>
void model_data<NUM>::perform_inf_start1(double* _pt, int _iCpt, double _dv)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	NUM * zzR = globalR;
	NUM * zzY = globalY;

	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	} else if (_iCpt==0){
		zzIRate[_iCpt] += _dv;
	}
}

template <typename NUM>
void model_data<NUM>::perform_inf_end1(double* _pt, int _iCpt, double _dv)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	NUM * zzR = globalR;
	NUM * zzY = globalY;

	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	} else if (_iCpt==0){
		zzIRate[_iCpt] += - _dv;
	}
}
template <typename NUM>
void model_data<NUM>::schedule_bolus2(int _mark, double _tRel, int _iCpt, double _amt, int _addl)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	double _tlag = 0;
	double _rate = 1, _dt = 1, _bio = 1;
	if (0){
	}
}

template <typename NUM>
void model_data<NUM>::schedule_iv2(int _mark, double _tRel, int _iCpt, double _amt, double _rate, int _addl)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	double _tlag = 0, _bio = 1, _dt;
	if (_rate == 0) return;
	_dt = _amt / _rate;
	if (0){
	}
}

template <typename NUM>
void model_data<NUM>::perform_bolus2(double* _pt, int _iCpt, double _dv, double _rate)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	NUM * zzR = globalR;
	NUM * zzY = globalY;

	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	}
}

template <typename NUM>
void model_data<NUM>::perform_inf_start2(double* _pt, int _iCpt, double _dv)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	NUM * zzR = globalR;
	NUM * zzY = globalY;

	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	}
}

template <typename NUM>
void model_data<NUM>::perform_inf_end2(double* _pt, int _iCpt, double _dv)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	NUM * zzR = globalR;
	NUM * zzY = globalY;

	long _neq = 0;
	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	}
}
template <typename NUM>
void model_data<NUM>::add_to_cpt(int _iCpt, double _dv)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	NUM * zzY = globalY;

	if (*_piState==ODESTATE_UNDERWAY) *_piState = ODESTATE_CHANGE;
	if (0){
	} else if (_iCpt==0){
		zzY[_iCpt] += _dv;
	} else if (_iCpt==1){
		zzY[_iCpt] += _dv;
	}
}

template <typename NUM>
int model_data<NUM>::get_var_value(const char* _nm, double* _pv)
{
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	NUM * zzY = globalY;

	eval_sparams(zzY);

    if (_nm[0]=='t' && strcmp(_nm, "t")==0) {
        *_pv = t;
        return 1;
    }
    if (_nm[0]=='W' && strcmp(_nm, "WT")==0) {
        *_pv = double(WT);
        return 1;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaDose")==0) {
        *_pv = double(AaDose);
        return 1;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaInfDose")==0) {
        *_pv = double(AaInfDose);
        return 1;
    }
    if (_nm[0]=='A' && strcmp(_nm, "AaInfRate")==0) {
        *_pv = double(AaInfRate);
        return 1;
    }
    if (_nm[0]=='A' && strcmp(_nm, "Aa")==0) {
        *_pv = double(Aa);
        return 1;
    }
    if (_nm[0]=='A' && strcmp(_nm, "A1")==0) {
        *_pv = double(A1);
        return 1;
    }
    if (_nm[0]=='C' && strcmp(_nm, "CMultStdev")==0) {
        *_pv = double(CMultStdev);
        return 1;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvVmax")==0) {
        *_pv = double(tvVmax);
        return 1;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvKm")==0) {
        *_pv = double(tvKm);
        return 1;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvV")==0) {
        *_pv = double(tvV);
        return 1;
    }
    if (_nm[0]=='t' && strcmp(_nm, "tvKa")==0) {
        *_pv = double(tvKa);
        return 1;
    }
    if (_nm[0]=='n' && strcmp(_nm, "nVmax")==0) {
        *_pv = double(nVmax);
        return 1;
    }
    if (_nm[0]=='V' && strcmp(_nm, "Vmax")==0) {
        *_pv = double(Vmax);
        return 1;
    }
    if (_nm[0]=='K' && strcmp(_nm, "Km")==0) {
        *_pv = double(Km);
        return 1;
    }
    if (_nm[0]=='V' && strcmp(_nm, "V")==0) {
        *_pv = double(V);
        return 1;
    }
    if (_nm[0]=='K' && strcmp(_nm, "Ka")==0) {
        *_pv = double(Ka);
        return 1;
    }
    if (_nm[0]=='C' && strcmp(_nm, "C")==0) {
        *_pv = double(C);
        return 1;
    }
    if (_nm[0]=='t' && strcmp(_nm, "t")==0) {
        *_pv = double(t);
        return 1;
    }
	return 0;
}

template <typename NUM>
void model_data<NUM>::_set_var_value(const char* _nm, double* _pv)
{
	double* _fixef = _args[0];
	NUM * _ranef5 = ranef5;

	NUM * zzY = globalY;

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

template <typename NUM>
void model_data<NUM>::set_var_value(const char* _nm, double* _pv)
{
	NUM * zzY = globalY;

    eval_sparams(zzY);

    _set_var_value(_nm, _pv);

    eval_sparams(zzY);
}

int IsSaemRanef(int * idx)
{
	static const int is_saem[] = {1};
    return is_saem[* idx];
}

void GetSaemCMatrixSize(int* nBaseEta, int* nExtraEta, int* nTheta){
	*nBaseEta = *nTheta = *nExtraEta = 0;

	(*nTheta)++; // tvVmax
	(*nBaseEta)++; // nVmax

	if ((-1 == -1)||_ISENABLED(-1)){(*nExtraEta)++; (*nTheta)++;} // bare theta CMultStdev

	if ((-1 == -1)||_ISENABLED(-1)){(*nExtraEta)++; (*nTheta)++;} // bare theta tvKm

	if ((-1 == -1)||_ISENABLED(-1)){(*nExtraEta)++; (*nTheta)++;} // bare theta tvV

	if ((-1 == -1)||_ISENABLED(-1)){(*nExtraEta)++; (*nTheta)++;} // bare theta tvKa

}
void GetSaemCMatrix(double* daCMat, int* theta2eta){
	double* _fixef = _args[0];
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

	// for additional bare thetas
	{
		daCMat[_iEta + _jTheta * _nEta] = 1; // extra eta for bare theta CMultStdev
		theta2eta[_jTheta] = -1;
		_iEta++;
		_jTheta++;
	}

	{
		daCMat[_iEta + _jTheta * _nEta] = 1; // extra eta for bare theta tvKm
		theta2eta[_jTheta] = -1;
		_iEta++;
		_jTheta++;
	}

	{
		daCMat[_iEta + _jTheta * _nEta] = 1; // extra eta for bare theta tvV
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

	return 0;
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

	// for additional bare thetas
	{
		iaCColToITheta[_jCTheta] = _iaIFixefErrToITheta[0]; // CMultStdev
		baCColIsLog[_jCTheta] = 0;
		_jCTheta++;
	}

	{
		iaCColToITheta[_jCTheta] = _iaIFixefErrToITheta[2]; // tvKm
		baCColIsLog[_jCTheta] = 0;
		_jCTheta++;
	}

	{
		iaCColToITheta[_jCTheta] = _iaIFixefErrToITheta[3]; // tvV
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
template <typename NUM>
void model_data<NUM>::end_obs_actions()
{
	NUM * zzY = globalY;

}
template <typename NUM>
void model_data<NUM>::new_get_pred(int* pnObs, double _t, const char* _nm, double _v1, bool doActions /*= true*/)
{
	NUM * zzY = globalY;

	double _dVFCorrectionFactor = 1;
	const double t = _t;
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	NUM _ypred(0); // predicted value
	NUM _var(1);
	NUM _sig(1); // current standard deviation of Epsilon
	double _std = 1; // current standard deviation of Epsilon
	NUM _sqrtvf(1); // dObs / dEpsilon (additive = 1, proportional = _ypred)
	double _vf = 1; // _sqrtvf*_sqrtvf (additive = 1, proportional = _ypred^2)

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
			RecordObservation(t, curTimeLastDose, iCurWhichDose, iCurWhichReset, logfl((double)_ypred), logfl(_v1), 1, 0);
		} else {
			RecordObservation(t, curTimeLastDose, iCurWhichDose, iCurWhichReset, (double)_ypred, _v1, (double)_var, 0);
		}

		(*pnObs)++; // count all observations

        return;
    }
}

const char * GetGradientErrorMessage()
{
	return "";
}

// jacobian info: 2=can matexp, 1=can stiffjac, 0=no jac, -1=no odes
void GetJacobianInfo(int* bUsingSyntheticGradients, int* _piJacobianInfo){
	if (*bUsingSyntheticGradients){
		* _piJacobianInfo = 0;
	} else {
		* _piJacobianInfo = 1;
	}
}
void GetMaxODELevel(int* bUsingSyntheticGradients, int* _piMaxODELevel){
	if (*bUsingSyntheticGradients){
		* _piMaxODELevel = ODELEVEL_STIFF;
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

template <typename NUM>
void model_data<NUM>::get_ll(double _dLL[], double _dGrad[], int* pnObs, double _t, const char* _nm, double _v1, int bBql)
{
	NUM * zzY = globalY;

	NUM _dLL1(0);
	double _probge[12];

    const double t = _t;
    double* _fixef = _args[0];
    const NUM * _ranef5 = ranef5;
    NUM _ypred(0), _var(1), _sig(1), _sqrtvf(1);
    double _std(1);
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
		_dLL[0] += (double)_dLL1;

		copy_grad(_dLL1, _dGrad);

		(*pnObs)++;
		if (_pcb){
			if (ALLOWLOGTRANSFORM && bEnableLogTransform){
				_pcb->OnObserve(curTime, curTimeLastDose, _nm, -1, logfl((double)_ypred), 0, logfl(_v1), (double)_sig, (double)_sqrtvf, bBql);
			} else {
				_pcb->OnObserve(curTime, curTimeLastDose, _nm, -1, (double)_ypred, 0, _v1, (double)_sig, (double)_sqrtvf, bBql);
			}
		}

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
	if (_iError != _iLastFreeError){
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
template <typename NUM>
void model_data<NUM>::get_obs_vf(int* _iWhichObs, double* _vf)
{
	#define t (_t + 0)
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	NUM * zzY = globalY;

	*_vf = 1;
	if (0);
	else if (*_iWhichObs == 0){
		double _std = _fixef[5];
		*_vf = (double)(sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))) * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))));
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

int GetVcv(const char **& names, const double *& vals, const int *& sizes)
{
    return 0;
}

template <typename NUM>
void model_data<NUM>::simulate(double _t, const char* _nm, double _v1, int bBql, int bSampleEpsilon, int bForVPC)
{
	NUM * zzY = globalY;

	double _probge[12];

    const double t = _t;
    double* _fixef = _args[0];
    const NUM * _ranef5 = ranef5;
    NUM _ypred(0), _var(1), _sig(1), _sqrtvf(1);
    double _std(1);
    int _i = 0;
	eval_sparams(zzY);

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
		if (!isfinite(_sig)){
			SetErrorMsg("Error (Model.cpp): non-finite _sig");
		}
		if (_sig == 0){
			SetErrorMsg("Error (Model.cpp): zero _sig");
		}
		_var = _sig * _sig;
		if (_pcb != 0){
			double _eps = (bSampleEpsilon ? _pcb->GetNormal() : 0); // get a N(0,1)
			NUM _ypred1(0);
			_eps *= _fixef[5];
			_ypred1 = (C + (_eps * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0))))));
			if (ALLOWLOGTRANSFORM && bEnableLogTransform){
				_pcb->OnObserve(curTime, curTimeLastDose, _nm, -1, logfl((double)_ypred), logfl((double)_ypred1), logfl(_v1), (double)_sig, (double)_sqrtvf, bBql);
			} else {
				_pcb->OnObserve(curTime, curTimeLastDose, _nm, -1, (double)_ypred, (double)_ypred1, _v1, (double)_sig, (double)_sqrtvf, bBql);
			}
		}

        return;
    }
}

void GetNumCategories(const char* _nm, int* pNCat){
	*pNCat = 0;
}

template <typename NUM>
void model_data<NUM>::simulate_for_tables(double _t, const char* _nm, double* _pObs, double* _pPred)
{
	NUM * zzY = globalY;

	const double t = _t;
	double* _fixef = _args[0];
	const NUM * _ranef5 = ranef5;
	double _eps = 0;
	eval_sparams(zzY);
	if (_pcb == 0) return;

    if (_nm[0]=='C' && strcmp(_nm, "CObs")==0){
        double _std = 0;
        _std = _fixef[5];
        // MD 151026: Tentative change: Turn off log transformation for sim tables
        if (FALSE && ALLOWLOGTRANSFORM && bEnableLogTransform){
            *_pPred = logfl((double)(C + (_eps * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))))));
            _eps = _pcb->GetNormal() * _fixef[5]; // get a N(0,stdev)
            *_pObs = logfl((double)(C + (_eps * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0)))))));
        } else {
            *_pPred = (double)(C + (_eps * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0))))));
            _eps = _pcb->GetNormal() * _fixef[5]; // get a N(0,stdev)
            *_pObs = (double)(C + (_eps * sqrt((1.0 + (pow(C, 2.0) * pow((CMultStdev / sigma()), 2.0))))));
        }

        return;
    }
}


