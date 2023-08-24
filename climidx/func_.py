"""
>--#########################################################################--<
>--------------------Calculations of Some Climate Indices---------------------<
>------------------------------by Changgui Lin--------------------------------<
>----------------------------------at SMHI------------------------------------<
>--------------------------------07 Feb 2020--------------------------created-<
>--------------------------------07 Feb 2020-------------------------modified-<
>--#########################################################################--<
* mEffPr_                     : cube, pr, evspsbl, climatology (season, year)
* mPrrn_                      : cube, pr, prsn, climatology (season, year)
* mDTR_                       : cube, tasmax, tasmin
* dHumiWarmDays_              : 1d arr, hurs, tas
* dHumiWarmDays_cube          : cube, hurs, tas
* dPr7_                       : 1d arr, pr
* dPr7_cube                   : cube, pr
* dLongestDryDays_            : 1d arr, pr
* dDryDays_                   : cube, pr
* dExtrPrDays_                : cube, pr
* dWarmDays_                  : cube, tasmax
* dCalmDays_                  : cube, sfcWind
* dConWarmDays_               : 1d arr, tasmax
* dConCalmDays_               : 1d arr, sfcWind
* dColdDays_                  : cube, tasmax
* dDegDay20_                  : 1d arr, tas
* dDegDay8_vegSeason_         : 1d arr, tas
* dDegDay17_                  : 1d arr, tas
* dEndSpringFrost_            : 1d arr, tasmin
* dFirstDayWithoutFrost_      : 1d arr, tasmin
* dFrostDays_                 : cube, tasmin
* dMinusDays_                 : cube, tas
* dFreezingDays_              : cube, tas
* dTropicNights_              : cube, tasmin
* dZeroCrossingDays_          : 1d arr, tasmax, tasmin
* dZeroCrossingDays_cube      : cube, tasmax, tasmin
* dStartEndVegSeason_         : 1d arr, tas, thr=5/2
* dWindyDays_                 : cube, wsgsmax
* dFreezRainDays_             : 1d arr, pr, ps, tas, ta925, ta850, ta700,
                                                     hus925, hus850, hus700
* dSncDays_                   : cube, snc
* dSndGT15Days_               : cube, snd
* dSndLE10Days_               : cube, snd
* dSndGT10LE20Days_           : cube, snd
* rho_fr_t_q_p_               : nd array, t, q, p
* dRhoPL_                     : 1d arr, ta, hus
* dRhoPS_                     : 1d arr, tas, huss, ps
* dDegDay_                    : cube, tas/tasmax/tasmin
* dPRSN_fr_PR_T_              : cube, pr, tas
* dPRRN_fr_PR_T_              : cube, pr, tas
* dRainOnSnow_                : cube, prrn, snc, snw
* dColdRainDays_              : cube, pr, tas
* dWarmSnowDays_              : cube, pr, tas
* dColdPRRNdays_              : cube, pr, prsn, tas
* dWarmPRSNdays_              : cube, prsn, tas
* dColdRainWarmSnowDays_      : cube, pr, tas
"""


import numpy as np
import dask.array as da
import iris
from iris.cube import Cube as _Cube
from uuuu import (
        consecutive_, rMEAN1d_,                                                #ffff
        pSTAT_cube, extract_byAxes_, rm_t_aux_cube, pst_                       #cccc
        )


__all__ = ['mEffPr_',
           'mPrrn_',
           'mDTR_',
           'dHumiWarmDays_',
           'dHumiWarmDays_cube',
           'dPr7_',
           'dPr7_cube',
           'dLongestDryDays_',
           'dDryDays_',
           'dExtrPrDays_',
           'dWarmDays_',
           'dCalmDays_',
           'dConWarmDays_',
           'dConCalmDays_',
           'dColdDays_',
           'dDegDay20_',
           'dDegDay8_vegSeason_',
           'dDegDay17_',
           'dEndSpringFrost_',
           'dFirstDayWithoutFrost_',
           'dFrostDays_',
           'dMinusDays_',
           'dFreezingDays_',
           'dTropicNights_',
           'dZeroCrossingDays_',
           'dZeroCrossingDays_cube',
           'dStartEndVegSeason_',
           'dWindyDays_',
           'dFreezRainDays_',
           'dSncDays_',
           'dSndGT15Days_',
           'dSndLE10Days_',
           'dSndGT10LE20Days_',
           'rho_fr_t_q_p_',
           'dRhoPL_',
           'dRhoPS_',
           'dDegDay_',
           'dPRSN_fr_PR_T_',
           'dPRRN_fr_PR_T_',
           'dRainOnSnow_',
           'dColdRainDays_',
           'dWarmSnowDays_',
           'dColdPRRNdays_',
           'dWarmPRSNdays_',
           'dColdRainWarmSnowDays_',
           'ws_cube',
           'extract_season_cube']


K0 = 273.15


def ws_data_func(u_data, v_data):
    return da.sqrt( u_data**2 + v_data**2 )

def ws_units_func(u_cube, v_cube):
    if u_cube.units != getattr(v_cube, 'units', u_cube.units):
        raise ValueError("units do not match")
    return u_cube.units

ws_cube = iris.analysis.maths.IFunc(ws_data_func, ws_units_func)


def _rl(x):
    return x.data if isinstance(x, _Cube) else x


def mEffPr_(cPr, cET, freq):
    c = cPr.copy(cPr.data - cET.data)
    return pSTAT_cube(c, *freq)


def mPrrn_(cPr, cPrsni, freq):
    c = cPr.copy(cPr.data - cPrsn.data)
    return pSTAT_cube(c, *freq)


def mDTR_(cTasmax, cTasmin, freq):
    c = cTasmax.copy(cTasmax.data - cTasmin.data)
    return pSTAT_cube(c, *freq)


def dHumiWarmDays_(rh, tas, yrs):
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = np.sum(np.logical_and(_rl(rh[yrs == yr]) > 90,
                                     _rl(tas[yrs == yr]) > K0 + 10))
    return o


def dHumiWarmDays_cube(cHurs, cTas, freq):
    c = cTas.copy(np.logical_and(cHurs.data > 90, cTas.data > K0 + 10))
    c = pSTAT_cube(c, *freq, stat='COUNT', function=lambda cell: cell)
    return c


def dPr7_(pr, yrs):
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = np.max(rMEAN1d_(_rl(pr[yrs == yr]), 7))
    return o


def dPr7_cube(cPr, freq):
    yr_doy_cube(cPr)
    c = cPr.rolling_window('time', iris.analysis.MEAN, 7)
    c = extract_byAxes_(c, 'time', c.coord('year').points%1 == 0)
    rm_t_aux_cube(c, keep='year')
    c = pSTAT_cube(c, *freq, stat='MAX')
    return c


def dLongestDryDays_(pr, yrs):
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = consecutive_(_rl(pr[yrs == yr]),
                            lambda x: x >= 1. / (24 * 3600))
    return o


def dDryDays_(cPr, freq):
    c = pSTAT_cube(cPr, *freq,
                   stat='COUNT',
                   function=lambda cell: cell < 1. / (24 * 3600))
    return c


def dExtrPrDays_(cPr, freq, thr=10.):
    o = pSTAT_cube(cPr, *freq,
                   stat='COUNT',
                   function=lambda cell: cell > thr / (24 * 3600))
    return o


def dWarmDays_(cTasmax, freq):
    c = pSTAT_cube(cTasmax, *freq,
                   stat='COUNT',
                   function=lambda cell: cell > K0 + 20)
    return c


def dCalmDays_(cSfcWind, freq):
    c = pSTAT_cube(cSfcWind, *freq,
                   stat='COUNT',
                   function=lambda cell: cell < 2)
    return c


def dConCalmDays_(sfcw, yrs):
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = consecutive_(_rl(sfcw[yrs == yr]),
                            lambda x: x >= 2)
    return o


def dConWarmDays_(tx, yrs):
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = consecutive_(_rl(tx[yrs == yr]),
                            lambda x: x <= K0 + 20)
    return o


def dColdDays_(cTasmax, freq):
    c = pSTAT_cube(cTasmax, *freq,
                   stat='COUNT',
                   function=lambda cell: cell < K0 - 7)
    return c


def dDegDay20_(tas, yrs):
    t0 = _rl(tas.copy())
    t0 -= K0 + 20
    t0[t0 < 0] = 0
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = np.sum(t0[yrs == yr])
    return o


def dDegDay8_vegSeason_(tas, yrs):
    t0 = _rl(tas.copy())
    t0 = rMEAN1d_(t0, 4, 'same')
    t1 = _rl(tas.copy())
    t1 -= K0 + 8
    t1[t1 < 0] = 0
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = np.sum(t1[np.logical_and(yrs == yr, t0 > K0 + 5)])
    return o


def dDegDay17_(tas, yrs):
    t0 = _rl(tas.copy())
    t0 = K0 + 17 - t0
    t0[t0 < 0] = 0
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = np.sum(t0[yrs == yr])
    return o


def dFirstDayWithoutFrost_(tn, yrs, doy):
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        ts = _rl(tn[yrs == yr])
        ydoy = doy[yrs == yr]
        try:
            tmp = ydoy[ts >= K0][0]
        except IndexError:
            tmp = np.nan
        o[i] = tmp
    return o


def dEndSpringFrost_(tn, yrs, doy):
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        ind = np.logical_and(yrs == yr, doy<=213)
        ts = _rl(tn[ind])
        ydoy = doy[ind]
        try:
            tmp = ydoy[ts < K0][-1]
        except IndexError:
            tmp = np.nan
        o[i] = tmp
    return o


def dMinusDays_(cTas, freq):
    c = pSTAT_cube(cTas, *freq,
                   stat='COUNT',
                   function=lambda cell: cell < K0)
    return c


def dFreezingDays_(cTasmax, freq):
    c = pSTAT_cube(cTasmax, *freq,
                   stat='COUNT',
                   function=lambda cell: cell < K0)
    return c


def dFrostDays_(cTasmin, freq):
    c = pSTAT_cube(cTasmin, *freq,
                   stat='COUNT',
                   function=lambda cell: cell < K0)
    return c


def dTropicNights_(cTasmin, freq):
    c = pSTAT_cube(cTasmin, *freq,
                   stat='COUNT',
                   function=lambda cell: cell > K0 + 17)
    return c


def dZeroCrossingDays_cube(cTasmax, cTasmin, freq):
    ind = np.logical_and(cTasmax.data > K0, cTasmin.data < K0)
    c = pSTAT_cube(cTasmax.copy(ind), *freq,
                   stat='COUNT',
                   function=lambda cell: cell)
    return c


def dZeroCrossingDays_(tx, tn, yrs):
    tt = np.logical_and(_rl(tx) > K0, _rl(tn) < K0)
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = np.sum(tt[yrs == yr])
    return o


def dStartEndVegSeason_(tas, yrs, doy, thr=5):
    uYrs = np.unique(yrs)
    o0 = np.ma.zeros((len(uYrs),))
    o1 = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        ydoy = doy[yrs == yr]
        ts = rMEAN1d_(_rl(tas[yrs == yr]), 4)
        if np.all(ts < K0 + thr):
            o0[i], o1[i] = ydoy[0], ydoy[0]
        else:
            inds = np.where(ts >= K0 + thr)[0] + 3
            o0[i], o1[i] = ydoy[inds[0]], ydoy[inds[-1]]
    return (o0, o1)


def dWindyDays_(cWsgsmax, freq):
    c = pSTAT_cube(cWsgsmax, *freq,
                   stat='COUNT',
                   function=lambda cell: cell > 21)
    return c


def dMarkFreezRain_(pr, ps, tas, t925, t850, t700, q925, q850, q700,
                    pr_min=1.8e-5, t_cold=273.24, t_melt=272.51,
                    rh_melt=89., h_cold=6900):
    import atmos
    aa = _rl(pr) > pr_min
    bb = _rl(tas) <= t_cold
    ta = np.vstack((_rl(t925), _rl(t850), _rl(t700)))
    qa = np.vstack((_rl(q925), _rl(q850), _rl(q700)))
    prl = np.asarray([92500, 85000, 70000]).reshape((3, 1))
    rh = atmos.calculate('RH', qv=qa, T=ta, p=prl)
    cc = ta > t_melt
    dd = rh >= rh_melt
    ee = (_rl(ps) - h_cold - prl) >= 0
    ff = np.any(np.logical_and(cc, dd, ee), axis=0)
    return np.logical_and(aa, bb, ff)


def dFreezRainDays_(pr, ps, tas, t925, t850, t700, q925, q850, q700, yrs,
                    pr_min=1.8e-5, t_cold=273.24, t_melt=272.51,
                    rh_melt=89., h_cold=6900):
    tt = dMarkFreezRain_(pr, ps, tas, t925, t850, t700, q925, q850, q700,
                         pr_min=pr_min, t_cold=t_cold, t_melt=t_melt,
                         rh_melt=rh_melt, h_cold=h_cold)
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = np.sum(tt[yrs == yr])
    return o


def dSncDays_(cSnc, freq, thr=.0):
    c = pSTAT_cube(cSnc, *freq,
                   stat='COUNT',
                   function=lambda cell: cell > thr)
    sthr = 'days with snow cover over {:g}%'.format(thr)
    pst_(c, sthr, 'days')
    return c


def dSndGT15Days_(cSnd, freq):
    c = pSTAT_cube(cSnd, *freq,
                   stat='COUNT',
                   function=lambda cell: cell > .15)
    return c


def dSndLE10Days_(cSnd, freq):
    c = pSTAT_cube(cSnd, *freq,
                   stat='COUNT',
                   function=lambda cell: np.logical_and(cell > 0.,
                                                        cell <= .1))
    return c


def dSndGT10LE20Days_(cSnd, freq):
    c = pSTAT_cube(cSnd, *freq,
                   stat='COUNT',
                   function=lambda cell: np.logical_and(cell > .1,
                                                        cell <= .2))
    return c


def rho_fr_t_q_p_(t, q, p):
    import atmos
    return  atmos.calculate('rho', T=_rl(t), qv=_rl(q), p=_rl(p))


def dRhoPL_(ta_pl, hus_pl, yrs, pl=92500):
    tt = rho_fr_t_q_p_(ta_pl, hus_pl, pl)
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = np.mean(tt[yrs == yr])
    return o


def dRhoPS_(tas, huss, ps, yrs):
    tt = rho_fr_t_q_p_(tas, huss, ps)
    uYrs = np.unique(yrs)
    o = np.ma.zeros((len(uYrs),))
    for i, yr in enumerate(uYrs):
        o[i] = np.mean(tt[yrs == yr])
    return o


def dDegDay_(cT, freq, thr=0., left=False):
    c = cT.copy()
    c.data -= thr + K0
    if left:
        c.data *= -1
    c.data[c.data < 0] = 0
    c = pSTAT_cube(c, *freq, stat='SUM')
    lg_ = 'lt' if left else 'gt'
    vn_ = 'dd_{:g}' if left else 'dd{:g}_'
    pst_(c, 'degree day {}{:g}'.format(lg_, thr),
         var_name=vn_.format(thr))
    return c


def dRainOnSnow_(cPrrn, cSnc, freq,
                 cSnw=None, thr_r=5., thr_c=25., thr_w=.003, attr=None):
    if cSnw:
        c = cPrrn.copy(np.logical_and(cPrrn.data > thr_r / 3600 / 24,
                                      cSnc.data > thr_c, cSnw.data > thr_w))
    else:
        c = cPrrn.copy(np.logical_and(cPrrn.data > thr_r / 3600 / 24,
                                      cSnc.data > thr_c))
    c = pSTAT_cube(c, *freq, stat='COUNT', function=lambda cell: cell)
    if attr:
        pst_(c, attrU=dict(ROS=attr))
    return c


def dPRSN_fr_PR_T_(cPr, cTas, thr=.58):
    data = cPr.data.copy()
    data[cTas.data >= K0 + thr] = 0.
    c = cPr.copy(data)
    pst_(c, 'snowfall_flux',
         var_name='prsn',
         attrU=dict(PRSN='PRSN from p & t (.lt. {:g})'.format(thr)))
    return c


def dPRRN_fr_PR_T_(cPr, cTas, thr=.58):
    data = cPr.data.copy()
    data[cTas.data < K0 + thr] = 0.
    c = cPr.copy(data)
    pst_(c, 'rainfall_flux',
         var_name='prrn',
         attrU=dict(PRRN='PRRN from p & t (.ge. {:g})'.format(thr)))
    return c


def dColdRainWarmSnowDays_(cPr, cTas, freq):
    ind = np.logical_and(cTas.data <= K0 + 2, cTas.data >= K0 - 2)
    ind0 = np.logical_and(ind, cPr.data > 1. / 24 / 3600)
    c = pSTAT_cube(cPr.copy(ind0), *freq,
                   stat='COUNT',
                   function=lambda cell: cell)
    return c


def dColdRainDays_(cPr, cTas, freq, thr_t=.58, thr_pr=1.):
    ind = np.logical_and(cTas.data <= K0 + 2, cTas.data >= K0 + thr_t)
    ind0 = np.logical_and(ind, cPr.data > thr_pr / 24 / 3600)
    c = pSTAT_cube(cPr.copy(ind0), *freq,
                   stat='COUNT',
                   function=lambda cell: cell)
    pst_(c, attrU=dict(PRRN='p & t .ge. {:g}'.format(thr_t)))
    return c


def dWarmSnowDays_(cPr, cTas, freq, thr_t=.58, thr_pr=1.):
    ind = np.logical_and(cTas.data >= K0 - 2, cTas.data < K0 + thr_t)
    ind0 = np.logical_and(ind, cPr.data > thr_pr / 24 / 3600)
    c = pSTAT_cube(cPr.copy(ind0), *freq,
                   stat='COUNT',
                   function=lambda cell: cell)
    pst_(c, attrU=dict(PRSN='p & t .lt. {:g}'.format(thr_t)))
    return c


def dColdPRRNdays_(cPr, cPrsn, cTas, freq, thr_pr=1.):
    ind = cTas.data <= K0 + 2
    ind0 = np.logical_and(ind, cPr.data - cPrsn.data > thr_pr / 24 / 3600)
    c = pSTAT_cube(cPr.copy(ind0), *freq,
                   stat='COUNT',
                   function=lambda cell: cell)
    return c


def dWarmPRSNdays_(cPrsn, cTas, freq, thr_pr=1.):
    ind = cTas.data >= K0 - 2
    ind0 = np.logical_and(ind, cPrsn.data > thr_pr / 24 / 3600)
    c = pSTAT_cube(cPrsn.copy(ind0), *freq,
                   stat='COUNT',
                   function=lambda cell: cell)
    return c
