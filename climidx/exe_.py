import numpy as np
import iris
from abc import abstractmethod

from uuuu import (
        isIter_, isMyIter_, isMonth_, ismono_, isSeason_, ll_, l__, mmmN_,     #ffff
        prg_,                                                                  #ffff
        pSTAT_cube, axT_, seasonyr_cube, initAnnualCube_, nTslice_cube,        #cccc
        concat_cube_, yr_doy_cube, ax_fn_mp_                                   #cccc
        )

from .io_ import _SV, _TO1, _get_freq, inDF
from .group_ import idxL
from .func_ import rho_fr_t_q_p_, ws_cube


def _tt(cube, y0y1=None, mmm=None):
    yr_doy_cube(cube)
    if mmm and isSeason_(mmm):
        seasonyr_cube(cube, mmm)
        yrs = cube.coord('seasonyr').points
        y_y_ = y0y1 if y0y1 else yrs[[0, -1]]
        y_y_ = y_y_ if ismono_(mmmN_(mmm)) else [y_y_[0] + 1, y_y_[-1] - 1]
    else:
        yrs = cube.coord('year').points
        y_y_ = y0y1 if y0y1 else yrs[[0, -1]]
    doy = cube.coord('doy').points
    return (axT_(cube), y_y_, yrs, doy)


def _szG(cL):
    if not isMyIter_(cL):
        return np.prod(cL.shape) * 8 * 1.e-9
    else:
        return np.prod(cL[0].shape) * len(cL) * 8 * 1.e-9


def _f_n(func, cL, *args, xm=80, **kwargs):
    if _szG(cL) < xm:
        func(cL, *args, **kwargs)
    else:
        n = int(np.ceil(_szG(cL) / xm))
        cLL = [nTslice_cube(i, n) for i in cL]
        nn = len(cLL[0])
        ll_('n_slice = {}'.format(nn))
        t000 = l__('loop nTslice')
        for i in range(nn):
            ###H
            #if i < 5:
            #    continue
            ###I
            func([ii[i].copy() for ii in cLL], *args,
                 _nm='__{}'.format(i), **kwargs)
            ll_(prg_(i, nn), t000)


def _afm_n(cL, ax, func, out, *args, npr=32, xm=80, **kwargs):
    if _szG(cL) < xm:
        ax_fn_mp_([i.data for i in cL] if isMyIter_(cL) else cL.data,
                  ax, func, out, *args, npr=npr, **kwargs)
    else:
        n = int(np.ceil(_szG(cL) / xm))
        cLL = [nTslice_cube(i, n) for i in cL] if isMyIter_(cL) else \
              nTslice_cube(cL, n)
        outL = [nTslice_cube(i, n) for i in out] if isMyIter_(out) else \
               nTslice_cube(out, n)
        nn = len(cLL[0]) if isMyIter_(cL) else len(cLL)
        t000 = l__('loop nTslice')
        for i in range(nn):
            ax_fn_mp_([ii[i].copy().data for ii in cLL] if isMyIter_(cL) else
                      cLL[i].data,
                      ax, func,
                      [ii[i] for ii in outL] if isMyIter_(out) else outL[i],
                      *args, npr=npr, **kwargs)
            ll_(prg_(i, nn), t000)
        out_ = [concat_cube_(iris.cube.CubeList(i)) for i in outL]\
               if isMyIter_(out) else concat_cube_(iris.cube.CubeList(outL))
        return out_


def _exe(idxD,                                                                 # dictionary for climate indices
         dn=None, gwl=None, po_=None, il_=None, user_cfg=None,                 # dgpic_
         y0y1=None,                                                            # required by _d2
         **c_kwArgs):                                                          # cubes

    dgpic_ = (dn, gwl, po_, il_, user_cfg)                                     # required by _SV, _TO1
    if any([i is None for i in dgpic_[:-1]]):                                  # check mandotory input arguments
        raise ValueError("inputs of 'dn', 'gwl', 'po_', 'il_' are mandotory!")
    os.makedirs(po_, exist_ok=True)                                            # directory for output

    iSKIP = []                                                                 # index to be skiped as already computed

    def _inCL(i_):
        vL_ = idxD[i_][2] if idxD[i_][5][1] is None else idxD[i_][5][1]
        if (all("c_{}".format(i) in c_kwArgs for i in vL_) and
            all(c_kwArgs["c_{}".format(i)] is not None for i in vL)):
            return tuple(c_kwArgs["c_{}".format(i)] for i in VL)

    def _iexe(i_):
        func_ = eval('_d{}'.format(idxD[i_][5][0]))
        cL = _inCL(i_)
        if cL is not None:
            func_(i_, cL)

    def _d0(i_, cL_, freq=None, _nm=None):
        freq = freq if freq else _get_freq(idxD, i_, dgpic_[4])
        t000 = l__('{} {}'.format(idxD[i_][0], i_))
        o = pSTAT_cube(cL_[0], *freq,
                       stat=idxD[i_][3] if idxD[i_][3] else 'MEAN',
                       )
        _SV(idxD, i_, o, dgpic_, freq=freq, _nm=_nm)                           # _SV(idxD_, i_, o, dgpic_, freq=None, _nm=None)
        ll_(i_, t000)

    def _d1(i_, cL_,
            fA_=(), fK_={}, freq=None):
        if isIter_(i_):                                                        # preprocessing regarding required index
            ik_ = i_[0]
            lmsg = '/'.join(('{}',) * len(i_)) + ' {}'
            vids = [idxD_[i][0] for i in i_]
            lmsg_ = lmsg.format(*vids, i_[0][:8])
        else:
            ik_ = i_
            lmsg_ = '{} {}'.format(idxD_[i_][0], i_)

        fA_ = fA_ if fA_ else idxD_[i_][5][2]
        fK_ = fK_ if fK_ else idxD_[i_][5][3]
        freq = freq if freq else _get_freq(idxD_, ik_, dgpic_[4])

        t000 = l__(lmsg_)
        o = idxD_[ik_][3](*cL_, freq, *fA_, **fK_)
        if not isIter_(i_) or len(i_) == 1:
            _SV(idxD_, ik_, o, dgpic_, freq=freq)
        else:
            for i, ii in zip(i_, o):
                _SV(idxD_, i, ii, dgpic_, freq=freq)
        ll_(lmsg_, t000)

    def _d2(i_, CL_,
            cK_={}, fA_=(), fK_={}, freq=None, out=False):
        if isIter_(i_):                                                        # preprocessing regarding required index
            ik_ = i_[0]
            lmsg = '/'.join(('{}',) * len(i_)) + ' {}'
            vids = [idxD[i][0] for i in i_]
            lmsg_ = lmsg.format(*vids, i_[0][:8])
        else:
            ik_ = i_
            lmsg_ = '{} {}'.format(idxD[i_][0], i_)

        fA_ = fA_ if fA_ else idxD_[i_][5][2]
        fK_ = fK_ if fK_ else idxD_[i_][5][3]
        cK_ = cK_ if cK_ else idxD_[i_][5][4]
        freq = freq if freq else _get_freq(idxD, ik_, dgpic_[4])
        freq = ('djf', 'mam', 'jja', 'son') if freq == 'season' else freq

        def _inic(y_y_, icK_, mmm=None):                                       # create output cubes
            icK_cp = icK_.copy()
            if mmm:
                icK_cp.update(dict(mmm=mmm))
            return initAnnualCube_(cL_[0], y_y_, **icK_cp)

        def _cube_cp(extract_func, mmm):                                       # if extraction regarding 'freq'
            return tuple(extract_func(i, mmm) for i in cL_)

        def _run_single_freq(x):                                               # knowing x
            mmm_ = None
            if isMonth_(x):
                cube_cp = _cube_cp(extract_month_cube, x)
                mmm_ = x
            elif isSeason_(x):
                cube_cp = _cube_cp(extract_season_cube, x)
                mmm_ = x
                cube_cp = _cube_cp(extract_season_cube, x)
                mmm_ = x
            elif ff == 'year':
                cube_cp = cL_
            else:
                emsg = "currently only 'year'/month/season acceptable!"
                raise Exception(emsg)
            ax_t, y_y_, tyrs, tdoy = _tt(cube_cp[0], y0y1, mmm_)               # axis, y0y1, years of data points, doy of data points
            if isinstance(cK_, dict):
                o = _inic(y_y_, cK_, mmm_)
            else:
                o = iris.cube.CubeList([_inic(y_y_, i, mmm_) for i in cK_])
            fA = ()                                                            # aditional arguments to index function
            if 'y' in fA_:
                fA += (tyrs,)
            if 'd' in fA_:
                fA += (tdoy,)
            o_ = _afm_n(cube_cp, ax_t, idxD[ik_][3], o, *fA, **fK_)
            o = o_ if o_ else o
            if not isIter_(i_) or len(i_) == 1:
                _SV(idxD, ik_, o, dgpic_, freq=(x,))
            else:
                for i, ii in zip(i_, o):
                    _SV(idxD, i, ii, dgpic_, freq=(x,))
            return o

        t000 = l__(lmsg_)
        ooo = []
        for ff in freq:
            tmp = _run_single_freq(ff)
            ooo.append(tmp)
            ll_('   >>>>{}'.format(ff), t000)
        ll_(lmsg_, t000)
        if out:
            return ooo[0] if len(ooo) == 1 else ooo

    def _d3(i_, CL_):
        VegSeason5 = ['VegSeasonDayStart-5',
                      'VegSeasonDayEnd-5',
                      'VegSeasonLength-5']
        VegSeason2 = ['VegSeasonDayStart-2',
                      'VegSeasonDayEnd-2',
                      'VegSeasonLength-2']
        if i_ in VegSeason5:
            VegSeason = VegSeason5
        elif i_ in VegSeason2:
            VegSeason = VegSeason2
        iSKIP += VegSeason
        o = _d2(VegSeason[:2], CL_, out=True)
        ik_ = VegSeason[-1]
        t000 = l__('{} {} ... predata'.format(idxD[ik_][0], ik_))
        o = o[1] - o[0]
        _SV(idxD, ik_, o, dgpic_)
        ll_(ik_, t000)

    def _rho_ps_p(cL, p, i_, _nm=None):
        t000 = l__('{} {} ... predata'.format(idxD[i_][0], i_))
        o = cL[0].copy(rho_fr_t_q_p_(cL[0].data, cL[1].data, p))
        if len(cL) > 2:
            o = iris.util.mask_cube(o, cL[2].data < p)
        pst_(o, 'air density', 'kg m-3', 'rho')
        ll_('{} {} ... predata'.format(idxD[i_][0], i_), t000)
        _d0(i_, (o,), _nm=_nm)

    def _d4(i_, CL_):
        t000 = l__('{} {}'.format(idxD[i_][0], i_))
        p = idxD_[i_][5][2][0]
        c_ps = c_kwArgs['c_ps'] if 'c_ps' in c_kwArgs else None
        cL = (*CL_, c_ps) if c_ps else CL_
        _f_n(_rho_ps_p, cL, p, i_)
        _TO1(idxD, i_, dgpic_)                                                 #_TO1(idxD_, i_, dgpic_, freq=None)
        ll_(i_, t000)

    def _rho_ps(cL, i_, _nm=None):
        t000 = l__('{} {} ... predata'.format(idxD[i_][0], i_))
        o = cL[0].copy(rho_fr_t_q_p_(cL[0].data, cL[1].data, cL[2].data))
        pst_(o, 'air density', 'kg m-3', 'rho')
        ll_('{} {} ... predata'.format(idxD[i_][0], i_), t000)
        _d0(i_, (o,), _nm=_nm)

    def _d5(i_, CL_):
        t000 = l__('{} {}'.format(idxD[i_][0], i_))
        _f_n(_rho_ps, cL, i_)
        _TO1(idxD, i_, dgpic_)
        ll_(i_, t000)

    def _d6(i_, CL_):
        cL, o = cL_[:-1], cL_[-1]
        fK_ = idxD_[i_][5][3]
        fK_.update(dict(cSnw=o))
        _d1(i_, cL, fK_=fK_)

    def _d7(i_, CL_):
        o = CL_[0].copy(CL_[0].data / CL_[1].data)
        _d0(i_, (o,))

    def _wind_msk(cL, p, i_, _nm=None):
        if len(cL) > 1:
            iris.util.mask_cube(cL[0], cL[1].data < p)
        _d0(i_, cL[:1], _nm=_nm)

    def _d8(i_, CL_):
        t000 = l__('{} {}'.format(idxD[i_][0], i_))
        p = idxD_[i_][5][2][0]
        c_ps = c_kwArgs['c_ps'] if 'c_ps' in c_kwArgs else None
        o = ws_cube(*CL_)
        cL = (o, c_ps) if c_ps else (o,)
        _f_n(_wind_msk, cL, p, i_)
        _TO1(idxD, i_, dgpic_)
        ll_(i_, t000)

    for i_ in il_:
        if i_ in iSKIP:
            continue
        _iexe(i_)


class _CLIMIDX:

    def __init__(self,
                 oDir,                                                         # dg[p]ic_ in _exe
                 user_cfg=None,                                                # dgpi[c]_ in _exe
                 folder='cordex', regD=None, varD=None,                        # inDF
                 idxD=None, iL=[], gMthd=None, vnX=9                           # idxL
                 ):
                #iDir, tInt,                                                   # inDF
                #DSnm=None, GWL=None,                                          # [dg]pic_ in _exe
                #y0y1=None                                                     # _d2 in _exe & inDF
        self.oDir = oDir
        self.folder = folder
        self.user_cfg = user_cfg
        iK_ = dict(iL=iL, gMthd=gMthd, vnX=vnX)
        if idxD:
            iK_.update(dict(idxD=idxD))
        self._idxL = idxL(**iK_)
        self.idxD = self._idxL.idxD
        iK_ = dict(folder=folder, regD=regD)
        if varD:
            iK_.update(dict(varD=varD))
        self._inDF = inDF(**iK_)

    @abstractmethod
    def dgyit_(self):                                                          # self.DSnm, self.GWL, self.y0y1, self.iDir, self.tInt
        pass

    def chk_dgyit_(self):
        try:
            self.DSnm, self.GWL, self.y0y1, self.iDir, self.tInt
        except:
            emsg = "CALL FAILED!!! Run self.dgyit_() before the call!"
            raise Exception(emsg)

    def __call__(self, ccc=None):

        def _xyz(vd):
            if isinstance(vd, tuple):
                il_ = vd[0]
                ccc = self._inDF(vd, self.iDir, self.tInt, period=self.y0y1)
                _exe(self.idxD,                                                # _exe(idxD,
                     dn=self.DSnm,                                             #      dn=None, gwl=None, po_=None, il_=None, user_cfg=None,
                     gwl=self.GWL,                                             #      y0y1=None,    
                     po_=self.oDir,                                            #      **c_kwArgs)
                     il_=il_,
                     user_cfg=self.user_cfg,
                     y0y1=self.y0y1,
                     **ccc
                     )
            elif isinstance(vd, list):
                for i in vd:
                    _xyz(i)
            else:
                emsg = "check self.VD!"
                raise Exception(emsg)

        self.chk_dgyit_()
        if ccc:
            _exe(self.idxD,
                 dn=self.DSnm,
                 gwl=self.GWL,
                 po_=self.oDir,
                 il_=self.iL,
                 user_cfg=self.user_cfg,
                 y0y1=self.y0y1,
                 **ccc)
        else:
            vD = self._idxL(self.tInt)
            _xyz(vD)
