import numpy as np
import iris
import iris.coord_categorisation as ica
import glob
import os

from uuuu import (
        isIter_, ll_, l__, ouniqL_, pure_fn_, p_deoverlap_, p_least_,          #ffff
        schF_keys_,                                                            #ffff
        concat_cube_, pSTAT_cube, extract_byAxes_, cubesv_                     #cccc
        )


from .diction_ import v__


_opj = os.path.join


def _cmip6_dir_finfo(
        idir,
        var='*', freq='*', gcm='*', exp='*', rip='*', gl='*', ver='*', p=None,
        ext='.nc'
        ):
    """
    Purpose: get data info from a cmip5 data directory
    """
    varL, freqL, gcmL, expL, ripL, glL, pL = [[] for _ in range(7)]
    s = '_'
    if isinstance(idir, str):
        idir = [idir]
    files = []
    for i in idir:
        if p is None:
            tmp = s.join([var, '*' + freq, gcm, exp, rip, gl]) + '*' + ext
        else:
            tmp = s.join([var, '*' + freq, gcm, exp, rip, gl, p]) + ext
        files += glob.glob(_opj(i, var, gl, ver, tmp))
    files.sort()
    for f in files:
        s_s = pure_fn_(f).split(s)
        varL.append(s_s[0])
        gcmL.append(s_s[2])
        expL.append(s_s[3])
        ripL.append(s_s[4])
        glL.append(s_s[5])
        if s_s[1].upper() != 'FX':
            freqL.append(s_s[1])
            pL.append(s_s[6])
        else:
            freqL.append(s_s[1][-3:])
    return {'var': ouniqL_(varL), 'freq': ouniqL_(freqL),
            'gcm': ouniqL_(gcmL), 'exp': ouniqL_(expL),
            'rip': ouniqL_(ripL), 'gl':ouniqL_(glL),
            'p': ouniqL_(pL), 'fn': files}


def _cmip5_dir_finfo(
        idir,
        var='*', freq='*', gcm='*', exp='*', rip='*', p=None,
        ext='.nc'
        ):
    """
    Purpose: get data info from a cmip5 data directory
    """
    varL, freqL, gcmL, expL, ripL, pL = [[] for _ in range(6)]
    s = '_'
    if isinstance(idir, str):
        idir = [idir]
    files = []
    for i in idir:
        if p is None:
            tmp = s.join([var, freq, gcm, exp, rip]) + '*' + ext
        else:
            tmp = s.join([var, freq, gcm, exp, rip, p]) + ext
        files += glob.glob(_opj(i, tmp))
    files.sort()
    for f in files:
        s_s = pure_fn_(f).split(s)
        varL.append(s_s[0])
        freqL.append(s_s[1])
        gcmL.append(s_s[2])
        expL.append(s_s[3])
        ripL.append(s_s[4])
        if s_s[1].upper() != 'FX':
            pL.append(s_s[5])
    return {'var': ouniqL_(varL), 'freq': ouniqL_(freqL),
            'gcm': ouniqL_(gcmL), 'exp': ouniqL_(expL),
            'rip': ouniqL_(ripL), 'p': ouniqL_(pL),
            'fn': files}


def _cordex_dir_finfo(
        idir, var='*', dm='*', gcm='*', exp='*', rip='*',
        rcm='*', ver='*', freq='*', p=None,
        ext='.nc'
        ):
    """
    Purpose: get data info from a cordex data directory
    """
    varL, dmL, gcmL, expL, ripL, rcmL, verL, freqL, pL = [[] for _ in range(9)]
    s = '_'
    if isinstance(idir, str):
        idir = [idir]
    files = []
    for i in idir:
        if p is None:
            tmp = s.join([var, dm, gcm, exp, rip, rcm, ver, freq]) + '*' + ext
        else:
            tmp = s.join([var, dm, gcm, exp, rip, rcm, ver, freq, p]) + ext
        files += glob.glob(_opj(i, tmp))
    files.sort()
    for f in files:
        s_s = pure_fn_(f).split(s)
        varL.append(s_s[0])
        dmL.append(s_s[1])
        gcmL.append(s_s[2])
        expL.append(s_s[3])
        ripL.append(s_s[4])
        rcmL.append(s_s[5])
        verL.append(s_s[6])
        freqL.append(s_s[7])
        if s_s[7].upper() != 'FX':
            pL.append(s_s[8])
    return {'var': ouniqL_(varL), 'dm': ouniqL_(dmL),
            'gcm': ouniqL_(gcmL), 'exp': ouniqL_(expL),
            'rip': ouniqL_(ripL), 'rcm': ouniqL_(rcmL),
            'ver': ouniqL_(verL), 'freq': ouniqL_(freqL),
            'p': ouniqL_(pL), 'fn': files}


def _min_fselect_(dir_finfo, period=None):
    """
    Purpose: select fewest files from finfo of a cmip5/cordex data directory
    """
    pp = p_deoverlap_(p_least_(dir_finfo['p'], *period) if period else
                      dir_finfo['p'])
    fn = [f for f in dir_finfo['fn'] if any(i in f for i in pp)]
    dir_finfo.update({'p': pp, 'fn': fn})



def _cmip6_dir_cubeL(
        idir,
        var='*', freq='*', gcm='*', exp='*', rip='*', gl='gr', ver='*', p='*',
        ext='.nc', period=None, ifconcat=False
        ):
    """
    Purpose: load cube list from a cmip5 data directory
    """
    #warnings.filterwarnings("ignore", category=UserWarning)
    info = _cmip6_dir_finfo(idir, var=var, freq=freq, gcm=gcm, exp=exp,
                           rip=rip, gl=gl, ver=ver, p=p, ext=ext)
    if (len(info['gcm']) * len(info['var']) * len(info['freq'])
        * len(info['rip'])) > 1 and ifconcat:
        raise Exception("no idea how to organize kArgs!")
    _min_fselect_(info, period)
    cubeL = iris.load(info['fn'])
    if var != '*':
        varCstr_ = iris.Constraint(cube_func=lambda c: c.var_name == var)
        cubeL = cubeL.extract(varCstr_)
    if len(cubeL) == 0:
        return None
    elif ifconcat:
        try:
            cube = concat_cube_(cubeL, atol=1e-5)
            p = sorted(info['p'])
            return {'cube': cube,
                    'p': '-'.join((p[0].split('-')[0], p[-1].split('-')[-1]))}
        except:
            ll_('_cmip6_dir_cubeL: concat_cube_ error, return None instead')
            return {'cube': None,
                    'p': '-'.join((p[0].split('-')[0], p[-1].split('-')[-1]))}
    else:
        return {'cube': cubeL, 'p': info['p']}


def _cmip5_dir_cubeL(
        idir,
        var='*', freq='*', gcm='*', exp='*', rip='*', p='*',
        ext='.nc', period=None, ifconcat=False
        ):
    """
    Purpose: load cube list from a cmip5 data directory
    """
    #warnings.filterwarnings("ignore", category=UserWarning)
    info = _cmip5_dir_finfo(idir, var=var, freq=freq, gcm=gcm, exp=exp,
                           rip=rip, p=p, ext=ext)
    if (len(info['gcm']) * len(info['var']) * len(info['freq'])
        * len(info['rip'])) > 1 and ifconcat:
        raise Exception("no idea how to organize kArgs!")
    _min_fselect_(info, period)
    cubeL = iris.load(info['fn'])
    if var != '*':
        varCstr_ = iris.Constraint(cube_func=lambda c: c.var_name == var)
        cubeL = cubeL.extract(varCstr_)
    if len(cubeL) == 0:
        return None
    elif ifconcat:
        try:
            cube = concat_cube_(cubeL, atol=1e-5)
            p = sorted(info['p'])
            return {'cube': cube,
                    'p': '-'.join((p[0].split('-')[0], p[-1].split('-')[-1]))}
        except:
            ll_('_cmip5_dir_cubeL: concat_cube_ error, return None instead')
            return {'cube': None,
                    'p': '-'.join((p[0].split('-')[0], p[-1].split('-')[-1]))}
    else:
        return {'cube': cubeL, 'p': info['p']}


def _norcp_dir_cubeL(
        idir,
        var='*', dm='*', gcm='*', exp='*', rcm='*', freq='*', ver='*', rip='*',
        p='*',
        ext='.nc', period=None, ifconcat=False
        ):
    """
    Purpose: load cube list from a cordex data directory
    """
    if isinstance(idir, str):
        idir_ = _opj(idir, var)
    else:
        idir_ = [_opj(i, var) for i in idir]
    return cordex_dir_cubeL(idir_, var=var, dm=dm, gcm=gcm, exp=exp, rcm=rcm,
                            freq=freq, ver=ver, rip=rip, p=p, ext=ext,
                            period=period, ifconcat=ifconcat)


def _cordex_dir_cubeL(
        idir,
        var='*', dm='*', gcm='*', exp='*', rcm='*', freq='*', ver='*', rip='*',
        p='*',
        ext='.nc', period=None, ifconcat=False):
    """
    Purpose: load cube list from a cordex data directory
    """
    #warnings.filterwarnings("ignore", category=UserWarning)
    info = _cordex_dir_finfo(idir, var=var, dm=dm, gcm=gcm, exp=exp,
                            rip=rip, rcm=rcm, p=p, ext=ext)
    if ((len(info['dm']) * len(info['gcm']) * len(info['var'])
         * len(info['rcm']) * len(info['freq']) * len(info['rip'])) > 1
        and ifconcat):
        raise Exception("no idea how to organize kArgs!")
    _min_fselect_(info, period)
    cubeL = iris.load(info['fn'])
    if var != '*':
        varCstr_ = iris.Constraint(cube_func=lambda c: c.var_name == var)
        cubeL = cubeL.extract(varCstr_)
    if len(cubeL) == 0:
        return None
    elif ifconcat:
        try:
            cube = concat_cube_(cubeL)
            p = sorted(info['p'])
            return {'cube': cube,
                    'p': '-'.join((p[0].split('-')[0], p[-1].split('-')[-1]))}
        except:
            ll_('_cordex_dir_cubeL : concat_cube_ error, return None instead')
            return {'cube': None,
                    'p': '-'.join((p[0].split('-')[0], p[-1].split('-')[-1]))}
    else:
        return {'cube': cubeL, 'p': info['p']}


def _xx(iDir, tInt):
    if isinstance(iDir, str):
        return _opj(iDir, '*' + tInt)
    elif isinstance(iDir, (tuple, list, set, np.ndarray)):
        return [_xx(i, tInt) for i in iDir]
    else:
        raise ValueError("'iDir' must be str type or array-like of str")


def _to_xhr(cube, x=6, valid=True):
    nh = 24 / x
    ica.add_categorised_coord(cube, 'xxx', 'time',
                              lambda coord, v: np.ceil(v / nh) * nh)
    o = cube.aggregated_by('xxx',iris.analysis.MEAN)
    if valid:
        dpo = np.diff(o.coord('xxx').points)
        ddpo = np.diff(np.where(dpo != 0)[0])[1]
        sss = None if dpo[ddpo-1] != 0 else 1
        eee = None if dpo[-ddpo] != 0 else -1
    else:
        sss, eee = None, None
    cube.remove_coord('xxx')
    return extract_byAxes_(o, 'time', np.s_[sss:eee])


def _xxx(cube, tInt0, tInt1):
    tInts = ['mon', 'day', '6hr', '3hr', '1hr']
    if any(i not in tInts for i in [tInt0, tInt1]):
        raise ValueError("unknown tIntuency names!")
    if tInts.index(tInt1) > tInts.index(tInt0):
        raise ValueError("cannot convert from low to high tIntuency!")
    if tInt1 == 'mon':
        return pSTAT_cube(cube, 'MEAN', 'month')
    elif tInt1 == 'day':
        return pSTAT_cube(cube, 'MEAN', 'day')
    elif tInt1 == '6hr':
        if cube.cell_methods and cube.cell_methods[0].method.upper() == 'MEAN':
            return _to_xhr(cube)
        else:
            nn = 2 if tInt0 == '3hr' else 6
            return extract_byAxes_(cube, 'time', np.s_[::nn])
    else:
        if cube.cell_methods and cube.cell_methods[0].method.upper() == 'MEAN':# average
            return _to_xhr(cube, x=3)
        else:                                                                  # instantaneous
            return extract_byAxes_(cube, 'time', np.s_[::3])


def _oe(iDir, tInt, folder='cordex', reg_d=None, **kwargs):
    tInts = ['mon', 'day', '6hr', '3hr', '1hr']
    _rf = eval('_{}_dir_cubeL'.format(folder))
    o = None
    e = None
    for f_ in tInts[tInts.index(tInt):]:
        p_ = _xx(iDir, f_)
        o = _rf(p_, ifconcat=True, **kwargs)
        if o:
            o = o['cube']
            if o:
                if reg_d:
                    o = intersection_(o, **reg_d)
                if 'period' in kwargs and kwargs['period'] is not None:
                    o = extract_period_cube(o, *kwargs['period'], yy=True)
                    e = None if o else 'yye'
                break
            else:
                e = 'cce'
    if f_ != freq and isinstance(o, iris.cube.Cube):
        o = _xxx(o, f_, freq)
    return (o, e)


def _dd(var, ccc, varD_):
    if 'c_{}'.format(var) not in ccc and var in varD_:
        func_ = varD_[var][0]
        if isinstance(func_, str):
            o = eval(func_)
        else:
            c_ = (ccc['_{}'.format(i)] for i in varD_[var][1])
            fA_ = varD_[var][2][0]
            fK_ = varD_[var][2][1]
            o = func_(*c_, *fA_, **fK_)
            func__ = varD_[var][3]
            if func__ is not None and 'c_ps' in ccc:
                o = func__(o, ccc['c_ps'].data < varD_[var][4])
        return o


class inDF:
    def __init__(self,
                 folder='cordex', regD=None,
                 varD=v__):
        self.iDir = iDir
        self.tInt = tInt
        self.varD = varD
        self.kA0 = dict(folder=folder, reg_d=regD)

    def __call__(self, vD, iDir, tInt, period=None):
        """
        vD fmt: ([idx0, ...], {c_var0: 'var0', ...})
        """
        ccc = dict()                                                           # derive first-hand variable cubes from data files
        t000 = l__("_oe(): {} variables".format(len(vD[1].keys())))
        ll_(', '.join(vD[1].keys()))
        ll_(', '.join(vD[0]))
        tmp, tmp_, tmp__ = [], [], []
        wsyes = False                                                          # SURFACE WIND SPEED
        for kk in vD[1].keys():
            if wsyes and kk in ('c_uas', 'c_vas'):                             # SURFACE WIND SPEED
                continue                                                       # SURFACE WIND SPEED
            cc, ee = _oe(iDir, tInt, var=vd[1][kk], period=period, **self.kA0)
            if cc:
                tmp.append(kk)
                ccc.update({kk: cc})
                if kk == 'c_sfcWind':                                          # SURFACE WIND SPEED
                    wsyes = True                                               # SURFACE WIND SPEED
            if ee == 'cce':                                                    # error in concate_cube_ call
                tmp_.append(kk)
            elif ee == 'yye':                                                  # error in data temporal coverage
                tmp__.append(kk)
        if len(tmp__) > 0:
            ll_('YYE: {}'.format(', '.join(tmp__)))
        if len(tmp_) > 0:
            ll_('CCE: {}'.format(', '.join(tmp_)))
        if len(tmp) > 0:
            ll_(', '.join(tmp))
        ll_("_oe()", t000)

        ddd = dict()                                                           # derive secondary variable cubes
        t000 = l__("_dd(): {} variables".format(len(vD[2].keys())))
        ll_(', '.join(vD[1].keys()))
        ll_(', '.join(vD[0]))
        for kk in vD[2].keys():
            cc = _dd(kk, ccc, self.varD)
            if cc:
                ddd.update({kk: cc})
        ccc.update(ddd)
        return ccc


def _meta(idxD_, i_, cube):
    pK_ = dict(var_name=re.sub('\W', '_', i_).lower())
    if idxD_[i_][6]:
        pK_.update(idxD_[i_][6])
    pst_(cube, **pK_)


def _dn(*terms):
    return '_'.join(i for i in terms if i)


def _fnfmt(idir, *terms):
    nm_ = '{}.nc'.format('_'.join(i for i in terms if i))
    return _opj(idir, nm_)


def _get_freq(idxD_, i_, user_cfg=None):
    ii = i_[0] if isIter_(i_) else i_
    if user_cfg:
        try:
            o = user_cfg['freq_cfg'][ii]
        except:
            o = idxD_[ii][4]
        return o
    else:
        return idxD_[ii][4]


def _SV(idxD_, i_, o, dgpic_, freq=None, _nm=None):
    dn, gwl, po_ = dgpic_[:3]

    def _fn(fff):
        return _fnfmt(po_, i_, dn, fff, gwl, _nm)

    freq = freq if freq else _get_freq(idxD_, i_, dgpic_[4])
    assert isIter_(freq), "type {!r} not acceptable here!".format(type(freq))
    if len(freq) == 1:
        _meta(idxD_, i_, o)
        cubesv_(o, _fn(freq[0]))
    else:
        for oo, ff in zip(o, freq):
            _meta(idxD_, i_, oo)
            cubesv_(oo, _fn(ff))


def _TO1(idxD_, i_, dgpic_, freq=None):
    dn, gwl, po_ = dgpic_[:3]
    def _fns(fff):
        return schF_keys_(
                po_,
                '_'.join((i for i in (i_, dn, fff, gwl, '__*') if i)))

    def _fn(fff):
        return _fnfmt(po_, i_, dn, fff, gwl)

    freq = freq if freq else _get_freq(idxD_, i_, dgpic_[4])
    if len(freq) == 1:
        fns = _fns(freq[0])
        if fns:
            o = iris.load(fns)
            cubesv_(concat_cube_(o), _fn(freq[0]))
            for i in fns:
                os.remove(i)
    else:
        for ff in freq:
            fns = _fns(ff)
            if fns:
                o = iris.load(fns)
                cubesv_(concat_cube_(o), _fn(ff))
                for i in fns:
                    os.remove(i)
