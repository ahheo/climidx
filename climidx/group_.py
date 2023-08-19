import os

from uuuu import ouniqL_, ss_fr_sl_, nli_, flt_l

from .diction_ import i__


def _vvv(a, vnX=9):                                                            # for memory limit, up to vnX (default=9) cubes
    a = sorted(a, key=lambda x:(len(x[1]), sorted(x[1][0])), reverse=True)     # sort indices according to the number of variables used
    while len(a) > 0:
        aa, bb = [a[0][0]], a[0][1]
        del(a[0])
        if len(a) > 0:
            for i in a.copy():
                tmp = ouniqL_(bb + i[1])
                if len(tmp) <= vnX:
                    aa.append(i[0])
                    bb = tmp
                    a.remove(i)
        yield (aa, bb)


def _vv(idxD_, iL_, tInt, gMthd=None, vnX=9):
    if tInt not in ('mon',) and gMthd is None:
        tmp = [(i, idxD_[i][2]) for i in iL_ if idxD_[i][1] == tInt]
        return list(_vvv(tmp, vnX))
    else:
        tmp = [idxD_[i][2] for i in iL_ if idxD_[i][1] == tInt]
        tmp_ = set(flt_l(tmp))
        if tInt not in ('mon',) and len(tmp_) > vnX:
            if gMthd == 'v':                                                   # according to shared variables
                tmp__ = ss_fr_sl_(list(map(set, tmp)))
                return [(
                    [i_ for i_ in iL_ if all(ii in i for ii in idxD_[i_][2])],
                    i
                    ) for i in tmp__]
            elif gMthd == 'i':                                                 # half & half
                iL__ = [iL_[:(len(iL_)//2)], iL_[(len(iL_)//2):]]
                tmp__ = [_vv(idxD_, i, tInt, gMthd=gMthd, vnX=vnX)
                         for i in iL__]
                return nli_(tmp__)
        else:
            return ([i for i in iL_ if idxD_[i][1] == tInt], tmp_)


def _vd2(iL_, idxD_):
    return ouniqL_(
        flt_l(idxD_[i][5][1] for i in iL_ if idxD_[i][5][1] is not None)
        )
    #bb = ouniqL_(flt_l(idxD_[i][2] for i in iL_))
    #return [i for i in aa if i not in bb]


def _vd_fr_vv(vl, idxD_):
    if isinstance(vl, tuple):
        vd = dict()
        for i in vl[1]:
            vd.update({'c_' + i: i})
        vd2 = dict()
        for i in _vd2(vl[0], idxD_):
            vd2.update({'c_' + i: i})
        return (vl[0], vd, vd2)
    elif isinstance(vl, list):
        return [_vd_fr_vv(i, idxD_) for i in vl]
    else:
        raise Exception('check input!')


class idxL:

    def __init__(self, idxD=i__, iL=[], gMthd=None, vnX=9):
        if not isinstance(idxD, dict):
            emsg = "type should be dict!"
            raise TypeError(emsg)
        self.idxD = idxD
        if iL:
            self.iL = iL
            if any(i not in self.idxD.keys() for i in self.iL):                # check selected indices all defined in self.idxD
                emsg = "specified indices contains undefined ones!"
                raise ValueError(emsg)
        else:
            self.iL = list(self.idxD.keys())
        self.gMthd = gMthd
        self.vnX = vnX

    def __call__(self, tInt):
        return _vd_fr_vv(
                _vv(
                    self.idxD, self.iL, tInt, self.gMthd, self.vnX
                    ),
                self.idxD)
