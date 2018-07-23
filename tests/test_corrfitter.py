# Copyright (c) 2017-18 G. Peter Lepage.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version (see <http://www.gnu.org/licenses/>).
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

from __future__ import print_function   # makes this work for python2 and 3

import inspect
import os
import unittest

import numpy as np
import gvar as gv
from corrfitter import *

try:
    import h5py
    NO_H5PY = False
except:
    NO_H5PY = True

class test_corr2(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_init(self):
        " Corr2 "
        m = Corr2(
            'tag', a='a', b='b', dE='dE', tmin=2, tp=6, otherdata='xtag',
            reverseddata='revtag'
            )
        self.assertTrue(
            m.datatag == 'tag' and m.a == ('a',None) and m.b == ('b',None)
            and m.dE == ('dE',None) and np.all(m.tdata == [0,1,2,3,4,5])
            and np.all(m.tfit == [2,3]) and m.tp == 6
            )
        self.assertEqual(m.otherdata, ['xtag'])
        self.assertEqual(m.reverseddata, ['revtag'])

        # osc
        m = Corr2(
            'tag', a=('a','ao'), b=('b','bo'), dE=('dE','dEo'),
            tmin=2, tp=6, otherdata=['xtag'], reverseddata=['revtag']
            )
        self.assertTrue(
            m.datatag == 'tag' and m.a == ('a','ao') and m.b == ('b','bo')
            and m.dE == ('dE','dEo') and np.all(m.tdata == [0,1,2,3,4,5])
            and np.all(m.tfit == [2,3]) and m.tp == 6
            )
        self.assertEqual(m.otherdata, ['xtag'])
        self.assertEqual(m.reverseddata, ['revtag'])

        # consistency check
        with self.assertRaises(ValueError):
            m = Corr2(
            'tag', a=('a','ao'), b='b', dE=('dE','dEo'),
            tmin=2, tp=6, otherdata=['xtag'], reverseddata=['revtag']
            )

    def test_init_tfit(self):
        " Corr2(...tfit=...) "
        m = Corr2(
            'tag', a='a', b='b', dE='dE', tmin=2, tdata=[1,2,3,4,5]
            )
        self.assertTrue(
            np.all(m.tdata == [1,2,3,4,5]) and
            np.all(m.tfit == [2,3,4,5])
            )

        m = Corr2(
            'tag', a='a', b='b', dE='dE', tmin=2, tmax=5
            )
        self.assertTrue(
            np.all(m.tdata == [0,1,2,3,4,5]) and
            np.all(m.tfit == [2,3,4,5])
            )

        m = Corr2(
            'tag', a='a', b='b', dE='dE', tmin=1, tmax=2, tp=6
            )
        self.assertTrue(
            np.all(m.tdata == [0,1,2,3,4,5]) and
            np.all(m.tfit == [1,2])
            )

        m = Corr2(
            'tag', a='a', b='b', dE='dE', tfit=[2,3,4], tp=6
            )
        self.assertTrue(
            np.all(m.tdata == [0,1,2,3,4,5]) and
            np.all(m.tfit == [2,3])
            )

        m = Corr2(
            'tag', a='a', b='b', dE='dE', tfit=[2,3,4], tp=7
            )
        self.assertTrue(
            np.all(m.tdata == [0,1,2,3,4,5,6]) and
            np.all(m.tfit == [2,3])
            )

        m = Corr2(
            'tag', a='a', b='b', dE='dE', tfit=[1,3], tdata=[0,1,3,5]
            )
        self.assertTrue(
            np.all(m.tdata == [0,1,3,5]) and
            np.all(m.tfit == [1,3])
            )

        m = Corr2(
            'tag', a='a', b='b', dE='dE', tfit=[1,2,3,4], tdata=[0,1,2,3,4,5],
            tp=-6
            )
        self.assertTrue(
            m.tp == -6 and
            np.all(m.tdata == [0,1,2,3,4,5]) and
            np.all(m.tfit == [1,2,3])
            )

    def test_buildprior(self):
        " Corr2.buildprior "
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=6)
        prior = gv.BufferDict()
        prior['a'] = ['1(1)', '2(1)', '3(1)', '4(1)']
        prior['b'] = ['5(1)', '6(1)', '7(1)', '8(1)']
        prior['dE'] = ['9(1)', '10(1)', '11(1)', '12(1)']
        prior = gv.gvar(prior)
        # basic behavior
        mprior = m.buildprior(prior)
        for k in prior:
            self.assertEqual(str(prior[k]), str(mprior[k]))
        # marginalization
        mprior = m.buildprior(prior, nterm=3)
        for k in prior:
            self.assertEqual(str(prior[k][:3]), str(mprior[k]))
        # extend=True
        prior = gv.BufferDict()
        prior['log(a)'] = ['1(1)', '2(1)', '3(1)', '4(1)']
        prior['log(b)'] = ['5(1)', '6(1)', '7(1)', '8(1)']
        prior['log(dE)'] = ['9(1)', '10(1)', '11(1)', '12(1)']
        prior = gv.gvar(prior)
        mprior = m.buildprior(prior)
        for k in prior:
            self.assertEqual(str(prior[k]), str(mprior[k]))
        # oscillating pieces
        m = Corr2(
            'tag', a=('a','ao'), b=('b','bo'), dE=('dE','dEo'), tmin=1, tp=6
            )
        prior = gv.BufferDict()
        prior['a'] = ['1(1)', '2(1)', '3(1)', '4(1)']
        prior['ao'] = ['1(1)', '2(1)', '3(1)', '4(1)']
        prior['b'] = ['5(1)', '6(1)', '7(1)', '8(1)']
        prior['bo'] = ['5(1)', '6(1)', '7(1)', '8(1)']
        prior['dE'] = ['9(1)', '10(1)', '11(1)', '12(1)']
        prior['dEo'] = ['9(1)', '10(1)', '11(1)', '12(1)']
        prior = gv.gvar(prior)
        for k in ['ao', 'bo', 'dEo']:
            prior[k] *= 10
        mprior = m.buildprior(prior)
        for k in prior:
            self.assertEqual(str(prior[k]), str(mprior[k]))

    def test_builddata(self):
        " Corr2.builddata "
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=6)
        data = gv.BufferDict()
        data['tag'] = ['1.0(1)', '2.0(1)', '3.0(1)', '4.0(1)', '5.0(1)', '6.0(1)']
        data['otag'] = ['2.0(1)', '3.0(1)', '4.0(1)', '5.0(1)', '6.0(1)', '7.0(1)']
        data['xtag'] = ['2.0(1)', '3.0(1)', '4.0(1)', '5.0(1)', '6.0(1)', '7.0(1)']
        data = gv.gvar(data)
        mdata = m.builddata(data)
        self.assertEqual(str(mdata),'[4.000(71) 4.000(71) 4.00(10)]')

        # reverse=True
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=6, reverse=True)
        mdata = m.builddata(data)
        self.assertEqual(str(mdata),'[4.000(71) 4.000(71) 4.00(10)]')

        # anti-periodic
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=-6)
        rdata = gv.gvar(data)
        rdata['tag'][-2:] *= -1
        mdata = m.builddata(rdata)
        self.assertEqual(str(mdata),'[4.000(71) 4.000(71) 4.00(10)]')

        # otherdata
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=6,  otherdata='otag')
        mdata = m.builddata(data)
        self.assertEqual(str(mdata),'[4.500(50) 4.500(50) 4.500(71)]')

        # reverseddata
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=6,  reverseddata='otag')
        mdata = m.builddata(data)
        self.assertEqual(str(mdata),'[4.500(50) 4.500(50) 4.500(71)]')

        # non-periodic
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tmax=4)
        mdata = m.builddata(data)
        self.assertEqual(str(mdata), str(data['tag'][1:5]))

        # non-periodic, reverse=True
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tmax=4, reverse=True)
        mdata = m.builddata(data)
        self.assertEqual(str(mdata[1:]), str(data['tag'][4:1:-1]))
        self.assertEqual(str(mdata[0]), str(data['tag'][5]))

    def test_builddataset(self):
        " Corr2.builddata "
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=6)
        dataset = gv.BufferDict()
        dataset['tag'] = np.array(
            [[1., 2., 3., 4., 5., 6.], [2., 3., 4., 5., 6., 7.]]
            )
        dataset['otag'] = np.array(
            [[1., 2., 3., 4., 5., 6.], [2., 3., 4., 5., 6., 7.]]
            )[:,::-1]
        dataset['xtag'] = np.array(
            [[1., 2., 3., 4., 5., 6.], [2., 3., 4., 5., 6., 7.]]
            ) * 10.
        mdataset = m.builddataset(dataset)
        self.assertEqual(
            str(mdataset),str(np.array([[4., 4, 4], [5, 5., 5.]]))
            )

        # reverse=True
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=6, reverse=True)
        mdataset = m.builddataset(dataset)
        self.assertEqual(
            str(mdataset),str(np.array([[4., 4, 4], [5, 5., 5.]]))
            )

        # anti-periodic
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=-6)
        rdataset = gv.BufferDict()
        rdataset['tag'] = np.array(
            [[1., 2., 3., 4., 5., 6.], [2., 3., 4., 5., 6., 7.]]
            )
        rdataset['xtag'] = np.array(
            [[1., 2., 3., 4., 5., 6.], [2., 3., 4., 5., 6., 7.]]
            ) * 10.
        rdataset['tag'][:, -2:] *= -1
        mdataset = m.builddataset(rdataset)
        self.assertEqual(
            str(mdataset),str(np.array([[4., 4, 4], [5, 5., 5.]]))
            )

        # otherdata
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=6,  otherdata='otag')
        mdataset = m.builddataset(dataset)
        self.assertEqual(
            str(mdataset),str(np.array([[3.5, 3.5, 3.5], [4.5, 4.5, 4.5]]))
            )

        # reverseddata
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=6,  reverseddata='otag')
        mdataset = m.builddataset(dataset)
        self.assertEqual(
            str(mdataset),str(np.array([[3.5, 3.5, 3.5], [4.5, 4.5, 4.5]]))
            )

        # non-periodic
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tmax=4)
        mdataset = m.builddataset(dataset)
        self.assertEqual(
            str(mdataset),str(dataset['tag'][:,1:5])
            )

        # non-periodic, reverse=True
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tmax=4, reverse=True)
        mdataset = m.builddataset(dataset)
        self.assertEqual(
            str(mdataset[:, 1:]),str(dataset['tag'][:,4:1:-1])
            )
        self.assertEqual(
            str(mdataset[:, 0]),str(dataset['tag'][:,5])
            )

    def test_fitfcn(self):
        " Corr2.fitfcn "
        def fcn(p, t, tp=1000, pfac=1., osc=False):
            t = np.array(t, dtype=int)
            E = np.cumsum(p['dE'])
            ans = np.sum(
                p['a'][:, None] * p['b'][:, None] * (
                    np.exp(- E[:, None] * t[None, :]) +
                    pfac * np.exp(- E[:, None] * (tp - t[None, :]))
                    ),
                axis=0
                )
            if osc:
                return ans * (-1) ** t
            else:
                return ans
        p = gv.BufferDict()
        p['a'] = [1., 2.]
        p['b'] = [3., 4.]
        p['dE'] = [.1, .2]
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=4)
        np.testing.assert_almost_equal(
            fcn(p, t=[1,2.], tp=4, pfac=+1 ), m.fitfcn(p)
            )
        np.testing.assert_almost_equal(
            fcn(p, t=m.tdata, tp=4, pfac=+1 ), m.fitfcn(p, t=m.tdata)
            )

        # anti-periodic
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tp=-4)
        np.testing.assert_almost_equal(
            fcn(p, t=[1,2.], tp=4, pfac=-1 ), m.fitfcn(p)
            )
        np.testing.assert_almost_equal(
            fcn(p, t=m.tdata, tp=4, pfac=-1 ), m.fitfcn(p, t=m.tdata)
            )

        # non-periodic
        m = Corr2('tag', a='a', b='b', dE='dE', tmin=1, tmax=4)
        np.testing.assert_almost_equal(
            fcn(p, t=[1.,2.,3.,4.], pfac=0 ), m.fitfcn(p)
            )

        # oscillating
        po = gv.BufferDict()
        po['a'] = [1.5, 2.5]
        po['b'] = [3.5, 4.5]
        po['dE'] = [.15, .25]
        for k in po:
            p[k + 'o'] = po[k]
        m = Corr2(
            'tag', a=('a','ao'), b=('b','bo'), dE=('dE','dEo'),
            tmin=1, tp=4, s=(1., -1.),
            )
        np.testing.assert_almost_equal(
            fcn(p, t=[1,2], tp=4) - fcn(po, t=[1,2], tp=4, osc=True),
             m.fitfcn(p)
            )

class test_corr3(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_init(self):
        " Corr3 "
        m = Corr3(
            'tag', a='a', dEa='dEa', b='b', dEb='dEb', Vnn='Vnn',
            reverse=True, symmetric_V=True, tmin=2, T=5
            )
        self.assertTrue(
            m.datatag == 'tag' and m.a == ('a',None) and m.b == ('b',None)
            and m.dEa == ('dEa',None) and m.dEb == ('dEb',None)
            and m.V == [['Vnn', None], [None, None]]
            and m.reverse == True and m.symmetric_V == True
            and np.all(m.tdata == [0,1,2,3,4,5]) and np.all(m.tfit == [2,3])
            and m.T == 5 and m.otherdata == [] and
            m.reverseddata == []
            )

        # tfit
        m = Corr3(
            'tag', a='a', dEa='dEa', b='b', dEb='dEb', Vnn='Vnn',
            tfit=[1,3,4], T=5, reverse=True,
            )
        self.assertTrue(
            m.reverse == True and m.symmetric_V == False
            and m.V == [['Vnn', None], [None, None]]
            )
        np.testing.assert_almost_equal([1,3,4], m.tfit)

        # tdata missing t's
        m = Corr3(
            'tag', a='a', dEa='dEa', b='b', dEb='dEb',
            Vnn='Vnn', tmin=2, tdata=[1,2,3,5], T=5,
            )
        self.assertTrue(
            m.reverse == False and m.symmetric_V == False
            and m.V == [['Vnn', None], [None, None]]
            )
        np.testing.assert_almost_equal([2,3], m.tfit)

        # otherdata
        m = Corr3(
            'tag', a='a', dEa='dEa', b='b', dEb='dEb', Vnn='Vnn',
            reverse=True, symmetric_V=True, tmin=2, T=5,
            otherdata='otag',
            )
        self.assertEqual(m.otherdata, ['otag'])
        self.assertEqual(m.reverseddata, [])

        # reverseddata
        m = Corr3(
            'tag', a='a', dEa='dEa', b='b', dEb='dEb', Vnn='Vnn',
            reverse=True, symmetric_V=True, tmin=2, T=5,
            reverseddata='rtag',
            )
        self.assertEqual(m.otherdata, [])
        self.assertEqual(m.reverseddata, ['rtag'])

    def test_init_osc(self):
        " Corr3 "
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo',
            reverse=True, symmetric_V=False, tmin=2, T=5
            )
        self.assertTrue(
            m.datatag == 'tag' and m.a == ('a','ao') and m.b == ('b','bo')
            and m.dEa == ('dEa','dEao') and m.dEb == ('dEb','dEbo')
            and m.V == [['Vnn', 'Vno'], ['Von', 'Voo']]
            and m.reverse == True and m.symmetric_V == False
            and np.all(m.tdata == [0,1,2,3,4,5]) and np.all(m.tfit == [2,3])
            and m.T == 5 and m.otherdata == [] and
            m.reverseddata == []
            )

    def test_init_consistency(self):
        " Corr3 consistency checks "
        args = dict(
            datatag='tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo',
            reverse=False, symmetric_V=False, tmin=2, T=5
            )
        m = Corr3(**args)

        # propagators
        with self.assertRaises(ValueError):
            nargs = dict(args)
            nargs['a'] = 'a'
            m = Corr3(**nargs)
        with self.assertRaises(ValueError):
            nargs = dict(args)
            nargs['dEb'] = 'dEb'
            m = Corr3(**nargs)
        with self.assertRaises(ValueError):
            nargs = dict(args)
            nargs['Voo'] = None
            m = Corr3(**nargs)

        # V[i][j]
        nargs = dict(args)
        nargs['a'] = 'a'
        nargs['dEa'] = 'dEa'
        nargs['Von'] = None
        nargs['Voo'] = None
        m = Corr3(**nargs)
        with self.assertRaises(ValueError):
            nargs['Vno'] = None
            m = Corr3(**nargs)
        with self.assertRaises(ValueError):
            nargs['Von'] = 'Von'
            m = Corr3(**nargs)
        nargs['Von'] = None
        with self.assertRaises(ValueError):
            nargs['Voo'] = 'Voo'
            m = Corr3(**nargs)

        # V[i][j]
        nargs = dict(args)
        nargs['b'] = 'b'
        nargs['dEb'] = 'dEb'
        nargs['Vno'] = None
        nargs['Voo'] = None
        m = Corr3(**nargs)
        with self.assertRaises(ValueError):
            nargs['Von'] = None
            m = Corr3(**nargs)
        with self.assertRaises(ValueError):
            nargs['Vno'] = 'Vno'
            m = Corr3(**nargs)
        nargs['Vno'] = None
        with self.assertRaises(ValueError):
            nargs['Voo'] = 'Voo'
            m = Corr3(**nargs)

        # V[i][j] symmetric_V=True
        nargs = dict(args)
        nargs['symmetric_V'] = True
        del nargs['Von']
        m = Corr3(**nargs)
        with self.assertRaises(ValueError):
            nargs['Vno'] = None
            m = Corr3(**nargs)
        nargs['Vno'] = 'Vno'
        with self.assertRaises(ValueError):
            nargs['a'] = 'a'
            nargs['dEa'] = 'dEa'
            m = Corr3(**nargs)
        nargs['b'] = 'b'
        nargs['dEb'] = 'dEb'
        nargs['Vno'] = None
        nargs['Voo'] = None
        m = Corr3(**nargs)
        with self.assertRaises(ValueError):
            nargs['Vno'] = 'Vno'
            m = Corr3(**nargs)

    def test_builddata(self):
        data = dict(
            tag = ['1(1)', '2(1)', '3(1)', '4(1)'],
            rtag = ['1(1)', '2(1)', '3(1)', '4(1)'],
            otag = ['4(1)', '3(1)', '2(1)', '1(1)'],
            )
        data = gv.gvar(data)

        # tmin
        m = Corr3(
            'tag', tmin=1, T=3,
            a='a', b='b', dEa='dEa', dEb='dEb', Vnn='Vnn',
            )
        self.assertEqual(str(m.builddata(data)), str(data['tag'][1:-1]))

        # reverse=True
        m = Corr3(
            'tag', tmin=1, T=3,
            a='a', b='b', dEa='dEa', dEb='dEb', Vnn='Vnn', reverse=True
            )
        self.assertEqual(str(m.builddata(data)), str(data['tag'][::-1][1:-1]))

        # reverseddata
        m = Corr3(
            'tag', tmin=1, T=3,
            a='a', b='b', dEa='dEa', dEb='dEb', Vnn='Vnn', reverseddata='rtag',
            )
        self.assertEqual(str(m.builddata(data)), '[2.50(71) 2.50(71)]')

        # otherdata
        m = Corr3(
            'tag', tmin=1, T=3,
            a='a', b='b', dEa='dEa', dEb='dEb', Vnn='Vnn', otherdata='otag',
            )
        self.assertEqual(str(m.builddata(data)), '[2.50(71) 2.50(71)]')

    def test_builddataset(self):
        dataset = dict(
            tag = [[1., 2., 3., 4.], [3., 4., 5., 6.]],
            otag = [[1., 2., 3., 4.], [3., 4., 5., 6.]],
            rtag = [[4., 3., 2., 1.], [6., 5., 4., 3.]]
            )

        # tmin
        m = Corr3(
            'tag', tmin=1, T=3,
            a='a', b='b', dEa='dEa', dEb='dEb', Vnn='Vnn',
            )
        self.assertEqual(
            str(m.builddataset(dataset)).replace(' ',''),
            '[[ 2.  3.]\n [ 4.  5.]]'.replace(' ',''),
            )

        # reverse=True
        m = Corr3(
            'tag', tmin=1, T=3,
            a='a', b='b', dEa='dEa', dEb='dEb', Vnn='Vnn', reverse=True
            )
        self.assertEqual(
            str(m.builddataset(dataset)).replace(' ',''),
            '[[ 3.  2.]\n [ 5.  4.]]'.replace(' ','')
            )

        # reverseddata
        m = Corr3(
            'tag', tmin=1, T=3,
            a='a', b='b', dEa='dEa', dEb='dEb', Vnn='Vnn', reverseddata='rtag',
            )
        self.assertEqual(
            str(m.builddataset(dataset)).replace(' ',''),
            '[[ 2.  3.]\n [ 4.  5.]]'.replace(' ','')
            )

        # otherdata
        m = Corr3(
            'tag', tmin=1, T=3,
            a='a', b='b', dEa='dEa', dEb='dEb', Vnn='Vnn', otherdata='otag',
            )
        self.assertEqual(
            str(m.builddataset(dataset)).replace(' ',''),
            '[[ 2.  3.]\n [ 4.  5.]]'.replace(' ','')
            )

    def test_buildprior(self):
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo',
            tmin=2, T=5,
            )
        prior = dict(
            a=['1(1)', '2(1)', '3(1)'],
            b=['3(1)', '4(1)', '5(1)'],
            dEa=['0.1(1)', '0.2(1)', '0.3(1)'],
            dEb=['0.3(1)', '0.4(1)', '0.5(1)'],
            Vnn=3 * [['10(1)', '20(1)', '30(1)']],
            Voo=2* [['40(1)', '50(1)']],
            Von=2 * [['1(1)', '2(1)', '3(1)']],
            Vno=3 * [['3(1)', '4(1)']],
            dummy=['1(10)']
            )
        prior = gv.gvar(prior)
        for k in ['a', 'b', 'dEa', 'dEb']:
            prior[k + 'o'] = 10 * prior[k][:2]
        mprior = m.buildprior(prior)
        for k in prior:
            if k in ['dummy']:
                continue
            self.assertEqual(str(prior[k]), str(mprior[k]))
        self.assertTrue('dummy' not in mprior)

        # marginalization
        mprior = m.buildprior(prior, nterm=(2,1))
        for k in prior:
            if k[0] == 'V' or k in ['dummy']:
                continue
            if k[-1] == 'o':
                self.assertEqual(str(prior[k][:1]), str(mprior[k]))
            else:
                self.assertEqual(str(prior[k][:2]), str(mprior[k]))
            nterm = [2,1]
            for i in range(2):
                for j in range(2):
                    Vij = m.V[i][j]
                    self.assertEqual(
                        str(mprior[Vij]),
                        str(prior[Vij][:nterm[i], :nterm[j]])
                        )
        self.assertTrue('dummy' not in mprior)

    def test_buildprior_sym(self):
        " Corr3.buildprior(...symmetric_V=True) "
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Voo='Voo',
            reverse=False, symmetric_V=True, tmin=2, T=5,
            )
        prior = dict(
            a=['1(1)', '2(1)', '3(1)'],
            b=['3(1)', '4(1)', '5(1)'],
            dEa=['0.1(1)', '0.2(1)', '0.3(1)'],
            dEb=['0.3(1)', '0.4(1)', '0.5(1)'],
            Vnn=['10(1)', '20(1)', '30(1)', '100(1)', '200(1)', '300(1)'],
            Voo=['40(1)', '50(1)', '60(1)', '40(1)', '50(1)', '60(1)'],
            Vno=3 * [['3(1)', '4(1)', '5(1)']],
            Von=3 * [['30(1)', '40(1)', '50(1)']],
            dummy=['1(10)']
            )
        prior = gv.gvar(prior)
        for k in ['a', 'b', 'dEa', 'dEb']:
            prior[k + 'o'] = 10 * prior[k]
        mprior = m.buildprior(prior)
        for k in prior:
            if k in ['dummy', 'Von']:
                continue
            self.assertEqual(str(prior[k]), str(mprior[k]))
        self.assertTrue('dummy' not in mprior)
        self.assertTrue('Von' not in mprior)

        # use Von instead of Vno with reverse=True
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Voo='Voo',
            reverse=True, symmetric_V=True, tmin=2, T=5,
            )
        mprior = m.buildprior(prior)
        for k in prior:
            if k in ['dummy', 'Von']:
                continue
            self.assertEqual(str(prior[k]), str(mprior[k]))
        self.assertTrue('dummy' not in mprior)
        self.assertTrue('Von' not in mprior)
        # check marginalization with reverse=True
        mprior = m.buildprior(prior, nterm=(2,1))

        # symmetric plus marginalization
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Voo='Voo',
            reverse=False, symmetric_V=True, tmin=2, T=5,
            )
        mprior = m.buildprior(prior, nterm=(2,1))
        for k in prior:
            if k[0] == 'V' or k in ['dummy', 'Von']:
                continue
            if k[-1] != 'o':
                self.assertEqual(str(prior[k][:2]), str(mprior[k]))
            else:
                self.assertEqual(str(prior[k][:1]), str(mprior[k]))
        self.assertEqual(
            str(mprior['Vnn']),
            str(np.array([prior['Vnn'][0], prior['Vnn'][1], prior['Vnn'][3]]))
            )
        self.assertEqual(
            str(mprior['Voo']),
            str(np.array([prior['Voo'][0]]))
            )
        self.assertEqual(
            str(mprior['Vno']),
            str(prior['Vno'][:2, :1])
            )
        self.assertTrue('Von' not in mprior)
        self.assertTrue('dummy' not in mprior)

    def test_buildprior_fitfcn(self):
        " Corr3.buildprior consistent with Corr3.fitfcn "
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo',
            tmin=2, T=5,
            )
        prior = dict(
            a=['1(1)', '2(1)', '3(1)'],
            b=['3(1)', '4(1)', '5(1)'],
            dEa=['0.1(1)', '0.2(1)', '0.3(1)'],
            dEb=['0.3(1)', '0.4(1)', '0.5(1)'],
            Vnn=3 * [['10(1)', '20(1)', '30(1)']],
            Voo=2* [['40(1)', '50(1)']],
            Von=2 * [['1(1)', '2(1)', '3(1)']],
            Vno=3 * [['3(1)', '4(1)']],
            dummy=['1(10)']
            )
        prior = gv.gvar(prior)
        for k in ['a', 'b', 'dEa', 'dEb']:
            prior[k + 'o'] = 10 * prior[k][:2]
        mprior = m.buildprior(prior)
        m.fitfcn(mprior)
        m.fitfcn(mprior, t=[1,2])

        # marginalize
        mprior = m.buildprior(prior, nterm=(2,1))
        m.fitfcn(mprior)
        m.fitfcn(mprior, t=[1,2])

        # symmetric_V=True
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Voo='Voo',
            reverse=False, symmetric_V=True, tmin=2, T=5,
            )
        prior = dict(
            a=['1(1)', '2(1)', '3(1)'],
            b=['3(1)', '4(1)', '5(1)'],
            dEa=['0.1(1)', '0.2(1)', '0.3(1)'],
            dEb=['0.3(1)', '0.4(1)', '0.5(1)'],
            Vnn=['10(1)', '20(1)', '30(1)', '100(1)', '200(1)', '300(1)'],
            Voo=['40(1)', '50(1)', '60(1)', '40(1)', '50(1)', '60(1)'],
            Vno=3 * [['3(1)', '4(1)', '5(1)']],
            Von=3 * [['30(1)', '40(1)', '50(1)']],
            dummy=['1(10)']
            )
        prior = gv.gvar(prior)
        for k in ['a', 'b', 'dEa', 'dEb']:
            prior[k + 'o'] = 10 * prior[k]
        mprior = m.buildprior(prior)
        m.fitfcn(mprior)
        m.fitfcn(mprior, t=[1,2])

        # symmetric_V=True, reverse=True
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Von='Von', Voo='Voo',
            reverse=True, symmetric_V=True, tmin=2, T=5,
            )
        prior = dict(
            a=['1(1)', '2(1)', '3(1)'],
            b=['3(1)', '4(1)', '5(1)'],
            dEa=['0.1(1)', '0.2(1)', '0.3(1)'],
            dEb=['0.3(1)', '0.4(1)', '0.5(1)'],
            Vnn=['10(1)', '20(1)', '30(1)', '100(1)', '200(1)', '300(1)'],
            Voo=['40(1)', '50(1)', '60(1)', '40(1)', '50(1)', '60(1)'],
            Vno=3 * [['3(1)', '4(1)', '5(1)']],
            Von=3 * [['30(1)', '40(1)', '50(1)']],
            dummy=['1(10)']
            )
        prior = gv.gvar(prior)
        for k in ['a', 'b', 'dEa', 'dEb']:
            prior[k + 'o'] = 10 * prior[k]
        mprior = m.buildprior(prior)
        m.fitfcn(mprior)
        m.fitfcn(mprior, t=[1,2])

        # marginalize
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Voo='Voo',
            reverse=False, symmetric_V=True, tmin=2, T=5,
            )
        mprior = m.buildprior(prior, nterm=(2,1))
        m.fitfcn(mprior)
        m.fitfcn(mprior, t=[1,2])


    def test_fitfcn(self):
        " Corr3.fitfcn "
        def prop(a, dE, t, osc=False):
            if a is None:
                return None
            t = np.array(t, dtype=int)
            E = np.cumsum(dE)
            ans = a[:, None] * np.exp(- E[:, None] * t[None, :])
            if osc:
                return ans * (-1) ** t[None, :]
            else:
                return ans
        na, nao, nb, nbo = 5, 4, 3, 2
        p = dict(
            a=na * ['1.0(1)'], ao=nao * ['2.0(1)'], b=nb * ['1.5(1)'],
            bo= nbo * ['2.5(1)'],
            dEa=na * ['0.1(1)'], dEb=nb * ['0.2(1)'], dEao=nao * ['0.15(10)'],
            dEbo=nbo * ['0.25(10)']
            )
        p['Vnn'] = gv.gvar(na * [ nb * ['1(1)']])
        p['Vno'] = gv.gvar(na * [nbo * ['2(1)']])
        p['Voo'] = gv.gvar(nao * [nbo * ['3(1)']])
        p['Von'] = gv.gvar(nao * [nb * ['4(1)']])
        p = gv.gvar(p)
        m = Corr3(
            'tag', a='a', dEa='dEa', b='b', dEb='dEb', Vnn='Vnn',
            tmin=1, T=5,
            )
        mprior = m.buildprior(p)
        def G(p, t=m.tfit, T=m.T):
            return np.sum(
                prop(p['a'], p['dEa'], t) *
                np.dot(p['Vnn'], prop(p['b'], p['dEb'], T-t)),
                axis=0
                )
        self.assertEqual(str(G(p)), str(m.fitfcn(p)))

        # osc
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo',
            tmin=2, T=5,
            )
        mprior = m.buildprior(p)
        def G(p, t=m.tfit, T=m.T, sa=m.sa, sb=m.sb):
            aprop = (
                prop(p['a'], p['dEa'], t) * sa[0],
                prop(p['ao'], p['dEao'], t, osc=True) * sa[1]
                )
            bprop = (
                prop(p['b'], p['dEb'], T-t) * sb[0],
                prop(p['bo'], p['dEbo'], T-t, osc=True) * sb[1]
                )
            V = [[p['Vnn'], p['Vno']], [p['Von'], p['Voo']]]
            ans = 0
            for i in range(2):
                for j in range(2):
                    ans += np.sum(
                        aprop[i] * np.dot(V[i][j], bprop[j]),
                        axis=0
                        )
            return ans
        self.assertEqual(str(G(p)), str(m.fitfcn(p)))

        # osc, reverse=True (has no effect on model, just on data)
        m = Corr3(
            'tag', b=('b', 'bo'), dEb=('dEb', 'dEbo'), a=('a', 'ao'),
            dEa=('dEa', 'dEao'), Vnn='Vnn', Von='Von', Vno='Vno', Voo='Voo',
            tmin=0, T=5, reverse=True
            )
        mprior = m.buildprior(p)
        self.assertEqual(str(G(p, m.tfit)), str(m.fitfcn(p)))

    def test_fitfcn_symm(self):
        " Corr3.fitfcn with symmetric_V=True "
        def prop(a, dE, t, osc=False):
            if a is None:
                return None
            t = np.array(t, dtype=int)
            E = np.cumsum(dE)
            ans = a[:, None] * np.exp(- E[:, None] * t[None, :])
            if osc:
                return ans * (-1) ** t[None, :]
            else:
                return ans
        na, nao, nb, nbo = 5, 4, 5, 4
        p = dict(
            a=na * ['1.0(1)'], ao=nao * ['2.0(1)'], b=nb * ['1.5(1)'],
            bo= nbo * ['2.5(1)'],
            dEa=na * ['0.1(1)'], dEb=nb * ['0.2(1)'], dEao=nao * ['0.15(10)'],
            dEbo=nbo * ['0.25(10)']
            )
        p['Vnn'] = gv.gvar((na * (na + 1)) // 2 * ['1(1)'])
        p['Vno'] = gv.gvar(na * [nbo * ['2(1)']])
        p['Voo'] = gv.gvar((nao * (nao + 1)) // 2 * ['3(1)'])
        p['Von'] = gv.gvar(nao * [nb * ['4(1)']])
        p = gv.gvar(p)
        m = Corr3(
            'tag', a='a', dEa='dEa', b='b', dEb='dEb', Vnn='Vnn',
            tmin=1, T=5, symmetric_V=True
            )
        mprior = m.buildprior(p)
        def build_V(V, n=na):
            ans = np.empty((n,n), dtype=V.dtype)
            for i in range(n):
                for j in range(i, n):
                    ans[i, j] = V[i * n + j - (i * (i + 1)) // 2]
                    if i != j:
                        ans[j, i] = ans[i, j]
            return ans
        def G(p, t=m.tfit, T=m.T):
            V = build_V(p['Vnn'], n=na)
            return np.sum(
                prop(p['a'], p['dEa'], t) *
                np.dot(V, prop(p['b'], p['dEb'], T-t)),
                axis=0
                )
        self.assertEqual(str(G(p)), str(m.fitfcn(p)))

        # osc
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Voo='Voo',
            tmin=1, T=5, symmetric_V=True,
            )
        mprior = m.buildprior(p)
        def G(p, t=m.tfit, T=m.T, sa=m.sa, sb=m.sb):
            aprop = (
                prop(p['a'], p['dEa'], t) * sa[0],
                prop(p['ao'], p['dEao'], t, osc=True) * sa[1]
                )
            bprop = (
                prop(p['b'], p['dEb'], T-t) * sb[0],
                prop(p['bo'], p['dEbo'], T-t, osc=True) * sb[1]
                )
            V = [[build_V(p['Vnn'], na), p['Vno']], [p['Vno'].T, build_V(p['Voo'],nao)]]
            ans = 0
            for i in range(2):
                for j in range(2):
                    ans += np.sum(
                        aprop[i] * np.dot(V[i][j], bprop[j]),
                        axis=0
                        )
            return ans
        self.assertEqual(str(G(p)), str(m.fitfcn(p)))

        # osc, reverse=True (no effect on model, just on data)
        m = Corr3(
            'tag', a=('a', 'ao'), dEa=('dEa', 'dEao'), b=('b', 'bo'),
            dEb=('dEb', 'dEbo'), Vnn='Vnn', Vno='Vno', Voo='Voo',
            tmin=0, T=5, symmetric_V=True, reverse=True,
            )
        mprior = m.buildprior(p)
        self.assertEqual(str(G(p,t=m.tfit)), str(m.fitfcn(p)))

class test_corrfitter(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_nterm(self):
        gv.ranseed(1)
        prior = dict(
            a=['1.00(1)', '1.10(1)'], ao=['1.20(1)', '1.30(1)'],
            b=['1.40(1)', '1.50(1)'], bo=['1.60(1)', '1.70(1)'],
            dE=['0.100(1)', '0.110(1)'],
            )
        prior = gv.gvar(prior)
        models = [
            Corr2('a', a='a', b='a', dE='dE', tmin=1, tmax=4),
            Corr2('b', a='b', b='b', dE='dE', tmin=1, tmax=4)
            ]
        fitter = CorrFitter(models=models, nterm=13)
        self.assertEqual(fitter.mopt, 13)
        kargs, okargs = fitter.set(nterm=12)
        self.assertEqual(fitter.mopt, 12)
        fitter.set(**okargs)
        self.assertEqual(fitter.mopt, 13)

    def test_lsqfit_2pt(self):
        " CorrFitter.lsqfit and chained_lsqfit 2pt amplitudes "
        gv.ranseed(1)
        prior = dict(
            a=['1.00(1)', '1.10(1)'], ao=['1.20(1)', '1.30(1)'],
            b=['1.40(1)', '1.50(1)'], bo=['1.60(1)', '1.70(1)'],
            dE=['0.100(1)', '0.110(1)'],
            )
        prior = gv.gvar(prior)
        models = [
            Corr2('a', a='a', b='a', dE='dE', tmin=1, tmax=4),
            Corr2('b', a='b', b='b', dE='dE', tmin=1, tmax=4)
            ]

        # fitter and fake data
        fitter = CorrFitter(models=models)
        p = fitter.buildprior(prior)
        fitfcn = fitter.buildfitfcn()
        pdata = gv.make_fake_data(fitfcn(p))

        # fit with lsqfit
        fit = fitter.lsqfit(pdata=pdata, prior=prior)
        chi2 = gv.chi2(fit.p, p)
        self.assertTrue(chi2.Q > 0.01)

        # fit with chained_lsqfit
        fit = fitter.chained_lsqfit(pdata=pdata, prior=prior)
        chi2 = gv.chi2(fit.p, p)
        self.assertTrue(chi2.Q > 0.01)

    def test_lsqfit_3pt(self):
        " CorrFitter.lsqfit and chained_lsqfit 2pt+3pt amplitudes "
        gv.ranseed(1)
        prior = dict(
            a=['1.00(1)', '1.10(1)'],
            b=['1.40(1)', '1.50(1)'],
            dEa=['0.100(1)', '0.110(1)'],
            dEb=['0.140(1)', '0.150(1)'],
            Vnn=2 * [2 * ['2.00(1)']], Vno=2 * [2 * ['2.10(1)']],
            Von=2 * [2 * ['2.20(1)']], Voo=2 * [2 * ['2.30(1)']],
            )
        prior = gv.gvar(prior)
        models = [
            Corr2('a', a='a', b='a', dE='dEa', tmin=1, tmax=4),
            Corr2('b', a='b', b='b', dE='dEb', tmin=1, tmax=4),
            Corr3(
                'ab', a='a', b='b', dEa='dEa', dEb='dEb', Vnn='Vnn',
                tmin=1, T=4,
                ),
            ]
        # fitter and fake data
        fitter = CorrFitter(models=models)
        p = fitter.buildprior(prior)
        fitfcn = fitter.buildfitfcn()
        pdata = gv.make_fake_data(fitfcn(p))

        # fit with lsqfit
        fit = fitter.lsqfit(pdata=pdata, prior=prior)
        chi2 = gv.chi2(fit.p, p)
        self.assertTrue(chi2.Q > 0.2)

        # fit with chained_lsqfit
        fit = fitter.chained_lsqfit(pdata=pdata, prior=prior)
        chi2 = gv.chi2(fit.p, p)
        self.assertTrue(chi2.Q > 0.2)

    def test_lsqfit_3pt_symm(self):
        " CorrFitter.lsqfit and chained_lsqfit 2pt+3pt amplitudes symm_V "
        gv.ranseed(1)
        prior = dict(
            a=['1.00(1)', '1.10(1)'],
            b=['1.40(1)', '1.50(1)'],
            dEa=['0.100(1)', '0.110(1)'],
            dEb=['0.140(1)', '0.150(1)'],
            Vnn=3 * ['2.00(1)'], Vno=2 * [2 * ['2.10(1)']],
            Von=2 * [2 * ['2.20(1)']], Voo=3 * ['2.30(1)'],
            )
        prior = gv.gvar(prior)
        models = [
            Corr2('a', a='a', b='a', dE='dEa', tmin=1, tmax=4),
            Corr2('b', a='b', b='b', dE='dEb', tmin=1, tmax=4),
            Corr3(
                'ab', a='a', b='b', dEa='dEa', dEb='dEb', Vnn='Vnn',
                tmin=1, T=4, symmetric_V=True
                ),
            ]
        # fitter and fake data
        fitter = CorrFitter(models=models)
        p = fitter.buildprior(prior)
        fitfcn = fitter.buildfitfcn()
        pdata = gv.make_fake_data(fitfcn(p))

        # fit with lsqfit
        fit = fitter.lsqfit(pdata=pdata, prior=prior)
        chi2 = gv.chi2(fit.p, p)
        self.assertTrue(chi2.Q > 0.2)

        # fit with chained_lsqfit
        fit = fitter.chained_lsqfit(pdata=pdata, prior=prior)
        chi2 = gv.chi2(fit.p, p)
        self.assertTrue(chi2.Q > 0.2)

    def test_lsqfit_3pt_osc(self):
        " CorrFitter.lsqfit and chained_lsqfit 2pt+3pt amplitudes "
        gv.ranseed(1)
        prior = dict(
            a=['1.00(1)', '1.10(1)'], ao=['2.00(1)', '2.10(1)'],
            b=['1.40(1)', '1.50(1)'], bo=['2.40(1)', '2.50(1)'],
            dEa=['0.100(1)', '0.110(1)'], dEao=['0.200(1)', '0.210(1)'],
            dEb=['0.140(1)', '0.150(1)'], dEbo=['0.240(1)', '0.250(1)'],
            Vnn=2 * [2 * ['2.00(1)']], Vno=2 * [2 * ['2.10(1)']],
            Von=2 * [2 * ['2.20(1)']], Voo=2 * [2 * ['2.30(1)']],
            )
        prior = gv.gvar(prior)
        models = [
            Corr2(
                'a', a=('a','ao'), b=('a','ao'), dE=('dEa','dEao'),
                tmin=1, tmax=4
                ),
            Corr2(
                'b', a=('b','bo'), b=('b','bo'), dE=('dEb','dEbo'),
                tmin=1, tmax=4
                ),
            Corr3(
                'ab', a=('a','ao'), b=('b','bo'), dEa=('dEa','dEao'),
                dEb=('dEb','dEbo'), Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo',
                tmin=1, T=4,
                ),
            ]
        # fitter and fake data
        fitter = CorrFitter(models=models)
        p = fitter.buildprior(prior)
        fitfcn = fitter.buildfitfcn()
        pdata = gv.make_fake_data(fitfcn(p))

        # fit with lsqfit
        fit = fitter.lsqfit(pdata=pdata, prior=prior)
        chi2 = gv.chi2(fit.p, p)
        self.assertTrue(chi2.Q > 0.1)

        # fit with chained_lsqfit
        fit = fitter.chained_lsqfit(pdata=pdata, prior=prior)
        chi2 = gv.chi2(fit.p, p)
        self.assertTrue(chi2.Q > 0.1)


    def test_lsqfit_3pt_osc_symm(self):
        " CorrFitter.lsqfit and chained_lsqfit 2pt+3pt amplitudes symm_V=True "
        gv.ranseed(1)
        prior = dict(
            a=['1.00(1)', '1.10(1)'], ao=['2.00(1)', '2.10(1)'],
            b=['1.40(1)', '1.50(1)'], bo=['2.40(1)', '2.50(1)'],
            dEa=['0.100(1)', '0.110(1)'], dEao=['0.200(1)', '0.210(1)'],
            dEb=['0.140(1)', '0.150(1)'], dEbo=['0.240(1)', '0.250(1)'],
            Vnn=3 * ['2.00(1)'], Vno=2 * [2 * ['2.10(1)']],
            Von=2 * [2 * ['2.20(1)']], Voo=3 * ['2.30(1)'],
            )
        prior = gv.gvar(prior)
        models = [
            Corr2(
                'a', a=('a','ao'), b=('a','ao'), dE=('dEa','dEao'),
                tmin=1, tmax=4
                ),
            Corr2(
                'b', a=('b','bo'), b=('b','bo'), dE=('dEb','dEbo'),
                tmin=1, tmax=4
                ),
            Corr3(
                'ab', a=('a','ao'), b=('b','bo'), dEa=('dEa','dEao'),
                dEb=('dEb','dEbo'), Vnn='Vnn', Vno='Vno', Voo='Voo',
                tmin=1, T=4, symmetric_V=True
                ),
            ]
        # fitter and fake data
        fitter = CorrFitter(models=models)
        p = fitter.buildprior(prior)
        fitfcn = fitter.buildfitfcn()
        pdata = gv.make_fake_data(fitfcn(p))

        # fit with lsqfit
        fit = fitter.lsqfit(pdata=pdata, prior=prior)
        chi2 = gv.chi2(fit.p, p)
        self.assertTrue(chi2.Q > 0.1)

        # fit with chained_lsqfit
        fit = fitter.chained_lsqfit(pdata=pdata, prior=prior)
        chi2 = gv.chi2(fit.p, p)
        self.assertTrue(chi2.Q > 0.1)

class test_eigenbasis(unittest.TestCase):
    def setUp(self):
        self.u = np.array([[1., 0.5], [1., 2.]])
        self.E = np.array([1., 3.])
        self.f = np.array([1., 1.])

    def make_G(self, tdata, keyfmt, srcs, f=None):
        G = collections.OrderedDict()
        tdata = np.array(tdata)
        if f is None:
            f = self.f
        else:
            f = np.array(f)
        for i, s1 in enumerate(srcs):
            for j, s2 in enumerate(srcs):
                key = keyfmt.format(s1=s1, s2=s2)
                G[key] = 0
                for n in range(2):
                    G[key] += (f[n] ** tdata) * (
                        self.u[n, i] * self.u[n, j] * gv.exp(-self.E[n] * tdata)
                        * gv.gvar('1.00(1)')
                        )
        return G

    def test_apply(self):
        " EigenBasis EigenBasis.apply EigenBasis.unapply "
        for tdata in [
            [1., 2., 3., 4.],
            [2., 4., 6., 8.],
            [0, 1., 2.],
            ]:
            tdata = np.array(tdata)
            G = self.make_G(tdata, keyfmt='{s1}{s2}', srcs='ab')
            basis = EigenBasis(
                data=G, keyfmt='{s1}{s2}', srcs='ab',
                t=2, tdata=tdata,
                )
            np.testing.assert_allclose(basis.E, self.E)
            newG = basis.apply(G, '{s1}{s2}')
            newG_mean = gv.mean(newG)
            np.testing.assert_allclose(newG_mean['00'], gv.exp(-self.E[0] * tdata))
            np.testing.assert_allclose(newG_mean['11'], gv.exp(-self.E[1] * tdata))
            np.testing.assert_allclose(newG_mean['01'], 0, atol=1e-10)
            np.testing.assert_allclose(newG_mean['10'], 0, atol=1e-10)
            oldG = basis.unapply(newG, '{s1}{s2}')
            for k in ['aa', 'ab', 'ba', 'bb']:
                np.testing.assert_allclose(gv.mean(oldG[k] - G[k]), 0, atol=1e-10)
                np.testing.assert_allclose(gv.sdev(oldG[k] - G[k]), 0, atol=1e-10)

    def test_svd(self):
        " EigenBasis.svd "
        tdata = [1,2,3,4]
        G = self.make_G(tdata, keyfmt='{s1}{s2}', srcs='ab')
        basis = EigenBasis(
            data=G, keyfmt='{s1}{s2}', srcs='ab',
            t=2, tdata=tdata,
            )
        Gsvd = basis.svd(G, svdcut=0.9)
        self.assertEqual(basis.svdn, 15)
        self.assertEqual(str(basis.svdcorrection), '0.000(30)')
        for k in G:
            np.testing.assert_allclose(gv.mean(G[k]), gv.mean(Gsvd[k]))
            self.assertTrue(np.all(gv.sdev(Gsvd[k]) > gv.sdev(G[k])))

    def test_make_prior(self):
        " EigenBasis.make_prior "
        tdata = np.arange(4.)
        datafmt = 'G.{s1}.{s2}'
        srcs = 'ab'
        G = self.make_G(tdata=tdata, keyfmt=datafmt, srcs=srcs)
        basis = EigenBasis(
            data=G, keyfmt=datafmt, srcs=srcs, t=(1,2),
            )
        nterm = 4
        ampl = '1.0(2)', '0.03(20)', '0.2(2.0)'
        dEfac = '1(1)'
        one, small, big = ampl

        # case 1 - canonical
        prior = basis.make_prior(
            nterm=nterm, keyfmt='a.{s1}', ampl=ampl, dEfac=dEfac, eig_srcs=True
            )
        dE = gv.gvar(nterm * [dEfac]) * (self.E[1] - self.E[0])
        dE[0] += self.E[0] - dE[0].mean
        self.assertEqual(str(gv.exp(prior['log(a.dE)'])), str(dE))
        a0 = gv.gvar(nterm * [big])
        a0[0] = gv.gvar(one)
        a0[1] = gv.gvar(small)
        self.assertEqual(str(prior['a.0']), str(a0))
        a1 = gv.gvar(nterm * [big])
        a1[0] = gv.gvar(small)
        a1[1] = gv.gvar(one)
        self.assertEqual(str(prior['a.1']), str(a1))

        # default values
        ampl = '1.0(3)', '0.03(10)', '0.2(1.0)'
        dEfac = '1(1)'
        one, small, big = ampl

        # case 2 - omit state
        prior = basis.make_prior(nterm=nterm, keyfmt='a.{s1}', states=[0], eig_srcs=True)
        a0 = gv.gvar(nterm * [big])
        a0[0] = gv.gvar(one)
        self.assertEqual(str(prior['a.0']), str(a0))
        a1 = gv.gvar(nterm * [big])
        a1[0] = gv.gvar(small)
        self.assertEqual(str(prior['a.1']), str(a1))

        # case 3 - swap states
        prior = basis.make_prior(nterm=nterm, keyfmt='a.{s1}', states=[1, 0], eig_srcs=True)
        dE = gv.gvar(nterm * ['1(1)']) * (self.E[1] - self.E[0])
        dE[0] += self.E[0] - dE[0].mean
        self.assertEqual(str(gv.exp(prior['log(a.dE)'])), str(dE))
        a1 = gv.gvar(nterm * [big])
        a1[0] = gv.gvar(one)
        a1[1] = gv.gvar(small)
        self.assertEqual(str(prior['a.1']), str(a1))
        a0 = gv.gvar(nterm * [big])
        a0[0] = gv.gvar(small)
        a0[1] = gv.gvar(one)
        self.assertEqual(str(prior['a.0']), str(a0))

    def test_make_prior_osc_osc(self):
        " EigenBasis.make_prior "
        tdata = np.arange(4.)
        datafmt = 'G.{s1}.{s2}'
        srcs = 'ab'
        G = self.make_G(tdata=tdata, keyfmt=datafmt, srcs=srcs)
        basis = EigenBasis(
            data=G, keyfmt=datafmt, srcs=srcs, t=(1,2),
            )
        nterm = 4
        ampl = '1.0(2)', '0.03(20)', '0.2(2.0)'
        dEfac = '1(1)'
        one, small, big = ampl

        # case 1 - canonical
        prior = basis.make_prior(
            nterm=nterm, keyfmt='a.{s1}', ampl=ampl, dEfac=dEfac, eig_srcs=True
            )
        dE = gv.gvar(nterm * [dEfac]) * (self.E[1] - self.E[0])
        dE[0] += self.E[0] - dE[0].mean
        self.assertEqual(str(gv.exp(prior['log(a.dE)'])), str(dE))
        a0 = gv.gvar(nterm * [big])
        a0[0] = gv.gvar(one)
        a0[1] = gv.gvar(small)
        self.assertEqual(str(prior['a.0']), str(a0))
        a1 = gv.gvar(nterm * [big])
        a1[0] = gv.gvar(small)
        a1[1] = gv.gvar(one)
        self.assertEqual(str(prior['a.1']), str(a1))

        # default values
        ampl = '1.0(3)', '0.03(10)', '0.2(1.0)'
        dEfac = '1(1)'
        one, small, big = ampl

        # case 2 - omit state
        prior = basis.make_prior(nterm=nterm, keyfmt='a.{s1}', states=[0], eig_srcs=True)
        a0 = gv.gvar(nterm * [big])
        a0[0] = gv.gvar(one)
        self.assertEqual(str(prior['a.0']), str(a0))
        a1 = gv.gvar([small] + (nterm - 1) * [big])
        self.assertEqual(str(prior['a.1']), str(a1))

        # case 3 - swap states
        prior = basis.make_prior(nterm=nterm, keyfmt='a.{s1}', states=[1, 0], eig_srcs=True)
        dE = gv.gvar(nterm * ['1(1)']) * (self.E[1] - self.E[0])
        dE[0] += self.E[0] - dE[0].mean
        self.assertEqual(str(gv.exp(prior['log(a.dE)'])), str(dE))
        a1 = gv.gvar(nterm * [big])
        a1[0] = gv.gvar(one)
        a1[1] = gv.gvar(small)
        self.assertEqual(str(prior['a.1']), str(a1))
        a0 = gv.gvar(nterm * [big])
        a0[0] = gv.gvar(small)
        a0[1] = gv.gvar(one)
        self.assertEqual(str(prior['a.0']), str(a0))

    def test_make_prior_osc(self):
        " EigenBasis.make_prior "
        osc = True
        tdata = np.arange(4.)
        datafmt = 'G.{s1}.{s2}'
        srcs = 'ab'
        G = self.make_G(tdata=tdata, keyfmt=datafmt, srcs=srcs, f=[-1., 1.])
        basis = EigenBasis(
            data=G, keyfmt=datafmt, srcs=srcs, t=(1,2), osc=osc
            )
        nterm = 4
        ampl = '1.0(2)', '0.03(20)', '0.2(2.0)'
        dEfac = '1(1)'
        one, small, big = ampl

        # case 1 - canonical
        prior = basis.make_prior(
            nterm=nterm, keyfmt=('a.{s1}','ao.{s1}'),
            ampl=ampl, dEfac=dEfac, eig_srcs=True, #states=[0,1,2]
            )
        self.assertEqual(
            str(gv.exp(prior['log(a.dE)'])),
            str(self.E[1] * np.array(gv.gvar(nterm * [dEfac])))
            )
        self.assertEqual(
            str(gv.exp(prior['log(ao.dE)'])),
            str(self.E[0] * np.array(gv.gvar(nterm * [dEfac])))
            )
        a0 = gv.gvar(nterm * [big])
        a0[0] = gv.gvar(one)
        self.assertEqual(str(prior['a.0']), str(a0))
        self.assertEqual(str(prior['ao.0']), str(a0))
        a1 = gv.gvar(nterm * [big])
        a1[0] = gv.gvar(small)
        self.assertEqual(str(prior['a.1']), str(a1))
        self.assertEqual(str(prior['ao.1']), str(a1))

class test_read_dataset(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_text(self):
        # make text file
        s = [1., 2., 3., 4.]
        v = list(np.array([[10.,11.], [12., 13.], [14., 15.], [16., 17.]]))
        ref_dset = dict(s=s, v=v)
        filetext = '\n'.join([
            '# comment',
            's 1.', 's 2.', 's 3.', 's 4.',
            'v 10. 11.', 'v 12. 13.', 'v 14. 15.', 'v 16. 17.'
            ])
        with open('test-gvar.txt', 'w') as textfile:
            textfile.write(filetext)
        # everything
        dset = read_dataset('test-gvar.txt')
        self.assertEqual(list(dset.keys()), ['s', 'v'])
        for k in dset:
            self.assertEqual(str(dset[k]), str(ref_dset[k]))
        # s only
        dset = read_dataset('test-gvar.txt', grep='[^v]')
        self.assertEqual(list(dset.keys()), ['s'])
        for k in ['s']:
            self.assertEqual(str(dset[k]), str(ref_dset[k]))
        # v only
        dset = read_dataset('test-gvar.txt', keys=['v'])
        self.assertEqual(list(dset.keys()), ['v'])
        for k in ['v']:
            self.assertEqual(str(dset[k]), str(ref_dset[k]))
        # binsize=2
        dset = read_dataset('test-gvar.txt', binsize=2)
        self.assertEqual(list(dset.keys()), ['s', 'v'])
        self.assertEqual(dset['s'], [1.5, 3.5])
        self.assertEqual(
            str(dset['v']),
            str([np.array([11., 12.]), np.array([15., 16.])])
            )
        os.remove('test-gvar.txt')

    @unittest.skipIf(NO_H5PY,"skipping test_hdf5 --- no h5py modules")
    def test_hdf5(self):
        # make hdf5 file
        s = [1., 2., 3., 4.]
        v = list(np.array([[10.,11.], [12., 13.], [14., 15.], [16., 17.]]))
        ref_dset = dict(s=s, v=v)
        with h5py.File('test-gvar.h5', 'w') as h5file:
            h5file['/run1/s'] = s
            h5file['/run2/v'] = v
        # everything
        dset = read_dataset('test-gvar.h5', h5group=['/run1', '/run2'])
        self.assertEqual(list(dset.keys()), ['s', 'v'])
        for k in dset:
            self.assertEqual(str(dset[k]), str(ref_dset[k]))
        # s only
        dset = read_dataset('test-gvar.h5', h5group=['/run1', '/run2'], grep='[^v]')
        self.assertEqual(list(dset.keys()), ['s'])
        for k in ['s']:
            self.assertEqual(str(dset[k]), str(ref_dset[k]))
        # v only
        dset = read_dataset('test-gvar.h5', h5group=['/run1', '/run2'], keys=['v'])
        self.assertEqual(list(dset.keys()), ['v'])
        for k in ['v']:
            self.assertEqual(str(dset[k]), str(ref_dset[k]))
        # binsize=2
        dset = read_dataset('test-gvar.h5', h5group=['/run1', '/run2'], binsize=2)
        self.assertEqual(list(dset.keys()), ['s', 'v'])
        self.assertEqual(dset['s'], [1.5, 3.5])
        self.assertEqual(
            str(dset['v']),
            str([np.array([11., 12.]), np.array([15., 16.])])
            )
        os.remove('test-gvar.h5')

class test_fastfit(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_periodic(self):
        gv.ranseed(1234)
        m = Corr2('tag', a='a', b='a', dE='dE', tmin=1, tp=32)
        prior = gv.gvar(dict(
            a=['2.00(1)', '2.00(1)'], dE=['0.50(1)', '0.50(1)'],
            ))
        G = gv.make_fake_data(m.fitfcn(prior, t=m.tdata), 0.1)
        fit = fastfit(G, ampl='4(2)', dE='0.5(5)', tmin=5, tp=32)
        self.assertTrue(abs(fit.E.mean - 0.5) < 5 * fit.E.sdev)
        self.assertTrue(abs(fit.ampl.mean - 4) < 5 * fit.ampl.sdev)

    def test_antiperiodic(self):
        gv.ranseed(1234)
        m = Corr2('tag', a='a', b='a', dE='dE', tmin=1, tp=-31)
        prior = gv.gvar(dict(
            a=['2.00(1)', '2.00(1)'], dE=['0.50(1)', '0.50(1)'],
            ))
        G = gv.make_fake_data(m.fitfcn(prior, t=m.tdata), 0.1)
        fit = fastfit(G, ampl='4(2)', dE='0.5(5)', tmin=5, tp=-31)
        self.assertTrue(abs(fit.E.mean - 0.5) < 5 * fit.E.sdev)
        self.assertTrue(abs(fit.ampl.mean - 4) < 5 * fit.ampl.sdev)

    def test_nonperiodic(self):
        gv.ranseed(12345)
        m = Corr2('tag', a='a', b='a', dE='dE', tmin=1, tmax=16)
        prior = gv.gvar(dict(
            a=['2.00(1)', '2.00(1)'], dE=['0.50(1)', '0.50(1)'],
            ))
        G = gv.make_fake_data(m.fitfcn(prior, t=m.tdata), 0.1)
        fit = fastfit(G, ampl='4(2)', dE='0.5(5)', tmin=6)
        self.assertTrue(abs(fit.E.mean - 0.5) < 5 * fit.E.sdev)
        self.assertTrue(abs(fit.ampl.mean - 4) < 5 * fit.ampl.sdev)

    def test_periodic_osc(self):
        gv.ranseed(1234)
        m = Corr2(
            'tag', a=(None, 'a'), b=(None, 'a'), dE=(None, 'dE'),
            tmin=1, tp=32, s=(0,1),
            )
        prior = gv.gvar(dict(
            a=['2.00(1)', '2.00(1)'], dE=['0.50(1)', '0.50(1)'],
            ))
        G = gv.make_fake_data(m.fitfcn(prior, t=m.tdata), 0.1)
        fit = fastfit(
            G, ampl=(None, '4(2)'), dE=(None, '0.5(5)'), tmin=5, tp=32,
            osc=True, s=(0, 1)
            )
        self.assertTrue(abs(fit.E.mean - 0.5) < 5 * fit.E.sdev)
        self.assertTrue(abs(fit.ampl.mean - 4) < 5 * fit.ampl.sdev)

if __name__ == '__main__':
    unittest.main()
