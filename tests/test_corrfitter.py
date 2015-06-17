#!/usr/bin/env python
# encoding: utf-8
"""
test_corrfitter.py
"""
# Copyright (c) 2012-2014 G. Peter Lepage.
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

import unittest
import inspect
import numpy as np
import gvar as gv
import gvar.dataset as ds
from corrfitter import *

PRINT_FITS = False      # print lots of fit info while doing tests
DISPLAY_PLOTS = False   # display plots for some fits
NSIG = 6.0              # number of sigmas allowed before signalling an error
NTERM = 3               # number of terms (ie, len(dE)) in correlators
SVDCUT = 1e-4           # svd cut used for all fits -- to minimize roundoff problems
FAST = False            # skips bootstrap tests if True

class ArrayTests(object):
    def __init__(self):
        pass
    ##
    def assert_gvclose(self,x,y,rtol=1e-5,atol=1e-8,prt=False):
        """ asserts that the means and sdevs of all x and y are close """
        if hasattr(x,'keys') and hasattr(y,'keys'): 
            if sorted(x.keys())==sorted(y.keys()):
                for k in x:
                    self.assert_gvclose(x[k],y[k],rtol=rtol,atol=atol)
                return
            else:
                raise ValueError("x and y have mismatched keys")
        self.assertSequenceEqual(np.shape(x),np.shape(y))
        x = np.asarray(x).flat
        y = np.asarray(y).flat
        if prt:
            print(np.array(x))
            print(np.array(y))
        for xi,yi in zip(x,y):
            self.assertGreater(atol+rtol*abs(yi.mean),abs(xi.mean-yi.mean))
            self.assertGreater(10*(atol+rtol*abs(yi.sdev)),abs(xi.sdev-yi.sdev))
    ##
    def assert_arraysclose(self,x,y,rtol=1e-5,prt=False):
        self.assertSequenceEqual(np.shape(x),np.shape(y))
        x = np.array(x).flatten()
        y = np.array(y).flatten()
        max_val = max(np.abs(list(x)+list(y)))
        max_rdiff = max(np.abs(x-y))/max_val
        if prt:
            print(x)
            print(y)
            print(max_val,max_rdiff,rtol)
        self.assertAlmostEqual(max_rdiff,0.0,delta=rtol)
    ##
    def assert_arraysequal(self,x,y):
        self.assertSequenceEqual(np.shape(x),np.shape(y))
        x = [float(xi) for xi in np.array(x).flatten()]
        y = [float(yi) for yi in np.array(y).flatten()]
        self.assertSequenceEqual(x,y)
    ##
##

class FitTests(object):
    def __init__(self):
        pass
    ##
    def assert_fitsagree(self, fitp1, fitp2):
        """ fit outputs fitp1 and fitp2 agree with each other (approximately)"""
        chi2 = 0
        dof = 0
        for k in fitp1:
            p1 = fitp1[k].flat
            p2 = fitp2[k].flat
            for p1i, p2i in zip(p1,p2):
                delta_mean = abs(abs(p1i.mean) - abs(p2i.mean))
                sdev = p1i.sdev
                self.assertLess(delta_mean, NSIG * sdev)
                delta_sdev = abs(p1i.sdev - p2i.sdev)
                self.assertLess(delta_sdev, NSIG * sdev)
    ##
    def assert_fitclose(self, fitp, p):
        """ GVars in fitp agree within NSIG sigma with parameters p. """
        chi2 = 0
        dof = 0
        for k in fitp:
            for fpi,pi in zip(fitp[k].flat,p[k].flat):
                delta_mean = abs(abs(fpi.mean)-abs(pi))
                sdev = fpi.sdev
                self.assertLess(delta_mean, NSIG * sdev)
                dof += 1
                chi2 += delta_mean**2/sdev**2
        if PRINT_FITS:
            print("param chi2/dof = %.2g [%d]\n" % (chi2/dof, dof))
    ##
##
    
class test_corr2(unittest.TestCase, FitTests, ArrayTests):
    def setUp(self):
        ## prior ##
        self.prior = gv.BufferDict()
        nt = NTERM
        self.prior['a'] = gv.gvar(nt*["0.50(1)"])
        self.prior['ao'] = gv.gvar(nt*["0.250(5)"])
        self.prior['logb'] = gv.log(gv.gvar(nt*["0.60(1)"]))
        self.prior['bo'] = gv.gvar(nt*["0.30(1)"])
        self.prior['logdE'] = gv.log(gv.gvar(nt*["0.50(1)"]))
        self.prior['logdEo'] = gv.log(gv.gvar(nt*["0.60(1)"]))
        ##
        ## actual parameters, time ranges, corr counter ##
        self.p = next(gv.raniter(self.prior))
        for x in ['b', 'dE', 'dEo']:
            self.p[x] = gv.exp(self.p['log' + x])
        self.tp = 10.
        self.tdata = np.arange(self.tp)
        self.tfit = self.tdata[1:]
        self.ncorr = 0
        ##
        self.ran = gv.gvar(0,1)
    ##
    def tearDown(self):
        del self.prior
        del self.p
        del self.tdata
        del self.tfit
        del self.tp
        del self.ncorr
    ##
    def getdoc(self):
        """ get __doc__ string for current function """
        frame = inspect.currentframe()
        caller_frame = inspect.getouterframes(frame)[1][0]
        caller_name = inspect.getframeinfo(caller_frame).function
        caller_func = eval("test_corr2."+caller_name)
        return caller_func.__doc__
    ##
    def mkcorr(self,a, b, dE, tp=None, othertags=[], s=1.):
        ans = Corr2(datatag=self.ncorr, a=a, b=b, dE=dE, tdata=self.tdata,
            tfit=self.tfit, tp=tp, s=s, othertags=othertags)
        self.ncorr += 1
        return ans
    ##
    def dofit(self, models, data_models=None, nterm=None, ratio=True):
        fitter = CorrFitter(models=models, nterm=nterm, ratio=ratio)
        if data_models is None:
            data_models = models
        data = make_data(models=data_models, p=self.p)
        if PRINT_FITS:
            print("Data:\n", data,"\n")
        fit = fitter.lsqfit(data=data, prior=self.prior, debug=True,
            print_fit=PRINT_FITS, svdcut=SVDCUT)
        #
        if PRINT_FITS:
            print("Corr2 exact parameter values:")
            for k in fit.p:
                print("%10s:"%str(k),self.p[k])
            print()
        self.assert_fitclose(fit.p,self.p)
        self.data = data
        self.fit = fit
        return fitter
    ##
    def dofit_chd(
        self, models, data_models=None, nterm=None, ratio=True,
        parallel=False, flat=False, fast=True):
        fitter = CorrFitter(models=models, nterm=nterm, ratio=ratio)
        if data_models is None:
            data_models = fitter.flat_models
        data = make_data(models=data_models, p=self.p)
        if PRINT_FITS:
            print("Data:\n", data,"\n")
        fit = fitter.chained_lsqfit(
            data=data, prior=self.prior, debug=True,
            print_fit=PRINT_FITS, svdcut=SVDCUT,
            parallel=parallel, flat=flat, fast=fast
            )
        #
        if PRINT_FITS:
            print("Corr2 exact parameter values:")
            for k in fit.p:
                print("%10s:"%str(k),self.p[k])
            print()
        self.assert_fitclose(fit.p, self.p)
        self.data = data
        self.fit = fit
        return fitter
    ##
    def dofastfit(self, model, osc=False):
        data = make_data(models=[model], p=self.p)
        fit = fastfit(data=data, prior=self.prior, model=model, osc=osc, 
                       svdcut=SVDCUT)
        Eeff = fit.E
        Aeff = fit.ampl
        i = 0 if not osc else 1
        Etrue = unpack(self.p, model.dE)[i][0]
        Atrue = unpack(self.p, model.a)[i][0] * unpack(self.p, model.b)[i][0]
        self.assertLess(abs(Eeff.mean - Etrue), NSIG * Eeff.sdev)
        self.assertLess(abs(Aeff.mean - Atrue), NSIG * Aeff.sdev)
        if PRINT_FITS:
            print("Eeff chi2/dof = %.2g [1]\n" % ((Eeff.mean - Etrue)**2/Eeff.sdev**2))
            print("Aeff chi2/dof = %.2g [1]\n" % ((Aeff.mean - Atrue)**2/Aeff.sdev**2))

    def test_read_dataset(self):
        " read_dataset "
        # format 1
        dset = read_dataset('format1-file')
        np.testing.assert_allclose( 
            dset['aa'], 
            [np.array([ 1.237,  0.912,  0.471]), np.array([ 1.035,  0.851,  0.426])]
            )
        np.testing.assert_allclose( 
            dset['bb'], 
            [np.array([ 3.214,  0.535,  0.125]), np.array([ 2.951,  0.625,  0.091])]
            )
        # format 2
        dset = read_dataset(dict(aa='format2-file.aa', bb='format2-file.bb'))
        np.testing.assert_allclose( 
            dset['aa'], 
            [np.array([ 1.237,  0.912,  0.471]), np.array([ 1.035,  0.851,  0.426])]
            )
        np.testing.assert_allclose( 
            dset['bb'], 
            [np.array([ 3.214,  0.535,  0.125]), np.array([ 2.951,  0.625,  0.091])]
            )

    def test_simulation(self):
        """ CorrFitter.simulated_data_iter """
        models = [ self.mkcorr(a="a", b="a", dE="logdE", tp=None) ]
        fitter = self.dofit(models)
        data = self.data
        diter = gv.BufferDict()
        k = list(data.keys())[0]
        # make n config dataset corresponding to data
        n = 100
        diter = gv.raniter(
            g = gv.gvar(gv.mean(self.data[k]), gv.evalcov(self.data[k]) * n),
            n = n
            )
        dataset = gv.dataset.Dataset()
        for d in diter:
            dataset.append(k, d)
        pexact = fitter.fit.pmean
        covexact = gv.evalcov(gv.dataset.avg_data(dataset)[k])
        for sdata in fitter.simulated_data_iter(n=2, dataset=dataset):
            sfit = fitter.lsqfit(
                data=sdata, prior=self.prior, p0=pexact, print_fit=False
                )
            diff = dict()
            for i in ['a', 'logdE']:
                diff[i] = sfit.p[i][0] - pexact[i][0]
            c2 = gv.chi2(diff)
            self.assertLess(c2/c2.dof, 15.)
            self.assert_arraysclose(gv.evalcov(sdata[k]), covexact)

    def test_periodic(self):
        """ corr2 -- periodic correlator """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a="a", b="a", dE="logdE", tp=self.tp) ]
        fitter = self.dofit(models)
        fitter = self.dofit_chd(models)
        if DISPLAY_PLOTS:
            fitter.display_plots()
    ##
    def test_fastfit_periodic(self):
        """ corr2 -- fastfit(periodic) """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        self.dofastfit(self.mkcorr(a="a", b="a", dE="logdE", tp=self.tp))
    ##        
    def test_lognormal(self):
        """ corr2 -- log normal parameters """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a="logb", b="logb", dE="logdE", tp=self.tp) ]
        self.dofit(models)
        self.dofit_chd(models)
    ##
    def test_nonperiodic(self):
        """ corr2 -- non-periodic correlator"""
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a="a", b="a", dE="logdE", tp=None) ]
        self.dofit(models)
        self.dofit_chd(models) 
    ##
    def test_fastfit_nonperiodic(self):
        """ corr2 -- fastfit(non-periodic) """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        self.dofastfit(self.mkcorr(a="a", b="a", dE="logdE", tp=None))
    ##        
    @unittest.skipIf(FAST,"skipping test_bootstrap for speed")
    def test_bootstrap(self):
        """ corr2 -- bootstrap """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a="a", b="a", dE="logdE", tp=None) ]
        fitter = self.dofit(models)
        bsdata = ds.Dataset()
        for fit in fitter.bootstrap_iter(n=40):
            bsdata.append(fit.pmean)
        bsfit = ds.avg_data(bsdata,bstrap=True)
        self.assert_fitclose(bsfit, self.p)
        self.assert_fitsagree(fitter.fit.p, bsfit)
    ##
    def test_antiperiodic(self):
        """ corr2 -- anti-periodic correlator """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a="a", b="a", dE="logdE", tp=-self.tp) ]
        self.dofit(models)
        self.dofit_chd(models)
    ##
    def test_fastfit_antiperiodic(self):
        """ corr2 -- fastfit(anti-periodic) """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        self.dofastfit(self.mkcorr(a="a", b="a", dE="logdE", tp=-self.tp))
    ##        
    def test_matrix1(self):
        """ corr2 -- 2x2 matrix fit (use othertags) """
        global NSIG
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a="a", b="a", dE="logdE"),
                   self.mkcorr(a="logb", b="logb", dE="logdE"),
                   self.mkcorr(a="a", b="logb", dE="logdE", othertags=[3]),
                   self.mkcorr(a="logb", b="a", dE="logdE", othertags=[2])]
        self.dofit(models=models[:-1], data_models=models)
        NSIG *= 1.5
        self.dofit_chd(models)
        NSIG /= 1.5
        self.assertEqual(models[2].all_datatags, [2, 3])
        self.assertEqual(models[3].all_datatags, [3, 2])
    ##
    def test_matrix2(self):
        """ corr2 -- 2x2 matrix fit (without othertags) """
        global NSIG
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a="a", b="a", dE="logdE"),
                   self.mkcorr(a="logb", b="logb", dE="logdE"),
                   self.mkcorr(a="a", b="logb", dE="logdE"),
                   self.mkcorr(a="logb", b="a", dE="logdE")]
        self.dofit(models)
        NSIG *= 2.
        self.dofit_chd(models)
        NSIG /= 2.
    ##
    def test_chained(self):
        """ test chained fit variations """
        global NSIG, PRINT_FITS
        models =[ 
                self.mkcorr(a='a', b='a', dE='logdE'),
                [
                self.mkcorr(a='a', b='b', dE='logdE'),
                self.mkcorr(a='b', b='a', dE='logdE')
                ],
                self.mkcorr(a='b', b='b', dE='logdE')
                ] 
        # PRINT_FITS = True
        NSIG *= 2.   # do this because so many fits
        self.dofit_chd(models, parallel=True, flat=False, fast=True)
        self.dofit_chd(models, parallel=True, flat=False, fast=False)
        self.dofit_chd(models, parallel=True, flat=True, fast=True)
        self.dofit_chd(models, parallel=True, flat=True, fast=False)
        self.dofit_chd(models, parallel=False, flat=False, fast=True)
        self.dofit_chd(models, parallel=False, flat=False, fast=False)
        self.dofit_chd(models, parallel=False, flat=True, fast=True)
        self.dofit_chd(models, parallel=False, flat=True, fast=False)

        self.dofit_chd(models, parallel=True, flat=False, fast=True, nterm=2)
        self.dofit_chd(models, parallel=True, flat=False, fast=False, nterm=2)
        self.dofit_chd(models, parallel=True, flat=True, fast=True, nterm=2)
        self.dofit_chd(models, parallel=True, flat=True, fast=False, nterm=2)
        self.dofit_chd(models, parallel=False, flat=False, fast=True, nterm=2)
        self.dofit_chd(models, parallel=False, flat=False, fast=False, nterm=2)
        self.dofit_chd(models, parallel=False, flat=True, fast=True, nterm=2)
        self.dofit_chd(models, parallel=False, flat=True, fast=False, nterm=2)
        NSIG /= 2.
        # PRINT_FITS = False

    def test_marginalization(self):
        """ corr2 -- marginalization (2x2 matrix)"""
        global NSIG
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a="a", b="a", dE="logdE"),
                   self.mkcorr(a="logb", b="logb", dE="logdE"),
                   self.mkcorr(a="a", b="logb", dE="logdE"),
                   self.mkcorr(a="logb", b="a", dE="logdE")]
        self.dofit(models, nterm=1, ratio=True)
        NSIG *= 1.5
        self.dofit_chd(models)
        NSIG /= 1.5
    ##
    def test_oscillating1(self):
        """ corr2 -- oscillating part """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a=("a","ao"), b=("a","ao"), 
                   dE=("logdE","logdEo")) ]
        self.dofit(models)
        self.dofit_chd(models)
    ##
    def test_oscillating2(self):
        """ corr2 -- oscillating part (1x2 matrix)"""
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a=("a","ao"), b=("a","ao"), 
                               dE=("logdE","logdEo")),
                   self.mkcorr(a=("a","ao"), b=("logb","bo"), 
                              dE=("logdE","logdEo")) ]
        self.dofit(models)
        self.dofit_chd(models)
    ##
    def test_marginalization2(self):
        """ corr2 -- marginalization (1x2 matrix fit w. osc.)"""
        if  PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a=("a","ao"), b=("a","ao"), 
                               dE=("logdE","logdEo"), s=(1,-1)),
                   self.mkcorr(a=("a","ao"), b=("logb","bo"), 
                              dE=("logdE","logdEo"), s=(1.,-1.)) ]
        self.dofit(models, nterm=(2,2), ratio=False)
        self.dofit_chd(models, nterm=(2,2), ratio=False)
        # N.B. setting ratio=True causes failures every 150 runs or so
    ##
    def test_oscillating3(self):
        """ corr2 -- oscillating part (only) """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a=(None,"a"), b=(None,"a"), 
                   dE=(None,"logdE"), s=(0,-1.)) ]
        self.dofit(models)
        self.dofit_chd(models)
    ##
    def test_fastfit_oscillating(self):
        """ corr2 -- fastfit(oscillating (only)) """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        self.dofastfit(self.mkcorr(a=(None,"a"), b=(None,"a"), 
                    dE=(None,"logdE"), s=(0,-1.)), osc=True)
    ##        
    def test_s1(self):
        """ corr2 -- s parameter #1"""
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a=("a","ao"), b=("a","ao"), 
                   dE=("logdE","logdEo"), s=(-1,1)) ]
        self.dofit(models)
        self.dofit_chd(models)
    ##
    def test_s2(self):
        """ corr2 -- s parameter #2"""
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ self.mkcorr(a=("a","ao"), b=("a","ao"), 
                   dE=("logdE","logdEo"), s=(1,1)) ]
        self.dofit(models)
        self.dofit_chd(models)
   ##
##
        
class test_corr3(unittest.TestCase, FitTests, ArrayTests):
    def setUp(self):
        ## prior ##
        self.prior = gv.BufferDict()
        nt = NTERM
        self.prior['a'] = gv.gvar(nt*["0.50(1)"])
        self.prior['ao'] = gv.gvar(nt*["0.250(5)"])
        self.prior['logb'] = gv.log(gv.gvar(nt*["0.60(1)"]))
        self.prior['bo'] = gv.gvar(nt*["0.30(1)"])
        self.prior['logdEa'] = gv.log(gv.gvar(nt*["0.50(1)"]))
        self.prior['logdEao'] = gv.log(gv.gvar(nt*["0.60(1)"]))
        self.prior['logdEb'] = gv.log(gv.gvar(nt*["0.45(1)"]))
        self.prior['logdEbo'] = gv.log(gv.gvar(nt*["0.65(1)"]))
        self.prior['Vnn'] = gv.gvar(nt*[nt*["2.00(1)"]])
        self.prior['Vno'] = gv.gvar(nt*[nt*["1.00(1)"]])
        self.prior['Von'] = gv.gvar(nt*[nt*["1.00(1)"]])
        self.prior['Voo'] = gv.gvar(nt*[nt*["2.00(1)"]])
        nsym = int(nt*(nt+1)/2)
        self.prior['Vnn_sym'] = gv.gvar(nsym*["2.00(1)"])
        self.prior['Voo_sym'] = gv.gvar(nsym*["2.00(1)"])
        ##
        ## actual parameters, time ranges, corr counter ##
        self.p = next(gv.raniter(self.prior))
        for x in ['b', 'dEa', 'dEao', 'dEb', 'dEbo']:
            self.p[x] = gv.exp(self.p['log' + x])
        self.T = 18.
        self.tdata = np.arange(self.T)
        self.tfit = self.tdata[1:]
        self.ncorr = 0
        ##
        self.ran = gv.gvar(0,1)
    ##
    def tearDown(self):
        del self.prior
        del self.p
        del self.T
        del self.tdata
        del self.tfit
        del self.ncorr
    ##
    def getdoc(self):
        """ get __doc__ string for current function """
        frame = inspect.currentframe()
        caller_frame = inspect.getouterframes(frame)[1][0]
        caller_name = inspect.getframeinfo(caller_frame).function
        caller_func = eval("test_corr3."+caller_name)
        return caller_func.__doc__
    ##
    def dofit(self, models, data_models=None, nterm=None, ratio=True):
        fitter = CorrFitter(models=models, nterm=nterm, ratio=ratio)
        if data_models is None:
            data_models = models
        data = make_data(models=data_models, p=self.p)
        if PRINT_FITS:
            print("Data:\n", data,"\n")
        fit = fitter.lsqfit(data=data, prior=self.prior, debug=True,
            print_fit=PRINT_FITS, svdcut=SVDCUT)
        #
        if PRINT_FITS:
            print("corr2/3 exact parameter values:")
            for k in fit.p:
                print("%10s:"%str(k),self.p[k])
            print()
        self.assert_fitclose(fit.p,self.p)
        return fitter
    ##
    def dofit_chd(self, models, data_models=None, nterm=None, ratio=True):
        fitter = CorrFitter(models=models, nterm=nterm, ratio=ratio)
        if data_models is None:
            data_models = models
        data = make_data(models=data_models, p=self.p)
        if PRINT_FITS:
            print("Data:\n", data,"\n")
        fit = fitter.chained_lsqfit(data=data, prior=self.prior, debug=True,
            print_fit=PRINT_FITS, svdcut=SVDCUT)
        #
        if PRINT_FITS:
            print("corr2/3 exact parameter values:")
            for k in fit.p:
                print("%10s:"%str(k),self.p[k])
            print()
        self.assert_fitclose(fit.p,self.p)
        return fitter
    ##
    def mkcorr2(self,a, b, dE, tp=None, othertags=[], s=1.):
        ans = Corr2(datatag=self.ncorr, a=a, b=b, dE=dE, tdata=self.tdata,
            tfit=self.tfit, tp=tp, s=s, othertags=othertags)
        self.ncorr += 1
        return ans
    ##
    def mkcorr3(self, a, b, dEa, dEb, Vnn, Vno=None, Von=None, Voo=None, #):
                symmetric_V=False, transpose_V=False, tpa=None, tpb=None):
        ans = Corr3(datatag=self.ncorr, a=a, b=b, dEa=dEa, dEb=dEb, Vnn=Vnn,
                    Vno=Vno, Von=Von, Voo=Voo, symmetric_V=symmetric_V,
                    transpose_V=transpose_V, tpa=tpa, tpb=tpb, T=self.T,
                    tdata=self.tdata, tfit=self.tfit)
        self.ncorr += 1
        return ans
    ##
    def test_symmetric(self):
        """ corr3 -- symmetric V"""
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ 
            self.mkcorr2(a="a", b="a", dE="logdEa", tp=None),
            self.mkcorr3(a="a", b="a", dEa="logdEa", dEb="logdEa", 
                         symmetric_V=True, Vnn="Vnn_sym")
        ]
        self.dofit(models)
        self.dofit_chd(models)
    ##
    def test_nonsymmetric(self):
        """ corr3 -- non-symmetric V"""
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ 
            self.mkcorr2(a="a", b="a", dE="logdEa"),
            self.mkcorr2(a="logb", b="logb", dE="logdEb"),
            self.mkcorr3(a="a", b="logb", dEa="logdEa", dEb="logdEb", 
                         Vnn="Vnn")
        ]
        self.dofit(models)
        self.dofit_chd(models)
    ##
    def test_transpose(self):
        """ corr3 -- transpose V"""
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ 
            self.mkcorr2(a="a", b="a", dE="logdEa"),
            self.mkcorr2(a="logb", b="logb", dE="logdEb"),
            self.mkcorr3(a="a", b="logb", dEa="logdEa", dEb="logdEb", 
                         Vnn="Vnn"),
            self.mkcorr3(b="a", a="logb", dEb="logdEa", dEa="logdEb", 
                         Vnn="Vnn", transpose_V=True)
            
        ]
        self.dofit(models)
        self.dofit_chd(models)
    ##
    def test_symmetric_osc(self):
        """ corr3 -- symmetric V with osc"""
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ 
            self.mkcorr2(a=("a", "ao"), b=("a", "ao"), 
                         dE=("logdEa", "logdEao")),
            self.mkcorr3(a=("a", "ao"), b=("a", "ao"), 
                         dEa=("logdEa", "logdEao"), 
                         dEb=("logdEa", "logdEao"), 
                         symmetric_V=True, Vnn="Vnn_sym",
                         Von="Von", Voo="Voo_sym")
        ]
        self.dofit(models)
        self.dofit_chd(models)
    ##
    def test_nonsymmetric_osc(self):
        """ corr3 -- non-symmetric V with osc"""
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ 
            self.mkcorr2(a=("a", "ao"), b=("a", "ao"), 
                         dE=("logdEa", "logdEao")),
            self.mkcorr2(a=("logb", "bo"), b=("logb", "bo"), 
                         dE=("logdEb", "logdEbo")),
            self.mkcorr3(a=("a", "ao"), b=("logb", "bo"), 
                         dEa=("logdEa", "logdEao"), 
                         dEb=("logdEb", "logdEbo"), 
                         Vnn="Vnn", Von="Von", Vno="Vno", Voo="Voo")
        ]
        fitter = self.dofit(models)
        if DISPLAY_PLOTS:
            fitter.display_plots()
    ##
    def test_transpose_osc(self):
        """ corr3 -- transpose V with osc"""
        global NSIG
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ 
            self.mkcorr2(a=("a", "ao"), b=("a", "ao"), 
                         dE=("logdEa", "logdEao")),
            self.mkcorr2(a=("logb", "bo"), b=("logb", "bo"), 
                         dE=("logdEb", "logdEbo")),
            self.mkcorr3(a=("a", "ao"), b=("logb", "bo"), 
                         dEa=("logdEa", "logdEao"), 
                         dEb=("logdEb", "logdEbo"), 
                         Vnn="Vnn"),
            self.mkcorr3(b=("a", "ao"), a=("logb", "bo"), 
                         dEb=("logdEa", "logdEao"), 
                         dEa=("logdEb", "logdEbo"), 
                         Vnn="Vnn", Vno="Vno", Von="Von", Voo="Voo",
                         transpose_V=True)
            
        ]
        self.dofit(models)
        NSIG *= 2.
        self.dofit_chd(models)
        NSIG /= 2.
    ##
    def test_periodic(self):
        """ corr3 -- periodic correlators """
        global NSIG
        if PRINT_FITS:
            print("======== " + self.getdoc())
        tp = 3*self.T
        models = [ 
            self.mkcorr2(a=("a", "ao"), b=("a", "ao"), 
                         dE=("logdEa", "logdEao"), tp=tp),
            self.mkcorr2(a=("logb", "bo"), b=("logb", "bo"), 
                         dE=("logdEb", "logdEbo"), tp=tp),
            self.mkcorr3(a=("a", "ao"), b=("logb", "bo"), 
                         dEa=("logdEa", "logdEao"), 
                         dEb=("logdEb", "logdEbo"), 
                         Vnn="Vnn", Von="Von", Vno="Vno", Voo="Voo",
                         tpa=tp, tpb=tp)
        ]
        NSIG *= 2.
        self.dofit(models)
        self.dofit_chd(models)
        NSIG /= 2.
    ##
    def test_antiperiodic(self):
        """ corr3 -- anti-periodic correlators """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        tp = -3*self.T 
        models = [ 
            self.mkcorr2(a=("a", "ao"), b=("a", "ao"), 
                         dE=("logdEa", "logdEao"), tp=tp),
            self.mkcorr2(a=("logb", "bo"), b=("logb", "bo"), 
                         dE=("logdEb", "logdEbo"), tp=tp),
            self.mkcorr3(a=("a", "ao"), b=("logb", "bo"), 
                         dEa=("logdEa", "logdEao"), 
                         dEb=("logdEb", "logdEbo"), 
                         Vnn="Vnn", Von="Von", Vno="Vno", Voo="Voo",
                         tpa=tp, tpb=tp)
        ]
        self.dofit(models)
        self.dofit_chd(models)
    ##
    @unittest.skipIf(FAST,"skipping test_bootstrap for speed")
    def test_bootstrap(self):
        """ corr3 -- bootstrap """
        if PRINT_FITS:
            print("======== " + self.getdoc())
        tp = 2*self.T
        models = [ 
            self.mkcorr2(a=("a", "ao"), b=("a", "ao"), 
                         dE=("logdEa", "logdEao"), tp=tp),
            self.mkcorr2(a=("logb", "bo"), b=("logb", "bo"), 
                         dE=("logdEb", "logdEbo"), tp=tp),
            self.mkcorr3(a=("a", "ao"), b=("logb", "bo"), 
                         dEa=("logdEa", "logdEao"), 
                         dEb=("logdEb", "logdEbo"), 
                         Vnn="Vnn", Von="Von", Vno="Vno", Voo="Voo",
                         tpa=tp, tpb=tp)
        ]
        self.dofit(models)
        fitter = self.dofit(models)
        bsdata = ds.Dataset()
        for fit in fitter.bootstrap_iter(n=50):
            bsdata.append(fit.pmean)
        bsfit = ds.avg_data(bsdata,bstrap=True)
        self.assert_fitclose(bsfit, self.p)
        self.assert_fitsagree(fitter.fit.p, bsfit)
    ##
    def test_marginalization(self):
        """ corr3 -- marginalization """
        global PRINT_FITS, NSIG
        if PRINT_FITS:
            print("======== " + self.getdoc())
        models = [ 
            self.mkcorr2(a=("a", "ao"), b=("a", "ao"), 
                         dE=("logdEa", "logdEao")),
            self.mkcorr3(a=("a", "ao"), b=("a", "ao"), 
                         dEa=("logdEa", "logdEao"), 
                         dEb=("logdEa", "logdEao"), 
                         Vnn="Vnn_sym" ,symmetric_V=True)
        ]
        NSIG *= 2.
        self.dofit(models, nterm=(NTERM-1, NTERM-1), ratio=True)
        NSIG /= 2.


class test_corrbasis(unittest.TestCase, FitTests, ArrayTests):
    def setUp(self):
        self.u = np.array([[1., 0.5], [1., 2.]])
        self.E = np.array([1., 3.])

    def make_G(self, tdata, keyfmt, srcs):
        G = collections.OrderedDict()
        tdata = np.array(tdata)
        for i, s1 in enumerate(srcs):
            for j, s2 in enumerate(srcs):
                key = keyfmt.format(s1=s1, s2=s2)
                G[key] = 0
                for n in range(2):
                    G[key] += (
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
        self.assert_gvclose(gv.exp(prior['log(a.dE)']), dE)
        a0 = gv.gvar(nterm * [big])
        a0[0] = gv.gvar(one)
        a0[1] = gv.gvar(small)
        self.assert_gvclose(prior['a.0'], a0)
        a1 = gv.gvar(nterm * [big])
        a1[0] = gv.gvar(small)
        a1[1] = gv.gvar(one)
        self.assert_gvclose(prior['a.1'], a1)

        # default values
        ampl = '1.0(3)', '0.03(10)', '0.2(1.0)'
        dEfac = '1(1)'
        one, small, big = ampl

        # case 2 - omit state
        prior = basis.make_prior(nterm=nterm, keyfmt='a.{s1}', states=[0], eig_srcs=True)
        a0 = gv.gvar(nterm * [big])
        a0[0] = gv.gvar(one)
        self.assert_gvclose(prior['a.0'], a0)
        a1 = gv.gvar(nterm * [big])
        self.assert_gvclose(prior['a.1'], a1)

        # case 3 - swap states
        prior = basis.make_prior(nterm=nterm, keyfmt='a.{s1}', states=[1, 0], eig_srcs=True)
        dE = gv.gvar(nterm * ['1(1)']) * (self.E[1] - self.E[0])
        dE[0] += self.E[0] - dE[0].mean
        self.assert_gvclose(gv.exp(prior['log(a.dE)']), dE)
        a1 = gv.gvar(nterm * [big])
        a1[0] = gv.gvar(one)
        a1[1] = gv.gvar(small)
        self.assert_gvclose(prior['a.1'], a1)
        a0 = gv.gvar(nterm * [big])
        a0[0] = gv.gvar(small)
        a0[1] = gv.gvar(one)
        self.assert_gvclose(prior['a.0'], a0)

def unpack(p,k):
    """ unpacks parameter k in dictionary p """
    ans = []
    for ki in k:
        if ki is None:
            ans.append(None)
        elif len(ki) > 3 and ki[:3] == 'log':
            ans.append(gv.exp(p[ki]))
        else:
            ans.append(p[ki])
    return ans
##

def make_f(tp):
    if tp is None:
        def f(E, t):
            return gv.exp(-E*t)
        ##
    elif tp >= 0:
        def f(E, t, tp=tp):
            return gv.exp(-E*t) + gv.exp(-E*(tp-t))
        ##
    else:
        def f(E, t, tp=tp):
            return gv.exp(-E*t) - gv.exp(-E*(-tp-t))
        ##
    return f
##

def corr3(p, m):
    """ build a corr3 -- p=param, m=model """
    def make_prop(p, t, a, dE, s, tp):
        f = make_f(tp)
        a = unpack(p, a)
        dE = unpack(p, dE)
        s = (s[0], 0.0 if s[1] == 0.0 else s[1]*(-1)**t)
        prop = []
        for ai, dEi, si in zip(a, dE, s):
            if ai is None or dEi is None:
                prop.append(None)
                continue
            ans = []
            sumdE = 0.0
            for aij, dEij in zip(ai, dEi):
                sumdE += dEij
                ans.append(si * aij * f(sumdE, t))
            prop.append(ans)
        return prop
    ##
    t = np.array(m.tdata)
    aprop = make_prop(p=p, t=t, a=m.a, dE=m.dEa, s=m.sa, tp=m.tpa)
    bprop = make_prop(p=p, t=m.T-t, a=m.b, dE=m.dEb, s=m.sb, tp=m.tpb)
    ans = 0.0
    for i, (apropi, Vi) in enumerate(zip(aprop, m.V)):
        if apropi is None:
            continue
        for j, (bpropj, Vij) in enumerate(zip(bprop,Vi)):
            if bpropj is None or Vij is None:
                continue
            V = gv.exp(p[Vij]) if Vij[:3] == 'log' else p[Vij]
            if i == j and m.symmetric_V:
                na = len(apropi)
                nb = len(bpropj)
                assert na == nb
                iterV = iter(V)
                V = np.empty((na,nb), dtype=V.dtype)
                for k in range(na):
                    for l in range(k,nb):
                        V[k, l] = next(iterV)
                        if k != l:
                            V[l, k] = V[k, l]
            if m.transpose_V or (i > j and m.symmetric_V):
                V = V.T
            for ak, Vk in zip(apropi, V):
                V_b = 0.0
                for bl, Vkl in zip(bpropj, Vk):
                    V_b += Vkl*bl
                ans += ak * V_b
    return ans
##   
    
def corr2(p, m):
    """ build a corr2 -- p=param, m=model"""
    f = make_f(m.tp)
    t = np.array(m.tdata)
    a = unpack(p, m.a)
    b = unpack(p, m.b)
    dE = unpack(p, m.dE)
    s = (m.s[0], m.s[1]*(-1)**t)
    ans = 0.0
    for ai, bi, dEi, si in zip(a, b, dE, s):
        # sum over normal and oscillating piece
        if ai is None or bi is None or dEi is None or si is None:
            continue
        E = np.array([sum(dEi[:n+1]) for n in range(len(dEi))])
        ans += si*np.array([np.sum(ai * bi * f(E,tj)) for tj in t])
    return ans
##

def add_noise(data,frac):
    """ add noise to correlators in list corrlist; frac = rel. size """
    global_noise = gv.gvar(1,frac)
    ans = gv.BufferDict()
    for k in data:
        ## add: a) uncorr. noise (smear for zeros); b) corr. noise ##
        corr = data[k]
        dcorr = np.abs(corr*frac)
        dcorr[1:-1] = (dcorr[1:-1] + dcorr[2:] + dcorr[:-2])/3.
        dcorr = gv.gvar(np.zeros(dcorr.shape),dcorr)
        dcorr = next(gv.bootstrap_iter(dcorr))
        ans[k] = (corr + dcorr)*global_noise
        ##
    return ans
##

def make_data(models,p):
    data = gv.BufferDict()
    for m in models:
        k = m.datatag
        data[k] = corr2(p=p, m=m) if isinstance(m, Corr2) else corr3(p=p, m=m)
    data = add_noise(data=data, frac=0.00001)
    return data
##

if __name__ == '__main__':
    unittest.main()

"""
Design Notes:
=============

* Script can run 1000 times without a failure but the marginalization tests
  can fail now and then (1 in 200 runs?). The bootstrap tests also fails
  occasionally, as can some of the corr3 tests. Other tests run many 1000s
  of times without failure.
  
* Corr2 tests fail completely if one tries to fit with fewer than the
  correct number of exponentials, in either the normal or oscillating
  parts. This is because statistical errors added to the correlators are
  tiny, making every piece important. For example, fitting in
  test_oscillating1() fails if one uses 2 instead of 3 oscillating terms in
  the fit (having generated data with 3 terms): chi**2/dof is around 50 and
  fit results deviate from the correct results by as much as 30 sigma or
  more. This is not realistic but it is quite useful for testing.
  
* The last bullet also means that the marginalization tests, which use only
  a single term in the fit (to fit data that uses 3 terms), are quite
  non-trivial --- as one would hope.

* Priors on energies and amplitudes are very tight in order to avoid the
  usual pathologies (eg, amplitude goes to zero and the energy is
  unconstrained). The point is to get fits that are certain to work. The
  small errors make things more gaussian and so more likely consistent with
  the assumptions underlying the fitting.


"""