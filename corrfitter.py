""" 
This module contains tools that facilitate least-squares fits, as functions
of time ``t``, of simulation (or other statistical) data for 2-point and
3-point correlators of the form::
        
    Gab(t)    =  <b(t) a(0)>
    Gavb(t,T) =  <b(T) V(t) a(0)>
        
where ``T > t > 0``. Each correlator is modeled using |Corr2| for 2-point
correlators, or |Corr3| for 3-point correlators in terms of amplitudes for
each source ``a``, sink ``b``, and vertex ``V``, and the energies
associated with each intermediate state. The amplitudes and energies are
adjusted in the least-squares fit to reproduce the data; they are defined
in a shared prior (typically a dictionary).
        
An object of type |CorrFitter| describes a collection of correlators and is
used to fit multiple models to data simultaneously. Fitting multiple
correlators simultaneously is important if there are statistical
correlations between the correlators. Any number of correlators may be
described and fit by a single |CorrFitter| object.

Typical code has the structure ::

    from corrfitter import CorrFitter
        
    data = make_data('mcfile')          # user-supplied routine
    models = make_models()              # user-supplied routine
    N = 4                               # number of terms in fit functions
    prior = make_prior(N)               # user-supplied routine
    fitter = CorrFitter(models=models)
    fit = fitter.lsqfit(data=data, prior=prior)  # do the fit
    print_results(fit, prior, data)     # user-supplied routine

where ``make_data`` assembles the correlator data to be fit, 
``make_prior`` defines Bayesian priors for the fit parameters, 
``fitter.lsqfit(...)`` does the fit, and ``print_results`` 
writes out the rsults. Sample code for each of these routines is
given in the tutorial documentation for |CorrFitter|.
"""

# Created by G. Peter Lepage, Cornell University, on 2010-11-26.
# Copyright (c) 2010-2013 G. Peter Lepage.
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

import lsqfit
import gvar as _gvar
import numpy
import math
import collections
import copy
__version__ = '3.4.2'

class BaseModel(object):
    """ Base class for correlator models. 
    
   Derived classes must define methods ``fitfcn``, ``buildprior``, and 
   ``builddata``, all of which are described below. In addition they
   can have attributes:
   
   .. attribute:: datatag
       
       |CorrFitter| builds fit data for the correlator by extracting the
       data in an input |Dataset| labelled by string ``datatag``. This
       label is stored in the ``BaseModel`` and must be passed to its
       constructor.
        
    .. attribute:: _abscissa
        
        (Optional) Array of abscissa values used in plots of the data and
        fit corresponding to the model. Plots are not made for a model that
        doesn't specify this attribute.
    """
    def __init__(self, datatag):
        super(BaseModel, self).__init__()
        self.datatag = datatag  # label
        self._abscissa = None   # for plots
    
    def fitfcn(self, p, nterm=None):
        """ Compute fit function fit parameters ``p`` using ``nterm`` terms. "
            
        :param p: Dictionary of parameter values.
        :type p: dictionary
        :param nterm: Restricts the number of non-oscillating terms in the
            fit function to ``nterm[0]`` and oscillating terms to
            ``nterm[1]``. Setting either (or both) to ``None`` implies that
            all terms in the prior are used.
        :type nterm: tuple of ``None`` or integers
        """
        raise NotImplementedError("fitfcn not defined")
    
    def builddata(self, data):
        """ Construct fit data. 
            
        Format of output must be same as format for fitfcn output.
        
        :param data: Dataset containing correlator data 
            (see ``gvar.dataset``).
        :type data: dictionary
        """
        raise NotImplementedError("builddata not defined")
    
    def buildprior(self, prior, nterm=None):
        """ Extract fit prior from ``prior``; resizing as needed. 
            
        If ``nterm`` is not ``None``, the sizes of the priors may need
        adjusting so that they correspond to the values specified in 
        ``nterm`` (for normal and oscillating pieces).
        
        :param prior: Dictionary containing *a priori* estimates of the 
            fit parameters.
        :type prior: dictionary
        :param nterm: Restricts the number of non-oscillating terms in the
            fit function to ``nterm[0]`` and oscillating terms to
            ``nterm[1]``. Setting either (or both) to ``None`` implies that
            all terms in the prior are used.
        :type nterm: tuple of ``None`` or integers        
        """
        raise NotImplementedError("buildprior not defined")

    @staticmethod
    def buildp0(p0, prior):
        if p0 is None or isinstance(p0, str):
            return p0
        new_p0 = _gvar.BufferDict()
        for i in prior:
            if i in p0:
                new_p0[i] = p0[i]
            else:
                raise ValueError("p0 missing key = " + str(i))
        return new_p0
    
    @staticmethod
    def _param(p, default=None):
        """ Parse fit-parameter label --- utility function. """
        if not isinstance(p, tuple):
            p = (p, default)
        ans = [BaseModel._paramkey(i) for i in p]
        return tuple(ans)

    _paramkey = staticmethod(lsqfit.transform_p.paramkey)
    _priorkey = staticmethod(lsqfit.transform_p.priorkey)

    @staticmethod
    def _transform_prior(prior):
        return lsqfit.transform_p(prior.keys()).transform(prior)

    @staticmethod
    def _dE(dE, logdE):
        """ Parse fit-parameter label for E differences --- utility fcn. """
        if dE is None:
            assert logdE is not None, "must specify dE"
            return BaseModel._param(logdE, None)
        else:
            return BaseModel._param(dE, None)
                  

class Corr2(BaseModel):
    """ Two-point correlators ``Gab(t) = <b(t) a(0)>``.
        
    |Corr2| models the ``t`` dependence of a 2-point correlator ``Gab(t)``
    using ::
        
        Gab(t) = sn * sum_i an[i]*bn[i] * fn(En[i], t)
               + so * sum_i ao[i]*bo[i] * fo(Eo[i], t)
        
    where ``sn`` and ``so`` are typically ``-1``, ``0``, or ``1`` and ::
        
        fn(E, t) =  exp(-E*t) + exp(-E*(tp-t)) # tp>0 -- periodic
               or   exp(-E*t) - exp(-E*(-tp-t))# tp<0 -- anti-periodic
               or   exp(-E*t)                  # if tp is None (nonperiodic)
        
        fo(E, t) = (-1)**t * fn(E, t)
        
    The fit parameters for the non-oscillating piece of ``Gab`` (first term)
    are ``an[i]``, ``bn[i]``, and ``dEn[i]`` where::
        
        dEn[0] = En[0] > 0
        dEn[i] = En[i]-En[i-1] > 0     (for i>0)
        
    and therefore ``En[i] = sum_j=0..i dEn[j]``. The fit parameters for
    the oscillating pied are defined analogously: ``ao[i]``, ``bo[i]``,
    and ``dEo[i]``.
        
    The fit parameters are specified by the keys corresponding to these
    parameters in a dictionary of priors supplied by |CorrFitter|. The keys
    are strings and are also used to access fit results. Any key that
    begins with "log" is assumed to refer to the logarithm of the parameter
    in question (that is, the exponential of the fit-parameter is used in
    the formula for ``Gab(t)``.) This is useful for forcing ``an``, ``bn``
    and/or ``dE`` to be positive.
        
    When ``tp is not None`` and positive, the correlator is assumed to be
    symmetrical about ``tp/2``, with ``Gab(t)=Gab(tp-t)``. Data from
    ``t>tp/2`` is averaged with the corresponding data from ``t<tp/2``
    before fitting. When ``tp`` is negative, the correlator is assumed to
    be anti-symetrical about ``-tp/2``.
        
    :param datatag: Key used to access correlator data in the input data 
        dictionary (see |CorrFitter|). ``data[self.datatag]`` is (1-d) 
        array containing the correlator values (|GVar|\s) if ``data`` is the 
        input data.
    :type datatag: string
    :param a: Key identifying the fit parameters for the source amplitudes
        ``an`` in the dictionary of priors provided by |CorrFitter|; or a
        two-tuple of keys for the source amplitudes ``(an, ao)``. The
        corresponding values in the dictionary of priors are (1-d) arrays
        of prior values with one term for each ``an[i]`` or ``ao[i]``.
        Replacing either key by ``None`` causes the corresponding term to
        be dropped from the fit function. These keys are used to label the
        corresponding parameter arrays in the fit results as well as in the
        prior. 
    :type a: string, or two-tuple of strings and/or ``None``
    :param b: Same as ``self.a`` but for the sinks ``(bn, bo)`` instead of
        the sources ``(an, ao)``.
    :type b: string, or two-tuple of strings and/or ``None``
    :param dE: Key identifying the fit parameters for the energy 
        differences ``dEn`` in the dictionary of priors provided by
        |CorrFitter|; or a two-tuple of keys for the energy differences
        ``(dEn, dEo)``. The corresponding values in the dictionary of priors
        are (1-d) arrays of prior values with one term for each ``dEn[i]``
        or ``dEo[i]``. Replacing either key by ``None`` causes the
        corresponding term to be dropped from the fit function. These keys
        are used to label the corresponding parameter arrays in the fit
        results as well as in the prior. 
    :type dE: string, or two-tuple of strings and/or ``None``
    :param s: Overall factor ``sn``, or two-tuple of overall factors 
        ``(sn, so)``. 
    :type s: number or two-tuple of numbers
    :param tdata: The ``t``\s corresponding to data entries in the input
        data. Note that ``len(self.tdata) == len(data[self.datatag])`` is
        required if ``data`` is the input data dictionary.
    :type tdata: list of integers
    :param tfit: List of ``t``\s to use in the fit. Only data with these
        ``t``\s (all of which should be in ``tdata``) is used in the fit.
    :type tfit: list of integers
    :param tp: If not ``None`` and positive, the correlator is assumed to 
        be periodic with ``Gab(t)=Gab(tp-t)``. If negative, the correlator
        is assumed to be anti-periodic with ``Gab(t)=-Gab(-tp-t)``. Setting
        ``tp=None`` implies that the correlator is not periodic, but rather
        continues to fall exponentially as ``t`` is increased indefinitely.
    :type tp: integer or ``None``
    :param othertags: List of additional data tags for data to be
        averaged with the ``self.datatag`` data before fitting.
    :type othertags: sequence of strings
    """
    def __init__(self, datatag, tdata, tfit, a, b, dE=None, logdE=None,    #):
                    s=1.0, tp=None, othertags=None):
        super(Corr2, self).__init__(datatag)
        self.a = self._param(a)
        self.b = self._param(b)
        self.dE = self._dE(dE, logdE)
        self.tdata = list(tdata)
        self.tp = tp
        self.s = self._param(s, -1.)
        # verify and compress tfit 
        ntfit = []
        for t in tfit:
            if tp is None:
                assert t in tdata, ("tfit incompatible with tdata: "
                                  +str(tfit)+" "+str(tdata))
                ntfit.append(t)
            else:
                t1, t2 = sorted([t, abs(tp)-t])
                if t1 in ntfit or t2 in ntfit:
                    continue
                assert (t >= 0 and t < abs(tp)), "illegal t in tfit: "+str(t)
                if t1 in tdata:
                    ntfit.append(t1)
                elif t2 in tdata:
                    ntfit.append(t2)
                else:
                    raise ValueError("tfit incompatible with tdata: "
                                      +str(tfit)+" "+str(tdata))
        
        self.tfit = numpy.array(ntfit)
        self.othertags = othertags
        self._abscissa = self.tfit
    
    def buildprior(self, prior, nterm):
        """ Create fit prior by extracting relevant pieces of ``prior``. 

        Priors for the fit parameters, as specificied by ``self.a`` etc., 
        are copied from ``prior`` into a new dictionary for use by the
        fitter. If a key ``"XX"`` cannot be found in ``prior``, the
        ``buildprior`` looks for one of ``"logXX"``, ``"log(XX)"``, 
        ``"sqrtXX"``, or ``"sqrt(XX)"`` and includes the corresponding
        prior instead.

        The number of terms kept in each part of the fit can be 
        specified using ``nterm = (n, no)`` where ``n`` is the 
        number of non-oscillating terms and ``no`` is the number 
        of oscillating terms. Setting ``nterm = None`` keeps 
        all terms.
        """

        newprior = _gvar.BufferDict()
        for ai, bi, dEi, ntermi in zip(self.a, self.b, self.dE, nterm):
            for x in [ai, bi, dEi]:
                if x is None:
                    continue
                x = self._priorkey(prior, x)
                newprior[x] = prior[x][None:ntermi]
        return newprior        
    
    def builddata(self, data):
        """ Assemble fit data from dictionary ``data``. 
            
        Extracts parts of array ``data[self.datatag]`` that are needed for
        the fit, as specified by ``self.tp`` and ``self.tfit``. The entries
        in the (1-D) array ``data[self.datatag]`` are assumed to be
        |GVar|\s and correspond to the ``t``s in ``self.tdata``.
        """
        tags = [self.datatag]
        if self.othertags is not None:
            tags.extend(self.othertags)
        tdata = self.tdata
        tp = self.tp
        if tp is not None:
            pfac = math.copysign(1,tp)
            tp = abs(tp)
        ans = []
        for tag in tags:
            odata = data[tag]
            ndata = [] 
            for t in self.tfit:
                idt = tdata.index(t)
                if tp is None or tp-t not in tdata or t == tp-t:
                    ndata.append(odata[idt])
                else:
                    ndata.append(lsqfit.wavg([odata[idt],
                                              pfac*odata[tdata.index(tp-t)]]))
            ans.append(ndata)
        fdata = numpy.array(ans[0]) if len(ans) == 1 else lsqfit.wavg(ans)
        return fdata 
    
    def fitfcn(self, p, nterm=None, t=None):
        """ Return fit function for parameters ``p``. """
        if t is None:
            t = self.tfit
        if self.tp is None:
            tp_t = None
        elif self.tp >= 0:
            tp_t = self.tp - t
            pfac = 1
        else:
            tp_t = -self.tp - t
            pfac = -1
        if nterm is None:
            nterm = (None, None)
        ofac = (None if self.s[0] == 0.0 else self.s[0],
                (None if self.s[1] == 0.0 else self.s[1]*(-1)**t))
        ans = 0.0
        for _ai, _bi, _dEi, ofaci, ntermi in zip(self.a, self.b, 
                                              self.dE, ofac, nterm):
            if _ai is None or _bi is None or _dEi is None or ofaci is None:
                continue
            if ntermi is not None:
                if ntermi == 0:
                    continue
                ai = p[_ai][:ntermi]
                bi = p[_bi][:ntermi]
                dEi = p[_dEi][:ntermi]
            else:
                ai = p[_ai]
                bi = p[_bi]
                dEi = p[_dEi]
            if tp_t is None:
                exp_t = _gvar.exp(-t)
                for aij, bij, sumdE in zip(ai, bi, numpy.cumsum(dEi)):
                    ans += ofaci * aij * bij * exp_t ** sumdE
            else:
                exp_t = _gvar.exp(-t)
                exp_tp_t = _gvar.exp(-tp_t)
                for aij, bij, sumdE in zip(ai, bi, numpy.cumsum(dEi)):
                    ans += ofaci * aij * bij * (exp_t ** sumdE + pfac * exp_tp_t ** sumdE)
        return ans    


class Corr3(BaseModel):
    """ Three-point correlators ``Gavb(t, T) = <b(T) V(t) a(0)>``.
        
    |Corr3| models the ``t`` dependence of a 3-point correlator
    ``Gavb(t, T)`` using ::
        
        Gavb(t, T) = 
         sum_i,j san*an[i]*fn(Ean[i],t)*Vnn[i,j]*sbn*bn[j]*fn(Ebn[j],T-t)
        +sum_i,j san*an[i]*fn(Ean[i],t)*Vno[i,j]*sbo*bo[j]*fo(Ebo[j],T-t)
        +sum_i,j sao*ao[i]*fo(Eao[i],t)*Von[i,j]*sbn*bn[j]*fn(Ebn[j],T-t)
        +sum_i,j sao*ao[i]*fo(Eao[i],t)*Voo[i,j]*sbo*bo[j]*fo(Ebo[j],T-t)
       
    where ::
        
        fn(E, t) =  exp(-E*t) + exp(-E*(tp-t)) # tp>0 -- periodic
               or   exp(-E*t) - exp(-E*(-tp-t))# tp<0 -- anti-periodic
               or   exp(-E*t)                  # if tp is None (nonperiodic)
        
        fo(E, t) = (-1)**t * fn(E, t)
        
    The fit parameters for the non-oscillating piece of ``Gavb`` (first term)
    are ``Vnn[i,j]``, ``an[i]``, ``bn[j]``, ``dEan[i]`` and ``dEbn[j]`` where,
    for example::
        
        dEan[0] = Ean[0] > 0
        dEan[i] = Ean[i]-Ean[i-1] > 0     (for i>0)
        
    and therefore ``Ean[i] = sum_j=0..i dEan[j]``. The parameters for the
    other terms are similarly defined.
        
    :param datatag: Tag used to label correlator in the input |Dataset|.
    :type datatag: string
    :param a: Key identifying the fit parameters for the source amplitudes
        ``an``, for ``a->V``, in the dictionary of priors provided by
        |CorrFitter|; or a two-tuple of keys for the source amplitudes
        ``(an, ao)``. The corresponding values in the dictionary of priors
        are (1-d) arrays of prior values with one term for each ``an[i]``
        or ``ao[i]``. Replacing either key by ``None`` causes the
        corresponding term to be dropped from the fit function. These keys
        are used to label the corresponding parameter arrays in the fit
        results as well as in the prior. 
    :type a: string, or two-tuple of strings or ``None``
    :param b: Same as ``self.a`` except for sink amplitudes ``(bn, bo)`` 
        for ``V->b`` rather than for ``(an, ao)``.
    :type b: string, or two-tuple of strings or ``None``
    :param dEa: Fit-parameter label for ``a->V`` intermediate-state energy 
        differences ``dEan``, or two-tuple of labels for the differences
        ``(dEan,dEao)``. Each label represents an array of energy differences.
        Replacing either label by ``None`` causes the corresponding term in
        the correlator function to be dropped. These keys
        are used to label the corresponding parameter arrays in the fit
        results as well as in the prior. 
    :type dEa: string, or two-tuple of strings or ``None``
    :param dEb: Fit-parameter label for ``V->b`` intermediate-state energy 
        differences ``dEbn``, or two-tuple of labels for the differences
        ``(dEbn,dEbo)``. Each label represents an array of energy differences.
        Replacing either label by ``None`` causes the corresponding term in
        the correlator function to be dropped. These keys
        are used to label the corresponding parameter arrays in the fit
        results as well as in the prior. 
    :type dEb: string, or two-tuple of strings or ``None``
    :param sa: Overall factor ``san`` for the non-oscillating ``a->V`` terms 
        in the correlator, or two-tuple containing the overall factors
        ``(san,sao)`` for the non-oscillating and oscillating terms.
    :type sa: number, or two-tuple of numbers
    :param sb: Overall factor ``sbn`` for the non-oscillating ``V->b`` terms 
        in the correlator, or two-tuple containing the overall factors
        ``(sbn,sbo)`` for the non-oscillating and oscillating terms.
    :type sb: number, or two-tuple of numbers
    :param Vnn: Fit-parameter label for the matrix of current matrix 
        elements ``Vnn[i,j]`` connecting non-oscillating states. Labels that
        begin with "log" indicate that the corresponding matrix elements are
        replaced by their exponentials; these parameters are logarithms of the
        corresponding matrix elements, which must then be positive.
    :type Vnn: string or ``None``
    :param Vno: Fit-parameter label for the matrix of current matrix 
        elements ``Vno[i,j]`` connecting non-oscillating to oscillating
        states. Labels that begin with "log" indicate that the corresponding
        matrix elements are replaced by their exponentials; these parameters
        are logarithms of the corresponding matrix elements, which must then
        be positive.
    :type Vno: string or ``None``
    :param Von: Fit-parameter label for the matrix of current matrix 
        elements ``Von[i,j]`` connecting oscillating to non-oscillating 
        states. Labels that begin with "log" indicate that the corresponding
        matrix elements are replaced by their exponentials; these parameters
        are logarithms of the corresponding matrix elements, which must then
        be positive.
    :type Von: string or ``None``
    :param Voo: Fit-parameter label for the matrix of current matrix 
        elements ``Voo[i,j]`` connecting oscillating states. Labels that begin
        with "log" indicate that the corresponding matrix elements are
        replaced by their exponentials; these parameters are logarithms of the
        corresponding matrix elements, which must then be positive.
    :type Voo: string or ``None``
    :param transpose_V: If ``True``, the transpose ``V[j,i]`` is used in
        place of ``V[i,j]`` for each current matrix element in the fit 
        function. This is useful for doing simultaneous fits to 
        ``a->V->b`` and ``b->V->a``, where the current matrix elements
        for one are the transposes of those for the other. Default value 
        is ``False``.
    :type transpose_V: boolean
    :param symmetric_V: If ``True``, the fit function for ``a->V->b`` is 
        unchanged (symmetrical) under the the interchange of ``a`` and
        ``b``. Then ``Vnn`` and ``Voo`` are square, symmetric matrices
        with ``V[i,j]=V[j,i]`` and their priors are one-dimensional arrays
        containing only elements ``V[i,j]`` with ``j>=i`` in the following
        layout::
        
            [V[0,0],V[0,1],V[0,2]...V[0,N],
                    V[1,1],V[1,2]...V[1,N],
                           V[2,2]...V[2,N],
                                 .
                                  .
                                   .
                                    V[N,N]]
                                    
        Furthermore the matrix specified for ``Von`` is transposed before
        being used by the fitter; normally the matrix specified for ``Von``
        is the same as the matrix specified for ``Vno`` when the amplitude
        is symmetrical. Default value is ``False``.
    :type symmetric_V: boolean
    :param tdata: The ``t``\s corresponding to data entries in the input
        |Dataset|.
    :type tdata: list of integers
    :param tfit: List of ``t``\s to use in the fit. Only data with these
        ``t``\s (all of which should be in ``tdata``) is used in the fit.
    :type tfit: list of integers
    :param tpa: If not ``None`` and positive, the ``a->V`` correlator is 
        assumed to be periodic with period ``tpa``. If negative, the
        correlator is anti-periodic with period ``-tpa``. Setting
        ``tpa=None`` implies that the correlators are not periodic.
    :type tpa: integer or ``None``
    :param tpb: If not ``None`` and positive, the ``V->b`` correlator is 
        assumed to be periodic with period ``tpb``. If negative, the
        correlator is periodic with period ``-tpb``. Setting ``tpb=None``
        implies that the correlators are not periodic.
    :type tpb: integer or ``None``
    """
    def __init__(self, datatag, T, tdata, tfit,          #):
                 Vnn, a, b, dEa=None, dEb=None, logdEa=None, logdEb=None, 
                 sa=1., sb=1.,
                 Vno=None, Von=None, Voo=None, transpose_V=False,
                 symmetric_V=False, tpa=None, tpb=None,
                 othertags=None):
        super(Corr3, self).__init__(datatag)
        self.a = self._param(a)
        self.dEa = self._dE(dEa, logdEa)
        self.sa = self._param(sa, -1.)
        self.b = self._param(b)
        self.dEb = self._dE(dEb, logdEb)
        self.sb = self._param(sb, -1.)
        self.V = [[Vnn, Vno], [Von, Voo]]
        self.transpose_V = transpose_V
        self.symmetric_V = symmetric_V
        self.T = T
        self.tdata = list(tdata)
        self.tpa = tpa
        self.tpb = tpb
        # validate tfit 
        ntfit = []
        for t in tfit:
            if t >= 0 and t <= T:
                ntfit.append(t)
        self.tfit = numpy.array(ntfit)
        self._abscissa = self.tfit
        self.othertags = othertags
    
    def buildprior(self, prior, nterm):
        """ Create fit prior by extracting relevant pieces of ``prior``. 

        Priors for the fit parameters, as specificied by ``self.a`` etc., 
        are copied from ``prior`` into a new dictionary for use by the
        fitter. If a key ``"XX"`` cannot be found in ``prior``, the
        ``buildprior`` looks for one of ``"logXX"``, ``"log(XX)"``, 
        ``"sqrtXX"``, or ``"sqrt(XX)"`` and includes the corresponding
        prior instead.

        The number of terms kept in each part of the fit can be 
        specified using ``nterm = (n, no)`` where ``n`` is the 
        number of non-oscillating terms and ``no`` is the number 
        of oscillating terms. Setting ``nterm = None`` keeps 
        all terms.
        """
        def resize_sym(Vii, nterm):
            N = int(numpy.round((((8*len(Vii)+1)**0.5 - 1.)/2.)))
            ans = []
            iterV = iter(Vii)
            for i in range(N):
                for j in range(i, N):
                    v = next(iterV)
                    if ((nterm[0] is None or i < nterm[0])
                        and (nterm[1] is None or j < nterm[1])):
                        ans.append(v)
            return numpy.array(ans)
        ans = _gvar.BufferDict()
        for x in [self.a, self.dEa, self.b, self.dEb]:
            for xi, ntermi in zip(x, nterm):
                if xi is not None:
                    xi = self._priorkey(prior, xi)
                    ans[xi] = prior[xi][None:ntermi]
        for i in range(2):
            for j in range(2):
                vij = self.V[i][j]
                if vij is None:
                    continue
                vij = self._priorkey(prior, vij)
                if i == j and self.symmetric_V:
                    ans[vij] = resize_sym(prior[vij], nterm)
                else:
                    ans[vij] = prior[vij][slice(None, nterm[i]), 
                                          slice(None, nterm[j])]
        return ans         
    
    def builddata(self, data):
        """ Assemble fit data from dictionary ``data``. 
            
        Extracts parts of array ``data[self.datatag]`` that are needed for
        the fit, as specified by ``self.tfit``. The entries in the (1-D)
        array ``data[self.datatag]`` are assumed to be |GVar|\s and
        correspond to the ``t``s in ``self.tdata``.
        """
        tags = [self.datatag]
        if self.othertags is not None:
            tags.extend(self.othertags)
        ans = []
        for tag in tags:
            odata = data[tag]
            tdata = self.tdata
            ndata = []
            for t in self.tfit:
                idt = tdata.index(t)
                ndata.append(odata[idt])
            ans.append(ndata)
        return numpy.array(ans[0]) if len(ans) == 1 else lsqfit.wavg(ans)
    
    def fitfcn(self, p, nterm=None, t=None):
        """ Return fit function for parameters ``p``. """
        # setup
        if t is None:
            t = self.tfit
        ta = t
        tb = self.T - t
        if self.tpa is None:
            tp_ta = None
        elif self.tpa >= 0:
            tp_ta = self.tpa - ta
            pafac = 1
        else:
            tp_ta = -self.tpa - ta
            pafac = -1
        #
        if self.tpb is None:
            tp_tb = None
        elif self.tpb >= 0:
            tp_tb = self.tpb - tb
            pbfac = 1
        else:
            tp_tb = -self.tpb - tb
            pbfac = -1
        if nterm is None:
            nterm = (None, None)

        # initial and final propagators  
        aprop = []  # aprop[i][j] i= n or o; j=excitation level
        ofac = (self.sa[0], (0.0 if self.sa[1] == 0.0 else self.sa[1]*(-1)**ta))
        for _ai, _dEai, ofaci, ntermai in zip(self.a, self.dEa, ofac, nterm):
            if _ai is None:
                aprop.append(None)
                continue
            ans = []
            if ntermai is None:
                ai =  p[_ai] 
                dEai = p[_dEai] 
            else:
                if ntermai <= 0:
                    aprop.append(None)
                    continue
                ai =  p[_ai][:ntermai]
                dEai = p[_dEai][:ntermai] 
            if tp_ta is None:
                exp_ta = _gvar.exp(-ta)
                ans = [
                    ofaci * aij * exp_ta ** sumdE 
                    for aij, sumdE in zip(ai, numpy.cumsum(dEai))
                    ]
            else:
                exp_ta = _gvar.exp(-ta)
                exp_tp_ta = _gvar.exp(-tp_ta)
                ans = [
                    ofaci * aij * (exp_ta ** sumdE + pafac * exp_tp_ta ** sumdE) 
                    for aij, sumdE in zip(ai, numpy.cumsum(dEai))
                    ]
            aprop.append(ans)
        bprop = []
        ofac = (self.sb[0], (0.0 if self.sb[1] == 0.0 else self.sb[1]*(-1)**tb))
        for _bi, _dEbi, ofaci, ntermbi in zip(self.b, self.dEb, ofac, nterm):
            if _bi is None:
                bprop.append(None)
                continue
            ans = []
            if ntermbi is None:
                bi = p[_bi] 
                dEbi = p[_dEbi]
            else:
                if ntermbi <= 0:
                    bprop.append(None)
                    continue
                bi = p[_bi][:ntermbi] 
                dEbi = p[_dEbi][:ntermbi] 
            if tp_tb is None:
                exp_tb = _gvar.exp(-tb)
                ans = [
                    ofaci * bij * exp_tb ** sumdE 
                    for bij, sumdE in zip(bi, numpy.cumsum(dEbi))
                    ]
            else:
                exp_tb = _gvar.exp(-tb)
                exp_tp_tb = _gvar.exp(-tp_tb)
                ans = [
                    ofaci * bij * (exp_tb ** sumdE + pbfac * exp_tp_tb ** sumdE) 
                    for bij, sumdE in zip(bi, numpy.cumsum(dEbi))
                    ]
            bprop.append(ans)
        
        # combine propagators with vertices 
        ans = 0.0
        for i, (apropi, Vi) in enumerate(zip(aprop, self.V)):
            if apropi is None:
                continue
            for j, (bpropj, Vij) in enumerate(zip(bprop, Vi)):
                if bpropj is None or Vij is None:
                    continue
                V = p[Vij]
                if i == j and self.symmetric_V:
                    # unpack symmetric matrix V 
                    na = len(apropi)
                    nb = len(bpropj)
                    assert na == nb, \
                        "Vnn and Voo must be square matrices if symmetric"
                    iterV = iter(V)
                    V = numpy.empty((na, nb), dtype=V.dtype)
                    for k in range(na):
                        for l in range(k, nb):
                            V[k, l] = next(iterV)
                            if k != l:
                                V[l, k] = V[k, l]
                    
                if self.transpose_V or (i>j and self.symmetric_V):
                    V = V.T
                for ak, Vk in zip(apropi, V):
                    acc = 0.0
                    for bl, Vkl in zip(bpropj, Vk):
                        acc += Vkl*bl
                    ans += ak*acc
        return ans                
            
class CorrFitter(object):
    """ Nonlinear least-squares fitter for a collection of correlators. 
        
    :param models: Correlator models used to fit statistical input data.
    :type models: list of correlator models or a single model
    :param svdcut: If ``svdcut`` is positive, eigenvalues ``ev[i]`` of the 
        (rescaled) data covariance matrix that are smaller than
        ``svdcut*max(ev)`` are replaced by ``svdcut*max(ev)`` in the
        covariance matrix. If ``svdcut`` is negative, eigenvalues less than
        ``|svdcut|*max(ev)`` are set to zero in the covariance matrix. The
        covariance matrix is left unchanged if ``svdcut`` is set equal to
        ``None`` (default). If ``svdcut`` is a 2-tuple, *svd* cuts are applied
        to both the correlator data (``svdcut[0]``) and to the prior
        (``svdcut[1]``).
    :type svdcut: number or ``None`` or 2-tuple
    :param svdnum: At most ``svdnum`` eigenmodes are retained in the 
        (rescaled) data covariance matrix; the modes with the smallest
        eigenvalues are discarded. ``svdnum`` is ignored if it is set to
        ``None``. If ``svdnum`` is a 2-tuple, *svd* cuts are applied to both
        the correlator data (``svdnum[0]``) and to the prior (``svdnum[1]``).
    :type svdnum: integer or ``None`` or 2-tuple
    :param tol: Tolerance used in :func:`lsqfit.nonlinear_fit` for the 
        least-squares fits (default=1e-10).
    :type tol: positive number less than 1
    :param maxit: Maximum number of iterations to use in least-squares fit 
        (default=500).
    :type maxit: integer
    :param nterm: Number of terms fit in the non-oscillating parts of fit 
        functions; or two-tuple of numbers indicating how many terms to fit
        for each of the non-oscillating and oscillating pieces in fits. If set
        to ``None``, the number is specified by the number of parameters in
        the prior.
    :type nterm: number or ``None``; or two-tuple of numbers or ``None``
    :param ratio: If ``True`` (the default), use ratio corrections for fit 
        data when the prior specifies more terms than are used in the fit. If
        ``False``, use difference corrections (see implementation notes,
        above).
    :type ratio: boolean
    """
    def __init__(self, models, svdcut=(1e-15, 1e-15), svdnum=None, tol=1e-10,
                maxit=500, nterm=None, ratio=True): # mc=None, ratio=True):
        super(CorrFitter, self).__init__()
        self.models = [models] if isinstance(models, BaseModel) else models
        self.svdcut = svdcut
        self.svdnum = svdnum
        self.tol = tol
        self.maxit = maxit
        self.fit = None
        self.dset = None
        # self.mc = mc
        self.ratio = ratio
        self.nterm = nterm if isinstance(nterm, tuple) else (nterm, None)
    
    def buildfitfcn(self, priorkeys):
        " Create fit function, with support for log-normal,... priors. "

        @lsqfit.transform_p(priorkeys, pindex=0, pkey='p')
        def _fitfcn(
            p, nterm=None, default_nterm=self.nterm, models=self.models
            ):
            """ Composite fit function. 
                
            :param p: Fit parameters.
            :type p: dict-like
            :param nterm: Number of terms fit in the non-oscillating parts of fit 
                functions; or two-tuple of numbers indicating how many terms to
                fit for each of the non-oscillating and oscillating pieces in
                fits. If set to ``None``, the number is specified by the number of
                parameters in the prior.
            :type nterm: number or ``None``; or two-tuple of numbers or ``None``
            :returns: A dictionary containing the fit function results for 
                parameters ``p`` from each model, indexed using the models' 
                ``datatag``\s.
            """
            ans = _gvar.BufferDict()
            if nterm is None:
                nterm = default_nterm
            for m in models:
                ans[m.datatag] = m.fitfcn(p, nterm=nterm)
            return ans

        return _fitfcn

    def builddata(self, data, prior, nterm=None):
        """ Build fit data, corrected for marginalized terms. """
        fitdata = _gvar.BufferDict()
        for m in self.models:
            fitdata[m.datatag] = m.builddata(data)
        # remove marginal fit parameters 
        if nterm is None:
            nterm = self.nterm
        # use priors to remove marginal parameters 
        if nterm == (None, None):
            return fitdata
        else:
            fitfcn = self.buildfitfcn(prior.keys())
            ftrunc = fitfcn(prior, nterm=nterm)
            fall = fitfcn(prior, nterm=(None, None))
            for m in self.models:
                if not self.ratio:
                    diff = ftrunc[m.datatag] - fall[m.datatag]
                    fitdata[m.datatag] += diff 
                else:
                    ii = (_gvar.mean(fall[m.datatag]) != 0.0)
                    ratio = ftrunc[m.datatag][ii]/fall[m.datatag][ii]
                    fitdata[m.datatag][ii] *= ratio 
        return fitdata
    
    def buildprior(self, prior, nterm=None):
        """ Build correctly sized prior for fit. """
        tmp = _gvar.BufferDict()
        if nterm is None:
            nterm = self.nterm
        for m in self.models:
            tmp.update(m.buildprior(prior=prior, nterm=nterm))
        # restore order of keys --- same as prior 
        ans = _gvar.BufferDict()
        for k in prior:
            if k in tmp:
                ans[k] = tmp[k]
        return ans
    
    def lsqfit(
        self, data, prior, p0=None, print_fit=True, nterm=None,
        svdcut=None, svdnum=None, tol=None, maxit=None, 
        aux_param={},
        **args
        ):
        """ Compute least-squares fit of the correlator models to data.
            
        :param data: Input data. The ``datatag``\s from the 
            correlator models are used as data labels, with 
            ``data[datatag]`` being a 1-d array of |GVar|\s 
            corresponding to correlator values.
        :type data: dictionary
        :param prior: Bayesian prior for the fit parameters used in the 
            correlator models.
        :type prior: dictionary
        :param p0: A dictionary, indexed by parameter labels, containing 
            initial values for the parameters in the fit. Setting
            ``p0=None`` implies that initial values are extracted from the
            prior. Setting ``p0="filename"`` causes the fitter to look in
            the file with name ``"filename"`` for initial values and to
            write out best-fit parameter values after the fit (for the next
            call to ``self.lsqfit()``).
        :param print_fit: Print fit information to standard output if 
            ``True``; otherwise print nothing.
        :param aux_param: A dictionary containing priors for 
            additional (auxilary) fit parameters. Where ``aux_param``
            and ``prior`` share keys, the entry in ``aux_param`` is
            used.
        :type aux_param: dictionary

        The following parameters overwrite the values specified in the
        |CorrFitter| constructor when set to anything other than ``None``:
        ``nterm``, ``svdcut``, ``svdnum``, ``tol``, and ``maxit``. Any
        further keyword arguments are passed on to
        :func:`lsqfit.nonlinear_fit`, which does the fit.
        """
        # setup
        if svdcut is None:
            svdcut = self.svdcut
        if svdnum is None:
            svdnum = self.svdnum
        if maxit is None:
            maxit = self.maxit
        if tol is None:
            tol = self.tol
        if nterm is None:
            nterm = self.nterm
        self.last_prior = prior
        self.aux_param = aux_param

        # do the fit and print results
        data = self.builddata(data=data, prior=prior, nterm=nterm)
        prior = self.buildprior(prior, nterm=nterm)
        fitfcn = self.buildfitfcn(prior.keys())
        for k in self.aux_param:
            prior[k] = self.aux_param[k]
        self.fit = lsqfit.nonlinear_fit( #
            data=data, p0=p0, fcn=fitfcn, prior=prior, 
            svdcut=svdcut, svdnum=svdnum, reltol=tol, 
            abstol=tol, maxit=maxit, **args
            )
        if print_fit:
            print(self.fit.format())
        return self.fit

    def chained_lsqfit(
        self, data, prior, p0=None, print_fit=True, nterm=None,
        svdcut=None, svdnum=None, tol=None, maxit=None, 
        aux_param={},
        **args
        ):
        """ Compute chained least-squares fit of the correlator models to data.
            
        A *chained* fit fits data for each model sequentially, using the
        best-fit parameters of one fit as priors for fit parameters in the 
        next fit. Results from the individual fits can be found in dictionary
        ``self.chained_fits``, which is indexed by the ``m.datatag``\s for 
        the models ``m``.

        :param data: Input data. The ``datatag``\s from the 
            correlator models are used as data labels, with 
            ``data[datatag]`` being a 1-d array of |GVar|\s 
            corresponding to correlator values.
        :type data: dictionary
        :param prior: Bayesian prior for the fit parameters used in the 
            correlator models.
        :type prior: dictionary
        :param p0: A dictionary, indexed by parameter labels, containing 
            initial values for the parameters in the fit. Setting
            ``p0=None`` implies that initial values are extracted from the
            prior. Setting ``p0="filename"`` causes the fitter to look in
            the file with name ``"filename"`` for initial values and to
            write out best-fit parameter values after the fit (for the next
            call to ``self.lsqfit()``).
        :param print_fit: Print fit information to standard output if 
            ``True``; otherwise print nothing.
        :param aux_param: A dictionary containing priors for 
            additional (auxilary) fit parameters. Where ``aux_param``
            and ``prior`` share keys, the entry in ``aux_param`` is
            used.
        :type aux_param: dictionary

        The following parameters overwrite the values specified in the
        |CorrFitter| constructor when set to anything other than ``None``:
        ``nterm``, ``svdcut``, ``svdnum``, ``tol``, and ``maxit``. Any
        further keyword arguments are passed on to
        :func:`lsqfit.nonlinear_fit`, which does the fit.
        """        # setup
        if svdcut is None:
            svdcut = self.svdcut
        if svdnum is None:
            svdnum = self.svdnum
        if maxit is None:
            maxit = self.maxit
        if tol is None:
            tol = self.tol
        if nterm is None:
            nterm = self.nterm
        self.last_prior = prior
        self.aux_param = collections.OrderedDict(aux_param)

        # prepare data, do fits
        data = self.builddata(data=data, prior=prior, nterm=nterm)
        fits = collections.OrderedDict()
        aux_param = collections.OrderedDict(self.aux_param)
        p0file = p0 if isinstance(p0, str) else None
        if p0file is not None:
            try:
                p0 = lsqfit.nonlinear_fit.load_parameters(p0file)
            except (IOError, EOFError):
                p0 = None
        for m in self.models:
            m_prior = m.buildprior(prior=prior, nterm=nterm)
            for i in aux_param:
                m_prior[i] = aux_param[i]
            m_p0 = m.buildp0(p0, m_prior)

            @lsqfit.transform_p(m_prior, 0, 'p')
            def m_fitfcn(
                p, nterm=None, default_nterm=self.nterm, fitfcn = m.fitfcn
                ):
                if nterm is None:
                    nterm = default_nterm
                return fitfcn(p, nterm=nterm)
            lastfit = lsqfit.nonlinear_fit(
                data=data[m.datatag], fcn=m_fitfcn, prior=m_prior, 
                p0=m_p0, svdcut=svdcut, svdnum=svdnum, reltol=tol, 
                abstol=tol, maxit=maxit, **args
                )
            fits[m.datatag] = lastfit
            # print ('----------', m.datatag)
            # print (fits[-1].format())
            aux_param.update(lastfit.p)
        # print ('------------ all together')
        self.chained_fits = fits
        self.fit = copy.copy(lastfit)
        chi2 = 0.0
        dof = 0.0
        logGBF = 0.0
        nit = 0
        time = 0.0
        svdcorrection = collections.OrderedDict()
        all_y = _gvar.BufferDict()
        for key in self.chained_fits:
            f = self.chained_fits[key]
            all_y[key] = f.y
            chi2 += f.chi2
            dof += f.dof
            logGBF += f.logGBF
            nit += f.nit
            time += f.time
            for i in f.svdcorrection:
                if i not in svdcorrection or svdcorrection[i] is None:
                    svdcorrection[i] = f.svdcorrection[i]
                elif f.svdcorrection[i] is not None:
                    svdcorrection[i] = numpy.concatenate((
                        svdcorrection[i], f.svdcorrection[i]
                        ))
        self.fit.chi2 = chi2
        self.fit.dof = dof
        self.fit.logGBF = logGBF
        self.fit.Q = lsqfit.gammaQ(self.fit.dof/2., self.fit.chi2/2.) 
        self.fit.nit = int(nit / len(self.chained_fits))
        self.fit.time = time
        self.fit.y = all_y
        self.fit.data = data
        self.fit.prior = self.buildprior(prior, nterm=nterm)
        self.fit.fcn = self.buildfitfcn(prior.keys())
        for k in self.aux_param:
            self.fit.prior[k] = self.aux_param[k]
        self.fit.svdcorrection = svdcorrection
        if p0file is not None:
            self.fit.dump_pmean(p0file)
        if print_fit:
            print(self.fit.format())
        return self.fit

    def bootstrap_iter(self, datalist=None, n=None):
        """ Iterator that creates bootstrap copies of a |CorrFitter| fit using 
        bootstrap data from list ``data_list``.
            
        A bootstrap analysis is a robust technique for estimating means and
        standard deviations of arbitrary functions of the fit parameters.
        This method creates an interator that implements such an analysis
        of list (or iterator) ``datalist``, which contains bootstrap
        copies of the original data set. Each ``data_list[i]`` is a different
        ``data`` input for ``self.lsqfit()`` (that is, a dictionary containing
        fit data). The iterator works its way through the data sets in
        ``data_list``, fitting the next data set on each iteration and
        returning the resulting :class:`lsqfit.LSQFit` fit object. Typical
        usage, for an |CorrFitter| object named ``fitter``, would be::
            
            for fit in fitter.bootstrap_iter(datalist):
                ... analyze fit parameters in fit.p ...
                    
        :param data_list: Collection of bootstrap ``data`` sets for fitter. If
                ``None``, the data_list is generated internally using the 
                means and standard deviations of the fit data (assuming
                gaussian statistics).
        :type data_list: sequence or iterator or ``None``
        :param n: Maximum number of iterations if ``n`` is not ``None``;
                otherwise there is no maximum.
        :type n: integer
        :returns: Iterator that returns a :class:`lsqfit.LSQFit` object 
                containing results from the fit to the next data set in
                ``data_list``.
        """
        if datalist is not None:
            datalist = (self.builddata(d, self.last_prior) 
                        for d in datalist)
        for bs_fit in self.fit.bootstrap_iter(n, datalist=datalist):
            yield bs_fit
    
    def collect_fitresults(self):
        """ Collect results from last fit for plots, tables etc.
            
        :returns: A dictionary with one entry per correlator model,
            containing ``(t,G,dG,Gth,dGth)`` --- arrays containing::
                
                t       = times
                G(t)    = data averages for correlator at times t
                dG(t)   = uncertainties in G(t)
                Gth(t)  = fit function for G(t) with best-fit parameters
                dGth(t) = uncertainties in Gth(t)
        """
        corr = self.fit.data
        corrth = self.fit.fcn(self.fit.p)
        ans = {}
        keys = []
        for m in self.models:
            tag = m.datatag
            x = m._abscissa
            if x is None:
                continue
            c = corr[tag]
            cth = corrth[tag]
            keys.append(tag)
            ans[tag] = (x, _gvar.mean(c), _gvar.sdev(c),
                        _gvar.mean(cth), _gvar.sdev(cth))
        self.keys = keys
        return ans
    
    def display_plots(self):
        """ Show plots of data/fit-function for each correlator.
            
        Assumes :mod:`matplotlib` is installed (to make the plots). Plots
        are shown for one correlator at a time. Press key ``n`` to see the
        next correlator; press key ``p`` to see the previous one; press key
        ``q`` to quit the plot and return control to the calling program;
        press a digit to go directly to one of the first ten plots. Zoom,
        pan and save using the window controls.
        """
        import matplotlib.pyplot as plt
        data = self.collect_fitresults()
        keys = self.keys
        # keys = data.keys()
        # keys.sort()
        fig = plt.figure()
        idx = [0]
        def plotdata(idx, fig=fig, keys=keys, data=data):
            if idx[0] >= len(keys):
                idx[0] = len(keys)-1
            elif idx[0] < 0:
                idx[0] = 0
            i = idx[0]
            k = keys[i]
            t, g, dg, gth, dgth = data[k]
            fig.clear()
            plt.title("%d) %s   (press 'n', 'p', 'q' or a digit)"
                        % (i, k))
            ax = fig.add_subplot(111)
            ax.set_ylabel(str(k)+' / '+'fit')
            ax.set_xlim(min(t)-1,max(t)+1)
            # ax.set_xlabel('t')
            ii = (gth != 0.0)       # check for exact zeros (eg, antiperiodic)
            ax.errorbar(t[ii], g[ii]/gth[ii], dg[ii]/gth[ii], fmt='o')
            ax.plot(t, numpy.ones(len(t), float), 'r-')
            ax.plot(t[ii], 1+dgth[ii]/gth[ii], 'r--')
            ax.plot(t[ii], 1-dgth[ii]/gth[ii], 'r--') 
            plt.draw()
            # fig.draw()
        
        def onpress(event, idx=idx):
            try:    # digit?
                idx[0] = int(event.key)
            except ValueError:
                if event.key == 'n':
                    idx[0] += 1
                elif event.key == 'p':
                    idx[0] -= 1
                elif event.key == 'q':
                    plt.close()
                    return
                else:
                    return
            plotdata(idx)
        
        fig.canvas.mpl_connect('key_press_event', onpress)
        plotdata(idx)
        plt.show()
    


class fastfit(object):
    """ Fast fit for the leading component of a :mod:`Corr2`.
        
    This function class estimates ``En[0]`` and ``an[0]*bn[0]`` in a two-point 
    correlator::
        
        Gab(t) = sn * sum_i an[i]*bn[i] * fn(En[i], t)
               + so * sum_i ao[i]*bo[i] * fo(Eo[i], t)
        
    where ``sn`` and ``so`` are typically ``-1``, ``0``, or ``1`` and ::
        
        fn(E, t) =  exp(-E*t) + exp(-E*(tp-t)) # tp>0 -- periodic
               or   exp(-E*t) - exp(-E*(-tp-t))# tp<0 -- anti-periodic
               or   exp(-E*t)                  # if tp is None (nonperiodic)
        
        fo(E, t) = (-1)**t * fn(E, t)
        
    The correlator is specified by ``model``, and ``prior`` is used to 
    remove (marginalize) all terms other than the ``En[0]`` term
    from the data. This gives a *corrected* correlator ``Gc(t)`` that 
    includes uncertainties due to the terms removed. Estimates of ``En[0]`` 
    are given by::
        
        Eeff(t) = arccosh(0.5*(Gc(t+1)+Gc(t-1))/Gc(t)),
        
    The final estimate is the weighted average ``Eeff_avg`` of the
    ``Eeff(t)``\s for different ``t``\s. Similarly, an estimate for the
    product of amplitutes, ``an[0]*bn[0]`` is obtained from the weighted
    average of ::
        
        Aeff(t) = Gc(t)/fn(Eeff_avg, t). 
        
    If ``osc=True``, an estimate is returned for ``Eo[0]`` rather
    than ``En[0]``, and ``ao[0]*bo[0]`` rather than ``an[0]*bn[0]``. 
    These estimates are most reliable when ``Eo[0]`` is smaller than
    ``En[0]`` (and so dominates at large ``t``).
        
    The results of the fast fit are stored and returned in an object of type 
    :class:`corrfitter.fastfit` with the following attributies:
        
    .. attribute:: E
        
        Estimate of ``En[0]`` (or ``Eo[0]`` if ``osc==True``) computed
        from the weighted average of ``Eeff(t)`` for ``t``\s in
        ``model.tfit``. The prior is also included in the weighted average.
        
    .. attribute:: ampl
        
        Estimate of ``an[0]*bn[0]`` (or ``ao[0]*bo[0]`` if ``osc==True``)
        computed from the weighted average of ``Aeff(t)`` for ``t``\s in
        ``model.tfit[1:-1]``. The prior is also included in the weighted 
        average.
 
    .. attribute:: chi2
        
        ``chi[0]`` is the ``chi**2`` for the weighted average of
        ``Eeff(t)``\s; ``chi[1]`` is the same for the ``Aeff(t)``\s.
            
    .. attribute:: dof
        
        ``dof[0]`` is the effective number of degrees of freedom in the
        weighted average of ``Eeff(t)``\s; ``dof[1]`` is the same for the
        ``Aeff(t)``\s.
            
    .. attribute:: Q
        
        ``Q[0]`` is the quality factor `Q` for the weighted average of
        ``Eeff(t)``\s; ``Q[1]`` is the same for the ``Aeff(t)``\s.
            
    .. attribute:: Elist
        
        List of ``Eeff(t)``\s used in the weighted average to estimate
        ``E``.
            
    .. attribute:: ampllist
            
        List of ``Aeff(t)``\s used in the weighted average to estimate
        ``ampl``.
        
    :param data: Input data. The ``datatag`` from the correlator model is
        used as a data key, with ``data[datatag]`` being a 1-d array of
        |GVar|\s corresponding to the correlator values.
    :type data: dictionary
    :param prior: Bayesian prior for the fit parameters in the 
        correlator model.
    :type prior: dictionary
    :param model: Correlator model for correlator of interest. The ``t``\s
        in ``model.tfit`` must be consecutive.
    :type model: Corr2
    :param osc: If ``True``, extract results for the leading oscillating
        term in the correlator (``Eo[0]``); otherwise ignore.
    :type osc: Bool
        
    In addition an *svd* cut can be specified, as in |CorrFitter|, using
    parameters ``svdcut`` and ``svdnum``. Also the type of marginalization
    use can be specified with parameter ``ratio`` (see |CorrFitter|).
    """
    def __init__(self, data, prior, model, svdcut=None, svdnum=None, #):
                 ratio=True, osc=False):
        assert isinstance(model, Corr2), "model must be type Corr2"
        prior = _gvar.BufferDict(prior)
        t = model.tfit
        assert numpy.all(t[1:] == t[:-1]+1), "model.tfit must contain consecutive ts"
        self.osc = osc
        if osc:
            # capture leading oscillating part 
            nterm = (0, 1)
            Gfac = (-1)**t
            a = model.a[1]
            b = model.b[1]
            dE = model.dE[1]
        else:
            # capture leading non-oscillating part 
            nterm = (1, 0)
            Gfac = 1.
            a = model.a[0]
            b = model.b[0]
            dE = model.dE[0]
        # extract relevant data 
        fitter = CorrFitter(models=[model], svdcut=svdcut, svdnum=svdnum,
                            ratio=ratio, nterm=nterm)
        G = fitter.builddata(data=data, prior=prior)[model.datatag] * Gfac
        
        # transform prior
        prior = BaseModel._transform_prior(prior)
        # compute priors for answers 
        E_prior =  prior[dE][0]
        ampl_prior =  prior[a][0]
        ampl_prior *=  prior[b][0]
        
        # compute E 
        self.Elist = _gvar.arccosh(0.5*(G[2:]+G[:-2])/G[1:-1])
        Elist = self.Elist.tolist() + [E_prior]
        self.E = lsqfit.wavg(Elist, svdcut=svdcut, svdnum=svdnum)
        self.chi2 = lsqfit.wavg.chi2
        self.dof = lsqfit.wavg.dof
        self.Q = lsqfit.wavg.Q
        
        # compute amplitude 
        p = {}
        p[a] =  numpy.array([1.])
        p[b] =  numpy.array([1.])
        p[dE] = numpy.array([self.E])

        G0 = model.fitfcn(p, nterm=nterm) * Gfac
        self.G = G
        self.G0 = G0
        ii = slice(1, -1)
        self.ampllist = G[ii]/G0[ii]
        amplist = self.ampllist.tolist() + [ampl_prior]
        self.ampl = lsqfit.wavg(self.ampllist, svdcut=svdcut, svdnum=svdnum)
        self.chi2 = numpy.array([self.chi2, lsqfit.wavg.chi2])
        self.dof = numpy.array([self.dof, lsqfit.wavg.dof])
        self.Q = numpy.array([self.Q, lsqfit.wavg.Q])
    
       

