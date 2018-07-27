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
# Copyright (c) 2010-2018 G. Peter Lepage.
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

import collections
import fileinput
import math
import re
import time
import warnings

import lsqfit
import gvar as _gvar
import numpy

try:
    import scipy
    import scipy.linalg
except ImportError:
    scipy = None

def _parse_param(p,default=None):
    " Parse fit-parameter label "
    return p if isinstance(p, tuple) else (p, default)

class Corr2(lsqfit.MultiFitterModel):
    """ Two-point correlators ``Gab(t) = <b(t) a(0)>``.

    |Corr2| models the ``t`` dependence of a 2-point correlator ``Gab(t)``
    using ::

        Gab(t) = sn * sum_i an[i] * bn[i] * fn(En[i], t)
               + so * sum_i ao[i] * bo[i] * fo(Eo[i], t)

    where ``sn`` and ``so`` are typically ``-1``, ``0``, or ``1`` and ::

        fn(E, t) =  exp(-E*t) + exp(-E*(tp-t)) # tp>0 -- periodic
               or   exp(-E*t) - exp(-E*(-tp-t))# tp<0 -- anti-periodic
               or   exp(-E*t)                  # if tp is None (nonperiodic)

        fo(E, t) = (-1)**t * fn(E, t)

    The fit parameters for the non-oscillating piece of ``Gab`` (first term)
    are ``an[i]``, ``bn[i]``, and ``dEn[i]`` where::

        dEn[0] = En[0]
        dEn[i] = En[i]-En[i-1] > 0     (for i>0)

    and therefore ``En[i] = sum_j=0..i dEn[j]``. The fit parameters for
    the oscillating piece are defined analogously: ``ao[i]``, ``bo[i]``,
    and ``dEo[i]``.

    The fit parameters are specified by the keys corresponding to these
    parameters in a dictionary of priors supplied to |CorrFitter|. The keys
    are strings and are also used to access fit results. A log-normal
    prior can be specified for a parameter by including an entry for
    ``log(c)`` in the prior, rather than for ``c`` itself. See the
    :mod:`lsqfit` documentation for information about other distributions
    that are available. Values for both ``log(c)`` and ``c`` are
    included in the parameter dictionary. Log-normal distributions
    are  useful for forcing ``an``, ``bn`` and/or ``dE`` to be positive.

    When ``tp is not None`` and positive, the correlator is assumed to be
    symmetrical about ``tp/2``, with ``Gab(t)=Gab(tp-t)``. Data from
    ``t>tp/2`` is averaged with the corresponding data from ``t<tp/2``
    before fitting. When ``tp`` is negative, the correlator is assumed to
    be anti-symetrical about ``-tp/2``.

    Args:
        datatag (str): Key used to access correlator data in the input data
            dictionary (see |CorrFitter|): ``data[self.datatag]`` is a (1-d)
            array containing the correlator values (|GVar|\s).
        a (str or tuple): Key identifying the fit parameters for the source
            amplitudes ``an`` in the dictionary of priors provided to
            |CorrFitter|; or a two-tuple of keys for the source amplitudes
            ``(an, ao)``. The corresponding values in the dictionary of priors
            are (1-d) arrays of prior values with one term for each ``an[i]``
            or ``ao[i]``. Replacing either key by ``None`` causes the
            corresponding term to be dropped from the fit function. These keys
            are used to label  the corresponding parameter arrays in the fit
            results as well as  in the prior.
        b (str or tuple): Same as ``self.a`` but for the sinks ``(bn, bo)``
            instead of the sources ``(an, ao)``.
        dE (str): Key identifying the fit parameters for the energy
            differences ``dEn`` in the dictionary of priors provided by
            |CorrFitter|; or a two-tuple of keys for the energy differences
            ``(dEn, dEo)``. The corresponding values in the dictionary of
            priors are (1-d) arrays of prior values with one term for each
            ``dEn[i]`` or ``dEo[i]``. Replacing either key by ``None`` causes
            the corresponding term to be dropped from the fit function. These
            keys are used to label the corresponding parameter arrays in the
            fit results as well as in the prior.
        s (float or tuple): Overall factor ``sn`` for non-oscillating part
            of fit function, or two-tuple of overall factors ``(sn, so)``
            for both pieces.
        tdata (list of ints): The ``t``\s corresponding to data
            entries in the input data. Note that ``len(self.tdata)``
            should equal ``len(data[self.datatag])``. If ``tdata`` is
            omitted, ``tdata=numpy.arange(tp)`` is assumed, or
            ``tdata=numpy.arange(tmax)`` if ``tp`` is not specified.
        tfit (list of ints): List of ``t``\s to use in the fit. Only data
            with these ``t``\s (all of which should be in ``tdata``) is  used
            in the fit. If ``tfit`` is omitted, it is assumed to be all ``t``
            values from ``tdata`` that are larger than or equal to ``tmin``
            (if specified) and smaller than or equal to ``tmax`` (if
            specified).
        tp (int or ``None``): If ``tp`` is positive, the correlator
            is assumed to be periodic with ``Gab(t)=Gab(tp-t)``.
            If negative, the correlator is assumed to be anti-periodic
            with ``Gab(t)=-Gab(-tp-t)``. Setting ``tp=None`` implies that
            the correlator is not periodic, but rather continues
            to fall exponentially as ``t`` is increased indefinitely.
        tmin (int or ``None``): If ``tfit`` is omitted, it is assumed
            to be all ``t`` values from ``tdata`` that are larger than or
            equal to ``tmin`` and smaller than or equal to ``tmax``
            (if specified). ``tmin`` is ignored if ``tfit`` is specified.
        tmax (int or ``None``): If ``tfit`` is omitted, it is assumed
            to be all ``t`` values from ``tdata`` that are larger than or
            equal to ``tmin`` (if specified) and smaller than or
            equal to ``tmax``. ``tmin`` is ignored if ``tfit`` and ``tdata``
            are specified.
        ncg (int): Width of bins used to coarse-grain the correlator before
            fitting. Each bin of ``ncg`` correlator values is replaced by
            its average. Default is ``ncg=1`` (ie, no coarse-graining).
        reverse (bool): If ``True``, the data associated with ``self.datatag``
            is time-reversed (``data -> [data[0], data[-1], data[-2]...data[1]]``).
            Ignored otherwise.
        otherdata (str or list or ``None``): Data tag or list of data tags for
            additional data that are averaged with the ``self.datatag``
            data before fitting. This is useful including data
            from correlators with the source and sink interchanged.
            Default is ``None``.
        reverseddata (str or list or ``None``): Data tag or list of data tags
            for data that is time-reversed and then averaged with
            the ``self.datatag`` data before fitting. Default is ``None``.
    """
    def __init__(
        self, datatag, a, b, dE, tmin=None, tmax=None, tdata=None, tfit=None,
        tp=None, s=1.0, ncg=1, reverse=False, otherdata=None, reverseddata=None,
        othertags=None  # backwards compatibility
        ):
        super(Corr2, self).__init__(datatag, ncg)
        # othertags is the old name for otherdata
        if othertags is not None and otherdata is None:
            otherdata = othertags
        if isinstance(otherdata, str):
            otherdata = [otherdata]
        elif otherdata is None:
            otherdata = []
        if isinstance(reverseddata, str):
            reverseddata = [reverseddata]
        elif reverseddata is None:
            reverseddata = []
        self.otherdata = list(otherdata)
        self.reverseddata = list(reverseddata)
        self.reverse = reverse
        self.a = _parse_param(a)
        self.b = _parse_param(b)
        self.dE = _parse_param(dE)
        self.s = _parse_param(s, -1.)
        self.tp = tp

        # check consistency
        for x in zip(self.a, self.b, self.dE):
            x = set(x)
            if None in x and len(x) > 1:
                raise ValueError(
                    'inconsistent a, b and dE for ' + str(self.datatag)
                    )

        # figure out tdata
        if tdata is None:
            # assume tdata starts at 0
            if tp is not None:
                self.tdata = numpy.arange(abs(self.tp))
            elif tmax is not None:
                self.tdata = numpy.arange(tmax + 1)
            else:
                raise ValueError(
                    'need to specify tdata or tp or tmax for '
                    + str(self.datatag)
                    )
        else:
            self.tdata = tdata

        # figure out tfit
        if tfit is None:
            if tmin is None:
                tmin = numpy.min(self.tdata)
            if tmax is None:
                tmax = numpy.max(self.tdata)
            if self.tp is not None and tmax > abs(self.tp) // 2:
                tmax = abs(self.tp) // 2
            tfit = []
            for t in self.tdata:
                if t >= tmin and t <= tmax:
                    tfit.append(t)
            tfit = numpy.sort(numpy.array(tfit))
        new_tfit = []
        for t in tfit:
            if self.tp is None:
                if t not in self.tdata:
                    raise ValueError(
                        'tfit not contained in tdata for ' + str(self.datatag)
                        )
                new_tfit.append(t)
            else:
                t1, t2 = numpy.sort([t, abs(self.tp) - t])
                if t1 in new_tfit or t2 in new_tfit:
                    continue
                if t1 in self.tdata:
                    new_tfit.append(t1)
                elif t2 in self.tdata:
                    new_tfit.append(t2)
                else:
                    raise ValueError(
                        'tfit not contained in tdata for ' + str(self.datatag)
                        )
        self.tfit = numpy.array(new_tfit)
        if len(self.tfit) == 0:
            raise ValueError("empty tfit for " + str(self.datatag))

    def __str__(self):
        ans = "{c.datatag}[a={c.a}"
        for f in ['b', 'dE', 's', 'tp']:
            ans += ', ' + f + '={c.' + f +'}'
        ans += ', tfit=[{t1}...{t2}]]'
        return ans.format(c=self, t1=self.tfit[0], t2=self.tfit[-1])

    def buildprior(self, prior, nterm=None, mopt=None):
        """ Create fit prior by extracting relevant pieces from ``prior``.

        This routine selects the entries in dictionary ``prior``
        corresponding to the model's fit parameters. If ``nterm`` is
        not ``None``, it also adjusts the number of terms that are
        retained.

        Args:
            prior (dictionary): Dictionary containing priors for fit
                parameters.
            nterm (``None`` or int or two-tuple): Setting ``nterm=(n,no)``
                restricts the number of terms to ``n`` in the
                non-oscillating part and ``no`` in the oscillating part
                of the fit function. Replacing either or both by
                ``None`` keeps all terms, as does setting ``nterm=None``.
                This optional argument is used to implement
                marginalization.
        """
        # mopt=nterm = #terms to keep when marginalizing (None => all)
        # N.B. extend=True is built in but must allow argument
        if nterm is None:
            nterm = mopt
        nterm = _parse_param(nterm, None)
        newprior = _gvar.BufferDict()
        for ai, bi, dEi, ntermi in zip(self.a, self.b, self.dE, nterm):
            len_x = None
            for x in [ai, bi, dEi]:
                if x is None:
                    continue
                x = _gvar.dictkey(prior, x)
                newprior[x] = prior[x][None:ntermi]
                if len_x is None:
                    len_x = len(newprior[x])
                elif len(newprior[x]) != len_x:
                    raise ValueError(
                        'length mismatch between a, b, and dE for '
                         + str(self.datatag)
                         )
        return newprior

    def builddataset(self, dataset):
        """ Assemble fit data from data set dictionary ``dataset``. """
        tdata = list(self.tdata)
        tp = self.tp
        if tp is not None:
            pfac = math.copysign(1, tp)
            tp = abs(tp)
        def collect_data(odata):
            odata = numpy.asarray(odata)
            ndata = []
            for t in self.tfit:
                idt = tdata.index(t)
                if tp is None or tp-t not in tdata or t == tp-t:
                    ndata.append(odata[:, idt])
                else:
                    idt_r = tdata.index(tp-t)
                    ndata.append((odata[:, idt] + pfac * odata[:, idt_r]) / 2.)
            return numpy.transpose(ndata)
        ans = []
        if self.reverse:
            otherdata = self.otherdata
            reverseddata = [self.datatag] + self.reverseddata
        else:
            otherdata = [self.datatag] + self.otherdata
            reverseddata = self.reverseddata
        for tag in otherdata:
            ans.append(collect_data(dataset[tag]))
        for tag in reverseddata:
            rdset = numpy.empty(dataset[tag].shape, dtype=float)
            rdset[:, 0] = dataset[tag][:, 0]
            rdset[:, 1:] = dataset[tag][:, -1:0:-1]
            ans.append(collect_data(rdset))
        return (
            numpy.array(ans[0]) if len(ans) == 1 else
            numpy.average(ans, axis=0)
            )

    def builddata(self, data):
        """ Assemble fit data from dictionary ``data``. """
        tdata = list(self.tdata)
        tp = self.tp
        if tp is not None:
            pfac = math.copysign(1, tp)
            tp = abs(tp)
        def collect_data(odata):
            ndata = []
            for t in self.tfit:
                idt = tdata.index(t)
                if tp is None or tp-t not in tdata or t == tp-t:
                    ndata.append(odata[idt])
                else:
                    idt_r = tdata.index(tp-t)
                    ndata.append(lsqfit.wavg([odata[idt], pfac*odata[idt_r]]))
            return ndata
        if self.reverse:
            otherdata = self.otherdata
            reverseddata = [self.datatag] + self.reverseddata
        else:
            otherdata = [self.datatag] + self.otherdata
            reverseddata = self.reverseddata
        ans = []
        for tag in otherdata:
            ans.append(collect_data(data[tag]))
        for tag in reverseddata:
            rdata = numpy.empty(data[tag].shape, dtype=data[tag].dtype)
            rdata[0] = data[tag][0]
            rdata[1:] = data[tag][-1:0:-1]
            ans.append(collect_data(rdata))
        return numpy.array(ans[0]) if len(ans) == 1 else lsqfit.wavg(ans)

    def fitfcn(self, p, t=None):
        """ Return fit function for parameters ``p``. """
        # setup
        if t is None:
            t = self.tfit
        else:
            t = numpy.asarray(t)
            if len(t.shape) != 1:
                raise ValueError('t must be 1-d array')
        if self.tp is None:
            tp_t = None
        elif self.tp >= 0:
            tp_t = self.tp - t
            pfac = 1
        else:
            tp_t = -self.tp - t
            pfac = -1
        ofac = (None if self.s[0] == 0.0 else self.s[0],
                (None if self.s[1] == 0.0 else self.s[1]*(-1)**t))

        # calculate function
        ans = 0.0
        for ai, bi, dEi, ofaci in zip(self.a, self.b, self.dE, ofac):
            if ai is None or ofaci is None:
                continue
            if tp_t is None:
                sumdE = numpy.cumsum(p[dEi])
                ans += ofaci * numpy.sum(
                    (p[ai] * p[bi])[None, :]
                    * _gvar.exp(-t[:, None] * sumdE[None, :]),
                    axis=1,
                    )
            else:
                sumdE = numpy.cumsum(p[dEi])
                ans += ofaci * numpy.sum(
                    (p[ai] * p[bi])[None, :] * (
                        _gvar.exp(-t[:, None] * sumdE[None, :])
                        + pfac * _gvar.exp(-tp_t[:, None] * sumdE[None, :])
                        ),
                    axis=1,
                    )
        return ans

class Corr3(lsqfit.MultiFitterModel):
    """ Three-point correlators ``Gavb(t, T) = <b(T) V(t) a(0)>``.

    |Corr3| models the ``t`` dependence of a 3-point correlator
    ``Gavb(t, T)`` using ::

        Gavb(t, T) =
         sum_i,j san * an[i] * fn(Ean[i],t) * Vnn[i,j] * sbn * bn[j] * fn(Ebn[j],T-t)
        +sum_i,j san * an[i] * fn(Ean[i],t) * Vno[i,j] * sbo * bo[j] * fo(Ebo[j],T-t)
        +sum_i,j sao * ao[i] * fo(Eao[i],t) * Von[i,j] * sbn * bn[j] * fn(Ebn[j],T-t)
        +sum_i,j sao * ao[i] * fo(Eao[i],t) * Voo[i,j] * sbo * bo[j] * fo(Ebo[j],T-t)

    where ::

        fn(E, t) =  exp(-E*t)
        fo(E, t) = (-1)**t * exp(-E*t)

    The fit parameters for the non-oscillating piece of ``Gavb`` (first term)
    are ``Vnn[i,j]``, ``an[i]``, ``bn[j]``, ``dEan[i]`` and ``dEbn[j]`` where,
    for example::

        dEan[0] = Ean[0]
        dEan[i] = Ean[i] - Ean[i-1] > 0     (for i>0)

    and therefore ``Ean[i] = sum_j=0..i dEan[j]``. The parameters for the
    other terms are similarly defined.

    Args:
        datatag (str): Tag used to label correlator in the input data.
        a (str or tuple): Key identifying the fit parameters for the source
            amplitudes ``an``, for ``a->V``, in the dictionary of priors
            provided to |CorrFitter|; or a two-tuple of keys for the source
            amplitudes ``(an, ao)``. The corresponding values in the
            dictionary of priors are (1-d) arrays of prior values with one
            term for each ``an[i]`` or ``ao[i]``. Replacing either key by
            ``None`` causes the corresponding term to be dropped from the fit
            function. These keys are used to label the corresponding parameter
            arrays in the fit results as well as in the prior.
        b (str or tuple): Same as ``self.a`` but  for the ``V->b`` sink
            amplitudes ``(bn, bo)``.
        dEa (str or tuple): Fit-parameter label for ``a->V``
            intermediate-state energy differences ``dEan``, or two-tuple of
            labels for the differences ``(dEan, dEao)``. Each label represents
            an array of energy differences. Replacing either label by ``None``
            causes the corresponding term in the correlator function to be
            dropped. These keys are used to label the corresponding parameter
            arrays in the fit results as well as in the prior.
        dEb (str or tuple): Same as ``self.dEa`` but for ``V->b`` sink
            energies ``(dEbn, dEbo)``.
        sa (float or tuple): Overall factor ``san`` for the non-oscillating
            ``a->V`` terms in the correlator, or two-tuple containing
            the overall factors ``(san, sao)`` for the non-oscillating and
            oscillating terms. Default is ``(1,-1)``.
        sb (float or tuple): Same as ``self.sa`` but for ``V->b`` sink
            overall factors ``(sbn, sbo)``.
        Vnn (str or ``None``): Fit-parameter label for the matrix of current
            matrix elements ``Vnn[i,j]`` connecting non-oscillating states.
            The matrix must be square and symmetric if ``symmetric_V=True``,
            and only the elements ``V[i,j]`` for ``j>=i`` are specified, using
            a 1-d array ``V_sym`` with the following layout::

                [V[0,0],V[0,1],V[0,2]...V[0,N],
                        V[1,1],V[1,2]...V[1,N],
                               V[2,2]...V[2,N],
                                     .
                                      .
                                       .
                                        V[N,N]]

            Note that ``V[i,j] = V_symm[i*N + j - i * (i+1) / 2]`` for
            ``j>=i``. Set ``Vnn=None`` to omit it.
        Vno (str or ``None``): Fit-parameter label for the matrix of current
            matrix elements ``Vno[i,j]`` connecting non-oscillating to
            oscillating states. Only one of ``Von`` and ``Vno`` can be
            specified if ``symmetric_V=True``; the other is defined to be its
            transform. Set ``Vno=None`` to omit it.
        Von (str or ``None``): Fit-parameter label for the matrix of current
            matrix elements ``Vno[i,j]`` connecting oscillating to non-
            oscillating states. Only one of ``Von`` and ``Vno`` can be
            specified if ``symmetric_V=True``; the other is defined to be its
            transform. Set ``Von=None`` to omit it.
        Voo (str or ``None``): Fit-parameter label for the matrix of current
            matrix elements ``Voo[i,j]`` connecting oscillating states. The
            matrix must be square and symmetric if ``symmetric_V=True``, and
            only the elements ``V[i,j]`` for ``j>=i`` are specified, using a
            1-d array ``V_sym`` with the following layout::

                [V[0,0],V[0,1],V[0,2]...V[0,N],
                        V[1,1],V[1,2]...V[1,N],
                               V[2,2]...V[2,N],
                                     .
                                      .
                                       .
                                        V[N,N]]

            Note that ``V[i,j] = V_symm[i*N + j - i * (i+1) / 2]`` for
            ``j>=i``. Set ``Voo=None`` to omit it.
        reverse (bool): If ``True``, the data associated with ``self.datatag``
            is time-reversed before fitting (interchanging ``t=0`` with
            ``t=T``). This is useful for doing simultaneous fits to
            ``a->V->b`` and ``b->V->a``, where one is time-reversed relative
            to the other: *e.g.*, ::

                models = [ ...
                    Corr3(
                        datatag='a->V->b', tmin=3, T=15,
                        a=('a', 'ao'), dEa=('dEa', 'dEao'),
                        b=('b', 'bo'), dEb=('dEb', 'dEbo'),
                        Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo',
                        ),
                    Corr3(
                        datatag='b->V->a', tmin=3, T=15,
                        a=('a', 'ao'), dEa=('dEa', 'dEao'),
                        b=('b', 'bo'), dEb=('dEb', 'dEbo'),
                        Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo',
                        reverse=True,
                        ),
                    ...
                    ]

            Another (faster) strategy for such situations is to average
            data from the second process with that from the  first, before
            fitting, using keyword ``reverseddata``. Default is ``False``.
        symmetric_V (bool): If ``True``, the fit function for ``a->V->b`` is
            unchanged (symmetrical) under the the interchange of ``a`` and
            ``b``. Then ``Vnn`` and ``Voo`` are square, symmetric matrices
            and their priors are one-dimensional arrays containing only
            elements ``V[i,j]`` with ``j>=i``, as discussed above.
            Only one of ``Von`` and ``Vno`` can be specified if
            ``symmetric_V=True``; the other is defined to be its transform.
        T (int): Separation between source and sink.
        tdata (list of ints): The ``t``\s corresponding to data entries
            in the input data. If omitted, is assumed equal to
            ``numpy.arange(T + 1)``.
        tfit (list of ints): List of ``t``\s to use in the fit. Only data
            with these ``t``\s (all of which should be in ``tdata``)
            is used in the fit. If omitted, is assumed equal to
            ``numpy.arange(tmin, T - tmin + 1)``.
        tmin (int or ``None``): If ``tfit`` is omitted, it is set equal
            to ``numpy.arange(tmin, T - tmin + 1)``. ``tmin`` is ignored
            if ``tfit`` is specified.
        ncg (int): Width of bins used to coarse-grain the correlator before
            fitting. Each bin of ``ncg`` correlator values is replaced by
            its average. Default is ``ncg=1`` (ie, no coarse-graining).
        reverseddata (str or list or ``None``): Data tag or list of data tags
            for additional data that are time-reversed and then averaged with
            the ``self.datatag`` data before fitting. This is useful for
            folding data from ``b->V->a`` into a fit for ``a->V->b``:
            *e.g.*, ::

                Corr3(
                    datatag='a->V->b',
                    a=('a', 'ao'), dEa=('dEa', 'dEao'),
                    b=('b', 'bo'), dEb=('dEb', 'dEbo'),
                    Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo',
                    tmin=3, T=15, reverseddata='b->V->a'
                    ),

            This is faster than using a separate model with ``reverse=True``.
            Default is ``None``.
        otherdata (str or list or ``None``): Data tag or list of data tags
            for additional data that are averaged with the ``self.datatag``
            data before fitting. Default is ``None``.
    """
    def __init__(
        self, datatag, T=None, tdata=None, tfit=None, tmin=None,
        a=None, b=None, dEa=None, dEb=None, sa=1., sb=1.,
        Vnn=None, Vno=None, Von=None, Voo=None,
        reverse=False, symmetric_V=False, transpose_V=None, ncg=1,
        reverseddata=None, otherdata=None,
        tpa=None, tpb=None, othertags=None # backwards compatibility; tp's ignored
        ):
        super(Corr3, self).__init__(datatag, ncg)
        if tpa is not None or tpb is not None:
            warnings.warn('parameters tpa and tbp are ignored (obsolete)')
        # othertags is the old name for otherdata
        if othertags is not None and otherdata is None:
            otherdata = othertags
        if isinstance(otherdata, str):
            otherdata = [otherdata]
        elif otherdata is None:
            otherdata = []
        if isinstance(reverseddata, str):
            reverseddata = [reverseddata]
        elif reverseddata is None:
            reverseddata = []
        self.reverseddata = list(reverseddata)
        self.reverse = reverse
        self.otherdata = list(otherdata)
        self.a = _parse_param(a)
        self.dEa = _parse_param(dEa)
        self.sa = _parse_param(sa, -1.)
        self.b = _parse_param(b)
        self.dEb = _parse_param(dEb)
        self.sb = _parse_param(sb, -1.)
        self.transpose_V = transpose_V
        if transpose_V is not None:
            warnings.warn("'transpose_V' is deprecated; use 'reverse' keyword instead")
        self.symmetric_V = symmetric_V
        if self.symmetric_V:
            # use transpose of Vno (if present) for Von (or vice versa)
            if Vno is not None:
                Von = None
        self.V = [[Vnn, Vno], [Von, Voo]]
        self.T = abs(T)

        # consistency checks
        for x in zip(self.a, self.dEa):
            x = set(x)
            if None in x and len(x) > 1:
                raise ValueError(
                    'inconsistent a and dEa for ' + str(self.datatag)
                    )
        for x in zip(self.b, self.dEb):
            x = set(x)
            if None in x and len(x) > 1:
                raise ValueError(
                    'inconsistent b and dEb for ' + str(self.datatag)
                    )

        for i in range(2):
            for j in range(2):
                if self.V[i][j] is not None:
                    if None in [self.a[i], self.b[j]]:
                        raise ValueError(
                            'inconsistent a, b and V for '
                            + str(self.datatag)
                            )
                else:
                    if self.symmetric_V and i != j:
                        if self.V[j][i] is not None:
                            continue
                    if None not in [self.a[i], self.b[j]]:
                        raise ValueError(
                            'inconsistent a, b and V for ' + str(self.datatag)
                            )

        # tdata amd tfit
        if tdata is None:
            tdata = numpy.arange(self.T + 1)
        self.tdata = tdata
        if tfit is None:
            if tmin is None:
                self.tfit = numpy.array(self.tdata)
            else:
                new_tfit = []
                tmin = abs(tmin)
                for t in self.tdata:
                    if t < tmin or self.T - t < tmin:
                        continue
                    new_tfit.append(t)
                self.tfit = numpy.array(new_tfit)
        else:
            new_tfit = []
            for t in tfit:
                if t in self.tdata:
                    new_tfit.append(t)
                else:
                    raise ValueError(
                        'tfit not contained in tdata for ' + str(self.datatag)
                        )
        self.tfit = numpy.sort(numpy.array(new_tfit))
        if len(self.tfit) == 0:
            raise ValueError("empty tfit for " + str(self.datatag))

    def buildprior(self, prior, mopt=None, nterm=None):
        if nterm is None:
            nterm = mopt
        nterm = _parse_param(nterm, None)
        def resize_sym(Vii, ntermi):
            # N = size of Vii; ntermi is new max. dimension
            N = int(numpy.round((((8. * len(Vii) + 1.) ** 0.5 - 1.) / 2.)))
            if ntermi is None or N == ntermi:
                return Vii
            ans = []
            iterV = iter(Vii)
            for i in range(N):
                for j in range(i, N):
                    v = next(iterV)
                    if j < ntermi:
                        ans.append(v)
            return numpy.array(ans)
        ans = _gvar.BufferDict()

        # unpack propagator parameters
        for x in [self.a, self.dEa, self.b, self.dEb]:
            for xi, ntermi in zip(x, nterm):
                if xi is not None:
                    xi = _gvar.dictkey(prior, xi)
                    ans[xi] = prior[xi][None:ntermi]

        # i,j range from n to o
        for i in range(2):
            for j in range(2):
                vij = self.V[i][j]
                if vij is None:
                    continue
                vij = _gvar.dictkey(prior, vij)
                if i == j and self.symmetric_V:
                    ans[vij] = (
                        prior[vij] if nterm[i] is None else
                        resize_sym(prior[vij], nterm[i])
                        )
                else:
                    if self.transpose_V:
                        ans[vij] = prior[vij][None:nterm[j], None:nterm[i]]
                    else:
                        ans[vij] = prior[vij][None:nterm[i], None:nterm[j]]

        # verify dimensions
        for ai, dEai in zip(self.a, self.dEa):
            if ai is None:
                continue
            ai, dEai = _gvar.get_dictkeys(prior, [ai, dEai])
            if len(ans[ai]) != len(ans[dEai]):
                raise ValueError(
                    'length mismatch between a and dEa for '
                    + str(self.datatag)
                    )
        for bj, dEbj in zip(self.b, self.dEb):
            if bj is None or dEbj is None:
                continue
            bj, dEbj = _gvar.get_dictkeys(prior, [bj, dEbj])
            if len(ans[bj]) != len(ans[dEbj]):
                raise ValueError(
                    'length mismatch between b and dEb for '
                    + str(self.datatag)
                    )
        for i in range(2):
            for j in range(2):
                Vij = self.V[i][j]
                if Vij is None:
                    continue
                ai, bj, Vij = _gvar.get_dictkeys(
                    prior, [self.a[i], self.b[j], Vij]
                    )
                if i == j and self.symmetric_V:
                    N = ans[ai].shape[0]
                    if ans[bj].shape[0] != N:
                        raise ValueError(
                            'length mismatch between a, b, and V for '
                            + str(self.datatag)
                            )
                    if len(ans[Vij].shape) != 1:
                        raise ValueError(
                            'symmetric_V=True => Vnn, Voo = 1-d arrays for '
                            + str(self.datatag)
                            )
                    if ans[Vij].shape[0] !=  (N * (N+1)) / 2:
                        raise ValueError(
                            'length mismatch between a, b, and V for '
                            + str(self.datatag)
                            )
                else:
                    ai, bj, Vij = _gvar.get_dictkeys(
                        prior, [self.a[i], self.b[j], Vij]
                        )
                    Vij_shape = (
                        ans[Vij].shape[::-1] if self.transpose_V else
                        ans[Vij].shape
                        )
                    if ans[ai].shape[0] != Vij_shape[0]:
                        raise ValueError(
                            'length mismatch between a and V for '
                            + str(self.datatag)
                            )
                    elif ans[bj].shape[0] != Vij_shape[1]:
                        raise ValueError(
                            'length mismatch between b and V for '
                            + str(self.datatag)
                            )
        return ans

    def builddataset(self, dataset):
        tdata = list(self.tdata)
        def collect_data(odata):
            odata = numpy.asarray(odata)
            ndata = []
            for t in self.tfit:
                ndata.append(odata[:, tdata.index(t)])
            return numpy.transpose(ndata)
        if self.reverse:
            ans = [collect_data(numpy.asarray(dataset[self.datatag])[:, ::-1])]
        else:
            ans = [collect_data(dataset[self.datatag])]
        for tag in self.otherdata:
            ans.append(collect_data(dataset[tag]))
        for tag in self.reverseddata:
            ans.append(collect_data(numpy.asarray(dataset[tag])[:, ::-1]))
        return (
            numpy.array(ans[0]) if len(ans) == 1 else
            numpy.average(ans, axis=0)
            )

    def builddata(self, data):
        tdata = list(self.tdata)
        def collect_data(odata):
            ndata = []
            for t in self.tfit:
                ndata.append(odata[tdata.index(t)])
            return ndata
        ans = []
        if self.reverse:
            ans = [collect_data(data[self.datatag][::-1])]
        else:
            ans = [collect_data(data[self.datatag])]
        for tag in self.otherdata:
            ans.append(collect_data(data[tag]))
        for tag in self.reverseddata:
            ans.append(collect_data(data[tag][::-1]))
        return numpy.array(ans[0]) if len(ans) == 1 else lsqfit.wavg(ans)

    def fitfcn(self, p, t=None):
        # setup
        if t is None:
            t = self.tfit
        ta = numpy.asarray(t)
        if len(ta.shape) != 1:
            raise ValueError('t must be 1-d array')
        tb = self.T - ta

        # initial propagators
        # aprop[i][j][t] where i = n or o and j = excitation level
        aprop = []
        ofac = (
            None if self.sa[0] == 0.0 else self.sa[0],
            None if self.sa[1] == 0.0 else self.sa[1] * (-1)**ta[None, :]
            )
        for ai, dEai, ofaci in zip(self.a, self.dEa, ofac):
            if ai is None or ofaci is None:
                aprop.append(None)
                continue
            sumdE = numpy.cumsum(p[dEai])
            aprop.append(
                ofaci * p[ai][:, None] * _gvar.exp(-ta[None, :] * sumdE[:, None])
                )
        aprop = numpy.array(aprop)

        # final propagators
        # bprop[i][j][t] where i = n or o and j = excitation level
        bprop = []
        ofac = (
            None if self.sb[0] == 0.0 else self.sb[0],
            None if self.sb[1] == 0.0 else self.sb[1] * (-1)**tb[None, :]
            )
        for bj, dEbj, ofacj in zip(self.b, self.dEb, ofac):
            if bj is None or ofacj is None:
                bprop.append(None)
                continue
            sumdE = numpy.cumsum(p[dEbj])
            bprop.append(
                ofacj * p[bj][:, None] * _gvar.exp(-tb[None, :] * sumdE[:, None])
                )
        bprop = numpy.array(bprop)

        # combine with vertices
        # combine propagators with vertices
        ans = 0.0
        for i, (apropi, Vi) in enumerate(zip(aprop, self.V)):
            if apropi is None:
                continue
            for j, (bpropj, Vij) in enumerate(zip(bprop, Vi)):
                if bpropj is None:
                    continue
                elif i != j and self.symmetric_V:
                    # Von is Vno.T or vice versa
                    if Vij is None:
                        Vij = self.V[j][i]
                        if Vij is None:
                            continue
                        V = p[Vij].T
                    else:
                        V = p[Vij]
                elif Vij is None:
                    continue
                else:
                    V = p[Vij]
                if i == j and self.symmetric_V:
                    # unpack symmetric matrix V (assumed square)
                    na = len(apropi)
                    iterV = iter(V)
                    V = numpy.empty((na, na), dtype=V.dtype)
                    for k in range(na):
                        for l in range(k, na):
                            V[k, l] = next(iterV)
                            if k != l:
                                V[l, k] = V[k, l]
                if self.transpose_V:
                    V = V.T
                if min(len(apropi), len(V), len(bpropj)) != 0:
                    ans += numpy.sum(apropi * numpy.dot(V, bpropj), axis=0)
        return ans

class CorrFitter(lsqfit.MultiFitter):
    """ Nonlinear least-squares fitter for a collection of correlator models.

    Args:
        models: List of models, derived from :mod:`lsqfit.MultiFitterModel`,
            to be fit to the data. Individual models in the list can
            be replaced by lists of models or tuples of models; see below.
        nterm (tuple or int or None): Number of terms fit in the
            non-oscillating part of fit functions; or a two-tuple of
            numbers indicating how many terms to fit in each of the
            non-oscillating and oscillating parts. Terms  omitted from the
            fit are marginalized (*i.e.*, included as corrections to the
            fit data). If set to ``None``, all parameters in the
            prior are fit, and none are marginalized.
        mopt: Marginalization options; alias for ``nterm``.
        ratio (bool): If ``True``, implement marginalization using
            ratios: ``data_marg = data * fitfcn(prior_marg) / fitfcn(prior)``.
            If ``False`` (default), implement using differences:
            ``data_marg = data + (fitfcn(prior_marg) - fitfcn(prior))``.
        fast (bool): Setting ``fast=True`` (default) strips any variable
            not required by the fit from the prior. This speeds
            fits but loses information about correlations between
            variables in the fit and those that are not.
        fitterargs: Additional arguments for the :class:`lsqfit.nonlinear_fit`,
            such as ``tol``, ``maxit``, ``svdcut``, ``fitter``, etc., as needed.
    """
    def __init__(self, models, **kargs):
        if 'ratio' not in kargs:
            kargs['ratio'] = False
        if 'fast' not in kargs:
            kargs['fast'] = True
        if 'nterm' in kargs:
            kargs['mopt'] = kargs['nterm']
        super(CorrFitter, self).__init__(models=models, **kargs)
        # replace nterm by mopt
        for tasktype, taskdata in self.tasklist:
            if tasktype == 'update-kargs' and 'nterm' in taskdata:
                taskdata['mopt'] = taskdata['nterm']
                del taskdata['nterm']

    def lsqfit(self, data=None, prior=None, pdata=None, p0=None, **kargs):
        """ Compute least-squares fit of models to data.

        :meth:`CorrFitter.lsqfit` fits all of the models together, in
        a single fit. It returns the :class:`lsqfit.nonlinear_fit` object
        from the fit.

        See documentation for :mod:`lsqfit.MultiFitter.lsqfit` for
        more information.

        Args:
            data: Input data. One of ``data`` or ``pdata`` must be
                specified but not both.
            pdata: Input data that has been processed by the
                models using :meth:`CorrFitter.process_data` or
                :meth:`CorrFitter.process_dataset`. One of
                ``data`` or ``pdata`` must be  specified but not both.
                ``pdata`` is obtained from ``data`` by collecting the output
                from ``m.builddata(data)`` for each model ``m`` and storing it
                in a dictionary with key ``m.datatag``.
            prior: Bayesian prior for fit parameters used by the models.
            p0: Dictionary , indexed by parameter labels, containing
                initial values for the parameters in the fit. Setting
                ``p0=None`` implies that initial values are extracted from the
                prior. Setting ``p0="filename"`` causes the fitter to look in
                the file with name ``"filename"`` for initial values and to
                write out best-fit parameter values after the fit (for the
                next call to ``self.lsqfit()``).
            wavg_svdcut (float):  SVD cut used in weighted averages for
                parallel fits.
            kargs: Arguments that override parameters specified when
                the :class:`CorrFitter` was created. Can also include
                additional arguments to be passed through to
                the :mod:`lsqfit` fitter.
        """
        if 'nterm' in kargs:
            kargs['mopt'] = kargs['nterm']
            del kargs['nterm']
        return super(CorrFitter, self).lsqfit(
            data=data, prior=prior, pdata=pdata, p0=p0, **kargs
            )

    def chained_lsqfit(
        self, data=None, prior=None, pdata=None, p0=None, **kargs
        ):
        """ Compute chained least-squares fit of models to data.

        :meth:`CorrFitter.chained_lsqfit` fits the models specified
        in ``models`` one at a time, in sequence, with the fit output
        from one being fed into the prior for the next.

        See documentation for :mod:`lsqfit.MultiFitter.chained_lsqfit` for
        (much) more information.

        Args:
            data: Input data. One of ``data`` or ``pdata`` must be
                specified but not both. ``pdata`` is obtained from ``data``
                by collecting the output from ``m.builddata(data)``
                for each model ``m`` and storing it in a dictionary
                with key ``m.datatag``.
            pdata: Input data that has been processed by the
                models using :meth:`CorrFitter.process_data` or
                :meth:`CorrFitter.process_dataset`. One of
                ``data`` or ``pdata`` must be  specified but not both.
            prior: Bayesian prior for fit parameters used by the models.
            p0: Dictionary , indexed by parameter labels, containing
                initial values for the parameters in the fit. Setting
                ``p0=None`` implies that initial values are extracted from the
                prior. Setting ``p0="filename"`` causes the fitter to look in
                the file with name ``"filename"`` for initial values and to
                write out best-fit parameter values after the fit (for the
                next call to ``self.lsqfit()``).
            kargs: Arguments that override parameters specified when
                the :class:`CorrFitter` was created. Can also include
                additional arguments to be passed through to
                the :mod:`lsqfit` fitter.
        """
        if 'nterm' in kargs:
            kargs['mopt'] = kargs['nterm']
            del kargs['nterm']
        return super(CorrFitter, self).chained_lsqfit(
            data=data, prior=prior, pdata=pdata, p0=p0, **kargs
            )

    def set(self, **kargs):
        """ Reset default keyword parameters.

        Assigns new default values from dictionary ``kargs`` to the fitter's
        keyword parameters. Keywords for the underlying :mod:`lsqfit` fitters
        can also be  included (or grouped together in dictionary
        ``fitterargs``).

        Returns tuple ``(kargs, oldkargs)`` where ``kargs`` is a dictionary
        containing all :class:`lsqfit.MultiFitter` keywords after they have
        been updated, and ``oldkargs`` contains the  original values for these
        keywords. Use ``fitter.set(**oldkargs)`` to restore the original
        values.
        """
        if 'nterm' in kargs:
            kargs['mopt'] = kargs['nterm']
            del kargs['nterm']
        return super(CorrFitter, self).set(**kargs)

    def display_plots(self, save=False, view='ratio'):
        """ Displays correlator plots.

        Deprecated. Use ``fit.show_plots(save, view)`` instead.
        See documentation for :mod:`lsqfit.MultiFitter` for more
        information.
        """
        warnings.warn(
            'display_plots is deprecated; use ``fit.show_plots(...) instead',
            DeprecationWarning,
            )
        self.fit.show_plots(save=save, view=view)

    @staticmethod
    def read_dataset(
        inputfiles, grep=None, keys=None, h5group='/', binsize=1,
        tcol=0, Gcol=1,
        ):
        """ Read correlator Monte Carlo data from files into a :class:`gvar.dataset.Dataset`.

        Three files formats are supported by :func:`read_dataset`, depending
        upon ``inputfiles``.

        If ``inputfiles`` is a string ending in ``'.h5'``, it is assumed to
        be the name of a file in hpf5 format. The file is opened as
        ``h5file`` and all hpf5 datasets in ``h5file[h5group]`` are
        collected into a dictionary and returned.

        The second file format is the text-file format supported by
        :class:`gvar.dataset.Dataset`: each line consists of a  tag or key
        identifying a correlator followed by data corresponding to  a single
        Monte Carlo measurement of the correlator. This format is assumed if
        ``inputfiles`` is a filename or a list of filenames. It allows a
        single file to contain an arbitrary number of measurements for an
        arbitrary number of different correlators. The data can also be spread
        over multiple files. A typical file might look like ::

            # this is a comment; it is ignored
            aa 1.237 0.912 0.471
            bb 3.214 0.535 0.125
            aa 1.035 0.851 0.426
            bb 2.951 0.625 0.091
            ...

        which describes two correlators, ``aa`` and ``bb``, each having
        three different ``t`` values.

        The third file format is assumed when ``inputfiles`` is a dictionary. The
        dictionary's keys and values identify the (one-dimensional) correlators
        and the files containing their Monte Carlo data, respectively. So the
        data for correlators ``aa`` and ``bb`` above are in separate files::

            fileinputs = dict(aa='aafile', bb='bbfile')

        Each line in these data files consists of an index ``t`` value followed by
        the corresponding value for correlator ``G(t)``.  The ``t``\s increase
        from line to line up to their maximum value,  at which point they repeat.
        The ``aafile`` file for correlator ``aa`` above  would look like::

            # this is a comment; it is ignored
            1 1.237
            2 0.912
            3 0.471
            1 1.035
            2 0.851
            3 0.426
            ...

        The columns in these files containing ``t`` and ``G(t)`` are
        assumed to be columns 0 and 1, respectively. These can be changed
        by setting arguments ``tcol`` and ``Gcol``, respectively.

        ``corrfitter.process_dataset`` supports keywords ``binsize``,
        ``grep`` and ``keys``. If ``binsize`` is greater than one,
        random samples are binned with bins of size ``binsize``.
        If ``grep`` is not ``None``, only keys that match or partially
        match regular expression ``grep`` are retained; others are ignored.
        If ``keys`` is not ``None``, only keys that are in list ``keys``
        are retained; others are discarded.
        """
        if not hasattr(inputfiles, 'keys'):
            # inputfiles is a filename or list of filenames (or files)
            if h5group == [] or h5group is None:   # not needed after next gvar update
                h5group = '/'
            dset = _gvar.dataset.Dataset(
                inputfiles, binsize=binsize, grep=grep,
                h5group=h5group, keys=keys,
                )
        else:
            # inputfiles is a dictionary
            # files are in t-G format
            if grep:
                grep = re.compile(grep)
            dset = _gvar.dataset.Dataset()
            for k in inputfiles:
                if keys and k not in keys:
                    continue
                if grep and grep.search(k) is None:
                    continue
                tlast = - float('inf')
                G = []
                for line in fileinput.input([inputfiles[k]]):
                    f = line.split()
                    if f[0][0] == '#':
                        # comment
                        continue
                    t = eval(f[tcol])
                    if t <= tlast:
                        dset.append(k, numpy.array(G))
                        G = [eval(f[Gcol])]
                    else:
                        G.append(eval(f[Gcol]))
                    tlast = t
                dset.append(k, numpy.array(G))
            if binsize > 1:
                dset = _gvar.dataset.bin_data(dset, binsize=binsize)
        return dset

    def simulated_pdata_iter(self, n, dataset, p_exact=None, pexact=None, rescale=1.):
        """ Create iterator that returns simulated fit pdata from ``dataset``.

        Simulated fit data has the same covariance matrix as
        ``pdata=self.process_dataset(dataset)``, but mean values that
        fluctuate randomly, from copy to copy, around
        the value of the fitter's fit function evaluated at ``p=p_exact``.
        The fluctuations are generated from bootstrap copies
        of ``dataset``.

        The best-fit results from a fit to such simulated copies of ``pdata``
        should agree with the numbers in ``p_exact`` to within the errors
        specified by the fits (to the simulated data) --- ``p_exact`` gives the
        "correct" values for the parameters. Knowing the correct value for
        each fit parameter ahead of a fit allows us to test the reliability of
        the fit's error estimates and to explore the impact of various fit
        options (*e.g.*, ``fitter.chained_fit`` versus ``fitter.lsqfit``,
        choice of SVD cuts, omission of select models, etc.)

        Typically one need examine only a few simulated fits in order
        to evaluate fit reliability, since we know the correct values
        for the parameters ahead of time. Consequently this method is
        much faster than traditional bootstrap analyses.

        ``p_exact`` is usually taken from the last fit done by the fitter
        (``self.fit.pmean``) unless overridden in the function call.
        Typical usage is as follows::

            dataset = gvar.dataset.Dataset(...)
            data = gvar.dataset.avg_data(dataset)
            ...
            fit = fitter.lsqfit(data=data, ...)
            ...
            for spdata in fitter.simulated_pdata_iter(n=4, dataset):
                # redo fit 4 times with different simulated data each time
                # here p_exact=fit.pmean is set implicitly
                sfit = fitter.lsqfit(pdata=spdata, ...)
                ... check that sfit.p (or functions of it) agrees ...
                ... with p_exact=fit.pmean to within sfit.p's errors      ...

        Args:
            n (int): Maximum number of simulated data sets made by iterator.
            dataset (dictionary): Dataset containing Monte Carlo copies of
                the correlators. ``dataset[datatag]`` is a two-dimensional
                array for the correlator corresponding to ``datatag``,
                where the first index labels the Monte Carlo copy
                and the second index labels time.
            p_exact (dictionary or ``None``): Correct parameter values for
                fits to the simulated data --- fit results should agree
                with ``p_exact`` to within errors. If ``None``, uses
                ``self.fit.pmean`` from the last fit.
            rescale (float): Rescale errors in simulated data by ``rescale``
                (*i.e.*, multiply covariance matrix by ``rescale ** 2``).
                Default is one, which implies no rescaling.
        """
        if pexact is not None:   # for legacy code
            p_exact = pexact
        if p_exact is None:
            if self.fit is None:
                raise ValueError('must specify p_exact')
            p_exact = self.fit.pmean
        pdata = CorrFitter.process_dataset(dataset, self.models)
        pdata_mean = _gvar.mean(pdata)
        pdata_cov = _gvar.evalcov(pdata.buf)
        del pdata
        fcn_mean = self.buildfitfcn()(p_exact)
        if rescale is None or rescale == 1.:
            rescale = None
            correction = _gvar.BufferDict(
                fcn_mean, buf=fcn_mean.buf - pdata_mean.buf
                )
            del fcn_mean
        else:
            pdata_cov *= rescale ** 2
        for bs_dataset in _gvar.dataset.bootstrap_iter(dataset, n=n):
            bs_mean = CorrFitter.process_dataset(
                bs_dataset, self.models, noerror=True
                )
            ans = _gvar.BufferDict()
            if rescale is None:
                for k in bs_mean:
                    ans[k] = bs_mean[k] + correction[k]
            else:
                for k in fcn_mean:
                    ans[k] = (
                        fcn_mean[k] + (bs_mean[k] - pdata_mean[k]) * rescale
                        )
            yield _gvar.BufferDict(
                ans, buf=_gvar.gvar(ans.buf, pdata_cov)
                )


# short cuts
process_data = CorrFitter.process_data
process_dataset = CorrFitter.process_dataset
read_dataset = CorrFitter.read_dataset

class EigenBasis(object):
    """ Eigen-basis of correlator sources/sinks.

    Given :math:`N` sources/sinks and the :math:`N \\times N`
    matrix :math:`G_{ij}(t)` of 2-point correlators created from every
    combination of source and sink, we can define a new basis
    of sources that makes the matrix correlator approximately
    diagonal. Each source in the new basis is associated with
    an eigenvector :math:`v^{(a)}` defined by the matrix equation

    .. math::

        G(t_1) v^{(a)} = \lambda^{(a)}(t_1-t_0) G(t_0) v^{(a)},


    for some :math:`t_0, t_1`. As :math:`t_0, t_1` increase, fewer
    and fewer states couple to :math:`G(t)`. In the limit where
    only :math:`N` states couple, the correlator

    .. math::

        \overline{G}_{ab}(t) \equiv v^{(a)T} G(t) v^{(b)}

    becomes diagonal, and each diagonal element couples to
    only a single state.

    In practice, this condition is only approximate: that is,
    :math:`\overline G(t)`  is approximately diagonal, with diagonal elements
    that overlap strongly with the lowest lying states, but
    somewhat with other states. These new sources are nevertheless useful
    for fits because there is an obvious prior for their amplitudes: ``prior[a][b]``
    approximately equal to one when ``b==a``, approximately zero
    when ``b!=a`` and ``b<N``, and order one otherwise.

    Such a prior can significantly enhance the stability of a multi-source fit,
    making it easier to extract reliable results for excited states.
    It encodes the fact that only a small number of states couple strongly
    to :math:`G(t)` by time :math:`t_0`, without being overly
    prescriptive about what their energies are. We can easilty project our
    correlator onto the new eigen-basis (using :func:`EigenBasis.apply`)
    in order to use this prior, but this is unnecessary.
    :func:`EigenBasis.make_prior` creates a prior of this type in the
    eigen-basis and then transforms it back to the original basis, thereby
    creating an equivalent prior for the amplitudes corresponding
    to the original sources.

    Typical usage is straightforward. For example, ::

        basis = EigenBasis(
            data,                           # data dictionary
            keyfmt='G.{s1}.{s2}',           # key format for dictionary entries
            srcs=['local', 'smeared'],      # names of sources/sinks
            t=(5, 7),                       # t0, t1 used for diagonalization
            )
        prior = basis.make_prior(nterm=4, keyfmt='m.{s1}')

    creates an *eigen-prior* that is optimized for fitting the
    2-by-2 matrix correlator given by ::

        [[data['G.local.local'],     data['G.local.smeared']]
         [data['G.smeared.local'],   data['G.smeared.smeared']]]

    where ``data`` is a dictionary containing all the correlators. Parameter
    ``t`` specifies the times used for the diagonalization: :math:`t_0=5`
    and :math:`t_1=7`. Parameter ``nterm`` specifies the number of terms
    in the fit. ``basis.make_prior(...)`` creates priors
    ``prior['m.local']`` and ``prior['m.smeared']``
    for the amplitudes corresponding to the local and smeared source,
    and a prior ``prior[log(m.dE)]`` for the logarithm of the
    energy differences between successive levels.

    The amplitudes ``prior['m.local']`` and ``prior['m.smeared']`` are
    complicated, with strong correlations between local and smeared entries
    for the same state. Projecting the prior unto the eigen-basis,
    however, reveals its underlying structure::

        p_eig = basis.apply(prior)

    implies ::

        p_eig['m.0'] = [1.0(3), 0.0(1), 0(1), 0(1)]
        p_eig['m.1'] = [0.0(1), 1.0(3), 0(1), 0(1)]

    where the different entries are now uncorrelated.  This structure
    registers our  expectation that the ``'m.0'`` source in the eigen-basis
    overlaps strongly with the  ground state, but almost not at all with the
    first excited state; and vice versa for the ``'m.1'`` source.
    Amplitude ``p_eig`` is noncommittal about  higher states. This structure
    is built into ``prior['m.local']`` and ``prior['smeared']``.

    It is easy to check that fit results are  consistent with the underlying
    prior. This can be done by projecting the  best-fit parameters unto the
    eigen-basis using  ``p_eig = basis.apply(fit.p)``. Alternatively, a
    table listing the  amplitudes in the new eigen-basis, together with the
    energies,  is printed by::

        print(basis.tabulate(fit.p, keyfmt='m.{s1}', eig_srcs=True))

    The prior can be adjusted, if needed, using the ``dEfac``, ``ampl``, and
    ``states`` arguments in :func:`EigenBasis.make_prior`.

    :func:`EigenBasis.tabulate` is also useful for printing the amplitudes
    for the original sources::

        print(basis.tabulate(fit.p, keyfmt='m.{s1}'))

    |EigenBasis| requires the scipy library in Python.

    Args:
        data: Dictionary containing the matrix correlator
            using the original basis of sources and sinks.

        keyfmt: Format string used to generate the keys
            in dictionary ``data`` corresponding to different
            components of the matrix of correlators. The
            key for :math:`G_{ij}` is assumed to be
            ``keyfmt.format(s1=i, s2=j)`` where ``i`` and ``j``
            are drawn from the list of sources, ``srcs``.

        srcs: List of source names used with ``keyfmt``
            to create the keys for finding correlator
            components :math:`G_{ij}` in the data dictionary.

        t: ``t=(t0, t1)`` specifies the ``t`` values
            used to diagonalize the correlation function.
            Larger ``t`` values are better than smaller ones,
            but only if the statistics are adequate.
            When fitting staggered-quark correlators, with oscillating
            components, choose ``t`` values where
            the oscillating pieces are positive (typically odd ``t``).
            If only one ``t`` is given, ``t=t0``, then ``t1=t0+2``
            is used with it. Fits that use |EigenBasis| typically
            depend only weakly on the choice of ``t``.

        tdata: Array containing the times for which there is
            correlator data. ``tdata`` is set equal to
            ``numpy.arange(len(G_ij))`` if it is not specified
            (or equals ``None``).

    The interface for :class:`EigenBasis` is experimental.
    It may change in the near future, as experience
    accumulates from its use.
    """
    def __init__(self, data, srcs, t, keyfmt='{s1}.{s2}', tdata=None, osc=False):
        if keyfmt is None:
            keyfmt = '{s1}.{s2}'
        if scipy is None:
            raise ImportError('need scipy for EigenBasis')
        self.keyfmt = keyfmt
        self.srcs = srcs
        self.eig_srcs = [str(s) for s in numpy.arange(len(self.srcs))]
        self.svdcorrection = _gvar.gvar('0(0)')
        self.svdn = 0
        self.osc = osc
        G = EigenBasis.assemble_data(
            data=data,
            keys=EigenBasis.generate_keys(keyfmt, self.srcs),
            )
        # diagonalize
        if tdata is None:
            tdata = numpy.arange(G.shape[2])
        else:
            tdata = numpy.asarray(tdata)
        self.tdata = tdata
        # find elements in G for diagonalization
        if isinstance(t, tuple):
            i0 = numpy.where(tdata==t[0])[0]
            i1 = numpy.where(tdata==t[1])[0]
            if len(i0) == 0:
                raise ValueError('t[0] not in tdata: ' + str(t[0]))
            if len(i1) == 0:
                raise ValueError('t[1] not in tdata: ' + str(t[1]))
            i0 = i0[0]
            i1 = i1[0]
        else:
            i0 = numpy.where(tdata==t)[0]
            if len(i0) == 0:
                raise ValueError('t not in tdata: ' + str(t))
            i0 = i0[0]
            i1 = i0 + 1 if (i0 + 1) < len(tdata) else i0 - 1
            if i1 < 0:
                raise ValueError('tdata too short: ' + str(tdata))
        t0 = tdata[i0]
        t1 = tdata[i1]
        if self.osc and (t1 - t0) % 2 == 0:
            raise ValueError('t1-t0 must be odd of osc==True')
        self.t = (t0, t1)

        G0 = numpy.empty(G.shape[:-1], float)
        G1 = numpy.empty(G.shape[:-1], float)
        for i in range(G.shape[0]):
            G0[i, i] = G[i, i, i0].mean
            G1[i, i] = G[i, i, i1].mean
            for j in range(i + 1, G.shape[1]):
                G0[i, j] = lsqfit.wavg([G[i, j, i0], G[j, i, i0]], tol=1e-10).mean
                G1[i, j] = lsqfit.wavg([G[i, j, i1], G[j, i, i1]], tol=1e-10).mean
                G0[j, i] = G0[i, j]
                G1[j, i] = G1[i, j]

        if not self.osc:
            w, v = scipy.linalg.eigh(G1, G0, eigvals_only=False)
            v = numpy.array(v.T)
            E = numpy.log(w) / (t0 - t1)
            # set normalization so vn.G.vn = exp(-En*t)
            for i, Ei in enumerate(E):
                v[i] *= _gvar.exp(-E[i] * t0 / 2.)
            E, v = zip(*sorted(zip(E,v)))
            self.v = numpy.array(v)
            self.E = numpy.array(E)
            self.v_inv = numpy.linalg.inv(self.v)
            self.neig = (len(self.E), 0)
        else:
            try:
                wall, vall = scipy.linalg.eigh(G1, G0)
            except numpy.linalg.linalg.LinAlgError:
                try:
                    wall, vall = scipy.linalg.eigh(G0, G1)
                    wall = 1. / wall
                except numpy.linalg.linalg.LinAlgError:
                    raise ValueError('cannot diagonalize G(t)')
            # wall = numpy.array([wi.real for wi in wall])
            vall = numpy.array(vall.T)
            w = wall[wall > 0]
            v = vall[wall > 0]
            E = numpy.log(w) / (t0 - t1)
            wo = wall[wall < 0]
            vo = vall[wall < 0]
            Eo = numpy.log(-wo) / (t0 - t1)
            # set normalization so vn.G.vn = exp(-En*t)
            for i, Ei in enumerate(E):
                v[i] *= numpy.exp(-Ei * t0 / 2.)
            for i, Eoi in enumerate(Eo):
                vo[i] *= numpy.exp(-Eoi * t0 / 2.)
            if len(E) > 1:
                E, v = zip(*sorted(zip(E,v)))
            if len(Eo) > 1:
                Eo, vo = zip(*sorted(zip(Eo,vo)))
            self.v = numpy.array(numpy.concatenate((v, vo)))
            self.E = numpy.array(numpy.concatenate((E, Eo)))
            self.v_inv = numpy.linalg.inv(self.v)
            self.neig = (len(E), len(Eo))

    def make_prior(
        self, nterm,
        keyfmt='{s1}',
        dEfac='1(1)',
        ampl=('1.0(3)', '0.03(10)', '0.2(1.0)'),
        states=None,
        eig_srcs=False,
        ):
        """ Create prior from eigen-basis.

        Args:
            nterm (int): Number of terms in fit function.

            keyfmt (str): Format string usded to generate keys for
                amplitudes and energies in the prior (a dictionary):
                keys are obtained from ``keyfmt.format(s1=a)`` where
                ``a`` is one of the original sources, ``self.srcs``,
                if ``eig_srcs=False`` (default), or one of the
                eigen-sources, ``self.eig_srcs``, if ``eig_srcs=True``.
                The key for the energy differences is generated by
                ``'log({})'.format(keyfmt.format(s1='dE'))``. The default
                is ``keyfmt={s1}``.

            dEfac (str or :class:`gvar.GVar`): A string or :class:`gvar.GVar`
                from which the priors for energy differences ``dE[i]`` are
                constructed. The mean value for ``dE[0]`` is set equal to the
                lowest energy obtained from the diagonalization. The mean
                values for the other ``dE[i]``\s are set equal to the
                difference between the lowest two energies from the
                diagonalization (or to the lowest energy if there is only
                one). These central values are then multiplied by
                ``gvar.gvar(dEfac)``. The default value, `1(1)`, sets the
                width equal to the mean value. The prior is the logarithm of
                the resulting values.

            ampl (tuple): A 3-tuple of strings or :class:`gvar.GVar`\s from
                which priors are contructed for amplitudes corresponding to
                the eigen-sources. ``gvar.gvar(ampl[0])`` is used for
                for source components where the overlap with a particular
                state is expected to be large; ``1.0(3)`` is the default value.
                ``gvar.gvar(ampl[1])`` is used for states that are expected
                to have little overlap with the source; ``0.03(10)`` is
                the default value. ``gvar.gvar(ampl[2])`` is used where
                there is nothing known about the overlap of a state with
                the source; ``0(1)`` is the default value.

            states (list): A list of the states in the correlator corresponding
                to successive eigen-sources, where ``states[i]`` is the state
                corresponding to ``i``-th source. The correspondence between
                sources and states is strong for the first sources, but
                can decay for subsequent sources, depending upon the quality
                of the data being used and the ``t`` values used in the
                diagonalization. In such situations one might specify
                fewer states than there are sources by making the length
                of ``states`` smaller than the number of sources. Setting
                ``states=[]`` assigns broad priors to the every component
                of every source. Parameter ``states`` can also be
                used to deal with situations where the order
                of later sources is not aligned with that of the
                actual states: for example, ``states=[0,1,3]`` connects the
                eigen-sources with the first, second and fourth states in the
                correlator. The default value, ``states=[0, 1 ... N-1]`` where
                ``N`` is the number of sources, assumes that sources
                and states are aligned.

            eig_srcs (bool): Amplitudes for the eigen-sources are
                tabulated if ``eig_srcs=True``; otherwise amplitudes
                for the original basis of sources are tabulated (default).
        """
        def _make_prior(states, keyfmt, E, skip):
            # nsrcs can be larger than nstates
            nsrcs = len(states)
            if numpy.any(states >= nterm):
                n = min(numpy.nonzero(states >= nterm)[0])
                states = states[:n]
            if len(states) > nterm:
                states = states[:nterm]
            one, small, big = ampl
            prior = collections.OrderedDict()
            # amplitudes
            nstates = len(states)
            nterm_ext = nterm + len(skip)
            for i, k in enumerate(EigenBasis.generate_keys(keyfmt, srcs=self.eig_srcs)):
                prior[k] = _gvar.gvar(nterm_ext * [big])
                if True: # i < nsrcs:
                    prior[k][states] = _gvar.gvar(nstates * [small])
                if i < nstates:
                    prior[k][states[i]] = _gvar.gvar(one)
            if not eig_srcs:
                prior = self.unapply(prior, keyfmt=keyfmt)
            if len(skip) > 0:
                idx = numpy.array([i for i in range(nterm_ext) if i not in skip])
                for k in prior:
                    prior[k] = numpy.array(prior[k][idx])
            # energies
            if len(E) == 0:
                E = self.E
            if len(E) >= 2:
                de = E[1] - E[0]
                e0 = E[0]
            else:
                de = E[0]
                e0 = de
            dE = _gvar.gvar(nterm * [str(dEfac)]) * de
            dE[0] = dE[0] + e0 - dE[0].mean
            prior['log({})'.format(keyfmt.format(s1='dE'))] = _gvar.log(dE)
            return prior
        # main body
        if keyfmt is None:
            keyfmt = ('{s1}', 'o.{s1}')
        elif isinstance(keyfmt, str):
            keyfmt = (keyfmt, 'o.{s1}')
        if states is None:
            states = (numpy.arange(self.neig[0]), numpy.arange(self.neig[1]))
        elif isinstance(states, tuple):
            states[0] = numpy.asarray(states[0], int)[:self.neig[0]]
            states[1] = numpy.asarray(states[1], int)[:self.neig[1]]
        else:
            states = numpy.asarray(states, int)
            states = (
                numpy.array(states[:self.neig[0]]),
                numpy.array(states[:self.neig[1]])
                )
        prior = _make_prior(
            states=states[0], keyfmt=keyfmt[0], E=self.E[:self.neig[0]],
            skip=range(self.neig[0], self.neig[0] + self.neig[1]),
            )
        if self.osc:
            oprior =  _make_prior(
                states=states[1] + self.neig[0], keyfmt=keyfmt[1], E=self.E[self.neig[0]:],
                skip=range(0, self.neig[0]),
                )
            prior.update(oprior)
        return prior

    def tabulate(self, p, keyfmt='{s1}', nterm=None, nsrcs=None, eig_srcs=False, indent=4 * ' '):
        """ Create table containing energies and amplitudes for ``nterm`` states.

        Given a correlator-fit result ``fit`` and a corresponding
        :class:`EigenBasis` object ``basis``, a table listing the energies
        and amplitudes for the first ``N`` states in correlators can be printed
        using ::

            print basis.tabulate(fit.p)

        where ``N`` is the number of sources and ``basis`` is an
        :class:`EigenBasis` object. The amplitudes are tabulated for the
        original sources unless parameter ``eig_srcs=True``, in which
        case the amplitudes are projected onto the the eigen-basis
        defined by ``basis``.

        Args:
            p: Dictionary containing parameters values.
            keyfmt: Parameters are ``p[k]`` where keys ``k`` are
                obtained from ``keyfmt.format(s1=s)`` where ``s`` is one
                of the original sources (``basis.srcs``) or one of the
                eigen-sources (``basis.eig_srcs``). The
                default definition is ``'{s1}'``.
            nterm: The number of states from the fit  tabulated.
                The default sets ``nterm`` equal to the number of
                sources in the basis.
            nsrcs: The number of sources tabulated. The default
                causes all sources to be tabulated.
            eig_srcs: Amplitudes for the eigen-sources are
                tabulated if ``eig_srcs=True``; otherwise amplitudes
                for the original basis of sources are tabulated (default).
            indent: A string prepended to each line of the table.
                Default is ``4 * ' '``.
        """

        if keyfmt is None:
            keyfmt = '{s1}'
        if nsrcs is None:
            nsrcs = len(self.srcs)
        if nterm is None:
            nterm = nsrcs
        try:
            dE = p[keyfmt.format(s1='dE')]
        except KeyError:
            logdE = p['log({})'.format(keyfmt.format(s1='dE'))]
            dE = numpy.exp(logdE)
        E = numpy.cumsum(dE)
        if nterm > len(E):
            nterm = len(E)
        if not eig_srcs:
            # tabulate amplitudes for original sources
            srcs = self.srcs[:nsrcs]
            sample_key = keyfmt.format(s1=self.srcs[0])
            if nsrcs > 0 and sample_key not in p:
                p = self.unapply(p, keyfmt=keyfmt)
        else:
            srcs = self.eig_srcs[:nsrcs]
            sample_key = keyfmt.format(s1=self.eig_srcs[0])
            if nsrcs > 0 and sample_key not in p:
                p = self.apply(p, keyfmt=keyfmt)
        keys = EigenBasis.generate_keys(keyfmt=keyfmt, srcs=srcs)
        all_lines = []
        fsize = numpy.zeros(len(keys) + 2, int)
        for l in range(nterm):
            line = [str(l), str(E[l])]
            for k in keys:
                if p[k][l] > 0:
                    line.append(' ' + str(p[k][l]))
                else:
                    line.append(str(p[k][l]))
            for i, f in enumerate(line):
                if len(f) > fsize[i]:
                    fsize[i] = len(f)
            all_lines.append(line)
        fmt = indent
        for fs in fsize:
            fmt += ' {{:{}}}'.format(fs)
        fmt += '\n'
        # header
        line = ['', 'E']
        for s in srcs:
            line.append(' ' + str(s))
        ans = fmt.format(*line)
        ans += indent + ((sum(fsize) + len(fsize)) * '=') + '\n'
        for line in all_lines:
            ans += fmt.format(*line)
        return ans

    def svd(self, data, keyfmt=None, svdcut=1e-15):
        """ Apply SVD cut to data in the eigen-basis.

        The SVD cut is applied to ``data[k]`` where key ``k``
        equals ``keyfmt.format(s1=s1)`` for vector data,
        or ``keyfmt.format(s1=s1, s2=s2)`` for matrix data with sources
        ``s1`` and ``s2`` drawn from ``self.srcs``. The data
        are transformed to the eigen-basis of sources/sinks before
        the cut is applied and then transformed back to the
        original basis of sources. Results are returned in
        a dictionary containing the modified correlators.

        If ``keyfmt`` is a list of formats, the SVD cut is
        applied to the collection of data formed from each
        format. The defaul value for ``keyfmt`` is
        ``self.keyfmt``.
        """
        if keyfmt is None:
            keyfmt = [self.keyfmt]
        elif isinstance(keyfmt, str):
            keyfmt = [keyfmt]
        newdata = self.apply(data, keyfmt=keyfmt)
        newdata = _gvar.svd(newdata, svdcut=svdcut)
        self.svdcorrection = numpy.sum(_gvar.svd.correction)
        self.svdn = _gvar.svd.nmod
        newdata = self.unapply(newdata, keyfmt=keyfmt)
        for i in data:
            if i not in newdata:
                newdata[i] = data[i]
        return newdata

    def apply(self, data, keyfmt='{s1}'):
        """ Transform ``data`` to the eigen-basis.

        The data to be transformed is ``data[k]`` where key ``k``
        equals ``keyfmt.format(s1=s1)`` for vector data,
        or ``keyfmt.format(s1=s1, s2=s2)`` for matrix data with sources
        ``s1`` and ``s2`` drawn from ``self.srcs``.
        A dictionary containing the transformed data is returned
        using the same keys but with the sources replaced by ``'0',
        '1' ...`` (from ``basis.eig_srcs``).

        If ``keyfmt`` is an array of formats, the transformation is
        applied for each format and a dictionary containing all of the
        results is returned. This is useful when the same sources
        and sinks are used for different types of correlators (e.g.,
        in both two-point and three-point correlators).
        """
        if isinstance(keyfmt, str):
            keyfmt = [keyfmt]
        newdata = collections.OrderedDict()
        oldsrcs = self.srcs
        newsrcs = self.eig_srcs
        for fmt in keyfmt:
            newdata.update(EigenBasis._apply(
                data=data, keyfmt=fmt, oldsrcs=oldsrcs, newsrcs=newsrcs,
                v=self.v
                ))
        return newdata

    def unapply(self, data, keyfmt='{s1}'):
        """ Transform ``data`` from the eigen-basis to the original basis.

        The data to be transformed is ``data[k]`` where key ``k``
        equals ``keyfmt.format(s1=s1)`` for vector data,
        or ``keyfmt.format(s1=s1, s2=s2)`` for matrix data with sources
        ``s1`` and ``s2`` drawn from ``self.eig_srcs``. A dictionary
        containing the transformed data is returned using the same keys but
        with the original sources (from ``self.srcs``).

        If ``keyfmt`` is an array of formats, the transformation is
        applied for each format and a dictionary containing all of the
        results is returned. This is useful when the same sources
        and sinks are used for different types of correlators (e.g.,
        in both two-point and three-point correlators).
        """
        if isinstance(keyfmt, str):
            keyfmt = [keyfmt]
        newdata = collections.OrderedDict()
        oldsrcs = self.eig_srcs
        newsrcs = self.srcs
        for fmt in keyfmt:
            newdata.update(EigenBasis._apply(
                data=data, keyfmt=fmt, oldsrcs=oldsrcs, newsrcs=newsrcs,
                v=self.v_inv
                ))
        return newdata

    @staticmethod
    def _apply(data, keyfmt, oldsrcs, newsrcs, v):
        """ Transform data from one basis to another. """
        if len(oldsrcs) != len(newsrcs):
            raise ValueError(
                'wrong number of sources: {o} vs {n}'.format(o=len(oldsrcs), n=len(newsrcs))
                )
        G = EigenBasis.assemble_data(
            data=data,
            keys=EigenBasis.generate_keys(keyfmt, oldsrcs),
            )
        newkeys = EigenBasis.generate_keys(keyfmt, newsrcs)
        newdata = collections.OrderedDict()
        for idx in numpy.ndindex(newkeys.shape):
            tmp = G
            for i in idx:
                tmp = v[i].dot(tmp)
            newdata[newkeys[idx]] = tmp
        return newdata

    @staticmethod
    def generate_keys(keyfmt, srcs):
        """ Generate a list of keys from ``keyfmt`` and ``srcs``. """
        if '{s1}' in keyfmt and '{s2}' in keyfmt:
            nsrc = len(srcs)
            ans = numpy.empty((nsrc, nsrc), object)
            for i, j in numpy.ndindex(nsrc, nsrc):
                ans[i, j] = keyfmt.format(s1=srcs[i], s2=srcs[j])
        elif '{s1}' in keyfmt:
            ans = numpy.array([keyfmt.format(s1=s) for s in srcs], object)
        else:
            raise ValueError('poorly formed keyfmt: ' + str(keyfmt))
        return ans

    @staticmethod
    def assemble_data(data, keys):
        """ Extract array from ``data`` specified by array ``keys``. """
        keys = numpy.asarray(keys, object)
        data_shape = numpy.shape(data[keys.flat[0]])
        data_type = data[keys.flat[0]].dtype
        ans = numpy.empty(keys.shape + data_shape, data_type)
        for idx in numpy.ndindex(keys.shape):
            ans[idx] = data[keys[idx]]
        return ans

class fastfit(object):
    """ Fast fit of a two-point correlator.

    This function class estimates ``E=En[0]`` and ``ampl=an[0]*bn[0]``
    for a two-point correlator modeled by ::

        Gab(t) = sn * sum_i an[i]*bn[i] * fn(En[i], t)
               + so * sum_i ao[i]*bo[i] * fo(Eo[i], t)

    where ``(sn, so)`` is typically ``(1, -1)`` and ::

        fn(E, t) =  exp(-E*t) + exp(-E*(tp-t)) # tp>0 -- periodic
               or   exp(-E*t) - exp(-E*(-tp-t))# tp<0 -- anti-periodic
               or   exp(-E*t)                  # if tp is None (nonperiodic)

        fo(E, t) = (-1)**t * fn(E, t)

    Prior estimates for the amplitudes and energies of excited states are
    used to remove (that is, marginalize) their contributions to give
    a *corrected* correlator ``Gc(t)`` that
    includes uncertainties due to the terms removed. Estimates of ``E``
    are given by::

        Eeff(t) = arccosh(0.5 * (Gc(t+1) + Gc(t-1)) / Gc(t)),

    The final estimate is the weighted average ``Eeff_avg`` of the
    ``Eeff(t)``\s for different ``t``\s. Similarly, an estimate for the
    amplitude ``ampl`` is obtained from the weighted
    average of ::

        Aeff(t) = Gc(t) / fn(Eeff_avg, t).

    If ``osc=True``, an estimate is returned for ``Eo[0]`` rather
    than ``En[0]``, and ``ao[0]*bo[0]`` rather than ``an[0]*bn[0]``.
    These estimates are reliable when ``Eo[0]`` is smaller than
    ``En[0]`` (and so dominates at large ``t``), but probably not
    otherwise.

    Examples:
        The following code examines a periodic correlator (period 64) at large
        times (``t >= tmin``), where estimates for excited states
        don't matter much:

            >>> import corrfitter as cf
            >>> print(G)
            [0.305808(29) 0.079613(24) ... ]
            >>> fit = cf.fastfit(G, tmin=24, tp=64)
            >>> print('E =', fit.E, ' ampl =', fit.ampl)
            E = 0.41618(13)  ampl = 0.047686(95)

        Smaller ``tmin`` values can be used if (somewhat) realistic priors
        are provided for the amplitudes and energy gaps:

            >>> fit = cf.fastfit(G, ampl='0(1)', dE='0.5(5)', tmin=3, tp=64)
            >>> print('E =', fit.E, ' ampl =', fit.ampl)
            E = 0.41624(11)  ampl = 0.047704(71)

        The result here is roughly the same as from the larger ``tmin``, but
        this would not be true for a correlator whose signal to noise ratio
        falls quickly with increasing time.

        :class:`corrfitter.fastfit` estimates the amplitude and energy at
        all times larger than ``tmin`` and then averages to get its final
        results. The chi-squared of the average (*e.g.*, ``fit.E.chi2``)
        gives an indication of the consistency of the estimates from different
        times. The chi-squared per degree of freedom is printed out for both
        the energy and the amplitude using ::

            >>> print(fit)
            E: 0.41624(11) ampl: 0.047704(71) chi2/dof [dof]: 0.9 0.8 [57] Q: 0.8 0.9

        Large values for ``chi2/dof`` indicate an unreliable results. In
        such cases the priors should be adjusted, and/or ``tmin`` increased,
        and/or an SVD cut introduced. The averages in the example above
        have good values for ``chi2/dof``.

    Parameters:
        G: An array of |GVar|\s containing the two-point correlator. ``G[j]``
            is assumed to correspond to time ``t=j``, where ``j=0...``.
        ampl: A |GVar| or its string representation giving an estimate for
            the amplitudes of the ground state and the excited states.  Use
            ``ampl=(ampln, amplo)`` when the  correlator contains oscillating
            states; ``ampln`` is the  estimate for non-oscillating states, and
            ``amplo`` for  oscillating states; setting one or the other to
            ``None``  causes the corresponding terms to be dropped.  Default
            value is ``'0(1)'``.
        dE: A |GVar| or its string representation giving an estimate for the
            energy separation between successive states. This estimate is
            also used to provide an estimate for the lowest energy
            when parameter ``E`` is not specified. Use  ``dE=(dEn, dEo)``
            when the correlator contains oscillating states: ``dEn`` is the
            estimate for non-oscillating states, and ``dEo`` for
            oscillating states; setting one or the other to ``None``
            causes the corresponding terms to be dropped.
            Default value is ``'1(1)'``.
        E: A |GVar| or its string representation giving an estimate for the
            energy of the lowest-lying state. Use  ``E=(En, Eo)``
            when the correlator contains oscillating states: ``En`` is the
            estimate for the lowest non-oscillating state, and ``Eo`` for
            lowest oscillating state. Setting ``E=None`` causes
            ``E`` to be set equal to ``dE``. Default value is ``None``.
        s: A tuple containing overall factors ``(sn, so)`` multiplying
            contributions from the normal and oscillating states.
            Default is ``(1,-1)``.
        tp (int or None): When not ``None``, the correlator is periodic
            with period ``tp`` when ``tp>0``, or anti-periodic with
            period ``-tp`` when ``tp<0``. Setting ``tp=None`` implies
            that the correlator is neither periodic nor anti-periodic.
            Default is ``None``.
        tmin (int): Only ``G(t)`` with ``t >= tmin`` are used. Default
            value is ``6``.
        svdcut (float or None): SVD cut used in the weighted average
            of results from different times. (See the
            :class:`corrfitter.CorrFitter` documentation for a discussion
            of SVD cuts.) Default is ``1e-6``.
        osc (bool): Set ``osc=True`` if the lowest-lying state is an
            oscillating state. Default is ``False``.

    Note that specifying a single |GVar| ``g`` (as opposed to a tuple) for any
    of parameters  ``ampl``, ``dE``, or ``E`` is equivalent to specifying the
    tuple  ``(g, None)`` when ``osc=False``, or the tuple ``(None, g)`` when
    ``osc=True``. A similar rule applies to parameter ``s``.

    :class:`corrfitter.fastfit` objects have the following attributes:

    Attributes:
        E: Energy of the lowest-lying state (|GVar|).
        ampl: Amplitude of the lowest-lying state (|GVar|).

    Both ``E`` and ``ampl`` are obtained by averaging results calculated
    for each time larger than ``tmin``. These are averaged to produce
    a final result. The consistency among results from different times
    is measured by the chi-squared of the average. Each of ``E`` and ``ampl``
    has the following extra attributes:

    Attributes:
        chi2: chi-squared for the weighted average.
        dof: The effective number of degrees of freedom in the weighted
            average.
        Q: The probability that the chi-squared could have been larger,
            by chance, assuming that the data are all Gaussain and consistent
            with each other. Values smaller than 0.05 or 0.1 suggest
            inconsistency. (Also called the *p-factor*.)

    An easy way to inspect these attributes is to print the fit object ``fit``
    using ``print(fit)``, which lists the values of the energy and amplitude,
    the ``chi2/dof`` for each of these, the number of degrees of freedom,
    and the ``Q`` for each.
    """
    def __init__(self, G, ampl='0(1)', dE='1(1)', E=None, s=(1,-1),
        tp=None, tmin=6, svdcut=1e-6, osc=False, nterm=10,
        ):
        import lsqfit
        if not isinstance(s, tuple):
            s = (0, s) if osc else (s, 0)
        def build(x, x0=(None, None)):
            if not isinstance(x, tuple):
                x = (None, x) if osc else (x, None)
            if not isinstance(x0, tuple):
                x0 = (None, x0) if osc else (x0, None)
            return (build_prior(x[0], x0[0]), build_prior(x[1], x0[1]))
        def build_prior(x, x0):
            if x is None:
                return x
            x = _gvar.gvar(x)
            dx = 0 if abs(x.mean) > 0.1 * x.sdev else 0.2 * x.sdev
            xmean = x.mean
            xsdev = x.sdev
            first_x = x if x0 is None else _gvar.gvar(x0)
            return (
                [first_x + dx] +
                [_gvar.gvar(xmean + dx, xsdev) for i in range(nterm - 1)]
                )
        a, ao = build(ampl)
        dE, dEo = build(dE, E)
        s, so = s
        if tp is None:
            def g(E, t):
                return _gvar.exp(-E * t)
            t = numpy.arange(0, len(G))[tmin:]
            G = G[tmin:]
        elif tp > 0:
            def g(E, t):
                return _gvar.exp(-E * t) + _gvar.exp(-E * (tp - t))
            tmax = -tmin + 1 if tmin > 1 else None
            t = numpy.arange(0, len(G))[tmin:tmax]
            G = G[tmin:tmax]
        elif tp < 0:
            def g(E, t):
                return _gvar.exp(-E * t) - _gvar.exp(-E * (-tp - t))
            tmid = int((-tp + 1) // 2)
            G0 = G[0]
            G = numpy.array(
                [G[0]] +
                list(lsqfit.wavg([G[1:tmid], -G[-1:-tmid:-1]], svdcut=svdcut))
                )
            t = numpy.arange(0, len(G))[tmin:]
            G = G[tmin:]
        else:
            raise ValueError('bad tp')
        if len(t) == 0:
            raise ValueError('tmin too large; not t values left')
        if osc:
            G *= (-1) ** t * so
            a, ao = ao, a
            dE, dEo = dEo, dE
            s, so = so, s
        dG = 0.
        E = numpy.cumsum(dE)
        for aj, Ej in list(zip(a, E))[1:]:
            dG += s * aj * g(Ej, t)
        if ao is not None and dEo is not None:
            Eo = numpy.cumsum(dEo)
            for aj, Ej in zip(ao, Eo):
                dG += so * aj * g(Ej, t) * (-1) ** t
        G = G - dG
        ratio = lsqfit.wavg(
            (0.5 * (G[2:] + G[:-2]) / G[1:-1]).tolist() + [_gvar.cosh(E[0])],
            svdcut=svdcut,
            )
        if ratio >= 1:
            self.E = type(ratio)(_gvar.arccosh(ratio), ratio.fit)
        else:
            raise RuntimeError(
                "can't estimate energy: cosh(E) = {}".format(ratio)
                )
        self.ampl = lsqfit.wavg(
            (G / g(self.E, t) / s).tolist() +  [a[0]],
            svdcut=svdcut,
            )

    def __str__(self):
        return (
            "E: {} ampl: {} chi2/dof [dof]: {:.1f} {:.1f} [{}] "
            "Q: {:.1f} {:.1f}"
            ).format(
            self.E, self.ampl, self.E.chi2 / self.E.dof,
            self.ampl.chi2 / self.ampl.dof, self.E.dof, self.E.Q, self.ampl.Q,
            )

