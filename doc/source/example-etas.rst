Annotated Example: Two-Point Correlator
============================================

.. |etas| replace:: :math:`\eta_s`
.. |Ds| replace:: :math:`D_s`
.. |CorrFitter| replace:: :class:`corrfitter.CorrFitter`
.. |Corr2| replace:: :class:`corrfitter.Corr2`
.. |Corr3| replace:: :class:`corrfitter.Corr3`
.. |Dataset| replace:: :class:`gvar.dataset.Dataset`
.. |GVar| replace:: :class:`gvar.GVar`

.. |~| unicode:: U+00A0
   :trim:

Introduction
-------------
The simplest use of :mod:`corrfitter` is calculating the
amplitude and energy of the ground state in a single
two-point correlator. Here we analyze an |etas| propagator where the
source and sink are the same.

The one slightly non-obvious aspect of this fit is its use of
log-normal priors for the energy differences ``dE`` between successive
states in the correlator. As discussed in :ref:`positive-parameters`,
this choice imposes an order on the states in relation to the
fit parameters by forcing all ``dE`` values to be positive.
Any such restriction helps stabilize a fit, improving both
efficiency and the final results.

Another design option that  helps stabilize the fit is to do a series of
fits, with increasing  number ``N`` of states in the fit function, where
the results  from the ``N-1`` fit are used by the fitter as the  starting
point (``p0``) for the ``N`` fit. The initial fits are bad, but this procedure
helps guide the fit parameters towards sensible values as the number of states
increases. See :ref:`faster-fits` for more discussion.

The source code (``etas.py``) and
data file (``etas-Ds.data``) are included with the :mod:`corrfitter`
distribution, in the ``examples/`` directory.
The data are from the HPQCD collaboration.


Code
------------------
Following the template outlined in :ref:`basic-fits`, the entire
code is:

.. literalinclude:: examples/etas.py
    :lines: 1-


Here the Monte Carlo data are read by ``make_data('etas.data')``  from file
``etas.data``. This file contains 225 lines,  each
with 64 numbers, of the form::

    etas    0.305044    0.0789607   0.0331313 ...
    etas    0.306573    0.0802435   0.0340765 ...
    ...

Each line is a different Monte Carlo estimate of the |etas|
correlator for t=0...63. The mean values and covariance matrix
are computed for the 64 elements of the correlator using
:func:`gvar.dataset.avg_data`, and the result is stored in
``data['etas']``, which is an array of Gaussian random
variables (objects of type :class:`gvar.GVar`).

A |CorrFitter| object, ``fitter``, is created for a single two-point
correlator from a list of models created by ``make_models()``. There
is only one model in the list because there is only one correlator.
It is a :class:`Corr2` object which specifies that:
the key (``datatag``) for extracting the correlator from the data dictionary is
``'etas'``; the propagator is periodic with period 64; each correlator
contains data for t values ranging from 0 to 63; only values
greater than or equal to 5 and less than 64-5 are fit; the
source and sink amplitudes are the same and labeled by ``'a'`` in
the prior; and the energy differences between successive states
are labeled ``'dE'`` in the prior.

Fits are tried with ``N`` states in the fit function, where ``N``
varies from 2 to 5. Usually ``N=2`` is too small, resulting in
a poor fit. Here we will find that results have converged by
``N=3``.

A prior, containing *a priori* estimates for the fit parameters,
is contructed for each ``N`` by ``make_prior(N)``. The amplitude priors,
``prior['a'][i]``,
are assumed to be 0±1, while the differences between successive
energies are taken to be, roughly, 0.5±0.5. These are broad priors,
based upon preliminary fits of the data.
We want to use log-normal
statistics for the energy differences, to guarantee that
they are positive (and the states ordered, in order of
increasing energy), so we use ``prior['logdE']`` for
the logarithms of the differences --- instead of
``prior['dE']`` for the differences themselves --- and take the logarithm
of the prior.

The fit is done by ``fitter.lsqfit()`` and
``print_results(fit)`` prints results for the first two states after
each fit (that is, for each ``N``). Note how results from
the fit to ``N`` terms is used as the starting point for the
fit with ``N+1`` terms, via parameter ``p0``. As mentioned above,
this speeds up the larger fits and also helps to stabilize them.

Results
--------
The output from this fit code is:

.. literalinclude:: examples/etas.out
    :lines: 1-

These fits are very fast --- a small fraction of a second each on a laptop.
Fit results converge by ``N=3`` states. The amplitudes and energy differences
for states above the first three are essentially identical to the prior
values; the Monte Carlo data are not sufficiently accurate to add any new
information about these levels. The fits for ``N>=3`` are excellent, with
chi-square per degree of freedom (``chi2/dof``) of 0.68. There are only
28 degrees of freedom here because the fitter, taking advantage of the
periodicity, folded the data about the midpoint in ``t`` and averaged,
before fitting. The ground state energy and amplitude are determined to
a part in 1,000 or better.

Correlated Data?
----------------------

It is worth checking whether the initial Monte Carlo data has correlations
from sample to sample, since such correlations lead to
underestimated fit errors. One approach is verify that results are
unchanged when the input data are binned. To bin the data we use ::

    def make_data(filename):
        """ Read data, compute averages/covariance matrix for G(t). """
        return gv.dataset.avg_data(cf.read_dataset(filename, binsize=2))

which averages successive samples (bins of 2). Binned data give the following
results from the last iteration and summary:

.. literalinclude:: examples/etas-binned.out
    :lines: 42-65

These agree pretty well with the previous results, suggesting that
correlations are not a problem.

Binning should have no significant effect on results if there are
no correlations, provided the total number of samples after binning
is sufficiently large (*e.g.*, more than 100—200). Strong
correlations cause error estimates to grow with increased
bin size (like the square root of ``binsize``). Binning reduces
correlations; data should be binned
with increasing bin sizes until fit error estimates stop growing.

Fast Fit and Effective Mass
----------------------------

The last two lines in the ``main()`` function of the code illustrate the  use
of :class:`corrfitter.fastfit` to get a very fast results for the lowest-energy
state. As discussed in :ref:`very-fast-fits`, :class:`corrfitter.fastfit`
provides  an alternative to the multi-exponential fits discussed above when
only  the lowest-energy parameters are needed. The method used is similar to a
traditional effective mass analysis except that estimates for contributions
from excited states are generated from priors and removed from the correlator
before determining the effective mass. This allows the code to use much
smaller ``t`` values than in the traditional approach, thereby obtaining
results that rival the multi-exponential fits.

In this example, :class:`corrfitter.fastfit` is used to analyze the two-point
correlator stored in array ``data['etas']``. The amplitudes for different
states are estimated to have size 0±1, while the spacings between energies
(and  between the first state and 0) are estimated to be 0.5±0.5. The code
averages results form all ``t`` values down to ``tmin=3``. Setting ``tp=64``
indicates that the correlator is periodic with period 64.

The last line of the output summarizes the results of the fast fit. The
energy and amplitude are almost identical to what was obtained from the
multi-exponential fits (note that ``fastfit.ampl`` is the same  as
``fit.a[0]**2``, which has value 0.047681(79)). :class:`corrfitter.fastfit`
estimates the energy and amplitude for each ``t`` greater than ``tmin``, and
then  averages the results. The consistency of results from different ``t``\s
is  measured by the chi-squared of the averages. The chi-squared per degree of
freedom is reported here to be 0.8 for the ``E`` average and 0.9 for the
``ampl`` average, indicating that there is good agreement between different
``t``\s.

While a fast fit is easier to set up, multi-exponential fits are  usually more
robust, and provide more detailed information about the fit. One  use for fast
fits is to estimate the sizes of parameters for use in designing  the priors
for a multi-exponential fit. There are often situations where *a priori*
knowledge about fit parameters is sketchy, especially for amplitudes. A  fast
fit to data at large ``t`` can quickly generate estimates for both  amplitudes
and energies, from which it is then easy to construct priors. In the code
above, for example, we could replace ``make_prior(N)``
by ``alt_make_prior(N, data['etas'])`` where::

    def alt_make_prior(N, G):
        fastfit = cf.fastfit(G=G, tmin=24, tp=64)
        da = 2 * fastfit.ampl.mean ** 0.5
        dE = 2 * fastfit.E.mean
        prior = collections.OrderedDict()
        prior['a'] = [gv.gvar(0, da) for i in range(N)]
        prior['log(dE)'] = gv.log([gv.gvar(dE, dE) for i in range(N)])
        return prior

This code does a fast fit using data from very large ``t``, where priors for
the excited states are unimportant. It then uses the results to create priors
for the amplitudes and energy differences for all states, assuming that the
ground state values are either larger, or smaller by no more than roughly a
factor of two. This customized prior gives results that are almost identical
to what was obtained using the original prior, above (in part because
the original prior is pretty sensible to begin with).

Designing a prior using :class:`corrfitter.fastfit` would be even more
useful when multiple sources and sinks are involved, as in a matrix fit.

