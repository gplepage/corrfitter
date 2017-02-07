Annotated Example: Matrix Correlator
=====================================

.. |etab| replace:: :math:`\eta_b`
.. |CorrFitter| replace:: :class:`corrfitter.CorrFitter`
.. |Corr2| replace:: :class:`corrfitter.Corr2`
.. |EigenBasis| replace:: :class:`corrfitter.EigenBasis`
.. |Dataset| replace:: :class:`gvar.dataset.Dataset`
.. |GVar| replace:: :class:`gvar.GVar`
.. |~| unicode:: U+00A0
   :trim:


Introduction
------------
Matrix correlators, built from multiple sources and sinks,
greatly improve results for the excited states in the correlators.
Here we analyze |etab| correlators using 4 sources/sinks using
a prior designed by |EigenBasis|.

A major challenge when studying excited states with multi-source
fits is the appearance in the fit of spurious states, with amplitudes that
are essentially zero, between the real states in the correlator.
These states contribute little to the correlators, because of their
vanishing amplitudes, but they usually have a strong negative impact
on the errors of states just below them and above them. They also can cause
the fitter to stall, taking 1000s of iterations to change nothing other
than the parameters of the spurious state.
|EigenBasis| addresses this problem by creating a prior that discourages
spurious states. It encodes the fact that only a small number of states
still couple to the matrix correlator by moderate values of ``t``, and
therefore, that there exist linear combinations of the sources that
couple strongly to individual low-lying states but not the others. This
leaves little room for spurious low-lying states.

This example involves only two-point correlators. Priors generated
by |EigenBasis| are also useful in fits to three-point correlators,
with multiple eigen-bases if different types of hadron are involved.

The source code (``etab.py``) and data file
(``etab.data``) are included with the :mod:`corrfitter` distribution,
in the ``examples/`` directory.
The data are from the HPQCD collaboration.

Code
-----------------
The main method follows the template in :ref:`basic-fits`, but
modified to handle the |EigenBasis| object ``basis``:

.. literalinclude:: examples/etab.py
    :lines: 1-8, 13-26

The eigen-basis is created by ``make_data('etab.h5')``:

.. literalinclude:: examples/etab.py
    :lines: 28-34

It reads Monte Carlo
data from file ``'etab.h5'``, which is an hdf5-format file that
contains 16 hdf5 datasets of interest::

    >>> import h5py
    >>> dset = h5py.File('etab.h5', 'r')
    >>> for v in dset.values():
    ...     print(v)
    <HDF5 dataset "1s0.dd": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.de": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.dg": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.dl": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.ed": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.ee": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.eg": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.el": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.gd": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.ge": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.gg": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.gl": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.ld": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.le": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.lg": shape (113, 23), type "<f8">
    <HDF5 dataset "1s0.ll": shape (113, 23), type "<f8">
    ...

Each of these contains 113 Monte Carlo samples for a
different correlator evaluated at times ``t=1,2...23``::

    >>> print(dset['1s0.ll'][:, :])
    [[ 0.360641  0.202182  0.134458 ...,  0.00118724  0.00091401  0.00070451]
     [ 0.365291  0.210573  0.143632 ...,  0.00133125  0.0010258   0.00078879]
     [ 0.362848  0.210732  0.143037 ...,  0.00130843  0.00101     0.00077627]
     ...,
     [ 0.364053  0.209284  0.141544 ...,  0.00120762  0.00093206  0.0007196 ]
     [ 0.365236  0.209835  0.139057 ...,  0.00116917  0.00089636  0.00069045]
     [ 0.362479  0.20709   0.136687 ...,  0.00106393  0.0008269   0.00064707]]

The sixteen different correlators have tags given by::

    '1s0.{s1}{s2}'.format(s1=s1, s2=s2)

where ``s1`` and ``s2`` are  drawn from the list
``['l', 'g', 'd', 'e']`` which labels the sources and sinks used
to create the correlators.

The data are read in, and their means and covariance matrix computed
using :func:`gvar.dataset.avg_data`. |EigenBasis| then creates an
eigen-basis by solving a generalized eigenvalue problem involving the
matrices of correlators at ``t=1`` and ``t=2``. (One might want
larger ``t`` values generally, but these data are too noisy.)
The eigenanalysis constructs a set of eigen-sources
that are linear combinations of the original sources chosen
so that each eigen-source overlaps strongly with one of the
lowest four states in the correlator, and weakly with all the others.
This eigen-basis is used later to construct the prior.

A correlator fitter, called ``fitter``, is created from the list of correlator
models returned by ``make_models(basis)``:

.. literalinclude:: examples/etab.py
    :lines: 36-48

There is one model for each correlator to be fit, so 16 in all. The keys
(``datatag``) for the correlator data are constructed from information
stored in the ``basis``. Each correlator
has data (``tdata``) for ``t=1...23``.
We fit all ``t`` values (``tfit``) for the diagonal elements of the
matrix correlator, but only about half the ``t`` values for other
correlators --- the information at large ``t`` is highly correlated between
different correlators, and therefore somewhat redundant. The
amplitudes are labeled by ``'etab.l'``, ``'etab.g'``, ``'etab.d'``, and
``'etab.e'`` in the prior. The energy differences are labeled by ``'etab.dE'``.

We try fits with ``N=1,2..9`` terms in the fit function. The number of
terms is encoded in the prior, which is constructed by
``make_prior(N, basis)``:

.. literalinclude:: examples/etab.py
    :lines: 50-51

The prior looks complicated ::

    k            prior[k]
    ------------ ---------------------------------------------------------------
    etab.l       [-0.51(17),  -0.45(16),   0.50(17),   0.108(94), -0.06(87) ...]
    etab.g       [-0.88(27),   0.104(97),  0.051(92), -0.004(90), -0.13(90) ...]
    etab.d       [-0.243(86), -0.50(15),  -0.268(91), -0.161(71), -0.21(60) ...]
    etab.e       [-0.211(83), -0.42(13),  -0.30(10),   0.26(10),  -0.12(61) ...]

    log(etab.dE) [-1.3(2.3),  -0.5(1.0),  -0.5(1.0),  -0.5(1.0),  -0.5(1.0) ...]



but its underlying structure becomes clear if we project it unto the
eigen-basis using ``prior_eig = basis.apply(prior, keyfmt='etab.{s1}')``::

    k         prior_eig[k]
    ------    ------------------------------------------------------------
    etab.0    [1.00(30),    .03(10),   .03(10),   .03(10),   .2(1.0) ... ]
    etag.1    [ .03(10),   1.00(30),   .03(10),   .03(10),   .2(1.0) ... ]
    etab.2    [ .03(10),    .03(10),  1.00(30),   .03(10),   .2(1.0) ... ]
    etab.3    [ .03(10),    .03(10),   .03(10),  1.00(30),   .2(1.0) ... ]

The *a priori* expectation built into ``prior_eig`` (and therefore ``prior``) is
that the ground state  overlaps strongly with the first source in the
eigen-basis, and weakly  with the other three. Similarly the first excited state
overlaps strongly with the second eigen-source, but none of the others. And so
on. The fifth  and higher excited states can overlap with every eigen-source.
The priors for the energy differences between successive levels are based
upon the energies obtained from the eigenanalysis (``basis.E``): the ``dE``
prior for the ground state is taken to be ``E0(E0)``, where ``E0=basis.E[0]``,
while for the other states it equals ``dE1(dE1)``,  where
``dE1=basis.E[1]-basis.E[0]``. The prior specifies log-normal  statistics for
``etab.dE`` and so replaces it by ``log(etab.dE)``.

The fit is done by ``fitter.lsqfit(...)``. An SVD cut is needed
(``svdcut=0.0004``) because the data are highly correlated.  The data are also
over-binned, to keep down the size of the  :mod:`corrfitter` distribution, and
this mandates an SVD cut, as well. Fit stability is sometimes improved by
applying the SVD cut to the data in the  eigen-basis rather than in the
original source basis. Here this could be done by replacing the last line of
``make_data(...)`` with::

        return basis.svd(data, svdcut=5e-3), basis

and dropping the ``svdcut=0.0004`` from ``fitter.lsqfit(...)``. The size
of the ``svdcut`` is usually different.


Final results are printed out by ``print_results(...)``
after the last fit is finished:

.. literalinclude:: examples/etab.py
    :lines: 53-67

This method first writes out two tables listing energies and amplitudes for
the first 4 states in the correlator. The first table shows results for the
original sources, while the second is for the eigen-sources. The correlators
are from NRQCD so only energy differences are physical. The energy differences
for each of the first two excited states relative to the ground states are
stored in dictionary ``outputs``. These are in lattice units. ``outputs``
also contains the ratio of ``3s-1s`` difference to the ``2s-1s`` difference,
and here the lattice spacing cancels out. The code automatically handles
statistical correlations between different energies as it does the arithmetic
for ``outputs`` --- the fit results are all |GVar|\s. The ``outputs``
are tabulated using :func:`gvar.fmt_values`. An error budget is also
produced, using :func:`gvar.fmt_errorbudget`, showing how much error
for each quantity comes from uncertainties in the prior and data, and
from uncertainties introduced by the SVD cut.

Finally plots showing the data divided by the fit for each correlator are
displayed (optionally).

Results
----------
Running the code produces the following output for the last fit (``N=7``):

.. literalinclude:: examples/etab.out
    :lines: 43-93

This is a good fit, with a chi-squared per degree of freedom of 1.0 for
260 degrees of freedom (the number of data points fit); the *Q* or *p*-value
is 0.51. This fit required 97 iterations, but took only a few seconds
on a laptop. The results are almost identical to those from ``N=6`` and ``N=8``.

The final energies and amplitudes for the original sources are listed as

.. literalinclude:: examples/etab.out
    :lines: 97-102

while for the eigen-sources they are

.. literalinclude:: examples/etab.out
    :lines: 104-109

The latter shows that the eigen-sources align quite well with the first
four states, as hoped. The errors, especially for the first three states,
are much smaller than the prior errors, which indicates strong signals for
these states.

Finally values and an error budget are presented for the ``2s-1s`` and
``3s-1s`` energy differences (in lattice units) and the ratio of the two:

.. literalinclude:: examples/etab.out
    :lines: 111-123

The first excited state is obviously more accurately determined than
the second state, but the fit improves our knowledge of both.
The energies for the fifth and higher states merely echo the
*a priori* information in the prior --- the data are not sufficiently
accurate to add much new information to what was in the prior. The prior
is less important for the three quantities tabulated here. The dominant
source of error in each case comes from either the SVD cut or the statistical
errors in the data.

Summary plots showing the data divided by the fit as a function of ``t``
for each of the 16 |~| correlators is shown below:

=================================   =================================   =================================   =================================
=================================   =================================   =================================   =================================
.. image:: examples/etab.1s0.ll.*   .. image:: examples/etab.1s0.lg.*   .. image:: examples/etab.1s0.ld.*   .. image:: examples/etab.1s0.le.*
.. image:: examples/etab.1s0.gl.*   .. image:: examples/etab.1s0.gg.*   .. image:: examples/etab.1s0.gd.*   .. image:: examples/etab.1s0.ge.*
.. image:: examples/etab.1s0.dl.*   .. image:: examples/etab.1s0.dg.*   .. image:: examples/etab.1s0.dd.*   .. image:: examples/etab.1s0.de.*
.. image:: examples/etab.1s0.el.*   .. image:: examples/etab.1s0.eg.*   .. image:: examples/etab.1s0.ed.*   .. image:: examples/etab.1s0.ee.*
=================================   =================================   =================================   =================================

These plots are displayed by the code above if flag
``DISPLAYPLOTS = True`` is set at the beginning of the code. The points with
error bars are correlator data points; the fit result is 1.0 in these
plots; and the dashed lines show the uncertainty in the fit function
values for the best-fit parameters. Fit and data agree well for
all correlators and all ``t`` values. As expected, strong correlations
exist between points with near-by ``t``\s.

Fit Stability
---------------
It is a good idea in fits like this one to test the stability of the
results to significant changes in the prior. This is especially true for
quantities like the ``3s-1s`` splitting that involve more highly excited
states. The default prior in effect assigns each of the four sources in the
new basis to one of the four states in the correlator with the lowest
energies. Typically the actual correspondence between source and low-energy
state weakens as the energy increases. So an obvious test is to rerun the
fit but with a prior that associates states with only three of the sources,
leaving the fourth source unconstrained. This is done by replacing

.. literalinclude:: examples/etab.py
    :lines: 50-51

with

.. literalinclude:: examples/etab-stab.py
    :lines: 50-51

in the code. The ``states`` option in the second ``basis.make_prior(...)``
assigns the three lowest lying states (in order of increasing energy)
to the first three eigen-sources, but leaves the fourth and higher states
unassigned. The prior for the amplitudes projected onto the eigen-basis
then becomes ::


    k         prior_eig[k]
    ------    ----------------------------------------------------------
    etab.0    [1.00(30),    .03(10),   .03(10),  .2(1.0),  .2(1.0) ... ]
    etab.1    [ .03(10),   1.00(30),   .03(10),  .2(1.0),  .2(1.0) ... ]
    etab.2    [ .03(10),    .03(10),  1.00(30),  .2(1.0),  .2(1.0) ... ]
    etab.3    [ .2(1.0),    .2(1.0),   .2(1.0),  .2(1.0),  .2(1.0) ... ]


where now no strong assumption is made about the overlaps of the first three
eigen-sources with the fourth state, or about the overlap of the fourth source
with any state. Running with this (more conservative) prior gives
the following results for the last fit and summary:

.. literalinclude:: examples/etab-stab.out
    :lines: 43-123

The energies and amplitudes for the first three states are almost unchanged,
which gives us confidence in the original results.
Results for the fourth and higher states have larger errors, as expected.

Note that while the chi-squared value for this last fit is almost identical to
that in the  original fit, the Bayes Factor (from ``logGBF``) is
exp(2175.1-2160.9)=1,469,000 times larger for the original fit. The Bayes Factor
gives us a sense of which prior the data prefer. Specifically it says that
our Monte Carlo data are 1,469,000 times more likely to have come from a model
with the original prior than from one with the more conservative prior. This
further reinforces our confidence in the original results.

Alternative Organization
-------------------------

Our |etab| fit uses the |EigenBasis| to construct a special prior for the fit,
but leaves the correlators unchanged. An alternative approach is  to project
the correlators unto the eigen-basis, and then to fit them with a fit function
defined directly in terms of the eigen-basis. This approach is conceptually
identical to that above, but in practice it gives somewhat different
results since the SVD cut  enters differently. A version of the code above
that uses this approach is:

.. literalinclude:: examples/etab-alt.py
    :lines: 1-8, 13-

Only five lines in this code differ from the original:
the SVD cut is different (``#1``); the data are projected onto the
eigen-basis (``#2``); the models are defined in terms
of the eigen-sources (``#3`` and ``#4``); and the prior is defined
for the eigen-sources (``#5``).

The fit results from the new code are very similar to before; there is
little difference between the two approaches in this case:

.. literalinclude:: examples/etab-alt.out
    :lines: 43-123
