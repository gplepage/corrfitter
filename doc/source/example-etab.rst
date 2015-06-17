Annotated Example: Matrix Correlator
=====================================

.. |etab| replace:: :math:`\eta_b`
.. |CorrFitter| replace:: :class:`corrfitter.CorrFitter`
.. |Corr2| replace:: :class:`corrfitter.Corr2`
.. |EigenBasis| replace:: :class:`corrfitter.EigenBasis`
.. |Dataset| replace:: :class:`gvar.dataset.Dataset`
.. |GVar| replace:: :class:`gvar.GVar`

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
    :lines: 1-25

The eigen-basis is created by ``make_data('etab.data')``:

.. literalinclude:: examples/etab.py
    :lines: 27-33

It reads Monte Carlo
data from file ``'etab.data'``, which has the following format::

    1s0.dd 0.143715 0.0588148 0.030329 ...
    1s0.dd 0.129067 0.0538892 0.0426075 ...
    ...
    1s0.de 0.120838 0.0510525 0.0306193 ...
    1s0.de 0.110661 0.0610127 0.0369064 ...
    ...
    1s0.dg 0.115676 0.0973511 0.0795044 ...
    1s0.dg 0.123485 0.10629 0.0885328  ...
    ...

There are 113 lines in the file for each distinct tag ``1s0.dd``, ``1s0.de`` ..., 
each with 23 numbers. Each line is a separate Monte Carlo estimate 
of the correlator identified by the tag for ``t=1...23``. There are
sixteen different correlators in all, with tags given by::

    '1s0.{s1}{s2}'.format(s1=s1, s2=s2)

where ``s1`` and ``s2`` are  drawn from the source list 
``['l', 'g', 'd', 'e']``.

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
    :lines: 35-47

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
    :lines: 49-50

The prior looks complicated ::

    k            prior[k]
    ------------ ---------------------------------------------------------------
    etab.l       [-0.51(17),  -0.45(16),   0.49(17),   0.103(93), -0.07(87) ...]
    etab.g       [-0.88(27),   0.107(98),  0.051(93), -0.004(90), -0.13(90) ...]
    etab.d       [-0.244(86), -0.50(15),  -0.274(92), -0.166(71), -0.22(60) ...]
    etab.e       [-0.212(84), -0.41(13),  -0.31(10),   0.27(10),  -0.12(62) ...]

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
    :lines: 52-66

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
Running the code produces the following output for the last fit (``N=9``):

.. literalinclude:: examples/etab.out
    :lines: 253-305

This is a good fit, with a chi-squared per degree of freedom of 0.99 for 
260 degrees of freedom (the number of data points fit); the *Q* or *p*-value
is 0.52. This fit required 146 iterations, but took only a few seconds
on a laptop. The results are almost identical to those from ``N=7`` and ``N=8``.

The final energies and amplitudes for the original sources are listed as

.. literalinclude:: examples/etab.out
    :lines: 309-314

while for the eigen-sources they are

.. literalinclude:: examples/etab.out
    :lines: 316-321

The latter shows that the eigen-sources align quite well with the first
four states, as hoped. The errors, especially for the first three states,
are much smaller than the prior errors, which indicates strong signals for 
these states.

Finally values and an error budget are presented for the ``2s-1s`` and 
``3s-1s`` energy differences (in lattice units) and the ratio of the two:

.. literalinclude:: examples/etab.out
    :lines: 323-335

The first excited state is obviously more accurately determined than
the second state, but the fit improves our knowledge of both. 
The energies for the fifth and higher states merely echo the 
*a priori* information in the prior --- the data are not sufficiently 
accurate to add much new information to what was in the prior. The prior
is less important for the three quantities tabulated here. The dominant
source of error in each case comes from either the SVD cut or the statistical
errors in the data.

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
    :lines: 49-50

with 

.. literalinclude:: examples/etab-stab.py
    :lines: 49-50

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
    :lines: 253-

The energies and amplitudes for the first three states are almost unchanged,
which gives us confidence in the original results.
Results for the fourth and higher states have larger errors, as expected.

Note that while the chi-squared value for this last fit is almost identical to
that in the  original fit, the Bayes Factor (from ``logGBF``) is
exp(2171-2158.6)=240,000 times larger for the original fit. The Bayes Factor
gives us a sense of which prior the data prefer. Specifically it says that
our Monte Carlo data are 240,000 times more likely to have come from a model
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

Only five lines in this code differ from the original:
the SVD cut is different (``#1``); the data are projected onto the
eigen-basis (``#2``); the models are defined in terms 
of the eigen-sources (``#3`` and ``#4``); and the prior is defined
for the eigen-sources (``#5``). 

The fit results from the new code are very similar to before; there is 
little difference between the two approaches in this case:

.. literalinclude:: examples/etab-alt.out
    :lines: 253-
