Annotated Example: Transition Form Factor and Mixing
=====================================================

.. |etas| replace:: :math:`\eta_s`
.. |Ds| replace:: :math:`D_s`
.. |CorrFitter| replace:: :class:`corrfitter.CorrFitter`
.. |Corr2| replace:: :class:`corrfitter.Corr2`
.. |Corr3| replace:: :class:`corrfitter.Corr3`
.. |Dataset| replace:: :class:`gvar.dataset.Dataset`
.. |GVar| replace:: :class:`gvar.GVar`
.. |chi2| replace:: :math:`\chi^2`
.. |~| unicode:: U+00A0
   :trim:

Introduction
--------------
Here we describe a complete Python code that uses :mod:`corrfitter`
to calculate the transition matrix element or form factor from
an |etas| meson to a |Ds| meson, together with the masses and amplitudes
of these mesons. A very similar code, for (speculative) |Ds|-|Ds| mixing,
is described at the end.

The form factor example combines data from two-point correlators, for the
amplitudes and energies, with data from three-point correlators,
for the transition matrix element. We fit all of the correlators
together, in a single fit, in order to capture correlations between
the various output parameters. The correlations are built into
the output parameters and consequently are reflected in any
arithmetic combination of parameters --- no bootstrap is needed to
calculate correlations or their impact on quantities derived
from the fit parameters. The best-fit parameters (in ``fit.p``)
are objects of type |GVar|.

Staggered quarks are used in
this simulation, so the |Ds| has oscillating components as well as
normal components in its correlators.

The source codes (``etas-Ds.py``, ``Ds-Ds.py``) and data files
(``etas-Ds.h5``, ``Ds-Ds.h5``) are included with
the :mod:`corrfitter` distribution, in the ``examples/`` directory.
The data are from the HPQCD collaboration.

Code
-----------------
The ``main`` method for the form-factor code follows the pattern described
in :ref:`basic-fits`:

.. literalinclude:: examples/etas-Ds.py
    :lines: 1-23

The Monte Carlo data are in a file named ``'etas-Ds.h5'``. We are doing
four fits, with 1, 2, 3, and 4 terms in the fit function. Each fit starts its
minimization at point ``p0``, which is set equal to the mean values of the
best-fit parameters from the previous fit (``p0 = fit.pmean``). This reduces
the  number of iterations needed for convergence in the ``N = 4`` fit, for
example, from 162 to 45. It also makes multi-term fits more stable.

After the fit, plots of the fit data divided by the fit
are displayed by ``fit.show_plots()``, provided
:mod:`matplotlib` is installed. A plot is made for each correlator, and the
ratios should equal one to within errors. To
move from one plot to the next press "n" on the keyboard; to move to a
previous plot press "p"; to cycle through different views of the data and fit
press "v"; and to quit the plots press "q".

We now look at each other major routine in turn.

a) make_data
________________
Method ``make_data('etas-Ds.h5')`` reads in the Monte Carlo data, averages
it, and formats it for use by |CorrFitter|:

.. literalinclude:: examples/etas-Ds.py
    :lines: 60-63

The data file ``etas-Ds.h5`` is in hdf5 format. It contains four datasets::

  >>> for v in dset.values():
  ...     print(v)
  <HDF5 dataset "3ptT15": shape (225, 16), type "<f8">
  <HDF5 dataset "3ptT16": shape (225, 17), type "<f8">
  <HDF5 dataset "Ds": shape (225, 64), type "<f8">
  <HDF5 dataset "etas": shape (225, 64), type "<f8">

Each corresponds to Monte Carlo data for a single correlator, which
is packaged as a two-dimensional :mod:`numpy` array whose
first index labels the Monte Carlo sample, and whose second index labels
time. For example, ::

  >>> print(dset['etas'][:, :])
  [[ 0.305044   0.0789607  0.0331313 ...,  0.0164646  0.0332153  0.0791385]
   [ 0.306573   0.0802435  0.0340765 ...,  0.0170088  0.034013   0.0801528]
   [ 0.306194   0.0800234  0.0338007 ...,  0.0168862  0.0337728  0.0799462]
   ...,
   [ 0.305955   0.0797565  0.0335741 ...,  0.0167847  0.0336077  0.0796961]
   [ 0.305661   0.0793606  0.0333133 ...,  0.0165365  0.0333934  0.0792943]
   [ 0.305365   0.079379   0.033445  ...,  0.0164506  0.0332284  0.0792884]]

is data for a two-point correlator describing the |etas|  meson. Each
of the 225 lines is a different Monte Carlo sample for the correlator, and
has 64 entries corresponding to ``t=0,1...63``. Note the periodicity
in this data.

Function ``gv.dataset.avg_data(dset)`` averages over the Monte Carlo samples
for all the correlators to compute their means and covariance matrix.
The end result is a dictionary whose keys are the keys used to
label the hdf5 datasets: for example, ::

    >>> data = make_data('etas-Ds.h5')
    >>> print(data['etas'])
    [0.305808(29) 0.079613(24) 0.033539(17) ... 0.079621(24)]
    >>> print(data['Ds'])
    [0.2307150(73) 0.0446523(32) 0.0089923(15) ... 0.0446527(32)]
    >>> print(data['3ptT16'])
    [1.4583(21)e-10 3.3639(44)e-10 ... 0.000023155(30)]

Here each entry in ``data`` is an array of |GVar|\s representing Monte
Carlo averages for the corresponding correlator at different times. This
is the format needed by |CorrFitter|. Note that the different correlators
are correlated with each other: for example, ::

    >>> print(gv.evalcorr([data['etas'][0], data['Ds'][0]]))
    [[ 1.          0.96432174]
     [ 0.96432174  1.        ]]

shows a 96% correlation between the ``t=0`` values in the |etas| and
|Ds| correlators.

b) make_models
__________________
Method ``make_models()`` specifies the theoretical models that will be used
to fit the data:

.. literalinclude:: examples/etas-Ds.py
    :lines: 65-92

Four models are specified, one for each correlator to be fit. The first two
are for the |etas| and |Ds| two-point correlators, corresponding to
entries in the data dictionary with keys ``'etas'`` and ``'Ds'``,
respectively.
These are periodic propagators, with period 64 (``tp``), and we want to
omit the first and last 5 (``tmin``) time steps in the correlator.
Labels for the fit parameters corresponding to the
sources (and sinks) are specified for each, ``'etas:a'`` and ``'Ds:a'``, as
are labels for the energy differences, ``'etas:dE'`` and ``'Ds:dE'``.  The
|Ds| propagator also has an oscillating piece because this data comes from
a staggered-quark analysis. Sources/sinks and energy differences are
specified for these as well: ``'Dso:a'`` and ``'Dso:dE'``.

Finally three-point models are specified for the data corresponding to
data-dictionary keys ``'3ptT15'`` and ``'3ptT16'``. These share several
parameters with the two-point correlators, but introduce new parameters
for the transition matrix elements: ``'Vnn'`` connecting normal states, and
``'Vno'`` connecting normal states with oscillating states.

c) make_prior
_________________
Method ``make_prior(N)`` creates *a priori* estimates for each fit
parameter, to be used as  priors in the fitter:

.. literalinclude:: examples/etas-Ds.py
    :lines: 94-117

Parameter ``N`` specifies how many terms are kept in the fit functions. The
priors are stored in a dictionary ``prior``. Each entry is an array, of
length ``N``, with one entry for each term in the fit function.
Each entry is a Gaussian random
variable, an object of type  |GVar|. Here we use the fact that
``gvar.gvar()`` can make a list of |GVar|\s from a list of strings of the form
``'0.1(1)'``: for example, ::

    >>> print(gv.gvar(['1(2)', '3(2)']))
    [1.0(2.0) 3.0(2.0)]

In this particular fit, we can assume that all the sinks/sources
are positive, and we can require that the energy differences be positive. To
force positivity, we use log-normal distributions for these parameters by
defining priors for ``'log(etas:a)'``, ``'log(etas:dE)'`` ... rather than
``'etas:a'``,  ``'etas:dE'`` ... (see :ref:`positive-parameters`). The *a
priori* values for these fit parameters are the logarithms of the values for
the parameters themselves: for example, each ``'etas:a'`` has prior ``0.3(3)``,
while the actual fit parameters, ``log(etas:a)``, have priors
``log(0.3(3)) = -1.2(1.0)``.

We override the default priors for the ground-state energies in each case.
This is not unusual since ``dE[0]``, unlike the other ``dE``\s,  is an energy,
not an energy difference. For the oscillating |Ds| state, we require
that its mass be ``0.3(3)`` larger than the |Ds| mass. One could  put
more precise information into the priors if that made sense given the goals
of the simulation. For example, if the main objective is a value for ``Vnn``,
one might include fairly exact information about the |Ds| and
|etas| masses in the prior, using results from experiment or from
earlier simulations. This would make no sense, however, if the goal is to
verify that simulations gives correct masses.

Note, finally, that a statement like ::

    prior['Vnn'] = gv.gvar(N * [N* ['0(1)']])       # correct

is *not* the same as ::

    prior['Vnn'] = N * [N * [gv.gvar('0(1)')]]      # wrong

The former creates ``N ** 2`` independent |GVar|\s, with one for each element
of ``Vnn``; it is one of the most succinct ways of creating a large number of
|GVar|\s. The latter creates only a single |GVar| and uses it repeatedly for
every element ``Vnn``, thereby forcing every element of ``Vnn``  to be equal
to every other element when fitting (since the difference between any two of
their priors is ``0±0``); it is almost certainly not what is desired.
Usually one wants to create the array of strings first, and then convert it to
|GVar|\s using ``gvar.gvar()``.

d) print_results
_____________________
Method ``print_results(fit, prior, data)`` reports on the best-fit values
for the fit parameters from the last fit:

.. literalinclude:: examples/etas-Ds.py
    :lines: 119-159

The best-fit parameter values are stored in dictionary ``p=fit.p``,
as are the exponentials of the log-normal parameters.
We also turn energy differences into energies using :mod:`numpy`'s cummulative
sum function :func:`numpy.cumsum`. The final output is:

.. literalinclude:: examples/etas-Ds.out
    :lines: 112-122

Finally we  create an error budget for the |etas|
and |Ds| masses, and for the ground-state transition amplitude
``Vnn``. The quantities of interest are specified in dictionary
``outputs``. For the error budget, we need another dictionary, ``inputs``,
specifying various inputs to the calculation, here the Monte Carlo data,
the SVD corrections, and the
priors. Each of these inputs
contributes to the errors in the  final results, as detailed in the
error budget:

.. literalinclude:: examples/etas-Ds.out
    :lines: 124-143

The error budget shows, for example, that the largest sources of uncertainty
in every quantity are the statistical errors in the input data.

e) SVD Cut
_____________________
The fits need an SVD cut, ``svdcut=8e-5``, because there are
only 225 |~| random samples for
each of 69 |~| data points. Following the recipe in :ref:`svd-cuts`,
we determine the value of the SVD cut using a separate
script consisting of

.. literalinclude:: examples/etas-Ds-svdcut.py
    :lines: 3-10

together with the ``make_models()`` method described above.
Running this script outputs a single line, telling us
to use ``svdcut=8e-5``::

    svdcut = 7.897481075407619e-05

It also displays a plot that shows how eigenvalues
of the data's correlation matrix that are below the SVD cutoff (dotted
red line) are significantly underestimated:

.. image:: examples/etas-Ds-svdcut.png
    :width: 70%

The SVD cut adds additional uncertainty to the data to increase
these eigenvalues so they do not cause trouble with the fit.

Results
--------------
The output from running the code is as follows:

.. literalinclude:: examples/etas-Ds.out
    :lines: 1-143

Note:

- This is a relatively simple fit, taking only a second or so on a
  laptop.

- Fits with only one or two terms in the fit function are poor, with
  ``chi2/dof``\s that are significantly larger than one.

- Fits with three terms work well, and adding futher terms has almost no
  impact. The |chi2| does not improve and parameters for the
  added terms differ little from their prior values (since the data are
  not sufficiently accurate to add new information).

- The quality of the fit is confirmed by the fit plots displayed at the
  end (press the 'n' and 'p' keys to cycle through the various plots,
  the 'v' key to cycle through different views of the data and fit,
  and the 'q' key to quit the plot). The plot for the |Ds| correlator,
  for example, shows correlator data divided by fit result as a
  function of ``t``:

  .. image:: examples/Ds.png
      :width: 70%

  The points with error bars are the correlator data points; the fit
  result is 1.0 in this plot, of course, and the shaded band shows the
  uncertainty in the fit function evaluated with the best-fit parameters.
  Fit and data agree to within errors. Note how the fit-function errors
  (the shaded band) track the data errors. In general the fit function
  is at least as accurate as the data. It can be much more accurate,
  for example, when the data errors grow rapidly with ``t``.

- In many applications precision can be improved by factors of 2—3 |~| or more
  by using multiple sources and sinks for the correlators. The code here
  is easily generalized to handle such a situation: each
  :class:`corrfitter.Corr2` and :class:`corrfitter.Corr3` in ``make_models()``
  is replicated with various different combinations of sources and
  sinks (one entry for each combination).

Testing the Fit
----------------
The |chi2| for the fit above is low for a fit to 69 |~| data points.
This would normally be evidence of a good fit, but, as
discussed in :ref:`goodness-of-fit`,
our SVD cut means that we need
to refit with added noise if we want to use |chi2| as a measure of
fit quality.
We implement this test by adding the following code at the end of the
``main()`` method  above:

.. literalinclude:: examples/etas-Ds.py
    :lines: 25-36, 51-58

This reruns the fit but with random noise (associated with the SVD cut
and priors) added to the data. The result is still a good fit, but with
a signficantly higher |chi2|, as expected:

.. literalinclude:: examples/etas-Ds.out
    :lines: 148-157

The noisy fit also agrees with the original fit about the most important
parameters from the fit. This test provides further evidence that
our fit is good.

A more complex test of the fitting protocol is obtained by using
simulated fits: see :ref:`simulated-fits`. We do this by adding

.. literalinclude:: examples/etas-Ds.py
    :lines: 38-50

to the end of the ``main()`` method. This code does ``n=2`` |~| simulations
of the full fit, using the means ``fit.pmean`` from the last fit
as ``p_exact``.
The code compares fit results with ``p_exact`` in each case,
and computes the |chi2| of the difference between the leading parameters
and ``p_exact``. The output is:

.. literalinclude:: examples/etas-Ds.out
    :lines: 159-177

This again confirms that the fit is working well.


Variation: Marginalization
--------------------------
Marginalization (see :ref:`marginalization`) can speed up fits like
this one. To use an 8-term fit function, while tuning parameters for only
``N`` terms, we change only four lines in the main program:

.. literalinclude:: examples/etas-Ds-marginalize.py
    :lines: 12-26

The first modification (``#1``) sets the prior to eight terms,
no matter what value ``N`` has.
The second modification (``#2``)
limits the
fits to ``N=1,2``, because that is all that will be needed to get good
values for the leading term.
The third (``#3``) tells ``fitter.lsqfit`` to fit parameters from
only the first ``N`` terms in the fit function; parts of the prior that are
not being fit are incorporated (*marginalized*) into the fit data.
The last modification (``#4``) changes what is printed out.
The output
shows that
results for the leading term have converged by ``N=2`` (and even ``N=1`` is
pretty good):

.. literalinclude:: examples/etas-Ds-marginalize.out
    :lines: 1-98

The tests applied to the first fit can be used here as well. For example,
setting ``add_svdnoise=True`` and ``add_priornoise=True`` in the
fit results in

.. literalinclude:: examples/etas-Ds-marginalize.out
    :lines: 103-112

This suggests a good fit. The results are consistent with the original
fit.


Variation: Chained Fit
------------------------
Chained fits are used if ``fitter.lsqfit(...)``
is replaced by ``fitter.chained_lsqfit(...)`` in ``main()``. Following
the advice at the end of :ref:`chained-fits`, we combine chained
fits with marginalization. Three parts of our original code need
modifications:

.. literalinclude:: examples/etas-Ds-chained.py
    :lines: 14-34

The first modification (``#1``) replaces the original list of models with
a structured list that instructs the (chained) fitter sequentially to:

  a) fit the ``etas`` 2-point correlator described in ``models[0]`` (``#1b``);

  b) fit the ``Ds`` 2-point correlator described in ``models[1]``  (``#1b``);

  c) reset fit parameters ``nterm=(2, 1)`` (marginalize all but the
     first three states) and ``svdcut=6.3e-5`` for
     subsequent fits (``#1c``);

  d) fit simultaneously the two 3-point correlators described
     in ``(models[2],models[3])``  (``#1d``).

The second modification (``#2``) replaces ``lsqfit`` by ``chained_lsqfit``.
It also removes the SVD cutfrom the 2-point fits;
as discussed above (``1c``), the SVD cut is reintroduced for the 3-point
fits.

The output for ``N=4`` terms is substantially shorter than for our
original code:

.. literalinclude:: examples/etas-Ds-chained.out
    :lines: 52-153

Again the results agree well with the original fit.
Fit results are listed from each step in the chain: first just
the ``etas`` 2-point correlator, then the ``Ds`` 2-point
correlator, and finally a combined fit of both 3-point
correlators. One might try less marginalization (e.g.,
``nterm=(1,1)``) to check that results are stable. Also one might
test the fits, as above.

Chained fits are particularly useful for very large data sets
(much larger than this one). Also marginalizing extraneous variables in the
3-point fits can make fitting more robust (because it is simpler).


Mixing
----------
Code to analyze |Ds|-|Ds| mixing is very similar to the code above for
a transition form factor. The ``main()`` and ``make_data()`` functions are
identical, except that here data are read from file ``'Ds-Ds.h5'`` and
the appropriate SVD cut is ``svdcut=0.003`` (see :ref:`svd-cuts`).
We need models for the two-point |Ds| correlator, and for two three-point
correlators describing the |Ds| to |Ds| transition:

.. literalinclude:: examples/Ds-Ds.py
    :lines: 68-90

The initial and final states in the three-point correlators are the same
here so we set parameter ``symmetricV=True`` in :class:`corrfitter.Corr3`.

The prior is also similar to the previous case:

.. literalinclude:: examples/Ds-Ds.py
    :lines: 92-111

We use log-normal distributions for the energy differences, and
amplitudes.
We store only the upper triangular parts of the ``Vnn`` and ``Voo`` matrices
since they are symmetrical (because ``symmetricV=True`` is set).

A minimal ``print_results()`` function is:

.. literalinclude:: examples/Ds-Ds.py
    :lines: 113-129

Running the mixing code gives the following output for ``N=4``:

.. literalinclude:: examples/Ds-Ds.out
    :lines: 22-112

The fits for individual correlators look good:

================================  ===================================== =====================================
================================  ===================================== =====================================
.. image:: examples/Ds-Ds.Ds.png  .. image:: examples/Ds-Ds.DsDsT15.png .. image:: examples/Ds-Ds.DsDsT18.png
================================  ===================================== =====================================

