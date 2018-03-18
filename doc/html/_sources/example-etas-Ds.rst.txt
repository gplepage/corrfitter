Annotated Example: Transition Form Factor and Mixing
=====================================================

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
    :lines: 1-24

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
previous plot press "p"; to quit the plots press "q".

We now look at each other major routine in turn.

a) make_data
________________
Method ``make_data('etas-Ds.h5')`` reads in the Monte Carlo data, averages
it, and formats it for use by |CorrFitter|:

.. literalinclude:: examples/etas-Ds.py
    :lines: 51-54

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
We also introduce an SVD cut (see :ref:`svd-cuts`) to account for the fact
that we have only 225 Monte Carlo samples for each piece of data.
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
    :lines: 56-83

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
    :lines: 85-108

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
    :lines: 110-153

The best-fit parameter values are stored in dictionary ``p=fit.p``,
as are the exponentials of the log-normal parameters.
We also turn energy differences into energies using :mod:`numpy`'s cummulative
sum function :func:`numpy.cumsum`. The final output is:

.. literalinclude:: examples/etas-Ds.out
    :lines: 112-123

Finally we  create an error budget for the |etas|
and |Ds| masses, for the mass difference between the |Ds| and its
opposite-parity partner, and for the ground-state transition amplitudes
``Vnn`` and ``Vno``. The quantities of interest are specified in dictionary
``outputs``. For the error budget, we need another dictionary, ``inputs``,
specifying various inputs to the calculation, here the Monte Carlo data and the
priors. Each of these inputs
contributes to the errors in the final results, as detailed in the
error budget:

.. literalinclude:: examples/etas-Ds.out
    :lines: 125-146

The error budget shows, for example, that the largest sources of uncertainty
in every quantity are the statistical errors in the input data.

Results
--------------
The output from running the code is as follows:

.. literalinclude:: examples/etas-Ds.out
    :lines: 1-146

Note:

- This is a relatively simple fit, taking only a couple of seconds on a
  laptop.

- Fits with only one or two terms in the fit function are poor, with
  ``chi2/dof``\s that are significantly larger than one.

- Fits with three terms work well, and adding futher terms has almost no
  impact. The chi-squared does not improve and parameters for the
  added terms differ little from their prior values (since the data are
  not sufficiently accurate to add new information).

- The quality of the fit is confirmed by the fit plots displayed at the
  end (press the 'n' and 'p' keys to cycle through the various plots,
  and the 'q' key to quit the plot). The plot for the |Ds| correlator,
  for example, shows correlator data divided by fit result as a
  function of ``t``:

  .. image:: examples/Ds.*
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

Variation: Marginalization
--------------------------
Marginalization (see :ref:`marginalized-fits`) can speed up fits like
this one. To use an 8-term fit function, while tuning parameters for only
``N`` terms, we change only four lines in the main program:

.. literalinclude:: examples/etas-Ds-marginalize.py
    :lines: 14-28

The first modification (``#1``)
limits the
fits to ``N=1,2``, because that is all that will be needed to get good
values for the leading term.
The second modification (``#2``) sets the prior to eight terms, no matter what value
``N`` has. The third (``#3``) tells ``fitter.lsqfit`` to fit parameters from
only the first ``N`` terms in the fit function; parts of the prior that are
not being fit are incorporated (*marginalized*) into the fit data.
The last modification (``#4``) changes what is printed out.
The output
shows that
results for the leading term have converged by ``N=2`` (and even ``N=1`` is
pretty good):

.. literalinclude:: examples/etas-Ds-marginalize.out
    :lines: 1-101


Variation: Chained Fit
------------------------
Chained fits (see :ref:`chained-fits`) are used if ``fitter.lsqfit(...)``
is replaced by ``fitter.chained_lsqfit(...)`` in ``main()``. The results
are about the same: for example,

.. literalinclude:: examples/etas-Ds-chained.out
    :lines: 353-387

Chained fits are particularly useful for very large data sets
(much larger than this one).

Test the Analysis
---------------------
We can test our analysis by adding
``test_fit(fitter, 'etas-Ds.h5')`` to the ``main``
program, where:

.. literalinclude:: examples/etas-Ds.py
    :lines: 27-49

This code does ``n=2`` simulations of the full fit, using the means of fit
results from the last fit done by ``fitter`` as ``pexact``.
The code compares fit results woth ``pexact`` in each case,
and computes the chi-squared of the difference between the leading
parameters and ``pexact``. The output is:

.. literalinclude:: examples/etas-Ds.out
    :lines: 151-

This shows that the fit is working well.

Other options are easily checked. For example,
only one line need be changed in ``test_fit`` in order to test
a marginalized fit::

    sfit = fitter.lsqfit(pdata=spdata, prior=prior, p0=pexact, nterm=(2,2))

Running this code gives:

.. literalinclude:: examples/etas-Ds-marginalize.out
    :lines: 106-

This is also fine and confirms that ``nterm=(2,2)`` marginalized fits
are a useful, faster substitute for full fits in this case.

Mixing
----------
Code to analyze |Ds|-|Ds| mixing is very similar to the code above for
a transition form factor. The ``main()`` and ``make_data()`` functions are
identical, except that here data are read from file ``'Ds-Ds.h5'`` and
the appropriate SVD cut is ``svdcut=0.014`` (see :ref:`svd-cuts`).
We need models for the two-point |Ds| correlator, and for two three-point
correlators describing the |Ds| to |Ds| transition:

.. literalinclude:: examples/Ds-Ds.py
    :lines: 33-55

The initial and final states in the three-point correlators are the same
here so we set parameter ``symmetricV=True`` in :class:`corrfitter.Corr3`.

The prior is also similar to the previous case:

.. literalinclude:: examples/Ds-Ds.py
    :lines: 57-76

We use log-normal distributions for the energy differences, and
amplitudes.
We store only the upper triangular parts of the ``Vnn`` and ``Voo`` matrices
since they are symmetrical (because ``symmetricV=True`` is set).

A minimal ``print_results()`` function is:

.. literalinclude:: examples/Ds-Ds.py
    :lines: 78-94

Running the mixing code gives the following output:

.. literalinclude:: examples/Ds-Ds.out
    :lines: 1-

The fits for individual correlators look good:

==============================  =================================== ===================================
==============================  =================================== ===================================
.. image:: examples/Ds-Ds.Ds.*  .. image:: examples/Ds-Ds.DsDsT15.* .. image:: examples/Ds-Ds.DsDsT18.*
==============================  =================================== ===================================

