Annotated Example: Transition Form Factor
============================================

.. |etas| replace:: :math:`\eta_s`
.. |Ds| replace:: :math:`D_s`
.. |CorrFitter| replace:: :class:`corrfitter.CorrFitter`
.. |Corr2| replace:: :class:`corrfitter.Corr2`
.. |Corr3| replace:: :class:`corrfitter.Corr3`
.. |Dataset| replace:: :class:`gvar.dataset.Dataset`
.. |GVar| replace:: :class:`gvar.GVar`

Introduction
--------------
Here we describe a complete Python code that uses :mod:`corrfitter` 
to calculate the transition matrix element or form factor from 
an |etas| meson to a |Ds| meson, together with the masses and amplitudes
of these mesons. A very similar code could be used to calculate 
mixing amplitudes, such as for *B* mixing.

This example combines data from two-point correlators, for the
amplitudes and energies, with data from three-point correlators, 
for the transition matrix element. We fit all of the correlators 
together, in a single fit, in order to capture correlations between
the various output parameters. The correlations are built into
the output parameters and consequently are reflected in any 
arithmetic combination of parameters --- no bootstrap is needed to 
calculate correlations or their impact on quantities derived 
from the fit parameters. The best-fit parameters (in ``fit.p`` and
``fit.transformed_p``) are objects of type |GVar|.

Staggered quarks are used in 
this simulation, so the |Ds| has oscillating components as well as 
normal components in its correlators.

The source code (``etas-Ds.py``) and data file 
(``etas-Ds.data``) are included with the :mod:`corrfitter` distribution,
in the ``examples/`` directory. The data are from the HPQCD collaboration.

Code
-----------------
The ``main`` method for this code follows the pattern described 
in :ref:`basic-fits`:

.. literalinclude:: examples/etas-Ds.py
    :lines: 1-25, 158-

The raw Monte Carlo data is in a file named ``'etas-Ds.data'``. We are doing
four fits, with 1, 2, 3, and 4 terms in the fit function. Each fit starts its
minimization at point ``p0``, which is set equal to the mean values of the
best-fit parameters from the previous fit (``p0 = fit.pmean``). This reduces
the  number of iterations needed for convergence in the ``N = 4`` fit, for
example, from 217 to 24. It also makes multi-term fits more stable.

The last line of ``main()`` 
displays plots of the fit data divided by the fit, provided 
:mod:`matplotlib` is installed. A plot is made for each correlator, and the
ratios should equal one to within errors. To 
move from one plot to the next press "n" on the keyboard; to move to a 
previous plot press "p"; to quit the plots press "q".

We now look at each other major routine in turn.

a) make_data
________________
Method ``make_data('etas-Ds.data')`` reads in the Monte Carlo data, averages
it, and formats it for use by |CorrFitter|. The data file (``'eta-Ds.data'``)
contains 225 lines,  each with 64 numbers on it, of the form::

    etas    0.305044    0.0789607   0.0331313 ...
    etas    0.306573    0.0802435   0.0340765 ...
    ...

Each of these lines is a single Monte Carlo estimate for the |etas| 
correlator on a lattice with 64 lattice points in the ``t`` direction; 
there are 225 Monte Carlo estimates in all. The same file also contains
225 lines describing the |Ds| meson correlator::

    Ds      0.230503    0.0445531   0.00895383 ...
    Ds      0.230947    0.0447479   0.00904294 ...
    ...

And it contains 225 lines each giving the 3-point amplitude for 
:math:`\eta_s \to D_s`
where the source and sink are separated by 15 and 16 time steps on the
lattice::

    3ptT15  4.63494e-10     1.11333e-09     2.46993e-09 ...
    3ptT15  4.85637e-10     1.15445e-09     2.59419e-09 ...
    ...

    3ptT16  1.42457e-10     3.27314e-10     7.61508e-10 ...
    3ptT16  1.47582e-10     3.4255e-10      7.95205e-10 ...
    ...

The first, second, third, *etc.* lines for each label come from the first, 
second, third, *etc.* Monte Carlo iterations, respectively; this 
synchronization allows
the code to compute correlations between different 
types of data.

Function :func:`corrfitter.read_dataset` is designed to read files 
in this format (among others). We use it to read the data, and
:func:`gvar.dataset.avg_data` to compute the means and covariance 
matrix of the data:

.. literalinclude:: examples/etas-Ds.py
    :lines: 49-51

This routine returns a dictionary whose keys are the strings used to label the
individual lines in ``etas-Ds.data``: for example, ::

    >>> data = make_data('etas-Ds.data')
    >>> print(data['Ds'])
    [0.2307150(73) 0.0446523(32) 0.0089923(15) ... 0.0446527(32)]
    >>> print(data['3ptT16'])
    [1.4583(21)e-10 3.3639(44)e-10 ... 0.000023155(30)]

Here each entry in ``data`` is an array of |GVar|\s representing the Monte
Carlo estimates (mean and covariance) for the corresponding correlator. This 
is the format needed by |CorrFitter|.

b) make_models
__________________
Method ``make_models()`` specifies the theoretical models that will be used 
to fit the data:

.. literalinclude:: examples/etas-Ds.py
    :lines: 53-85

Four models are specified, one for each correlator to be fit. The first two
are for the |etas| and |Ds| two-point correlators, corresponding to
entries in the data dictionary with keys ``'etas'`` and ``'Ds'``,
respectively. 
These are periodic propagators, with period 64 (``tp``), and we want to
omit the first and last 5 (``tmin``) time steps in the correlator. The
``t``\s to be fit are listed in ``tfit``, while the ``t``\s contained in the
data are in ``tdata``. Labels for the fit parameters corresponding to the
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
parameter, to be used as priors in the fitter:

.. literalinclude:: examples/etas-Ds.py
    :lines: 87-110

Parameter ``N`` specifies how many terms are kept in the fit functions. The
priors are specified in a dictionary ``prior``. Each entry is an array, of
length ``N``, with one entry for each term. Each entry is a Gaussian random
variable, specified by an object of type  |GVar|. Here we use the fact that
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
    :lines: 112-156

The best-fit parameter values are stored in dictionary ``p=fit.transformed_p``,
as are the exponentials of the log-normal parameters.
We also turn energy differences into energies using :mod:`numpy`'s cummulative
sum function :func:`numpy.cumsum`. The final output is:

.. literalinclude:: examples/etas-Ds.out
    :lines: 157-168

Finally we  create an error budget for the |etas|
and |Ds| masses, for the mass difference between the |Ds| and its
opposite-parity partner, and for the ground-state transition amplitudes
``Vnn`` and ``Vno``. The quantities of interest are specified in dictionary
``outputs``. For the error budget, we need another dictionary, ``inputs``,
specifying various inputs to the calculation: the Monte Carlo data, the
priors, and the results from any *svd* cuts (none here). Each of these inputs
contributes to the errors in the final results, as detailed in the
error budget:

.. literalinclude:: examples/etas-Ds.out
    :lines: 170-191

The error budget shows, for example, that the largest sources of uncertainty
in every quantity are the statistical errors in the input data. 

Results
--------------
The output from running the code is as follows:

.. literalinclude:: examples/etas-Ds.out
    :lines: 1-191

Note:

- This is a relatively simple fit, taking only a couple of seconds on a 
  laptop.

- Fits with only one or two terms in the fit function are poor, with 
  ``chi2/dof``\s that are significantly larger than one.

- Fits with three terms work well, and adding futher terms has almost no 
  impact. The ``chi**2`` does not improve and parameters for the
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
  result is 1.0 in this plot, of course, and the dashed lines show the 
  uncertainty in the fit function evaluated with the best-fit parameters.
  Fit and data agree to within errors. Note how the fit-function errors
  (the dashed lines) track the data errors. In general the fit function
  is at least as accurate as the data. It can be much more accurate, 
  for example, when the data errors grow rapidly with ``t``.


Variation: Marginalization
--------------------------
Marginalization (see :ref:`marginalized-fits`) can speed up fits like 
this one. To use an 8-term fit function, while tuning parameters for only
``N`` terms, we change only four lines in the main program:

.. literalinclude:: examples/etas-Ds-marginalize.py
    :lines: 14-27

The first modification (``#1``) is in the definition of ``fitter``, where we add
an extra argument to tell |CorrFitter| what kind of marginalization 
to use (that is, not the ratio method). The second modification (``#2``) 
limits the
fits to ``N=1,2``, because that is all that will be needed to get good
values for the leading term.
The third modification (``#3``) sets the prior to eight terms, no matter what value
``N`` has. The last (``#4``) tells ``fitter.lsqfit`` to fit parameters from 
only the first ``N`` terms in the fit function; parts of the prior that are
not being fit are incorporated (*marginalized*) into the fit data. The output 
shows that
results for the leading term have converged by ``N=2`` (and even ``N=1`` isn't
so bad):

.. literalinclude:: examples/etas-Ds-marginalize.out
    :lines: 1-81


Variation: Chained Fit
------------------------
Chained fits (see :ref:`chained-fits`) are used if ``fitter.lsqfit(...)`` 
is replaced by ``fitter.chained_lsqfit(...)`` in ``main()``. The results
are about the same: for example,

.. literalinclude:: examples/etas-Ds-chained.out
    :lines: 170-175

We obtain more or less the same results,

.. literalinclude:: examples/etas-Ds-chained.out
    :lines: 273-278


if we polish the final results from the chained fit using 
a final call to ``fitter.lsqfit`` (see :ref:`chained-fits`):

.. literalinclude:: examples/etas-Ds-chained.py
    :lines: 30

Another variation is to replace the last line (``return models``) 
in ``make_models()`` by::

    return [models[:2]] + models[2:]

This causes the two 2-point correlators (``models[:2]``) to be fit
in parallel, which makes sense since they share no parameters.
The result of the (parallel) fit of the 2-point correlators is used
as a prior for the chained fits of the 3-point correlators (``models[2:]``).
The fit results are mostly unchanged, although the polishing fit 
is significantly faster (more than 2x) in this case:

.. literalinclude:: examples/etas-Ds-chained.out
    :lines: 441-446

Test the Analysis
---------------------
We can test our analysis by adding 
``test_fit(fitter, 'etas-Ds.data')`` to the ``main`` 
program, where:

.. literalinclude:: examples/etas-Ds.py
    :lines: 28-47

This code does ``n=2`` simulations of the full fit, using the means of fit 
results from the last fit done by ``fitter`` as ``pexact``. 
The code prints out each fit,
and for each it computes the ``chi**2`` of the difference between the leading
parameters and ``pexact``. The output is:

.. literalinclude:: examples/etas-Ds.out
    :lines: 196-

This shows that the fit is working well, at least for the leading 
parameter for each key. 

Other options are easily checked. For example,
only one line need be changed in ``test_fit`` in order to test 
a marginalized fit::

    sfit = fitter.lsqfit(data=sdata, prior=prior, p0=pexact, nterm=(2,2))

Running this code gives:

.. literalinclude:: examples/etas-Ds-marginalize.out
    :lines: 86-

This is also fine and confirms that ``nterm=(2,2)`` marginalized fits 
are a useful, faster substitute for full fits. Indeed the simulation 
suggests that the marginalized fit is somewhat more accurate
than the original fit for the oscillating-state parameters (``Vno``, 
``log(Dso:a)``, ``log(Dso:dE)`` --- compare the simulated results with
the ``nterm=4`` results from the original fit, as these were used to 
define ``pexact``).