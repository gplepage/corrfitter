:mod:`corrfitter` - Least-Squares Fit to Correlators
=====================================================

.. moduleauthor:: G.P. Lepage <g.p.lepage@cornell.edu>

.. |CorrFitter| replace:: :class:`corrfitter.CorrFitter`
.. |Corr2| replace:: :class:`corrfitter.Corr2`
.. |Corr3| replace:: :class:`corrfitter.Corr3`
.. |BaseModel| replace:: :class:`corrfitter.BaseModel`
.. |Dataset| replace:: :class:`gvar.dataset.Dataset`
.. |GVar| replace:: :class:`gvar.GVar`

Introduction 
------------------
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
    
        
Basic Fits
----------
To illustrate, consider data for two 2-point correlators: ``Gaa`` with the
same source and sink (``a``), and ``Gab`` which has source ``a`` and
(different) sink ``b``. The data are contained in a dictionary ``data``,
where ``data['Gaa']`` and ``data['Gab']`` are one-dimensional arrays
containing values for ``Gaa(t)`` and ``Gab(t)``, respectively, with
``t=0,1,2...63``. Each array element in ``data['Gaa']`` and ``data['Gab']``
is a gaussian random variable of type |GVar|, and specifies the mean and
standard deviation for the corresponding data point::
        
    >>> print data['Gaa']
    [0.159791 +- 4.13311e-06 0.0542088 +- 3.06973e-06 ... ]
    >>> print data['Gab']
    [0.156145 +- 1.83572e-05 0.102335 +- 1.5199e-05 ... ]
        
|GVar|\s can also capture any statistical correlations between different
pieces of data.
        
We want to fit this data to the following formulas::
        
    Gaa(t,N) = sum_i=0..N-1  a[i]**2 * exp(-E[i]*t)
    Gab(t,N) = sum_i=0..N-1  a[i]*b[i] * exp(-E[i]*t)
        
Our goal is to find values for the amplitudes, ``a[i]`` and ``b[i]``, and the
energies, ``E[i]``, so that these formulas reproduce the average values for
``Gaa(t,N)`` and ``Gab(t,N)`` that come from the data, to within the data's
statistical errors. We use the same ``a[i]``\s and ``E[i]``\s in both
formulas. The fit parameters used by the fitter are the ``a[i]``\s and
``b[i]``\s, as well as the differences ``dE[i]=E[i]-E[i-1]`` for ``i>0`` and
``dE[0]=E[0]``. The energy differences are usually positive by construction
(see below) and are easily converted back to energies using::
        
    E[i] = sum_j=0..i dE[j]
        
A typical code has the following structure::
        
    from corrfitter import CorrFitter
        
    data = make_data('mcfile')          # user-supplied routine
    models = make_models()              # user-supplied routine
    N = 4                               # number of terms in fit functions
    prior = make_prior(N)               # user-supplied routine
    fitter = CorrFitter(models=models)
    fit = fitter.lsqfit(data=data, prior=prior)  # do the fit
    print_results(fit, prior, data)     # user-supplied routine
        
We discuss each user-supplied routine in turn.
        
``make_data('mcfile')``
____________________________
``make_data('mcfile')`` creates the dictionary containing the data that is to
be fit. Typically such data comes from a Monte Carlo simulation. Imagine that
the simulation creates a file called ``'mcfile'`` with layout ::
        
     # first correlator: each line has Gaa(t) for t=0,1,2...63
     Gaa  0.159774739530e+00 0.541793561501e-01 ...
     Gaa  0.159751906801e+00 0.542054488624e-01 ...
     Gaa  ...
     .
     .
     .
     # second correlator: each line has Gab(t) for t=0,1,2...63
     Gab  0.155764170032e+00 0.102268808986e+00 ...
     Gab  0.156248435021e+00 0.102341455176e+00 ...
     Gab  ...
     .
     .
     .
        
where each line is one Monte Carlo measurement for one or the other
correlator, as indicated by the tags at the start of each line. (Lines for
``Gab`` may be interspersed with lines for ``Gaa`` since every line has a
tag.) The data can be analyzed using the :mod:`gvar.dataset` module::
        
    import gvar
        
    def make_data(filename):
        dset = gvar.dataset.Dataset(filename)
        return gvar.dataset.avg_data(dset)
        
This reads the data from file into a dataset object (type
:class:`gvar.dataset.Dataset`) and then computes averages for each
correlator and ``t``, together with a covariance matrix for the set of
averages. Thus ``data = make_data('mcfile')`` creates a dictionary where
``data['Gaa']`` is a 1-d array of |GVar|\s obtained by averaging over the
``Gaa`` data in the ``'mcfile'``, and ``data['Gab']`` is a similar array
for the ``Gab`` correlator.
        
``make_models()``
____________________
``make_models()`` identifies which correlators in the fit data are to be fit,
and specifies theoretical models (that is, fit functions) for these
correlators::
        
    from corrfitter import Corr2
        
    def make_models():
        models = [ Corr2(datatag='Gaa', tdata=range(64), tfit=range(64),
                        a='a', b='a', dE='dE'),
                        
                   Corr2(datatag='Gab', tdata=range(64), tfit=range(64),
                        a='a', b='b', dE='dE')
                 ]
        return models
        
For each correlator, we specify: the key used in the input data dictionary
``data`` for that correlator (``datatag``); the values of ``t`` for which
results are given in the input data (``tdata``); the values of ``t`` to
keep for fits (``tfit``, here the same as the range in the input data, but
could be any subset); and fit-parameter labels for the source (``a``) and
sink (``b``) amplitudes, and for the intermediate energy-differences
(``dE``). Fit-parameter labels identify the parts of the prior,
discussed below, corresponding to the actual fit parameters (the labels are
dictionary keys). Here the two models, for ``Gaa`` and ``Gab``, are
identical except for the data tags and the sinks. ``make_models()`` returns
a list of models; the only parts of the input fit data that are fit are
those for which a model is specified in ``make_models()``.
        
Note that if there is data for ``Gba(t,N)`` in addition to ``Gab(t,N)``, and
``Gba = Gab``, then the (weighted) average of the two data sets will be
fit if ``models[1]`` is replace by::
        
    Corr2(datatag='Gab', tdata=range(64), tfit=range(64),
         a=('a', None), b=('b', None), dE=('dE', None),
         othertags=['Gba'])
        
The additional argument ``othertags`` lists other data tags that correspond
to the same physical quantity; the data for all equivalent data tags is
averaged before fitting (using :func:`lsqfit.wavg`). Alternatively (and
equivalently) one could add a third ``Corr2`` to ``models`` for ``Gba``,
but it is more efficient to combine it with ``Gab`` in this way if they are
equivalent.
        
``make_prior(N)``
____________________
This routine defines the fit parameters that correspond to each fit-parameter
label used in ``make_models()`` above. It also assigns *a priori* values to
each parameter, expressed in terms of gaussian random variables (|GVar|\s),
with a mean and standard deviation. The prior is built using class
:class:`gvar.BufferDict`::
        
    import lsqfit
    import gvar
        
    def make_prior(N):
        prior = gvar.BufferDict()       # prior = {}  works too
        prior['a'] = [gvar.gvar(0.1, 0.5) for i in range(N)]
        prior['b'] = [gvar.gvar(1., 5.) for i in range(N)]
        prior['dE'] = [gvar.gvar(0.25, 0.25) for i in range(N)]
        return prior
        
(:class:`gvar.BufferDict` can be replaced by an ordinary Python dictionary;
it is used here because it remembers the order in which the keys are added.)
``make_prior(N)`` associates arrays of ``N`` gaussian random variables
(|GVar|\s) with each fit-parameter label, enough for ``N`` terms in the fit
function. These are the *a priori* values for the fit parameters, and they
can be retrieved using the label: setting ``prior=make_prior(N)``, for
example, implies that ``prior['a'][i]``, ``prior['b'][i]`` and
``prior['dE'][i]`` are the *a priori* values for ``a[i]``, ``b[i]`` and
``dE[i]`` in the fit functions (see above). The *a priori* value for each
``a[i]`` here is set to ``0.1+-0.5``, while that for each ``b[i]`` is
``1+-5``::
        
    >>> print prior['a']
    [0.1 +- 0.5 0.1 +- 0.5 0.1 +- 0.5 0.1 +- 0.5]
    >>> print prior['b']
    [1 +- 5 1 +- 5 1 +- 5 1 +- 5]
        
Similarly the *a priori* value for each energy difference is ``0.25+-0.25``.
(See the :mod:`lsqfit` documentation for further information on priors.)
        

``print_results(fit,prior,data)``
_________________________________
The actual fit is done by ``fit=fitter.lsqfit(...)``, which also prints out
a summary of the fit results (this output can be suppressed if desired).
Further results are reported by ``print_results(fit, prior, data)``: for
example, ::
        
    def print_results(fit, prior, data): 
        a = fit.p['a']                              # array of a[i]s
        b = fit.p['b']                              # array of b[i]s
        dE = fit.p['dE']                            # array of dE[i]s
        E = [sum(dE[:i+1]) for i in range(len(dE))] # array of E[i]s
        print 'Best fit values:
        print '     a[0] =',a[0]                
        print '     b[0] =',b[0]              
        print '     E[0] =',E[0]
        print 'b[0]/a[0] =',b[0]/a[0]
        outputs = {'E0':E[0], 'a0':a[0], 'b0':b[0], 'b0/a0':b[0]/a[0]}
        inputs = {'a'=prior['a'], 'b'=prior['b'], 'dE'=prior['dE'],
                  'data'=[data[k] for k in data])
        print fit.fmt_errorbudget(outputs, inputs)
        
The best-fit values from the fit are contained in ``fit.p`` and are accessed
using the labels defined in the prior and the |Corr2| models. Variables like
``a[0]`` and ``E[0]`` are |GVar| objects that contain means and standard
deviations, as well as information about any correlations that might exist
between different variables (which is relevant for computing functions of the
parameters, like ``b[0]/a[0]`` in this example). 
        
The last line of ``print_results(fit,prior,data)`` prints an error budget for
each of the best-fit results for ``a[0]``, ``b[0]``, ``E[0]`` and
``b[0]/a[0]``, which are identified in the print output by the labels
``'a0'``, ``'b0'``, ``'E0'`` and ``'b0/a0'``, respectively. The error for any
fit result comes from uncertainties in the inputs --- in particular, from the
fit data and the priors. The error budget breaks the total error for a
result down into the components coming from each source. Here the sources are
the *a priori* errors in the priors for the ``'a'`` amplitudes, the ``'b'``
amplitudes, and the ``'dE'`` energy differences, as well as the errors in
the fit data ``data``. These sources are labeled in the print output by
``'a'``, ``'b'``, ``'dE'``, and ``'data'``, respectively. (See the
:mod:`gvar`/:mod:`lsqfit` tutorial for further details on partial standard
deviations and :func:`gvar.fmt_errorbudget`.)
                
Plots of the fit data divided by the fit function, for each correlator, are
displayed by calling ``fitter.display_plots()`` provided the :mod:`matplotlib`
module is present.

.. _faster-fits:   
     
Faster Fits
-----------
Good fits often require fit functions with several exponentials and many
parameters. Such fits can be costly. One strategy that can speed things up is
to use fits with fewer terms to generate estimates for the most important
parameters. These estimates are then used as starting values for the full
fit. The smaller fit is usually faster, because it has fewer parameters, but
the fit is not adequate (because there are too few parameters). Fitting the
full fit function is usually faster given reasonable starting estimates, from
the smaller fit, for the most important parameters. Continuing with the
example from the previous section, the code ::
        
    data = make_data('mcfile')  
    fitter = CorrFitter(models=make_models())
    p0 = None
    for N in [1,2,3,4,5,6,7,8]:
        prior = make_prior(N)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0)
        print_results(fit, prior, data)
        p0 = fit.pmean
            
does fits using fit functions with ``N=1...8`` terms. Parameter mean-values
``fit.pmean`` from the fit with ``N`` exponentials are used as starting values
``p0`` for the fit with ``N+1`` exponentials, hopefully reducing the time
required to find the best fit for ``N+1``.

Faster Fits --- Postive Parameters
----------------------------------
Priors used in |CorrFitter| assign an *a priori* gaussian/normal
distribution to each parameter. It is possible instead to assign a log-normal
distribution, which forces the corresponding parameter to be positive. This is
done by choosing a label in the prior that begins with "log": for example,
``'logdE'`` instead of ``'dE'``. The fitter implements the log-normal
distribution by using the parameter's logarithm, instead of the parameter
itself, as the new fit parameter; the logarithm has a gaussian/normal
distribution. The original parameter is recovered inside the fit function (and
elsewhere) by taking the  exponential of the new fit parameter.

Using log-normal distributions where possible can significantly improve the
stability of a fit. This is because otherwise the fit function typically has
many symmetries that lead to large numbers of equivalent but different best
fits. For example, the fit functions ``Gaa(t,N)`` and ``Gab(t,N)`` above are
unchanged by exchanging ``a[i]``, ``b[i]`` and ``E[i]`` with ``a[j]``,
``b[j]`` and ``E[j]`` for any ``i`` and ``j``. We can remove this degeneracy
by using a log-normal distribution for the ``dE[i]``\s since this guarantees
that all ``dE[i]``\s are positive, and therefore that ``E[0],E[1],E[2]...``
are ordered (in decreasing order of importance to the fit at large ``t``).
        
Another symmetry of ``Gaa`` and ``Gab``, which leaves both fit functions
unchanged, is replacing ``a[i],b[i]`` by ``-a[i],-b[i]``. Yet another is to
add a new term to the fit functions with ``a[k],b[k],dE[k]`` where ``a[k]=0``
and the other two have arbitrary values. Both of these symmetries can be
removed by using a log-normal distribution for the ``a[i]`` priors, thereby
forcing all ``a[i]>0``.
        
The log-normal distributions for the ``a[i]`` and ``dE[i]`` are introduced
into the code example above by changing the corresponding labels in
``make_prior(N)``,  and taking logarithms of the corresponding prior values::
        
    from gvar import log                        # numpy.log() works too
        
    def make_models():                          # same as before
        models = [ Corr2(datatag='Gaa', tdata=range(64), tfit=range(64),
                         a='a', b='a', dE='dE'),
                        
                   Corr2(datatag='Gab', tdata=range(64), tfit=range(64),
                         a='a', b='b', dE='dE')
                 ]
        return models
        
    def make_prior(N):
        prior = gvar.BufferDict()               # prior = {}  works too
        prior['loga'] = [log(gvar.gvar(0.1, 0.5)) for i in range(N)]
        prior['b'] = [gvar.gvar(1., 5.) for i in range(N)]
        prior['logdE'] = [log(gvar.gvar(0.25, 0.25)) for i in range(N)]
        return prior
        
This replaces the original fit parameters, ``a[i]`` and ``dE[i]``, by new fit
parameters, ``log(a[i])`` and ``log(dE[i])``. The *a priori* distributions for
the logarithms are gaussian/normal, with priors of ``log(0.1+-0.5)`` and
``log(0.25+-0.25)`` for the ``log(a)``\s and ``log(dE)``\s respectively. 

Note that the labels are unchanged here in ``make_models()``. It is
unnecessary to change labels in the models; |CorrFitter| will automatically
connect the  modified terms in the prior with the appropriate terms in the
models. This allows one to switch back and forth between log-normal and normal
distributions without changing the models --- only the names in the prior
need be changed. |CorrFitter| also supports "sqrt-normal" distributions,
which are indicated by ``"sqrt"`` at the start of a parameter-name in the
prior; the actual parameter in the fit function is the square of this fit-
parameter, and so is again positive.

Note also that only three lines in ``print_results(fit,prior,data)``, above,
would change if we had used log-normal priors for ``a`` and ``dE``::

    from gvar import exp                        # numpy.exp() works too
    ...
    a = exp(fit.p['loga'])                      # array of a[i]s
    ...
    dE = exp(fit.p['logdE'])                    # array of dE[i]s
    ...
    inputs = {'loga':prior['loga'], 'b':prior['b'], 'logdE':fit.prior['logdE'],
              'data':[data[k] for k in data]}
    ...

Faster Fits --- Chained Fits
----------------------------
Large complicated fits, where lots of models and data are fit simultaneously,
can  take a very long time. This is especially true if there are strong
correlations in the data. Such correlations can also cause  problems from
numerical roundoff errors when the inverse of the data's covariance matrix is
computed for the ``chi**2`` function, requiring large *svd* cuts which can
degrade precision (see below). An alternative approach, that often leads to
more robust and faster fits, is to used *chained* fits.  In a chained fit,
each model is fit by itself in sequence, but with the best-fit parameters from
each fit serving as priors for fit parameters in  the next fit. All parameters
from one fit become fit parameters in the next, including those parameters
that are  not explicitly needed by the next fit (since they may be correlated
with the input data for the next fit or with its priors). Statistical 
correlations between data/priors from different models are preserved 
throughout (this is essential to the method).

The results from a chained fit are identical to a standard simultaneous fit in
the limit of large statistics (that is, in the gaussian limit), but a  chained
fit never involves fitting more than a single correlator at a time. Single-
correlator fits are generally much faster than simultaneous multi-correlator
fits, and roundoff errors (and therefore *svd* cuts) are much less of a 
problem. Consequently chained fits can be more accurate in practice than
conventional simultaneous fits, especially for high-statistics data.

Converting to chained fits is trivial: simply replace ``fitter.lsqfit(...)``
by ``fitter.chained_lsqfit(...)``. The output from this function represents
the results for the entire chain of fits, and so can be used in exactly the
same way as the  output from :func:`fitter.lsqfit` (and is usually quite
similar, to within statistical errors). Results from the different links in
the chain --- that is, from the fits for individual models --- can be accessed
after the fit using ``fitter.chained_fits[datatag]`` where ``datatag`` is the
data tag for the model of interest.


Faster Fits --- Marginalization
-------------------------------       
Often we care only about parameters in the leading term of the fit function,
or just a few of the leading terms. The non-leading terms are needed for a
good fit, but we are uninterested in the values of their parameters. In such
cases the non-leading terms can be absorbed into the fit data, leaving behind
only the leading terms to be fit (to the modified fit data) --- non-leading
parameters are, in effect, integrated out of the analysis, or *marginalized*.
The errors in the modified data are adjusted to account for uncertainties in
the marginalized terms, as specified by their priors. The resulting fit
function has many fewer parameters, and so the fit can be much faster.
        
Continuing with the example in :ref:`faster-fits`, imagine that ``Nmax=8``
terms are needed to get a good fit, but we only care about parameter values
for the first couple of terms. The code from that section can be rearranged to
fit only the leading ``N`` terms where ``N<Nmax``, while incorporating
(marginalizing) the remaining, non-leading terms as corrections to the data::
        
    Nmax = 8
    data = make_data('mcfile')  
    prior = make_prior(Nmax)        # build priors for Nmax terms
    models = make_models()
    p0 = None
    for N in [1,2,3]:
        fitter = CorrFitter(models=models, nterm=N)  # fit only N terms
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0)
        print_results(fit, prior, data)
        p0 = fit.pmean
        
Here the ``nterm`` parameter in |CorrFitter| specifies how many terms are used
in the fit functions. The prior specifies ``Nmax`` terms in all, but only
parameters in ``nterm=N`` terms are varied in the fit. The remaining terms
specified by the prior are automatically incorporated into the fit data by
|CorrFitter|.
        
Remarkably this method is often as accurate with ``N=1`` or ``2`` as a full
``Nmax``-term fit with the original fit data; but it is much faster. If this
is not the case, check for singular priors, where the mean is much smaller
than the standard deviation. These can lead to singularities in the covariance
matrix for the corrected fit data. Such priors are easily fixed: for example,
use ``gvar.gvar(0.1,1.)`` rather than ``gvar.gvar(0.0,1.)`` or
``gvar.gvar(0.001,1.)``. In some situations an *svd* cut (see below) can also
help.
      

Variations
----------
Any 2-point correlator can be turned into a periodic function of ``t`` by
specifying the period through parameter ``tp``. Doing so causes the 
replacement (for ``tp>0``) ::
        
    exp(-E[i]*t)   ->   exp(-E[i]*t) + exp(-E[i]*(tp-t))
        
in the fit function. If ``tp`` is negative, the function is replaced by 
an anti-periodic function with period ``abs(tp)`` and (for ``tp<0``)::
        
    exp(-E[i]*t)   ->   exp(-E[i]*t) - exp(-E[i]*(abs(tp)-t))
        
Also (or alternatively) oscillating terms can be added to the fit by
modifying parameter ``s`` and by specifying sources, sinks and energies for
the oscillating pieces. For example, one might want to replace the sum of
exponentials with two sums ::
        
    sum_i a[i]**2 * exp(-E[i]*t) - sum_i ao[i]**2 (-1)**t * exp(-Eo[i]*t)
        
in a (nonperiodic) fit function. Then an appropriate model
would be, for example, ::
        
    Corr2(datatag='Gaa', tdata=range(64), tfit=range(64),
          a=('a','ao'), b=('a','ao'), dE=('logdE','logdEo'), s=(1,-1))
        
where ``ao`` and ``dEo`` refer to additional fit parameters describing
the oscillating component. In general parameters for amplitudes and
energies can be tuples with two components: the first describing normal
states, and the second describing oscillating states. To omit one or the
other, put ``None`` in place of a label. Parameter ``s[0]`` is an overall
factor multiplying the non-oscillating terms, and ``s[1]`` is the
corresponding factor for the oscillating terms.
        
An *svd* cut can be applied to the covariance matrix for the data by
specifying parameters ``svdcut`` and/or ``svdnum``. (See documentation for
:mod:`lsqfit`; it is useful to set ``svdnum`` equal to the number of
measurements used to determine the covariance matrix for ``G(t)`` since
that is the largest number of eigenmodes possible in the covariance
matrix.)
  
Very Fast (But Limited) Fits
-----------------------------
At large ``t``, correlators are dominated by the term with the smallest
``E``, and often it is only the parameters in that leading term that are
needed. In such cases there is a very fast analysis that is often almost
as accurate as a full fit. An example is::
    
    from corrfitter import fastfit
    
    data = make_data('mcfile')    # user-supplied routine - fit data
    N = 10                        # number of terms in fit functions
    prior = make_prior(N)         # user-supplied routine - fit prior
    model = Corr2(a=..., b=..., ...) # create model describing correlator
    fit = fastfit(data=data, prior=prior, model=model)
    print('E[0] =', fit.E)                  # E[0]
    print('a[0]*b[0] =', fit.ampl)          # a[0]*b[0]
    print('chi2/dof =', fit.chi2/fit.dof)   # good fit if of order 1 or less
    print('Q =', fit.Q)             # good fit if Q bigger than about 0.1
    
:class:`fastfit` estimates ``E[0]`` by using the prior, in effect, to
remove (*marginalize*) all terms from the correlator other than the
``E[0]`` term: so the data ``Gdata(t)`` for the correlator is replaced by,
for example, ::

    Gdata(t) - sum_i=1..N-1  a[i]*b[i] * exp(-E[i]*t)

where ``a[i]``, ``b[i]``, and ``E[i]`` for ``i>0`` are replaced by their
values in the prior. The modified prior is then fit by a single term,
``a[0] * b[0] * exp(-E[0]*t)``, which means that a fit is not necessary
(since the functional form is so simple). It is important to check the
``chi**2`` of the fit, to make sure the fit is good. If it is not, try
restricting ``model.tfit`` to larger ``t``\s (:class:`fastfit` averages
estimates from all ``t``\s in ``model.tfit``).
    
The marginalization of terms with larger ``E``\s allows the code to use
information from much smaller ``t``\s than otherwise, increasing precision.
It also quantifies the uncertainty caused by the existence of these terms.
This simple analysis is a special case of the more general marginalization
strategy discussed in :ref:`faster-fits`, above.
    
3-Point Correlators
-------------------
Correlators ``Gavb(t,T) = <b(T) V(t) a(0)>`` can also be included in fits
as functions of ``t``. In the illustration above, for example, we might
consider additional Monte Carlo data describing a form factor with the 
same intermediate states before and after ``V(t)``. Assuming the data is
tagged by ``aVbT15`` and describes ``T=15``, the corresponding entry in the
collection of models might then be::
    
    Corr3(datatag="aVbT15", T=15, tdata=range(16), tfit=range(16),
        Vnn='Vnn',                # parameters for V
        a='a', dEa='dE',          # parameters for a->V
        b='b', dEb='dE',          # parameters for V->b
        )
    
This models the Monte Carlo data for the 3-point function using the
following formula::
    
    sum_i,j a[i] * exp(-Ea[i]*t) * Vnn[i,j] * b[j] * exp(-Eb[j]*t)
    
where the ``Vnn[i,j]``\s are new fit parameters related to ``a->V->b`` form
factors. Obviously multiple values of ``T`` can be studied by including
multiple |Corr3| models, one for each value of ``T``. Either or both of the
initial and final states can have oscillating components (include ``sa``
and/or ``sb``), or can be periodic (include ``tpa`` and/or ``tpb``). If
there are oscillating states then additional ``V``\s must be specified:
``Vno`` connecting a normal state to an oscillating state, ``Von``
connecting oscillating to normal states, and ``Voo`` connecting oscillating
to oscillating states.
    
There are two cases that require special treatment. One is when
simultaneous fits are made to ``a->V->b`` and ``b->V->a``. Then the
``Vnn``, ``Vno``, *etc.* for ``b->V->a`` are the (matrix) transposes of
the the same matrices for ``a->V->b``. In this case the models for the two
would look something like::
    
    models = [ 
        ...
        Corr3(datatag="aVbT15", T=15, tdata=range(16), tfit=range(16),
            Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo',
            a=('a','ao'), dEa=('dE','dEo'), sa=(1,-1), # a->V
            b=('b','bo'), dEb=('dE','dEo'), sb=(1,-1)  # V->b
            ),
        Corr3(datatag="bVaT15", T=15, tdata=range(16), tfit=range(16),
            Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo', transpose_V=True,
            a=('b','bo'), dEa=('dE','dEo'), sa=(1,-1), # b->V
            b=('a','ao'), dEb=('dE','dEo'), sb=(1,-1)  # V->a
            ),
        ...
    ]
    
The same ``V``\s are specified for the second correlator, but setting
``transpose_V=True`` means that the transpose of each matrix is used
in the fit for that correlator.
    
The second special case is for fits to ``a->V->a`` where source and sink
are the same. In that case, ``Vnn`` and ``Voo`` are symmetric matrices, and
``Von`` is the transpose of ``Vno``. The model for such a case would look
like::
    
    Corr3(datatag="aVbT15", T=15, tdata=range(16), tfit=range(16),
        Vnn='Vnn', Vno='Vno', Von='Vno', Voo='Voo', symmetric_V=True,
        a=('a','ao'), dEa=('dE', 'dEo'), sa=(1, -1), # a->V
        b=('a','ao'), dEb=('dE', 'dEo'), sb=(1, -1)  # V->a
        )
    
Here ``Vno`` and ``Von`` are set equal to the same matrix, but specifying
``symmetric_V=True`` implies that the transpose will be used for ``Von``.
Furthermore ``Vnn`` and ``Voo`` are symmetric matrices when
``symmetric_V==True`` and so only the upper part of each matrix is needed.
In this case ``Vnn`` and ``Voo`` are treated as one-dimensional arrays with
``N(N+1)/2`` elements corresponding to the upper parts of each matrix,
where ``N`` is the number of exponentials (that is, the number of
``a[i]``\s).
    
Bootstrap Analyses
------------------
A *bootstrap analysis* gives more robust error estimates for fit parameters
and functions of fit parameters than the conventional fit when errors are
large, or fluctuations are non-gaussian. A typical code looks something like::
    
    import gvar as gv
    import gvar.dataset as ds
    from corrfitter import CorrFitter
    # fit
    dset = ds.Dataset('mcfile')
    data = ds.avg_data(dset)            # create fit data
    fitter = Corrfitter(models=make_models())
    N = 4                               # number of terms in fit function
    prior = make_prior(N)
    fit = fitter.lsqfit(prior=prior, data=data)  # do standard fit
    print 'Fit results:'
    print 'a',exp(fit.p['loga'])        # fit results for 'a' amplitudes
    print 'dE',exp(fit.p['logdE'])      # fit results for 'dE' energies
    ...
    ...
    # bootstrap analysis
    print 'Bootstrap fit results:'
    nbootstrap = 10                     # number of bootstrap iterations
    bs_datalist = (ds.avg_data(d) for d in dset.bootstrap_iter(nbootstrap))
    bs = ds.Dataset()                   # bootstrap output stored in bs
    for bs_fit in fitter.bootstrap_iter(bs_datalist): # bs_fit = lsqfit output
        p = bs_fit.pmean    # best fit values for current bootstrap iteration
        bs.append('a', exp(p['loga']))  # collect bootstrap results for a[i]
        bs.append('dE', exp(p['logdE']))# collect results for dE[i]
        ...                             # include other functions of p 
        ...
    bs = ds.avg_data(bs, bstrap=True)   # medians + error estimate
    print 'a', bs['a']                  # bootstrap result for 'a' amplitudes
    print 'dE', bs['dE']                # bootstrap result for 'dE' energies
    ....
        
This code first prints out the standard fit results for the ``'a'`` amplitudes
and ``'dE'`` energies. It then makes ``10`` bootstrap copies of the original
input data, and fits each using the best-fit parameters from the original fit
as the starting point for the bootstrap fit. The variation in the best-fit
parameters from fit to fit is an indication of the uncertainty in those
parameters. This example uses a :class:`gvar.dataset.Dataset` object ``bs`` to
accumulate the results from each bootstrap fit, which are computed using the
best-fit values of the parameters (ignoring their standard deviations). Other
functions of the fit parameters could be included as well. At the end
``avg_data(bs, bstrap=True)`` finds median values for each quantity in
``bs``, as well as a robust estimate of the uncertainty (to within 30% since
``nbootstrap`` is only ``10``).
    
The list of bootstrap data sets ``bs_datalist`` can be omitted in this example
in situations where the input data has high statistics. Then the bootstrap
copies are generated internally by :func:`fitter.bootstrap_iter()` from the
means and covariance matrix of the input data (assuming gaussian statistics).
    
    
New Models
----------
Classes to describe new models are usually derived from |BaseModel|. These
can be for fitting new types of correlators. They can also be used in other
ways --- for example, to add constraints. Imagine a situation where one
wants to constrain the third energy (``E2``) in a fit to be ``0.60(1)``.
This can be accomplished by adding ``E2_Constraint()`` to the list of
models in |Corrfitter| where::
    
    import gvar
    import corrfitter
    
    class E2_Constraint(corrfitter.BaseModel):
        def __init__(self):
            super(E2_Constraint,self).__init__('E2-constraint') # data tag
        
        def fitfcn(self,x,p):
            dE = gvar.exp(p['logdE'])
            return sum(dE[:3])              # E2 formula in terms of p
        
        def builddata(self,d):
            return gvar.gvar(0.6,0.01)      # E2 value
            
        def buildprior(self, prior, nterm):
            return {}
        
Any number of constraints like this can be added to the list of models.
    
Note that this constraint could instead be built into the priors for
``logdE`` by introducing correlations between different parameters.
    
Implementation
--------------
|CorrFitter| allows models to specify how many exponentials to include in the
fit function (using parameters ``nterm``, ``nterma`` and ``ntermb``). If that
number is less than the number of exponentials specified by the prior, the
extra terms are incorporated into the fit data before fitting. The default
procedure is to multiply the data by ``G(t,p,N)/G(t,p,max(N,Nmax))`` where:
``G(p,t,N)`` is the fit function with ``N`` terms for parameters ``p`` and
time ``t``; ``N`` is the number of exponentials specified in the models;
``Nmax`` is the number of exponentials specified in the prior; and here
parameters ``p`` are set equal to their values in the prior (correlated
|GVar|\s).
    
An alternative implementation for the data correction is to add
``G(t,p,N)-G(t,p,max(N,Nmax))`` to the data. This implementation is selected
when parameter ``ratio`` in |CorrFitter| is set to ``False``. Results are
similar to the other implementation, though perhaps a little less robust.

Background information on the some of the fitting strategies used by
|corrfitter| can be found by doing web searches for "hep-lat/0110175" and
"arXiv:1111.1363". These are two papers by G.P. Lepage and collaborators 
whose published versions are: G.P. Lepage et al, Nucl.Phys.Proc.Suppl. 
106 (2002) 12-20; and K. Hornbostel et al, Phys.Rev. D85 (2012) 031504.


Correlator Model Objects
-------------------------
Correlator objects describe theoretical models that are fit to
correlator data by varying the models' parameters. 

A model object's parameters are specified through priors for the fit. A
model assigns labels to each of its parameters (or arrays of related
parameters), and these labels are used to identify the corresponding
parameters in the prior. Parameters can be shared by more than one model
object.

A model object also specifies the data that it is to model. The data is 
identified by the data tag that labels it in the input file or |Dataset|.

.. autoclass:: corrfitter.Corr2(datatag, tdata, tfit, a, b, dE, s=1.0, tp=None, othertags=None)
   :members:
   
.. autoclass:: corrfitter.Corr3(datatag, T, tdata, tfit, Vnn, a, b, dEa, dEb, sa=1.0, sb=1.0, Vno=None, Von=None, Voo=None, transpose_V=False, symmetric_V=False, tpa=None, tpb=None, othertags=None)
   :members:
   
.. autoclass:: corrfitter.BaseModel
   :members:


|CorrFitter| Objects
---------------------
|CorrFitter| objects are wrappers for :func:`lsqfit.nonlinear_fit()` which
is used to fit a collection of models to a collection of Monte Carlo data.

.. autoclass:: corrfitter.CorrFitter
   :members:

Fast Fit Objects
-----------------

.. autoclass:: corrfitter.fastfit
