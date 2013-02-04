""" Introduction 
------------
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
        
The priors assign an *a priori* gaussian or normal distribution to each
parameter. It is possible instead to assign a log-normal distribution, which
forces the parameter to be positive. This is done by choosing a label that
begins with "log": for example, ``'logdE'`` instead of ``'dE'``. The fitter
implements the log-normal distribution by using the parameter's logarithm,
instead of the parameter itself, as a new fit parameter; the logarithm has a
gaussian/normal distribution. The original parameter is recovered by taking
the exponential of the new fit parameter.
        
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
into the code example by changing the corresponding labels in
``make_models()`` and ``make_prior(N)``, and taking logarithms of the
corresponding prior values::
        
    from gvar import log                        # numpy.log() works too
        
    def make_models():
        models = [ Corr2(datatag='Gaa', tdata=range(64), tfit=range(64),
                         a='loga', b='loga', dE='logdE'),
                        
                   Corr2(datatag='Gab', tdata=range(64), tfit=range(64),
                         a='loga', b='b', dE='logdE')
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
        
Note that only three lines in ``print_results(fit,prior,data)`` would change
if we had used log-normal priors for ``a`` and ``dE``, as discussed in the
previous section::
        
    from gvar import exp                        # numpy.exp() works too
    ...
    a = exp(fit.p['loga'])                      # array of a[i]s
    ...
    dE = exp(fit.p['logdE'])                    # array of dE[i]s
    ...
    inputs = {'loga':prior['loga'], 'b':prior['b'], 'logdE':fit.prior['logdE'],
              'data':[data[k] for k in data]}
    ...
        
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
        
Often we care only about parameters in the leading term of the fit function,
or just a few of the leading terms. The non-leading terms are needed for a
good fit, but we are uninterested in the values of their parameters. In such
cases the non-leading terms can be absorbed into the fit data, leaving behind
only the leading terms to be fit (to the modified fit data) --- non-leading
parameters are, in effect, integrated out of the analysis, or *marginalized*.
The errors in the modified data are adjusted to account for uncertainties in
the marginalized terms, as specified by their priors. The resulting fit
function has many fewer parameters, and so the fit can be much faster.
        
Continuing with the example above, imagine that ``Nmax=8`` terms are needed
to get a good fit, but we only care about parameter values for the first
couple of terms. The code above can be rearranged to fit only the leading
``N`` terms where ``N<Nmax``, while incorporating (marginalizing) the
remaining, non-leading terms as corrections to the data::
        
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
        
Remarkably this method is usually as accurate with ``N=1`` or ``2`` as a full
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
        Vnn='Vnn',                  # parameters for V
        a='a',dEa='logdE',          # parameters for a->V
        b='b',dEb='logdE',          # parameters for V->b
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
            a=('a','ao'), dEa=('logdE','logdEo'), sa=(1,-1), # a->V
            b=('b','bo'), dEb=('logdE','logdEo'), sb=(1,-1)  # V->b
            ),
        Corr3(datatag="bVaT15", T=15, tdata=range(16), tfit=range(16),
            Vnn='Vnn', Vno='Vno', Von='Von', Voo='Voo', transpose_V=True,
            a=('b','bo'), dEa=('logdE','logdEo'), sa=(1,-1), # b->V
            b=('a','ao'), dEb=('logdE','logdEo'), sb=(1,-1)  # V->a
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
        a=('a','ao'), dEa=('logdE', 'logdEo'), sa=(1, -1), # a->V
        b=('a','ao'), dEb=('logdE', 'logdEo'), sb=(1, -1)  # V->a
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
"""

# Created by G. Peter Lepage, Cornell University, on 2010-11-26.
# Copyright (c) 2010-2012 G. Peter Lepage.
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
__version__ = '3.2.5'

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
    ##
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
    ##
    def builddata(self, data):
        """ Construct fit data. 
            
        Format of output must be same as format for fitfcn output.
        
        :param data: Dataset containing correlator data 
            (see ``gvar.dataset``).
        :type data: dictionary
        """
        raise NotImplementedError("builddata not defined")
    ##
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
    ##
    def _param(self, p, default=None):
        """ Parse fit-parameter label --- utility function. """
        if isinstance(p, tuple):
            return p
        else:
            return (p, default)
    ##
    def _dE(self, dE, logdE):
        """ Parse fit-parameter label for E differences --- utility fcn. """
        if dE is None:
            assert logdE is not None, "must specify dE"
            logdE = self._param(logdE, None)
            for x in logdE:
                assert (x is None or x[:3] == 'log'), \
                    "logdE label must begin with 'log'"
            return logdE
        else:
            return self._param(dE, None)
    ##              
##

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
        prior. Fit parameters with labels that begin with "log" are
        replaced by their exponentials in the fit function; and in the
        prior, therefore, these parameters would be logarithms of the
        corresponding amplitudes (and the amplitudes must then be positive.)
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
        results as well as in the prior. Fit parameters with labels that
        begin with "log" are replaced by their exponentials in the fit
        function; and in the prior, therefore, these parameters would be
        logarithms of the corresponding amplitudes (and the amplitudes must
        then be positive.)
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
        ## verify and compress tfit ##
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
        ##
        self.tfit = numpy.array(ntfit)
        self.othertags = othertags
        self._abscissa = self.tfit
    ##
    def buildprior(self, prior, nterm):
        newprior = _gvar.BufferDict()
        for ai, bi, dEi, ntermi in zip(self.a, self.b, self.dE, nterm):
            for x in [ai, bi, dEi]:
                if x is None:
                    continue
                newprior[x] = prior[x][None:ntermi]
        return newprior        
    ##
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
    ##
    def fitfcn(self, p, nterm=None):
        """ Return fit function for parameters ``p``. (Ingores ``x``.) """
        ans = 0.0
        t = self.tfit
        if self.tp is None:
            tp_t = None
        elif self.tp >= 0:
            tp_t = self.tp - t
            pfac = 1
        else:
            tp_t = -self.tp - t
            pfac = -1
        # tp_t = None if self.tp is None else self.tp-t
        if nterm is None:
            nterm = (None, None)
        ofac = (None if self.s[0] == 0.0 else self.s[0],
                (None if self.s[1] == 0.0 else self.s[1]*(-1)**t))
        for ai, bi, dEi, ofaci, ntermi in zip(self.a, self.b, 
                                              self.dE, ofac, nterm):
            if ai is None or bi is None or dEi is None or ofaci is None:
                continue
            if ntermi is not None:
                if ntermi == 0:
                    continue
                ai = (p[ai][:ntermi] if ai[:3] != 'log' 
                        else _gvar.exp(p[ai][:ntermi]))
                bi = (p[bi][:ntermi] if bi[:3] != 'log' 
                        else _gvar.exp(p[bi][:ntermi]))
                dEi = (p[dEi][:ntermi] if dEi[:3] != 'log' 
                        else _gvar.exp(p[dEi][:ntermi]))
            else:   
                ai = p[ai] if ai[:3] != 'log' else _gvar.exp(p[ai])
                bi = p[bi] if bi[:3] != 'log' else _gvar.exp(p[bi])
                dEi = p[dEi] if dEi[:3] != 'log' else _gvar.exp(p[dEi])
            sumdE = 0.0
            for a, b, dE in zip(ai, bi, dEi):
                sumdE += dE
                ans += ofaci*a*b*((_gvar.exp(-sumdE*t)) if tp_t is None else 
                              (_gvar.exp(-sumdE*t)+pfac*_gvar.exp(-sumdE*tp_t)) )
        return ans
    ##
##

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
        results as well as in the prior. Fit parameters with labels that
        begin with "log" are replaced by their exponentials in the fit
        function; and in the prior, therefore, these parameters would be
        logarithms of the corresponding amplitudes (and the amplitudes must
        then be positive.)
    :type a: string, or two-tuple of strings or ``None``
    :param b: Same as ``self.a`` except for sink amplitudes ``(bn, bo)`` 
        for ``V->b`` rather than for ``(an, ao)``.
    :type b: string, or two-tuple of strings or ``None``
    :param dEa: Fit-parameter label for ``a->V`` intermediate-state energy 
        differences ``dEan``, or two-tuple of labels for the differences
        ``(dEan,dEao)``. Each label represents an array of energy differences.
        Replacing either label by ``None`` causes the corresponding term in
        the correlator function to be dropped. Fit parameters with labels that
        begin with "log" are replaced by their exponentials in the fit
        function; here, therefore, these parameters would be logarithms of the
        corresponding energy differences, which differences must then be
        positive.
    :type dEa: string, or two-tuple of strings or ``None``
    :param dEb: Fit-parameter label for ``V->b`` intermediate-state energy 
        differences ``dEbn``, or two-tuple of labels for the differences
        ``(dEbn,dEbo)``. Each label represents an array of energy differences.
        Replacing either label by ``None`` causes the corresponding term in
        the correlator function to be dropped. Fit parameters with labels that
        begin with "log" are replaced by their exponentials in the fit
        function; here, therefore, these parameters would be logarithms of the
        corresponding energy differences, which differences must then be
        positive.
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
        ## verify tfit ##
        ntfit = []
        for t in tfit:
            if t >= 0 and t <= T:
                ntfit.append(t)
        ##
        self.tfit = numpy.array(ntfit)
        self._abscissa = self.tfit
        self.othertags = othertags
    ##
    def buildprior(self, prior, nterm):
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
        ##
        ans = _gvar.BufferDict()
        for x in [self.a, self.dEa, self.b, self.dEb]:
            for xi, ntermi in zip(x, nterm):
                if xi is not None:
                    ans[xi] = prior[xi][None:ntermi]
        for i in range(2):
            for j in range(2):
                vij = self.V[i][j]
                if vij is None:
                    continue
                if i == j and self.symmetric_V:
                    ans[vij] = resize_sym(prior[vij], nterm)
                else:
                    ans[vij] = prior[vij][slice(None, nterm[i]), 
                                          slice(None, nterm[j])]
        return ans         
    ##
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
    ##
    def fitfcn(self, p, nterm=None):
        """ Return fit function for parameters ``p``. """
        ta = self.tfit
        tb = self.T-self.tfit
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
        # tp_ta = None if self.tpa is None else self.tpa-ta
        # tp_tb = None if self.tpb is None else self.tpb-tb
        if nterm is None:
            nterm = (None, None)
        ## initial and final propagators ## 
        aprop = []  # aprop[i][j] i= n or o; j=excitation level
        ofac = (self.sa[0], (0.0 if self.sa[1] == 0.0 else self.sa[1]*(-1)**ta))
        for ai, dEai, ofaci, ntermai in zip(self.a, self.dEa, ofac, nterm):
            if ai is None:
                aprop.append(None)
                continue
            ans = []
            sumdE = 0.0
            if ntermai is None:
                ai =  p[ai] if ai[:3] != 'log' else _gvar.exp(p[ai])
                dEai = p[dEai] if dEai[:3] != 'log' else _gvar.exp(p[dEai])
            else:
                if ntermai <= 0:
                    aprop.append(None)
                    continue
                ai =  (p[ai][:ntermai] if ai[:3] != 'log' 
                            else _gvar.exp(p[ai][:ntermai]))
                dEai = (p[dEai][:ntermai] if dEai[:3] != 'log' 
                            else _gvar.exp(p[dEai][:ntermai]))
            for a, dE in zip(ai, dEai):
                sumdE += dE
                ans.append(ofaci*a*((_gvar.exp(-sumdE*ta)) 
                        if tp_ta is None else 
                        (_gvar.exp(-sumdE*ta)+pafac*_gvar.exp(-sumdE*tp_ta))))
            aprop.append(ans)
        bprop = []
        ofac = (self.sb[0], (0.0 if self.sb[1] == 0.0 else self.sb[1]*(-1)**tb))
        for bi, dEbi, ofaci, ntermbi in zip(self.b, self.dEb, ofac, nterm):
            if bi is None:
                bprop.append(None)
                continue
            ans = []
            sumdE = 0.0
            if ntermbi is None:
                bi = p[bi] if bi[:3] != 'log' else _gvar.exp(p[bi])
                dEbi = p[dEbi] if dEbi[:3] != 'log' else _gvar.exp(p[dEbi])
            else:
                if ntermbi <= 0:
                    bprop.append(None)
                    continue
                bi = (p[bi][:ntermbi] if bi[:3] != 'log' 
                            else _gvar.exp(p[bi][:ntermbi]))
                dEbi = (p[dEbi][:ntermbi] if dEbi[:3] != 'log' 
                            else _gvar.exp(p[dEbi][:ntermbi]))
            for b, dE in zip(bi, dEbi):
                sumdE += dE
                ans.append(ofaci*b*((_gvar.exp(-sumdE*tb)) 
                        if tp_tb is None else 
                        (_gvar.exp(-sumdE*tb)+pbfac*_gvar.exp(-sumdE*tp_tb))))
            bprop.append(ans)
        ##
        ## combine propagators with vertices ##
        ans = 0.0
        for i, (apropi, Vi) in enumerate(zip(aprop, self.V)):
            if apropi is None:
                continue
            for j, (bpropj, Vij) in enumerate(zip(bprop, Vi)):
                if bpropj is None or Vij is None:
                    continue
                V = _gvar.exp(p[Vij]) if Vij[:3] == 'log' else p[Vij]
                if i == j and self.symmetric_V:
                    ## unpack symmetric matrix V ##
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
                    ##
                if self.transpose_V or (i>j and self.symmetric_V):
                    V = V.T
                for ak, Vk in zip(apropi, V):
                    acc = 0.0
                    for bl, Vkl in zip(bpropj, Vk):
                        acc += Vkl*bl
                    ans += ak*acc
        ##
        return ans                
    ##
##                    
            
class CorrFitter(object):
    """ Nonlinear least-squares fitter for a collection of correlators. 
        
    :param models: Correlator models used to fit statistical input data.
    :type models: list of correlator models
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
    def __init__(self, models, svdcut=None, svdnum=None, tol=1e-10, #):
                maxit=500, nterm=None, ratio=True): # mc=None, ratio=True):
        super(CorrFitter, self).__init__()
        self.models = models
        self.svdcut = svdcut
        self.svdnum = svdnum
        self.tol = tol
        self.maxit = maxit
        self.fit = None
        self.dset = None
        # self.mc = mc
        self.ratio = ratio
        self.nterm = nterm if isinstance(nterm, tuple) else (nterm, None)
    ##
    def fitfcn(self, p, nterm=None):
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
        ans = {}
        if nterm is None:
            nterm = self.nterm
        for m in self.models:
            ans[m.datatag] = m.fitfcn(p, nterm=nterm)
        return ans
    ## 
    def builddata(self, data, prior, nterm=None):
        """ Build fit data, corrected for marginalized terms. """
        fitdata = _gvar.BufferDict()
        for m in self.models:
            fitdata[m.datatag] = m.builddata(data)
        ## remove marginal fit parameters ##
        if nterm is None:
            nterm = self.nterm
        ## use priors to remove marginal parameters ##
        if nterm == (None, None):
            return fitdata
        else:
            for m in self.models:
                ftrunc = m.fitfcn(prior, nterm=nterm)
                fall = m.fitfcn(prior, nterm=None)
                if not self.ratio:
                    diff = ftrunc-fall
                    fitdata[m.datatag] += diff 
                else:
                    ii = (_gvar.mean(fall) != 0.0)
                    ratio = ftrunc[ii]/fall[ii]
                    fitdata[m.datatag][ii] *= ratio 
        ##
        ##
        return fitdata
    ##
    def buildprior(self, prior, nterm=None):
        """ Build correctly sized prior for fit. """
        tmp = _gvar.BufferDict()
        if nterm is None:
            nterm = self.nterm
        for m in self.models:
            tmp.update(m.buildprior(prior=prior, nterm=nterm))
        ## restore order of keys --- same as prior ##
        ans = _gvar.BufferDict()
        for k in prior:
            if k in tmp:
                ans[k] = tmp[k]
        ##
        return ans
    ##
    def lsqfit(self, data, prior, p0=None, print_fit=True, nterm=None,  #):
            svdcut=None, svdnum=None, tol=None, maxit=None, **args):
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
            
        The following parameters overwrite the values specified in the
        |CorrFitter| constructor when set to anything other than ``None``:
        ``nterm``, ``svdcut``, ``svdnum``, ``tol``, and ``maxit``. Any
        further keyword arguments are passed on to
        :func:`lsqfit.nonlinear_fit`, which does the fit.
        """
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
        data = self.builddata(data=data, prior=prior, nterm=nterm)
        prior = self.buildprior(prior, nterm=nterm)
        self.fit = lsqfit.nonlinear_fit( #
            data=data, p0=p0, fcn=self.fitfcn,
            prior=prior, svdcut=svdcut, svdnum=svdnum,
            reltol=tol, abstol=tol, maxit=maxit, **args)
        if print_fit:
            print(self.fit.format())
        return self.fit
    ##
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
    ##
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
    ##
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
        ##
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
        ##
        fig.canvas.mpl_connect('key_press_event', onpress)
        plotdata(idx)
        plt.show()
    ##
##

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
        t = model.tfit
        assert numpy.all(t[1:] == t[:-1]+1), "model.tfit must contain consecutive ts"
        self.osc = osc
        if osc:
            ## capture leading oscillating part ##
            nterm = (0, 1)
            Gfac = (-1)**t
            a = model.a[1]
            b = model.b[1]
            dE = model.dE[1]
            ##
        else:
            ## capture leading non-oscillating part ##
            nterm = (1, 0)
            Gfac = 1.
            a = model.a[0]
            b = model.b[0]
            dE = model.dE[0]
            ##
        ## compute priors for answers ##
        E_prior = numpy.exp(prior[dE][0]) if dE[:3] == "log" else prior[dE][0]
        ampl_prior = numpy.exp(prior[a][0]) if a[:3] == "log" else prior[a][0]
        ampl_prior *= numpy.exp(prior[b][0]) if b[:3] == "log" else prior[b][0]
        ##
        ## extract relevant data ##
        fitter = CorrFitter(models=[model], svdcut=svdcut, svdnum=svdnum,
                            ratio=ratio, nterm=nterm)
        G = fitter.builddata(data=data, prior=prior)[model.datatag] * Gfac
        ##
        ## compute E ##
        self.Elist = _gvar.arccosh(0.5*(G[2:]+G[:-2])/G[1:-1])
        Elist = self.Elist.tolist() + [E_prior]
        self.E = lsqfit.wavg(Elist, svdcut=svdcut, svdnum=svdnum)
        self.chi2 = lsqfit.wavg.chi2
        self.dof = lsqfit.wavg.dof
        self.Q = lsqfit.wavg.Q
        ##
        ## compute amplitude ##
        p = {}
        p[a] = numpy.array([0.]) if a[:3] == "log" else numpy.array([1.])
        p[b] = numpy.array([0.]) if b[:3] == "log" else numpy.array([1.])
        p[dE] = numpy.log([self.E]) if dE[:3] == "log" else numpy.array([self.E])
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
    ##
##       

