:mod:`corrfitter` - Least-Squares Fit to Correlators
=====================================================

.. module:: corrfitter
   :synopsis: Least-Squares Fit to Correlators.

.. moduleauthor:: G.P. Lepage <g.p.lepage@cornell.edu>

.. |CorrFitter| replace:: :class:`corrfitter.CorrFitter`
.. |Corr2| replace:: :class:`corrfitter.Corr2`
.. |Corr3| replace:: :class:`corrfitter.Corr3`
.. |EigenBasis| replace:: :class:`corrfitter.EigenBasis`
.. |BaseModel| replace:: :class:`corrfitter.BaseModel`
.. |Dataset| replace:: :class:`gvar.dataset.Dataset`
.. |GVar| replace:: :class:`gvar.GVar`
.. |chi2| replace:: :math:`\chi^2`

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

.. autoclass:: corrfitter.Corr2(datatag, a, b, dE, s=1.0, tp=None, tmin=None, tmax=None, tdata=None, tfit=None, reverse=False, reverseddata=[], otherdata=[])
   :members:

.. autoclass:: corrfitter.Corr3(datatag, T, Vnn, a, b, dEa, dEb, sa=1.0, sb=1.0, Vno=None, Von=None, Voo=None, tdata=None, tfit=None, tmin=None. reverse=False, symmetric_V=False, reverseddata=[], otherdata=[])
   :members:


|CorrFitter| Objects
---------------------
|CorrFitter| objects are wrappers for :func:`lsqfit.nonlinear_fit()` which
is used to fit a collection of models to a collection of Monte Carlo data.

.. autoclass:: corrfitter.CorrFitter
   :members:

:class:`corrfitter.EigenBasis` Objects
--------------------------------------
:class:`corrfitter.EigenBasis` objects are useful for analyzing two-point and
three-point correlators with multiplle sources and sinks.
The current interface for :class:`EigenBasis` is experimental.
It may change in the near future, as experience
accumulates from its use.


.. autoclass:: corrfitter.EigenBasis

    In addition to ``keyfmt``, ``srcs``, ``t`` and ``tdata`` above,
    the main attributes are:

    .. attribute:: E

        Array of approximate energies obtained from the eigenanalysis.

    .. attribute:: eig_srcs

        List of labels for the sources in the eigen-basis: ``'0'``, ``'1'`` ...

    .. attribute:: svdcorrection

        The sum of the SVD corrections added to the data by the last call
        to :func:`EigenBasis.svd`.

    .. attribute:: svdn

        The number of degrees of freedom modified by the SVD correction
        in the last call to :func:`EigenBasis.svd`.

    .. attribute:: v

        ``v[a]`` is the eigenvector corresponding to source ``a``
        in the new basis, where ``a=0,1...``.


    .. attribute:: v_inv

        ``v_inv[i]`` is the inverse-eigenvector for transforming from the
        new basis back to the original basis.

    The main methods are:

    .. automethod:: apply

    .. automethod:: make_prior

    .. automethod:: svd

    .. automethod:: tabulate

    .. automethod:: unapply



Fast Fit Objects
-----------------

.. autoclass:: corrfitter.fastfit
