corrfitter
----------
This module contains tools for least-squares fits, as functions
of time ``t``, of simulation (or other statistical) data for 2-point and
3-point correlators of the form::

    Gab(t)    =  <b(t) a(0)>
    Gavb(t,T) =  <b(T) V(t) a(0)>

Each correlator is modeled using ``Corr2`` for 2-point correlators, or
``Corr3`` for 3-point correlators in terms of amplitudes for each source
``a``, sink ``b``, and vertex ``V``, and the energies associated with each
intermediate state. The amplitudes and energies are adjusted in the
least-squares fit to reproduce the data; they are specified in a shared prior
(typically a dictionary).

An object of type ``CorrFitter`` describes a collection of correlators and is
used to fit multiple models to data simultaneously. Any number of
correlators may be described and fit by a single ``CorrFitter`` object.
``CorrFitter`` objects can also be used to to extract the appropriate fit
data from ``Dataset`` objects.

This module has been used extensively for analyzing results from lattice QCD
simulations. Documentation for this module can be found by looking at
``doc/html/index.html`` or <https://corrfitter.readthedocs.io>. The ``.py``
files in directory ``examples/`` also illustrate how to use corrfitter; see
``examples/README`` for more information.


Information on how to install the library is in the file INSTALLATION.
To test the module try 'make tests'.

This module requires the ``lsqfit`` and ``gvar`` Python packages.

Background information on the some of the fitting strategies used by
``corrfitter`` can be found by doing a web searches for "hep-lat/0110175",
"arXiv:1111.1363", and :arXiv:1406.2279" (appendix). These are papers by
G.P. Lepage and collaborators whose published versions are:
G.P. Lepage et al, Nucl.Phys.Proc.Suppl. 106 (2002) 12-20;
K. Hornbostel et al, Phys.Rev. D85 (2012) 031504; and
C.M. Bouchard et al, Phys.Rev. D90 (2014) 054506.

``corrfitter`` version numbers have the form ``major.minor.patch`` where
incompatible changes are signaled by incrementing the ``major`` version
number, the ``minor`` number signals new features, and  the ``patch`` number
signals bug fixes.

| Created by G. Peter Lepage (Cornell University) on 2008-02-12.
| Copyright (c) 2008-18 G. Peter Lepage.
