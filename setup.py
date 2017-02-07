from distutils.core import setup
import corrfitter

setup(name='corrfitter',
    version=corrfitter.__version__,
    description='Utilities for fitting correlators in lattice QCD.',
    author='G. Peter Lepage, Cornell University',
    author_email='g.p.lepage@cornell.edu',
    license='GPLv3+',
    py_modules=['corrfitter'],
    requires=["lsqfit (>=9.1.1)", 'numpy (>=1.7)', 'gvar (>=8.2)'],
    install_requires=['lsqfit>=9.1.1', 'gvar>=8.2', 'numpy>=1.7'],
    platforms="Any",
    url="https://github.com/gplepage/corrfitter.git",
    long_description="""\
    This module contains tools that facilitate least-squares fits, as functions
    of time ``t``, of simulation (or other statistical) data for 2-point and
    3-point correlators of the form::

        Gab(t)    =  <b(t) a(0)>
        Gavb(t,T) =  <b(T) V(t) a(0)>

    Each correlator is modeled using |Corr2| for 2-point correlators, or
    |Corr3| for 3-point correlators in terms of amplitudes for each source
    ``a``, sink ``b``, and vertex ``V``, and the energies associated with each
    intermediate state. The amplitudes and energies are adjusted in the
    least-squares fit to reproduce the data; they are specified in a shared prior
    (typically a dictionary).

    An object of type |CorrFitter| describes a collection of correlators and is
    used to fit multiple models to data simultaneously. Any number of
    correlators may be described and fit by a single |CorrFitter| object.
    |CorrFitter| objects can also be used to to extract the appropriate fit
    data from |Dataset| objects.

    This module has been used extensively for analyzing results from lattice
    QCD simulations.
    """
    ,
    classifiers = [                     #
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Topic :: Scientific/Engineering :: Physics'
        ]

)