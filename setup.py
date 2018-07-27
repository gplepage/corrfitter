from distutils.command.build_py import build_py as _build_py
from distutils.core import setup

CORRFITTER_VERSION = '8.0.1'

class build_py(_build_py):
    def run(self):
        """ Append version number to corrfitter.py """
        with open('src/corrfitter.py', 'a') as cffile:
            cffile.write("\n__version__ = '%s'\n" % CORRFITTER_VERSION)
        _build_py.run(self)

# pypi
with open('README.rst', 'r') as file:
    long_description = file.read()

setup(name='corrfitter',
    version=CORRFITTER_VERSION,
    description='Utilities for fitting correlators in lattice QCD.',
    author='G. Peter Lepage, Cornell University',
    author_email='g.p.lepage@cornell.edu',
    license='GPLv3+',
    packages={''},
    package_dir={'':'src'},
    cmdclass={'build_py': build_py},
    requires=["lsqfit (>=11.0)", 'numpy (>=1.7)', 'gvar (>=9.0)'],
    install_requires=['lsqfit>=10.0', 'gvar>=9.0', 'numpy>=1.7'],
    platforms="Any",
    url="https://github.com/gplepage/corrfitter.git",
    long_description=long_description,
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
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Physics'
        ]

)