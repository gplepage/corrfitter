# from distutils.command.build_py import build_py
from distutils.core import setup

CORRFITTER_VERSION = open('src/corrfitter/_version.py', 'r').readlines()[0].split("'")[1]

# pypi
with open('README.rst', 'r') as file:
    long_description = file.read()

setup(name='corrfitter',
    version=CORRFITTER_VERSION,
    description='Utilities for fitting correlators in lattice QCD.',
    author='G. Peter Lepage, Cornell University',
    author_email='g.p.lepage@cornell.edu',
    license='GPLv3+',
    packages={'corrfitter'},
    package_dir={'corrfitter':'src/corrfitter'},
    # cmdclass={'build_py': build_py},
    requires=["lsqfit (>=11.6)", 'numpy (>=1.7)', 'gvar (>=11.6)'],
    install_requires=['lsqfit>=11.6', 'gvar>=11.6', 'numpy>=1.7'],
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