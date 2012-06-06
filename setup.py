from distutils.core import setup
import corrfitter

setup(name='corrfitter',
    version=corrfitter.__version__,
    description='Utilities for fitting correlators.',
    author='G. Peter Lepage, Cornell University',
    author_email='g.p.lepage@cornell.edu',
    license='GPLv3',
    py_modules = ['corrfitter']
)