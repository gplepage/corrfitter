from __future__ import print_function   # makes this work for python2 and 3

import gvar as gv
import numpy as np
import collections
from corrfitter import CorrFitter, Corr2, Corr3

DISPLAYPLOTS = False         # display plots at end of fitting
try:
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

def main():
    data = make_data('etas-Ds.h5')
    models = make_models()                                                     # 1
    models = [models[0], models[1], dict(nterm=(1,0)), (models[2], models[3])] # 1
    fitter = CorrFitter(models=models)                                         # 1
    p0 = None
    for N in [1, 2, 3, 4]:
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N)
        fit = fitter.chained_lsqfit(data=data, prior=prior, p0=p0)             # 2
        print(fit.formatall(pstyle=None if N < 4 else 'v'))                    # 3
        p0 = fit.pmean
    print_results(fit, prior, data)
    if DISPLAYPLOTS:
        fit.show_plots()

import warnings
def print_results(fit, prior, data):
    # remove warnings caused by absence of Vno
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        etas_print_results(fit, prior, data)

# reuse code from etas-Ds.py
import importlib
etas_Ds = importlib.import_module('etas-Ds')
make_data = etas_Ds.make_data
make_models = etas_Ds.make_models
make_prior = etas_Ds.make_prior
etas_print_results = etas_Ds.print_results


if __name__ == '__main__':
    main()

"""
              metas: 0.41622(12)
                mDs: 1.20162(19)
           mDso-mDs: 0.264(16)
                Vnn: 0.76634(98)
                Vno: 0.0(1.0)

Partial % Errors:
                  metas       mDs  mDso-mDs       Vnn       Vno
---------------------------------------------------------------
  statistics:      0.03      0.02      4.52      0.12       nan
 log(etas:a):      0.00      0.00      0.00      0.00       nan
log(etas:dE):      0.00      0.00      0.00      0.00       nan
   log(Ds:a):      0.00      0.00      0.33      0.00       nan
  log(Ds:dE):      0.00      0.00      0.20      0.00       nan
  log(Dso:a):      0.00      0.00      2.35      0.00       nan
 log(Dso:dE):      0.00      0.00      2.92      0.00       nan
         Vnn:      0.00      0.00      0.01      0.04       nan
         Vno:      0.00      0.00      0.00      0.02       inf
---------------------------------------------------------------
       total:      0.03      0.02      5.88      0.13       inf
"""