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
    models = make_models()
    prior = make_prior(8)
    fitter = CorrFitter(models=make_models())
    p0 = None
    for N in [1, 2, 3, 4]:
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N)
        fit = fitter.chained_lsqfit(data=data, prior=prior, p0=p0)
        print(fit.formatall(pstyle=None if N < 4 else 'v'))
        p0 = fit.pmean
    print_results(fit, prior, data)
    if DISPLAYPLOTS:
            fitter.display_plots()
    # fit with structured models
    # print('--- refit with structured models')
    # models = [models[:2]] + models[2:]
    # fitter = CorrFitter(models=models)
    # fit = fitter.chained_lsqfit(data=data, prior=prior, p0=p0)
    # print_results(fit, prior, data)



# reuse code from etas-Ds.py
import importlib
etas_Ds = importlib.import_module('etas-Ds')
make_data = etas_Ds.make_data
make_models = etas_Ds.make_models
make_prior = etas_Ds.make_prior
print_results = etas_Ds.print_results


if __name__ == '__main__':
    main()
