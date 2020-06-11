from __future__ import print_function   # makes this work for python2 and 3



import gvar as gv
import corrfitter as cf

DISPLAYPLOTS = False         # display plots at end of fitting
try:
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

def main():
    data = make_data('etas-Ds.h5')
    models = make_models()                                              # 1a
    models = [
      models[0], models[1],                                             # 1b
      dict(nterm=(2, 1), svdcut=6.3e-5),                                # 1c
      (models[2], models[3])                                            # 1d
      ]
    fitter = cf.CorrFitter(models=models)                               # 1e
    p0 = None
    for N in [1, 2, 3, 4]:
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N)
        fit = fitter.chained_lsqfit(data=data, prior=prior, p0=p0)      # 2
        print(fit.format(pstyle=None if N < 4 else 'm'))
        p0 = fit.pmean
    print_results(fit, prior, data)
    if DISPLAYPLOTS:
        fit.show_plots()

    # check fit quality by adding noise
    print('\n==================== add svd, prior noise')
    noisy_fit = fitter.chained_lsqfit(
        data=data, prior=prior, p0=fit.pmean, svdcut=SVDCUT,
        noise=True,
        )
    print(noisy_fit.format(pstyle=None))
    p = key_parameters(fit.p)
    noisy_p = key_parameters(noisy_fit.p)
    print('      fit:', p)
    print('noisy fit:', noisy_p)
    print('          ', gv.fmt_chi2(gv.chi2(p - noisy_p)))

    # simulated fit
    for sim_pdata in fitter.simulated_pdata_iter(
        n=2, dataset=cf.read_dataset('etas-Ds.h5'), p_exact=fit.pmean
        ):
        print('\n==================== simulation')
        sim_fit = fitter.chained_lsqfit(
            pdata=sim_pdata, prior=prior, p0=fit.pmean, svdcut=SVDCUT,
            )
        print(sim_fit.format(pstyle=None))
        p = key_parameters(fit.pmean)
        sim_p = key_parameters(sim_fit.p)
        print('simulated - exact:', sim_p - p)
        print('          ', gv.fmt_chi2(gv.chi2(p - sim_p)))

import importlib
import sys
if sys.version_info > (2,):
    etas_Ds = importlib.import_module('etas-Ds')
else:
    etas_Ds = importlib.__import__('etas-Ds')
make_data = etas_Ds.make_data
make_models = etas_Ds.make_models
make_prior = etas_Ds.make_prior
print_results = etas_Ds.print_results
key_parameters = etas_Ds.key_parameters
SVDCUT = 1e-12 # etas_Ds.SVDCUT

# import importlib
# etas_Ds = importlib.import_module('etas-Ds')
# make_data = etas_Ds.make_data
# make_models = etas_Ds.make_models
# make_prior = etas_Ds.make_prior
# etas_print_results = etas_Ds.print_results


if __name__ == '__main__':
    gv.ranseed(123456)
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