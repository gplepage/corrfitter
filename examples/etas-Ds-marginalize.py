from __future__ import print_function   # makes this work for python2 and 3

import gvar as gv
import corrfitter as cf

DISPLAYPLOTS = False         # display plots at end of fitting?
try:
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

def main():
    data = make_data('etas-Ds.h5')
    fitter = cf.CorrFitter(models=make_models())
    p0 = None
    prior = make_prior(8)                                               # 1
    for N in [1, 2]:                                                    # 2
        print(30 * '=', 'nterm =', N)
        fit = fitter.lsqfit(
            data=data, prior=prior, p0=p0, nterm=(N, N), svdcut=SVDCUT  # 3
            )
        print(fit)                                                      # 4
        p0 = fit.pmean
    print_results(fit, prior, data)
    if DISPLAYPLOTS:
        fit.show_plots()

    # check fit quality by adding noise
    print('\n==================== add svd, prior noise')
    noisy_fit = fitter.lsqfit(
        data=data, prior=prior, p0=fit.pmean, svdcut=SVDCUT, nterm=(N, N),
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
        sim_fit = fitter.lsqfit(
            pdata=sim_pdata, prior=prior, p0=fit.pmean, svdcut=SVDCUT,
            nterm=(N, N),
            )
        print(sim_fit.format(pstyle=None))
        p = key_parameters(fit.pmean)
        sim_p = key_parameters(sim_fit.p)
        print('simulated - exact:', sim_p - p)
        print('          ', gv.fmt_chi2(gv.chi2(p - sim_p)))

# reuse code from etas-Ds.py
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
SVDCUT = etas_Ds.SVDCUT

if __name__ == '__main__':
    gv.ranseed(123456)
    main()