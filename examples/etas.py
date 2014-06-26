from __future__ import print_function   # makes this work for python2 and 3

import collections
import gvar as gv
import numpy as np
from corrfitter import CorrFitter, Corr2, read_dataset

def main():
    data = make_data(filename='etas-Ds.data')
    fitter = CorrFitter(models=make_models())
    p0 = None
    for N in range(2, 6):                                   
        print(30 * '=', 'nterm =', N)
        prior = make_prior(N)
        fit = fitter.lsqfit(data=data, prior=prior, p0=p0)  
        p0 = fit.pmean    
        print_results(fit)                                  

def make_data(filename):
    """ Read data, compute averages/covariance matrix for G(t). """
    return gv.dataset.avg_data(read_dataset(filename))

def make_models():
    """ Create corrfitter model for G(t). """
    corr = Corr2(
        datatag='etas', tp=64, tdata=range(64), tfit=range(5, 64-5), 
        a='a', b='a', dE='dE'
        )
    return [corr]

def make_prior(N):
    """ Create prior for N-state fit. """
    prior = collections.OrderedDict()
    prior['a'] = gv.gvar(N * ['0(1)'])
    prior['logdE'] = gv.log(gv.gvar(N * ['0.5(5)']))                       
    return prior

def print_results(fit):
        p = fit.transformed_p                               
        E = np.cumsum(p['dE'])                              
        a = p['a']                                          
        print('{:2}  {:15}  {:15}'.format('E', E[0], E[1])) 
        print('{:2}  {:15}  {:15}\n'.format('a', a[0], a[1]))

if __name__ == '__main__':
    main()