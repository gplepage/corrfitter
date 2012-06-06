from corrfitter import Corr2,CorrFitter
from lsqfit import wavg
from gvar import gvar,log,exp,evalcov,mean,sdev,BufferDict
from gvar.dataset import Dataset,avg_data
import lsqfit
import gvar as gd
import time
import numpy
from numpy import linalg

lsqfit.nonlinear_fit.fmt_parameter = '%7.3f +- %7.3f'

NTERM_PRIOR = 8             # number of terms in prior
MC = None
RATIO = True

USE_MARGINALIZATION = False  # run mode
    
DISPLAY_PLOTS = False       # display plots at end
try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

def main():
    dfile = 'ups.bin24'     # input data file
    svdcut = 1e-3          # needed even without marginalization
    svdnum = 113            # number of samples in ups.bin24
    data = avg_data(Dataset(dfile))
    # for k in data:
    #     print k,gd.is_independent(data[k])[0]
    # c = gd.evalcov(data)
    models = make_models()
    if USE_MARGINALIZATION:
        prior = make_prior(NTERM_PRIOR)
    starttime = time.clock()
    p0 = None
    ntermlist = [4] if USE_MARGINALIZATION else [8]
    for nterm in ntermlist:
        if not USE_MARGINALIZATION:
            prior = make_prior(nterm)
        fitter = CorrFitter(models=models,svdcut=svdcut,
                            nterm=nterm if USE_MARGINALIZATION else None,
                            mc=MC,ratio=RATIO,maxit=10000)
        print 30*'=','nterm =',nterm,'  nterm_prior =',len(prior['logdE']),
        print '   MC =',MC,RATIO
        fit = fitter.lsqfit(prior=prior,data=data,p0=p0,svdnum=svdnum)
        # print fit
        print_results(fitter,prior,data)
        if nterm==1:
            print_effmass(fitter,prior,data)
        if not isinstance(p0,str):      # not using a file?
            print '\nusing last fit\n'
            p0 = fit.pmean
        print '\nt=%.1f ====' % (time.clock()-starttime)
        print
    if DISPLAY_PLOTS:
        fitter.display_plots()
##

def print_effmass(fitter,prior,data):
    meff = []
    xdata = fitter.fit.data[1]
    print '\neff mass:'
    for m in fitter.models:
        G = xdata[m.datatag]
        meff.append(wavg(log(G[:-1]/G[1:])))
        print '%10s:'%m.datatag,meff[-1],'  chi2/dof = %.2f  Q = %.2f'%(wavg.chi2/wavg.dof,wavg.Q)
    print '\n%10s:'%'all not d',wavg(meff[:-4]),'  chi2/dof = %.2f  Q = %.2f'%(wavg.chi2/wavg.dof,wavg.Q)
##

def print_results(fitter,prior,data):
    fit = fitter.fit
    p = fit.p
    dE = exp(p['logdE'])
    E = [sum(dE[:i+1])-(dE[0] if i!=0 else 0.0) for i in range(len(dE))]
    print '   E =',fmtlist(E),'  (En-E0 for n>0)'
    if 'logl' in p:
        print '   l =',fmtlist(exp(p['logl']))
    for k in ['g','d','e']:
        if k in p:
            print '   %s ='%k,fmtlist(p[k])
    print
    ## error budget ##
    inputs = dict(stat=[data[k] for k in data])
    inputs.update(prior)
    outputs = dict(E0=E[0])
    if len(E)>1:
        outputs['dE1'] = E[1]
    print fit.fmt_partialsdev(outputs,inputs,ndigit=3)
    ##
    print '\ncorrelators:\n   ',[m.datatag for m in fitter.models]
##
        
def make_prior(nterm):
    prior = BufferDict()
    c0 = 0.1
    dc = 1.0
    prior.add('logl',[log(gvar(0.4,0.4)) for i in range(nterm)])
    prior.add('g',numpy.array([gvar(c0,dc) for i in range(nterm)]))
    prior.add('e',[gvar(c0,dc) for i in range(nterm)])
    prior.add('d',[gvar(c0,dc) for i in range(nterm)])
    prior.add('logdE',log([gvar(0.5,0.5) for i in range(nterm)]))
    prior['logdE'][0] = log(gvar(.256,0.1))
    return prior
##

def make_models():
    # tfit = [1,2,3,5,8,13,21]    # fibonacci numbers
    tfit = range(1,24)          # full range
    models = [
    Corr2('1s0.ll',a='logl',b='logl',dE='logdE',tdata=range(1,24),tfit=tfit)
    ,
    Corr2('1s0.lg',a='logl',b='g',dE='logdE',
        tdata=range(1,24),tfit=tfit,othertags=['1s0.gl'])
    ,
    Corr2('1s0.gg',a='g',b='g',dE='logdE',tdata=range(1,24),tfit=tfit)
    ,
    Corr2('1s0.le',a='logl',b='e',dE='logdE',
        tdata=range(1,24),tfit=tfit,othertags=['1s0.el'])
    ,
    Corr2('1s0.ge',a='g',b='e',dE='logdE',
        tdata=range(1,24),tfit=tfit,othertags=['1s0.eg'])
    ,
    Corr2('1s0.ee',a='e',b='e',dE='logdE',tdata=range(1,24),tfit=tfit)
    ,
    Corr2('1s0.ld',a='logl',b='d',dE='logdE',
        tdata=range(1,24),tfit=tfit,othertags=['1s0.dl'])
    ,
    Corr2('1s0.gd',a='g',b='d',dE='logdE',
        tdata=range(1,24),tfit=tfit,othertags=['1s0.dg'])
    ,
    Corr2('1s0.ed',a='e',b='d',dE='logdE',
        tdata=range(1,24),tfit=tfit,othertags=['1s0.de'])
    ,
    Corr2('1s0.dd',a='d',b='d',dE='logdE',tdata=range(1,24),tfit=tfit)
    ]
    return models[:]
##

def fmtlist(x):
    """ Make neat lists for printing. """
    return '  '.join([xi.fmt() for xi in x])
##
    
if __name__ == '__main__':
    if True:
        main()
    else:
        import hotshot, hotshot.stats
        prof = hotshot.Profile("lsqfit-test1.prof")
        prof.runcall(main)
        prof.close()
        stats = hotshot.stats.load("lsqfit-test1.prof")
        stats.strip_dirs()
        stats.sort_stats('time', 'calls')
        stats.print_stats(40)
