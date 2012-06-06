#!/usr/bin/env python
# encoding: utf-8
"""
test-1.py  --- standard usage (including error budgets and plots)

Created by Peter Lepage on 2010-11-26.
Copyright (c) 2010/2011 Cornell University. All rights reserved.
"""

import os
from corrfitter import Corr2,Corr3,CorrFitter
from gvar import gvar,log,exp,BufferDict
from gvar.dataset import Dataset,avg_data
from numpy import array,arange
import lsqfit
##
lsqfit.nonlinear_fit.fmt_parameter = '%7.3f +- %7.3f'

DISPLAYPLOTS = True         # display plots at end of fitting
try: 
    import matplotlib
except ImportError:
    DISPLAYPLOTS = False

TEST = True         # testing mode?

def main():
    dfile = "coarse_3pt_etas_Ds.bint8"   # data file
    data = avg_data(Dataset(dfile)) # compute avg and cov for data
    pfile = "test-1.p"                   # last fit stored in file test-1.p
    fitter = CorrFitter(models=build_models())
    nexp_list = NEXP_TEST if TEST else [2,3,4,5,6]
    p0 = P0_TEST if TEST else pfile
    for nexp in nexp_list:
        print '========================== nexp =',nexp
        prior = build_prior(nexp)
        fit = fitter.lsqfit(data=data,prior=prior,p0=p0) # pfile)
        print_results(fit,prior,data)
        print '\n\n'
    if DISPLAYPLOTS:
        fitter.display_plots()
##
      
def print_results(fit,prior,data):
    """ print out additional results from the fit """
    ## etas parameters ##
    dEetas = exp(fit.p['log(etas:dE)'])
    print 'etas:dE =',fmtlist(dEetas[:3])
    print ' etas:E =',fmtlist([sum(dEetas[:i+1]) for i in range(3)])
    aetas = exp(fit.p['log(etas:a)'])
    print ' etas:a =',fmtlist(aetas[:3])
    print
    ##
    ## Ds parameters ##
    dEDs = exp(fit.p['log(Ds:dE)'])
    print 'Ds:dE =',fmtlist(dEDs[:3])
    print ' Ds:E =',fmtlist([sum(dEDs[:i+1]) for i in range(3)])
    aDs = exp(fit.p['log(Ds:a)'])
    print ' Ds:a =',fmtlist(aDs[:3])
    print
    ##
    ## vertices ##
    print 'etas->V->Ds  =',fit.p['Vnn'][0,0].fmt()
    print 'etas->V->Dso =',fit.p['Vno'][0,0].fmt()
    print
    ##
    ## error budget ##
    outputs = dict(Vnn=fit.p['Vnn'][0,0],Vno=fit.p['Vno'][0,0],
                 Eetas=dEetas[0],EDs=dEDs[0])
    inputs = {'stat.':[data[k] for k in data]} # statistical errors in data
    inputs.update(prior)                       # all entries in prior
    print fit.fmt_partialsdev(outputs,inputs,ndigit=3)
    ##
##

def fmtlist(x):
    """ Make neat lists for printing. """
    return '  '.join([xi.fmt() for xi in x])
##
    
def build_prior(nexp):
    """ build prior """
    prior = BufferDict()
    ## etas ##
    prior.add('log(etas:a)',[log(gvar(0.3,0.3)) for i in range(nexp)])
    prior.add('log(etas:dE)',[log(gvar(0.5,0.5)) for i in range(nexp)])
    prior['log(etas:dE)'][0] = log(gvar(0.4,0.2))
    ##
    ## Ds ##
    prior.add('log(Ds:a)',[log(gvar(0.3,0.3)) for i in range(nexp)])
    prior.add('log(Ds:dE)',[log(gvar(0.5,0.5)) for i in range(nexp)])
    prior['log(Ds:dE)'][0] = log(gvar(1.2,0.2))
    prior.add('log(Ds:ao)',[log(gvar(0.1,0.1)) for i in range(nexp)])
    prior.add('log(Ds:dEo)',[log(gvar(0.5,0.5)) for i in range(nexp)])
    prior['log(Ds:dEo)'][0] = log(exp(prior['log(Ds:dE)'][0])+gvar(0.3,0.3))
    ##
    ## V ##
    prior.add('Vnn',[[gvar(0.01,1.0) for i in range(nexp)] for j in range(nexp)])
    prior.add('Vno',[[gvar(0.01,1.0) for i in range(nexp)] for j in range(nexp)])
    return prior
##

def build_models():
    """ build models """
    tmin = 5
    tdata = range(64)
    tfit = range(tmin,64-tmin) # all ts
    
    tp = 64 # periodic
    models = [
        Corr2(datatag='etas',tp=tp,tdata=tdata,tfit=tfit,
            a='log(etas:a)',b='log(etas:a)',dE='log(etas:dE)'),
            
        Corr2(datatag='Ds',tp=tp,tdata=tdata,tfit=tfit,
            a=('log(Ds:a)','log(Ds:ao)'),b=('log(Ds:a)','log(Ds:ao)'),
            dE=('log(Ds:dE)','log(Ds:dEo)'),s=(1.,-1.)),

        Corr3(datatag='3ptT15',tdata=range(16),T=15,tfit=range(tmin,16-tmin),
            a='log(etas:a)',dEa='log(etas:dE)',tpa=tp,
            b=('log(Ds:a)','log(Ds:ao)'),dEb=('log(Ds:dE)','log(Ds:dEo)'),tpb=tp,sb=(1,-1.),
            Vnn='Vnn',Vno='Vno',transpose_V=False),
            
        Corr3(datatag='3ptT16',tdata=range(17),T=16,tfit=range(tmin,17-tmin),
            a='log(etas:a)',dEa='log(etas:dE)',tpa=tp,
            b=('log(Ds:a)','log(Ds:ao)'),dEb=('log(Ds:dE)','log(Ds:dEo)'),tpb=tp,sb=(1,-1.),
            Vnn='Vnn',Vno='Vno',transpose_V=False)
    ]
    return models
##    

NEXP_TEST = [6]
P0_TEST = {'log(Ds:a)': [-1.5382029985053618, -1.2890389183673718,
-0.8266404705754684, -1.1281693592882478, -1.1942970642079684,
-1.2029575678081768], 'log(Ds:ao)': [-2.7187514959968135, -2.4930018835083185,
-2.153056999309124, -2.274349119294338, -2.3000607164691345,
-2.3023830129920464], 'Vno': [[-0.7767194750934922, 0.29893034835380455,
0.0018184442826528018, 0.021717846191188498, 0.011730183184535546,
0.010175999805737345], [0.21499388626034838, 0.08090442617659241,
0.004526958942468914, 0.009566576678925012, 0.009981460666758885,
0.009999494947514592], [-0.03929180298702926, 0.016001668510532977,
0.010207619966221364, 0.009993775842607086, 0.009999490734868723,
0.009999977541285848], [-0.0008212538513334825, 0.009446705426502066,
0.010011631248341204, 0.010000150342742585, 0.009999996724907924,
0.009999999715913033], [0.008707361530432445, 0.009892605564266226,
0.009999973201624102, 0.010000008998777822, 0.01000000009622645,
0.009999999998091305], [0.009873335420234476, 0.00998801754685697,
0.009999941109021849, 0.010000000137364686, 0.010000000005831669,
0.010000000000057843]], 'Vnn': [[0.7674463006011443, -0.5085825363348335,
0.45476872843648597, 0.07910330776006197, 0.01672199952261604,
0.010575060501324804], [0.0428274543801929, 0.2276246794008917,
0.024339418136236402, 0.011105918132437442, 0.010079750081074862,
0.01000540600309465], [-0.0789057702177111, 0.039088907022644885,
0.010402246249041023, 0.010015038267077255, 0.010001273028618115,
0.010000092479777158], [-0.04722804062576639, 0.01216744637810462,
0.01003297837759048, 0.010000246006116778, 0.010000007960775022,
0.010000000698523368], [0.00035115678498912886, 0.010139057728244761,
0.01000276274869412, 0.010000023410840267, 0.010000000149913561,
0.01000000000466716], [0.008870833254729392, 0.01000779450180237,
0.010000194918187728, 0.010000001937760121, 0.01000000001438937,
0.01000000000008824]], 'log(Ds:dEo)': [0.36852247027597856,
-1.450434416125243, -0.716397528343509, -0.7262479860745519,
-0.6965486249012215, -0.6934059166295101], 'log(etas:dE)':
[-0.8765386251062858, -0.44666827867143977, -0.7029624277436295,
-0.6648624116442549, -0.6824061676451147, -0.6916868235591809], 'log(Ds:dE)':
[0.18374400617661235, -0.6825021780114466, -0.669942670009535,
-0.7649157041649213, -0.7045972594319748, -0.6943318698851685], 'log(etas:a)':
[-1.5213670965733321, -1.5776623136419392, -1.2117105815007296,
-1.2436876034918443, -1.213887220332993, -1.205310268830514]}

# for k in P0_TEST:
#     P0_TEST[k] = array(P0_TEST[k])

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
        stats.print_stats(140)

