#!/usr/bin/env python
# encoding: utf-8
"""
dataset.py --- simplified replacement for old module; for legacy purposes only 
               (use gvar.dataset for new stuff).
"""
# Created by G. Peter Lepage, Cornell University, on 2012-05-22.
# Copyright (c) 2010-2012 G. Peter Lepage.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version (see <http://www.gnu.org/licenses/>).
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import gvar
import numpy

def flexargfcn(f):
    """ decorator for Dataset analysis functions """
    def newf(self,arg=None):
        if arg is not None:
            if (hasattr(arg,"__hash__") and arg.__hash__ is not None 
                and arg in self): 
                return f(self,self[arg])
            else:
                ans = {}
                for k in arg:
                    ans[k] = f(self,self[k]) 
                return ans
        else:
            ans = {}
            for k in self:
                ans[k] = f(self,self[k]) 
            return ans
    ##
    newf.__doc__ = f.__doc__
    return newf
##

class Dataset(gvar.dataset.Dataset):
    """ legacy class --- don't use otherwise """
    def __init__(self, *args,**kargs):
        if 'bstrap' in kargs:
            bstrap = kargs['bstrap']
            del kargs['bstrap']
        else:
            bstrap = False
        super(Dataset, self).__init__(*args,**kargs)
        self.bstrap = bstrap
    ##
    def copy(self,dlist):
        if type(dlist) not in [list,tuple]:
            dlist = [dlist]
        for d in dlist:
            self.extend(d)
    ##
    def nmeas(self):
        dlen = [len(self[k]) for k in self]
        return min(dlen),max(dlen)
    ##
    @flexargfcn
    def avg(self,x):
        return numpy.mean(x,axis=0)
    ##
    @flexargfcn
    def median(self,x):
        return gvar.dataset.avg_data(x,median=True,spread=self.bstrap)
    ##
    @flexargfcn
    def sdev(self,x):
        if self.bstrap:
            return numpy.std(xx,axis=0)
        else:
            return numpy.std(xx,axis=0)/math.sqrt(len(xx))
    ##
    @flexargfcn
    def cov(self,d):
        oldshape = d[0].shape
        d = [di.flatten() for di in d]
        # d = d.reshape((d.shape[0],-1)) # flatten all other dimensions
        dcov = numpy.cov(d,rowvar=False,bias=True)
        if self.bstrap:
            return dcov.reshape(oldshape+oldshape)
        else:
            return dcov.reshape(oldshape+oldshape)/float(len(d))
    ##     
    def gdev(self,arg=None):
        if arg is None:
            arg = self.keys()
            strip_ans = False
        elif (hasattr(arg,'__hash__') and arg.__hash__ is not None 
              and arg in self):
            arg = [arg]
            strip_ans = True
        if len(arg)==0:
            return dict()
        dd = dict()
        for k in arg:
            dd[k] = self[k]
        ans = gvar.dataset.avg_data(dd,spread=self.bstrap)
        return ans[arg[0]] if strip_ans else ans
    ##
    def bin(self,nbin=2):
        return gvar.bin_data(self,binsize=nbin)
    ##
    def bootstrap_iter(self,n=None):
        return gvar.dataset.bootstrap_iter(self,n)
    ##
    def assemble(self,template,newtag=''):
        ans = Dataset()
        ans[newtag] = self.arrayzip(template)
        return ans
##  

def tabulate_avg(avgout,format=(6,3)):
    """ Tabulates averages and standard deviations.
        
    tabulate_avg(...) creates a nicely formatted table displaying the
    output from functions like ``dataset.Dataset.gdev``. Here ``avgout`` is
    the output. Parameter ``format`` specifies the output format:
    ``format=(N,D)`` implies that format ``'%N.Df(%Dd)'`` is used to print
    ``avg,int(10**D * std_dev)``. The table is returned as a single string,
    for printing.
    """
    table = []
    output = avgout.items()
    output.sort()
    for tag,avsd in output:
        try:
            av = avsd.mean
            sd = avsd.sdev
        except AttributeError:
            av = gvar.mean(avsd)
            sd = gvar.sdev(avsd)
        lines = ''
        line = '%15s' % str(tag)
        try:
            sdfac = 10**format[1]
            fmt = (' %'+str(format[0])+'.'+str(format[1])+
                  'f(%'+str(format[1])+'d)')
            def avgfmt(av,sd,fmt=fmt,sdfac=sdfac):
                try:
                    return fmt % (av,int(sdfac*sd+0.5))
                except:
                    return (' %g (%.4g)' % (av,sd))
            ##
        except:
            def avgfmt(av,sd):
                return (' %g (%.4g)' % (av,sd))
            ##
        na = len(av)
        if len(sd)<na:
            na = len(sd)
        if na>=1:
            for i in xrange(na):
                if len(sd.shape)==2:
                    sdi = math.sqrt(sd[i][i])
                else:
                    sdi = sd[i]
                nextfield = avgfmt(av[i],sdi)
                if (len(nextfield)+len(line))>78:
                    lines = lines + line + '\n'
                    line = ''.ljust(15) + nextfield
                else:
                    line = line + nextfield
            table.append(lines + line +'\n')
    return '\n'.join(table)
##    
    

if __name__ == '__main__':
    import gvar
    gvar.ranseed((1950,1))
    r1 = gvar.gvar(8.,1.)
    r2 = gvar.gvar([-10.,-9.],[2.,3.])
    r3 = gvar.gvar([[0.,1.],[2.,3.]],[[1.,2.],[3.,4.]])
    r3_iter = gvar.raniter(r3)
    r2_iter = gvar.raniter(r2)
    N = 1001
    d = Dataset(bstrap=False)
    for x in range(N):
        d.append('x',r3_iter.next())
    for x in range(N):
        d.append('y',r2_iter.next())
    d2 = Dataset(bstrap=False)
    for x in range(N):
        d2.append('z',r1())
    d.copy(d2)
    med = d.gdev()
    for k in med:
        print( k,med[k])
    avg = d.avg()
    for k in avg:
        print( k,avg)
    print( d.nmeas())
    nd = d.assemble(['y','y'],'yy')
    print( nd.avg())
    nd = d.grep('x|y')
    print( nd.keys())
        
    output = """    
    y [-10.0069 +- 0.00632473 -8.98159 +- 0.00949937]
    x [[-0.00311314 +- 0.00316463 1.00027 +- 0.00634238]
     [2.00083 +- 0.00950991 2.99893 +- 0.0126476]]
    z 8.00641 +- 0.00315755
    """     # 15.5 sec with N = 100001