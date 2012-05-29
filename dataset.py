#!/usr/bin/env python
# encoding: utf-8
"""
dataset.py --- simplified replacement for old module; for legacy purposes only 
               (use gvar.Dataset for new stuff)

Created by Peter Lepage on 2012-05-22.
Copyright (c) 2012 Cornell University. All rights reserved.
"""

import gvar

class Dataset(gvar.Dataset):
    """ legacy class --- don't use otherwise """
    def __init__(self, *args,**kargs):
        super(Dataset, self).__init__(*args,**kargs)
    ##
    def gdev(self):
        return gvar.avg_data(self)
    ##
    def bin(self,nbin=2):
        return gvar.bin_data(self,binsize=nbin)
    ##
##  

if __name__ == '__main__':
    a = Dataset()
    a.append(x=1,y=[10,100])
    a.append(x=2,y=[20,200])
    a.append(x=3,y=[30,300])
    print (a.bin(3))

