from __future__ import print_function   # makes this work for python2 and 3

import gvar as gv
import corrfitter as cf

def main():
    dset = cf.read_dataset('etab.h5', grep='1s0')
    s = gv.dataset.svd_diagnosis(dset, models=make_models())
    print('svdcut =', s.svdcut)
    s.plot_ratio(show=True)

from etab import make_models

if __name__ == '__main__':
    main()
