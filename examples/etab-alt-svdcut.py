from __future__ import print_function   # makes this work for python2 and 3

import gvar as gv
import corrfitter as cf

def main():
    data, basis = make_data('etab.h5')
    s = gv.dataset.svd_diagnosis((data, 113), models=make_models())
    print('svdcut =', s.svdcut)
    s.plot_ratio(show=True)

import importlib
import sys
if sys.version_info > (2,):
    etab_alt = importlib.import_module('etab-alt')
else:
    etab_alt = importlib.__import__('etab-alt')
make_models = etab_alt.make_models
make_data = etab_alt.make_data

if __name__ == '__main__':
    main()
