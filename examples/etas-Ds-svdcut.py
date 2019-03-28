from __future__ import print_function   # makes this work for python2 and 3

import gvar as gv
import corrfitter as cf

def main():
    dset = cf.read_dataset('etas-Ds.h5')
    s = gv.dataset.svd_diagnosis(dset, models=make_models())
    print('svdcut =', s.svdcut)
    s.plot_ratio(show=True)
    # chained fit
    models = make_models()
    models = models[:2] + [tuple(models[2:])]
    for m in models:
        s = gv.dataset.svd_diagnosis(dset, models=[m])
        print('svdcut (chained) =', s.svdcut)
        s.plot_ratio(show=True)

import importlib
import sys
if sys.version_info > (2,):
    make_models = importlib.import_module('etas-Ds').make_models
else:
    make_models = importlib.__import__('etas-Ds').make_models

if __name__ == '__main__':
    gv.ranseed(123456)
    main()
    # if True:
    #     main()
    # else:
    #     import cProfile, pstats, StringIO
    #     pr = cProfile.Profile()
    #     pr.enable()
    #     main()
    #     pr.disable()
    #     s = StringIO.StringIO()
    #     sortby = 'tottime'
    #     ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    #     ps.print_stats()
    #     print (s.getvalue())
