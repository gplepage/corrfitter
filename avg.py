#! /usr/bin/env python
""" Average dataset files.

Usage:

    avg.py file1 file2 ... (for text files)

    avg.py file.h5 group1 group2 ... (for hdf5 file)

Returns a table showing averages and standard deviations for all
quantities in the dataset files. (See gvar.dataset.Dataset for
information on file formats.)
"""
from __future__ import print_function

import fileinput
import warnings
import sys

import gvar
from gvar.dataset import Dataset, avg_data

warnings.simplefilter("ignore")

if sys.argv and sys.argv[1][-3:] == '.h5':
    dset = Dataset(sys.argv[1], h5group=sys.argv[2:])
else:
    dset = Dataset(fileinput.FileInput(openhook=fileinput.hook_compressed))
ndset = Dataset()
# put number of items in key
for k in dset:
    ndset[str(k) + (' [{}]'.format(len(dset[k])))] = dset[k]
data = avg_data(ndset)
print(gvar.tabulate(data, headers=False))