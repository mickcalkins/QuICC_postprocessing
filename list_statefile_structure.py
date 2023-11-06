""" 

   This script will list the structure of a QuICC 'state' file

"""

from __future__ import print_function # this is for python 3 compatability

import numpy as np
import h5py

filename = 'state0003.hdf5'

f = h5py.File(filename, 'r')

print('Here is the structure of file', filename)

def printname(name):
    print(name)

f.visit(printname)

print('For file', filename, 'the time is', np.array( f.get('run/time') ) )

