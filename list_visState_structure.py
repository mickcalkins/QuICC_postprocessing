""" this file will list the structure of a QuICC 'visState' file
"""

from __future__ import print_function # this is for python 3 compatability

import numpy as np
import h5py

filename = 'visState0001.hdf5'

f = h5py.File(filename, 'r')

print("Here is the structure of file", filename)

def printname(name):
    print(name)

f.visit(printname)


    #Q = np.float( np.array( f.get('physical/chandrasekhar') ) )
print( np.array( f.get('run/time') ) )

