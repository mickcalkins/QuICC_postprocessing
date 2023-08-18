
'''This is a simple example that calls the read data subroutine'''

from __future__ import print_function # this is for python 3 compatability
from read_state import read_params
from read_state import read_truncation
from read_state import read_data

# provide the name of the file to be read
filename = 'state0001'

# import relevant parameters/quantities/arrays
Q, Ra, Pm, Pr = read_params(filename)
trunc_phys, trunc_spec, trunc_trans = read_truncation(filename)
print(trunc_phys, trunc_spec, trunc_trans)
#vel_pol, vel_tor = read_data(filename)


