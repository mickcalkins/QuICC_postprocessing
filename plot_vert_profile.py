""" 

    This file plots data from a QuICC visualization file 

"""

import numpy as np
import matplotlib.pyplot as plt
from read_visState import read_params
from read_visState import read_grid
from read_visState import read_data


# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000'

# import relevant parameters/quantities/arrays
Q, Ra, Pm, Pr = read_params(filename)
print 'Ra=',Ra
print 'Pm=',Pm
print 'Pr=',Pr
x, y, z = read_grid(filename)
temp, mean_temp, u, v, w, vort_x, vort_y, vort_z, bx, by, bz = read_data(filename)


# choose horizontal x/y location
x_loc = (x.size)/2
y_loc = (y.size)/2

# plot data 
data = temp
data_z = np.squeeze(np.squeeze(data[:, x_loc, y_loc]))

# rms of data
data2 = np.square(data)
mx = np.average(data2, axis=1)
mxy = np.average(mx, axis=1)
rms = np.sqrt(mxy)

zp = 0.5*(1. - z)

plt.plot(data_z, zp, 'ko-')

plt.show()

