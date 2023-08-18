""" 

    This file plots a horizontal slice of data from a QuICC visualization file.

"""

import numpy as np
import matplotlib.pyplot as plt
from read_visState import read_params_RBC
from read_visState import read_grid
from read_visState import read_data_RBC


# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000'

# import relevant parameters/quantities/arrays
Ra, Pr = read_params_RBC(filename)
x, y, z = read_grid(filename)
#temp, mean_temp, u, v, w, vort_x, vort_y, vort_z, bx, by, bz = read_data(filename)
temp, mean_temp, u, v, w, vort_x, vort_y, vort_z = read_data_RBC(filename)

# create mesh
X, Y = np.meshgrid(x, y)

#fig, ax = plt.subplots(nrows=1)
#fig, (ax0, ax1) = plt.subplots(ncols=2)

# plot data in the midplane and within thermal boundary layer
data = temp 
midz = (z.size)/2

bl_point = 10
data_slice_mid = np.squeeze(data[midz:midz+1, :, :])
data_slice_bl = np.squeeze(data[bl_point:bl_point+1, :, :])

#c0 = ax0.pcolormesh(X, Y, data_slice_mid, shading='gouraud')
#c1= ax1.pcolormesh(X, Y, data_slice_bl, shading='gouraud')

# set the limits of the plot to the limits of the data
#ax0.axis('equal')
#ax1.axis('equal')
#fig.colorbar(c0, ax=ax0)
#fig.colorbar(c1, ax=ax1)

#ax0.set_axis_off()
#ax1.set_axis_off()

plt.figure(1)
plt.pcolormesh(X, Y, data_slice_mid, shading='gouraud')
plt.axis('equal')
plt.axis('off')

plt.show()
#plt.savefig('temperature.png')
