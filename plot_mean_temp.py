''' 

This file plots the mean temperature data from a QuICC visualization file.

'''

import numpy as np
import matplotlib.pyplot as plt
from read_visState import read_params
from read_visState import read_grid
from read_visState import read_data


# provide the name of the file to be read
filename = 'visState0000'

# import relevant parameters/quantities/arrays
Q, Ra, Pm, Pr = read_params(filename)
x, y, z = read_grid(filename)
temp, mean_temp, u, v, w, vort_x, vort_y, vort_z, bx, by, bz = read_data(filename)

# take the horizontal average and plot data 
data = mean_temp 
data_z = np.squeeze(np.squeeze(data[:, 0, 0]))

t_mean = np.flipud( 0.5*(data_z + 1. - z) )

zp = 0.5*(1. - z)

plt.rcParams["mathtext.fontset"] = "cm"
#plt.rc('font',family='Times')
plt.rcParams.update({'font.size': 16})
#plt.plot(t_mean, zp, 'ko-')
plt.plot(t_mean, zp, 'k-')

plt.xlabel(r'$\overline{T}$')
plt.ylabel(r'$z$', rotation=0)

plt.show()
#plt.savefig("mean_temperature.eps", dpi=150)
