""" 

    This file plots an rms profile from a QuICC visualization file.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import read_visState 
import compute_stats


# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000.hdf5'

# import relevant parameters/quantities/arrays - model specific
Ra, Pr = read_visState.read_params_RBC(filename)
#print('Rayleigh number = ', '{:.2e}'.format(Ra))
#print('Prandtl number = ', '{:.2e}'.format(Pr))

# choose data
#entry = input("Enter the quantity to be plotted (temp, mean_temp, ux, uy, uz, vortx, vorty, vortz):  ")
#entry = input("Enter the quantity to be plotted (temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz):  ")
entry = 'temp'

x, y, z = read_visState.read_grid(filename)
temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz = read_visState.read_data_MC(filename)
#temp, mean_temp, ux, uy, uz, vortx, vorty, vortz = read_visState.read_data_RBC(filename)
data = temp

# get grid size
Nz, Nx, Ny = data.shape
#print('Grid size(Nx,Ny,Nz) =', Nx, Ny, Nz)

# choose horizontal x/y location
x_loc = int((x.size)/2)
y_loc = int((y.size)/2)

#print('Generating plot...')

# data at particular (x_loc, y_loc)
data_z = np.squeeze(np.squeeze(data[:, x_loc, y_loc]))

# average of data - organized as (z,x,y)
#data_ave_x = np.average(data, axis=1)
#data_ave = np.average(data_ave_x, axis=1)
#data_ave_mat = np.transpose( np.tile(data_ave, (Ny, Nx, 1)) )


# rms of data
#data_no_mean = data - data_ave_mat
#data2 = np.square(data_no_mean)
#mx = np.average(data2, axis=1)
#mxy = np.average(mx, axis=1)
#rms = np.sqrt(mxy)

data_ave, data_rms = compute_stats.stats(data)

zp = 0.5*(1. - z)

#plt.plot(data_ave, zp, 'ko-')
#plt.plot(data_z, zp, 'ko-')
plt.plot(data_rms, zp, 'ko-')
plt.xlabel(entry)
plt.ylabel(r'$z$', rotation=0)

savefile = 'rms_profile_' + entry + '.png'
plt.savefig(savefile)
#plt.show()

