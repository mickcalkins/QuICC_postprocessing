""" 

    This file plots quicc simulation data that is helpful for checking resolution. 

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import read_visState 
import compute_stats

# provide the name of the file to be read - no need for .hdf5
filename = 'visState0210'

# import relevant parameters/quantities/arrays - model specific
Ra, Pr = read_visState.read_params_RBC(filename)
#print('Rayleigh number = ', '{:.2e}'.format(Ra))
#print('Prandtl number = ', '{:.2e}'.format(Pr))


# grab the data
x, y, z = read_visState.read_grid(filename)
#temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz = read_visState.read_data_MC(filename)
temp, mean_temp, ux, uy, uz, vortx, vorty, vortz = read_visState.read_data_RBC(filename)


# rms of data - temperature and x-velocity
temp_ave, temp_rms = compute_stats.stats(temp)
ux_ave, ux_rms = compute_stats.stats(ux)

# rescaled z-grid
zp = 0.5*(1. - z)


# grab midpoints in z and x
midz = int((z.size)/2)
midx = int((x.size)/2)
midy = int((y.size)/2)

# define indices for extraction
bl_point = int(6)    # z-point within the boundary layer
x_point = midx  # x-point
y_point = midy  # y-point

# extract data at specified z/x locations
data = temp
data_midz = np.squeeze( data[midz, :, :] )
data_blz = np.squeeze( data[bl_point, :, :] )
spec_midz_y = np.fft.fft( np.squeeze( data_midz[x_point, :] ) )
spec_blz_y = np.fft.fft( np.squeeze( data_blz[x_point, :] ) )
spec_midz_x = np.fft.fft( np.squeeze( data_midz[:, y_point] ) )
spec_blz_x = np.fft.fft( np.squeeze( data_blz[:, y_point] ) )

#data = ux
#data_midz = np.squeeze( data[midz, :, :] )
#data_blz = np.squeeze( data[bl_point, :, :] )
#ux_spec_midz = np.fft.fft( np.squeeze( data_midz[x_point, :] ) )
#ux_spec_blz = np.fft.fft( np.squeeze( data_blz[x_point, :] ) )


# remove the aliased modes
de_size = spec_blz_x.size
true_size_x = np.int((de_size/2.)/1.5)

de_size = spec_blz_y.size
true_size_y = np.int((de_size/2.)/1.5)

#fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(12,6))
#fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(10,5))
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8))

axs[0, 0].plot(temp_rms, zp, 'ko-')
axs[0, 0].set_xlabel(r'$\theta_{rms}$')
axs[0, 0].set_ylabel(r'$z$')

axs[0, 1].plot(ux_rms, zp, 'ko-')
axs[0, 1].set_xlabel(r'$\theta_{rms}$')
#axs[0, 1].set_xlabel(r'$u_{rms}$')
axs[0, 1].set_ylabel(r'$z$')

axs[1, 0].plot(np.abs(spec_blz_x[0:true_size_x]), label=r'$|\widehat{\theta}_{bl}|$')
axs[1, 0].plot(np.abs(spec_midz_x[0:true_size_x]), label=r'$|\widehat{\theta}_{mid}|$')
axs[1, 0].set_xlabel(r'$N_x$')
axs[1, 0].set_ylabel(r'$|\widehat{\theta}|$')
axs[1, 0].set_xscale('log')
axs[1, 0].set_yscale('log')

#axs[1, 1].plot(np.abs(ux_spec_blz[0:true_size]), label=r'$|\widehat{u}_{bl}|$')
#axs[1, 1].plot(np.abs(ux_spec_midz[0:true_size]), label=r'$|\widehat{u}_{mid}|$')
axs[1, 1].plot(np.abs(spec_blz_y[0:true_size_y]), label=r'$|\widehat{\theta}_{bl}|$')
axs[1, 1].plot(np.abs(spec_midz_y[0:true_size_y]), label=r'$|\widehat{\theta}_{mid}|$')
axs[1, 1].set_xlabel(r'$N_y$')
axs[1, 1].set_ylabel(r'$|\widehat{\theta}|$')
#axs[1, 1].set_ylabel(r'$|\widehat{u}|$')
axs[1, 1].set_xscale('log')
axs[1, 1].set_yscale('log')


plt.savefig('resolution_check.png')

