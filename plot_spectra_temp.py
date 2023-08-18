""" 

    This file plots horizontal spectra from a QuICC visualization file 

"""
import sys
#from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import read_visState 

# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000.hdf5'

# import relevant parameters/quantities/arrays
#Q, Ra, Pm, Pr = read_params(filename)
#Q, Ra, Pm, Pr = read_visState.read_params_MC(filename)
Ra, Pr = read_visState.read_params_RBC(filename)

# choose data
#entry = input("Enter the quantity to be plotted (temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz):  ")
#entry = input("Enter the quantity to be plotted (temp, mean_temp, ux, uy, uz, vortx, vorty, vortz):  ")
entry = 'temp'

x, y, z = read_visState.read_grid(filename)
temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz = read_visState.read_data_MC(filename)
#temp, mean_temp, ux, uy, uz, vortx, vorty, vortz = read_visState.read_data_RBC(filename)
data = temp

# midpoints in z and x
midz = int((z.size)/2)
midx = int((x.size)/2)
midy = int((y.size)/2)

# define indices for extraction
bl_point = int(10)    # z-point within the boundary layer
x_point = midx  # x-point
y_point = midy  # y-point

# extract data at specified z/x locations
data_midz = np.squeeze( data[midz, :, :] )
data_blz = np.squeeze( data[bl_point, :, :] )
data_spec_midz_y = np.fft.fft( np.squeeze( data_midz[x_point, :] ) )
data_spec_blz_y = np.fft.fft( np.squeeze( data_blz[x_point, :] ) )
data_spec_midz_x = np.fft.fft( np.squeeze( data_midz[:, y_point] ) )
data_spec_blz_x = np.fft.fft( np.squeeze( data_blz[:, y_point] ) )


# remove the aliased modes
de_size_y = data_spec_blz_y.size
true_size_y = int((de_size_y/2.)/1.5)
de_size_x = data_spec_blz_x.size
true_size_x = int((de_size_x/2.)/1.5)

k_y = np.arange(0, true_size_y)
k_x = np.arange(0, true_size_x)

max_blz_x = max(np.abs(data_spec_blz_x[0:true_size_x]))
max_midz_x = max(np.abs(data_spec_midz_x[0:true_size_x]))
max_x = min(max_blz_x, max_midz_x)
max_x_plot = max(max_blz_x, max_midz_x)

max_blz_y = max(np.abs(data_spec_blz_y[0:true_size_y]))
max_midz_y = max(np.abs(data_spec_midz_y[0:true_size_y]))
max_y = min(max_blz_y, max_midz_y)
max_y_plot = max(max_blz_y, max_midz_y)

# make the plots
fig1,ax1 = plt.subplots(1,1)
plt.rcParams["mathtext.fontset"] = "cm"

#ax1.figure(1)
ax1.loglog(k_x+1, np.abs(data_spec_blz_x[0:true_size_x]), label=r'$|\widehat{f}_{bl}|$')
ax1.loglog(k_x+1, np.abs(data_spec_midz_x[0:true_size_x]), label=r'$|\widehat{f}_{mid}|$')
plt.xlabel(r'$k_x+1$', fontsize=18)
#ax1.set_ylim(1e-15*max_x, 2*max_x_plot)
ax1.set_ylim(1e-17, 2*max_x_plot)
legend = plt.legend(loc='best', shadow=False, ncol = 2, fontsize = 18, frameon=False)

ax1.tick_params(labelright=True, right=True)

ax1.axhline(y = 1e-3*max_x, color = 'r', linestyle = '--')

plt.tight_layout()
plt.savefig('spectra_x.png')

fig2,ax2 = plt.subplots(1,1)
plt.rcParams["mathtext.fontset"] = "cm"

#ax2.figure(2)
ax2.loglog(k_y+1, np.abs(data_spec_blz_y[0:true_size_y]),label=r'$|\widehat{f}_{bl}|$')
ax2.loglog(k_y+1, np.abs(data_spec_midz_y[0:true_size_y]),label=r'$|\widehat{f}_{mid}|$')
plt.xlabel(r'$k_y+1$', fontsize=18)
ax2.set_ylim(1e-6*max_y, 2*max_y_plot)
legend = plt.legend(loc='best', shadow=False, ncol = 2, fontsize = 18, frameon=False)

ax2.tick_params(labelright=True, right=True)

ax2.axhline(y = 1e-3*max_y, color = 'r', linestyle = '--')

plt.tight_layout()
plt.savefig('spectra_y.png')

#plt.show()

