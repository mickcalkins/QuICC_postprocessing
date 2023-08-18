""" 

    This file plots horizontal spectra from a QuICC visualization file 

"""
import sys
#from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import read_visState 

# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000'

# import relevant parameters/quantities/arrays
#Q, Ra, Pm, Pr = read_params(filename)
#Q, Ra, Pm, Pr = read_visState.read_params_MC(filename)
Ra, Pr = read_visState.read_params_RBC(filename)

# choose data
#entry = input("Enter the quantity to be plotted (temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz):  ")
entry = input("Enter the quantity to be plotted (temp, mean_temp, ux, uy, uz, vortx, vorty, vortz):  ")

print('You entered the following:  ', entry)

print('Reading the data...')
x, y, z = read_visState.read_grid(filename)
#temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz = read_visState.read_data_MC(filename)
temp, mean_temp, ux, uy, uz, vortx, vorty, vortz = read_visState.read_data_RBC(filename)

if not entry in ['temp', 'mean_temp', 'ux', 'uy', 'uz', 'vortx', 'vorty', 'vortz', 'bx', 'by', 'bz']:
    print('You did not enter a valid option, exiting program')
    sys.exit() 
elif entry == 'temp':
    data = temp
elif entry == 'mean_temp':
    data = mean_temp
elif entry == 'ux':
    data = ux 
elif entry == 'uy':
    data = uy 
elif entry == 'uz':
    data = uz 
elif entry == 'vortx':
    data = vortx 
elif entry == 'vorty':
    data = vorty 
elif entry == 'vortz':
    data = vortz 
elif entry == 'bx':
    data = bx 
elif entry == 'by':
    data = by 
elif entry == 'bz':
    data = bz 

# get grid size
print('Grid size =', data.shape)


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
true_size_y = np.int((de_size_y/2.)/1.5)
de_size_x = data_spec_blz_x.size
true_size_x = np.int((de_size_x/2.)/1.5)

k_y = np.arange(0, true_size_y)
k_x = np.arange(0, true_size_x)

# make the plots
fig,ax = plt.subplots(1,1)
plt.rcParams["mathtext.fontset"] = "cm"

max_blz_x = max(np.abs(data_spec_blz_x[0:true_size_x]))
max_midz_x = max(np.abs(data_spec_midz_x[0:true_size_x]))
max_x = max(max_blz_x, max_midz_x)

max_blz_y = max(np.abs(data_spec_blz_y[0:true_size_y]))
max_midz_y = max(np.abs(data_spec_midz_y[0:true_size_y]))
max_y = max(max_blz_y, max_midz_y)

plt.figure(1)
plt.loglog(k_x+1, np.abs(data_spec_blz_x[0:true_size_x]), label=r'$|\widehat{f}_{bl}|$')
plt.loglog(k_x+1, np.abs(data_spec_midz_x[0:true_size_x]), label=r'$|\widehat{f}_{mid}|$')
plt.xlabel(r'$k_x+1$', fontsize=18)
ax.set_ylim(1e-20, 2*max_x)
legend = plt.legend(loc='best', shadow=False, ncol = 2, fontsize = 18, frameon=False)

plt.tight_layout()
plt.savefig('spectra_x.png')

plt.figure(2)
plt.loglog(k_y+1, np.abs(data_spec_blz_y[0:true_size_y]),label=r'$|\widehat{f}_{bl}|$')
plt.loglog(k_y+1, np.abs(data_spec_midz_y[0:true_size_y]),label=r'$|\widehat{f}_{mid}|$')
plt.xlabel(r'$k_y+1$', fontsize=18)
ax.set_ylim(1e-17, 2*max_y)
legend = plt.legend(loc='best', shadow=False, ncol = 2, fontsize = 18, frameon=False)

plt.tight_layout()
plt.savefig('spectra_y.png')

#plt.show()

