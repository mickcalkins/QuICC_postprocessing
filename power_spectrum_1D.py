""" 

    This script computes a 1D power spectrum from QuICC data.

"""
import sys
import pandas as pd
#from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import read_visState 

# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000.hdf5'

# import relevant parameters/quantities/arrays
#Q, Ra, Pm, Pr = read_params(filename)
#Q, Ra, Pm, Pr = read_visState.read_params_MC(filename)
Q, Ra = read_visState.read_params_HMC(filename)
#Ra, Pr = read_visState.read_params_RBC(filename)

print('Reading the data...')
x, y, z_c = read_visState.read_grid(filename)
z = 0.5*(1.0 - z_c)
#temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz = read_visState.read_data_MC(filename)
temp, mean_temp, ux, uy, uz, vortx, vorty, vortz = read_visState.read_data_RBC(filename)

data = uz

# get grid size
Nz, Nx, Ny = data.shape
print('Grid size(Nx,Ny,Nz) =', Nx, Ny, Nz)

# midpoints in z and x
midz = int((z.size)/2)

data_spec_midz = np.fft.fft2( data[midz, :, :] )

data_spec = np.zeros(data.shape, dtype=complex)
power = np.zeros(data.shape, dtype=complex)

for i in range(Nz):
    data_spec[i,:,:] = np.fft.fft2( data[i, :, :] )
    power[i,:,:] = 0.5*data_spec[i,:,:]*np.conj( data_spec[i,:,:] )

power_ave = np.trapz(power, z, axis=0)
power_ave_1D = np.real(np.sum(power_ave, axis=1))

power_midz = 0.5*data_spec_midz*np.conj(data_spec_midz)
power_midz_1D = np.real(np.sum(power_midz, axis=1))

# remove the aliased modes
true_size_x = int(Nx/3)

#de_size_y = data_spec_blz_y.size
#true_size_y = np.int((de_size_y/2.)/1.5)
#de_size_x = data_spec_blz_x.size
#true_size_x = np.int((de_size_x/2.)/1.5)

#k_y = np.arange(0, true_size_y)
k_x = np.arange(0, true_size_x-1)

# write the data to file
dat = {'k_x': k_x, 'Power': power_ave_1D[0:true_size_x-1]}
df = pd.DataFrame(data=dat)
df.to_excel('Power_1D' + '_Q' + '{:.1e}'.format(Q) + '_Ra' + '{:.1e}'.format(Ra) + '.xlsx')


# make the plots
fig,ax = plt.subplots(1,1)
plt.rcParams["mathtext.fontset"] = "cm"


plt.figure(1)
#plt.loglog(k_x+1, np.abs(data_spec_blz_x[0:true_size_x]), label=r'$|\widehat{f}_{bl}|$')
#plt.loglog(k_x+1, power_midz_1D[0:true_size_x-1])
plt.loglog(k_x+1, power_ave_1D[0:true_size_x-1], 'r')
plt.xlabel(r'$k_x+1$', fontsize=18)
#ax.set_ylim(1e-20, 2*max_x)
#legend = plt.legend(loc='best', shadow=False, ncol = 2, fontsize = 18, frameon=False)

plt.tight_layout()
plt.show()
#plt.savefig('spectra_x.png')


