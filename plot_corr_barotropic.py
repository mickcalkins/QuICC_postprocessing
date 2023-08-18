""" 

    This file plots the QG barotropic streamfunction or vorticity
    from a QuICC visualization file.

"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from scipy import integrate
#from read_visState import read_params
from read_visState import read_grid
from read_visState import read_data_QG

# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000'

# import relevant parameters/quantities/arrays
#Q, Ra, Pm, Pr = read_params(filename)
#Ra, Pr, Pm = read_nondims(filename)
#print 'Rayleigh number = ', '{:.2e}'.format(Ra)
#print 'Prandtl number = ', '{:.2e}'.format(Pr)
#print 'Magnetic Prandtl number = ', '{:.2e}'.format(Pm)

# choose data
#data_entry = input("Enter the quantity to be plotted (stream, vortz):  ")

#print('You entered the following:  ', data_entry)


print('Reading the data...')
x, y, z_c = read_grid(filename)
temp, ux, uy, uz, stream, vortz = read_data_QG(filename)

z = 0.5*(1.0 - z_c)

# get grid size
print('Grid size =', ux.shape)
Nz = ux.shape[0]
Nx = ux.shape[1]
Ny = ux.shape[2]


# create mesh for plot
dim1, dim2 = np.meshgrid(x, y)

# for depth-varying conductivity
#z_bot = z[0:len(z)/2]
#z_top = z[len(z)/2:]


# compute the barotropic mode for both u and vort 
ux_ave = integrate.trapz(ux[:, :, :], z, axis=0)
vort_ave = integrate.trapz(vortz[:, :, :], z, axis=0)

ux_p = np.zeros((Nz, Nx, Ny))
vort_p = np.zeros((Nz, Nx, Ny))
corr_p = np.zeros((Nz, Nx, Ny))


# compute fluctuating quantities
for i in range(len(z)):
    ux_p[i, :, :] = ux[i, :, :] - ux_ave[:, :] 
    vort_p[i, :, :] = vortz[i, :, :] - vort_ave[:, :] 
    corr_p[i, :, :] = ux_p[i, :, :]*vort_p[i, :, :]

# average the correlation function
corr_ave = integrate.trapz(corr_p[:, :, :], z, axis=0)

#data = ux_p[10, :, :]
data = corr_ave

#data_bot = data_top

#sys.exit()

#fig, ax = plt.subplots(nrows=1)


print('Making the figure...')

#fig, (ax0, ax1) = plt.subplots(ncols=2)
fig, ax = plt.subplots(ncols=1)
num_bin = 1000

levels = MaxNLocator(nbins=num_bin).tick_values(data.min(), data.max())
cmap = plt.get_cmap('PiYG')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

baro = ax.pcolormesh(x, y, data, cmap=cmap, norm=norm, shading='gouraud')
fig.colorbar(baro, ax=ax, fraction=0.046, pad=0.04)
ax.set(aspect='equal')

ax.set_title(r'$\langle u^{\prime} \zeta^{\prime} \rangle$')

plt.tight_layout()


ax.axis('off')
#ax0.axis('off')
#ax1.axis('off')

#if data_entry == 'vortz':
#   plt.savefig('VortBarotropic.png')
#elif data_entry == 'stream':

plt.savefig('Correlation_R80.png')

plt.show()

sys.exit()
