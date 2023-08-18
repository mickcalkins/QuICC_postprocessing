""" 

    This file plots the QG barotropic streamfunction or vorticity
    from a QuICC visualization file.

"""

from __future__ import print_function

from six.moves import input
import sys, os
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
data_entry = input("Enter the quantity to be plotted (stream, vortz):  ")

print('You entered the following:  ', data_entry)


print('Reading the data...')
x, y, z_c = read_grid(filename)
temp, ux, uy, uz, stream, vortz = read_data_QG(filename)

z = 0.5*(1.0 - z_c)

if not data_entry in ['stream', 'vortz']:
    print('You did not enter a valid option, exiting program')
    sys.exit() 
elif data_entry == 'stream':
    data = stream 
elif data_entry == 'vortz':
    data = vortz 

# get grid size
print('Grid size =', data.shape)

# create mesh for plot
dim1, dim2 = np.meshgrid(x, y)

# for depth-varying conductivity
#z_bot = z[0:len(z)/2]
#z_top = z[len(z)/2:]



# integrate user-specified variable in z
#data_top = integrate.trapz(data[0:len(z)/2, :, :], z_top, axis=0)
#data_bot = integrate.trapz(data[len(z)/2:, :, :], z_bot, axis=0)
data = integrate.simps(data[:, :, :], z, axis=0)

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

#levels_bot = MaxNLocator(nbins=num_bin).tick_values(data_bot.min(), data_bot.max())
#cmap_bot = plt.get_cmap('PiYG')
#norm_bot = BoundaryNorm(levels_bot, ncolors=cmap_bot.N, clip=True)

#bot = ax0.pcolormesh(x, y, data_bot, cmap=cmap_bot, norm=norm_bot, shading='gouraud')
#fig.colorbar(bot, ax=ax0, fraction=0.046, pad=0.04)
#ax0.set(aspect='equal')

#levels_top = MaxNLocator(nbins=num_bin).tick_values(data_top.min(), data_top.max())
#cmap_top = plt.get_cmap('PiYG')
#norm_top = BoundaryNorm(levels_top, ncolors=cmap_top.N, clip=True)

#top = ax1.pcolormesh(x, y, data_top, cmap=cmap_top, norm=norm_top, shading='gouraud')
#fig.colorbar(top, ax=ax1, fraction=0.046, pad=0.04)
#ax1.set(aspect='equal')


ax.set_title('barotropic mode')
#ax0.set_title('bottom layer (conductor)')
#ax1.set_title('top layer (insulator)')

plt.tight_layout()


ax.axis('off')
#ax0.axis('off')
#ax1.axis('off')

if data_entry == 'vortz':
   plt.savefig('VortBarotropic.png')
elif data_entry == 'stream':
   plt.savefig('PsiBarotropic.png')


plt.show()

sys.exit()
