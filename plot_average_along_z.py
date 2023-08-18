""" 

    This file produces an xy visualization of a quantity averaged along
    the direction of the imposed magnetic field.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from scipy import integrate
import read_visState 

# angle of tilt - hard coded for now
angle_deg = 45.0
angle = angle_deg*(np.pi/180.) # radians
eta3 = np.sin(angle) # z tilt factor

# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000'

# import relevant parameters/quantities/arrays
Q, Ra, Pm, Pr = read_visState.read_params_MC(filename)
#print 'Rayleigh number = ', '{:.2e}'.format(Ra)
#print 'Prandtl number = ', '{:.2e}'.format(Pr)
#print 'Magnetic Prandtl number = ', '{:.2e}'.format(Pm)
# choose data

entry = input("Enter the quantity to be plotted (temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz):  ")

print('You entered the following:  ', entry)


print('Reading the data...')
x, y, zc = read_visState.read_grid(filename)
z = 0.5*(1. - zc)
temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz = read_visState.read_data_MC(filename)

if not entry in ['temp', 'mean_temp', 'ux', 'uy', 'uz', 'vortx', 'vorty', 'vortz', 'bx', 'by', 'bz']:
    print('You did not enter a valid quantity option, exiting program')
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
Nz, Nx, Ny = data.shape
print('Grid size(Nx,Ny,Nz) =', Nx, Ny, Nz)

# create mesh for plot
dim1, dim2 = np.meshgrid(x, y)


# integrate user-specified variable in z
data = integrate.trapz(data[0:Nz, :, :], z, axis=0)


#data_bot = data_top

#sys.exit()

#fig, ax = plt.subplots(nrows=1)


print('Making the figure...')

fig, ax0 = plt.subplots(ncols=1)
num_bin = 1000

levels = MaxNLocator(nbins=num_bin).tick_values(data.min(), data.max())
cmap = plt.get_cmap('PiYG')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

ave = ax0.pcolormesh(x, y, data, cmap=cmap, norm=norm, shading='gouraud')
fig.colorbar(ave, ax=ax0, fraction=0.046, pad=0.04)
ax0.set(aspect='equal')

ax0.set_title('depth-average')

plt.tight_layout()


ax0.axis('off')

#plt.savefig('DepthAverage.png')

plt.show()

