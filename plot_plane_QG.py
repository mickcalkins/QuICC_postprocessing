""" 

    This file plots data from a QuICC visualization file 

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from read_visState import read_params
from read_visState import read_grid
from read_visState import read_data_QG


# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000'

# import relevant parameters/quantities/arrays
Q, Ra, Pm, Pr = read_params(filename)
#Ra, Pr, Pm = read_nondims(filename)
#print 'Rayleigh number = ', '{:.2e}'.format(Ra)
#print 'Prandtl number = ', '{:.2e}'.format(Pr)
#print 'Magnetic Prandtl number = ', '{:.2e}'.format(Pm)

# choose data
data_entry = raw_input("Enter the quantity to be plotted (temp, ux, uy, uz, vortz):  ")

print 'You entered the following:  ', data_entry

# choose plane
plane_entry = raw_input("Enter the plane to be plotted (xy, xz):  ")

print 'You entered the following:  ', plane_entry

if not plane_entry in ['xy', 'xz']:
    print 'You did not enter a valid plane, exiting program'
    sys.exit() 

print 'Reading the data...'
x, y, z = read_grid(filename)
temp, ux, uy, uz, vortz = read_data_QG(filename)

if not data_entry in ['temp', 'ux', 'uy', 'uz', 'vortz']:
    print 'You did not enter a valid option, exiting program'
    sys.exit() 
elif data_entry == 'temp':
    data = temp
elif data_entry == 'ux':
    data = ux 
elif data_entry == 'uy':
    data = uy 
elif data_entry == 'uz':
    data = uz 
elif data_entry == 'vortz':
    data = vortz 

# get grid size
print 'Grid size =', data.shape

# create mesh
if plane_entry == 'xy':
   dim1, dim2 = np.meshgrid(x, y)
   #zp = (z.size)/2
   #zp = (z.size)/4
   #zp = 3*(z.size)/4
   #zp = 139
   zp = 1
   data_slice = np.squeeze(data[zp, :, :])
elif plane_entry == 'xz':
   dim1,dim2 = np.meshgrid(x, z)
   midy = (y.size)/2
   data_slice = np.squeeze(data[:, midy, :])

#fig, ax = plt.subplots(nrows=1)
#fig, (ax0, ax1) = plt.subplots(ncols=2)

# plot data in the midplane and within thermal boundary layer

#bl_point = 10
#data_slice_bl = np.squeeze(data[bl_point:bl_point+1, :, :])


print 'Making the figure...'

#midplane slice 
#c0 = ax0.pcolormesh(X, Y, data_slice_mid, shading='gouraud')
#ax0.axis('equal')
#fig.colorbar(c0, ax=ax0)
#ax0.set_axis_off()

#c1= ax1.pcolormesh(X, Y, data_slice_bl, shading='gouraud')
#ax1.axis('equal')
#fig.colorbar(c1, ax=ax1)

#ax1.set_axis_off()

plt.figure(1)

cmap = plt.get_cmap('PiYG')
#norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

plt.pcolormesh(dim1, dim2, data_slice, cmap=cmap, shading='gouraud')
#plt.pcolormesh(X, Y, data_slice_bl, shading='gouraud')
plt.axis('equal')
plt.axis('off')

plt.show()
#plt.savefig('temperature.png')

