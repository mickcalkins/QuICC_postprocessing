""" 

    This file plots data from a QuICC visualization file 

"""

from __future__ import print_function

from six.moves import input
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from read_visState import read_data_Beta


# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000'

# import relevant parameters/quantities/arrays
#Q, Ra, Pm, Pr = read_params(filename)
#Ra, Pr, Pm = read_nondims(filename)
#print 'Rayleigh number = ', '{:.2e}'.format(Ra)
#print 'Prandtl number = ', '{:.2e}'.format(Pr)
#print 'Magnetic Prandtl number = ', '{:.2e}'.format(Pm)

# choose data
data_entry = input("Enter the quantity to be plotted (temp, ux, uy, uz, stream):  ")

print('You entered the following:  ', data_entry)

# choose plane
plane_entry = input("Enter the plane to be plotted (xy, xz, yz):  ")

print('You entered the following:  ', plane_entry)

if not plane_entry in ['xy', 'xz', 'yz']:
    print('You did not enter a valid plane, exiting program')
    sys.exit() 

print('Reading the data...')
Ra, x, y, z, temp, ux, uy, uz, stream = read_data_Beta(filename)

if not data_entry in ['temp', 'ux', 'uy', 'uz', 'stream']:
    print('You did not enter a valid option, exiting program')
    sys.exit() 
elif data_entry == 'temp':
    data = temp
elif data_entry == 'ux':
    data = ux 
elif data_entry == 'uy':
    data = uy 
elif data_entry == 'uz':
    data = uz 
elif data_entry == 'stream':
    data = stream 

# get grid size
print('Data size =', data.shape)
print('x size =', x.shape)
print('y size =', y.shape)


# create mesh
if plane_entry == 'xy':
   #dim1, dim2 = np.meshgrid(x, y)
   dim1, dim2 = np.meshgrid(y, x)
   #zp = (z.size)/2
   #zp = (z.size)/4
   #zp = 3*(z.size)/4
   #zp = 139
   zp = 1
   data_slice = np.squeeze(data[zp, :, :])
   print('Data slice size =', data_slice.shape)
   print('dim1 size =', dim1.shape)
   print('dim2 size =', dim2.shape)
elif plane_entry == 'xz':
   dim1,dim2 = np.meshgrid(x, z)
   midy = (y.size)/2
   data_slice = np.squeeze(data[:, :, int(midy)])
   print('Data slice size =', data_slice.shape)
   print('dim1 size =', dim1.shape)
   print('dim2 size =', dim2.shape)
elif plane_entry == 'yz':
   dim1,dim2 = np.meshgrid(y, z)
   midx = (x.size)/2
   data_slice = np.squeeze(data[:, int(midx), :])

#fig, ax = plt.subplots(nrows=1)
#fig, (ax0, ax1) = plt.subplots(ncols=2)

# plot data in the midplane and within thermal boundary layer

#bl_point = 10
#data_slice_bl = np.squeeze(data[bl_point:bl_point+1, :, :])


print('Making the figure...')

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
#plt.axis('equal')
plt.axis('off')

plt.show()
#plt.savefig('temperature.png')

