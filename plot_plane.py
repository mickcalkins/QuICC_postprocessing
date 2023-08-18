""" 

    This file plots data from a QuICC visualization file for the plane-layer Navier-Stokes models.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import read_visState
from read_parameter_file import read_box
#from read_visState import read_params
#from read_visState import read_grid
#from read_visState import read_data


# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000.hdf5'

# import relevant parameters/quantities/arrays
#Q, Ra, Pm, Pr = read_visState.read_params_MC(filename)
Ra, Pr = read_visState.read_params_RBC(filename)
#print 'Rayleigh number = ', '{:.2e}'.format(Ra)
#print 'Prandtl number = ', '{:.2e}'.format(Pr)
#print 'Magnetic Prandtl number = ', '{:.2e}'.format(Pm)

# get aspect ratio
box2D, box3D, kc2D, kc3D = read_box('parameters.cfg')
Gamma = box2D*(2.0*np.pi/kc2D) 
print('Aspect ratio = ', Gamma)

# choose data 
entry = input("Enter the quantity to be plotted (temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz):  ")

print('You entered the following:  ', entry)


# choose plane
plane_entry = input("Enter the plane to be plotted (xy, xz, yz):  ")

print('You entered the following:  ', plane_entry)

if not plane_entry in ['xy', 'xz', 'yz']:
    print('You did not enter a valid plane, exiting program')
    sys.exit() 

print('Reading the data...')
x_c, y_c, z_c = read_visState.read_grid(filename)
z = 0.5*(1. - z_c)
x = np.linspace(0, Gamma, num=np.size(x_c), endpoint=False) 
y = np.linspace(0, Gamma, num=np.size(y_c), endpoint=False) 

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

print('Making the figure...')

plt.figure(1)

#cmap = plt.get_cmap('PiYG')
#cmap = plt.get_cmap('bwr')
#cmap = plt.get_cmap('seismic')
#cmap = plt.get_cmap('plasma')
cmap = plt.get_cmap('magma')

# create mesh
if plane_entry == 'xy':
   dim1, dim2 = np.meshgrid(x, y)
   zp = int((z.size)/2)
   #zp = int(8)
   data_slice = np.squeeze(data[zp, :, :])
   plt.pcolormesh(dim1, dim2, np.transpose(data_slice), cmap=cmap, shading='gouraud')
   #print(dim1.shape)
   #print(dim2.shape)
   #print(data_slice.shape)
elif plane_entry == 'xz':
   dim1,dim2 = np.meshgrid(x, z)
   midy = int((y.size)/2)
   data_slice = np.squeeze(data[:, :, midy])
elif plane_entry == 'yz':
   dim1,dim2 = np.meshgrid(y, z)
   midx = int((x.size)/2)
   data_slice = np.squeeze(data[:, midx, :])
   plt.pcolormesh(dim1, dim2, data_slice, cmap=cmap, shading='gouraud')




#norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

#plt.pcolormesh(dim1, dim2, np.transpose(data_slice), cmap=cmap, shading='gouraud')
#plt.pcolormesh(dim2, dim1, data_slice, shading='gouraud')
#plt.pcolormesh(X, Y, data_slice_bl, shading='gouraud')
plt.axis('equal')
plt.axis('off')

plt.savefig( entry + '_' + plane_entry + '_plane' + '.png', dpi=300, bbox_inches="tight")
plt.show()


