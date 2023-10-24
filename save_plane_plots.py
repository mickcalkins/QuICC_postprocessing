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

# just used to name file
run_entry = str(sys.argv[1])

# provide the name of the file to be read 
filename = 'visState0036.hdf5'

# import relevant parameters/quantities/arrays
#Q, Ra, Pm, Pr = read_visState.read_params_MC(filename)
#Ra, Pr = read_visState.read_params_RBC(filename)
#print 'Rayleigh number = ', '{:.2e}'.format(Ra)
#print 'Prandtl number = ', '{:.2e}'.format(Pr)
#print 'Magnetic Prandtl number = ', '{:.2e}'.format(Pm)

# get aspect ratio
#box2D, box3D, kc2D, kc3D = read_box('parameters.cfg')
Gamma = 6*(2.0*np.pi/2.22) 

# choose data 
#entry = input("Enter the quantity to be plotted (temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz):  ")

#print('You entered the following:  ', entry)


x_c, y_c, z_c = read_visState.read_grid(filename)
z = 0.5*(1. - z_c)
x = np.linspace(0, Gamma, num=np.size(x_c), endpoint=False) 
y = np.linspace(0, Gamma, num=np.size(y_c), endpoint=False) 

temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz = read_visState.read_data_MC(filename)


data = vortz

#cmap = plt.get_cmap('PiYG')
#cmap = plt.get_cmap('bwr')
#cmap = plt.get_cmap('seismic')
#cmap = plt.get_cmap('plasma')
cmap = plt.get_cmap('magma')

plt.figure(1)
dim1, dim2 = np.meshgrid(x, y)
zp = int((z.size)/2)
data_slice = np.squeeze(data[zp, :, :])
plt.pcolormesh(dim1, dim2, np.transpose(data_slice), cmap=cmap, shading='gouraud')
plt.axis('equal')
plt.axis('off')
plt.savefig( 'temp' + '_midplane_' + run_entry + '.png', dpi=300, bbox_inches="tight")

plt.figure(2)
zp = int(0)
data_slice = np.squeeze(data[zp, :, :])
plt.pcolormesh(dim1, dim2, np.transpose(data_slice), cmap=cmap, shading='gouraud')
plt.axis('equal')
plt.axis('off')
plt.savefig( 'temp' + '_bl_' + run_entry + '.png', dpi=300, bbox_inches="tight")



#if plane_entry == 'xz':
#dim1,dim2 = np.meshgrid(x, z)
#midy = int((y.size)/2)
#data_slice = np.squeeze(data[:, :, midy])

#if plane_entry == 'yz':
plt.figure(3)
dim1,dim2 = np.meshgrid(y, z)
midx = int((x.size)/2)
data_slice = np.squeeze(data[:, midx, :])
plt.pcolormesh(dim1, dim2, data_slice, cmap=cmap, shading='gouraud')
plt.axis('equal')
plt.axis('off')
plt.savefig( 'temp' + '_yz_' + run_entry + '.png', dpi=300, bbox_inches="tight")




