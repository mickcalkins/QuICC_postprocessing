""" 

    This file computes the cfl-based timestep requirement

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
filename = 'visState0001.hdf5'

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

print('Reading the data...')
#temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz = read_visState.read_data_MC(filename)
temp, mean_temp, ux, uy, uz, vortx, vorty, vortz = read_visState.read_data_RBC(filename)


# get grid size
Nz, Nx, Ny = uz.shape
print('Grid size(Nx,Ny,Nz) =', Nx, Ny, Nz)

# set up the grids
x_c, y_c, z_c = read_visState.read_grid(filename)
z = 0.5*(1. - z_c)
dz = np.diff(z, append=1)
x = np.linspace(0, Gamma, num=np.size(x_c), endpoint=False) 
dx = abs(x[1]-x[0])
y = np.linspace(0, Gamma, num=np.size(y_c), endpoint=False) 
dy = abs(y[1]-y[0])

dt_z = np.zeros((Nx,Ny,Nz))
dt_x = np.zeros((Nx,Ny,Nz))
dt_y = np.zeros((Nx,Ny,Nz))

for i in range(0, Nx):
    for j in range(0, Ny):
        for k in range(0, Nz):
            dt_z[i,j,k] = abs(dz[k]/uz[k,i,j])
            dt_x[i,j,k] = abs(dx/ux[k,i,j])
            dt_y[i,j,k] = abs(dy/uy[k,i,j])

dt_x_min = np.min(dt_x)
dt_y_min = np.min(dt_y)
dt_z_min = np.min(dt_z)

print('CFL-based timestep sizes:')
print('minimum dt_x =', '{:.2e}'.format(dt_x_min))
print('minimum dt_y =', '{:.2e}'.format(dt_y_min))
print('minimum dt_z =', '{:.2e}'.format(dt_z_min))


