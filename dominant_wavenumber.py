""" 

    This file finds the dominant wavenumber from a horizontal FFT.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import read_visState
from read_parameter_file import read_box
#from read_visState import read_params
#from read_visState import read_grid
#from read_visState import read_data

pi = np.pi

# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000'

# import relevant parameters/quantities/arrays
#Q, Ra, Pm, Pr = read_visState.read_params_MC(filename)
Ra, Pr = read_visState.read_params_RBC(filename)
#print 'Rayleigh number = ', '{:.2e}'.format(Ra)
#print 'Prandtl number = ', '{:.2e}'.format(Pr)
#print 'Magnetic Prandtl number = ', '{:.2e}'.format(Pm)

# get aspect ratio
box2D, box3D, kc2D, kc3D = read_box('parameters.cfg')
Gamma = box2D*(2.0*pi/kc2D) 

x_c, y_c, z_c = read_visState.read_grid(filename)
z = 0.5*(1. - z_c)
x = np.linspace(0, Gamma, num=np.size(x_c), endpoint=False) 
y = np.linspace(0, Gamma, num=np.size(y_c), endpoint=False) 
z_mid = int((z.size)/2)

temp, mean_temp, ux, uy, uz, vortx, vorty, vortz, bx, by, bz = read_visState.read_data_MC(filename)

data = temp

# get grid size
Nz, Nx, Ny = data.shape
print('Grid size(Nx,Ny,Nz) =', Nx, Ny, Nz)

kx_max = Nx/2
ky_max = Ny/2

kx = np.arange(0, kx_max+1, 1)
ky = np.arange(0, ky_max+1, 1)
x = np.arange(0, 2.0*pi, (2.0*pi/Nx))
y = np.arange(0, 2.0*pi, (2.0*pi/Ny))
X, Y = np.meshgrid(x, y)
Kx, Ky = np.meshgrid(kx, ky)

data_slice = np.transpose( np.squeeze(data[z_mid, :, :]) )

f_hat = np.fft.rfft2(data_slice)
Power = np.real( f_hat*np.conjugate(f_hat) )
Power_real = Power[0:int(ky_max)+1, 0:int(kx_max)+1]
result = np.where(Power_real == np.amax(Power_real))
listOfCoordinates = list(zip(result[0], result[1]))
print('dominant wavenumber [k_y, k_x] = ', listOfCoordinates )

plt.pcolormesh(Kx[0:10, 0:10], Ky[0:10, 0:10], Power_real[0:10, 0:10], shading='auto')
plt.xlabel(r'k_x')
plt.ylabel(r'k_y')


#norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

#plt.pcolormesh(dim1, dim2, np.transpose(data_slice), cmap=cmap, shading='gouraud')
#plt.pcolormesh(dim2, dim1, data_slice, shading='gouraud')
#plt.pcolormesh(X, Y, data_slice_bl, shading='gouraud')
plt.axis('equal')
#plt.axis('off')

plt.savefig( 'wavenumber' + '.png')
#plt.show()


