""" 

    This file plots data from a QuICC visualization file for the tilted magnetoconvection model.

    The data is remapped to a skewed coordinate system aligned with the imposed magnetic field.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import read_visState
#from read_visState import read_params
#from read_visState import read_grid
#from read_visState import read_data

# angle of tilt - hard coded for now
angle_deg = 45.0
angle = angle_deg*(np.pi/180.) # radians
eta3 = np.sin(angle) # z tilt factor

# provide the name of the file to be read - no need for .hdf5
filename = 'visState0000'

# import relevant parameters/quantities/arrays
Q, Ra, Pm, Pr = read_visState.read_params_MC(filename)
#Ra, Pr, Pm = read_nondims(filename)
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

# define tilted magnetic field coordinate
eta = (1./eta3)*z

# create mesh for x-eta plot
dim1,dim2 = np.meshgrid(x, eta)
midy = (y.size)/2
y_loc = int(midy)
data_slice = np.squeeze(data[:, :, y_loc])


# define wavenumber vector in x-direction
k_pos = np.linspace(0, int(Nx/2 - 1), int(Nx/2))
k_neg = np.linspace(int(-Nx/2), -1, int(Nx/2))
kx = np.concatenate((k_pos, k_neg), axis=0)

# convert an xz plane of data to skewed coordinates
data_skew = np.zeros( (Nz,Nx) ) 
data_hat = np.zeros( (Nz,Nx), dtype=np.complex_ ) 
data_hat_p = np.zeros( (Nz,Nx), dtype=np.complex_ ) 

for i in range(0, Nz):
    data_hat[i, :] = np.fft.fft( np.squeeze( data_slice[i, :] ) )

for i in range(0, Nz):
    for j in range(0, Nx):
        data_hat_p[i, j] = data_hat[i, j]*np.exp( -1.j*kx[j]*eta[i] ) 

for i in range(0, Nz):
    for j in range(0, Nx):
        data_skew[i, :] = np.real( np.fft.ifft( data_hat_p[i, :] ) )


#for i in range(0, Nz):
#    for j in range(0, Ny):
#        data_hat[:, i] = np.fft.fft(data_slice[:, i])



print('Making the figure...')

plt.figure(1)

cmap = plt.get_cmap('PiYG')
#norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

plt.pcolormesh(dim1, dim2, data_skew, cmap=cmap, shading='gouraud')
#plt.pcolormesh(X, Y, data_slice_bl, shading='gouraud')
plt.axis('equal')
plt.axis('off')

#plt.savefig('uz_skew.png')
plt.show()


