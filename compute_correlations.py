""" 
    
    This file computes correlations from a QuICC visualization file.

"""

from __future__ import print_function
import pandas as pd
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import read_visState
from read_parameter_file import read_box

# 2D/3D?
dimensions = str(sys.argv[1])

# search for visState files
path = "."
file_names = np.flipud([fn for fn in os.listdir(path)
        if fn.startswith('visState')])

num_files = len(file_names)

print('Found', num_files,'files:', file_names[:])

first_file = file_names[0]

# list parameters and find size of array
#Ra, Pr = read_visState.read_params_RBC(first_file)
Q, Ra = read_visState.read_params_HMC(first_file)
print('Rayleigh number = ', '{:.2e}'.format(Ra))
print('Chandrasekhar number = ', '{:.2e}'.format(Q))

# get aspect ratio
box2D, box3D, kc2D, kc3D = read_box('parameters.cfg')
Gamma = box2D*(2.0*np.pi/kc2D) 
#print('Aspect ratio = ', Gamma)

uz = read_visState.read_one_quantity(first_file)
data = uz
data_shape = data.shape
numz = data_shape[0]
numx = data_shape[1]
numy = data_shape[2]

# get grid size
#print( 'Grid size (z, x, y) =', data.shape )

x_c, y_c, z_c = read_visState.read_grid(first_file)
z = 0.5*(1. - z_c)
x = np.linspace(0, Gamma, num=np.size(x_c), endpoint=False) 
y = np.linspace(0, Gamma, num=np.size(y_c), endpoint=False) 
dx = x[1]-x[0]
dy = y[1]-y[0]

ndata = np.zeros(data_shape)
corr_x_yz = np.zeros(data_shape)
acorr_x_yz = np.zeros(data_shape)
corr_y_xz = np.zeros(data_shape)
acorr_y_xz = np.zeros(data_shape)
acorr_x_all = np.zeros([int(numx),int(num_files)])
acorr_x = np.zeros(int(numx))
acorr_y_all = np.zeros([int(numy),int(num_files)])
acorr_y = np.zeros(int(numy))

kk = 0

# loop through all files
for k in file_names:
    print('Reading file: ', k)
    data = read_visState.read_one_quantity(k)

    if dimensions == '3D':
       # compute autocorrelation in x
       # loop over each z- and y-value
       for i in range(numz):
           for j in range(numy):
               mean = np.mean(data[i,:,j])
               var = np.var(data[i,:,j])
               ndata[i,:,j] = data[i,:,j]-mean
               corr_x_yz[i,:,j] = signal.correlate(ndata[i,:,j], ndata[i,:,j], mode='full', method='fft')[len(ndata[i,:,j])-1:]
               acorr_x_yz[i,:,j] = corr_x_yz[i,:,j]/var/len(ndata[i,:,j])

       # average in y and z
       acorr_x_z = np.mean(acorr_x_yz,axis=2)
       acorr_x_avg = np.trapz(acorr_x_z, z, axis=0) 
       acorr_x = acorr_x_avg + acorr_x
       acorr_x_all[:,kk] = acorr_x_avg 

    # compute autocorrelation in y
    # loop over each z- and x-value
    for i in range(numz):
        for j in range(numx):
            mean = np.mean(data[i,j,:])
            var = np.var(data[i,j,:])
            ndata[i,j,:] = data[i,j,:]-mean
            corr_y_xz[i,j,:] = signal.correlate(ndata[i,j,:], ndata[i,j,:], mode='full', method='fft')[len(ndata[i,j,:])-1:]
            acorr_y_xz[i,j,:] = corr_y_xz[i,j,:]/var/len(ndata[i,j,:])

    # average in x and z
    acorr_y_z = np.mean(acorr_y_xz,axis=1)
    acorr_y_avg = np.trapz(acorr_y_z, z, axis=0) 
    acorr_y = acorr_y_avg + acorr_y
    acorr_y_all[:,kk] = acorr_y_avg 
    
    kk = kk + 1

# take average value
acorr_x = acorr_x/num_files
acorr_y = acorr_y/num_files

integral_length_y = np.trapz(y, acorr_y)
zero_length_index_y = np.where(acorr_y < 0)[0][0]
zero_length_y = y[zero_length_index_y-1] + 0.5*dy
print('y-correlation length =', integral_length_y )
print('y-zero crossing length =', zero_length_y )

# write the data to file
daty = {'y': y, 'auto_corr_y': acorr_y}
df_y = pd.DataFrame(data=daty)
df_y.to_excel('auto_corr_y' + '_Q' + '{:.1e}'.format(Q) + '_Ra' + '{:.1e}'.format(Ra) + '.xlsx')


if dimensions == '3D':
   integral_length_x = np.trapz(x, acorr_x)
   zero_length_index_x = np.where(acorr_x < 0)[0][0]
   zero_length_x = x[zero_length_index_x-1] + 0.5*dx
   print('x-correlation length =', integral_length_x )
   print('x-zero crossing length =', zero_length_x )
   datx = {'x': x, 'auto_corr_x': acorr_x}
   df_x = pd.DataFrame(data=datx)
   df_x.to_excel('auto_corr_x' + '_Q' + '{:.1e}'.format(Q) + '_Ra' + '{:.1e}'.format(Ra) + '.xlsx')


# plots here
plt.figure(1)
plt.plot(y, acorr_y_all, 'k--', y, acorr_y, 'r-')
plt.xlabel('y')
plt.ylabel(r'R_y')
plt.savefig('auto_correlation_y.png')

if dimensions == '3D':
   plt.figure(2)
   plt.plot(x, acorr_x_all, 'k--', x, acorr_x, 'r-')
   plt.xlabel('x')
   plt.ylabel(r'R_x')
   plt.savefig('auto_correlation_x.png')




