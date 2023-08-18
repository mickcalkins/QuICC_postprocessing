'''

This file generates images for depth-averaged flow from many visState files.

'''

from __future__ import print_function

import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from scipy import integrate
#from joblib import Parallel, delayed
#import multiprocessing
from read_visState import read_grid
from read_visState import read_data_QG


# build the list of files to read 
filelist = glob.glob('visState*')

first_file = 'visState0000'

# generate the grid
x, y, z_c = read_grid(first_file)
z = 0.5*(1.0 - z_c)

# create mesh for plot
dim1, dim2 = np.meshgrid(x, y)

    
# read first file
temp, ux, uy, uz, stream, vortz = read_data_QG(first_file)

# define quantity to be plotted
data = stream

# create image properties based on data from first file
fig, ax = plt.subplots(ncols=1)
num_bin = 200

levels = MaxNLocator(nbins=num_bin).tick_values(data.min(), data.max())
cmap = plt.get_cmap('PiYG')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

# loop over the file list
for infile in filelist:
    base = infile[:-5]  # strip off the extension .hdf5
    file_num = base[-4:]
    
    print(file_num)

    temp, ux, uy, uz, stream, vortz = read_data_QG(base)
    
    data = stream
    data = integrate.simps(data[:, :, :], z, axis=0)

    print('Making a figure...')

    baro = ax.pcolormesh(x, y, data, cmap=cmap, norm=norm, shading='gouraud')
    #fig.colorbar(baro, ax=ax, fraction=0.046, pad=0.04)
    ax.set(aspect='equal')
 
    plt.tight_layout()
    ax.axis('off')
    plt.savefig('Barotropic' + file_num + '.png')

