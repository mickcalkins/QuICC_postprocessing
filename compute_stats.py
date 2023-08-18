""" 
    This file computes basic z-dependent stats from QuICC plane layer simulation data. 

"""

import numpy as np
import h5py


def stats(data):
   
    Nz, Nx, Ny = data.shape

    # average of data
    data_ave_x = np.average(data, axis=1)
    data_ave = np.average(data_ave_x, axis=1)
    data_ave_mat = np.transpose( np.tile(data_ave, (Ny, Nx, 1)) )
    
    # rms of data
    data_no_ave = data - data_ave_mat
    data2 = np.square(data_no_ave)
    mx = np.average(data2, axis=1)
    mxy = np.average(mx, axis=1)
    data_rms = np.sqrt(mxy)

    return data_ave, data_rms


