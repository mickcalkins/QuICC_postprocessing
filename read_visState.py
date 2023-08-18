""" 
    This file is designed to read QuICC hdf5 files into python. 
    Before proceeding, use the 'list_visState_structure.py' file 
    to determine the layout/structure of the file since the 
    the structure will depend upon the model used.

"""

import numpy as np
import h5py


def read_params_RBC(filename):

    f = h5py.File(filename, 'r')

    # control parameters
    Ra = float( np.array( f.get('physical/rayleigh') ) )
    Pr = float( np.array( f.get('physical/prandtl') ) )

    return Ra, Pr

def read_params_HMC(filename):

    f = h5py.File(filename, 'r')

    # control parameters
    Q = float( np.array( f.get('physical/chandrasekhar') ) )
    Ra = float( np.array( f.get('physical/rayleigh') ) )

    return Q, Ra

def read_params_MC(filename):

    f = h5py.File(filename, 'r')

    # control parameters
    Q = float( np.array( f.get('physical/chandrasekhar') ) )
    Ra = float( np.array( f.get('physical/rayleigh') ) )
    Pm = float( np.array( f.get('physical/magnetic_prandtl') ) )
    Pr = float( np.array( f.get('physical/prandtl') ) )

    return Q, Ra, Pm, Pr

def read_grid(filename):
 
    f = h5py.File(filename, 'r')

    # x,y,z grids
    x = np.array( f.get('mesh/grid_x') )
    y = np.array( f.get('mesh/grid_y') )
    z = np.array( f.get('mesh/grid_z') )
     
    return x, y, z

def read_data_RBC(filename):
 
    f = h5py.File(filename, 'r')

    # variables
    temp = np.array( f.get('fluct_temperature/fluct_temperature') )
    mean_temp = np.array( f.get('mean_temperature/mean_temperature') )
    u = np.array( f.get('velocity/velocity_x') )
    v = np.array( f.get('velocity/velocity_y') )
    w = np.array( f.get('velocity/velocity_z') )
    vort_x = np.array( f.get('velocity_curl/velocity_curl_x') )
    vort_y = np.array( f.get('velocity_curl/velocity_curl_y') )
    vort_z = np.array( f.get('velocity_curl/velocity_curl_z') )
    
    return temp, mean_temp, u, v, w, vort_x, vort_y, vort_z

def read_one_quantity(filename):
 
    f = h5py.File(filename, 'r')

    # variables
    #temp = np.array( f.get('fluct_temperature/fluct_temperature') )
    #mean_temp = np.array( f.get('mean_temperature/mean_temperature') )
    #u = np.array( f.get('velocity/velocity_x') )
    #v = np.array( f.get('velocity/velocity_y') )
    data = np.array( f.get('velocity/velocity_z') )
    #vort_x = np.array( f.get('velocity_curl/velocity_curl_x') )
    #vort_y = np.array( f.get('velocity_curl/velocity_curl_y') )
    #vort_z = np.array( f.get('velocity_curl/velocity_curl_z') )
    #bx = np.array( f.get('magnetic/magnetic_x') )
    #by = np.array( f.get('magnetic/magnetic_y') )
    #bz = np.array( f.get('magnetic/magnetic_z') )
    
    return data


def read_data_MC(filename):
 
    f = h5py.File(filename, 'r')

    # variables
    temp = np.array( f.get('fluct_temperature/fluct_temperature') )
    mean_temp = np.array( f.get('mean_temperature/mean_temperature') )
    u = np.array( f.get('velocity/velocity_x') )
    v = np.array( f.get('velocity/velocity_y') )
    w = np.array( f.get('velocity/velocity_z') )
    vort_x = np.array( f.get('velocity_curl/velocity_curl_x') )
    vort_y = np.array( f.get('velocity_curl/velocity_curl_y') )
    vort_z = np.array( f.get('velocity_curl/velocity_curl_z') )
    bx = np.array( f.get('magnetic/magnetic_x') )
    by = np.array( f.get('magnetic/magnetic_y') )
    bz = np.array( f.get('magnetic/magnetic_z') )
    
    return temp, mean_temp, u, v, w, vort_x, vort_y, vort_z, bx, by, bz

def read_data_QG(filename):
 
    f = h5py.File(filename, 'r')

    # variables
    temp = np.array( f.get('temperature/temperature') )
    u = np.array( f.get('velocityx/velocityx') )
    v = np.array( f.get('velocityy/velocityy') )
    w = np.array( f.get('velocityz/velocityz') )
    stream = np.array( f.get('streamfunction/streamfunction') )
    vort_z = np.array( f.get('vorticityz/vorticityz') )
    
    return temp, u, v, w, stream, vort_z

def read_data_QGMHD(filename):
 
    f = h5py.File(filename, 'r')

    # variables
    bx = np.array( f.get('fbx/fbx') )
    by = np.array( f.get('fby/fby') )
    bz = np.array( f.get('fbz/fbz') )
    jz = np.array( f.get('fjz/fjz') )
    temp = np.array( f.get('temperature/temperature') )
    u = np.array( f.get('velocityx/velocityx') )
    v = np.array( f.get('velocityy/velocityy') )
    w = np.array( f.get('velocityz/velocityz') )
    stream = np.array( f.get('streamfunction/streamfunction') )
    vort_z = np.array( f.get('vorticityz/vorticityz') )
    
    return bx, by, bz, jz, temp, u, v, w, stream, vort_z

def read_data_Beta(filename):
 
    f = h5py.File(filename, 'r')

    # input parameters
    Ra = np.array( f.get('physical/rayleigh') )

    # grid
    x = np.array( f.get('mesh/grid_x') )
    y = np.array( f.get('mesh/grid_y') )
    z = np.array( f.get('mesh/grid_z') )

    # variables
    temp = np.array( f.get('temperature/temperature') )
    u = -1.0*np.array( f.get('streamfunction_grad/streamfunction_grad_y') )
    v = np.array( f.get('streamfunction_grad/streamfunction_grad_x') )
    w = np.array( f.get('velocityz/velocityz') )
    stream = np.array( f.get('streamfunction/streamfunction') )
    
    return Ra, x, y, z, temp, u, v, w, stream
