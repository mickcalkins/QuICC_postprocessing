""" 
    This file is designed to read QuICC 'state' hdf5 files into python. 
    Before proceeding, use the 'list_state_structure.py' file 
    to determine the layout/structure of the file since the 
    the structure will depend upon the model used.

    This version is for the following model:

        BoussinesqPlaneMCModel

"""

import numpy as np
import h5py

global file_ext
file_ext = '.hdf5'

def read_params(filename):

    f = h5py.File(filename + file_ext, 'r')

    # control parameters
    Q = np.float( np.array( f.get('physical/chandrasekhar') ) )
    Ra = np.float( np.array( f.get('physical/rayleigh') ) )
    Pm = np.float( np.array( f.get('physical/magnetic_prandtl') ) )
    Pr = np.float( np.array( f.get('physical/prandtl') ) )

    return Q, Ra, Pm, Pr

def read_truncation(filename):
 
    f = h5py.File(filename + file_ext, 'r')

    # truncation 
    trunc_phys = np.array( f.get('truncation/physical/dim1D') )
    trunc_spec = np.array( f.get('truncation/spectral/dim1D') )
    trunc_trans = np.array( f.get('truncation/transform/dim1D') )
     
    return trunc_phys, trunc_spec, trunc_trans 

def read_data(filename):
 
    f = h5py.File(filename + file_ext, 'r')

    # variables
    vel_pol = np.array( f.get('velocity/velocity_pol') )
    vel_tor = np.array( f.get('velocity/velocity_tor') )
        
    return vel_pol, vel_tor 
