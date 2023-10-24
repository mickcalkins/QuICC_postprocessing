import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import visStateMethods2 as vsm
import scipy.fftpack as sft
import read_parameter_file as rpf
import h5py as hf
import numpy.ma as ma
import time

def energy_spectra(directory, fn = 'All_Files', comp=0,interp = True):
#comp-- array to specify which element(s) to take spectra of 0-- full field, 1--barotropic, 2--baroclinic. Can choose more than one value ([0,1,2]) gives all values
#interp -- boolean to specify whether to interpolate between integer k.
    param = directory + 'parameters.cfg'
    box2D,kc2D = rpf.read_box(param)
    dir_contents = []
    print(box2D,kc2D)
    if fn=='All_Files':
        for filename in os.listdir(directory):
            if filename.startswith('visState0') and not(filename.startswith('visState0000')): 
                dir_contents.append(filename)
        dir_contents.sort()
        leading_str = ''
        #full_file = directory + filename
    else:
        dir_contents.append(fn)
        leading_str = 'Energy_Spectra_Files/'+fn[8:12]
    MYDIR = (directory + "Energy_Spectra_Files")
    CHECK_FOLDER = os.path.isdir(MYDIR)
    
    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print("created folder : ", MYDIR)
    
    else:
        print(MYDIR, "folder already exists.")
    
    print(dir_contents)
    LL = len(dir_contents)
    KE = np.empty(LL,dtype=object)
    Ia = np.empty(LL,dtype=object)
    Ta = np.empty(LL,dtype=object)
    is_long = np.zeros(LL)
    ts_long = np.zeros(LL)
    
    index = 0
    print('Beginning to search through hdf5')
    #File to track progress
    progress=open(leading_str+'output_energy_spectra_progress.txt','a') 
    tag_str = ['','barotropic','baroclinic'] 
    bt = comp 
    for full_file in dir_contents:
        print(index)
        with hf.File(full_file,'r') as f:
            #vx =f['/velocityx/velocityx'][:]
            #vy =f['/velocityy/velocityy'][:]
            #vz =f['/velocityz/velocityz'][:]
            vx =f['/velocity/velocity_x'][:]
            vy =f['/velocity/velocity_y'][:]
            vz =f['/velocity/velocity_z'][:]
        if not interp:
            (Ia[index],Ta[index],KE[index],is_long[index],ts_long[index]) = vsm.integralscale(vx, vy, vz,kc2D/box2D, bt = bt,interp = interp)
        else:
            (Ia[index],Ta[index],KE[index]) = vsm.integralscale(vx, vy, vz,kc2D/box2D, bt = bt,interp = interp)
        #Write to progress file
        progress.write('Wrote Energy Spectra for ' + full_file + ' ('+tag_str[bt]+')' + time.strftime("%H:%M:%S")+'\n' )
        progress.flush()
        
        index += 1    
   
    nz,nx,ny = np.shape(vx)
    
    
    z = vsm.gauss_cheb(nz)
    z = -0.5*z + 0.5
    kmax = len(KE[0]) #Number of modes. = maximum k value + 1
    da_KE = np.empty((LL,kmax),dtype=float)
    for index in range(LL):
        for kindex in range(kmax): #loop over each K to depth average
            da_KE[index,kindex] = vsm.depth_average(KE[index][kindex,:])
            #da_KE[index,kindex] = np.trapz(KE[index][kindex,:],z)
      
    sumtoavg = np.zeros(kmax)
    for index in range(LL):
        sumtoavg=sumtoavg + da_KE[index,:]
    avg = sumtoavg/LL        
      
        
    k = np.linspace(0,kmax-1,kmax,dtype = int)
    
    df = pd.DataFrame(list(zip(k, avg)), columns=['k', 'E_k'])
    df.to_csv(leading_str+'energy_spectra'+tag_str[bt]+'.csv', sep='\t', index=False)
    if LL > 0:
        ta_Ia = np.average(Ia)
        ta_Ta = np.average(Ta)

    data_Ia = 0.5*ta_Ia + 0.5*np.flip(ta_Ia,axis=0)
    data_Ta = 0.5*ta_Ta + 0.5*np.flip(ta_Ta,axis=0)
    mIs_v = vsm.depth_average(data_Ia)
    mTs_v = vsm.depth_average(data_Ta)
    #mIs_v = np.trapz(data_Ia,z)
    #mTs_v=np.trapz(data_Ta,z)
    if not interp:
        cms= ['Integral Scale','Taylor microscale','IS_Long','TS_Long']          
        LV = [[mIs_v,mTs_v,np.average(is_long),np.average(ts_long)]]
    else:
        cms= ['Integral Scale','Taylor microscale']          
        LV = [[mIs_v,mTs_v]]
    df2 = pd.DataFrame(LV,columns = cms)
    df2.to_csv(leading_str+'Length_Scale'+tag_str[bt]+'.txt', sep = ',', index=False)   
    
