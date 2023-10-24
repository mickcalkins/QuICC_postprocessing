import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import h5py as hf
import numpy.fft as ft
import read_parameter_file as rpf
import numpy.matlib as mb
import numpy.ma as ma
import scipy.fftpack as sft
from scipy import sparse

def integralscale(A,B,C,kcn=1.3048/10,bt =0,interp = True):
#Takes three physical space fields and converts them
#into the Fourier baroclinic spectrum
#A -- field in physical space (eg. vx)
#B -- field in physical space (eg. vy)
#C -- field in physcial space (eg. vz)
#kcn -- non-dimensional critical wavelength/#of critical wavenumbers
#bt -- tag to determine which component to evaluate. 0--full field, 1-- barotropic, 2-- baroclinic
#Returns
#Ia -- Integralscale (nz x 1 vector)
#Ta -- Taylor microscale (nz x 1 vector)
#CS -- Total energy spectra for each z layer(kx (Dim2D) x nz matrix) 

    #FFT A,B,C
    if bt==0:
        pass
    elif bt==1:
        A = cbp(A)
        B = cbp(B)
        C = cbp(C)
    elif bt==2:
        A = A-cbp(A)
        B = B-cbp(B)
        C = C - cbp(C)

    C_hat = ft.fft2(C)
    C2 = np.abs(C_hat)**2
    del C_hat
    
    A_hat = ft.fft(ft.fft(A,axis=1),axis=2)#,norm = 'forward')
    A2 = np.abs(A_hat)**2
   
    del A_hat
    
    B_hat = ft.fft2(B)#,norm = 'forward')
    B2 = np.abs(B_hat)**2
    
    del B_hat
    nz,nx,ny = np.shape(A)
    N = nx        #Number of horizontal points in physical space
    order = np.append(range(0,int(N/2+1)),range(-int(np.ceil(N/2))+1,0)) #for initializing grid
    KX,KY = np.meshgrid(order,order,indexing='ij') #K grid
    if not interp:
        is_long = []
        ts_long = []
        key_x = int(nx/2+1)
        kx_int = order#np.array(range(0,int(nx/2+1)))
        key_y = int(ny/2+1)
        ky_int =np.append(range(0,int(ny/2+1)),range(-int(np.ceil(ny/2))+1,0))# np.array(range(0,int(ny/2+1)))
       #Really, I should not multiply the 000 mode by 4
        E = A2 + B2 + C2
        K = np.sqrt(KX**2 + KY**2)
      #  E = 1*(A2[:,0:key_x,0:key_y] + B2[:,0:key_x,1:key_y] + C2[:,0:key_x,0:key_y])
      #  K = np.sqrt(KX[0:key_x,0:key_y]**2 + KY[0:key_x,0:key_y]**2)
        K[0,0] = 0.0001
        K= np.fft.fftshift(K,axes = (0,1))
        kx_int = np.fft.fftshift(kx_int)
        ky_int = np.fft.fftshift(ky_int)
        for E_perp in E:
            E_perp[0,0] = 0
            E_perp = np.fft.fftshift(E_perp,axes = (0,1))
            i_num = E_perp*(1/K)
            t_den = E_perp*(K*K)
            a = np.trapz([np.trapz(n,ky_int) for n in i_num],kx_int)
            b = np.trapz([np.trapz(n,ky_int) for n in E_perp],kx_int)
            is_long = np.append(is_long,a/b/kcn)
            c = np.trapz([np.trapz(n,ky_int) for n in t_den],kx_int)
            ts_long = np.append(ts_long,np.sqrt(b/c)/kcn)   
        print(is_long)
        del K
        del E,E_perp,i_num,t_den 

        isl = depth_average(is_long)
        tsl = depth_average(ts_long)   
    CS = np.zeros((N//2+1,nz)) #Initialize array for total kinetic energy
    k = np.sqrt(KX**2 + KY**2)
    k_trunc = np.floor(k)
    r = k -np.floor(k)
    PowerSpectrum = A2 + B2 + C2
    del A2
    del B2 
    del C2
    PS_interpolate = np.multiply(1.0-r,PowerSpectrum)
    PS_interpolate2 = r*PowerSpectrum
    
    P1 = PS_interpolate.reshape(*PS_interpolate.shape[:1],-1)         	
    P2 = PS_interpolate2.reshape(*PS_interpolate2.shape[:1],-1)
    mask1D = np.where(k_trunc==0,1,0)
    mask1D = mask1D.flatten()
    mask1D = mask1D[np.newaxis,:]
    mask = np.repeat(mask1D,nz,axis=0)
    mask = sparse.coo_matrix(mask)
    
    for kappa in range(N//3+1):
        A = mask.multiply(P1)
        #PS1_masked = PS1_masked.reshape(*PS1_masked.shape[:1],-1)

        sm = np.asarray(A.sum(axis=1)).flatten()

        CS[kappa,:] = CS[kappa,:] +sm# np.sum(PS1_masked,axis = (1))

        A = mask.multiply(P2)#np.multiply(mask,PS_interpolate2)
        
        CS[kappa+1,:] = np.asarray(A.sum(axis=1)).flatten()

        mask1D = np.where(k_trunc==kappa+1,1,0)

        mask1D = mask1D.flatten()

        mask = sparse.coo_matrix(mask1D)

        ONES = np.transpose(np.ones((1,nz)))

        mask = mask.multiply(ONES)
    
    CS = (1.0/2)*np.abs(CS)/(N**4) #Normalize



    kperp = np.linspace(0,int(N/2),int(N/2) + 1)
    kp = kperp[:,np.newaxis]
    kp = np.repeat(kp,nz,axis =1)

    k1 = 1

    #Initialize
    ET = np.empty((nz),dtype = float)
    I1 = np.empty(nz, dtype = float)
    T1 = np.empty(nz, dtype = float)
    ET = np.trapz(CS,kperp,axis = 0)
    I1 = np.trapz(CS[1:,:]/kp[1:,:],kperp[1:],axis = 0)/kcn
    T1 = np.trapz((CS*(kp**2)),kperp[:],axis = 0)*kcn*kcn
    
    Ia = I1/ET
    Ta = np.sqrt(ET/T1)
    if not interp:
        return(Ia,Ta,CS,isl,tsl)
    return(Ia,Ta,CS)

# Returns a Chebyshev grid on the interval [-1,1]
# with entry x[0] = 1 and x(end) = -1
def gauss_cheb(N):
    a = np.linspace(1,N,N)
    x = np.cos(np.pi*(2*a -1)/2/N)
    return x
def cbp(A,switch_on = True):
#Removes depth averaged component of A if switch_on==True
#returns an array of 0's the size of A if switch_on == False
    (nz,nx,ny) = np.shape(A)
    val = np.zeros((nz,nx,ny))
    if switch_on:
        ak = 1./(nz)*sft.dct(A,axis = 0,type = 2)
        ak[0,:,:] = 0.5*ak[0,:,:]
        weights = np.zeros((nz,nx,ny))
        ONES = np.ones((nx,ny))
        for i in range(nz):
            #print(i%2) 
            if i%2 ==0:
                weights[i,:,:] = 1/(1-i*i)*ONES
            weights[0,:,:] = ONES
        prod = weights * ak
        mat = np.sum(prod,axis=0)
        for i in range(nz):
            val[i,:,:] = mat
        return val
    else:
        return val

def depth_average(A,dct_type = 2):
#type 1 vs. type 2 is really a matter of what your collocation points are.
#If your collocation points are Gauss-Lobatto points (x_i = cos(pi i/N), i=0,...,N-1) you should use type 1.
#If your collocation points are Gauss points (x_i = cos(pi (2*i+1)/(2N)), i=0,...,N-1) you should use type 2.
#QuICC uses Gauss points
    nz = len(A)
    nz_hat = int(nz*2/3) - 1
    #val = np.zeros((nz,nx,ny))
    if dct_type == 2:
        ak = 1./(nz)*sft.dct(A,axis = 0,type=2)
    if dct_type ==1:
        ak = 1./(nz-1)*sft.dct(A,axis = 0,type =1)
    ak = ak[0:nz_hat+1]
    #ONES = np.ones((nx,ny))
    ak[0] = 0.5 * ak[0]
   # ak[-1,:,:] = 1 * ak[-1,:,:]
    weights = np.zeros(nz_hat+1)
    for i in range(1,nz_hat+1,1):
        #print(i%2) 
        if i%2 ==0:
            weights[i] = 1/(1-i*i)
        weights[0] =1# ONES
    prod = weights*ak
    #prod = 1./2*weights * ak
    mat = np.sum(prod,axis=0)
    return(mat)
