visStateMethods2.py
    Contains subroutines for calculating spectra and length scales.
   
     integralscale
        The main function is integralscale which is equivalent to Ming's matlab script of the same name. 
        One extra feature is that you can specify whether you want full field, depth-averaged, or depth-averaged removed using the 'bt' flag. 
        All of the crap with masked arrays and sparse arrays below line 100 is intended to reduce the memory requirements (similarly with the constant deleting of variables).
   
     gauss_cheb
        I believe that this is almost exactly the same as Ming's script, except that maybe it varies by a pi/2N. regardless, it returns the Chebyshev grid that QuICC uses for the F3DQG model. 

spectra.py

    energy_spectra
        Uses the ouput of integral scale to write arrays of spectra.
        #directory -- string to specify location of visStates
        #fn -- string to specify visState. default 'All_Files'. When set to 'All_Files' the routine will get spectra for each visState in the directory and time average the result. 'energy_spectra.csv' and 'Length_Scale.txt' will be written into the directory and contain only the time averaged results. Instead, you can opt to specify the full name of a visState (ie 'visState0001.hdf5') and the routine will calculate spectra and length scales for that file only. It will save the results in the directory Energy_Spectra_Files as (0001energy_spectra.csv/ 0001Length_Scale.txt).  
        #comp-- array to specify which element(s) to take spectra of 0-- full field, 1--barotropic, 2--baroclinic. default is 0 (full field) (This is just bt in visStateMethods2.integralscale)
        #interp -- default = True. boolean to specify whether to return length scales that are calulated with non-interpolated kinetic energy spectra. In my experience this doesn't really change much. 

            if False, then returns length scales calculated with interpolated and non-interpolated kinetic energy spectra. 
            if True, returns only length scales calculated with interpolated ke spectra.
            
test.py
    Just for testing. plop a parameters.cfg and visState file into this directory and check to see if it works. NOTE: the scripts ignore files named "visState0000.hdf5" Since the way QuICC makes visStates the 0000 one is almost always a duplicate (since you copy it).
