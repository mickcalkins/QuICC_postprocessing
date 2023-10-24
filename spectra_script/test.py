import spectra as pp

#calculates spectra/ls for visState0001.hdf5
#Saves results in Energy_Spectra_Files/0001{energy_spectra.csv,Length_Scale.txt}
pp.energy_spectra('./',fn ='visState0001.hdf5')
#calculates spectra/ls for all visStates
#Saves time averaged results in {energy_spectra.csv,Length_Scale.txt}
pp.energy_spectra('./')

