""" Reads in a QuICC kinetic_energy.dat file and plots

Usage: python plot_ke_quicc.py RRBC rotation

"""

import numpy as np
import matplotlib.pyplot as plt
import read_parameter_file 
import sys

model = str(sys.argv[1])
timescale = str(sys.argv[2])

# file name to read
filename = 'parameters.cfg'

# get run parameters 
if model == 'RBC':
   Ra, Pr = read_parameter_file.read_nondims(filename)
elif model == 'RRBC':
   Ra, Pr, Ek = read_parameter_file.read_nondims_RRBC(filename)

#print 'Rayleigh number = ', '{:.2e}'.format(Ra)
#print 'Prandtl number = ', '{:.2e}'.format(Pr)
#print 'Magnetic Prandtl number = ', '{:.2e}'.format(Pm)


# get normalization parameters from file
box2D, box3D, kc2D, kc3D = read_parameter_file.read_box(filename)

# define normalization constant
fac = (box2D/kc2D)**(2.)

# initialize
time = []
ke = []

# open file
#f = open('kinetic_energy.dat', 'r')
f = open('kinetic_energy_all.dat', 'r')


# load line by line
for line in f:
    line = line.strip()
    columns = line.split()
    time.append(columns[0])
    ke.append(columns[1])

# convert to array
time = np.asarray(time, dtype=np.double)
ke = fac*np.asarray(ke, dtype=np.double)

# compute Reynolds number
if timescale == 'viscous':
   Re = np.sqrt(2.*ke)
elif timescale == 'rotation':
   Re = (1.0/Ek)*np.sqrt(2.*ke)

#Rm = np.sqrt(2.*ke)
#Re = Rm/Pm 

# compute basic stats
#Rm_ave = np.mean(Rm)
#Rm_std = np.std(Rm)

Re_ave = round(np.mean(Re), 4)
Re_std = round(np.std(Re), 4)
Reynolds = 'Re = ' + str(Re_ave) + ' +/- ' + str(Re_std)

#print 'Time-averaged Magnetic Reynolds # =', Rm_ave, '+/-', Rm_std
print('Time-averaged Reynolds # =', Re_ave, '+/-', Re_std)


f.close()

#plt.rcParams["mathtext.fontset"] = "cm"
#plt.rcParams.update({'font.size': 10})

#plt.plot(time,ke)
plt.plot(time,Re)
#plt.plot(time,Re_ave*np.ones(len(time)), 'k--')
#plt.plot(time,Re_ave*np.ones(len(time))+Re_std, 'k-.')
#plt.plot(time,Re_ave*np.ones(len(time))-Re_std, 'k-.')
plt.xlabel(r'$t$')
plt.ylabel(r'$Re$', rotation=0)
#plt.ylabel(r'$KE$', rotation=0)

plt.annotate(Reynolds, xy=(0.5, 0.95), xycoords='axes fraction')

plt.savefig('Reynolds_all_vs_time.png')
#plt.show()

