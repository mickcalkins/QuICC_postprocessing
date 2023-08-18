""" 

Reads in a QuICC kinetic_energy.dat file and plots it.

"""

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from read_parameter_file import read_box
from read_parameter_file import read_nondims

# file name to read
filename = 'parameters.cfg'

# get run parameters
Ra, Pr = read_nondims(filename)
print('Rayleigh number = ', '{:.2e}'.format(Ra))
print('Prandtl number = ', '{:.2e}'.format(Pr))
#print 'Magnetic Prandtl number = ', '{:.2e}'.format(Pm)

# get normalization parameters from file
box2D, box3D, kc2D, kc3D = read_box(filename)

# define normalization constant
#fac = (box2D/kc2D)**(2.)
fac = (box2D*box3D)/(kc2D*kc3D)

# initialize
time = []
ke = []
ke1 = []
ke2 = []

# open file
f = open('kinetic_energy.dat', 'r')

# read and ignore header lines
header1 = f.readline()
header2 = f.readline()
header3 = f.readline()

# load line by line
for line in f:
    line = line.strip()
    columns = line.split()
    time.append(columns[0])
    ke.append(columns[1])
    ke1.append(columns[2])
    ke2.append(columns[3])

# convert to array
time = np.asarray(time, dtype=np.float)
ke = fac*np.asarray(ke, dtype=np.float)
ke1 = fac*np.asarray(ke1, dtype=np.float)
ke2 = fac*np.asarray(ke2, dtype=np.float)
#Rm = np.sqrt(2.*ke)
#Re = Rm/Pm 

# compute basic stats
#Rm_ave = np.mean(Rm)
#Rm_std = np.std(Rm)
#Re_ave = np.mean(Re)
#Re_std = np.std(Re)
ke_ave = np.mean(ke)
ke1_ave = np.mean(ke1)
ke2_ave = np.mean(ke2)

print('Average KE = ', ke_ave)
print('Average KE_1 = ', ke1_ave)
print('Average KE_2 = ', ke2_ave)
#print 'Time-averaged Magnetic Reynolds # =', Rm_ave, '+/-', Rm_std
#print 'Time-averaged Reynolds # =', Re_ave, '+/-', Re_std


f.close()

#plt.rcParams["mathtext.fontset"] = "cm"
#plt.rcParams.update({'font.size': 10})

plt.plot(time,ke)
#plt.plot(time,Re)
#plt.plot(time,Re_ave*np.ones(len(time)), 'k--')
#plt.plot(time,Re_ave*np.ones(len(time))+Re_std, 'k-.')
#plt.plot(time,Re_ave*np.ones(len(time))-Re_std, 'k-.')
plt.xlabel(r'$t$')
#plt.ylabel(r'$Re$', rotation=0)
plt.ylabel(r'$KE$', rotation=0)
plt.show()

