""" Reads in a QuICC kinetic_energy.dat file and plots"""

import numpy as np
import matplotlib.pyplot as plt
from read_parameter_file import read_box
from read_parameter_file import read_nondims
import sys

t_avg = float(sys.argv[1])

print('Data will be averaged after time =', t_avg)


# file name to read
filename = 'parameters.cfg'

# get normalization parameters from file
box2D, box3D, kc2D, kc3D = read_box(filename)

# define normalization constant
fac = (box2D/kc2D)**(2.)

# initialize
time = []
KE = []
Nu = []

# open file
f1 = open('kinetic_energy_all.dat', 'r')
f2 = open('nusselt_all.dat', 'r')
#f1 = open('kinetic_energy.dat', 'r')
#f2 = open('nusselt.dat', 'r')

# read and ignore header lines
header1 = f1.readline()
header2 = f1.readline()
header3 = f1.readline()

header1 = f2.readline()
header2 = f2.readline()
header3 = f2.readline()

# load line by line
for line in f1:
    line = line.strip()
    columns = line.split()
    time.append(columns[0])
    KE.append(columns[1])

f1.close()

for line in f2:
    line = line.strip()
    columns = line.split()
    time.append(columns[0])
    Nu.append(columns[1])

f2.close()


# convert to array
time = np.asarray(time, dtype=np.double)
KE = fac*np.asarray(KE, dtype=np.double)
Nu = np.asarray(Nu, dtype=np.double)
Re = np.sqrt(2.*KE)
#Rm = np.sqrt(2.*KE)
#Re = Rm/Pm 

t_ind = np.where(time > t_avg)[0][0]

print('first index of averaging =', t_ind)


# compute basic stats
#Rm_ave = np.mean(Rm)
#Rm_std = np.std(Rm)
Re_ave = round(np.mean(Re[t_ind:]), 4)
Re_std = round(np.std(Re[t_ind:]), 4)
Reynolds = 'Re = ' + str(Re_ave) + ' +/- ' + str(Re_std)
print(Reynolds)

Nu_ave = round(np.mean(Nu[t_ind:]), 4)
Nu_std = round(np.std(Nu[t_ind:]), 4)
Nusselt = 'Nu = ' + str(Nu_ave) + ' +/- ' + str(Nu_std)
print(Nusselt)




