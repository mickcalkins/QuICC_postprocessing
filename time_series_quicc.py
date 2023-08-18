""" 

This script makes several plots after a quicc simulation.

"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from read_parameter_file import read_box
from read_parameter_file import read_nondims
import sys

# file name to read
filename = 'parameters.cfg'

# get run parameters
Ra, Pr = read_nondims(filename)

# get normalization parameters from file
box2D, box3D, kc2D, kc3D = read_box(filename)

# define normalization constant
fac = (box2D/kc2D)**(2.)

# initialize
time_nu = []
Nu = []
time_ke = []
ke = []

# open file
f_nu = open('nusselt.dat', 'r')
f_ke = open('kinetic_energy.dat', 'r')

# read and ignore header lines
header1 = f_nu.readline()
header2 = f_nu.readline()
header3 = f_nu.readline()
header4 = f_ke.readline()
header5 = f_ke.readline()
header6 = f_ke.readline()

# load line by line
for line in f_nu:
    line = line.strip()
    columns = line.split()
    time_nu.append(columns[0])
    Nu.append(columns[1])

f_nu.close()

# load line by line
for line in f_ke:
    line = line.strip()
    columns = line.split()
    time_ke.append(columns[0])
    ke.append(columns[1])

f_ke.close()

# convert to arrays
time_nu = np.asarray(time_nu, dtype=np.float)
Nu = np.asarray(Nu, dtype=np.float)
time_ke = np.asarray(time_ke, dtype=np.float)
ke = fac*np.asarray(ke, dtype=np.float)
Re = np.sqrt(2.*ke)

Nu_ave = round(np.mean(Nu), 2)
Nu_std = round(np.std(Nu), 4)
Nusselt = 'Nu = ' + str(Nu_ave) + ' +/- ' + str(Nu_std)

Re_ave = round(np.mean(Re), 2)
Re_std = round(np.std(Re), 4)
Reynolds = 'Re = ' + str(Re_ave) + ' +/- ' + str(Re_std)


fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(12,6))

ax0.plot(time_nu, Nu)
ax0.set_xlabel('time')
ax0.set_ylabel('Nusselt')
ax0.annotate(Nusselt, xy=(0.5, 0.95), xycoords='axes fraction')

ax1.plot(time_ke, Re)
ax1.set_xlabel('time')
ax1.set_ylabel('Reynolds')
ax1.annotate(Reynolds, xy=(0.5, 0.95), xycoords='axes fraction')

#plt.plot(time_nu,Nu)
#plt.xlabel('time')
#plt.ylabel('Nusselt')


plt.savefig('time_series_quicc.png')
#plt.show()

