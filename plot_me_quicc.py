""" Reads in a QuICC magnetic_energy.dat file and plots"""

import numpy as np
import matplotlib.pyplot as plt
from read_parameter_file import read_box
from read_parameter_file import read_nondims

# file name to read
filename = 'parameters.cfg'

# get run parameters
Ra, Pr, Pm = read_nondims(filename)
print 'Rayleigh number = ', '{:.2e}'.format(Ra)
print 'Prandtl number = ', '{:.2e}'.format(Pr)
print 'Magnetic Prandtl number = ', '{:.2e}'.format(Pm)

# get normalization parameters from file
box2D, kc2D = read_box(filename)

# define normalization constant
fac = (box2D/kc2D)**(2.)

# initialize
time = []
ke = []
me = []

# open files
fk = open('kinetic_energy.dat', 'r')
fm = open('magnetic_energy.dat', 'r')

# read and ignore header lines
header1 = fm.readline()
header2 = fm.readline()
header3 = fm.readline()
header4 = fk.readline()
header5 = fk.readline()
header6 = fk.readline()

# load line by line
for line in fm:
    line = line.strip()
    columns = line.split()
    time.append(columns[0])
    me.append(columns[1])

# load line by line
for line in fk:
    line = line.strip()
    columns = line.split()
    ke.append(columns[1])

# convert to array
time = np.asarray(time, dtype=np.float)
#me = (1./np.sqrt(Pm))*fac*np.asarray(me, dtype=np.float)
me = Pm*fac*np.asarray(me, dtype=np.float)
ke = fac*np.asarray(ke, dtype=np.float)
ratio = me/ke
B = np.sqrt(2.*me)

# compute basic stats
B_ave = np.mean(B)
B_std = np.std(B)

print 'Time-averaged B-field =', B_ave, '+/-', B_std

fm.close()

plt.plot(time,me)
#plt.plot(time,ratio)
#plt.plot(time,B)
#plt.plot(time,B_ave*np.ones(len(time)), 'k--')
#plt.plot(time,B_ave*np.ones(len(time))+B_std, 'k-.')
#plt.plot(time,B_ave*np.ones(len(time))-B_std, 'k-.')
plt.xlabel(r'$t_{\eta}$')
#plt.ylabel(r'$B_{rms}$', rotation=0)
#plt.ylabel(r'$E_{mag}/E_{kin}$', rotation=90)
plt.ylabel(r'$E_{mag}$', rotation=90)

plt.show()

