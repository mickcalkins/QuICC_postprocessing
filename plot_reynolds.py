""" Reads in a QuICC kinetic_energy.dat file and plots

Usage: python plot_ke_quicc.py RRBC rotation

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
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
f = open('kinetic_energy.dat', 'r')
#f = open('kinetic_energy_all.dat', 'r')

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

# convert to array
time = np.asarray(time, dtype=float)
ke = fac*np.asarray(ke, dtype=float)

# compute Reynolds number
if timescale == 'viscous':
   Re = np.sqrt(2.*ke)
elif timescale == 'rotation':
   Re = (1.0/Ek)*np.sqrt(2.*ke)

# cumulative average
Re_series = pd.Series(Re)
windows = Re_series.expanding()
moving_averages = windows.mean()
Re_moving_averages = moving_averages.tolist()

# exponentially weighted average
moving_averages_exp = round(Re_series.ewm(alpha=0.005, adjust=False).mean(), 2)
Re_moving_averages_exp = moving_averages_exp.tolist()

Re_ave = round(np.mean(Re), 4)
Re_std = round(np.std(Re), 4)
Reynolds = 'Re = ' + str(Re_ave) + ' +/- ' + str(Re_std)

#print 'Time-averaged Magnetic Reynolds # =', Rm_ave, '+/-', Rm_std
print('Time-averaged Reynolds # =', Re_ave, '+/-', Re_std)


f.close()

#plt.rcParams["mathtext.fontset"] = "cm"
#plt.rcParams.update({'font.size': 10})

plt.plot(time, Re, label=r'Raw Data')
plt.plot(time, Re_moving_averages, '--', label=r'Cumulative Average')
#plt.plot(time, Re_moving_averages_exp, label=r'Exponential Average')
plt.xlabel(r'time')
plt.ylabel(r'Reynolds', rotation=90)

plt.legend()

plt.annotate(Reynolds, xy=(0.5, 0.95), xycoords='axes fraction')

plt.savefig('Reynolds_vs_time.png')
#plt.show()

