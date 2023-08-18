""" 

This file reads in a QuICC kinetic_energy.dat file and plots it.

"""

import numpy as np
import matplotlib.pyplot as plt
import read_parameter_file 
import read_visState 

# file name to read
filename = 'parameters.cfg'

# get run parameters
Ra, Pr = read_parameter_file.read_nondims(filename)
#Ra, Pr, Pm = read_nondims(filename)
#print 'Rayleigh number = ', '{:.2e}'.format(Ra)
#print 'Prandtl number = ', '{:.2e}'.format(Pr)
#print 'Magnetic Prandtl number = ', '{:.2e}'.format(Pm)

# get normalization parameters from file
box2D, box3D, kc2D, kc3D = read_parameter_file.read_box(filename)


# define normalization constant
fac = 0.5*(box2D/kc2D)**(2.)

# initialize
time = []
ke = []
ke_x = []
ke_y = []
ke_z = []

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
    ke_x.append(columns[2])
    ke_y.append(columns[3])
    ke_z.append(columns[4])

# convert to array
time = np.asarray(time, dtype=float)
#ke = fac*np.asarray(ke, dtype=float)
ke = fac*np.asarray(ke, dtype=float)
ke_x = fac*np.asarray(ke_x, dtype=float)
ke_y = fac*np.asarray(ke_y, dtype=float)
ke_z = fac*np.asarray(ke_z, dtype=float)
Re = np.sqrt(2.*ke)
Re_x = np.sqrt(2.*ke_x)
Re_y = np.sqrt(2.*ke_y)
Re_z = np.sqrt(2.*ke_z)

# compute basic stats
ke_ave = np.mean(ke)
ke_std = np.std(ke)

Re_ave = round(np.mean(Re), 4)
Re_std = round(np.std(Re), 4)


Re_x_ave = round(np.mean(Re_x), 4)
Re_x_std = round(np.std(Re_x), 4)

Re_y_ave = round(np.mean(Re_y), 4)
Re_y_std = round(np.std(Re_y), 4)

Re_z_ave = round(np.mean(Re_z), 4)
Re_z_std = round(np.std(Re_z), 4)

Reynolds = 'Re = ' + str(Re_ave) + ' +/- ' + str(Re_std)

print('Time-averaged Reynolds # =', Re_ave, '+/-', Re_std)
print('Time-averaged x-Reynolds # =', Re_x_ave, '+/-', Re_x_std)
print('Time-averaged y-Reynolds # =', Re_y_ave, '+/-', Re_y_std)
print('Time-averaged z-Reynolds # =', Re_z_ave, '+/-', Re_z_std)


f.close()

#plt.rcParams["mathtext.fontset"] = "cm"
#plt.rcParams.update({'font.size': 10})

plt.plot(time,Re, label=r'Re')
plt.plot(time,Re_x, label=r'$Re_{x}$')
plt.plot(time,Re_y, label=r'$Re_{y}$')
plt.plot(time,Re_z, label=r'$Re_{z}$')
#plt.plot(time,Re)
#plt.plot(time,Re_ave*np.ones(len(time)), 'k--')
#plt.plot(time,Re_ave*np.ones(len(time))+Re_std, 'k-.')
#plt.plot(time,Re_ave*np.ones(len(time))-Re_std, 'k-.')
plt.xlabel(r'$t$')
#plt.ylabel(r'$Re$', rotation=0)
plt.ylabel(r'$Re$', rotation=0)
legend = plt.legend(loc='best', ncol = 2) 

plt.annotate(Reynolds, xy=(0.5, 0.95), xycoords='axes fraction')

plt.savefig('Reynolds_vs_time.png')

#plt.show()

