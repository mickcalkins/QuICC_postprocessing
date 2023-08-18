""" 
    Computes Reynolds number for each velocity component.
    A time average is performed with the data.

    dependencies: none

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as Spline

#variable = str(sys.argv[1])

#filename = variable + '.dat'

# open files
f4 = open('velocityx_rms.dat', 'r')
f5 = open('velocityy_rms.dat', 'r')
f6 = open('velocityz_rms.dat', 'r')


# read and ignore header lines
header1 = f4.readline()
header2 = f5.readline()
header3 = f6.readline()

header1 = f4.readline()
header2 = f5.readline()
header3 = f6.readline()

header1 = f4.readline()
header2 = f5.readline()
header3 = f6.readline()

# initialize data objects as lists
u_x = []
u_y = []
u_z = []
time = []

# load line by line
# we separate time from data 

# load line by line
# we separate time from data 
for line in f4:
    line = line.strip()  # get rid of \n character
    rows = line.split()  # breaks data into separate rows
    time.append(rows[0]) # first data goes into time
    u_x.append(rows[-len(rows)+1:])    

f4.close()

for line in f5:
    line = line.strip()  # get rid of \n character
    rows = line.split()  # breaks data into separate rows
    time.append(rows[0]) # first data goes into time
    u_y.append(rows[-len(rows)+1:])    

f5.close()

for line in f6:
    line = line.strip()  # get rid of \n character
    rows = line.split()  # breaks data into separate rows
    time.append(rows[0]) # first data goes into time
    u_z.append(rows[-len(rows)+1:])    

f6.close()

# convert to array: data[0] = z grid
u_x = np.asarray(u_x, dtype=float)
u_y = np.asarray(u_y, dtype=float)
u_z = np.asarray(u_z, dtype=float)
z = u_x[0]

u_x_1 = np.flipud(np.mean(u_x[2:], axis=0))
u_y_1 = np.flipud(np.mean(u_y[2:], axis=0))
u_z_1 = np.flipud(np.mean(u_z[2:], axis=0))

u_x_2 = np.flipud(u_x_1)
ux = 0.5*(u_x_1 + u_x_2)
u_y_2 = np.flipud(u_y_1)
uy = 0.5*(u_y_1 + u_y_2)
u_z_2 = np.flipud(u_z_1)
uz = 0.5*(u_z_1 + u_z_2)

ux_avg = np.trapz(ux, z)
uy_avg = np.trapz(uy, z)
uz_avg = np.trapz(uz, z)

print('Reynolds-x = ', ux_avg)
print('Reynolds-y = ', uy_avg)
print('Reynolds-z = ', uz_avg)

sys.exit()

plt.plot(ux, z, 'k--')
plt.plot(uy, z, 'r--')
plt.plot(uz, z, 'b--')

#plt.plot(data_avg+data_rms, data[0], 'k--')
#plt.plot(data_avg-data_rms, data[0], 'k--')

#plt.xlabel(variable)
#plt.xlabel(r'$\overline{\mathbf{J}^2}$')
#plt.xlabel(r'$\overline{\bf{\zeta}^2}$')
plt.xlabel(r'$\theta_{rms}$')
plt.ylabel(r'$z$', rotation=0)

plt.tight_layout()

#plt.savefig('wrms_R80Q100.png')
#plt.savefig('wrms_R80Q01.png')
#plt.savefig('dissB_R80Q100.png')
#plt.savefig('dissB_R80Q01.png')
#plt.savefig('dissV_R80Q01.png')

plt.show()
