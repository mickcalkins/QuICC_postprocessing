""" 
    Computes magnetic Taylor microscale.
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
f1 = open('fjx_rms.dat', 'r')
f2 = open('fjy_rms.dat', 'r')
f3 = open('fjz_rms.dat', 'r')
f4 = open('fbx_rms.dat', 'r')
f5 = open('fby_rms.dat', 'r')
f6 = open('fbz_rms.dat', 'r')


# read and ignore header lines
header1 = f1.readline()
header2 = f1.readline()
header3 = f1.readline()

header1 = f2.readline()
header2 = f2.readline()
header3 = f2.readline()

header1 = f3.readline()
header2 = f3.readline()
header3 = f3.readline()

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
vort_x = []
vort_y = []
vort_z = []
u_x = []
u_y = []
u_z = []
time = []

# load line by line
# we separate time from data 
for line in f1:
    line = line.strip()  # get rid of \n character
    rows = line.split()  # breaks data into separate rows
    time.append(rows[0]) # first data goes into time
    vort_x.append(rows[-len(rows)+1:])    

f1.close()

for line in f2:
    line = line.strip()  # get rid of \n character
    rows = line.split()  # breaks data into separate rows
    time.append(rows[0]) # first data goes into time
    vort_y.append(rows[-len(rows)+1:])    

f2.close()

for line in f3:
    line = line.strip()  # get rid of \n character
    rows = line.split()  # breaks data into separate rows
    time.append(rows[0]) # first data goes into time
    vort_z.append(rows[-len(rows)+1:])    

f3.close()

# convert to array: data[0] = z grid
vort_x = np.asarray(vort_x, dtype=float)
vort_y = np.asarray(vort_y, dtype=float)
vort_z = np.asarray(vort_z, dtype=float)
time = time[1:len(time)] # strip out hash symbol

# find the number of times data was output
num_times = len(vort_x)-1
#print('number of output = ', num_times)

vort_x_avg = np.flipud(np.mean(vort_x[2:], axis=0))
vort_y_avg = np.flipud(np.mean(vort_y[2:], axis=0))
vort_z_avg = np.flipud(np.mean(vort_z[2:], axis=0))
enstrophy_1 = vort_x_avg**2 + vort_y_avg**2 + vort_z_avg**2
enstrophy_2 = np.flipud(enstrophy_1)
enstrophy = 0.5*(enstrophy_1 + enstrophy_2)


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

u_x_avg = np.flipud(np.mean(u_x[2:], axis=0))
u_y_avg = np.flipud(np.mean(u_y[2:], axis=0))
u_z_avg = np.flipud(np.mean(u_z[2:], axis=0))
ke_1 = u_x_avg**2 + u_y_avg**2 + u_z_avg**2
ke_2 = np.flipud(ke_1)
ke = 0.5*(ke_1 + ke_2)

z = u_x[0]

ke_int = np.trapz(ke, z)
enstrophy_int = np.trapz(enstrophy, z)

taylor = np.sqrt(ke_int/enstrophy_int)

print('Ohmic dissipation rate/Q = ', enstrophy_int)
print('Magnetic Taylor microscale = ', taylor)

sys.exit()

plt.plot(ke, z, 'k--')
plt.plot(ke_1, z, 'r--')
plt.plot(ke_2, z, 'b--')

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
