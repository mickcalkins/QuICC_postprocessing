""" 
    Computes thermal boundary layer thickness from rms temperature fluctuation. 
    A time average is performed with the data.

    dependencies: none

    usage example: 
    
    python find_tbl_thickness.py
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

#variable = str(sys.argv[1])

#filename = variable + '.dat'

# open file
f = open('temperature_rms.dat', 'r')

# read and ignore header lines
header1 = f.readline()
header2 = f.readline()
header3 = f.readline()

# initialize data objects as lists
data = []
time = []

# load line by line
# we separate time from data 
for line in f:
    line = line.strip()  # get rid of \n character
    rows = line.split()  # breaks data into separate rows
    time.append(rows[0]) # first data goes into time
    data.append(rows[-len(rows)+1:])    

f.close()

# convert to array: data[0] = z grid
data = np.asarray(data, dtype=float)
time = time[1:len(time)] # strip out hash symbol

# find the number of times data was output
num_times = len(data)-1
print('Number of writes = ', num_times)

temp_1 = np.flipud(np.mean(data[2:], axis=0))
temp_rms = np.flipud(np.std(data[2:], axis=0))
temp_2 = np.flipud(temp_1)
temp = 0.5*(temp_1 + temp_2)
#temp = np.insert(temp, 0, 0)
#temp = np.append(temp, 0)

z_out = data[0]
#z_out = np.insert(z_out, 0, 0)
#z_out = np.append(z_out, 1)

#f_interp = interpolate.interp1d(z_out, temp)

#fac = 20.0

#z_new = np.cos(np.pi*np.arange(0, fac*len(z_out)+1)/(fac*len(z_out)))
#z_new = 0.5*(1.0 - z_new)

#temp_new = f_interp(z_new)

#max_index = np.argmax(temp_new)

#z_max = z_new[max_index]

#if 0.5 < z_max < 1: 
#    tbl = 1.0 - z_max
#elif 0 < z_max < 0.5:
#    tbl = z_max 


#spl = Spline(z, temp_avg)

dtemp = np.gradient(temp,z_out)
ind_neg = np.where(dtemp < 0)[0][0]

slope = ( dtemp[ind_neg] - dtemp[ind_neg-1] )/( z_out[ind_neg] - z_out[ind_neg-1] )

z_cross = z_out[ind_neg-1] - (1/slope)*dtemp[ind_neg-1]

print('boundary layer thickness = ', z_cross)

sys.exit()

#plt.plot(np.flipud(data[2:].T), data[0], color='tab:gray')
plt.plot(temp, z_out, 'ko-')
plt.plot(temp+temp_rms, z_out, 'r--')
plt.plot(temp-temp_rms, z_out, 'r--')
#plt.plot(dtemp, z_out, 'o')

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
