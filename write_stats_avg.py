""" 
    Reads in a stats file and writes the time-averaged data to a file.

    You need to provide name of variable

    dependencies: none

    usage example: 
    
    python write_stats_avg.py velocityz_rms
"""

import sys
import numpy as np
import matplotlib.pyplot as plt


variable = str(sys.argv[1])

filename = variable + '.dat'

# open file
f = open(filename, 'r')

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
data = np.asarray(data, dtype=np.float)
time = time[1:len(time)] # strip out hash symbol

# find the number of times data was output
num_times = len(data)-1


data_avg = np.flipud(np.mean(data[2:], axis=0))
data_rms = np.flipud(np.std(data[2:], axis=0))

#sys.exit()

plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams.update({'font.size': 20})
plt.rcParams["font.family"] = "Times New Roman"

#plt.plot(np.flipud(data[2:].T), data[0], color='tab:gray')
plt.plot(data_avg, data[0], 'k-')
#plt.plot(data_avg+data_rms, data[0], 'k--')
#plt.plot(data_avg-data_rms, data[0], 'k--')

#plt.xlabel(variable)
#plt.xlabel(r'$\overline{\mathbf{J}^2}$')
#plt.xlabel(r'$\overline{\bf{\zeta}^2}$')
plt.xlabel(r'$w_{rms}$')
plt.ylabel(r'$z$', rotation=0)

plt.tight_layout()

#plt.savefig('wrms_R80Q100.png')
plt.savefig('wrms_R80Q01.png')
#plt.savefig('dissB_R80Q100.png')
#plt.savefig('dissB_R80Q01.png')
#plt.savefig('dissV_R80Q01.png')

plt.show()
