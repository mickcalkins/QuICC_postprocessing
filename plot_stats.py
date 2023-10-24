""" 
    Reads in a stats file and plots the data
    You need to provide name of variable

    dependencies: none

    usage example: 
    
    python plot_stats.py velocityz_rms
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
data = np.asarray(data, dtype=float)
time = time[1:len(time)] # strip out hash symbol

# find the number of times data was output
num_times = len(data)-1

print('number of records =', num_times)

#sys.exit()

data_avg = np.flipud(np.mean(data[1:], axis=0))
data_rms = np.flipud(np.std(data[1:], axis=0))


plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams.update({'font.size': 18})

plt.plot(np.flipud(data[2:].T), data[0], color='tab:gray')
plt.plot(data_avg, data[0], 'k-')
plt.plot(data_avg+data_rms, data[0], 'k--')
plt.plot(data_avg-data_rms, data[0], 'k--')

plt.xlabel(variable)
plt.ylabel(r'$z$', rotation=0)

plt.show()
