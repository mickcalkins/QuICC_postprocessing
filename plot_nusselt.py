""" 

Reads in a 'nusselt.dat' file and plots it.

"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


# initialize
time = []
Nu = []

# open file
f = open('nusselt.dat', 'r')
#f = open('nusselt_all.dat', 'r')

# read and ignore header lines
header1 = f.readline()
header2 = f.readline()
header3 = f.readline()

# load line by line
for line in f:
    line = line.strip()
    columns = line.split()
    time.append(columns[0])
    Nu.append(columns[1])

f.close()

# convert to array
time = np.asarray(time, dtype=np.double)
Nu = np.asarray(Nu, dtype=np.double)
Nu_series = pd.Series(Nu)
windows = Nu_series.expanding()
moving_averages = windows.mean()
Nu_moving_average = moving_averages.tolist()


Nu_ave = round(np.mean(Nu), 4)
Nu_std = round(np.std(Nu), 4)
Nusselt = 'Nu = ' + str(Nu_ave) + ' +/- ' + str(Nu_std)
print('Time-averaged Nusselt # =', Nu_ave, '+/-', Nu_std)

plt.plot(time,Nu)
plt.plot(time,Nu_moving_average, '--')

plt.xlabel('time')
plt.ylabel('Nusselt')

plt.annotate(Nusselt, xy=(0.5, 0.95), xycoords='axes fraction')

plt.savefig('Nusselt_vs_time.png')
#plt.show()

