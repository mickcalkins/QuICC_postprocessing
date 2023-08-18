""" 
    Reads in a stats file and plots the data.
    
    You need to provide name of variable.

    Dependencies: none

    Usage example: 
    
    python plot_stats.py velocityz_rms

"""

import sys
import numpy as np
import matplotlib.pyplot as plt



filename1 = 'dissB_avg_R80Q01.dat'
filename2 = 'dissB_avg_R80Q100.dat'
filename3 = 'velocityz_rms_R80Q01.dat'
filename4 = 'velocityz_rms_R80Q100.dat'

# open file
f1 = open(filename1, 'r')
f2 = open(filename2, 'r')
f3 = open(filename3, 'r')
f4 = open(filename4, 'r')

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
header2 = f4.readline()
header3 = f4.readline()

# initialize data objects as lists
data1 = []
time1 = []
data2 = []
time2 = []
data3 = []
time3 = []
data4 = []
time4 = []

# load line by line
# we separate time from data 
for line in f1:
    line = line.strip()  # get rid of \n character
    rows = line.split()  # breaks data into separate rows
    time1.append(rows[0]) # first data goes into time
    data1.append(rows[-len(rows)+1:])    
for line in f2:
    line = line.strip()  # get rid of \n character
    rows = line.split()  # breaks data into separate rows
    time2.append(rows[0]) # first data goes into time
    data2.append(rows[-len(rows)+1:])    
for line in f3:
    line = line.strip()  # get rid of \n character
    rows = line.split()  # breaks data into separate rows
    time3.append(rows[0]) # first data goes into time
    data3.append(rows[-len(rows)+1:])    
for line in f4:
    line = line.strip()  # get rid of \n character
    rows = line.split()  # breaks data into separate rows
    time4.append(rows[0]) # first data goes into time
    data4.append(rows[-len(rows)+1:])    

f1.close()
f2.close()
f3.close()
f4.close()


data1 = np.asarray(data1, dtype=np.float)
data2 = np.asarray(data2, dtype=np.float)
data3 = np.asarray(data3, dtype=np.float)
data4 = np.asarray(data4, dtype=np.float)

time = time1[1:len(time1)] # strip out hash symbol
num_times = len(data1)-1


data1_avg = np.flipud(np.mean(data1[2:], axis=0))
data2_avg = np.flipud(np.mean(data2[2:], axis=0))
data3_avg = np.flipud(np.mean(data3[2:], axis=0))
data4_avg = np.flipud(np.mean(data4[2:], axis=0))


plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams.update({'font.size': 20})
plt.rcParams["font.family"] = "Times New Roman"

#plt.plot(data1_avg, data1[0], 'k-', label=r'$\epsilon_B, \widetilde{Q}=0.1$')
#plt.plot(10.*data2_avg, data1[0], 'k--', label=r'$10 \times \epsilon_B, \widetilde{Q}=100$')

plt.plot(data3_avg, data1[0], 'k-', label=r'$\widetilde{Q}=0.1$')
plt.plot(data4_avg, data1[0], 'k--', label=r'$\widetilde{Q}=100$')

#plt.xlabel(r'$\epsilon_B$')
plt.xlabel(r'$w_{rms}$')
plt.ylabel(r'$z$', rotation=0)

legend = plt.legend(loc='best', shadow=False, ncol = 1, fontsize = 18, frameon=False) 

plt.tight_layout()

#plt.xticks(np.arange(0, 35, step=5.))
#plt.savefig('wrms_comparison.png')
#plt.savefig('wrms_comparison.png')
#plt.savefig('wrms_R80Q01.png')
#plt.savefig('dissB_R80Q100.png')
#plt.savefig('dissB_comparison.png')
#plt.savefig('dissV_R80Q01.png')

plt.show()
