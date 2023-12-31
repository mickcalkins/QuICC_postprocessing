

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
import sys


#df = pd.read_csv('energy_spectra.csv', delimiter='\t')
df = pd.read_csv('energy_spectrabaroclinic.csv', delimiter='\t')

k = df.k
E_k = df.E_k


#rc('text', usetex=True)
plt.rcParams.update({'axes.linewidth': 1.1})
plt.rcParams.update({'font.size': 18})



fig,ax = plt.subplots(1,1)

plt.plot(k[1:], E_k[1:], 'k')
plt.xlabel(r'$k$', labelpad=0)
plt.ylabel(r'$E_k$', rotation=0, labelpad=15)
plt.yscale('log')
plt.xscale('log')
#ax.set_ylim(1e-4, 1e2) 
#plt.tick_params(axis='both', which='minor', length=minor_tick_length, width=tick_width)
#plt.tick_params(axis='both', which='major', length=tick_length, width=tick_width)

plt.tight_layout()
ax.tick_params(direction='in', which='both', top=True, right=True)

plt.savefig("energy_spectra.png", format='png', bbox_inches="tight")

#sys.exit()
#plt.show()
