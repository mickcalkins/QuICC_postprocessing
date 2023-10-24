

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc


df = pd.read_csv('energy_spectra.csv', delimiter='\t')

k = df.k
E_k = df.E_k

rc('text', usetex=True)
plt.rcParams.update({'axes.linewidth': 1.1})
plt.rcParams.update({'font.size': 18})



fig,ax = plt.subplots(1,1)

plt.plot(k, E_k, 'k')
plt.xlabel(r'$k$', labelpad=0)
plt.ylabel(r'$E_k$', rotation=0, labelpad=15)
plt.yscale('log')
plt.xscale('log')
ax.set_ylim(1e-4, 1e2) 
#plt.tick_params(axis='both', which='minor', length=minor_tick_length, width=tick_width)
#plt.tick_params(axis='both', which='major', length=tick_length, width=tick_width)

plt.tight_layout()
ax.tick_params(direction='in', which='both', top=True, right=True)

plt.savefig("energy_spectra.pdf", format='pdf', bbox_inches="tight")

plt.show()
