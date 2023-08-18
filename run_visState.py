'''

This is a simple example that calls the read data subroutine.

'''

import read_visState as readvis

# provide the name of the file to be read
filename = 'visState0000'

# import relevant parameters/quantities/arrays

Ra, Pr = readvis.read_params_RBC(filename)

print 'Ra, Pr = ', Ra, Pr

x, y, z = readvis.read_grid(filename)

temp, mean_temp, u, v, w, vort_x, vort_y, vort_z = readvis.read_data_RBC(filename)

print 'temp', temp


