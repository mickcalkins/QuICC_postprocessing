
import numpy as np
import h5py

data = h5py.File('visState0000.hdf5', 'r')

names = list(data)

print('here is what is inside', names)

print('specific output', names[3]) 

# these are both groups according to tests below
fluct_group = data[names[0]]
mesh_group = data[names[3]]


print(mesh_group)

#print(mesh_group[0].value)

# get the actual data with value attribute
fluct = fluct_group[names[0]].value


# gives max value of array
print(np.amax(fluct))

# print dimensions of array
print(fluct.shape)

print(isinstance(mesh_group, h5py.Group))

# check out contents of an item
item = mesh_group 

isFile = isinstance(item, h5py.File)
isGroup = isinstance(item, h5py.Group)
isDataset = isinstance(item, h5py.Dataset)


# get information about item
#print('item.id',item.id)
#print('item.ref',item.ref)
#print('item.parent',item.parent)
#print('item.file',item.file)
#print('item.name',item.name)

# get the list of daughters in the group
print(mesh_group.items())
mesh_names = list(mesh_group.items())
print(list(mesh_group[mesh_names[0]].value))

exit()
# convert the group in dictionary and iterate over their key and values
#for key,val in dict(data).iteritems():
#    print key,val

# this prints attributes of the file
#print data.attrs.keys()

#n1 = data.get('fluct_temperature')
#print np.array(n1)


# this is the most helpful/easiest to see data layout
def printname(name):
    print(name)

f.visit(printname)
