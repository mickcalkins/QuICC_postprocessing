''' This code reads in QuICC parameters.cfg file and finds  
    certain parameters '''

#filename = 'parameters.cfg'

# open file
#global fh, filename
#fh = open(filename, 'r')

def read_box(filename):

    fh = open(filename, 'r')

    # names of things to find
    find_box2 = '/box2D'
    find_box3 = '/box3D'
    find_kc2 = '/kc2D'
    find_kc3 = '/kc3D'

    for num, line in enumerate(fh, 1):
        if find_box2 in line:
            data = line.lstrip()
            datanew = data.replace('<box2D>','')
            box2D = float(datanew.replace('</box2D>',''))
        if find_box3 in line:
            data = line.lstrip()
            datanew = data.replace('<box3D>','')
            box3D = float(datanew.replace('</box3D>',''))
        if find_kc2 in line:
            data = line.lstrip()
            datanew = data.replace('<kc2D>','')
            kc2D = float(datanew.replace('</kc2D>',''))
        if find_kc3 in line:
            data = line.lstrip()
            datanew = data.replace('<kc3D>','')
            kc3D = float(datanew.replace('</kc3D>',''))
    
    return box2D, box3D, kc2D, kc3D

def read_res(filename):

    fh = open(filename, 'r')

    # names of things to find
    find_dim2 = '/dim2D'
    find_dim3 = '/dim3D'

    for num, line in enumerate(fh, 1):
        if find_dim2 in line:
            data = line.lstrip()
            datanew = data.replace('<dim2D>','')
            dim2D = float(datanew.replace('</dim2D>',''))
        if find_dim3 in line:
            data = line.lstrip()
            datanew = data.replace('<dim3D>','')
            dim3D = float(datanew.replace('</dim3D>',''))
      
    return dim2D, dim3D

def read_nondims(filename):

    fh = open(filename, 'r')

    # names of things to find
    find_Ra = '/rayleigh'
    find_Pr = '/prandtl'
    #find_Pm = '/magnetic_prandtl'

    for num, line in enumerate(fh, 1):
        if find_Ra in line:
            data = line.lstrip()
            datanew = data.replace('<rayleigh>','')
            Ra = float(datanew.replace('</rayleigh>',''))
        if find_Pr in line:
            data = line.lstrip()
            datanew = data.replace('<prandtl>','')
            Pr = float(datanew.replace('</prandtl>',''))
        #if find_Pm in line:
            #data = line.lstrip()
            #datanew = data.replace('<magnetic_prandtl>','')
            #Pm = float(datanew.replace('</magnetic_prandtl>',''))
    
    return Ra, Pr

def read_nondims_RRBC(filename):

    fh = open(filename, 'r')

    # names of things to find
    find_Ra = '/rayleigh'
    find_Pr = '/prandtl'
    find_Ek = '/ekman'

    for num, line in enumerate(fh, 1):
        if find_Ra in line:
            data = line.lstrip()
            datanew = data.replace('<rayleigh>','')
            Ra = float(datanew.replace('</rayleigh>',''))
        if find_Pr in line:
            data = line.lstrip()
            datanew = data.replace('<prandtl>','')
            Pr = float(datanew.replace('</prandtl>',''))
        if find_Ek in line:
            data = line.lstrip()
            datanew = data.replace('<ekman>','')
            Ek = float(datanew.replace('</ekman>',''))
    
    return Ra, Pr, Ek
