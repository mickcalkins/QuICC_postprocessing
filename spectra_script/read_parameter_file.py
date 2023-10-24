''' This code reads in QuICC parameters.cfg file and finds  
    certain parameters '''

#filename = 'parameters.cfg'

# open file
#global fh, filename
#fh = open(filename, 'r')

def read_box(filename):

    fh = open(filename, 'r')

    # names of things to find
    find_box = '/box2D'
    find_kc = '/kc2D'

    for num, line in enumerate(fh, 1):
        if find_box in line:
            data = line.lstrip()
            datanew = data.replace('<box2D>','')
            box2D = float(datanew.replace('</box2D>',''))
        if find_kc in line:
            data = line.lstrip()
            datanew = data.replace('<kc2D>','')
            kc2D = float(datanew.replace('</kc2D>',''))
    
    return box2D, kc2D

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
