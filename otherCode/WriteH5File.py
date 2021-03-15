import datetime
import h5py  as h5
import numpy as np
import os



def sanityChecks(h_old, filename, columns, seeds):

    data = h_old[filename]

    if columns == ['All']:
        columns = list(data.keys())

    for column in columns:
        try:
            test = data[column]
        except:
            raise ValueError("column %s does not exist in %s"\
                             %(column, filename))
    
    if filename=='doubleCompactObjects':
        seedheader='seed'
    elif filename=='systems':
        seedheader='SEED'
    elif filename=='supernovae':
        seedheader='randomSeed'
    elif filename=='formationChannels':
        seedheader='m_randomSeed'
    elif filename=='RLOF':
        seedheader='randomSeed'
    elif filename=='commonEnvelopes':
        seedheader='randomSeed'

    else:
        seedheader='SEED'
        # print('filename = ', filename)

    seedData = data[seedheader][()]

    mask      = np.in1d(seeds, seedData)
    not_exist = np.logical_not(mask)
    if np.sum(not_exist) !=0:
        raise ValueError("seed(s) %s do not exist in %s"\
                         %(seeds[not_exist], filename))


def createDataInNewH5(h_old, h_new, filename, columns, seeds):

    h_new.create_group(filename)

    dataOld = h_old[filename]

    if filename=='doubleCompactObjects':
        seedheader='seed'
    elif filename=='systems':
        seedheader='SEED'
    elif filename=='supernovae':
        seedheader='randomSeed'
    elif filename=='formationChannels':
        seedheader='m_randomSeed'
    elif filename=='RLOF':
        seedheader='randomSeed'
    elif filename=='commonEnvelopes':
        seedheader='randomSeed'
    else:
        print('filename = ', filename)


    mask    = np.in1d(dataOld[seedheader][()],seeds )

    if columns == ['All']:
        columns = list(dataOld.keys())

    for column in columns:
        print(column)
        # unit       = dataOld[column].attrs['units']
        columnOld  = dataOld[column][()]
        data       = columnOld[mask]
        dataNew    = h_new[filename].create_dataset(column, data=data)
        # dataNew.attrs['units'] = unit
        #add attribute comment
    
def reduceH5(pathToOld = None, pathToNew = None, dictFiles=None, dictColumns=None, dictSeeds=None):


        if pathToOld is None:
            raise ValueError("pathToOld not given")
        if pathToNew is None:
            raise ValueError("pathToNew not given")
        if dictFiles is None:
            raise ValueError("dictionary of files not given is None")
        if dictColumns is None:
            raise ValueError("dictionary of columns not given is None")
        if dictColumns is None:
            raise ValueError("dictionary of seeds not given is None")

        if((len(dictFiles) != len(dictColumns)) or\
           (len(dictSeeds) != len(dictColumns)) or\
           (len(dictFiles) != len(dictSeeds))):
            raise ValueError("the 3 dictionaries are not of same lengt!")

        h_old = h5.File(pathToOld)
        h_new = h5.File(pathToNew, 'w')


        keys = dictFiles.keys()

        for key in keys:
            fileName = dictFiles[key]          
            sanityChecks(h_old, dictFiles[key], dictColumns[key], dictSeeds[key])
            createDataInNewH5(h_old, h_new, dictFiles[key], dictColumns[key], dictSeeds[key])
            
        h_old.close()
        h_new.close()
