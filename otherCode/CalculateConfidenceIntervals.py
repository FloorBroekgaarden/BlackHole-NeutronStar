# from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import time
import sys
import copy
#Quick fudge to make import from ../Scripts work
import sys
sys.path.append('../Scripts')

# from time import sleep
# from IPython.display import clear_output, display
import pandas as pd 
import string
# # sys.path.append(pathPostProcessing+'/2_CosmicIntegration')


# import ClassCosmicIntegrator  as CI #Given settings and redshifts returns rates (2D arrays) Loads the data
# import ClassEvents            as CE
# import ClassBayes             as CB
# import coencodeVarious        as CV
from PostProcessingScripts import * 


print(MSSFRnameslist)
DCOTypeList = ['BHBH', 'BHNS', 'NSNS']

NumberBPSmodels=15
alphabet = list(string.ascii_uppercase)
BPSnameslist = alphabet[:NumberBPSmodels]



        

def calculateConfidenceIntervals(BPSmodelName='A', MSSFRnameslist=MSSFRnameslist, DCOtype='BHNS', whichWeight='det'): 
    """ 
    Calculate confidence intervals  (distribution quantiles) summary 
    input:
    

    """


    # prepare DataFrame 
    xvarHeaders = ['Mass1', 'Mass2', 'tc',\
                   'log10(tc)', 'TotMass', 'ChirpMass', 'q', 'metallicitySystems', 'log10metallicitySystems', 'tdelay',\
                   'log10(tdelay)']

    xvarUnits = ['Msun', 'Msun', 'Gyr',\
                   'log10(Gyr)', 'Msun', 'Msun', '#', '#', '#', 'Gyr', 'log10(Gyr)']

    # quantiles that I want to know
    y_quantiles  =          [0.005,   0.05,   0.16,   0.25,   0.5,   0.75,   0.84,   0.95,  0.995]
    indexnames   = ['unit', '0.005', '0.05', '0.16', '0.25', '0.5', '0.75', '0.84', '0.95', '0.995']
    # nr of rows and columns that will be used:
    ncol_var = len(xvarHeaders)   
    ncol_MSSFR = len(MSSFRnameslist)
    ncol_Rate_det = 1

    nrows = len(y_quantiles) + 1 # +1 for units (see below)
    # store variables, and Observed and intrinsic rates for all MSSFR variations:
    ncol = ncol_var * (ncol_MSSFR) # for each MSSFR I want to give quantiles for each xparam 
    df_placeholder = np.zeros((nrows, ncol)) # will be filled in loop: 

    headernames=[]
    units=[]
    for ind_s, ss in enumerate(xvarHeaders):
        for ind_mssfr,  mssfr in enumerate(MSSFRnameslist):
            sss = ss + '_' + mssfr
            headernames.append(sss)
            units.append(xvarUnits[ind_s])

    # store dataFrame with zeros that we will fill on the go:
    dfw = pd.DataFrame(data=df_placeholder, columns=headernames, index=indexnames)   
    # add units at first row (index=0)
    dfw.iloc[0]=units        
        
                
        
        
    print('now at m=', BPSmodelName)










    if DCOtype=='BHNS':
        DCOname='BHNS'
    elif DCOtype=='BHBH':
        DCOname='BBH'
    elif DCOtype=='NSNS':
        DCOname='BNS'


    # constants
    Zsolar=0.0142
    nModels = 15
    # BPScolors       = sns.color_palette("husl", nModels)
    # lw = 3.5
    # Virgo         = '/Volumes/Virgo/DATA/BHNS/'
    # VirgoAllDCO = '/Volumes/Virgo/DATA/AllDCO/'
    # AndromedaBHNS = '/Volumes/Andromeda/DATA/BHNS/'
    # AndromedaAllDCO  = '/Volumes/Andromeda/DATA/AllDCO/'

    # alphabet = list(string.ascii_uppercase)
    # BPSnameslist = alphabet[:nModels]

    # BPSdir = ['fiducial/', 'fiducial/', 'alpha0_5/', 'alpha2_0/', 'unstableCaseBB/', 'rapid/', 'zeroBHkick/', 'massTransferEfficiencyFixed_0_25/', 'massTransferEfficiencyFixed_0_5/', 'massTransferEfficiencyFixed_0_75/', 'ccSNkick_100km_s/', 'ccSNkick_30km_s/']

    # dictBPSnameToDir   = dict(zip(BPSnameslist, BPSdir))    
    # dictBPSnameToColor = dict(zip(BPSnameslist, BPScolors))

    # path for files 
    path_dir = '/Volumes/Andromeda/DATA/AllDCO_bugfix/'
    nModels=15
    BPSnameslist = list(string.ascii_uppercase)[0:nModels]
    modelDirList = ['fiducial', 'massTransferEfficiencyFixed_0_25', 'massTransferEfficiencyFixed_0_5', 'massTransferEfficiencyFixed_0_75', \
                   'unstableCaseBB', 'alpha0_5', 'alpha2_0', 'fiducial', 'rapid', 'maxNSmass2_0', 'maxNSmass3_0', 'noPISN',  'ccSNkick_100km_s', 'ccSNkick_30km_s', 'noBHkick' ]

    alphabetDirDict =  {BPSnameslist[i]: modelDirList[i] for i in range(len(BPSnameslist))}



    #####



    path_ = path_dir
    path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOname +'_' + BPSmodelName + '.h5'
            



    # read in data 
    fdata = h5.File(path)



    # obtain BH and NS masses
    M1 = fdata['doubleCompactObjects']['M1'][...].squeeze()
    M2 = fdata['doubleCompactObjects']['M2'][...].squeeze()
    MBH, MNS = obtainM1BHandM2BHassymetric(M1, M2)
    del M1
    del M2

    tc = fdata['doubleCompactObjects']['tc'][...].squeeze() /1000.  # in Gyr
    metallicitySystems = fdata['doubleCompactObjects']['Metallicity1'][...].squeeze()
    tform = fdata['doubleCompactObjects']['tform'][...].squeeze() /1000.  # devide by 1000 to get units in Gyr 
    tdelay = tc + tform  # delay time 






    xvarlist = [MBH, MNS, tc, np.log10(tc), (MBH+MNS), chirpmass(MBH, MNS), (MBH/MNS), metallicitySystems, np.log10(metallicitySystems), tdelay, np.log10(tdelay)]
    # gives 90%, 99% and median

    # which weights do we want?
    if whichWeight =='det':
        fparam_name ='weights_detected'
    elif whichWeight=='z0':
        fparam_name = 'weights_intrinsic'



    for ind_mssfr, mssfr in enumerate(MSSFRnameslist):
        weightheader = 'w_' + mssfr
        weights = fdata[fparam_name][weightheader][...].squeeze()

        
        print('now at mssfr ', ind_mssfr+1, 'out of ', len(MSSFRnameslist))
        # should do bootstrapping here at some point.



        for ind_xvar, xvar in enumerate(xvarlist):

          Nrepeats=1
          boot_xvar = np.asarray(xvar)
          boot_index = np.arange(len(boot_xvar))

          bootstrap_array = np.zeros((Nrepeats, len(y_quantiles)))


          for ind_r, repeat in enumerate(range(Nrepeats)):
              # (bootstrap) re-sample a random set of samples with replacement from existing samples
              # do this by drawing random sample indecis nrs, each nr corresponds to a sample  
              boot_randindex = np.random.choice(boot_index, size=len(boot_index), replace=True, p=None)
              boot_randweight   = np.asarray(weights)[boot_randindex] # get the corresponding weights 
              boot_randsamples   = np.asarray(boot_xvar)[boot_randindex] # get the corresponding samples


              xvar_quantiles = weighted_quantile(values=boot_randsamples, quantiles=y_quantiles, \
                     sample_weight=boot_randweight)
              bootstrap_array[ind_r] = xvar_quantiles


          xvar_quantiles_bootstrapped = np.mean(bootstrap_array, axis=0)

          dfw_key = xvarHeaders[ind_xvar] + '_' + mssfr
          dfw[dfw_key][1:] = xvar_quantiles_bootstrapped

        
        
    dfCSVname = pathQuantiles + 'ConfidenceIntervals_model_' + BPSmodelName + '_' + DCOtype + '.csv' 
    dfw.to_csv(dfCSVname) 
    
    
    print()
    print('finished')

    fdata.close()
    
    return







pathQuantiles = '/Users/floorbroekgaarden/Projects/BlackHole-NeutronStar/csvFiles/distributionQuantiles/'
for DCOtype in ['BHNS']: #'BHBH', 'BHNS',
    print('---------------------------')
    print('---------------------------')
    print('---------------------------')
    print('at DCO model', DCOtype)
    for ind_m, BPSmodelName in enumerate(['A', 'B', 'C',  'D', 'E', 'F', 'G', 'H',  'I' ,'J', 'K', 'L', 'M', 'N', 'O']):
        print()
        print('------------------------------------------')
        calculateConfidenceIntervals(BPSmodelName=BPSmodelName, MSSFRnameslist=MSSFRnameslist, DCOtype=DCOtype, whichWeight='det')













# def calculateConfidenceIntervalsBoot(BPSmodelName='A', MSSFRnameslist=MSSFRnameslist, DCOtype='BHNS', whichWeight='det'): 
#     """ 
#     Calculate confidence intervals  
#     input:
    

#     """
    
#     # prepare DataFrame 
#     xvarHeaders = ['Mass1', 'Mass2', 'tc',\
#                    'log10(tc)', 'TotMass', 'ChirpMass', 'q', 'metallicitySystems']

#     xvarUnits = ['Msun', 'Msun', 'Gyr',\
#                    'Gyr', 'Msun', 'Msun', '#', '#']

#     # quantiles that I want to know
#     y_quantiles  =          [0.005,   0.05,   0.16,   0.25,   0.5,   0.75,   0.84,   0.95,  0.995]
#     indexnames   = ['unit', '0.005', '0.05', '0.16', '0.25', '0.5', '0.75', '0.84', '0.95', '0.995']
#     # nr of rows and columns that will be used:
#     ncol_var = len(xvarHeaders)   
#     ncol_MSSFR = len(MSSFRnameslist)
#     ncol_Rate_det = 1

#     nrows = len(y_quantiles) + 1 # +1 for units (see below)
#     # store variables, and Observed and intrinsic rates for all MSSFR variations:
#     ncol = ncol_var * (ncol_MSSFR) # for each MSSFR I want to give quantiles for each xparam 
#     df_placeholder = np.zeros((nrows, ncol)) # will be filled in loop: 

#     headernames=[]
#     units=[]
#     for ind_s, ss in enumerate(xvarHeaders):
#         for ind_mssfr,  mssfr in enumerate(MSSFRnameslist):
#             sss = ss + '_' + mssfr
#             headernames.append(sss)
#             units.append(xvarUnits[ind_s])

#     # store dataFrame with zeros that we will fill on the go:
#     dfw = pd.DataFrame(data=df_placeholder, columns=headernames, index=indexnames)   
#     # add units at first row (index=0)
#     dfw.iloc[0]=units        
        
                
        
        
#     print('now at m=', BPSmodelName)
#     # read in model
#     df = pd.read_csv('CompactData_model_' + BPSmodelName + '_' + DCOtype + '.csv', header=[0], skiprows=[1])


#     xvar1 = np.asarray(df['mass1'])
#     xvar2 = np.asarray(df['mass2'])
#     M1, M2 = obtainM1BHandM2BHassymetric(m1=xvar1, m2=xvar2)
#     tc = df['tc']/1000. # Gyr
#     metallicitySystems = df['metallicitySystems']


#     xvarlist = [M1, M2, tc, np.log10(tc), (M1+M2), chirpmass(M1, M2), (M1/M2), metallicitySystems]
#     # gives 90%, 99% and median

#     for ind_mssfr, mssfr in enumerate(MSSFRnameslist):
#         if whichWeight =='det':
#             dfkey_w = mssfr+'_R_det'
#         elif whichWeight=='z0':
#             dfkey_w = mssfr+'_R_z0'
#         weights = np.asarray(df[dfkey_w])
        
#         print('now at mssfr ', ind_mssfr+1, 'out of ', len(MSSFRnameslist))
#         # should do bootstrapping here at some point.



#         for ind_xvar, xvar in enumerate(xvarlist):

#           Nrepeats=100
#           boot_xvar = np.asarray(xvar)
#           boot_index = np.arange(len(boot_xvar))

#           bootstrap_array = np.zeros((Nrepeats, len(y_quantiles)))


#           for ind_r, repeat in enumerate(range(Nrepeats)):
#               # (bootstrap) re-sample a random set of samples with replacement from existing samples
#               # do this by drawing random sample indecis nrs, each nr corresponds to a sample  
#               boot_randindex = np.random.choice(boot_index, size=len(boot_index), replace=True, p=None)
#               boot_randweight   = np.asarray(weights)[boot_randindex] # get the corresponding weights 
#               boot_randsamples   = np.asarray(boot_xvar)[boot_randindex] # get the corresponding samples


#               xvar_quantiles = weighted_quantile(values=boot_randsamples, quantiles=y_quantiles, \
#                      sample_weight=boot_randweight)
#               bootstrap_array[ind_r] = xvar_quantiles


#           xvar_quantiles_bootstrapped = np.mean(bootstrap_array, axis=0)

#           dfw_key = xvarHeaders[ind_xvar] + '_' + mssfr
#           dfw[dfw_key][1:] = xvar_quantiles_bootstrapped

        
        
#     dfCSVname = 'ConfidenceIntervals_model_' + BPSmodelName + '_' + DCOtype + '_boot.csv' 
#     dfw.to_csv(dfCSVname) 
    
    
#     print()
#     print('finished')
    
#     return


# def calculateConfidenceIntervals(BPSmodelName='A', MSSFRnameslist=MSSFRnameslist, DCOtype='BHNS', whichWeight='det'): 
#     """ 
#     Calculate confidence intervals  
#     input:
    

#     """
    
#     # prepare DataFrame 
#     xvarHeaders = ['Mass1', 'Mass2', 'tc',\
#                    'log10(tc)', 'TotMass', 'ChirpMass', 'q', 'metallicitySystems']

#     xvarUnits = ['Msun', 'Msun', 'Gyr',\
#                    'Gyr', 'Msun', 'Msun', '#', '#']

#     # quantiles that I want to know
#     y_quantiles  =          [0.005,   0.05,   0.16,   0.25,   0.5,   0.75,   0.84,   0.95,  0.995]
#     indexnames   = ['unit', '0.005', '0.05', '0.16', '0.25', '0.5', '0.75', '0.84', '0.95', '0.995']
#     # nr of rows and columns that will be used:
#     ncol_var = len(xvarHeaders)   
#     ncol_MSSFR = len(MSSFRnameslist)
#     ncol_Rate_det = 1

#     nrows = len(y_quantiles) + 1 # +1 for units (see below)
#     # store variables, and Observed and intrinsic rates for all MSSFR variations:
#     ncol = ncol_var * (ncol_MSSFR) # for each MSSFR I want to give quantiles for each xparam 
#     df_placeholder = np.zeros((nrows, ncol)) # will be filled in loop: 

#     headernames=[]
#     units=[]
#     for ind_s, ss in enumerate(xvarHeaders):
#         for ind_mssfr,  mssfr in enumerate(MSSFRnameslist):
#             sss = ss + '_' + mssfr
#             headernames.append(sss)
#             units.append(xvarUnits[ind_s])

#     # store dataFrame with zeros that we will fill on the go:
#     dfw = pd.DataFrame(data=df_placeholder, columns=headernames, index=indexnames)   
#     # add units at first row (index=0)
#     dfw.iloc[0]=units        
        
                
        
        
#     print('now at m=', BPSmodelName)
#     # read in model
#     df = pd.read_csv('CompactData_model_' + BPSmodelName + '_' + DCOtype + '.csv', header=[0], skiprows=[1])


#     xvar1 = np.asarray(df['mass1'])
#     xvar2 = np.asarray(df['mass2'])
#     M1, M2 = obtainM1BHandM2BHassymetric(m1=xvar1, m2=xvar2)
#     tc = df['tc']/1000. # Gyr
#     metallicitySystems = df['metallicitySystems']


#     xvarlist = [M1, M2, tc, np.log10(tc), (M1+M2), chirpmass(M1, M2), (M1/M2), metallicitySystems]
#     # gives 90%, 99% and median

#     for ind_mssfr, mssfr in enumerate(MSSFRnameslist):
#         if whichWeight =='det':
#             dfkey_w = mssfr+'_R_det'
#         elif whichWeight=='z0':
#             dfkey_w = mssfr+'_R_z0'
#         weights = np.asarray(df[dfkey_w])
        
#         print('now at mssfr ', ind_mssfr+1, 'out of ', len(MSSFRnameslist))
#         # should do bootstrapping here at some point.



#         for ind_xvar, xvar in enumerate(xvarlist):

#           Nrepeats=1
#           boot_xvar = np.asarray(xvar)
#           boot_index = np.arange(len(boot_xvar))

#           bootstrap_array = np.zeros((Nrepeats, len(y_quantiles)))


#           for ind_r, repeat in enumerate(range(Nrepeats)):
#               # (bootstrap) re-sample a random set of samples with replacement from existing samples
#               # do this by drawing random sample indecis nrs, each nr corresponds to a sample  
#               boot_randindex = np.random.choice(boot_index, size=len(boot_index), replace=True, p=None)
#               boot_randweight   = np.asarray(weights)[boot_randindex] # get the corresponding weights 
#               boot_randsamples   = np.asarray(boot_xvar)[boot_randindex] # get the corresponding samples


#               xvar_quantiles = weighted_quantile(values=boot_randsamples, quantiles=y_quantiles, \
#                      sample_weight=boot_randweight)
#               bootstrap_array[ind_r] = xvar_quantiles


#           xvar_quantiles_bootstrapped = np.mean(bootstrap_array, axis=0)

#           dfw_key = xvarHeaders[ind_xvar] + '_' + mssfr
#           dfw[dfw_key][1:] = xvar_quantiles_bootstrapped

        
        
#     dfCSVname = 'ConfidenceIntervals_model_' + BPSmodelName + '_' + DCOtype + '.csv' 
#     dfw.to_csv(dfCSVname) 
    
    
#     print()
#     print('finished')
    
#     return
