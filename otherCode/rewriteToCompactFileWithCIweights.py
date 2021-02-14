#pathToData = './COMPASOutput.h5'
from __future__ import print_function
from __future__ import division

import sys

import h5py  as h5   #for handling data format
import numpy as np  #for array handling
import os           #For checking existence data
import WriteH5File

import h5py as h5
import sys

sys.path.append('/Users/floorbroekgaarden/Projects/BHNS_project/Scripts')
import ClassCosmicIntegrator  as CI #Given settings and redshifts returns rates (2D arrays) Loads the data
import coencodeVarious        as CV
from PostProcessingScripts import * 

def maskTargetDCOsSTROOPWAFEL(DCOtype, boolDCOmask, f): #, otherSelection, otherparam):
    """returns mask of DCOs of interest
    fxparam  is hdf5 keyname of file where variable for which you want to mask DCOs is in 
    DCOtype = 'BBH' / 'ALL' / 'BHNS' or 'BNS' 
    boolDCOmask = [Hubble, RLOF, Optimistic] # boolean values whether to mask mergers in a HUbble time, 
    binaries that have RLOFSecondaryAfterCEE = True, and Optimistic binaries (i.e. optimisticCEFlag == 0)
    pathToDirectory is pathname to Directory where _oratory & _sampling directories are
    """
    
    Hubble, RLOF, Optimistic = boolDCOmask
    

 
    
    fDCO = f['doubleCompactObjects']
    
    
    
    # mask binaries of given DCO type
    if (DCOtype == 'BNS') | (DCOtype=='NSNS'):
        mask0 = ((fDCO['stellarType1'][...] == 13) & (fDCO['stellarType2'][...] == 13))
    elif (DCOtype == 'BHNS') | (DCOtype == 'NSBH'):
        mask0 = ((fDCO['stellarType1'][...] == 13) & (fDCO['stellarType2'][...] == 14)) | \
            ((fDCO['stellarType1'][...] == 14) & (fDCO['stellarType2'][...] == 13) )          
    elif (DCOtype == 'BBH') | (DCOtype=='BHBH'):
        mask0 = ((fDCO['stellarType1'][...] == 14) & (fDCO['stellarType2'][...] == 14))
    elif (DCOtype == 'all') | (DCOtype == 'ALL') :
        mask0 = ((fDCO['stellarType1'][...] == 14) | (fDCO['stellarType1'][...] == 13))
    else:
        print('error: DCO type not known')
        
    # Hubble mask
    if Hubble:
        mask1 = (fDCO['mergesInHubbleTimeFlag'][...]==True) 
    elif not Hubble:
        mask1 = (fDCO['mergesInHubbleTimeFlag'][...]==True) |  (fDCO['mergesInHubbleTimeFlag'][...]==False) 
    # RLOF mask
    if RLOF:
        mask2 = (fDCO['RLOFSecondaryAfterCEE'][...]==False)
    elif not RLOF:
        mask2 = (fDCO['RLOFSecondaryAfterCEE'][...]==False) | (fDCO['RLOFSecondaryAfterCEE'][...]==True)
    # Optimistic mask :  if True mask systems that have optimistic CE flag ==1
    if Optimistic:
        mask3 = np.logical_not(fDCO["optimisticCEFlag"][...] == 1) + \
        np.logical_not(fDCO["optimisticCEFlag"][...] == 0) 

    elif not Optimistic:
        mask3 = np.logical_not(fDCO["optimisticCEFlag"][...] == 1)  
    
    # combine the different masks and the oratory and refinement masks
    combinedmask = mask0 * mask1 * mask2 * mask3
    combinedmask = combinedmask.squeeze()
    # if otherSelection =='UFD':
    #     KpcToKM = 3.086 * 10**(16) # kpc to km  
    #     MyrToYr = 1E6 # yrs
    #     YrToSec = 3.154 *1E7 #sec        
    #     UFD_epsilon = otherparam[0]
    #     UFD_Rvir = otherparam[1]
    #     Xbh1 = otherparam[2]
    #     Rns = otherparam[3]

    #     fSN = f['supernovae']
    #     seedsOfIntererst = fDCO['seed'][...].squeeze()
    #     seedsSN = fSN['randomSeed'][...].squeeze()
    #     bools = np.in1d(seedsSN, seedsOfIntererst)        
        
    #     tc  = fDCO['tc'][...].squeeze()
    #     vsys = fSN['systemicVelocity'][...].squeeze()[bools]
    #     vsysSN2 = vsys[1:][::2]
    #     traveldistance = tc * vsysSN2 *  MyrToYr * YrToSec
    #     radiusUFDgalaxy = UFD_epsilon * UFD_Rvir * KpcToKM
    #     maskCandidatesUFD = (traveldistance <= radiusUFDgalaxy) | ((vsysSN2 <= 44) & (tc * MyrToYr *YrToSec<= radiusUFDgalaxy)) 
        
    #     combinedmask = maskCandidatesUFD*combinedmask
    

    return combinedmask


def reduceH5file(pathToData, pathToDataWithoutCOMPASname, DCOtype, optimistic, BPSmodelName):
    """

    DCOtype in 'BBH', 'BNS', 'BHNS'
     """

    ff  = h5.File(pathToData) # + 'COMPASOutput.h5')
    # print(pathToData + 'COMPASOutput.h5')
    print("The main files I have at my disposal are:\n",list(ff.keys()))


    # print("The main files I have at my disposal are:\n",list(Data['formationChannels'].keys()))

	# Which Files do I want?
	# options: ['RLOF', 'XRayBinaries', 'commonEnvelopes', 'cppSource', 'doubleCompactObjects', 'formationChannels', 'pulsarEvolution', 'runtimes', 'supernovae', 'systems'])
    filesOfInterest   = {1:'doubleCompactObjects',2:'systems',\
                             3:'supernovae', 4:'formationChannels', 5:'RLOF', 6:'commonEnvelopes'} #{1:'doubleCompactObjects',2:'systems', 3:'supernovae'}

	# #Give a list of columns you want, if you want all, say ['All']
	# columnsOfInterest = {1:['All'],\
	#                      2:['All'],\
	#                      3:['SEED', 'MZAMS_1', 'MZAMS_2']}
    columnsOfInterest =	   {1:['COCoreMassDCOFormation1', 'COCoreMassDCOFormation2',  'ECSNPrimary', 'ECSNSecondary', \
                                                         'HeCoreMassDCOFormation1', 'HeCoreMassDCOFormation2', 'ID', 'M1', 'M1ZAMS', 'M2', 'M2ZAMS', 'Metallicity1', 'Metallicity2', 'PISNPrimary', 'PISNSecondary', 'PPISNPrimary', 'PPISNSecondary', \
                                                         'PrimaryMTCase', 'RL1to2PostCEE', 'RL1to2PreCEE', 'RL2to1PostCEE', 'RL2to1PreCEE', 'RLOFSecondaryAfterCEE', 'SecondaryMTCase', 'SemiMajorAxisPostCEE', 'SemiMajorAxisPreCEE', 'USSNPrimary', 'USSNSecondary',\
                                                          'coreMassDCOFormation1', 'coreMassDCOFormation2', 'doubleCommonEnvelopeFlag', 'drawnKick1', 'drawnKick2', 'eccentricityDCOFormation', 'eccentricityInitial', 'eccentricityPrior2ndSN',\
                                                            'kickDirectionPower',  'mergesInHubbleTimeFlag', 'optimisticCEFlag', 'phiSupernova1', 'phiSupernova2',   'recycledPrimary', 'recycledSecondary', 'relativeVelocity2ndSN', 'samplingPhase', 'seed', \
                                                            'separationDCOFormation', 'separationInitial', 'separationPrior2ndSN', 'sigmaKickBH', 'sigmaKickNS', 'stellarType1', 'stellarType2', 'tc', 'tform', 'thetaSupernova1', 'thetaSupernova2', 'totalMassDCOFormation1', 'totalMassDCOFormation2', 'weight'],\
                                                        2:['ID', 'Metallicity1', 'Metallicity2', 'SEED', 'disbound', 'eccentricity',  'mass1', 'mass2', 'meanAnomaly1', 'meanAnomaly2', 'omega1', 'omega2', 'phi1', 'phi2', 'samplingPhase', 'separation', 'stellar_merger', 'theta1', 'theta2', 'weight'],\
                                                        3:['MassCOCoreSN', 'MassCoreSN', 'MassStarCompanion', 'MassStarSN',  'Survived','drawnKickVelocity', 'eccentricityAfter', 'eccentricityBefore', 'experiencedRLOF', 'fallback', 'flagECSN', 'flagHpoorSN', 'flagHrichSN', 'flagPISN', 'flagPPISN', 'flagRLOFontoaNS', 'flagSN', 'flagUSSN', 'kickVelocity', \
                                                                  'phi', 'previousStellarTypeCompanion', 'previousStellarTypeSN', 'psi', 'randomSeed', 'runawayFlag', 'separationAfter', 'separationBefore', 'systemicVelocity', 'theta', 'time', 'uK', 'vRel', 'whichStar'],\
                                                          4:['All'], 5:['All'], 6:['All'] \
                                                             }

	# #example of the seeds dictionary the actual one will be defined later
	# seedsOfInterest   = {1:None,\
	#                      2:None,\
	#                      3:None}
    
    # seedsDCO = ff['doubleCompactObjects']['seed'][()]
    seedsSystems = ff['systems']['SEED'][...].squeeze()
    seedsSN = ff['supernovae']['randomSeed'][...].squeeze()


    seedsFC = ff['formationChannels']['m_randomSeed'][...].squeeze()
    seedsRLOF = ff['RLOF']['randomSeed'][...].squeeze()
    seedsCE = ff['commonEnvelopes']['randomSeed'][...].squeeze()



    boolDCOmask = []
    boolDCOmask.append(1)
    boolDCOmask.append(1)
    boolDCOmask.append(optimistic)
    print('boolDCOmask', boolDCOmask)

    maskDCO = maskTargetDCOsSTROOPWAFEL(DCOtype=DCOtype, boolDCOmask=boolDCOmask, f=ff) #,\
                                            # otherSelection=None, otherparam=None)


    #what are the seeds of the systems that form BBHs
    seedsDCO = ff['doubleCompactObjects']['seed'][...].squeeze()[maskDCO]
    seedsSystemsMask = np.in1d(seedsSystems,seedsDCO)
    seedsSNMask = np.in1d(seedsSN, seedsDCO)

    seedsFCMask = np.in1d(seedsFC, seedsDCO)
    seedsRLOFMask = np.in1d(seedsRLOF, seedsDCO)
    seedsCEMask = np.in1d(seedsCE, seedsDCO)

    seedsSystems = seedsSystems[seedsSystemsMask]
    seedsSN = seedsSN[seedsSNMask]

    seedsFC = seedsFC[seedsFCMask]
    seedsRLOF = seedsRLOF[seedsRLOFMask]
    seedsCE = seedsCE[seedsCEMask]




    # print((seedsDCO))
    # print((seedsSystems))
    # print((seedsSN))   
    # print(len(seedsDCO))
    # print(len(seedsSystems))
    # print(len(seedsSN))

    seedsOfInterest   = {1:seedsDCO,\
                          2:seedsSystems,\
                          3:seedsSN,\
                          4:seedsFC,\
                          5:seedsRLOF,\
                          6:seedsCE}


    pathToNewData = pathToDataWithoutCOMPASname + 'COMPASCompactOutput_' + DCOtype + '_' + BPSmodelName +'.h5'
    # pathToNewData = pathToNewDataWithoutCOMPASname  + 'COMPASCompactOutput_' + DCOtype + '_' + BPSmodelName +'.h5' ## TEMP  


    print(pathToNewData)
    WriteH5File.reduceH5(pathToOld = pathToData, pathToNew = pathToNewData,\
                     dictFiles=filesOfInterest, dictColumns=columnsOfInterest, dictSeeds=seedsOfInterest)






    ff.close() 




    print('-----------------------------')
    print('now doing the second part (adding Cosmic Integration Weights)')





    


    #path to the data
    pathCOMPASOutput = pathToDataWithoutCOMPASname
    modelname = BPSmodelName
    #Will only look at BBHs so might as well set everything
    minz = 0.
    if DCOtype=='BHNS':
        maxz = .50
        resz = 100 # change to 100 //floor 
    elif DCOtype=='BNS':
        maxz = .25
        resz = 100 # change to 100 //floor 
    elif DCOtype=='BBH': 
        maxz = 1.5
        resz = 250 # change to 100 //floor 
    Data = CI.CosmicIntegrator(COMPASpath = pathCOMPASOutput, DCOtypes=DCOtype,\
           minRedshift=minz,   maxRedshift=maxz, nrRedshiftBins=resz, optimistic=optimistic, Cosmology='WMAP', COMPASbinaryFraction=1.0)

    #I use the custom cosmology because this was the flatlambda prescription used before WMAP Stevenson et al 2019
    #Doesnt matter to much (between WMAP and 
    #this it is 22, and 22.7 per year) but to prevent redoing all the numbers in the tex for referee

    print(Data.COMPAS.mass1)
    print(len(Data.COMPAS.mass1))


    #######
    fdata = h5.File(pathToNewData)

    Rdet = fdata.create_group("weights_detected")
    Rz0 = fdata.create_group("weights_intrinsic")

    RdetPerRedshift = fdata.create_group("weights_detectedPerRedshift")
    RintPerRedshift = fdata.create_group("weights_intrinsicPerRedshift")



    GSMFs = ['Panter et al. (2004) Single', 'Furlong et al. (2015) Single', 'Furlong et al. (2015) Double']
    MZs   = [ 'Langer et al. (2006)'      , 'Langer et al. +offset (2006)', 'Ma et al. (2015)']
    SFRs  = ['Madau et al. (2014)'         ,'Strolger et al. (2004)',     'Madau et al. (2017)']


    # Neijssel:
             
    Data.MSSFR.Zprescription         = 'logNormal' 
    Data.MSSFR.SFRprescription       = 'Neijssel et al. (2019)'
    Data.MSSFR.logNormalPrescription = 'Neijssel Phenomenological'
    Data.MSSFR.GSMFprescription      = None
    Data.MSSFR.ZMprescription        = None
    Data.cosmologicalIntegration()


    weightsDet        =  np.sum(Data.PerSystemPerRedshift_ratesObserved*Data.COMPAS.weight, axis=0) # //floor weight
    ratesPerSystem_z0 = Data.PerSystemPerRedshift_ratesIntrinsic[0,:] * Data.COMPAS.weight

    mssfrName='000'


    # append total observed rate (integrated over redshift horizon LIGO/VIRGO/KAGRA)
    Rdetxyz = Rdet.create_dataset(u"w_%s"%mssfrName, data=weightsDet)
    Rdetxyz.attrs[u"units"] = u"yr^{-1}"
    # append total intrinsic rate at redshift 0 
    Rdetxyz = Rz0.create_dataset(u"w_%s"%mssfrName, data=ratesPerSystem_z0)
    Rdetxyz.attrs[u"units"] = u"Gpc^{-3}yr^{-1}"    


    Redshifts      = Data.Shell_centerRedshift
    for ind_zz in range(len(Data.PerSystemPerRedshift_ratesIntrinsic)):
            weightsDet_zz = \
            Data.PerSystemPerRedshift_ratesObserved[ind_zz,:]*Data.COMPAS.weight

            ratesPerSystem_zz =\
            Data.PerSystemPerRedshift_ratesIntrinsic[ind_zz,:] * Data.COMPAS.weight

            zz = Redshifts[ind_zz]
            zz=np.round(zz,4)


            # append total observed rate (integrated over redshift horizon LIGO/VIRGO/KAGRA)
            Rdetxyz = RdetPerRedshift.create_dataset(u"R_%s_z_%s"%(mssfrName,zz), data=np.sum(weightsDet_zz))
            Rdetxyz.attrs[u"units"] = u"yr^{-1}"
            # append total intrinsic rate at redshift 0 
            Rdetxyz = RintPerRedshift.create_dataset(u"R_%s_z_%s"%(mssfrName,zz), data=np.sum(ratesPerSystem_zz))
            Rdetxyz.attrs[u"units"] = u"Gpc^{-3}yr^{-1}"             

        
        
        
          



    for ind_GSMF, GSMF in enumerate(GSMFs):
        ind_y = ind_GSMF + 1
        for ind_MZ, MZ in enumerate(MZs):
            ind_z = ind_MZ +1
            for ind_SFR, SFR in enumerate(SFRs):
                ind_x = ind_SFR+1
                
                mssfrName = '%s%s%s'%(ind_x, ind_y, ind_z)   
                print('now at mssfr xyz=', mssfrName)
                     
                Data.MSSFR.Zprescription         = 'MZ_GSMF'
                Data.MSSFR.SFRprescription       = SFR
    #                 Data.MSSFR.logNormalPrescription = logNormal[nrL]
                Data.MSSFR.GSMFprescription      = GSMF
                Data.MSSFR.ZMprescription        = MZ
                Data.cosmologicalIntegration()
                

                weightsDet        =  np.sum(Data.PerSystemPerRedshift_ratesObserved*Data.COMPAS.weight, axis=0) # //floor weight
                ratesPerSystem_z0 = Data.PerSystemPerRedshift_ratesIntrinsic[0,:] * Data.COMPAS.weight




                # append total observed rate (integrated over redshift horizon LIGO/VIRGO/KAGRA)
                Rdetxyz = Rdet.create_dataset(u"w_%s"%mssfrName, data=weightsDet)
                Rdetxyz.attrs[u"units"] = u"yr^{-1}"
                # append total intrinsic rate at redshift 0 
                Rdetxyz = Rz0.create_dataset(u"w_%s"%mssfrName, data=ratesPerSystem_z0)
                Rdetxyz.attrs[u"units"] = u"Gpc^{-3}yr^{-1}"         
                
                
                
                Redshifts      = Data.Shell_centerRedshift
                for ind_zz in range(len(Data.PerSystemPerRedshift_ratesIntrinsic)):
                        weightsDet_zz = \
                        Data.PerSystemPerRedshift_ratesObserved[ind_zz,:]*Data.COMPAS.weight

                        ratesPerSystem_zz =\
                        Data.PerSystemPerRedshift_ratesIntrinsic[ind_zz,:] * Data.COMPAS.weight

                        zz = Redshifts[ind_zz]
                        zz=np.round(zz,3)

                        # append total observed rate (integrated over redshift horizon LIGO/VIRGO/KAGRA)
                        Rdetxyz = RdetPerRedshift.create_dataset(u"R_%s_z_%s"%(mssfrName,zz), data=np.sum(weightsDet_zz))
                        Rdetxyz.attrs[u"units"] = u"yr^{-1}"
                        # append total intrinsic rate at redshift 0 
                        Rdetxyz = RintPerRedshift.create_dataset(u"R_%s_z_%s"%(mssfrName,zz), data=np.sum(ratesPerSystem_zz))
                        Rdetxyz.attrs[u"units"] = u"Gpc^{-3}yr^{-1}"  
    
    # add SEEDS to each file                   
    # Rdetxyz = RdetPerRedshift.create_dataset(u"SEED", data=seedsDCO)
    # Rdetxyz.attrs[u"units"] = u"#"    
    # Rdetxyz = RintPerRedshift.create_dataset(u"SEED", data=seedsDCO)
    # Rdetxyz.attrs[u"units"] = u"#"   
    Rdetxyz = Rdet.create_dataset(u"SEED", data=seedsDCO)
    Rdetxyz.attrs[u"units"] = u"#"              
    Rdetxyz = Rz0.create_dataset(u"SEED", data=seedsDCO)
    Rdetxyz.attrs[u"units"] = u"#"


    fdata.close() 

    print()
    print('-----------------------------------------------')
    print('completed')
    print('SUCCESS!!!! YES!')







# for DCOtype in ['BNS']:

#     pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/fiducial/COMPASOutput.h5'
#     pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/fiducial/'
#     optimistic=int(0)
#     BPSmodelName = 'A'
#     print()
#     print('-------------------------------------------------------')
#     print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
#     print('which is in this directory', pathToDataWithoutCOMPASname)
#     print('optimistic = ', optimistic)

#     reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)








for DCOtype in ['BHNS', 'BNS']:
    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/massTransferEfficiencyFixed_0_25/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/massTransferEfficiencyFixed_0_25/'
    optimistic=int(0)
    BPSmodelName = 'B'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


for DCOtype in ['BHNS', 'BNS']:
    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/massTransferEfficiencyFixed_0_5/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/massTransferEfficiencyFixed_0_5/'
    optimistic=int(0)
    BPSmodelName = 'C'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


for DCOtype in ['BHNS', 'BNS']:
    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/massTransferEfficiencyFixed_0_75/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/massTransferEfficiencyFixed_0_75/'
    optimistic=int(0)
    BPSmodelName = 'D'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



for DCOtype in ['BHNS', 'BNS']:

    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/unstableCaseBB/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/unstableCaseBB/'
    optimistic=int(0)
    BPSmodelName = 'E'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


for DCOtype in ['BHNS', 'BNS']:

    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/alpha0_5/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/alpha0_5/'
    optimistic=int(0)
    BPSmodelName = 'F'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



for DCOtype in ['BHNS','BNS']:

    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/alpha2_0/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/alpha2_0/'
    optimistic=int(0)
    BPSmodelName = 'G'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



for DCOtype in ['BHNS', 'BNS']:

    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/fiducial/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/fiducial/'
    optimistic=int(1)
    BPSmodelName = 'H'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


for DCOtype in ['BHNS','BNS']:

    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/rapid/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/rapid/'
    optimistic=int(0)
    BPSmodelName = 'I'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



for DCOtype in ['BHNS','BNS']:

    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/maxNSmass2_0/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/maxNSmass2_0/'
    optimistic=int(0)
    BPSmodelName = 'J'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



for DCOtype in ['BHNS','BNS']:

    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/maxNSmass3_0/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/maxNSmass3_0/'
    optimistic=int(0)
    BPSmodelName = 'K'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)


for DCOtype in ['BHNS','BNS']:

    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/noPISN/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/noPISN/'
    optimistic=int(0)
    BPSmodelName = 'L'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)



for DCOtype in ['BHNS','BNS']:

    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/ccSNkick_100km_s/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/ccSNkick_100km_s/'
    optimistic=int(0)
    BPSmodelName = 'M'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)

for DCOtype in ['BHNS','BNS']:

    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/ccSNkick_30km_s/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/ccSNkick_30km_s/'
    optimistic=int(0)
    BPSmodelName = 'N'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)

for DCOtype in ['BHNS','BNS']:

    pathToData = '/Volumes/Andromeda/DATA/AllDCO_bugfix/noBHkick/COMPASOutput.h5'
    pathToDataWithoutCOMPASname = '/Volumes/Andromeda/DATA/AllDCO_bugfix/noBHkick/'
    optimistic=int(0)
    BPSmodelName = 'O'
    print()
    print('-------------------------------------------------------')
    print(' now at ', DCOtype, 'running for modelname', BPSmodelName)
    print('which is in this directory', pathToDataWithoutCOMPASname)
    print('optimistic = ', optimistic)

    reduceH5file(pathToData=pathToData, pathToDataWithoutCOMPASname=pathToDataWithoutCOMPASname, DCOtype=DCOtype, optimistic=optimistic, BPSmodelName=BPSmodelName)













# if __name__ == "__main__":
#     pathToData = (sys.argv[1])
#     pathToDataWithoutCOMPASname  = (sys.argv[2])
#     DCOtype = (sys.argv[3])
#     optimistic = int(sys.argv[4])
#     BPSmodelName = (sys.argv[5])

#     print(optimistic)
    
# #    print('test')
#     reduceH5file(pathToData,  pathToDataWithoutCOMPASname,  DCOtype, optimistic, BPSmodelName)


### EXAMPLE:
#python reduceHdf5DataMinimal.py '/Volumes/Virgo/DATA/BHNS/alpha_10/COMPASOutput.h5' '/Volumes/Virgo/DATA/BHNS/alpha_10/'  'BNS' 0 'A'
# python rewriteToCompactFileWithCIweights.py '/Volumes/Andromeda/DATA/AllDCO/fiducial/COMPASOutput.h5' '/Volumes/Andromeda/DATA/AllDCO/fiducial/'  'BNS' 1 'B'
