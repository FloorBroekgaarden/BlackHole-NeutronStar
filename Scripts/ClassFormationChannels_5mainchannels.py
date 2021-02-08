from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import time
import sys
import gc
import string

import ClassCOMPAS     as CC ###
from PostProcessingScripts import *
import math # for tan function and degrees function





dictDCOtypeDCOlabel = {'BBH':'BHBH', 'BNS':'NSNS', 'BHNS':'BHNS'}

def returnSeedsPercentageClassicAndOnlyStableMT(pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/fiducial/',\
                                types='BHNS',  withinHubbleTime=True, optimistic=False, \
                                binaryFraction=1):
    
    
    #The first step is to look at the number of systems interacting at this metallicity
    #at the time of simulating RLOF and CEE printing was not robust yet
    # to catch all the systems. Some MS and touching systems might
    #not be there but should not be too many

    # retrive hdf5 file 
    # f          = h5.File(pathCOMPASOutput+'COMPASOutput.h5')
    f          = h5.File(pathCOMPASOutput)

    # obtain the BHNS seeds: 
    Data            = CC.COMPASData(path=pathCOMPASOutput, lazyData=True, Mlower=5., \
                     Mupper=150., binaryFraction=binaryFraction)
    Data.setCOMPASDCOmask(types=types,  withinHubbleTime=withinHubbleTime, optimistic=optimistic)
    Data.setCOMPASData()
    
    SeedsBHNS    = Data.seeds


    
    allSeedsSystems  = f['systems']['SEED'][...].squeeze()
    weights    = f['systems']['weight'][...].squeeze()
    maskSystemsBHNS = np.in1d(allSeedsSystems, SeedsBHNS)
    seedsSystemBHNS = allSeedsSystems[maskSystemsBHNS]


    normalisation = float(np.sum(weights[maskSystemsBHNS]))
    # print('normalisation (weighted rate DCO) = ', normalisation)
    # print('# DCOs of interest = ', len(weights[maskSystemsBHNS]))


    #which seeds are interacting (Z doesnt matter)
    RLOFseeds  = np.unique(f['RLOF']['randomSeed'][...].squeeze())
    CEEseeds   = np.unique(f['commonEnvelopes']['randomSeed'][...].squeeze())
#     print RLOFseeds.shape, CEEseeds.shape
    InteractingSeeds = np.unique(np.concatenate((RLOFseeds, CEEseeds)))


    #Which of the Seeds in SeedsSystems are interacting?
    boolInteracting  = np.in1d(seedsSystemBHNS, InteractingSeeds)
    # print 'of all BHNS seeds %s percent experiencedMT'\
    #     %((np.sum(weights[maskSystemsBHNS][boolInteracting])/normalisation)*100)
    #from here on we use the seedsOf interest to do the next step in the formation channel
    seedsInterest = allSeedsSystems[maskSystemsBHNS][boolInteracting]    
    #these arrays are long delete and clear memory
    # del maskZ
    # del allSeedsSystems
    # del CEEseeds
    # del RLOFseeds
    # del boolInteracting
    # del InteractingSeeds
    #garbage collect forces deleted arrays to clear from memory
    gc.collect()

    
    
    # During this research the RLOF printing was very much in development
    # This means that some output such as the x-th moment of RLOF was 
    # not robust (yet). Hence we use a clever trick with the seeds
    # Credits  idea to J.W.Barrett, made it more general to get not only
    # firt RLOF per system but x-th. Snipper is in coencodeVarious.cp

    #carefull the cosmicInt simulation is large and this will 
    #take a lot of RAM-memory - be patient and hope (fingerscrossed)
    #array of 90.000.000 for allSeeds RLOF possibly more for RLOF
    #close other windows :O X'(

    #UPDATE we only want of metallicity Z
    RLOFseeds   = f['RLOF']['randomSeed'][...].squeeze()
    RLOF_Z_MASK = np.in1d(RLOFseeds,seedsInterest)

    # dont care about 4th and 5th moment etc will take longer to
    # calculate those
    nthMoment   = getXmomentOfMT(RLOFseeds, maxCounter=3)
    # print(len(RLOFseeds))


    gc.collect()
    
    ############################################
    #############   MASS TRANSFER   ############
    ############################################
    
    # percentage of systems post MS primaries
    # MS transfering onto MS secondary  (channel 1 and 2)
    donorPostMS = ((f['RLOF']['type1Prev'][...].squeeze() >=1) &\
                  (f['RLOF']['type1Prev'][...].squeeze() <7))
    companionMS = ((f['RLOF']['type2Prev'][...].squeeze() == 1) |\
                  (f['RLOF']['type2Prev'][...].squeeze() == 0))
    primary     = (f['RLOF']['flagRLOF1'][...].squeeze() == True)
    firstMoment = (nthMoment == 1)



    # # Check which RLOFseeds are BHNS systems
    # Data        = CC.COMPASData(path=pathCOMPASOutput, lazyData=True, Mlower=5., \
    #                  Mupper=150., binaryFraction=1)
    # Data.setCOMPASDCOmask(types='BHNS',  withinHubbleTime=True, optimistic=False)
    # Data.setCOMPASData()
    # SeedsHubble = Data.seeds[Data.Hubble==True]
    DCOmask     = np.in1d(RLOFseeds,  seedsSystemBHNS)



    maskInterest      = [donorPostMS & companionMS & firstMoment & primary & DCOmask]
    seedsRemain       = f['RLOF']['randomSeed'][...].squeeze()[maskInterest]
    # print(seedsRemain)



    # weights    = f['systems']['weight'][...].squeeze()
    # systemSeed = f['systems']['SEED'][...].squeeze()

    # maskSystemsBHNS 
    # DCOmaskSystems     = np.in1d(systemSeed, SeedsHubble)
    # normalisation = float(np.sum(weights[maskSystemsBHNS]))
    systemsOfInterest  = np.in1d(allSeedsSystems, seedsRemain)


    nrSystemsInterest = np.sum(weights[systemsOfInterest])
    # print "post MS primar onto MS secondary "
    # print "of interacting systems %s percentage"\
    #       %((nrSystemsInterest/normalisation)*100)
    # del weights
    # del systemSeed
    # del systemsOfInterest

    # del donorPostMS
    # del companionMS
    # del firstMoment
    # del primary
    # del maskInterest
    gc.collect()





    ######################
    ######################
    ######################




    # test the number of stable mass transfers.




    donorPostMS = ((f['commonEnvelopes']['stellarType1'][...].squeeze() >=1) &\
                  (f['commonEnvelopes']['stellarType1'][...].squeeze() <7))
    companionMS = (f['commonEnvelopes']['stellarType2'][...].squeeze() == 1) |\
                  (f['commonEnvelopes']['stellarType2'][...].squeeze() == 0)
    seedsCEE     = f['commonEnvelopes']['randomSeed'][...].squeeze()
    seeds        = np.in1d(f['commonEnvelopes']['randomSeed'][...].squeeze(), seedsRemain)

    DCOmaskCEE    = np.in1d(seedsCEE, seedsSystemBHNS)

    maskInterest = donorPostMS & companionMS & seeds & DCOmaskCEE

    seedsCEE    = f['commonEnvelopes']['randomSeed'][...].squeeze()\
                                    [maskInterest]


    seedsRemain = seedsRemain[np.logical_not(np.in1d(seedsRemain, seedsCEE))]

    weights    = f['systems']['weight'][...].squeeze()
    systemSeed = f['systems']['SEED'][...].squeeze()
    systemsOfInterest  = np.in1d(allSeedsSystems, seedsRemain)




    nrStable     = float(np.sum(weights[systemsOfInterest]))
    # print "percentage stable after 1st MT = %s"  %((nrStable/normalisation)*100)
#     del seedsCEE
#     del donorPostMS
#     del companionMS
#     del maskInterest
    gc.collect()    
    
    
    
    #########################
    ####### 2nd MT ##########
    ########################
    
    
    #nr systems second mass transfer 
    secondMoment = (nthMoment == 2)[RLOF_Z_MASK]
    secondary    = (f['RLOF']['flagRLOF2'][...].squeeze() == True)[RLOF_Z_MASK]
    seeds        = np.in1d(f['RLOF']['randomSeed'][...].squeeze()[RLOF_Z_MASK], seedsRemain)
    maskInterest = secondMoment & secondary & seeds
    seedsRemain  = f['RLOF']['randomSeed'][...].squeeze()[RLOF_Z_MASK][maskInterest]

    systemsOfInterest  = np.in1d(allSeedsSystems, seedsRemain)


    # print "percentage having a second MT from secondary  = %s"  %((np.sum(weights[systemsOfInterest])/normalisation)*100)
    del secondMoment
    del secondary
    del seeds
    gc.collect()
    # print seedsRemain[0]    

    
    #################. CE #############
    # test the number of stable mass transfers.
    donorUnstripped = (f['commonEnvelopes']['stellarType2'][...].squeeze() <7)
    #donor not allowed to be a HG star
    seeds           = np.in1d(f['commonEnvelopes']['randomSeed'][...].squeeze(), seedsRemain)
    maskInterest    = seeds & donorUnstripped 
    
    
    
    ######### CLASSIC #############
    seedsRemain_classic        = f['commonEnvelopes']['randomSeed'][...].squeeze()[maskInterest]
    systemsOfInterest_classic  = np.in1d(allSeedsSystems, seedsRemain_classic)
    # print "percentage unstable second MT =  %s"  %((np.sum(weights[systemsOfInterest_classic])/normalisation)*100)
    
    seedsChannel_classic = seedsRemain_classic
    percentageChannel_classic = ((np.sum(weights[systemsOfInterest_classic])/normalisation)*100)
    
    ######################
    
    
    
    ######### ONLY STABLE MT #############
    seedsCEE    = f['commonEnvelopes']['randomSeed'][...].squeeze()[maskInterest]
    seedsRemain_onlyStableMT = seedsRemain[np.logical_not(np.in1d(seedsRemain, seedsCEE))]


    systemsOfInterest_onlyStableMT  = np.in1d(allSeedsSystems, seedsRemain_onlyStableMT)
    # print "percentage unstable second MT =  %s"  %((np.sum(weights[systemsOfInterest_onlyStableMT])/normalisation)*100)
    
    seedsChannel_onlyStableMT = seedsRemain_onlyStableMT
    percentageChannel_onlyStableMT = ((np.sum(weights[systemsOfInterest_onlyStableMT])/normalisation)*100)    
    
    ######################
    print('percentage Classic        = ', percentageChannel_classic)
    print('percentage Only stable MT =  ', percentageChannel_onlyStableMT)
    
    ######### SAVE SEEDS AND PERCENTAGES #############
    seedsPercentageClassic = seedsChannel_classic, percentageChannel_classic
    seedsPercentageOnlyStableMT = seedsChannel_onlyStableMT, percentageChannel_onlyStableMT
    
    return seedsPercentageClassic, seedsPercentageOnlyStableMT








def returnSeedsPercentageDoubleCoreCEE(pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/fiducial/',\
                                types='BHNS',  withinHubbleTime=True, optimistic=False, \
                                binaryFraction=1):
    
    
    #The first step is to look at the number of systems interacting at this metallicity
    #at the time of simulating RLOF and CEE printing was not robust yet
    # to catch all the systems. Some MS and touching systems might
    #not be there but should not be too many

    # retrive hdf5 file 
    # f          = h5.File(pathCOMPASOutput+'COMPASOutput.h5')
    f          = h5.File(pathCOMPASOutput)

    # obtain the BHNS seeds: 
    Data            = CC.COMPASData(path=pathCOMPASOutput, lazyData=True, Mlower=5., \
                     Mupper=150., binaryFraction=binaryFraction)
    Data.setCOMPASDCOmask(types=types,  withinHubbleTime=withinHubbleTime, optimistic=optimistic)
    Data.setCOMPASData()
    
    SeedsBHNS    = Data.seeds


    
    allSeedsSystems  = f['systems']['SEED'][...].squeeze()
    weights    = f['systems']['weight'][...].squeeze()
    maskSystemsBHNS = np.in1d(allSeedsSystems, SeedsBHNS)
    seedsSystemBHNS = allSeedsSystems[maskSystemsBHNS]


    normalisation = float(np.sum(weights[maskSystemsBHNS]))
    # print('normalisation (weighted rate DCO) = ', normalisation)
    # print('# DCOs of interest = ', len(weights[maskSystemsBHNS]))


    #which seeds are interacting (Z doesnt matter)
    RLOFseeds  = np.unique(f['RLOF']['randomSeed'][...].squeeze())
    CEEseeds   = np.unique(f['commonEnvelopes']['randomSeed'][...].squeeze())
#     print RLOFseeds.shape, CEEseeds.shape
    InteractingSeeds = np.unique(np.concatenate((RLOFseeds, CEEseeds)))


    #Which of the Seeds in SeedsSystems are interacting?
    boolInteracting  = np.in1d(seedsSystemBHNS, InteractingSeeds)
    # print 'of all BHNS seeds %s percent experiencedMT'\
    #     %((np.sum(weights[maskSystemsBHNS][boolInteracting])/normalisation)*100)
    #from here on we use the seedsOf interest to do the next step in the formation channel
    seedsInterest = allSeedsSystems[maskSystemsBHNS][boolInteracting]    
    #these arrays are long delete and clear memory
    # del maskZ
    # del allSeedsSystems
    # del CEEseeds
    # del RLOFseeds
    # del boolInteracting
    # del InteractingSeeds
    #garbage collect forces deleted arrays to clear from memory
    gc.collect()

    
    
    # During this research the RLOF printing was very much in development
    # This means that some output such as the x-th moment of RLOF was 
    # not robust (yet). Hence we use a clever trick with the seeds
    # Credits  idea to J.W.Barrett, made it more general to get not only
    # firt RLOF per system but x-th. Snipper is in coencodeVarious.cp

    #carefull the cosmicInt simulation is large and this will 
    #take a lot of RAM-memory - be patient and hope (fingerscrossed)
    #array of 90.000.000 for allSeeds RLOF possibly more for RLOF
    #close other windows :O X'(

    #UPDATE we only want of metallicity Z
    RLOFseeds   = f['RLOF']['randomSeed'][...].squeeze()
    RLOF_Z_MASK = np.in1d(RLOFseeds,seedsInterest)

    # dont care about 4th and 5th moment etc will take longer to
    # calculate those
    nthMoment   = getXmomentOfMT(RLOFseeds, maxCounter=3)
    # print(len(RLOFseeds))


    gc.collect()
    
    ############################################
    #############   MASS TRANSFER   ############
    ############################################
    

    # percentage of systems post MS primaries
    # MS transfering onto MS secondary  (channel 1 and 2)
    primaryPostHG = (f['RLOF']['type1Prev'][...].squeeze() >2)[RLOF_Z_MASK]
    secondaryPostHG = (f['RLOF']['type2Prev'][...].squeeze() >2)[RLOF_Z_MASK]
    firstMoment = (nthMoment == 1)[RLOF_Z_MASK]

    DCOmask     = np.in1d(RLOFseeds[RLOF_Z_MASK],  seedsSystemBHNS)
    
    maskInterest      = [primaryPostHG & secondaryPostHG & firstMoment & DCOmask]
    seedsRemain       = f['RLOF']['randomSeed'][...].squeeze()[RLOF_Z_MASK][maskInterest]
    
    
    systemsOfInterest  = np.in1d(allSeedsSystems, seedsRemain)
    nrSystemsInterest = np.sum(weights[systemsOfInterest])
    

    testMoment   = getXmomentOfMT(seedsRemain, maxCounter=10)
    # print np.max(testMoment)


    # print "post MS primar onto MS secondary "
    # print "of interacting systems %s percentage"\
    #       %((nrSystemsInterest/normalisation)*100)

    gc.collect()
    
    
    
    # both stars are post MS
    primaryPostHG = (f['commonEnvelopes']['stellarType1'][...].squeeze() >2) &\
                    (f['commonEnvelopes']['stellarType1'][...].squeeze() <7)
    secondaryPostHG = (f['commonEnvelopes']['stellarType2'][...].squeeze()>2) &\
                    (f['commonEnvelopes']['stellarType2'][...].squeeze() <7)
    seeds        = np.in1d(f['commonEnvelopes']['randomSeed'][...].squeeze(), seedsRemain)
    
    
    seedsCEE     = f['commonEnvelopes']['randomSeed'][...].squeeze()
    DCOmaskCEE    = np.in1d(seedsCEE, seedsSystemBHNS)    
    
    maskInterest = primaryPostHG & secondaryPostHG & seeds & DCOmaskCEE
    
    seedsRemain    = f['commonEnvelopes']['randomSeed'][...].squeeze()[maskInterest]
    
    systemsOfInterest  = np.in1d(allSeedsSystems, seedsRemain)
    
    
    print("percentage double core CE channel = %s"%((np.sum(weights[systemsOfInterest])/normalisation)*100))
    del primaryPostHG
    del secondaryPostHG
    del seeds

    gc.collect()    
    
    
    
    
    
    seedsChannel = seedsRemain
    percentageChannel = ((np.sum(weights[systemsOfInterest])/normalisation)*100)
    return seedsChannel, percentageChannel







def returnSeedsPercentageSingleCoreCEE(pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/fiducial/',\
                                types='BHNS',  withinHubbleTime=True, optimistic=False, \
                                binaryFraction=1):
    
    
    #The first step is to look at the number of systems interacting at this metallicity
    #at the time of simulating RLOF and CEE printing was not robust yet
    # to catch all the systems. Some MS and touching systems might
    #not be there but should not be too many

    # retrive hdf5 file 
    # f          = h5.File(pathCOMPASOutput+'COMPASOutput.h5')
    f          = h5.File(pathCOMPASOutput)

    # obtain the BHNS seeds: 
    Data            = CC.COMPASData(path=pathCOMPASOutput, lazyData=True, Mlower=5., \
                     Mupper=150., binaryFraction=binaryFraction)
    Data.setCOMPASDCOmask(types=types,  withinHubbleTime=withinHubbleTime, optimistic=optimistic)
    Data.setCOMPASData()
    
    SeedsBHNS    = Data.seeds


    
    allSeedsSystems  = f['systems']['SEED'][...].squeeze()
    weights    = f['systems']['weight'][...].squeeze()
    maskSystemsBHNS = np.in1d(allSeedsSystems, SeedsBHNS)
    seedsSystemBHNS = allSeedsSystems[maskSystemsBHNS]


    normalisation = float(np.sum(weights[maskSystemsBHNS]))
    # print('normalisation (weighted rate DCO) = ', normalisation)
    # print('# DCOs of interest = ', len(weights[maskSystemsBHNS]))


    #which seeds are interacting (Z doesnt matter)
    RLOFseeds  = np.unique(f['RLOF']['randomSeed'][...].squeeze())
    CEEseeds   = np.unique(f['commonEnvelopes']['randomSeed'][...].squeeze())
#     print RLOFseeds.shape, CEEseeds.shape
    InteractingSeeds = np.unique(np.concatenate((RLOFseeds, CEEseeds)))


    #Which of the Seeds in SeedsSystems are interacting?
    boolInteracting  = np.in1d(seedsSystemBHNS, InteractingSeeds)
    # print 'of all BHNS seeds %s percent experiencedMT'\
    #     %((np.sum(weights[maskSystemsBHNS][boolInteracting])/normalisation)*100)
    #from here on we use the seedsOf interest to do the next step in the formation channel
    seedsInterest = allSeedsSystems[maskSystemsBHNS][boolInteracting]    
    #these arrays are long delete and clear memory
    # del maskZ
    # del allSeedsSystems
    # del CEEseeds
    # del RLOFseeds
    # del boolInteracting
    # del InteractingSeeds
    #garbage collect forces deleted arrays to clear from memory
    gc.collect()

    
    
    # During this research the RLOF printing was very much in development
    # This means that some output such as the x-th moment of RLOF was 
    # not robust (yet). Hence we use a clever trick with the seeds
    # Credits  idea to J.W.Barrett, made it more general to get not only
    # firt RLOF per system but x-th. Snipper is in coencodeVarious.cp

    #carefull the cosmicInt simulation is large and this will 
    #take a lot of RAM-memory - be patient and hope (fingerscrossed)
    #array of 90.000.000 for allSeeds RLOF possibly more for RLOF
    #close other windows :O X'(

    #UPDATE we only want of metallicity Z
    RLOFseeds   = f['RLOF']['randomSeed'][...].squeeze()
    RLOF_Z_MASK = np.in1d(RLOFseeds,seedsInterest)

    # dont care about 4th and 5th moment etc will take longer to
    # calculate those
    nthMoment   = getXmomentOfMT(RLOFseeds, maxCounter=3)
    # print(len(RLOFseeds))


    gc.collect()
    
    ############################################
    #############   MASS TRANSFER   ############
    ############################################
    

    # percentage of systems post MS primaries
    # MS transfering onto MS secondary  (channel 1 and 2)
    primaryPostHG = (f['RLOF']['type1Prev'][...].squeeze() >2)[RLOF_Z_MASK]
    secondaryMS   = ((f['RLOF']['type2Prev'][...].squeeze() == 1) |\
                                  (f['RLOF']['type2Prev'][...].squeeze() == 0))[RLOF_Z_MASK]
    
    firstMoment = (nthMoment == 1)[RLOF_Z_MASK]

    DCOmask     = np.in1d(RLOFseeds[RLOF_Z_MASK],  seedsSystemBHNS)
    
    maskInterest      = [primaryPostHG & secondaryMS & firstMoment & DCOmask]
    seedsRemain       = f['RLOF']['randomSeed'][...].squeeze()[RLOF_Z_MASK][maskInterest]
    
    
    systemsOfInterest  = np.in1d(allSeedsSystems, seedsRemain)
    nrSystemsInterest = np.sum(weights[systemsOfInterest])
    

    testMoment   = getXmomentOfMT(seedsRemain, maxCounter=10)
    # print np.max(testMoment)


    # print "post MS primar onto MS secondary "
    # print "of interacting systems %s percentage"\
    #       %((nrSystemsInterest/normalisation)*100)

    gc.collect()
    
    
    
    # primary post MS, secondary MS 
    primaryPostHG = (f['commonEnvelopes']['stellarType1'][...].squeeze() >2) &\
                    (f['commonEnvelopes']['stellarType1'][...].squeeze() <7)
    secondaryMS = f['commonEnvelopes']['stellarType2'][...].squeeze()<2
    seeds        = np.in1d(f['commonEnvelopes']['randomSeed'][...].squeeze(), seedsRemain)
    
    
    seedsCEE     = f['commonEnvelopes']['randomSeed'][...].squeeze()
    DCOmaskCEE    = np.in1d(seedsCEE, seedsSystemBHNS)    
    
    maskInterest = primaryPostHG & secondaryMS & seeds & DCOmaskCEE
    maskInterest = (primaryPostHG==1) & (secondaryMS==1) & (seeds==1) & (DCOmaskCEE==1)
    tempMask = (maskInterest==1)
    # print
    # print('---testting')
    seedsRemain    = f['commonEnvelopes']['randomSeed'][...].squeeze()[maskInterest]
    
    systemsOfInterest  = np.in1d(allSeedsSystems, seedsRemain)
    # print('old masking')
    # print(len(f['commonEnvelopes']['randomSeed'][...].squeeze()))
    # print(len(seedsRemain))
    # print('new masking')
    # # print len(f['commonEnvelopes']['randomSeed'][...].squeeze()[tempMask])
    # print(len(weights[systemsOfInterest]), '= systems of Interest')
    
    # print('end test')
    
    
    print("percentage single core CE channel = %s"%((np.sum(weights[systemsOfInterest])/normalisation)*100))
    del primaryPostHG
    del secondaryMS
    del seeds

    gc.collect()    
    
    
    
    
    
    seedsChannel = seedsRemain
    percentageChannel = ((np.sum(weights[systemsOfInterest])/normalisation)*100)
    return seedsChannel, percentageChannel








def returnSeedsPercentageOther(pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/fiducial/',\
                                types='BHNS',  withinHubbleTime=True, optimistic=False, \
                                binaryFraction=1, channelsSeedsList=[]):
    
    
    # f          = h5.File(pathCOMPASOutput+'COMPASOutput.h5')
    f          = h5.File(pathCOMPASOutput)

    # obtain the BHNS seeds: 
    Data            = CC.COMPASData(path=pathCOMPASOutput, lazyData=True, Mlower=5., \
                     Mupper=150., binaryFraction=binaryFraction)
    Data.setCOMPASDCOmask(types=types,  withinHubbleTime=withinHubbleTime, optimistic=optimistic)
    Data.setCOMPASData()
    

    mask_notInChannels = np.ones_like(Data.seeds)
    for seeds in channelsSeedsList:
        mask_notInChannels -= np.in1d(Data.seeds, np.array(seeds))

    maskk = (mask_notInChannels==1)

    seedsChannel = Data.seeds[maskk]
    percentageChannel = (np.sum(Data.weight[maskk])/np.sum(Data.weight))*100
    
    print('percentage other channel = ', percentageChannel)    
    return seedsChannel, percentageChannel



# seedsChannelSCCE, percentageChannelSCCE = returnSeedsPercentageOther()

# def createEmptyCSVplaceholder(DCOtype='BBH', nBPSmodels=15):


   


#     DCOname=dictDCOtypeDCOlabel[DCOtype]
#     BPSnameslist = list(string.ascii_uppercase)[0:nBPSmodels]   
#     channel_names = ['total', 'I_classic', 'II_only_stable_MT', 'III_single_core_CE', 'IV_double_core_CE', 'V_other']


#     NAMES = []
#     # stringgg = 'GW190814rate'

#     for ind_m, m_ in enumerate(BPSnameslist):
#         for ind_c, c_ in enumerate(channel_names):
#             str_ = m_ + ' ' + c_ + '  [Msun^{-1}]'

#             NAMES.append(str_)

            
            


#     datas=[]
#     nMetallicities = 53
#     Zlist=['0_0001','0_00011', '0_00012', '0_00014', '0_00016', '0_00017',\
#     '0_00019', '0_00022', '0_00024', '0_00027', '0_0003', '0_00034',\
#     '0_00037', '0_00042', '0_00047', '0_00052', '0_00058', '0_00065',\
#     '0_00073', '0_00081', '0_0009', '0_00101', '0_00113', '0_00126',\
#     '0_0014', '0_00157', '0_00175', '0_00195', '0_00218', '0_00243',\
#     '0_00272', '0_00303', '0_00339', '0_00378', '0_00422', '0_00471',\
#     '0_00526', '0_00587', '0_00655', '0_00732', '0_00817', '0_00912',\
#     '0_01018', '0_01137', '0_01269', '0_01416', '0_01581', '0_01765', '0_01971', '0_022', '0_0244', '0_02705', '0_03']

#     for i in range(len(NAMES)):
#         datas.append(np.zeros(nMetallicities))
#         # datas.append(np.zeros(nMetallicities))
        
            
#     df = pd.DataFrame(data=datas, index=NAMES, columns=Zlist).T
#     df.columns =   df.columns.map(str)
#     df.index.names = ['Z_i']
#     df.columns.names = ['model']

#         # print(df) 

#     df.to_csv('formationRatesTotalAndPerChannel_'+DCOname+ '_' +  '.csv')
#     return 


# INITIALIZE=False 

# if INITIALIZE == True:
#     createEmptyCSVplaceholder(DCOtype='BNS', nBPSmodels=15)
#     createEmptyCSVplaceholder(DCOtype='BNS', nBPSmodels=15)
#     createEmptyCSVplaceholder(DCOtype='BHNS', nBPSmodels=15)



# def writeFormationRatesAndChannelsToFile(DCOtype='BBH', \
#     pathCOMPASOutput='/Volumes/Andromeda/DATA/AllDCO_bugfix/',\
#      alphabetDirDict=[], nBPSmodels=15):
    
    
#     BPSnameslist = list(string.ascii_uppercase)[0:nBPSmodels]   
#     channel_names = ['total', 'I_classic', 'II_only_stable_MT', 'III_single_core_CE', 'IV_double_core_CE', 'V_other']
#     # temp = range(nModels+3)
#     DCOname=dictDCOtypeDCOlabel[DCOtype]
    

#     print('now at DCO type  ', DCOtype)
        
#     for ind_m, bps_model in enumerate(BPSnameslist):    
#         print()
#         print('now at model ', alphabetDirDict[bps_model])
            
#         # set always optimistic CE false, unless we are doing the optimistic variation
#         OPTIMISTIC=False
#         if bps_model=='H':
#             OPTIMISTIC=True
#             print('doing optimistic version of fiducial')
            
#         # path to datafile 
#         path = pathCOMPASOutput+alphabetDirDict[bps_model] + '/'

            
#         #But I want only within Hubble time 
#         Data            = CC.COMPASData(path=path, lazyData=True, Mlower=5., \
#                          Mupper=150, binaryFraction=1)
#         Data.setCOMPASDCOmask(types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC)
#         Data.setCOMPASData()
        
#         metallicities = Data.metallicitySystems
#         seeds    = Data.seeds[Data.Hubble==True]
#         weights = Data.weight
            
            
        
       

        

#         seedsPercentageClassic, seedsPercentageOnlyStableMT = returnSeedsPercentageClassicAndOnlyStableMT(pathCOMPASOutput=path,\
#                                         types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                         binaryFraction=1)
#         seedsClassic, percentageClassic = seedsPercentageClassic
#         seedsOnlyStableMT, percentageOnlyStableMT = seedsPercentageOnlyStableMT



#         seedsDoubleCE, percentageDoubleCE = returnSeedsPercentageDoubleCoreCEE(pathCOMPASOutput=path,\
#                                         types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                         binaryFraction=1)


#         seedsSingleCE, percentageSingleCE = returnSeedsPercentageSingleCoreCEE(pathCOMPASOutput=path,\
#                                         types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                         binaryFraction=1)



#         seedschannels = [seedsClassic, seedsOnlyStableMT, seedsSingleCE, seedsDoubleCE]


#         seedsOther, percentageOther = returnSeedsPercentageOther(pathCOMPASOutput=path,\
#                                         types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                         binaryFraction=1, channelsSeedsList=seedschannels)


#         seedschannels = [seedsClassic, seedsOnlyStableMT, seedsSingleCE, seedsDoubleCE, seedsOther]




#         dictChannelsBHNS = { 'classic':seedsClassic, \
#                             'immediate CE':seedsSingleCE,\
#                                  'stable B no CEE':seedsOnlyStableMT, \
#                              r'double-core CE':seedsDoubleCE,  \
#                                 'other':seedsOther\
#                                }



#         formationRateTotal           = np.zeros(len(Data.metallicityGrid))  
#         formationRateClassic         = np.zeros(len(Data.metallicityGrid)) 
#         formationRateOnlyStableMT    = np.zeros(len(Data.metallicityGrid)) 
#         formationRateSingleCE        = np.zeros(len(Data.metallicityGrid)) 
#         formationRateDoubleCE        = np.zeros(len(Data.metallicityGrid)) 
#         formationRateOther           = np.zeros(len(Data.metallicityGrid)) 


#         for nrZ, Z in enumerate(Data.metallicityGrid):
#             maskZ = (metallicities == Z)
#             formationRateTotal[nrZ] = np.sum(weights[maskZ]) # //floor weights


#             # mask different channels
#             InClassic       = np.in1d(seeds, seedsClassic)
#             InOnlyStableMT  = np.in1d(seeds, seedsOnlyStableMT)
#             InSingleCE      = np.in1d(seeds, seedsSingleCE)
#             InDoubleCE      = np.in1d(seeds, seedsDoubleCE)
#             InOther         = np.in1d(seeds, seedsOther)

#             maskClassic         = (metallicities == Z) & (InClassic==1)
#             maskOnlyStableMT    = (metallicities == Z) & (InOnlyStableMT==1)
#             maskSingleCE        = (metallicities == Z) & (InSingleCE==1)
#             maskDoubleCE        = (metallicities == Z) & (InDoubleCE==1)
#             maskOther           = (metallicities == Z) & (InOther==1)

#             formationRateClassic[nrZ]         = np.sum(weights[maskClassic])
#             formationRateOnlyStableMT[nrZ]    = np.sum(weights[maskOnlyStableMT])
#             formationRateSingleCE[nrZ]        = np.sum(weights[maskSingleCE]) 
#             formationRateDoubleCE[nrZ]        = np.sum(weights[maskDoubleCE])
#             formationRateOther[nrZ]           = np.sum(weights[maskOther])

#         formationRateTotal = np.divide(formationRateTotal, Data.totalMassEvolvedPerZ) + 0 #lowerY        
#         formationRateClassic = np.divide(formationRateClassic, Data.totalMassEvolvedPerZ)
#         formationRateOnlyStableMT = np.divide(formationRateOnlyStableMT, Data.totalMassEvolvedPerZ)
#         formationRateSingleCE = np.divide(formationRateSingleCE, Data.totalMassEvolvedPerZ)
#         formationRateDoubleCE = np.divide(formationRateDoubleCE, Data.totalMassEvolvedPerZ)
#         formationRateOther = np.divide(formationRateOther, Data.totalMassEvolvedPerZ)


#         df = pd.read_csv('formationRatesTotalAndPerChannel_'+DCOname+ '_' +  '.csv', index_col=0)
#         # namez0 = bps_model +' total  [Msun^{-1}]'
#         for ind_c, c_ in enumerate(channel_names):
#             str_ = bps_model + ' ' + c_ + '  [Msun^{-1}]'

#             # total rates 
#             if c_=='total':
#                 df[str_] = formationRateTotal 
#             elif c_=='I_classic':
#                 df[str_] = formationRateClassic
#             elif c_=='II_only_stable_MT':
#                 df[str_] = formationRateOnlyStableMT
#             elif c_=='III_single_core_CE':
#                 df[str_] = formationRateSingleCE
#             elif c_=='IV_double_core_CE':
#                 df[str_] = formationRateDoubleCE
#             elif c_=='V_other':
#                 df[str_] = formationRateOther


#         df.to_csv('formationRatesTotalAndPerChannel_'+DCOname+ '_' +  '.csv')


#     print('finished')

    # return
