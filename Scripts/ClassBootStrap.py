import numpy as np
import h5py as h5
import os.path
import sys 





#Need to do the sampling on the same level as were the montecarlo is done
#which is at the all initials stage PER METALLICITY

#Prefer not to store seeds of all initials in the class because they are large


#I want to directly compare to the number of mergers of interest that I had originally
#i.e. from initial sample to drawn sample of interest after COMPAS evolution

#Class is a bit of overkill but who knows for future stuff

class Bootstrap(object):


    def __init__(self, verbose = False, pathToFile = None, \
                        fileName = 'BootstrapSample.h5', nrBootstraps=5):


        self.nrBootstraps       = nrBootstraps
        self.pathToFile         = pathToFile
        self.fileName           = fileName

        self.verbose            = verbose



    def PrintVerbose(self, message):
        if self.verbose:
            print(message)

    def sampleFromAll(self, allSeeds, allMetallicities):
        self.PrintVerbose('start seeds sampling')

        seedsSample= np.zeros(len(allSeeds), dtype='int64')
        for nrZ, Z in enumerate(np.unique(allMetallicities)):
            boolZ  = (allMetallicities == Z)
            seedsZ = allSeeds[boolZ]
            seedsZresampled = np.random.choice(seedsZ, size=len(seedsZ),\
                                                  replace=True, p=None)
            seedsSample[boolZ] = seedsZresampled
        self.PrintVerbose('end seeds sampling')
        return seedsSample


    def rateMultiplier(self, seedsInterest, seedsSample):
        """ I want to associate which seeds I drew to seeds Interest
            difficulty is in the duplicate draws
        """
        self.PrintVerbose('multiplier start calculation')
        
        #count duplicates and sort seeds
        seedsSample, counts = np.unique(seedsSample, return_counts=True)
        #Only care about the seeds that make BBHs
        boolReduceCounts= np.in1d(seedsSample, seedsInterest)
        seedsSample     = seedsSample[boolReduceCounts]
        counts          = counts[boolReduceCounts]
        self.PrintVerbose('multiplier reduced sample to look at interest')

        #Indices on how to sort the BBH seeds
        sorter          = np.argsort(seedsInterest)
        InterestSorted  = seedsInterest[sorter]

        #Reduce to only have seeds we sample
        #so we can directly correlate with counts
        boolReduceMergers = np.in1d(InterestSorted, seedsSample)
        self.PrintVerbose('multiplier bool reduce seeds interest drawn')

        multiplyRates = np.zeros(len(seedsInterest))
        #The code below in words is:
        #If we would sort the rates of the BBH we sampled the counts are= counts
        multiplyRates[sorter[boolReduceMergers]] = counts
        #This way multiplyrates has the 
        #original order so we can directly associate 
        #with other arrays of seeds Interest like rates

        self.PrintVerbose('end multiplier calculation')
        return multiplyRates

    def createAndWriteToFile(self, allSeeds=None, allMetallicities=None, seedsInterest=None):

        #Does the file already exist?
        if os.path.isfile(self.pathToFile+self.fileName):
            print("File exist so I am done Already")
            bootstrapped = h5.File(self.pathToFile+self.fileName)
            self.nrBootstraps = len(bootstrapped.keys())
            print("sample size is %s" %(self.nrBootstraps))
        else:
            print("File does not exits depending on nrBootstraps might take some time")
            h5f = h5.File(self.pathToFile+self.fileName, 'a')
            counter = 0
            while counter < self.nrBootstraps:
                seedsSample   = self.sampleFromAll(allSeeds, allMetallicities)
                multiplyRates = self.rateMultiplier(seedsInterest, seedsSample)
                h5f.create_dataset(str(int(counter)), data=multiplyRates)
                counter+=1
                print "\r " + str(counter)+ '/' +str( self.nrBootstraps),
                sys.stdout.flush()
            h5f.close()

    def ReadSample(self, nrSample=None, maskDCO=None):

        bootstrapped = h5.File(self.pathToFile+self.fileName)
        return bootstrapped[str(int(nrSample))][...][maskDCO.squeeze()]


