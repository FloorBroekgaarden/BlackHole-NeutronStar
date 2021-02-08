import numpy as np
import h5py as h5 #for reading in data
from collections import Counter
import os
class FormationChannels(object):

    def __init__(self, path = 'setPath', deleteColumns=[], \
                              chunkSize=100000,startFile=0, endFile=None, verbose=False):

        #general settings data
        self.path   = path
        if os.path.isfile(path):
            self.h5file           = h5.File(path)
        else:
            raise ValueError("h5 file not found. Invalid path?")

        #settings for gathering the formation channels
        #see function returnUniqueChannelsCounts
        self.deleteColumns = deleteColumns
        self.chunkSize     = chunkSize 
        self.startFile     = startFile 
        self.endFile       = endFile


        self.verbose       = verbose
        self.seedsInterest = None 


        #Results 
        self.header         = None   #which part of formation channel header is used
        self.sortedChannels = None   #array with the sorted unique channels from formation channel
        self.sortedCounts   = None   #array with the counts of sorted unique channels
        self.sortedStrings  = []     #list with the string-format of the sorted unique channels
        self.sortedSeeds    = []

        self.weights = None # weights of samples  //floor 



    def returnUniqueChannelsCounts(self): # //floor added weights 
        """
        Sort the formation channels based on the columns of interest.
        Do this in chunks in case of big file sizes or small memory.
        
        Input:
        
        deleteColumns: List of strings taken from file.keys()
                       These, columns will not be taken into account when comparing
                       Formation Channels.
        pathToH5:      string to the H5 output directory
        chuckSize:     integer, how many rows to read in at once, the bigger the more
                       memory demanding
        start:         integer, index row from which to start slicing the data
        end:           integer, index row where to stop reading the file.
        booleanFilter: booleanArray, length data set formation channels, to subselect
                       systems of interest for example only the ones in the merger file
        
        Returns:
        
        header:         List of strings that returns the header of the uniqueChannels
                        this might not be same order as the header of the file.
        
        uniqueChannels: 2D numpy.array where every row is a formation channel
                        Note that it only contains the columns looked at.
                        This array is sorted where the most occuring channel 
                        is the first row and the least occuring channel the last.
        counts:         1D array same length as the number of rows of uniqueChannels
                        Each element is the number of counts corresponding to
                        the same row of uniqueChannels.    
        """
        
        #Read in the file
        form = self.h5file['formationChannels']
        
        #Reduce the keys to only the columns we want to 
        #always ignore seed since this is always different
        columns      = form.keys()
        i = columns.index("m_randomSeed")  
        del columns[i]
        #Now delete other columns to ignore
        deleteColumns = self.deleteColumns
        for string in deleteColumns:
            i = columns.index(string)  
            del columns[i]    
      

        #define chunks and start end
        chunkSize   = self.chunkSize
        startFile   = self.startFile
        startChunk  = self.startFile
        endChunk    = None


        #If end not supplied go over entire length file
        if self.endFile == None:
            #use any column to find the length
            endFile = len(form[columns[0]][...])
        else:
            endFile = self.endFile

        #create header jsut in case it goes through the file in a weird order
        header         = []
        for nrH, head in enumerate(columns):
            header.append(head)
        if self.verbose:
            print 'self.header \n Sorted which columns from formation channels we use'
        #overallCounter will be the dictionary tracking
        #keys   = channel
        #values = count
        overallCounter = Counter()
        if self.verbose:
            print 'Counting/sorting the unique Channels'
        while endChunk!= endFile:
            #Test that the last chunk ends at the proper index
            endChunk = startChunk + chunkSize

            if endChunk > endFile:
                endChunk  = endFile        
                chunkSize = endChunk - startChunk
            #temporary array to cast columns in
            #every row is a column of length chunk
            if self.booleanFilter is None:
                tempArray = np.zeros(shape=(len(columns), chunkSize))
            else:
                tempArray = np.zeros(shape=(len(columns), \
                                     np.sum(self.booleanFilter[startChunk:endChunk])))
            for nrC, column in enumerate(columns):
                #Cant slice beyond the end
                if self.booleanFilter is None:
                    tempArray[nrC] = form[column][startChunk:endChunk].squeeze()
                else:
                    tempArray[nrC]  = form[column][startChunk:endChunk]\
                                     [self.booleanFilter[startChunk:endChunk]].squeeze()
            #now the columns are the rows so to turn it into columns take the transpose
            tempArray = tempArray.T
            # cast them to tuples (which are hashable)


            # // floor added weights
            self.weights = self.h5file['systems']['weight'][...].squeeze()

            if self.booleanFilter is None:
                tempweights = self.weights[startChunk:endChunk]
            else:
                tempweights = self.weights[startChunk:endChunk]\
                [self.booleanFilter[startChunk:endChunk]]
            # // floor added weights
            
            # // floor added weights
###########
        # cast them to tuples (which are hashable)        
            tuples = [(tuple(row),tempweights[ind]) for ind,row in enumerate(tempArray)]
################


            # \\floor
            weightedcounts = {}
            for valuetemp in tuples:
                if valuetemp[0] in weightedcounts:
                    weightedcounts[valuetemp[0]] += valuetemp[1]
                else:
                    weightedcounts[valuetemp[0]] = valuetemp[1]
            # print Counter(weightedcounts)
            # \\floor
            overallCounter += Counter(weightedcounts) # // floor
            startChunk = endChunk #//floor 

        # // floor ## below 
        # if self.booleanFilter is None:
        #     if np.sum(overallCounter.values()) != (endFile - startFile):
        #         raise ValueError('nr of channels counted does not correspond to nr channels looked at')
        # else:
        #     if np.sum(overallCounter.values()) != (np.sum(self.booleanFilter)):
        #         raise ValueError('nr of channels counted does not correspond to nr channels looked at')
        # // floor above



        if self.verbose:
            print 'self.sortedChannels self.sortedCounts \n Combining the counts and channels in arrays'
        channels = overallCounter.keys()
        counts   = overallCounter.values()
        zippedChannelsCounts = zip(channels, counts)
        # in the zipped list, of everything in there sort by the second element of that thing    
        sortedChannelsCounts = sorted(zippedChannelsCounts, key=lambda keyValue: keyValue[1])

        sortedChannels, sortedCounts = zip(*sortedChannelsCounts)
     
        #the channels are a list of tuples I like arrays so
        sortedChannels = np.array([list(row) for row in sortedChannels])
        
        #lowest count is sorted first I want the highest count so inverse both
        sortedChannels = sortedChannels[::-1]
        sortedCounts   = sortedCounts[::-1]

        self.header         = header
        self.sortedChannels = sortedChannels
        self.sortedCounts   = sortedCounts

        print(np.sum(sortedCounts), 'sum sorted counts')

        if self.verbose:
            print 'Done'
        return



    def returnArrayRanks(self):

        #Read in the file
        fForm = self.h5file['formationChannels']
        #TODO update this function so it does it by chunks
        #Bottleneck is reading in the file by boolean slicing
        sortedTuples = [tuple(x) for x in self.sortedChannels]
        sortedRankes = range(len(sortedTuples))
        newDict = dict(zip(sortedTuples, sortedRankes))
        if self.verbose:
            print 'self.rankArray \n Creating column with rank of each system, patience required :)'
        tempArray = np.zeros(shape=(len(self.header), np.sum(self.booleanFilter)))
        for nrC, column in enumerate(self.header):
            tempArray[nrC] = fForm[column][...][self.booleanFilter].squeeze()
            #print nrC
        
        tuples = [tuple(x) for x in tempArray.T]
        rankArray = np.zeros(len(tuples))
        for nrR, row in enumerate(tuples):
            rankArray[nrR] = newDict[row]
        if self.verbose:
            print 'Done :D '

        self.rankArray = rankArray
        return

    def formationToString(self, header, row, rank, count):
        """
        This function translates the row into a human readable string
        given the specific header that was used to create the row.
        
        The code below is a bit elaborate but I hope
        it is clearer to see what is going on and easier to
        adapt if the header changes.
        """
        dictIdx = dict(zip(header, range(len(header))))
        dictRow = dict(zip(header, row))
        

        
        #The possible headers of which the value indicates time-ordering
        timeHeader = ['mt_primary_ep1','mt_primary_ep2', 'mt_primary_ep3',\
                     'mt_secondary_ep1','mt_secondary_ep2', 'mt_secondary_ep3',\
                     'SN_primary_type_1', 'SN_primary_type_2', 'SN_primary_type_3',\
                     'SN_secondary_type_1','SN_secondary_type_2', 'SN_secondary_type_3',\
                     'CEE','CEE_failed', 'CEE_wet']
        #Ok look up the columns we have used and the time counters
        headerTimes = []
        valueTimes = []
        for key in dictRow:
            if key in timeHeader:
                value = int(dictRow[key])
                #Only want info if something happened so nonzero
                if value != 0:
                    headerTimes.extend([key])
                    valueTimes.extend([value])
        
        #Put the header in order given the times in the values
        args = np.argsort(valueTimes)
        headerTimes = np.array(headerTimes)[args]
        
        #each of the header that have time info might have additional
        #info such as stellar types or wether it is the primary or secondary
        #that started the phase
        #This dictionary links the above list of timeHeader
        #to the name of the possible additional columns
        dictInfo = {'mt_primary_ep1':('mt_primary_ep1_K1','mt_primary_ep1_K2'),\
                    'mt_primary_ep2':('mt_primary_ep2_K1','mt_primary_ep2_K2'),\
                    'mt_primary_ep3':('mt_primary_ep3_K1','mt_primary_ep3_K2'),\
                    'mt_secondary_ep1':('mt_secondary_ep1_K1', 'mt_secondary_ep1_K2'),\
                    'mt_secondary_ep2':('mt_secondary_ep2_K1', 'mt_secondary_ep2_K2'),\
                    'mt_secondary_ep3':('mt_secondary_ep3_K1', 'mt_secondary_ep3_K2'),\
                    'CEE': ('CEE_instigator'),\
                    'CEE_failed':('CEE_failed_instigator'),\
                    'CEE_wet':('CEE_wet_instigator')}

        #look in our headerTimes and see what we have
        #then look if we asked for additional info and append that column name in place
        headerWithInfo = []
        for string in headerTimes:
            headerWithInfo.extend([string])
            try:
                addInfo = dictInfo[string]
                for addString in addInfo:
                    if addString in header:
                        headerWithInfo.extend([addString])
            except:
                pass
        
        #We also have some generic information which if present we just append to end
        genericInfo = ['stellar_type_K1','stellar_type_K2','merged_in_Hubble_time','binary_disbound']
        for generic in genericInfo:
            if generic in header:
                headerWithInfo.extend([generic])
        
        #ok so by now the headerWithInfo contains the a ordered list of strings (headernames)
        #which we just need to map to the index of the value and replace with the string we want.
        dictTypes = {0:'MS<0 ',1:'MS ', 2:'HG ', 3:'FGB ', 4:'CHeB ', 5:'EAGB ', 6:'TPAGB ',\
             7:'HeMS ',8:'HeHG ',9:'HeGB ', 10:'HeWD ', 11:'COWD ', 12:'ONeWD ',\
             13:'NS ', 14:'BH ', 15:'MR '}
        dictPS = {1:' Primary', 2:' Secondary'}
        dictBool = {0:'False', 1:'True'}
        dummy = ['']
        stellarTypeInfo = {'mt_primary_ep1_K1':(' P=', dictTypes),\
                           'mt_primary_ep1_K2':(' S=', dictTypes),\
                           'mt_primary_ep2_K1':(' P=', dictTypes),\
                           'mt_primary_ep2_K2':(' S=', dictTypes),\
                           'mt_primary_ep3_K1':(' P=', dictTypes),\
                           'mt_primary_ep3_K2':(' S=', dictTypes),\
                           'mt_secondary_ep1_K1':(' P=', dictTypes),\
                           'mt_secondary_ep1_K2':(' S=', dictTypes),\
                           'mt_secondary_ep2_K1':(' P=', dictTypes),\
                           'mt_secondary_ep2_K2':(' S=', dictTypes),\
                           'mt_secondary_ep3_K1':(' P=', dictTypes),\
                           'mt_secondary_ep3_K2':(' S=', dictTypes),\
                           'CEE_instigator':(' started by', dictPS),\
                           'CEE_failed_instigator':(' started by', dictPS),\
                           'CEE_wet_instigator':(' started by', dictPS)}
        
        superNovaInfo = {  'SN_primary_type_1':'ccSN primary ',\
                           'SN_primary_type_2':'ECSN primary ',\
                           'SN_primary_type_3':'USSN primary ',\
                           'SN_secondary_type_1':'ccSN secondary ',\
                           'SN_secondary_type_2':'ECSN secondary ',\
                           'SN_secondary_type_3':'USSN secondary '}
        
        additionalInfo = {'stellar_type_K1':(' last type Primary in Form=', dictTypes),\
                          'stellar_type_K2':(' last type Secondary in Form=', dictTypes),\
                          'merged_in_Hubble_time':(' merged within Hubble=', dictBool),\
                          'binary_disbound':(' binary disrupted=', dictBool)
                         }
        
        outputString =''
        for name in headerWithInfo:
            if name in stellarTypeInfo.keys():
                string, func = stellarTypeInfo[name]
                idxRow = int(dictIdx[name])
                value = int(row[idxRow])
                outputString += string + func[value]
            elif name in superNovaInfo.keys():
                outputString += ' ->' + superNovaInfo[name]
            elif name in additionalInfo.keys():
                string, func = additionalInfo[name]
                idxRow = int(dictIdx[name])
                value = int(row[idxRow])
                outputString +=' ->'+ string + func[value]
            else:
                outputString += ' ->' + name
        return outputString




    def seedsOfInterestDCO(self, types='BBH', withinHubbleTime=True, optimistic=False):
        #We do not want all the formation channels just the ones that form BBHs
        fDCO    = self.h5file['doubleCompactObjects']
        if types == 'BBH':
            maskTypes = (fDCO['stellarType1'][...].squeeze() == 14) &\
                        (fDCO['stellarType2'][...].squeeze() == 14)
        elif types == 'BNS':
            maskTypes = (fDCO['stellarType1'][...].squeeze() == 13) &\
                        (fDCO['stellarType2'][...].squeeze() == 13)
        elif types == 'BHNS':
            maskTypes = ((fDCO['stellarType1'][...].squeeze() == 14) &\
                        (fDCO['stellarType2'][...].squeeze() == 13)) |\
                        ((fDCO['stellarType1'][...].squeeze() == 13) &\
                        (fDCO['stellarType2'][...].squeeze() == 14))
                
        if withinHubbleTime == True:
            maskHubble = (fDCO['mergesInHubbleTimeFlag'][...].squeeze()==True)
        else:
            #Array where all are true
            maskHubble = np.ones(len(fDCO['mergesInHubbleTimeFlag'][...].squeeze()), dtype=bool)
                          
        if optimistic == True:
            #we do not care about the optimistic flag (both False and True allowed)
            #Array where all are true
            maskOptimistic = np.ones(len(fDCO['optimisticCEFlag'][...].squeeze()), dtype=bool)
        else:
            #optimistic scenario not allowed (pessimistic) hence the flag must be false
            #This removes systems with CEE from HG donors (no core envelope separation)
            maskOptimistic = fDCO['optimisticCEFlag'][...].squeeze() == False
                          
        #we never want in first timestep after CEE, because 
        #we define it as a system that should not have survived the CEE
        maskNoRLOFafterCEE =  (fDCO['RLOFSecondaryAfterCEE'][...].squeeze()==False)
                          

        maskDCO = maskTypes & maskHubble & maskOptimistic & maskNoRLOFafterCEE
        #what are the seeds of the systems that form the DCO?
        seedsOfInterest = fDCO['seed'][...].squeeze()[maskDCO]
        return seedsOfInterest

    def groupSeeds(self):
        if self.verbose:
            print 'self.sortedSeeds \n Creating array per channel with all seeds of that channel'
        fForm   = self.h5file['formationChannels']
        allSeeds= fForm['m_randomSeed'][...].squeeze()[self.booleanFilter]

        for nr in np.unique(self.rankArray):
            maskSeed = (self.rankArray == nr)
            self.sortedSeeds.append(allSeeds[maskSeed])
        return

    def formationChannelsSeeds(self, seeds=None, types=None, withinHubbleTime=True, optimistic=False):
        #We only want the formation cahnnels of the seeds above
        fForm   = self.h5file['formationChannels']
        if self.verbose:
            print 'self.booleanFilter \n Looking up the seeds/systems and creating mask for formation Channels. '
        if (seeds is None) & (types is None):
            raise ValueError("Need either seeds (array) or types (string) input")
        if (seeds is not None) & (types is not None):
            raise ValueError("Given both seeds and types , what do you want None?")
        elif types is None:
            self.booleanFilter = np.in1d(fForm['m_randomSeed'][...].squeeze(),seeds)
        elif (seeds is None) and (types in ['BBH', 'BHNS', 'BNS']):
            seeds = self.seedsOfInterestDCO(types=types, withinHubbleTime=withinHubbleTime,\
                                             optimistic=optimistic)
            self.booleanFilter = np.in1d(fForm['m_randomSeed'][...].squeeze(),seeds)
        else:
             raise ValueError("cannot create boolean filter of these types yet+"\
                            "currently only for 'BBH', 'BHNS', 'BNS'")
        if np.sum(self.booleanFilter) == 0:
            raise ValueError("These seeds do not exist in formation channel output")

        if self.verbose:
            print "Doing the formation channels for %s systems" %(np.sum(self.booleanFilter))
        self.returnUniqueChannelsCounts()
        self.returnArrayRanks()
        self.groupSeeds()
        if self.verbose:
            print 'self.sortedStrings \n Constructing human readable string for each of the unique Channels (Magic) '
        for nrC in np.unique(self.rankArray.astype(int)):
            count      = self.sortedCounts[nrC]
            row        = self.sortedChannels[nrC]
            stringChannel = self.formationToString(self.header, row, nrC, count)
            self.sortedStrings.append(self.formationToString(self.header, row, nrC, count))

        
        if self.verbose:
            print " :D :D \n finished in total we have %s channels for %s systems"\
                    %(len(self.sortedChannels), np.sum(self.sortedCounts))
        # if np.sum(self.sortedCounts) != np.sum(self.booleanFilter):  #//floor weights
        #     raise ValueError("number of counts %s ,  nr seeds %s, sum BooleanFilter %s,"\  #//floor weights
        #                      %(np.sum(self.sortedCounts), len(seeds), np.sum(self.booleanFilter))\  #//floor weights
        #                      + " something went wrong")
        return 

    def resetResults(self):
        self.header         = None   #which part of formation channel header is used
        self.sortedChannels = None   #array with the sorted unique channels from formation channel
        self.sortedCounts   = None   #array with the counts of sorted unique channels
        self.sortedStrings  = []     #list with the string-format of the sorted unique channels
        self.sortedSeeds    = []
        return

    def returnSeedsNotOrContainingString(self, test=None, string=None):
    
        seedsReturn = []

        for nr in range(len(self.sortedStrings)):
            channelString =  self.sortedStrings[nr]
            if test == 'include':
                if string in channelString:
                    seedsReturn.extend(self.sortedSeeds[nr])
            elif test == 'exclude':
                if string not in channelString:
                    seedsReturn.extend(self.sortedSeeds[nr])
            else:
                raise ValueError("test argument is either 'include' or 'exclude'")

        return seedsReturn

