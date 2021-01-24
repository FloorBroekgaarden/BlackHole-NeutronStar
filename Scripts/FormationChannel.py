import numpy as np
import h5py as h5 #for reading in data
from collections import Counter
from IPython.display import clear_output

def returnUniqueChannelsCounts(deleteColumns=[], pathToH5=None, \
                               chunkSize=100000,startFile=0, endFile=None,
                              booleanFilter = None, weights=None):
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
    f = h5.File(pathToH5+'COMPASOutput.h5')
    form = f['formationChannels']
    
    #Reduce the keys to only the columns we want to 
    #always ignore seed since this is always different
    columns      = form.keys()
    i = columns.index("m_randomSeed")  
    del columns[i]
    #Now delete other columns to ignore
    deleteColumns = deleteColumns
    for string in deleteColumns:
        i = columns.index(string)  
        del columns[i]    
  

    #define chunks and start end
    chunkSize   = chunkSize
    startFile   = startFile
    startChunk  = startFile
    endChunk    = None


    #If end not supplied go over entire length file
    if endFile == None:
        #use any column to find the length
        endFile = len(f['formationChannels'][columns[0]][...])
    else:
        endFile = endFile

    #create header jsut in case it goes through the file in a weird order
    header         = []
    for nrH, head in enumerate(columns):
        header.append(head)
        
    #overallCounter will be the dictionary tracking
    #keys   = channel
    #values = count
    overallCounter = Counter()
    sumtuples = 0
    print 'Counting/sorting the unique Channels'
    while endChunk!= endFile:
        # clear_output(wait=True)
        #Test that the last chunk ends at the proper index
        endChunk = startChunk + chunkSize
        print(str(np.round(endChunk/float(endFile), 3)*100)+'%')

        if endChunk > endFile:
            endChunk  = endFile        
            chunkSize = endChunk - startChunk
        #temporary array to cast columns in
        #every row is a column of length chunk
        if booleanFilter is None:
            tempArray = np.zeros(shape=(len(columns), chunkSize))
        else:
            tempArray = np.zeros(shape=(len(columns), np.sum(booleanFilter[startChunk:endChunk])))
        for nrC, column in enumerate(columns):
            #Cant slice beyond the end
            if booleanFilter is None:
                tempArray[nrC] = form[column][startChunk:endChunk].squeeze()
            else:
                tempArray[nrC]  = form[column][startChunk:endChunk]\
                                 [booleanFilter[startChunk:endChunk]].squeeze()
        #now the columns are the rows so to turn it into columns take the transpose
        tempArray = tempArray.T
        
        # \floor add weights:
        if weights is None:
            if booleanFilter is None:
                tempweights = np.zeros_like(chunkSize)
            else:
                tempweights = np.zeros_like(np.sum(booleanFilter[startChunk:endChunk]))
        else: # if weights are not all 1
            if booleanFilter is None:
                tempweights = weights[startChunk:endChunk]
            else:
                tempweights = weights[startChunk:endChunk]\
                [booleanFilter[startChunk:endChunk]]


        # cast them to tuples (which are hashable)        
        tuples = [(tuple(row),tempweights[ind]) for ind,row in enumerate(tempArray)]
        #count the things in this chunk give them as keys and values
        #and add them to the overall dictionary
        # print Counter(tuples), 'counter(tuples)'
        # print 'tuples' ,len(tuples)



        weightedcounts = {}
        for valuetemp in tuples:
            if valuetemp[0] in weightedcounts:
                weightedcounts[valuetemp[0]] += valuetemp[1]
            else:
                weightedcounts[valuetemp[0]] = valuetemp[1]
        # print Counter(weightedcounts)



        sumtuples += len(tuples)
        overallCounter += Counter(weightedcounts)
        startChunk = endChunk
        
    if booleanFilter is None :
        if np.sum(overallCounter.values()) != (endFile - startFile):
            raise ValueError('nr of channels counted does not correspond to nr channels looked at')
    else:
        if (np.sum(overallCounter.values()) != (np.sum(booleanFilter))) & (weights is None):
            raise ValueError('nr of channels counted does not correspond to nr channels looked at')
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
    print sumtuples, 'sumtuples!'
    return header, sortedChannels, sortedCounts



def returnArrayRanks(deleteColumns=[], pathToH5=None, chunkSize=100000,startFile=0,\
                     endFile=None, booleanFilter=None, saveArray=False, weights = None):

    header, sortedChannels, sortedCounts = returnUniqueChannelsCounts(deleteColumns=deleteColumns, \
                            pathToH5=pathToH5, chunkSize=chunkSize,startFile=startFile, endFile=endFile,\
                                                                     booleanFilter=booleanFilter, weights=weights)

    #Read in the file
    f = h5.File(pathToH5+'COMPASOutput.h5')
    fForm = f['formationChannels']
    #TODO update this function so it does it by chunks
    #Bottleneck is reading in the file by boolean slicing
    sortedTuples = [tuple(x) for x in sortedChannels]
    sortedRankes = range(len(sortedTuples))
    newDict = dict(zip(sortedTuples, sortedRankes))

    print 'Creating column with rank of each system, patience required :) '
    tempArray = np.zeros(shape=(len(header), np.sum(booleanFilter)))
    for nrC, column in enumerate(header):
        tempArray[nrC] = fForm[column][...][booleanFilter].squeeze()
        #print nrC
    
    tuples = [tuple(x) for x in tempArray.T]
    rankArray = np.zeros(len(tuples))
    for nrR, row in enumerate(tuples):
        rankArray[nrR] = newDict[row]

    print 'Done :D '
    return header, sortedChannels, sortedCounts, rankArray

def formationToString(header, row, rank, count):
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


def returnChannelsGivenMask(pathToH5, maskDCO, maxRank=6):

    f = h5.File(pathToH5+'COMPASOutput.h5')
    seeds     = f['doubleCompactObjects']['seed'][maskDCO]

    #Only want formataion channels of BBH seeds

    #Syntax -> f['formationChannels']['m_randomSeed'][...].squeeze()
    # of file formation channels -> of the column seeds -> read in all -> squeeze in 1D array
    bools = np.in1d(f['formationChannels']['m_randomSeed'][...].squeeze(),seeds)


    header, sortedChannels, sortedCounts, rankArray =\
    returnArrayRanks(deleteColumns=[], pathToH5=pathToH5, chunkSize=100000,\
                        startFile=0, endFile=None, booleanFilter=bools)

    #I only care about the top 5 channels. Will make all above 5 a single number (6) so
    #It is easy to mask
    mask = rankArray >=maxRank 
    rankArray[mask] = maxRank
    printer = np.round(np.sum(mask)/float(len(mask)), 3)*100
    print 'maxRank=%s, percentage systems in remainder = %s' %(maxRank, printer)
    return header, sortedChannels, sortedCounts, rankArray


def printChannels(header, sortedChannels, sortedCounts, \
                  rankArray, totalSystems, maxRank=1000):
    accounted = 0.
    for channel in range(maxRank-1):
        count      = sortedCounts[channel]
        row        = sortedChannels[channel]
        rank       = channel + 1
        percentage = count/float(totalSystems)*100
        accounted += percentage
        print '\n ----channel  '+str(rank)+' ---%='+str(percentage)+'----  count='+str(count)
        print(formationToString(header, row, rank, count))

    print '\n channels', maxRank-1, ' account for %', accounted ,  'of all binaries evolved'
