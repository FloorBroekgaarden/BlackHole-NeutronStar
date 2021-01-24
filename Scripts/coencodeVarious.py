import numpy as np
import matplotlib.pyplot as plt
#Set latex environment for plots/labels
import matplotlib
matplotlib.rc('font', **{'family': 'sans-serif'})#, 'sans-serif': ['Helvetica']})
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

"""
Created by Coen Neijssel during master thesis :D Winter 2015
Thourough clean up/spring cleaning and OOP!!!! during 1st year PhD on 30-31 May-2017

The goal is 2D binning, useful for imshow or colourcoded histograms
The default is set to give you the delay time distribution of your dataframe (or it will break :p)
"""
class binning2D(object):

    def __init__(self, verbose=True, binningY=False, weightedHist=False, weightsData=None, \
                 normalizePerBinWidthX =True, normalizeTotalNumber = False,\
                 arrayXdata=None, minValueX=0.01, maxValueX=15e4,  nrBinsX=200, logBinningX=True,\
                 arrayYdata=None, minValueY=0.01, maxValueY=15e4,  nrBinsY=10, logBinningY=True):


        """ General settings for Binning
        """
        self.verbose   = verbose
        self.binningY  = binningY

        self.weightedHist     = weightedHist
        self.nameColumnWeights= weightsData
        
        self.normalizePerBinWidthX    = normalizePerBinWidthX
        self.normalizeTotalNumber     = normalizeTotalNumber      
    
        """ Settings  bins X-axis direction
        """
        self.arrayXdata      = arrayXdata
        self.minValueX       = minValueX
        self.maxValueX       = maxValueX
        self.nrBinsX         = nrBinsX
        self.logBinningX     = logBinningX
        self.edgesX, self.binWidthsX = self.createBins(self.logBinningX, self.minValueX,\
                                                       self.maxValueX, self.nrBinsX)

        """ Settings  bins Y-axis direction
        """
        self.arrayYdata      = arrayYdata
        self.minValueY       = minValueY
        self.maxValueY       = maxValueY
        self.nrBinsY         = nrBinsY
        self.logBinningY     = logBinningY
        self.edgesY, self.binWidthsY = self.createBins(self.logBinningY, self.minValueY,\
                                                       self.maxValueY, self.nrBinsY)

        self.verbosePrintOptions()

    def createBins(self,logBinning, minValue, maxValue, nrBins):
        if logBinning:
            minValue = np.log10(minValue)
            maxValue = np.log10(maxValue)
            edges    = np.logspace(minValue, maxValue, nrBins+1, base =10.0)
        else:
            edges    = np.linspace(minValue, maxValue, nrBins+1)

        binWidths = (edges[1:] - edges[:-1])

        return edges, binWidths









    def binningX(self, data):

        nrLines  = len(data)
        array2D = np.zeros(shape=(self.nrBinsX, nrLines), dtype=bool)
        
        """ Bin the data and change the above Boolean accoringly if the line falls in a bin
            then that bin is assigned true. Basically I am binning lines not data, this is important
            since this means that all the information in each line is maintained. 
        """
        dataX = arrayXdata
       
        for i in range(len(self.edgesX)-1): #-1 because last edge does not relate to a bin
            array2D[i] = (dataX  <= self.edgesX[i+1] ) & (dataX > self.edgesX[i])

        """"This array now looks like

                       line 1   line 2  line 3  line 4   .... len(data)
                       ---------------------------------------------
            bin  1   |  False   True    False  False              |    <---
            bin  2   |  True    False   False  False              |       |
            bin  3   |  False   False   True   False              |       |
              .      |                                ..          |       |
              .      |                                  ..        |       | 
            nr_bins_h|                                            |       |
                     -----------------------------------------------      |
                                                                          |
                                                                          |
                This sum is total value in bin1 since True reads as 1 in np.sum()


             You could either sum only the True to see the  number of systems in each bin 
            (e.g. no weights) or you could multiply each element with its weight.

        """

        return array2D


    def binningYFunc(self, array, data):
        nrLines  = len(data)
        dataY    = arrayYdata
        array3D = np.zeros(shape=(self.nrBinsX, self.nrBinsY, nrLines), dtype=bool)
        for nr, line in enumerate(array):          #for every bin
            dat = line * dataY                          #looks if the line is in that bin since line is 
                                                      #True or False. If True then value continues for
                                                      # binning. If false then value is zero and will
                                                      # fall outside of bins.
            for j in range(len(self.edgesY)-1):           #bin them again with booleans
               array3D[nr][j] = (dat <= self.edgesY[j+1] ) & (dat > self.edgesY[j])
        return array3D

    def weightArray(self, array, data):
        weightedArray = np.zeros(shape=array.shape)
        if self.weightedHist:
            if self.binningY:
                weights = weightsData
                print(np.sum(weights))
                xs, ys, zs = array.shape
                for x in range(xs):
                    for y in range(ys):
                        weightedArray[(x,y)] = np.multiply(array[(x,y)],  weights)
                array = weightedArray
            else:
                weights = weightsData
                weightedArray = np.multiply(array, weights)
                array = weightedArray
        return array



    def normalizeBins(self, array):
        if self.normalizePerBinWidthX:
            array = array/self.binWidthsX
        if self.normalizeTotalNumber:
            array = array/float(np.sum(array) )
        return array


    def collapseTransposeArray(self, array):
        if self.binningY:
            array = np.sum(array, axis =2)
        else:
            array = np.sum(array, axis =1)

        #row is now a bin so Transpose that each column is a bin
        array = array.T
        return array


    def mainFunction(self, data):

        array = self.binningX(data)

        if self.binningY:
            array = self.binningYFunc(array, data)

        array = self.weightArray(array, data)
        array = self.collapseTransposeArray(array)
        array = self.normalizeBins(array)
        return array, self.edgesX, self.edgesY

    def verbosePrintOptions(self):
        if self.verbose:
            verboseText ='''\n \n            2D-Binning created by Coen Neijssel :D
                           \r         Options for settings and their default values\n
                           \r--------------------------------------------------------------
                           \r|                                                            |
                           \r|  verbose=True                 binningY=False               |
                           \r|  weightedHist=False           nameColumnWeights='weights'  |
                           \r|  normalizePerBinWidthX =True  normalizeTotalNumber = False |
                           \r|  nameColumnBinsX='tc'         minValueX=0.01               |
                           \r|  maxValueX=15e4               nrBinsX=200                  |
                           \r|  logBinningX=True             nameColumnBinsY='tc'         |
                           \r|  minValueY=0.01               maxValueY=15e4,              |
                           \r|  nrBinsY=10                   logBinningY=True             |
                           \r|                                                            |
                           \r--------------------------------------------------------------
                        '''
            print( verboseText)  #used \r to not print indentation)







def layoutAxes(ax, nameX='', nameY='', \
               labelSizeMajor = 10, fontsize = 25, second=False, labelpad=None):
    """
    Tiny code to do the layout for axes in matplotlib
    """
    tickLengthMajor = 10
    tickLengthMinor = 5
    tickWidthMajor  = 1.5
    tickWidthMinor  = 1.5
    
    #rc('axes', linewidth=2)
    #label1 always refers to first axis not the twin 
    if not second:
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
    if second:
        for tick in ax.xaxis.get_major_ticks():
            tick.label2.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label2.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.2)
    ax.tick_params(length=tickLengthMajor, width=tickWidthMajor, which='major')
    ax.tick_params(length=tickLengthMinor, width=tickWidthMinor, which='minor')
    ax.set_xlabel(nameX, fontsize=fontsize,labelpad=labelpad)#,fontweight='bold')
    ax.set_ylabel(nameY, fontsize=fontsize,labelpad=labelpad)#, fontweight='bold')    

    return ax





def stepifyMe(array, edgesX, edgesY, ax):
    """use the binedges to make a x coordinates for step histogram
       instead of lines
    """
    doubled = np.repeat(edgesX, 2)
    mask = np.ones(len(doubled), dtype=bool)
    mask[0] = False
    mask[-1] = False
    stepX = doubled[mask]
    array = np.repeat(array, 2, axis=ax) #also need twice the y values
    return array, stepX


def scatterBinner(arrays, axisLabels, booleans,  groupLabels, colours, \
    binsX, binsY, locLegend=1, saveFile=None, everyOther=1, lineWidth=0.,\
    logSide=False, logTop=False, logX=False, logY=False, dpi=300, dotsize=1.):
    """
    
    arrays:      contains two sub arrays of equal length where
                 array[0] goes on the x-axis and array[1] on the y-axis
                 note that they must be of same length for the scatter plot
                 
    axisLabels:  List of the label on the x-axis and y-axis
    
    booleans:    contains x-number of boolean series which groups the data
    
    groupLabels: The name to appear on the legend for each boolean Group
                 
    colours:     A list of colours to colour each boolean group by
                 for now only accepts arguments as r,g,b,c,k etc
       
    minX/maxX/resX: The min max and nr of bins to bin the data on the x-axis
    minY/maxY/resY: The min max and nr of bins to bin the data on the y-axis
    
    locLegend : location legend in scatter plot (top right=1 then counter clockwise)
    """

    
    #Create a list with the x,y data in subgroups
    paramX = []
    paramY = []
    for boolean in booleans:
        paramX.append(arrays[0][boolean])
        paramY.append(arrays[1][boolean])

    #create the boundaries and shapes of the subplots using gridspec
    fig = plt.figure(figsize=(15,10))
    gs1 = gridspec.GridSpec(60,60)
    gs1.update(left=0.05, right=0.48, wspace=0.0, hspace=0.0)
    
    ax2 = fig.add_subplot((gs1[20:60, 0:40]))
    #below I use shareX and shareY so the limits of the plot are the same
    #as the scatterplot (limits plot NOT the bins)
    ax1 = fig.add_subplot(gs1[0:20, 0:40], sharex=ax2)
    ax3 = fig.add_subplot(gs1[20:60, 40:60], sharey=ax2)



    #make scatterplot
    for nrL, label, in enumerate(groupLabels):
        ax2.scatter(paramX[nrL], paramY[nrL], \
                    c=colours[nrL],  alpha=0.3, lw=0., s=dotsize)
    ax2.legend(loc=1)
    ax2.set_xlabel(axisLabels[0])
    ax2.set_ylabel(axisLabels[1])

    
    #Use the inbuild stacked barchart form matplotlib but with our bins and colours
    n, bins, patches = ax3.hist(list(paramY), binsY, normed=False, stacked=True, \
                                    orientation='horizontal', color=colours, alpha=0.75, lw=lineWidth)
        
    n, bins, patches = ax1.hist(list(paramX), binsX, normed=False, stacked=True, \
                                    orientation='vertical', color=colours, alpha=0.75, lw=lineWidth)
    #I dont want ticks on the shared axis
    ax3.get_yaxis().set_visible(False)
    ax1.get_xaxis().set_visible(False)
    #Remove the first tick (lowest) so it doesnt overlap with the scatterplot
    ax1.set_yticks(ax1.get_yticks()[1:])
    ax3.set_xticks(ax3.get_xticks()[1:])

    #Of the hitgram on the right that is flipped flip the labels too
    for tick in ax3.get_xticklabels():
        tick.set_rotation(270)

    #make custom legend which I find prettier
    recs =[]
    classes = []
    for nrG, group in enumerate(groupLabels):
        if nrG % everyOther == 0:
            recs.append(mpatches.Rectangle((0,0),1,1,fc=colours[nrG],alpha=0.6, lw=0.0))
            classes.append(group)
    ax2.legend(recs,classes,loc=locLegend,  ncol=1, fancybox=False, shadow=False, prop={'size': 10}) 

    if logY == True:
        ax3.set_yscale('log')
        ax2.set_yscale('log')
    if logX == True:
        ax2.set_xscale('log')
        ax1.set_xscale('log')
    if logSide == True:
        ax3.set_xscale('log')
    if logTop  == True:
        ax1.set_yscale('log')


    if saveFile is None:
        plt.show()
    else:    
        plt.savefig(saveFile, bbox_inches='tight', dpi=dpi)


    return




def plot2Darray(binnedArray=None, binEdgesX=None, binEdgesY=None, \
                colourMap='inferno', offsetColours=0, step=True, alpha=0.7,\
                labelsX='', labelY='', labelL='',\
                locLegend=1, nrColLegend=1, everyOther=1,\
                axes=None,typePlot='fill', boxBelow=False, bbox=[0., 0., 10, 10,10,10]):
    
    #transform array so we have step
    if typePlot != 'density':
        centerY = (binEdgesY[:-1] + binEdgesY[1:]) / 2.
        if step:
            doubled = np.repeat(binEdgesX, 2)
            mask = np.ones(len(doubled), dtype=bool)
            mask[0] = False
            mask[-1] = False
            stepX = doubled[mask]
            binnedArray = np.repeat(binnedArray, 2, axis=1)
            
        else:
            center = (binEdgesX[:-1] + binEdgesX[1:]) / 2.
            stepX = center
    #Create custom colormap to index over
    cm = plt.get_cmap(colourMap)
    nColors = len(binnedArray)+offsetColours
    mycolors = [cm(x) for x in np.linspace(0,1 , nColors)] 
    
    
    if typePlot=='fill':
        lowerFill = np.zeros(len(binnedArray[0]))
        for nrR, row in enumerate(binnedArray):
            upperFill = lowerFill + row
            axes.fill_between(stepX, lowerFill, upperFill, color=mycolors[nrR], alpha=alpha, lw=0.)
            lowerFill = upperFill
    if typePlot=='lines':
        for nrR, row in enumerate(binnedArray):
            axes.plot(stepX, row, color=mycolors[nrR], alpha=alpha)
        pass
        
    #create legend:
    if typePlot=='lines' or typePlot=='fill':
        #make custom legend which I find prettier
        recs =[]
        classes = []
        counter = 0
        for nrR, row in enumerate(binnedArray):
            if np.sum(row)==0:
                pass
            else:
                if counter % everyOther ==0:
                    recs.append(mpatches.Rectangle((0,0),1,1,fc=mycolors[nrR],alpha=alpha, lw=0.0))
                    classes.append(np.round(binEdgesY[nrR], 2))
                counter+=1
    if not boxBelow:
        axes.legend(recs,classes,loc=locLegend,  ncol=nrColLegend, fancybox=False, shadow=False, prop={'size': 10},\
                bbox_to_anchor=(1.3, 1.), title=labelL)
    else:
        legend = axes.legend(recs,classes,bbox_to_anchor=(bbox[0], bbox[1]),ncol=bbox[2], title=labelL,\
                              prop={'size': bbox[3]} , fontsize=bbox[4])
        legend.get_title().set_fontsize(bbox[5])
    if typePlot=='density':
        pass
    axes.set_xlim(binEdgesX[0], binEdgesX[-1])
    return axes



def returnNormedPerBinPanel(arrays, binsX, binsY,\
                           NormTotal=False, NormEachBin=False, weights=None): 
    
    
    #First create the subgroups that we are going to stack
    booleans  = np.zeros(shape=(binsY[2], len(arrays[1])), dtype=bool)
    binEdgesY = np.linspace(binsY[0], binsY[1], binsY[2]+1)
    
    for nrY, _, in enumerate(binEdgesY[:-1]):
        boolean = (arrays[1] >= binEdgesY[nrY]) & (arrays[1] < binEdgesY[nrY+1])
        booleans[nrY] = boolean
        
    #Now create the 2D array which will have nrBins in X as columns
    #and each row is a subgroup from above binned in the x direction
    binnedArray = np.zeros(shape=(binsY[2], binsX[2]))
    binEdgesX = np.linspace(binsX[0], binsX[1], binsX[2]+1)
    
    for nrGroup, boolGroup in enumerate(booleans):
        dataGroup = arrays[0][boolGroup]
        for nrX, _, in enumerate(binEdgesX[:-1]):
            boolean = (dataGroup >= binEdgesX[nrX]) & (dataGroup < binEdgesX[nrX+1])
            if weights is None:
                binnedArray[nrGroup][nrX] =  np.sum(boolean)
            else:
                weightsGroup = weights[boolGroup]
                weightedBool = numpy.multiply(boolean,weightsGroup)
                binnedArray[nrGroup][nrX] =  np.sum(weightedBool)
            
    #Now we have a two-D array with the value in each bin(X,Y)
    #This can be directly plotted for imshow for a density map
    #Or we can normalize it in several manners.
    
    if NormTotal:
        total = float(np.sum(binnedArray))
        binnedArray = binnedArray/total
    #if NormEachGroup:
    #    pass
    if NormEachBin:
        sumPerBin = np.sum(binnedArray, axis=0).astype(float)
        mask = sumPerBin == 0
        sumPerBin[mask]=1 #if there is nothing in it we might as well divide by 1 to prevent error
        binnedArray = np.divide(binnedArray, sumPerBin)
    
    return binnedArray, binEdgesX, binEdgesY





def scatterPanel(fig, axes,  binnedArray, binEdgesX, binEdgesY, \
                 resolution=0, labels=[None, None], colourmap='jet', log=True, maskColour='k'):
    
    tinyX = (max(binEdgesX) - min(binEdgesX))/float(len(binEdgesX))
    tinyY = (max(binEdgesY) - min(binEdgesY))/float(len(binEdgesY))

    if log:
        mask = (binnedArray !=0)
        binnedArray[mask] = np.log10(binnedArray[mask])
    else:
        mask = binnedArray != 0

    centerX = (binEdgesX[1:] + binEdgesX[:-1])/2.
    centerY = (binEdgesY[1:] + binEdgesY[:-1])/2.
        
    #plotting
    vmin = np.min(binnedArray[mask])
    vmax = np.max(binnedArray[mask])
    for nrY, Y in enumerate(centerY):
        Y = np.ones(len(centerX))*Y
        #masked ones get colour
        
        maskInvert = np.logical_not(mask[nrY])
        axes.scatter(centerX[maskInvert], Y[maskInvert], c=maskColour, \
                     alpha=1., s=30, lw=0., marker=",",\
                     vmin=vmin, vmax=vmax, cmap=colourmap)
        sc = axes.scatter(centerX[mask[nrY]], Y[mask[nrY]], c=binnedArray[nrY][mask[nrY]], \
                     alpha=0.9, s=35, lw=0.,marker=",",\
                     vmin=vmin, vmax=vmax, cmap=colourmap)
    axes.set_xlim(min(centerX)-tinyX, max(centerX)+tinyX)
    axes.set_ylim(min(centerY)-tinyY, max(centerY)+tinyY)
    axes.set_xlabel(labels[0], fontsize=20)
    axes.set_ylabel(labels[1], fontsize=20)
    fig.colorbar(sc, ax=axes )
    return 

def oneD_stepBinning(axes, array, label, resolution=100):
    
    bins = np.linspace(min(array),max(array),resolution)
    values, edges = np.histogram(array, bins)
    hist = values/float(np.max(values))
    
    cumulative = np.cumsum(values/float(np.sum(values)))
    center = (edges[1:] + edges[:-1])/2.
    
    doubled = np.repeat(edges, 2)
    mask = np.ones(len(doubled), dtype=bool)
    mask[0] = False
    mask[-1] = False
    stepX = doubled[mask]
    hist = np.repeat(hist, 2, axis=0)
    
    axes.plot(stepX, hist, c='k', lw=2.)
    axes.plot(center, cumulative, c='k', linestyle=':', lw=2.)
    axes.set_xlabel(label, fontsize=20)
    return
    

def trianglePlot(fig, axes, Parameters, ParLabels, resolution=100, colourmap='jet'):

    for i in range(len(Parameters)):
        for j in range(len(Parameters)):
            if i > j:
                arrays = [Parameters[j], Parameters[i]]
                labels = [ParLabels[j] , ParLabels[i]]
                scatterPanel(fig, axes[i][j], arrays, resolution, labels=labels, colourmap=colourmap)
            elif i == j:
                oneD_stepBinning(axes[i][j], Parameters[i],  ParLabels[i],  resolution)
            else:
                fig.delaxes(axes[i][j])
            for tick in axes[i][j].xaxis.get_major_ticks():
                tick.label.set_fontsize(20)
            for tick in axes[i][j].yaxis.get_major_ticks():
                tick.label.set_fontsize(20)


    
def getXmomentOfMT(Seeds, maxCounter=10):
    #this function might become obsolete if we finetune the RLOF output
    #with help from idea Jim Barrett
    #to have a number for x-moment of RLOF

    #make seeds into 1D array and calculate difference, meaning everytime the next line
    #has same seed it will be zero else it will be more, except for the very first line.
    offsetIndices = np.diff(Seeds)
    #I dont care about the difference just that it is nonzero, make it all into 0-1s
    offsetIndices[offsetIndices>=1] = 1

    #Create am empty array to turn into a boolean slice, since the np.diff ommits
    #first line we add one to the length.
    indices       = np.zeros(len(offsetIndices)+1)
    indices[0]    = 1
    #Now an array with 0 and 1s where every one is the first line of a different seed.
    #This effectively is the first moment of mass transfer of the system.
    indices[1:]   = offsetIndices

    
    
    #so nr 1 is first moment,
    counter = 2
    while (0. in indices) and (counter <=maxCounter):
        #get indices
        indexFilled   = np.where(indices != 0)
        #add 1 essentially move one row down
        indexFilledTemp  = indexFilled[0] + 1
        #if not marked alreaydy i.e. in indexFilled 
        notMarked = np.logical_not(np.in1d(indexFilledTemp,indexFilled))
        #and if index not bigger than array
        notTooBig = indexFilledTemp < (len(indices) -1)
        #give me those indices
        indexFilledTemp = indexFilledTemp[notMarked & notTooBig]
        #and fill in the RLOF counter as anotehr moment of RLOF
        indices[np.array(indexFilledTemp,)] = counter
        counter+=1
    return indices


def drawConcentricFigure():
    #Having some fun drawing
    minRedshift   = 0
    maxRedshift   = 2.
    nRedshiftBins = 10.
    redshiftEdges = np.linspace(0,maxRedshift,nRedshiftBins+1) #The bin edges in redshift
    redshifts = 0.5*(redshiftEdges[:-1] + redshiftEdges[1:])      #Central value of each redshift bin
    angles = np.linspace(0,2*np.pi, 1e4)

    fig, axes = plt.subplots(1,1, figsize=(15,3.75))

    for nrBin, zBin in enumerate(redshiftEdges[:-1]):
        #Draw circles of at every redshift edge
        if nrBin == 0:
            label='shell Edges'
        else:
            label = None
        x = redshiftEdges[nrBin] * np.cos(angles)
        y = redshiftEdges[nrBin] * np.sin(angles)
        axes.plot(x,y,c='k', lw=2, linestyle=':', label=label)
        
        #for one shell draw a line thick enough
        #so it looks like the shell filled
        if nrBin == 5:
            x = redshifts[nrBin] * np.cos(angles)
            y = redshifts[nrBin] * np.sin(angles)
            axes.plot(x,y,c='b', lw=35, alpha=0.2, label='volume Shell')
        
    #redshift point at center/halfway point shell which we use for luminosity distance
    axes.scatter(redshifts, np.zeros(len(redshifts)), label='redshift at center shell')

    #draw arrow from detector to redshift for luminosity Distance
    axes.arrow(0, 0.05, redshifts[5]-0.05, 0.0, head_width=0.04, fc='k')
    axes.annotate('luminosity distance to shell', xy=(0.1, 0.1), xytext=(0.1, 0.1), fontsize=15)
    #Extra text
    axes.arrow(1.9, 0.3, 0., -0.22, head_width=0.04, fc='k')
    axes.annotate('calculate mergers at \n each of these points', xy=(1.9, 0.3), xytext=(1.9, 0.3), fontsize=15)
    #Central point to show location detector
    axes.scatter(0,0, c='r', label='detector', s=100)

    axes.get_yaxis().set_visible(False)
    axes.legend(loc=2, prop={'size':18})
    axes.set_ylim(-0.5, 0.5)
    nameX = r'$\rm redshift \ z$'
    nameY = r'$ $'
    axes= layoutAxes(axes, nameX=nameX, nameY=nameY)

    plt.show()


#SFR prescription for illustration
def SFR_Madau(z): #[Msun yr-1 Gpc-3]
    return 0.015* ((1+z)**2.7) / ( 1 + ((1+z)/2.9)**5.6) * 1e9 

def drawingConcept(x, y):
    fig,axes = plt.subplots(1,1, figsize=(30,9))
    redshifts = np.linspace(0,3,100)
    axes.plot(redshifts, SFR_Madau(redshifts)/float(1e8), c='k', linestyle=':', label='SFR Madau et al.', lw=2.)
    axes.set_ylabel('SFR*1e8 [Msun yr-1 Gpc-3]')
    axes.set_xlabel('redshift z')
    #Extra text
    axes.arrow(0.0+x, 0.3, 0.5, 0.0, head_width=0.05, fc='k')
    axes.arrow(0.1+x, 0.6, -0.07, -0.23, head_width=0.02, fc=None)
    axes.arrow(0.63+x, 0.35+y, -0.05, 0.1, head_width=0.02, fc=None)
    axes.annotate('delay time', xy=(0.1, 0.31), xytext=(0.1+x, 0.31), fontsize=25)
    axes.annotate('binary merges here', xy=(0.1, 0.6), xytext=(0.1+x, 0.6), fontsize=25)
    axes.annotate('binary born here with this SFR', xy=(0.63, 0.35), xytext=(0.63+x, 0.35+y), fontsize=25)
    nameX = r'$\rm redshift \ z$'
    nameY = r'$\rm SFR*1e8 [Msun yr-1 Gpc-3]$'
    axes= layoutAxes(axes, nameX=nameX, nameY=nameY)
    axes.legend(loc=2, prop={'size':25})
    plt.show()
