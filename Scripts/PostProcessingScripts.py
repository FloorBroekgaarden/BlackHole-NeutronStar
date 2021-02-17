from __future__ import division # in case this script is used in python 2 
import h5py as h5


import numpy as np
import string



import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec


import matplotlib.gridspec as gridspec
# import matplotlib
import matplotlib.pyplot as plt
# for e.g., minor ticks 
from matplotlib.ticker import (FormatStrFormatter,
                               AutoMinorLocator)
#Set latex environment for plots/labels
import matplotlib
matplotlib.rc('font', **{'family': 'sans-serif'})#, 'sans-serif': ['Helvetica']})
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']



from matplotlib.offsetbox import AnchoredText
from matplotlib import rc                                                                                                                                                                                                                    
from matplotlib import rcParams
import seaborn as sns

from astropy import units as u
from astropy import constants as const


from scipy.spatial.distance import cdist



rc('font', family='serif', weight = 'bold')
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
rc('axes', linewidth=2)

matplotlib.rcParams['xtick.major.size'] = 12
matplotlib.rcParams['ytick.major.size'] = 12
matplotlib.rcParams['xtick.minor.size'] = 8
matplotlib.rcParams['ytick.minor.size'] = 8
matplotlib.rcParams['font.weight']= 'bold'
matplotlib.rcParams.update({'font.weight': 'bold'})

fs = 24 # fontsize for plots
rc('axes', linewidth=2)


def layoutAxes(ax, nameX='', nameY='', \
               labelSizeMajor = 10, fontsize = 25, second=False, labelpad=None, setMinor=True):
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
    
    if setMinor==True:
        # add minor ticks:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

    return ax



bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.75) # for box around text in plot




        
dictChannelsBHNSListBolt = [r'\textbf{case A}', r'\textbf{case B}',  r'\textbf{case C}',\
                            r'\textbf{case B only stable}',\
                            r'\textbf{case B immediate CE}',
                            r'\textbf{case C immediate} \textbf{CE}',\
                             r'\textbf{double-core CE}', r'\textbf{other}']




zorderlist = { 'stable B':10, 'stable B no CEE':13, \
                    'case B immediate CE':12,'stable C':15,\
                    r'case C immediate CE':17, 
                    'stable A':14, \
                 r'double-core CE':11, 'other':16\
                 }








dictChannelsBHNSListBolt = [r'\textbf{(I) Classic}', \
                            r'\textbf{(II) Only stable mass transfer}',\
                            r'\textbf{(III) Single-core CE as first mass transfer}',\
                             r'\textbf{(IV) Double-core CE as first mass transfer}', r'\textbf{(V) Other}']




dictChannelsBHNSList = ['classic', \
                      'stable B no CEE', \
                    'immediate CE',\
                 r'double-core CE', 'other']




zorderlist = { 'classic':10, 'stable B no CEE':13, \
                    'immediate CE':12,\
                 r'double-core CE':11, 'other':16\
                 }


# default settings for labels and names of BPS models 
nModels=15
BPSnameslist = list(string.ascii_uppercase)[0:nModels]
modelDirList = ['fiducial', 'massTransferEfficiencyFixed_0_25', 'massTransferEfficiencyFixed_0_5', 'massTransferEfficiencyFixed_0_75', \
               'unstableCaseBB', 'alpha0_5', 'alpha2_0', 'fiducial', 'rapid', 'maxNSmass2_0', 'maxNSmass3_0', 'noPISN',  'ccSNkick_100km_s', 'ccSNkick_30km_s', 'noBHkick' ]

alphabetDirDict =  {BPSnameslist[i]: modelDirList[i] for i in range(len(BPSnameslist))}


physicalNamesBPSmodels = [r'\textbf{fiducial}',\
                           r'$\beta=0.25$', r'$\beta=0.5$',  r'$\beta=0.75$',r'\textbf{unstable case BB}',\
                           r'$\alpha_{\rm{CE}}=0.5$',  r'$\alpha_{\rm{CE}}=2$', r'\textbf{optimistic CE}',\
                          r'\textbf{rapid SN}', r'$\rm{max} \ m_{\rm{NS}}=2.0\,\rm{M}_{\odot}$', r'$\rm{max} \ m_{\rm{NS}}=3.0\,\rm{M}_{\odot}$',\
                          r'\textbf{no PISN}', r'\textbf{SN} '+ r'$\sigma_{\rm{rms}}^{\rm{1D}}=100\,\rm{km}\,\rm{s}^{-1}$',r'\textbf{SN} '+ r'$\sigma_{\rm{rms}}^{\rm{1D}}=30\,\rm{km}\,\rm{s}^{-1}$',\
                          r'\textbf{SN} '+ r'$v_{\rm{k,BH}}=0\,\rm{km}\,\rm{s}^{-1}$' ]



alphabetPhysicalNameDict =  {BPSnameslist[i]: physicalNamesBPSmodels[i] for i in range(len(BPSnameslist))}



physicalNamesBPSmodelsWithEnter = [r'\textbf{fiducial}',\
                           r'$\beta=0.25$', r'$\beta=0.5$',  r'$\beta=0.75$',r'\textbf{unstable}' + '\n'+ r'\textbf{case BB}',\
                           r'$\alpha_{\rm{CE}}=0.5$',  r'$\alpha_{\rm{CE}}=2$', r'\textbf{optimistic CE}',\
                          r'\textbf{rapid SN}', r'$\rm{max} \ m_{\rm{NS}}$' +'\n' + r'$2.0\,\rm{M}_{\odot}$', r'$\rm{max} \ m_{\rm{NS}}$' +'\n' + r'$3.0\,\rm{M}_{\odot}$',\
                          r'\textbf{no PISN}', r'\textbf{SN} '+ r'$\sigma_{\rm{rms}}^{\rm{1D}}$' +'\n' + r'$100\,\rm{km}\,\rm{s}^{-1}$',r'\textbf{SN} '+ r'$\sigma_{\rm{rms}}^{\rm{1D}}$' +'\n' + r'$30\,\rm{km}\,\rm{s}^{-1}$',\
                          r'\textbf{SN} '+ r'$v_{\rm{k,BH}}$' +'\n' + r'$0\,\rm{km}\,\rm{s}^{-1}$' ]

alphabetPhysicalNameDictWithEnter =  {BPSnameslist[i]: physicalNamesBPSmodelsWithEnter[i] for i in range(len(BPSnameslist))}



colorlist = [ '#118AB2', '#EF476F', '#FFD166', '#073B4C', 'gray']

def obtainDataSTROOPWAFEL(param, pathToDirectory):
    """returns for STROOPWAFEL (AIS) simulation the data of wanted variable
    combines the data from AIS_oratory and AIS_sampling 
    
    param = [xparam, fxparam] ,  are the name of the variable and hdf5 keyname where it is in
    e.g. param = ['M1', 'doubleCompactObjects'] (see also: print(list(f.keys())))
    pathToDirectory is pathname to Directory where AIS_oratory & AIS_sampling directories are
    """ 

    xparam, fxparam = param

    pathAIS = pathToDirectory +'/COMPASOutput.h5'    

    fAIS = h5.File(pathAIS)
      
        
    ##### get parameter from two directories and combine them ############
    xvalues         = fAIS[fxparam][xparam][...].squeeze()
    return   xvalues



def maskTargetDCOsSTROOPWAFEL(DCOtype, boolDCOmask, f, otherSelection, otherparam):
    """returns mask of DCOs of interest
    fxparam  is hdf5 keyname of file where variable for which you want to mask DCOs is in 
    DCOtype = 'BBH' / 'ALL' / 'BHNS' or 'BNS' 
    boolDCOmask = [Hubble, RLOF, Pessimistic] # boolean values whether to mask mergers in a HUbble time, 
    binaries that have RLOFSecondaryAfterCEE = True, and Pessimistic binaries (i.e. optimisticCEFlag == 0)
    pathToDirectory is pathname to Directory where _oratory & _sampling directories are
    """
    
    Hubble, RLOF, Pessimistic = boolDCOmask
    

 
    
    fDCO = f['doubleCompactObjects']
    
    
    
    
    # mask binaries of given DCO type
    if DCOtype == 'BNS':
        mask0 = ((fDCO['stellarType1'][...] == 13) & (fDCO['stellarType2'][...] == 13))
    elif (DCOtype == 'BHNS') | (DCOtype == 'NSBH'):
        mask0 = ((fDCO['stellarType1'][...] == 13) & (fDCO['stellarType2'][...] == 14)) | \
            ((fDCO['stellarType1'][...] == 14) & (fDCO['stellarType2'][...] == 13) )          
    elif DCOtype == 'BBH':
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
    # Pessimistic mask :  if True mask systems that have optimistic CE flag ==1
    if Pessimistic:
        mask3 = np.logical_not(fDCO["optimisticCEFlag"][...] == 1)
    elif not Pessimistic:
        mask3 = np.logical_not(fDCO["optimisticCEFlag"][...] == 1) + \
        np.logical_not(fDCO["optimisticCEFlag"][...] == 0)   
    
    # combine the different masks and the oratory and refinement masks
    combinedmask = mask0 * mask1 * mask2 * mask3
    combinedmask = combinedmask.squeeze()
    if otherSelection =='UFD':
        KpcToKM = 3.086 * 10**(16) # kpc to km  
        MyrToYr = 1E6 # yrs
        YrToSec = 3.154 *1E7 #sec        
        UFD_epsilon = otherparam[0]
        UFD_Rvir = otherparam[1]
        Xbh1 = otherparam[2]
        Rns = otherparam[3]

        fSN = f['supernovae']
        seedsOfIntererst = fDCO['seed'][...].squeeze()
        seedsSN = fSN['randomSeed'][...].squeeze()
        bools = np.in1d(seedsSN, seedsOfIntererst)        
        
        tc  = fDCO['tc'][...].squeeze()
        vsys = fSN['systemicVelocity'][...].squeeze()[bools]
        vsysSN2 = vsys[1:][::2]
        traveldistance = tc * vsysSN2 *  MyrToYr * YrToSec
        radiusUFDgalaxy = UFD_epsilon * UFD_Rvir * KpcToKM
        maskCandidatesUFD = (traveldistance <= radiusUFDgalaxy) | ((vsysSN2 <= 44) & (tc * MyrToYr *YrToSec<= radiusUFDgalaxy)) 
        
        combinedmask = maskCandidatesUFD*combinedmask
    

    
    return combinedmask






def obtainweightsSTROOPWAFEL(pathToDirectory):
    """returns weights for all DCOs and all systems for STROOPWAFEL
    pathToDirectory is pathname to Directory where AIS_oratory & AIS_sampling directories are 
    """
    
    pathAIS = pathToDirectory +'/COMPASOutput.h5'   # '/home/floor/Data_Thesis/bdMC/Z0_002'

    fAIS = h5.File(pathAIS)


    ##### get the DCO and all system weights  ############
    DCOsweights          = fAIS['doubleCompactObjects']['weight'][...].squeeze()


    
    systemsweights          = fAIS['systems']['weight'][...].squeeze()

    
    
    return DCOsweights, systemsweights


def chirpmass(m1, m2):
    numer = (m1*m2)**(3./5)
    denom = (m1+m2)**(1./5)
    
    return numer/denom




def obtainM1BHandM2BHassymetric(m1, m2):
    m1bh, m2bh = np.zeros_like(m1), np.zeros_like(m1)
    maskm1heavier = ( m1 >= m2)
    maskm2heavier = (m1 < m2)
    
    m1bh[maskm1heavier] = m1[maskm1heavier] 
    m1bh[maskm2heavier] = m2[maskm2heavier]
    m2bh[maskm1heavier] = m2[maskm1heavier]
    m2bh[maskm2heavier] = m1[maskm2heavier]
    
    return m1bh, m2bh # m1bh has all the heaviest systems









def getMaskBHNS(m1bh, m2bh):
    # add later on the 2nd explodes first 
    
    maskBHNS = m1bh >= m2bh # we have a BH=NS

    
    return maskBHNS



def below3Msun(m1bh):
    # add later on the 2nd explodes first 
    
    maskBHNS = m1bh <= 3 # we have a BH=NS

    
    return maskBHNS







class gaussian_kde(object):
    """Representation of a kernel-density estimate using Gaussian kernels.

    Kernel density estimation is a way to estimate the probability density
    function (PDF) of a random variable in a non-parametric way.
    `gaussian_kde` works for both uni-variate and multi-variate data.   It
    includes automatic bandwidth determination.  The estimation works best for
    a unimodal distribution; bimodal or multi-modal distributions tend to be
    oversmoothed.

    Parameters
    ----------
    dataset : array_like
        Datapoints to estimate from. In case of univariate data this is a 1-D
        array, otherwise a 2-D array with shape (# of dims, # of data).
    bw_method : str, scalar or callable, optional
        The method used to calculate the estimator bandwidth.  This can be
        'scott', 'silverman', a scalar constant or a callable.  If a scalar,
        this will be used directly as `kde.factor`.  If a callable, it should
        take a `gaussian_kde` instance as only parameter and return a scalar.
        If None (default), 'scott' is used.  See Notes for more details.
    weights : array_like, shape (n, ), optional, default: None
        An array of weights, of the same shape as `x`.  Each value in `x`
        only contributes its associated weight towards the bin count
        (instead of 1).

    Attributes
    ----------
    dataset : ndarray
        The dataset with which `gaussian_kde` was initialized.
    d : int
        Number of dimensions.
    n : int
        Number of datapoints.
    neff : float
        Effective sample size using Kish's approximation.
    factor : float
        The bandwidth factor, obtained from `kde.covariance_factor`, with which
        the covariance matrix is multiplied.
    covariance : ndarray
        The covariance matrix of `dataset`, scaled by the calculated bandwidth
        (`kde.factor`).
    inv_cov : ndarray
        The inverse of `covariance`.

    Methods
    -------
    kde.evaluate(points) : ndarray
        Evaluate the estimated pdf on a provided set of points.
    kde(points) : ndarray
        Same as kde.evaluate(points)
    kde.pdf(points) : ndarray
        Alias for ``kde.evaluate(points)``.
    kde.set_bandwidth(bw_method='scott') : None
        Computes the bandwidth, i.e. the coefficient that multiplies the data
        covariance matrix to obtain the kernel covariance matrix.
        .. versionadded:: 0.11.0
    kde.covariance_factor : float
        Computes the coefficient (`kde.factor`) that multiplies the data
        covariance matrix to obtain the kernel covariance matrix.
        The default is `scotts_factor`.  A subclass can overwrite this method
        to provide a different method, or set it through a call to
        `kde.set_bandwidth`.

    Notes
    -----
    Bandwidth selection strongly influences the estimate obtained from the KDE
    (much more so than the actual shape of the kernel).  Bandwidth selection
    can be done by a "rule of thumb", by cross-validation, by "plug-in
    methods" or by other means; see [3]_, [4]_ for reviews.  `gaussian_kde`
    uses a rule of thumb, the default is Scott's Rule.

    Scott's Rule [1]_, implemented as `scotts_factor`, is::

        n**(-1./(d+4)),

    with ``n`` the number of data points and ``d`` the number of dimensions.
    Silverman's Rule [2]_, implemented as `silverman_factor`, is::

        (n * (d + 2) / 4.)**(-1. / (d + 4)).

    Good general descriptions of kernel density estimation can be found in [1]_
    and [2]_, the mathematics for this multi-dimensional implementation can be
    found in [1]_.

    References
    ----------
    .. [1] D.W. Scott, "Multivariate Density Estimation: Theory, Practice, and
           Visualization", John Wiley & Sons, New York, Chicester, 1992.
    .. [2] B.W. Silverman, "Density Estimation for Statistics and Data
           Analysis", Vol. 26, Monographs on Statistics and Applied Probability,
           Chapman and Hall, London, 1986.
    .. [3] B.A. Turlach, "Bandwidth Selection in Kernel Density Estimation: A
           Review", CORE and Institut de Statistique, Vol. 19, pp. 1-33, 1993.
    .. [4] D.M. Bashtannyk and R.J. Hyndman, "Bandwidth selection for kernel
           conditional density estimation", Computational Statistics & Data
           Analysis, Vol. 36, pp. 279-298, 2001.

    Examples
    --------
    Generate some random two-dimensional data:

    >>> from scipy import stats
    >>> def measure(n):
    >>>     "Measurement model, return two coupled measurements."
    >>>     m1 = np.random.normal(size=n)
    >>>     m2 = np.random.normal(scale=0.5, size=n)
    >>>     return m1+m2, m1-m2

    >>> m1, m2 = measure(2000)
    >>> xmin = m1.min()
    >>> xmax = m1.max()
    >>> ymin = m2.min()
    >>> ymax = m2.max()

    Perform a kernel density estimate on the data:

    >>> X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    >>> positions = np.vstack([X.ravel(), Y.ravel()])
    >>> values = np.vstack([m1, m2])
    >>> kernel = stats.gaussian_kde(values)
    >>> Z = np.reshape(kernel(positions).T, X.shape)

    Plot the results:

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
    ...           extent=[xmin, xmax, ymin, ymax])
    >>> ax.plot(m1, m2, 'k.', markersize=2)
    >>> ax.set_xlim([xmin, xmax])
    >>> ax.set_ylim([ymin, ymax])
    >>> plt.show()

    """
    def __init__(self, dataset, bw_method=None, weights=None):
        self.dataset = np.atleast_2d(dataset)
        if not self.dataset.size > 1:
            raise ValueError("`dataset` input should have multiple elements.")
        self.d, self.n = self.dataset.shape
            
        if weights is not None:
            self.weights = weights / np.sum(weights)
        else:
            self.weights = np.ones(self.n) / self.n
            
        # Compute the effective sample size 
        # http://surveyanalysis.org/wiki/Design_Effects_and_Effective_Sample_Size#Kish.27s_approximate_formula_for_computing_effective_sample_size
        self.neff = 1.0 / np.sum(self.weights ** 2)

        self.set_bandwidth(bw_method=bw_method)

    def evaluate(self, points):
        """Evaluate the estimated pdf on a set of points.

        Parameters
        ----------
        points : (# of dimensions, # of points)-array
            Alternatively, a (# of dimensions,) vector can be passed in and
            treated as a single point.

        Returns
        -------
        values : (# of points,)-array
            The values at each point.

        Raises
        ------
        ValueError : if the dimensionality of the input points is different than
                     the dimensionality of the KDE.

        """
        points = np.atleast_2d(points)

        d, m = points.shape
        if d != self.d:
            if d == 1 and m == self.d:
                # points was passed in as a row vector
                points = np.reshape(points, (self.d, 1))
                m = 1
            else:
                msg = "points have dimension %s, dataset has dimension %s" % (d,
                    self.d)
                raise ValueError(msg)

        # compute the normalised residuals
        chi2 = cdist(points.T, self.dataset.T, 'mahalanobis', VI=self.inv_cov) ** 2
        # compute the pdf
        result = np.sum(np.exp(-.5 * chi2) * self.weights, axis=1) / self._norm_factor

        return result

    __call__ = evaluate

    def scotts_factor(self):
        return np.power(self.neff, -1./(self.d+4))

    def silverman_factor(self):
        return np.power(self.neff*(self.d+2.0)/4.0, -1./(self.d+4))

    #  Default method to calculate bandwidth, can be overwritten by subclass
    covariance_factor = scotts_factor

    def set_bandwidth(self, bw_method=None):
        """Compute the estimator bandwidth with given method.

        The new bandwidth calculated after a call to `set_bandwidth` is used
        for subsequent evaluations of the estimated density.

        Parameters
        ----------
        bw_method : str, scalar or callable, optional
            The method used to calculate the estimator bandwidth.  This can be
            'scott', 'silverman', a scalar constant or a callable.  If a
            scalar, this will be used directly as `kde.factor`.  If a callable,
            it should take a `gaussian_kde` instance as only parameter and
            return a scalar.  If None (default), nothing happens; the current
            `kde.covariance_factor` method is kept.

        Notes
        -----
        .. versionadded:: 0.11

        Examples
        --------
        >>> x1 = np.array([-7, -5, 1, 4, 5.])
        >>> kde = stats.gaussian_kde(x1)
        >>> xs = np.linspace(-10, 10, num=50)
        >>> y1 = kde(xs)
        >>> kde.set_bandwidth(bw_method='silverman')
        >>> y2 = kde(xs)
        >>> kde.set_bandwidth(bw_method=kde.factor / 3.)
        >>> y3 = kde(xs)

        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111)
        >>> ax.plot(x1, np.ones(x1.shape) / (4. * x1.size), 'bo',
        ...         label='Data points (rescaled)')
        >>> ax.plot(xs, y1, label='Scott (default)')
        >>> ax.plot(xs, y2, label='Silverman')
        >>> ax.plot(xs, y3, label='Const (1/3 * Silverman)')
        >>> ax.legend()
        >>> plt.show()

        """
        if bw_method is None:
            pass
        elif bw_method == 'scott':
            self.covariance_factor = self.scotts_factor
        elif bw_method == 'silverman':
            self.covariance_factor = self.silverman_factor
        elif np.isscalar(bw_method): # and not isinstance(bw_method, string_types):
            self._bw_method = 'use constant'
            self.covariance_factor = lambda: bw_method
        elif callable(bw_method):
            self._bw_method = bw_method
            self.covariance_factor = lambda: self._bw_method(self)
        else:
            msg = "`bw_method` should be 'scott', 'silverman', a scalar " \
                  "or a callable."
            raise ValueError(msg)

        self._compute_covariance()

    def _compute_covariance(self):
        """Computes the covariance matrix for each Gaussian kernel using
        covariance_factor().
        """
        self.factor = self.covariance_factor()
        # Cache covariance and inverse covariance of the data
        if not hasattr(self, '_data_inv_cov'):
            # Compute the mean and residuals
            _mean = np.sum(self.weights * self.dataset, axis=1)
            _residual = (self.dataset - _mean[:, None])
            # Compute the biased covariance
            self._data_covariance = np.atleast_2d(np.dot(_residual * self.weights, _residual.T))
            # Correct for bias (http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_covariance)
            self._data_covariance /= (1 - np.sum(self.weights ** 2))
            self._data_inv_cov = np.linalg.inv(self._data_covariance)

        self.covariance = self._data_covariance * self.factor**2
        self.inv_cov = self._data_inv_cov / self.factor**2
        self._norm_factor = np.sqrt(np.linalg.det(2*np.pi*self.covariance)) #* self.n

        
        

def lowess(x, y, f=2. / 3., iter=3):
    """lowess(x, y, f=2./3., iter=3) -> yest
    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.
    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.
    """
    n = len(x)
    r = int(ceil(f * n))
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
    w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
    w = (1 - w ** 3) ** 3
    yest = np.zeros(n)
    delta = np.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:, i]
            b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
            A = np.array([[np.sum(weights), np.sum(weights * x)],
                          [np.sum(weights * x), np.sum(weights * x * x)]])
            beta = linalg.solve(A, b)
            yest[i] = beta[0] + beta[1] * x[i]

        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta ** 2) ** 2

    return yest




# EM ejecta functions:




def calculateRisco(m_bhtemp, Xefftemp):
    # this is prograde orbit
    # see also https://duetosymmetry.com/tool/kerr-isco-calculator/

    # everything in cgs
    c = 2.99792458E10 #[cm s-1] 
    G = 6.67259E-8   
    Msun = 1.99E33 # gr
    Rsun = 6.96E10 # cm     
    
    factorFront =   ((G*m_bhtemp)/c**2) #m_bhtemp #s
    
    Z1 = 1 + (1 - Xefftemp**2)**(1/3) * ((1 + Xefftemp)**(1/3) + (1 - Xefftemp)**(1/3) )
    Z2 = np.sqrt((3* Xefftemp**2 + Z1**2))
    
    Risco = factorFront * (3 + Z2 - np.sqrt((3-Z1)*(3+Z1 +2*Z2)))
    return Risco



def calculateEjectedMassMerger(m_ns, r_ns, m_bh, Xeff ):
 # from 1807.00011, Eq 4 
    # returns M_rem in solar masses 
    # input r and m in solar masses and R sun. Xeff in [0,1] (symmetric) 
    # RNS in km
    
    
    # everything in cgs
    c = 2.99792458E10 #[cm s-1] 
    G = 6.67259E-8   
    Msun = 1.99E33 # gr
    Rsun = 6.96E10 # cm         
    
    
    # convert to cgs
    r_ns  = r_ns*0.1*10**6 #np.asarray([1.E6]* len(m_ns)) # to cm
    m_ns_cgs = Msun * m_ns
    m_bh_cgs = Msun * m_bh
    
    
    alpha, beta, gamma, delta = 0.406, 0.139, 0.255, 1.761
    C_NS = G * m_ns_cgs / (r_ns * c**2)
    
    R_isco = calculateRisco(m_bh_cgs, Xeff)
    
    R_isco_norm  = R_isco / (m_bh_cgs * (G/c**2)) 
    
    Q = m_bh_cgs / m_ns_cgs
    
    eta = Q / (1 + Q)**2
    
    FirstTerm  = alpha*(1 - 2*C_NS) / eta**(1/3)
    SecondTerm = beta* R_isco_norm * C_NS / eta 
    
    A = np.asarray(FirstTerm - SecondTerm + gamma)
    B = np.zeros_like(m_ns_cgs)
    
    Mrem_model = np.maximum(A,B)**(delta)
    
    Mrem_model /= Msun # in solar masses 
    
    # and the true M remnant mass (not normalized and in solar masses =)
    Mrem_solar = Mrem_model * m_ns_cgs  
    return Mrem_solar # in [Msun]






def convert_a_to_P_circular(separation, M1, M2):
    """calculate Period from separation
    separation is separation (needs to be given in astropy units)
    M1 and M2 are masses of the binary
    
    """
    G = const.G # [gr cm s^2]
    

    mu = G*(M1+M2)
    period = 2*np.pi * np.sqrt(separation**3/mu)
    
    
    
    return period        


Arrays_minNSmassEjecta_labels = [r'$(R_{\rm{NS}},\chi_{\rm{BH}})=11.5,0$',\
                                 r'$(R_{\rm{NS}},\chi_{\rm{BH}})=13,0$',\
                                 r'$(R_{\rm{NS}},\chi_{\rm{BH}})=11.5,0.5$', \
                                 r'$(R_{\rm{NS}},\chi_{\rm{BH}})=13,0.5$']







# SFRD options : 
GSMFs = ['Panter et al. (2004) Single', 'Furlong et al. (2015) Single', 'Furlong et al. (2015) Double']
MZs   = [ 'Langer et al. (2006)'      , 'Langer et al. +offset (2006)', 'Ma et al. (2015)']
SFRs  = ['Madau et al. (2014)'         ,'Strolger et al. (2004)',     'Madau et al. (2017)']




# define list with names of the MSSFR variations :-) 
MSSFRnameslist = []
MSSFRnameslist.append('000') # add phenomenological 

for ind_GSMF, GSMF in enumerate(GSMFs):
    ind_y = ind_GSMF + 1
    for ind_MZ, MZ in enumerate(MZs):
        ind_z = ind_MZ +1
        for ind_SFR, SFR in enumerate(SFRs):
            ind_x = ind_SFR+1
            
            
            
            
        

            MSSFRnameslist.append('%s%s%s'%(ind_x, ind_y, ind_z))
# print(modelnameslistA)

DCOTypeList = ['BHBH', 'BHNS', 'NSNS']

NumberBPSmodels=15
alphabet = list(string.ascii_uppercase)
BPSnameslist = alphabet[:NumberBPSmodels]


# GW stuff


NSNSrate0 = [320-240,320+490] # Gpc-3 yr-1 from: https://arxiv.org/pdf/2001.01761.pdf
BHBHrate0 = [23.9-8.6,23.9+14.9] # Gpc-3 yr-1 from: https://arxiv.org/pdf/1811.12907.pdf
BHNSrate0 = [0,610] # Gpc-3 yr-1 from: https://arxiv.org/pdf/1811.12907.pdf 



NSNSrate0 = [250,2810] # Gpc-3 yr-1 from: https://arxiv.org/pdf/2001.01761.pdf
BHBHrate0 = [9.7,101] # Gpc-3 yr-1 from: https://arxiv.org/pdf/1811.12907.pdf
BHNSrate0 = [0,610] # Gpc-3 yr-1 from: https://arxiv.org/pdf/1811.12907.pdf 


from math import log10, floor
def round_to_1(x):
    """ round to one significant digit"""
    return round(x, -int(floor(log10(abs(x)))))

def round_to_2(x):
    """ round to one significant digit"""
    return round(x, -int(floor(log10(abs(x))))+1)

# from https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy

def roundAndFormat1(xxx):
    """ changes numbers xxx into string percentages that are nice integer numbers 
    it is a hack to make sure that numbers like 65.241 become 65, whereas numbers like 6.7 become 7
    so its basically rounding to the nearest integer number and then formatting it nicely for legend
    """
    st = '{:0.2}'.format(xxx) # round
    st = (round(float(st),0))
    st = str(st)
    st = ('%f' %float(st)).rstrip('0').rstrip('.') # format
    return str(st)

def roundAndFormat(xxx):
    st = '{:0.2}'.format(xxx) # round
    st = ('%f' %float(st)).rstrip('0').rstrip('.') # format
    return str(st)





def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)








# GW DATA 
# array with gravitational wave estimates from LIGO-Virgo public data
#name m1 m1_plus m1_minus m2 m2_plus m2_minus spin spin_plus spin_minus chirp chirp_plus chirp_minus, massratio massratio_plus massratio_minus
# FROM https://www.gw-openscience.org/catalog/GWTC-1-confident/html/
# CATALOG PAPER 
gw150914 = [35.6, 4.8, 3,        30.6, 3, 4.4,       -0.01, 0.12, 0.13,     28.6, 1.7, -1.5,            0,0,0]                     
gw151012 = [23.3, 14, 5.5,       13.6, 4.1, 4.8,      0.04, 0.28, 0.19,     15.2, 2.1, -1.2,            0,0,0]                     
gw151226 = [13.7, 8.8, 3.2,      7.7, 2.2, 2.6,       0.18, 0.2, 0.12,      8.9, 0.3, -0.3,             0,0,0]           
gw170104 = [31, 7.2, 5.6,        20.1, 4.9, 4.5,     -0.04, 0.17, 0.2,      21.4, 2.2, -1.8,            0,0,0]       
gw170608 = [10.9 ,5.3 ,1.7,      7.6, 1.3, 2.1,       0.03, 0.19, 0.07,     7.9, 0.2, -0.2,             0,0,0]       
gw170729 = [50.6, 16.6, 10.2,    34.3, 9.1, 10.1,     0.36, 0.2, 0,         35.4, 6.5, -4.8,            0,0,0]       
gw170809 = [35.2 ,8.3, 6 ,       23.8 ,5.2, 5.1,      0.07 ,0.16 ,0.16,     24.9, 2.1, -1.7,            0,0,0]    
gw170814 = [30.7 ,5.7, 3,        25.3, 2.9, 4.1,      0.07, 0.12, 0.11,     24.1, 1.4, -1.1,            0,0,0]   
gw170817 = [1.46 ,0.12 ,0.1,     1.27 ,0.09 ,0.09,    0 ,   0.02 ,0.01,     1.186, 0.001, -0.001,       0,0,0]    
gw170818 = [35.5 ,7.5 ,4.7,      26.8, 4.3 ,5.2,     -0.09, 0.18, 0,        26.5, 2.1, -1.7,            0,0,0]   
gw170823 = [39.6 ,10 ,6.6,       29.4, 6.3, 7.1,      0.08, 0.2, 0.22,      29.2, 4.6, -3.6,            0,0,0]


#mass ratio GW170817 from https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.119.161101 
qrange_gw170817 = [0.7,1]
Zrange_gw170817 = np.asarray([0.2, 1])*0.0142 #from https://arxiv.org/pdf/1710.05861.pdf

GWdata = [gw150914, gw151012, gw151226, gw170104,gw170608,gw170729 , gw170809, gw170814, gw170817, gw170818, gw170823]







def QinBHspinmodel(separationPreSN2, M1, M2, maskNSBH):
    # returns spins from Qin et al + 2018 model 
    
    # start with all zeseparationDCOFormationarationDCOFormationBH spins
    BHspins = np.zeros_like(separationPreSN2)
    
    # now add spins for NS-BH following Qin et al 2018:
    # this is based on the separation prior to the second SN  
    PeriodPreSN2 = convert_a_to_P_circular(separation=separationPreSN2*u.Rsun, M1=M1*u.Msun, M2=M2*u.Msun)
    PeriodPreSN2 = PeriodPreSN2.to(u.d).value # in days 
    
    # only do non zero spins
    # first mask super tight NSBH that will get spin 1
    maskNSBHChi1 = (np.log10(PeriodPreSN2) < -0.3) & (maskNSBH ==1)
    BHspins[maskNSBHChi1] = np.ones(np.sum(maskNSBHChi1)) # fill with ones
#     print('#total, = ', len(maskNSBHChi1))
#     print('# with Chi = 1, = ', np.sum(maskNSBHChi1))
    
    # now the variable spin
    maskNSBHChi_var = (np.log10(PeriodPreSN2) > -0.3) &  (np.log10(PeriodPreSN2) < 0.3)  &(maskNSBH ==1)
    m_, c_ = -5./3, 0.5 # from Qin + 2018 
    spins_var =  m_ * np.log10(PeriodPreSN2[maskNSBHChi_var])  + c_   
    BHspins[maskNSBHChi_var] = spins_var
#     print('# with Chi var = ', np.sum(maskNSBHChi_var))
    
    return BHspins











def convert_a_to_P_circular(separation, M1, M2):
    """calculate Period from separation
    separation is separation (needs to be given in astropy units)
    M1 and M2 are masses of the binary
    
    """
    G = const.G # [gr cm s^2]
    

    mu = G*(M1+M2)
    period = 2*np.pi * np.sqrt(separation**3/mu)
    
    
    
    return period








    
def getXmomentOfMT(Seeds, maxCounter=10):
    # from Neijssel + 19 
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







### OLD CODE FOR COLORS 

# DEFINE COLOR LIST 

# colorlist = ['', '', '' , '', '', '', '', '']
# colorlist[7] = 'gray'

# for nrC in range(7):
#     if (nrC>0) & (nrC<4):
#         cm       = plt.get_cmap('viridis_r')
#         indd = nrC +2

#     if (nrC>=4) & (nrC<=5):
#         cm       = plt.get_cmap('plasma_r')
#         indd = nrC
#     if nrC>5 or nrC==0:
#         cm       = plt.get_cmap('plasma_r')
#         if nrC==0:
#             indd=0 
#         else:
#             indd=1
   
        
#     #             cmapCustom = matplotlib.colors.LinearSegmentedColormap.from_list("", [   "white", cm[i]])    
#     mycolors = [cm(x) for x in np.linspace(0,1 , (7))] 
#     colorlist[nrC] = mycolors[indd]


# colorlist = ['', '', '' , '', '', '', '', '']
# colorlist[7] = 'gray'

# for nrC in range(7):
#     if (nrC>=0) & (nrC<=2):
#         cm       = plt.get_cmap('plasma_r')
#         indd = nrC

#         mycolors = [cm(x) for x in np.linspace(0,1 , (6))] 
#         colorlist[nrC] = mycolors[indd] 

#     if nrC==3:
#         colorlist[nrC]='red'


#     if (nrC>=4) & (nrC<=5):
#         cm       = plt.get_cmap('plasma_r')
#         indd = nrC-1
#         mycolors = [cm(x) for x in np.linspace(0,1 , (6))] 
#         colorlist[nrC] = mycolors[indd] 

#     if nrC==6:
#         colorlist[nrC]='green'


# colorlist = ['', '', '' , '', '', '', '', '']
# colorlist[7] = 'gray'

# for nrC in range(7):
#     if (nrC>=0) & (nrC<=3):
#         cm       = plt.get_cmap('plasma')
#         indd = nrC+1

#         mycolors = [cm(x) for x in np.linspace(0,1 , (5))] 
#         colorlist[nrC] = mycolors[indd] 

# #     if nrC==5:
# #         colorlist[nrC]='lightblue'


#     if (nrC>=4) & (nrC<=6):
#         cm       = plt.get_cmap('viridis_r')
#         indd = nrC-3
#         mycolors = [cm(x) for x in np.linspace(0,1 , (5))] 
#         colorlist[nrC] = mycolors[indd] 

#     if nrC==6:
#         colorlist[nrC]='green'
        



    # #             cmapCustom = matplotlib.colors.LinearSegmentedColormap.from_list("", [   "white", cm[i]])    
    # mycolors = [cm(x) for x in np.linspace(0,1 , (7))] 
    # colorlist[nrC] = mycolors[indd]

# colorlist = ['', '', '' , '', '', '', '', '']
# colorlist[7] = 'gray'




# own definitions
# colorlist = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 'gold'] # colours of channels 
# colorlist = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 'deepskyblue', 'gold','#8c564b', 'gold'] # colours of channels 

# ChannelLabelList = ['channel 1','channel 2','channel 3','channel 4','channel 5','channel 6', 'channel 7' ]  # labels of channels  
# ChannelLabelListShort = ['1','2','3','4','1b','2b', '5' ] # shorter notation of ChannelLabelList




