import numpy as np
from scipy import signal as sig
from scipy.optimize import curve_fit,Bounds

"""Functions to fit a distribution of dihedral angles to a sum of Gaussian distributions 
   and returns the mean and sigma values for each Gaussian.

   Adapted from code by Alexandre Ferreira. 
"""

def _histogram(dihedral_angles: str, datacol: int=1, probcutoff: float=0.0075):
    """Compute the probability distribution given a timeseries of dihedral angles

       centers,histo (list,list) : centers of the bins and values of the probability distribution.

    :param dihedral_angles: text file with the timeseries. Angle values are assumed to be in the second column of the file.
    :type dihedral_angles: str
    :param datacol: column (0-indexed) where data is. Defaults to 1.
    :type datacol: int, optional
    :param probcutoff: cutoff value for the population of a bin to take into account in the final histogram. Defaults to 0.001
    :type probcutoff: float, optional
    :returns centers: centers of the bins of the probability distribution.
    :rtype centers: list
    :returns histo: values of the histogram of the probability distribution
    :rtype histo: list

    """
    
    commentchars=['#','@']
    l=-np.pi	      # left-hand edge value
    r=np.pi	      # right-hand edge value
    bw=5*np.pi/180    # 5 degrees bin

    nb=int((r-l)/bw)  # number of bins
    
    histo=[0 for i in range(nb)]
    centers=[l+bw*(i+0.5) for i in range(nb)]

    n=0
    
    with open(dihedral_angles,'r') as f:
        for line in f:
            if line[0] not in commentchars:
                cols=line.split()
                b=int((float(cols[datacol])-l)/bw)
                histo[b]+=1
                n+=1

    for b in range(nb):
        histo[b]/=n
        real_histo=histo
        histo[b]*=(histo[b]>probcutoff)

    # in order to be able to treat also the cases where the peak is located around +-pi [rad]
    # the histogram is replicated between -3*pi and 3*pi
    # this increases the number of peaks that are found at a very low computational cost.
    centers=[c-2*np.pi for c in centers]+centers+[c+2*np.pi for c in centers]
    histo=histo*3
        
    return centers,histo


def _sumGaussFunct(x: list, *parameters: list):
    """Sum of Gaussian functions
 
    :param x: x axis values
    :type x: list
    :param parameters: list of parameters of length 3*N where N is the number of Gaussians to estimate. Order of the parameters is A,mu,sigma.
    :type parameters: list
    :returns f: y axis values
    :rtype f: list

    """
    
    ngauss=int(len(parameters)/3)

    a=np.asarray(parameters[0:ngauss])
    mu=np.asarray(parameters[ngauss:2*ngauss])
    s=np.asarray(parameters[2*ngauss:3*ngauss])

    f=[ 0 for c in x ]
    
    for i in range(ngauss):

        a_val=a[i]
        mu_val=mu[i]
        s_val=s[i]

        for ic,c in enumerate(x):

            f[ic] += a_val*np.exp(-((c-mu_val)**2)/(2*s_val**2))

    return f

def gaussianFit(dihedral_angles: str, datacol: int=1, plotFit: bool=False):
    """Estimate fit parameters

    :param dihedral_angles: text file with the values of dihedral angles
    :type dihedral_angles: str
    :param datacol: column (0-indexed) where data is. Defaults to 1.
    :type datacol: int, optional
    :param plotFit: plot the original histogram and the fitted gaussians? Defaults to False.
    :type plotFit: bool, optional
    :returns centers: centers of the fitted distributions
    :rtype centers: list
    :returns sigma: sigma values of the fitted distributions
    :rtype sigma: list  

    """
    
    bincenters,histo=_histogram(dihedral_angles,datacol=datacol)

    cutoffdist=20
    
    # look for peaks with a inter-peak distance of at least cutoffdist
    peaks,properties=sig.find_peaks(histo, distance=cutoffdist)  

    ngauss=len(peaks)
    guess_a=[ histo[p] for p in peaks ]
    guess_mu=[ bincenters[p] for p in peaks ]
    guess_sigma=[5*np.pi/180 for p in peaks ]   # 5 degrees

    sumf = _sumGaussFunct(bincenters,guess_a,guess_mu,guess_sigma)

    #bounds=[ tuple(np.array([0]*ngauss+[-10]*ngauss*2).flatten()), tuple(np.array([10]*ngauss*3).flatten())  ]

    bounds=Bounds(lb=np.array([0]*ngauss+[-10]*ngauss*2).flatten(), ub=np.array([3]*ngauss+[10]*ngauss*2).flatten())

    
    fitparms=np.array(guess_a+guess_mu+guess_sigma)

    popt,pcov=curve_fit(_sumGaussFunct,bincenters,histo,p0=np.array([guess_a, guess_mu, guess_sigma]).flatten(),bounds=bounds)

    if (peaks[0]+(len(histo)-peaks[-1])<cutoffdist):

        if histo[peaks[0]]>=histo[peaks[-1]]:
            peaks=np.delete(peaks,-1)
            index=range(ngauss,len(popt),ngauss)
            popt=np.delete(popt,index)
        else:
            peaks=np.delete(peaks,0)
            index=range(0,len(popt),ngauss)
            popt=np.delete(popt,index)
            
        ngauss-=1

    temp_centers=popt[ngauss:2*ngauss]
    temp_sigma=popt[2*ngauss:3*ngauss]

    heights=[]
    centers=[]
    sigmas=[]

    for icenter,center in enumerate(temp_centers):
        if center>=-np.pi and center<=np.pi:
            heights.append(popt[0:ngauss][icenter])
            centers.append(center)
            sigmas.append(temp_sigma[icenter])

        
    fitFunct=_sumGaussFunct(bincenters,heights,centers,sigmas)

    originFunct=_sumGaussFunct(bincenters,[guess_a,guess_mu,guess_sigma])

    
    
        

    

    # for debug purposes it is possible to plot the original
    # histogram and the fitted gaussians
    if plotFit:
    
        import matplotlib.pyplot as plt

        fig,ax=plt.subplots()

        ax.set_xlim([-np.pi,np.pi])
        
        ax.plot(bincenters,histo,label='data',linewidth=2)
        ax.plot(bincenters,fitFunct,label='fit',linestyle='--')

        #ax.legend()
        
        plt.show()
    
    return centers,sigmas

