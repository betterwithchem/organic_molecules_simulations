import numpy as np
from scipy import signal as sig
from scipy.optimize import curve_fit

"""Fit a distribution of dihedral angles to a sum of Gaussian distributions 
   and returns the mean and sigma values for each Gaussian 
"""


def histogram(dihedral_angles: str, probcutoff=0.001):
    """Compute the probability distribution given a timeseries of dihedral angles
    Args:
       dihedral_angles (str) : text file with the timeseries. Angle values are assumed to be in the second column of the file.
       probcutoff (float, optional) : cutoff value for the population of a bin to take into account in the final histogram.

    Returns:
       centers,histo (list,list) : centers of the bins and values of the probability distribution.
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
                b=int((float(cols[1])-l)/bw)
                histo[b]+=1
                n+=1

    for b in range(nb):
        histo[b]/=n
        real_histo=histo
        histo[b]*=(histo[b]>probcutoff)
        
    return centers,histo


def sumGaussFunct(x,*args):
    """Sum of Gaussian functions
 
    Args:
       x (list): x axis values
       *args (list): list of parameters of length 3*N where N is the number of Gaussians to estimate. Order of the parameters is A,mu,sigma.

    Returns:
      f (list): y axis values
    """
    
    ngauss=int(len(args)/3)

    a=np.asarray(args[0:ngauss])
    mu=np.asarray(args[ngauss:2*ngauss])
    s=np.asarray(args[2*ngauss:3*ngauss])

    f=[ 0 for c in x ]
    
    for i in range(ngauss):

        a_val=a[i]
        mu_val=mu[i]
        s_val=s[i]

        for ic,c in enumerate(x):

            f[ic] += a_val*np.exp(-((c-mu_val)**2)/(2*s_val**2))

    return f



def gaussianFit(dihedral_angles: str, datacol=1):
    """Estimate fit parameters

    Args:
       dihedral_angles (str) : text file with the values of dihedral angles
       datacol (int, optional) : column (0-indexed) where data is. Defaults to 1.

    Returns:
       centers (list) : mean values of each Gaussian curve obtained from the fit.
       sigma (list) : standard deviation values of each Gaussian curve obtained from the fit.
    """
    
    bincenters,histo=histogram(dihedral_angles)

    cutoffdist=20
    
    # look for peaks with a inter-peak distance of at least cutoffdist
    peaks,properties=sig.find_peaks(histo, distance=cutoffdist)  

    ngauss=len(peaks)
    guess_a=[ histo[p] for p in peaks ]
    guess_mu=[ bincenters[p] for p in peaks ]
    guess_sigma=[5*np.pi/180 for p in peaks ]   # 5 degrees

    sumf = sumGaussFunct(bincenters,guess_a,guess_mu,guess_sigma)

    popt,pcov=curve_fit(sumGaussFunct,bincenters,histo,p0=[guess_a,guess_mu,guess_sigma])

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
    
    fitFunct=sumGaussFunct(bincenters,popt[0:ngauss],popt[ngauss:2*ngauss],popt[2*ngauss:3*ngauss])

    originFunct=sumGaussFunct(bincenters,[guess_a,guess_mu,guess_sigma])
    
    centers=[popt[ngauss:2*ngauss]]
    sigma=[popt[2*ngauss:3*ngauss]]
    
    
    return centers,sigma
