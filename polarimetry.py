#routines for producing polarimetric images, fits files and plots, including binning of data for clean presentation

import numpy as np
import math as mt

def calcQU(Oray,Eray,**kwargs):
    """
    Calculates linear Stokes parameter from observed O-ray and E-ray intensities.

    For use when each Stokes parameter is only measured at a single angle or when using the double difference method of polarimetric differential imaging
    """
    return Oray - Eray

def calcpQU(O0,E0,O45,E45,**kwargs):
    """
    Calculates fractional Stokes parameter based on PDI method to cancel instrumental and atmospheric effects
    """
    R=np.sqrt((O0-E0)/(O45-E45))
    return (R-1.)/(R+1.)

def polfrac(I,Q=None,U=None,V=None,**kwargs):
    ind=np.where(I > 0.)
    poldeg=np.zeros_like(I)
    if V is none:
        poldeg[ind]=(np.sqrt(Q[ind]*Q[ind] + U[ind]*U[ind]))/I[ind]
    else:
        poldeg[ind]=(np.sqrt(Q[ind]*Q[ind] + U[ind]*U[ind] + V[ind]*V[ind]))/I[ind]
    return poldeg

def polang(I,Q,U,**kwargs):
    gamma=np.zeros_like(I)
#    ind=np.where(Q > 0. & U > 0.)
    ind=np.where((Q[i,j] > 0.) & (U[i,j] > 0.))
    gamma[ind] = 0.5*np.arctan(U[ind]/Q[ind])
    ind=np.where( (Q[i,j] < 0.) & (U[i,j] < 0.))
    gamma[ind] = 0.5*(np.arctan(abs(U[ind])/abs(Q[ind])) + np.pi)
    ind=np.where( (Q[i,j] < 0.) & (U[i,j] > 0.))
    gamma[ind] = 0.5*np.arctan2(U[ind],Q[ind])
    ind=np.where((Q[i,j] > 0.) & (U[i,j] < 0.))
    gamma[ind] = 0.5*np.arctan2(U[ind],Q[ind]) + np.pi
    ind=np.where((Q[i,j]==0) & (U[i,j]==0))
    gamma[ind] = 0.    
    ind=np.where(( abs(Q[i,j]) == 0) & (U[i,j] < 0))
    gamma[ind] = 0.75*np.pi
    ind=np.where( (abs(Q[i,j]) == 0) & (U[i,j] > 0))
    gamma[ind] = 0.25*np.pi
    ind=np.where( (abs(U[i,j]) == 0) & (Q[i,j] < 0))
    gamma[ind] = 0.5*np.pi   
    ind=np.where( (abs(U[i,j]) == 0) & (Q[i,j] > 0))
    gamma[ind] = 0.     
    return np.degrees(gamma)

def poldegang(I,Q,U,noang=False,**kwargs): #for getting polarisation angle and degree at same time
    ind=np.where(I > 0.)
    poldeg=np.zeros_like(I)
    poldeg[ind]=(np.sqrt(Q[ind]*Q[ind] + U[ind]*U[ind]))/I[ind]
    gamma=np.zeros_like(I)
#    ind=np.where(Q > 0. & U > 0.)
    if (noang==False):
        ind=np.where((Q[i,j] > 0.) & (U[i,j] > 0.))
        gamma[ind] = 0.5*np.arctan(U[ind]/Q[ind])
        ind=np.where( (Q[i,j] < 0.) & (U[i,j] < 0.))
        gamma[ind] = 0.5*(np.arctan(abs(U[ind])/abs(Q[ind])) + np.pi)
        ind=np.where( (Q[i,j] < 0.) & (U[i,j] > 0.))
        gamma[ind] = 0.5*np.arctan2(U[ind],Q[ind])
        ind=np.where((Q[i,j] > 0.) & (U[i,j] < 0.))
        gamma[ind] = 0.5*np.arctan2(U[ind],Q[ind]) + np.pi
        ind=np.where((Q[i,j]==0) & (U[i,j]==0))
        gamma[ind] = 0.    
        ind=np.where(( abs(Q[i,j]) == 0) & (U[i,j] < 0))
        gamma[ind] = 0.75*np.pi
        ind=np.where( (abs(Q[i,j]) == 0) & (U[i,j] > 0))
        gamma[ind] = 0.25*np.pi
        ind=np.where( (abs(U[i,j]) == 0) & (Q[i,j] < 0))
        gamma[ind] = 0.5*np.pi   
        ind=np.where( (abs(U[i,j]) == 0) & (Q[i,j] > 0))
        gamma[ind] = 0.     
    return poldeg,np.degrees(gamma)

def rebin(arr,shape,**kwargs):
    sh=shape[0],arr.shape[0]//shape[0],shape[1],arr.shape[1]//shape[1]
    return arr.reshape(sh).mean(-1).mean(1) #switch to nanmean() if NaNs are present to average over non-NaN values

def imbin(image,nbins,**kwargs):
    if (mt.fmod(image.shape[0],nbins[0]) == 0) & (mt.fmod(image.shape[1],nbins[1]) == 0):
        #this case works
        binnedimage=rebin(image,nbins)
    else:
        raise binerror('For now, image.shape % nbins must be 0')
    return binnedimage
