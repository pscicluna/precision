from __future__ import print_function

#routines for dark and bias subtraction
import numpy as np
import astropy.io.fits as fits
from .utility import Observation

class Dark(object):
    ''' 
    Future class for making everything more OOP-sie
    '''
    def __init__(self):
        pass

    def __repr__(self):
        pass

    def __str__(self):
        pass

    

def loaddarks(darkframes=None,darkdir='',**kwargs):
    darks=np.array([])
    #print darkdir,darkframes
    if isinstance(darkframes,basestring):
        f=darkdir+darkframes+".fits"
        hdu=fits.open(f)
        #remember to add exposure time info!
        try:
            darks=np.r_[darks,hdu[0].data]
        except ValueError:
            darks=hdu[0].data
            darks=np.r_[darks,hdu[0].data/hdu[0].header['EXPTIME']] #counts/s ?
        hdu.close()
        #pass
    elif isinstance(darkframes,np.ndarray) or isinstance(darkframes,list):
        for f in darkframes:
            hdu=fits.open(f)
            #remember to add exposure time info!
            darks=np.r_[darks,hdu[0].data]
            darks=np.r_[darks,hdu[0].data/hdu[0].header['EXPTIME']] #counts/s ?
            hdu.close()
    return darks

def makemasterdark(darkframes=None,computevariance=True,computeRON=False,**kwargs):
    print(darkframes)
    #take many dark frames as input, then 
    if darkframes is not None:
        #check what has been passed
        if isinstance(darkframes,Observation):
            print(darkframes.obs_main[0][:][1],darkframes.datadir)
            darkframes=loaddarks(darkframes.obs_main[0][:][1],darkdir=darkframes.datadir)
            
            #for f in darkframes.obs_main:
            #    pass
            pass
        elif isinstance(darkframes[0],basestring):
            #contains strings, need to read files
            #masterdark=np.nanmedian(loaddarks(darkframes),axis=0)
            darkframes=loaddarks(darkframes)
        elif isinstance(darkframes[0],np.ndarray):
            #contains fits data already
            pass # horrible hack :)
        else:
            raise TypeError("You need to pass either data arrays or the names of files containing dark frames")
        pass        
    else:
        raise ValueError("Dark frames must be defined either as file names or data previously read in")
    masterdark=np.nanmedian(darkframes,axis=0) #check that this is robust for cubes :/
    #calculate variance
    if (computevariance):
        if len(darkframes.shape) == 3:
            darkvar=(np.pi/(2.*darkframes.shape[0]))*(np.std(darkframes,axis=0))**2.
        else:
            darkvar=masterdark
    # return master dark frame1
    if computeRON: #for testing only. Expensive to take so many differences and create so many subapertures. Based on algorithm described in NACO pipeline manual. Also provide same routine in separate routine if needed.
        pass
    else:
        RON = np.nan
    return masterdark,darkvar,RON

def darksubtract(scienceframes=None,darkframes=None,**kwargs):
    #performs dark subtraction from science frames
    if darkframes is not None:
        #check what has been passed
        if isinstance(darkframes[0],basestring):
            #contains strings, need to read files
            if isinstance(darkframes,[list,ndarray]):
                #many files, need to create master dark
                masterdark,darkvar=makemasterdark(darkframes)
            else:
                #only one file, treat it as master, read in
                masterdark=makemasterdark(darkframes,computevariance=False)[0]
                darkvar=np.sqrt(masterdark)
        elif isinstance(darkframes[0],ndarray):
            #contains fits data already
            if isinstance(darkframes,[list,ndarray]):
                #many arrays, need to create master dark
                masterdark=makemasterdark(darkframes)
            else:
                #only one array, treat it as master, performs subtraction
                masterdark=darkframes
                darkvar=np.sqrt(masterdark)
        else:
            raise TypeError("You need to pass either numpy data arrays or the names of files containing dark frames")
    else:
        raise ValueError("Dark frames must be defined either as file names or data previously read in")
    return scienceframes-masterdark,masterdark,darkvar,RON

def biassubtract():
    pass

def readnoise():
    pass
