#routines for flatfielding, etc
import numpy as np
import astropy.io.fits as fits
from astropy.stats import sigma_clip

def loadflats(flatframes=None,masterdark=None):
    flats=np.array([])
    for f in flatframes:
        hdu=fits.open(f)
        if masterdark is None: #check to see if dark subtraction has been requested
            masterdark = np.zeros_like(hdu[0].data)
        if len(hdu[0].data.shape) == 3): #it's a cube!
#            flat=np.zeros([hdu[0].data.shape])
            for i in range(hdu[0].data.shape[0]):
                median=np.nanmedian(hdu[0].data[i,:,:])
                flat=(hdu[0].data[i,:,:]-masterdark)/median
                flats=np.r_[flats,flat]
        else:
            median=np.nanmedian(hdu[0].data)
            flat=((hdu[0].data-masterdark)/median))
            flats=np.r_[flats,flat]#hdu[0].data] #(hdu[0].data/median) #needs revision in case of datacube
        hdu.close()
    return flats

def makemasterflat(flatframes=None,masterdark=None,computevar=True,method='med',computebadpix=True,sig=5.0):
    if flatframes is not None:
        #check what has been passed
        if isinstance(flatframes[0],basestring):
            #contains strings, need to read files
            flatframes=loadflats(flatframes,masterdark)
        elif isinstance(flatframes[0],np.ndarray):
            #contains fits data already, but needs to be normalised and dark-subtracted
            pass # horrible hack :)
        else:
            raise TypeError("You need to pass either data arrays or the names of files containing flat frames")
    else:
        raise ValueError("Flat frames must be defined either as file names or data previously read in")
    masterflat=np.nanmedian(flatframes,axis=0)
    if computevar:
        flatvar=(np.pi/(2.*flatframes.shape[0]))*(np.std(flatframes,axis=0))**2
    if computebadpix:
        badpixelmap=np.zeros_like(masterflat)
        badpixelmap[np.unravel_index(sigma_clip(masterflat,sig=5.0).mask,badpixelmap.shape)] = 1 #test this thoroughly!! 5-sigma outliers chosen as, given size of IRDIS detector, there should be ~1 false bad pixel if pixel values are normally distributed; exact value should be determined through experimentation
        if computevar:
            return masterflat,flatvar,badpixelmap
        else:
            return masterflat,badpixelmap
    else:
        if computevar:
            return masterflat,flatvar
        else:
            return masterflat


def flatfield(scienceframe=None,flatframes=None,darkframes=None,masterdark=None):
    if flatframes is not None:
        if isinstance(flatframes[0],basestring):
            if isinstance(flatframes,[list,ndarray]):
                masterflat,flatvar=makemasterflat(flatframes,masterdark)
            else:
                masterflat=matmasterflat(flatframes,masterdark,computevar=False)[0]
                flatvar=np.sqrt(masterflat)
        elif isinstance(flatframes[0],ndarray):
            if isinstance(flatframes,[list,ndarray]):
                masterflat,flatvar=makemasterflat(flatframes,masterdark)
            else:
                masterflat=flatframes-masterdark
                flatvar=np.sqrt(masterflat)
        else:
            raise TypeError("You need to pass either numpy data arrays or the names of files containing flat frames")
    else:
        raise ValueError("Flat frames must be defined either as file names or data previously read in")
    pass


