import numpy as np
import scipy.ndimage as simage
import astropy.io.fits as fits
import math as mt
from scipy.signal import medfilt2d
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
import datetime
from __init__ import __version__ #hack, replace later
#from instrument import Instrument

class Instrument(object):
    #base class to be inherited by instrument-specific classes with details and redution routines, possibly simulation routines too
    """Generic instrument class. Contains things that all routines will need (e.g. dark subtraction, flat fielding, sky subtraction, bad pixel identification etc) while specific reduction routines are provided in extension classes. As much as possible, uncertainty tracking should be included here. """
    def __init__(self,scifiles=None,badpixfile=None,skyfile=None,flatfiles=None,darkfiles=None,mode=None,optics=None,**kwargs):
        self.scifiles=scifiles
        self.badpixfile=badpixfile
        self.skyfile=skyfile
        self.flatfiles=flatfiles
        self.darkfiles=darkfiles
        self.mode=mode #instrument modes, dictionary
        self.optics=optics #instrument optics, dictionary

    def __repr__(self):
        raise NotImplementedError("Please implement this in your subclass")
        
    def __str__(self):
        raise NotImplementedError("Please implement this in your subclass")

    def setbadpixelmap(self,**kwargs):
        hdu=fits.open(badpixfile)
        self.badpixmap=hdu[0].data
    
    def setskyframe(self,**kwargs):
        hdu=fits.open(self.skyfile)
        self.skyframe=hdu[0].data
    
    def setbadpix(self,sciframe,badpixmap,**kwargs):
        #set all bad pixels to NaN
        sciframe[badpixmap==1]=float('Nan') #does this work?
    
    def badpixcorrect(self,data=None,**kwargs):
        if data is None:
            data=self.sciframes
            recopy=True
        else:
            recopy=False
        if (len(data.shape) == 3): #datacube
            for i in range(data.shape[0]):
                #make sure to test how many times the windowing should be repeated - adding many small windows doens't cost too much - needs much optimisation though!
                mask=self.badpixmap#np.logical_or(self.badpixmap,np.isnan(data[i,:,:]))
                #print mask.shape
                
                #data[i,self.badpixmap]=medfilt2d(data[i,:,:],25)[self.badpixmap]
                data[i,mask]=medfilt2d(data[i,:,:],9)[mask]
                data[i,mask]=medfilt2d(data[i,:,:],5)[mask]
                data[i,:,:] = interpolate_replace_nans(data[i,:,:],Gaussian2DKernel(stddev=1))
                #exit()
        else: #only one image
            #data[self.badpixmap]=medfilt2d(data,25)[self.badpixmap]
            data[self.badpixmap]=medfilt2d(data,5)[self.badpixmap]
        if recopy:
            self.sciframes=data
        else:
            return data

    def readobservations(self,**kwargs):
        #read in all the science frames and compile a list of them.
        self.obslist=[]
        for f in self.scifiles:
            self.obslist.append(fits.open(f))
    
    def skysub(self,sciframe,skyframe,**kwargs):
        sciframe-=skyframe
        return sciframe

    def align(self,**kwargs):
        #align exposures prior to medianing
#        self.alignedframes=simage.interpolation.shift(
        pass
    
    def derotate(self,**kwargs):
        pass

    def logger(self,msg,**kwargs):
        """
        Writes logging information to a file

        Inputs:
        msg      The message to be written - a string or sequence of strings
        logfile  A file-like object to which msg should be written

        Outputs:
        None
        """
        logtime = str(datetime.datetime.now())
        outstr = logtime+": "+msg
        self.logfile.write(outstr)
        #raise NotImplementedError("Pipeline logging is not implemented yet")

    def updateHeader(self,**kwargs):
        pass

#future classes, they're here as reminders for now, everything will get moved around later as structure starts to become clearer and I get data to work on
class IFS(Instrument):
    pass

class NACO(Instrument):
    pass

class NACOSpec(NACO):
    pass

class VISIR(Instrument):
    pass

class VISIRSpec(VISIR):
    pass

class HiCIAO(Instrument):
    pass

class VAMPIRES(Instrument):
    pass

