import numpy as np
import scipy.ndimage as simage
import astropy.io.fits as fits
import math as mt
#from instrument import Instrument

class Instrument(object):
    #base class to be inherited by instrument-specific classes with details and redution routines, possibly simulation routines too
    """Generic instrument class. Contains things that all routines will need (e.g. dark subtraction, flat fielding, sky subtraction, bad pixel identification etc) while specific reduction routines are provided in extension classes. As much as possible, uncertainty tracking should be included here. """
    def __init__(self,scifiles=None,badpixfile=None,skyfile=None,flatfiles=None,darkfiles=None,mode=None,optics=None):
        self.scifiles=scifiles
        self.badpixfile=badpixfile
        self.skyfile=skyfile
        self.flatfiles=flatfiles
        self.darkfiles=darkfiles
        self.mode=mode #instrument modes, dictionary
        self.optics=optics #instrument optics, dictionary

    def __str__(self):
        pass

    def setbadpixelmap(self):
        hdu=fits.open(badpixfile)
        self.badpixmap=hdu[0].data
    
    def setskyframe(self):
        hdu=fits.open(self.skyfile)
        self.skyframe=hdu[0].data
    
    def setbadpix(sciframe,badpixmap):
        #set all bad pixels to NaN
        sciframe[badpixmap==1]=float('Nan') #does this work?

    def readobservations(self):
        #read in all the science frames and compile a list of them.
        self.obslist=[]
        for f in self.scifiles:
            self.obslist.append(fits.open(f))
    
    def skysub(sciframe,skyframe):
        sciframe-=skyframe
        return sciframe

    def align(self):
        #align exposures prior to medianing
#        self.alignedframes=simage.interpolation.shift(
        pass
    
    def derotate(self):
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

