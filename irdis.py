import numpy as np
import scipy.ndimage as simage
import astropy.io.fits as fits
import math as mt
from instrument import Instrument
#import astropysics.ccd as ccd
import photutils as pu
import skimage as skim
from skimage.transform import hough_circle
from skimage.util import img_as_ubyte
import darkbias
import flatfield
#from skimage.feature import peak_local_max, canny


#need a list of filters, conoragraphs, etc
class Irdis(Instrument):
    def __init__(self,scifiles=None,badpixfile=None,skyfile=None,flatfiles=None,
                 darkfiles=None,mode=None,optics=None,options=None):
        self.scifiles=scifiles
        self.badpixfile=badpixfile
        self.skyfile=skyfile
        self.flatfiles=flatfiles
        self.darkfiles=darkfiles
        self.mode=mode 
        self.optics=optics
        self.options=options
    
    def findcentre(self,**kwargs):
        pass
        #how on earth do I do this? try finding the coronagraph?
        #in case of ALC look for central spot using photutils
        #else use scikit-image to find circles for CLC
 #       tempimg=img_as_byte(frame)
 #       edges=canny(tempimg)
 #       hough_radii=1 #needs to depend on coronagraph choice
 #       hough_res=hough_circle(edges,hough_radii)
        #now need to use info to align both cubes
        #return centres of circles

    def splitchannels(self,**kwargs):
        #call this after dark/flat/sky correction and centre finding but before anything else - will probably end up most used...
        return self.sciframes[:,:,0:1024],self.sciframes[:,:,1024:2048]

    def CI(self,**kwargs): #test on VY CMa
        #
        pass

    def SDI(self,**kwargs): #test on GD50, VY CMa
        pass

    def ADI(self,**kwargs): #test on GD50
        pass
    
    def LOCI(self,**kwargs):
        pass

    def DPI(self,**kwargs): #test on HR 3090
        pass

class IrdisLSS(Irdis):
    pass
