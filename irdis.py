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
from image_registration.register_images import register_images
#from skimage.feature import peak_local_max, canny


#need a list of filters, conoragraphs, etc
class Irdis(Instrument):
    def __init__(self,scifiles=None,badpixfile=None,skyfiles=None,
                 flatfiles=None,darkfiles=None,mode=None,optics=None,
                 options=None,**kwargs):
        self.scifiles=scifiles
        self.badpixfile=badpixfile
        self.skyfiles=skyfiles
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

    def splitchannels(self,data,**kwargs):
        #call this after dark/flat/sky correction and centre finding but before anything else - will probably end up most used...
        return data[:,:,0:1024],data[:,:,1024:2048]

    def CI(self,**kwargs): #test on VY CMa
        #dark/flat/sky correction (optional)

        #undither and concatenate exposure cubes

        #optional - derotate

        #split L and R channels and concatenate cubes

        #collapse cubes

        #flux
        pass

    def SDI(self,**kwargs): #test on GD50, VY CMa

        #dark/flat/sky correction (optional)

        #undither and concatenate exposure cubes

        #optional - derotate

        #split channels - create concatenated cube but also collapse cube of each channel and do subtraction
        pass

    def ADI(self,**kwargs): #test on GD50
        pass
    
    def LOCI(self,**kwargs):
        pass

    def DPI(self,**kwargs): #test on HR 3090
        #perform dark/flat/sky correction (optional)

        #undither

        #split L and R channels into O- and E-rays respectively
        Oray,Eray=self.splitchannels()
        
        #calculate Stokes parameter(s) the old fashioned way
        #concatenate channels and collapse for I 

        #O - E for Q/U

        #if double difference, do +Q - -Q

        #compute fractional pol

        #do ratio method too
        
        pass

    def reduce(self,**kwargs):
        self.masterdark,self.darkvar,self.RON=darkbias.makemasterdark(self.darkframes,**kwargs)
        
        self.masterflat,self.flatvar,self.badpixmap=flatfield.makemasterflat(self.flatframes,self.masterdark,**kwargs)
        
        for f in self.skyfiles:
            data,header=self.readdata(f)
            data=self.badpixcorrect(data)
            self.skyframes=np.r_[self.skyframes,data]
        #now median combine all sky frames
        if (len(self.skyframes.shape) == 3):
            self.mastersky=np.nanmedian(self.skyframes,axis=0)
            self.skyvar=(np.pi/(2.*self.skyframes.shape[0]))*(np.std(self.skyframes,axis=0))**2
            self.skyfiles=None
        else:
            self.mastersky=self.skyframes
            self.skyvar=self.mastersky

        self.medianNoDerot=np.array([])
        self.medianDerot=np.array([])
        self.varNoDerot=np.array([])
        self.varDerot=np.array([])
        for f in self.scifiles:
            #read data
            data,header=self.readdata(f)
            #intermediate processing
            data=self.badpixcorrect(data)
            data=self.skysub(data,self.mastersky)
            #split channels and align
            datal,datar=self.splitchannels(data)
            chanshift=register_images(datal[0],datar[0],usamp=10.)
            datar=simage.interpolation.shift(datar,chanshift)
            #data=np.r_[datal,datar]
            #datal=None
            #datar=None

            #----------------------------NOW THINGS ARE DIFFERENT DEPENDING ON METHOD!!--------------------------
            self.medianNoDerot=np.r_[self.medianNoDerot,
                                     [simage.interpolation.shift(np.nanmedian(datal,axis=0),
                                                                [header['HIERARCH ESO INS1 DITH POSX'],
                                                                 header['HIERARCH ESO INS1 DITH POSY']]),
                                      simage.interpolation.shift(np.nanmedian(datar,axis=0),
                                                                 [header['HIERARCH ESO INS1 DITH POSX'],
                                                                  header['HIERARCH ESO INS1 DITH POSY']])
                                      ]
                                     ] #make sure the sign is right here (by inspection!!)
            self.varNoDerot=np.r_[self.varNoDerot,
                                  [simage.interpolation.shift((np.pi/(2.*datal.shape[0]))*(np.std(datal,axis=0))**2,
                                                              [header['HIERARCH ESO INS1 DITH POSX'],
                                                               header['HIERARCH ESO INS1 DITH POSY']]),
                                   simage.interpolation.shift((np.pi/(2.*datar.shape[0]))*(np.std(datar,axis=0))**2,
                                                              [header['HIERARCH ESO INS1 DITH POSX'],
                                                               header['HIERARCH ESO INS1 DITH POSY']])
                                   ]
                                  ]
            #now derotate cube if pupil stabilised
            if self.rot=='PUPIL':
                pass
            #first find centre of rotation
                
            #then calculate rotation as a function of time

            #then derotate each frame of each half of the detector
            self.sciframes

        if self.mode=='CI':
            print 'This mode is not available yet'
        elif self.mode=='DPI':
            print 'This mode is not available yet'
        elif self.mode=='SDI':
            print 'This mode is not available yet'
        else:
            print 'IRDIS mode not recognised'

    def readdata(self,filename,**kwargs):
        hdu=fits.open(filename)
        cube=hdu[0].data #extract data itself
        #then extract important header info ? (might have already done this before, not sure about architecture yet)
        cube=cube-self.masterdark
        cube=cube/self.masterflat
        header=hdu[0].header
        extra='??'
        hdu.close()
        return cube,header#,extra

class IrdisLSS(Irdis):
    pass
