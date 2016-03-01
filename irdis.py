import numpy as np
import scipy.ndimage as simage
import astropy.io.fits as fits
import astropy.stats as astats
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
    
    def findcentre(self,data,**kwargs):
        mean,median,std=astats.sigma_clipped_stats(data,sigma=3.0)
        sources=pu.daofind(data[500:550,450:500] - median, # assumes datacube #[500:550,450:500]
                           fwhm = 3.0,
                           threshold=5.*std,
                           sharplo=0.3,
                           sharphi=0.5,
                           roundhi=0.3,
                           roundlo=-0.3
                           )
        sources.sort('peak') #take brightest source
        try:
            self.central=sources[:][-1]
        except:
            self.central=sources
#        self.central=central
        self.sources=sources
        #return central,sources
#        pass
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

    def derotate(self,data,centre,**kwargs):
        #pad image so that centre is at the centre of the array
        padX=[data.shape[1] - centre[0],centre[0]] #check which is row and column in image and centre routines!
        padX=[data.shape[0] - centre[1],centre[1]]
        datap=np.pad(data,[padY,padX],'constant')
        #derotate
        data=simage.interpolation.rotate(datap,angle,reshape=False)[padY[0]:-padY[1],padX[0]-padX[1]]
        pass

    def CI(self,**kwargs): #test on VY CMa
        self.finalNoDerot=np.mean(self.medianNoDerot,axis=0) #produces L- and R- channel images
        self.finalNoDerot=np.mean(self.finalNoDerot,axis=0) #combine both channels
        self.finalvarNoDerot=np.mean(self.varNoDerot,axis=0)         #small number approximation - 
        self.finalvarNoDerot=(np.mean(self.finalvarNoDerot,axis=0) / #variance on mean = mean of variances / N
                              (2*self.varNoDerot.shape[0]))          #or alternatively
                                                                     #uncertainty on mean = mean of uncertainties / sqrt(N)
        pass

    def SDI(self,**kwargs): #test on GD50, VY CMa
        self.CI()

        self.SDINoDerot=np.mean(self.medianNoDerot,axis=0) #produces L- and R- channel images
        self.SDINoDerot=self.SDINoDerot[1]-self.SDINoDerot[0] #check which way round this is supposed to be! and scaling!

        self.varSDINoDerot=np.mean(self.varNoDerot,axis=0) / self.varNoDerot.shape[0]
        self.varSDINoDerot=np.sum(self.varSDINoDerot, axis=0)
        
        pass

    def ADI(self,**kwargs): #test on GD50
        #classical ADI, take median non-derotated frame and subtract it from each derotated frame, then collapse the whole cube
        pass
    
    def LOCI(self,**kwargs):
        pass

    def DPI(self,**kwargs): #test on HR 3090
        self.CI()

        #split L and R channels into O- and E-rays respectively
        self.MeanOray=np.mean(self.medianNoDerot[:,0,:,:],axis=0)
        self.MeanEray=np.mean(self.medianNoDerot[:,1,:,:],axis=0)# self.splitchannels()
        
        #calculate Stokes parameter(s) the old fashioned way
        #concatenate channels and collapse for I 
        self.I=self.MeanOray + self.MeanEray

        #O - E for Q/U
        self.pol=self.MeanOray - self.MeanEray #this all needs updating to enable it to process an entire DPI observation, correctly interpreting what combinations of +/- Q/U to do

        #if double difference, do +Q - -Q

        #compute fractional pol

        #do ratio method too
        
        pass

    def reduce(self,**kwargs):
        if self.darkframes is not None:
            self.masterdark,self.darkvar,self.RON=darkbias.makemasterdark(self.darkframes,**kwargs)
        else:
            self.masterdark=None

        if self.flatframes is not None:        
            self.masterflat,self.flatvar,self.badpixmap=flatfield.makemasterflat(self.flatframes,self.masterdark,**kwargs)
        else:
            self.masterflat=None
        
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
        self.headers=np.array([])
        isci=0
        for f in self.scifiles:
            #read data
            data,header=self.readdata(f)
            if isci==0:
                #pull important info out of header from first science file
                self.rot=header['HIERARCH ESO INS4 COMB ROT']
                self.parangInit=header['HIERARCH ESO TEL PARANG START']
                self.optics={filt: [header['HIERARCH ESO INS1 FILT NO'],header['HIERARCH ESO INS1 FILT ID'],header['HIERARCH ESO INS1 FILT NAME']],opti:[header['HIERARCH ESO INS1 OPTI2 NO'],header['HIERARCH ESO INS1 OPTI2 ID'],header['HIERARCH ESO INS1 OPTI2 NAME']],stop: [header['HIERARCH ESO INS1 OPTI1 NO'],header['HIERARCH ESO INS1 OPTI1 ID'],header['HIERARCH ESO INS1 OPTI1 NAME']]}
                #coros and stops could be in IRDIS (INS1) or in CPI (INS4)
                pass
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

            self.headers=np.r_[self.headers,header]
            #now derotate cube if pupil stabilised
            if self.rot=='PUPIL':
                self.findcentre(self.medianNoDerot[isci][0])
#                pass
            #first find centre of rotation
                
            #then calculate rotation as a function of time
                self.parang=[header['HIERARCH ESO TEL PARANG START'],header['HIERARCH ESO TEL PARANG START']]
                self.pdelt=(self.parang[1]-self.parang[0])/data.shape[0]
                for i in range(data.shape[0]):
                    angle=-1.*self.parang[0]-self.parangInit + (i+0.5)*self.pdelt #rotation angle at centre of exposure relative to beginning of entire sequence - add absolute rotations as well!
                    #then derotate each frame of each half of the detector
                    self.derotate(datal[i,:,.],self.centre,angle) #simage.interpolation.rotate(datal[i,:,:],angle,reshape=False) #somehow I must be able to pass in the centre of rotation...I guess I could also shift it so that it is centred correctly first.
            self.sciframes
            isci+=1


        if self.mode=='CI':
            self.CI()
        elif self.mode=='DPI':
            self.DPI()
        elif self.mode=='SDI':
            self.SDI()
        else:
            print 'IRDIS mode not recognised'

    def readdata(self,filename,**kwargs):
        hdu=fits.open(filename)
        cube=hdu[0].data #extract data itself
        #then extract important header info ? (might have already done this before, not sure about architecture yet)
        if self.masterdark is not None:
            cube=cube-self.masterdark
        if self.masterflat is not None:
            cube=cube/self.masterflat
        header=hdu[0].header
        extra='??'
        hdu.close()
        return cube,header#,extra

class IrdisLSS(Irdis):
    pass
