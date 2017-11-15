import numpy as np
import scipy.ndimage as simage
import astropy.io.fits as fits
import astropy.stats as astats
import math as mt
from instrument import Instrument
#import astropysics.ccd as ccd
import photutils as pu
#import skimage as skim
#from skimage.transform import hough_circle
#from skimage.util import img_as_ubyte
import darkbias
import flatfield
from image_registration.register_images import register_images
#from skimage.feature import peak_local_max, canny
import matplotlib.pyplot as plt
import matplotlib.colors as cols
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
import skimage.transform as tf

#need a list of filters, conoragraphs, etc
class Irdis(Instrument):
    def __init__(self, observation,#,scifiles=None,badpixfile=None,skyfiles=None,
                 #flatfiles=None,darkfiles=None,mode=None,optics=None,
                 options=None, interProd=True, **kwargs):
        self.science=observation.obs_main
        self.datadir=observation.datadir
        self.mode=observation.obs_type
        self.darks=observation.cals['IRD_DARK']
        if 'CLI' in self.mode:
            self.flats=observation.cals['IRD_CLI_FLAT']
            try:
                self.sky=observation.cals['IRD_SCI_CLI_SKY']  #or does this need to be 'IRD_SKY_BG_RAW'? #Sky
            except KeyError:
                self.sky=observation.cals['IRD_SKY_BG_RAW']  #or does this need to be 'IRD_SKY_BG_RAW'? #Sky
        try:
            self.flux=observation.cals['IRD_FLUX_CALIB_CORO_RAW'] #flux
        except KeyError:
            self.flux=None #flux observations were not taken (e.g. this may be a reference psf)
#        self.scifiles=scifiles
#        self.badpixfile=badpixfile
#        self.skyfiles=skyfiles
#        self.flatfiles=flatfiles
#        self.flatcals
#        self.darkfiles=darkfiles
#        self.mode=mode 
#        self.optics=optics
        self.options=options
        self.interProd=interProd #write out intermediate products? (master dark, etc)
        #calibration angles
        self.pupilOffset=135.87 #Pupil offset angle
        self.tn=-1.764 #True North offset
        #astrometric calibration numbers - probably won't use these and just say that pipeline is not suitable for high-precision astrometry
        self.plateScale=12.251 #IRDIS plate scale (after distortion correction) of irdis detector
        self.Anamorph=1.0062 #anamophic distortion multiplier for Y direction of detector

    def __repr__(self):
        return "<IRDIS mode: %s science data: %s sky: %s flux: %s flats: %s darks: %s options: %s>" % (self.mode, self.science, self.sky, self.flux, self.flats, self.darks)

    def __str__(self):
        return "IRDIS mode: %s science data: %s sky: %s flux: %s flats: %s darks: %s options: %s" % (self.mode, self.science, self.sky, self.flux, self.flats, self.darks, self.options)

    def findSources(self,sourcefinder,data,window,guess,median,interact=True,**kwargs):
        sources=sourcefinder.find_stars(data[np.int(guess[0])-window:np.int(np.ceil(guess[0]))+window,#daofind(data[np.int(guess[0])-window:np.int(np.ceil(guess[0]))+window,
                                              np.int(guess[1])-window:np.int(np.ceil(guess[1]))+window] - median)
        sources.sort('flux')
        if interact:
#            print sources
            sources.pprint(max_lines=-1)
            plt.imshow(data[np.int(guess[0])-window:np.int(np.ceil(guess[0]))+window,
                            np.int(guess[1])-window:np.int(np.ceil(guess[1]))+window] - median,
                       origin='lower')
            plt.show(block=False)
            #try:
            a=raw_input('What is the ID of the source at the centre of rotation?\n')
            plt.close()
            #except ValueError:
            if int(a) < 0:
                #plt.close()
                sources,new_window=self.findSources(sourcefinder,data,int(1.5*window),guess,median)
                #sources.pprint(max_lines=-1)
                #plt.imshow(data[np.int(guess[0])-window:np.int(np.ceil(guess[0]))+window,
                #                np.int(guess[1])-window:np.int(np.ceil(guess[1]))+window] - median)
                #plt.show(block=False)
            else:
                b=(sources['id'] == int(a))
                return sources[:][b],window
        return sources,window#sources
    
    def findcentre(self,data,guess,interact=True,window=60,**kwargs):
        print data.shape,guess
        print np.int(guess[0])
        print np.int(np.ceil(guess[1]))
        #help(np.int(np.ceil(guess[1])) + window)
        #help(data)
        mean,median,std=astats.sigma_clipped_stats(data,sigma=3.0)
        sourcefinders=pu.DAOStarFinder( #data[500:550,450:500] - median, # assumes datacube #[500:550,450:500]
                           fwhm = 2.0,
                           threshold=std#,
                           #sharplo=0.3,
                           #sharphi=0.5,
                           #roundhi=0.3,
                           #roundlo=-0.3
                           )
        #sources
        #b=
        self.central,new_window=self.findSources(sourcefinders,data,window,guess,median,interact)
#        sources.sort('flux')#'peak') #take brightest source
        #try making this interactive?
        #if interact:
#            print sources
        #    sources.pprint(max_lines=-1)
        #    plt.imshow(data[np.int(guess[0])-window:np.int(np.ceil(guess[0]))+window,
        #                    np.int(guess[1])-window:np.int(np.ceil(guess[1]))+window] - median)
        #    plt.show(block=False)
            #try:
        #    a=raw_input('What is the ID of the source at the centre of rotation?\n')
            #except ValueError:
        #    if int(a) < 0:
        #        sources=self.findSources(sourcefinders,data,int(1.5*window),guess,median)
        #        sources.pprint(max_lines=-1)
        #        plt.imshow(data[np.int(guess[0])-window:np.int(np.ceil(guess[0]))+window,
        #                        np.int(guess[1])-window:np.int(np.ceil(guess[1]))+window] - median)
        #        plt.show(block=False)
                #try:
        #        a=raw_input('What is the ID of the source at the centre of rotation?\n')
                #            print a,sources.dtype
        #    b=(sources['id'] == int(a))
#            print b
#            print sources[b]
            #try:
                #self.central.append(sources[:][b])
            #except:
        #self.central=sources[:][b]
        #    plt.close()
        #try:
        #    self.central=sources[:][-1]
        #except:
        #    self.central=sources
#        self.central=central
        #self.sources=sources
        #return central,sources
        #print 
        self.centre=np.array([self.central['xcentroid'].data+(np.int(guess[1])-new_window),self.central['ycentroid'].data+(np.int(guess[0])-new_window)])
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

    def derotate(self,data,centre,angle,**kwargs):
        #pad image so that centre is at the centre of the array
        padX=[data.shape[1] - centre[0],centre[0]] #check which is row and column in image and centre routines! - this doesn't work for non-integer cases. Have to do some other shifting I guess
        padY=[data.shape[0] - centre[1],centre[1]]
        print padX,padY
        datap=np.pad(data,[padY,padX],'constant')
        #derotate
        data=simage.interpolation.rotate(datap,angle,reshape=False)[padY[0]:-padY[1],padX[0]-padX[1]]

    def CI(self,**kwargs): #test on VY CMa
        self.finalNoDerot=np.mean(self.medianNoDerot,axis=0) #produces L- and R- channel images
        #self.finalNoDerot=np.mean(self.finalNoDerot,axis=0) #combine both channels
        #self.finalvarNoDerot=np.mean(self.varNoDerot,axis=0)         #small number approximation - 
        self.finalvarNoDerot=(np.mean(self.varNoDerot,axis=0) / #variance on mean = mean of variances / N
                              (self.varNoDerot.shape[0]))#(2*self.varNoDerot.shape[0]))          #or alternatively
                                                                     #uncertainty on mean = mean of uncertainties / sqrt(N)
        #need to fix uncertainty tracking...current version should underestimate variances, only alteratives (apart from MC) appear to overestimate them
#        self.finalvarNoDerot+=
        #if self.rot=='PUPIL':        
        #    self.finalDerot=np.mean(self.medianDerot,axis=0) #produces L- and R- channel images
        #    self.finalDerot=np.mean(self.finalDerot,axis=0) #combine both channels
        #    self.finalVarDerot=np.mean(self.varDerot,axis=0)         #small number approximation - 
        #    self.finalVarDerot=(np.mean(self.finalVarDerot,axis=0) / #variance on mean = mean of variances / N
        #                          (2*self.varDerot.shape[0]))          #or alternatively
        return 
#        pass

    def SDI(self,**kwargs): #test on GD50, VY CMa
        self.CI()

        self.SDINoDerot=np.mean(self.medianNoDerot,axis=0) #produces L- and R- channel images
        self.SDINoDerot=self.SDINoDerot[1]-self.SDINoDerot[0] #check which way round this is supposed to be! and scaling!

        self.varSDINoDerot=np.mean(self.varNoDerot,axis=0) / self.varNoDerot.shape[0]
        self.varSDINoDerot=np.sum(self.varSDINoDerot, axis=0)
        
        if self.rot=='PUPIL':
            self.SDIDerot=np.mean(self.medianDerot,axis=0)
            self.SDIDerot=self.SDIDerot[1]-self.SDIDerot[0]
            self.varSDIDerot=np.mean(self.varDerot,axis=0) / self.varDerot.shape[0]
            self.varSDIDerot=np.sum(self.varSDIDerot,axis=0)
        return

    def ADI(self,**kwargs): #test on GD50
        #classical ADI, take median non-derotated frame and subtract it from each non-derotated frame, then derotate all frames and collapse the whole cube
        self.medianADI=np.array([])
        self.varADI=np.array([])
        isci=0
        for f in self.scifiles:
            #read data
            data,header=self.readdata(f)
            data=self.badpixcorrect(data)
            data=self.skysub(data,self.mastersky)
            #undither data
            data=simage.interpolation.shift(data,[header['HIERARCH ESO INS1 DITH POSX'], #check sign
                                                  header['HIERARCH ESO INS1 DITH POSY'] #of shifts
                                                  ]
                                            )
            #split channels and align
            datal,datar=self.splitchannels(data)
            chanshift=register_images(datal[0],datar[0],usamp=10.)
            datar=simage.interpolation.shift(datar,[0,chanshift[0],chanshift[1]])
            #subtract speckle image (non-derotated median frames)
            datal=datal-self.finalNoDerot
            datar=datar-self.finalNoDerot
            #derotate channels
            centreGuess=[header['HIERARCH ESO SEC CORO XC'],header['HIERARCH ESO SEC CORO XC']]
            self.findcentre(self.medianNoDerot[isci][0],centreGuess)
#            self.centre=[self.central['xcentroid'],self.central['ycentroid']]                
            #then calculate rotation as a function of time
            self.parang=[header['HIERARCH ESO TEL PARANG START'],header['HIERARCH ESO TEL PARANG START']]
            self.pdelt=(self.parang[1]-self.parang[0])/data.shape[0]
            for i in range(data.shape[0]):
                angle=-1.*self.parang[0]-self.parangInit + (i+0.5)*self.pdelt #rotation angle at centre of exposure relative to beginning of entire sequence - add absolute rotations as well!
                    #then derotate each frame of each half of the detector
                self.derotate(datal[i,:,:],self.centre,angle) #simage.interpolation.rotate(datal[i,:,:],angle,reshape=False) #somehow I must be able to pass in the centre of rotation...I guess I could also shift it so that it is centred correctly first.
                
            centreGuess=[header['HIERARCH ESO SEC CORO XC'],header['HIERARCH ESO SEC CORO XC']]
            self.findcentre(self.medianNoDerot[isci][1],centreGuess)
            self.centre=[self.central['xcentroid'],self.central['ycentroid']]
            #then calculate rotation as a function of time
            self.parang=[header['HIERARCH ESO TEL PARANG START'],header['HIERARCH ESO TEL PARANG START']]
            self.pdelt=(self.parang[1]-self.parang[0])/data.shape[0]
            for i in range(data.shape[0]):
                angle=-1.*self.parang[0]-self.parangInit + (i+0.5)*self.pdelt #rotation angle at centre of exposure relative to beginning of entire sequence - add absolute rotations as well!
                    #then derotate each frame of each half of the detector
                self.derotate(datar[i,:,:],self.centre,angle)
            #now build ADI medians for each dither position
            self.medianADI=np.r_[self.medianADI,
                                 [np.nanmedian(datal,axis=0),np.nanmedian(datar,axis=0)]
                                 ]
            self.varADI=np.r_[self.varADI,
                              [(np.pi/(2.*datal.shape[0]))*(np.std(datal,axis=0))**2,
                               (np.pi/(2.*datar.shape[0]))*(np.std(datar,axis=0))**2
                               ]
                              ]
            isci+=1
        #average all frames ('cADI') 
        self.finalADI=np.mean(self.medianADI,axis=0) #produces L- and R- channel images
        self.finalADI=np.mean(self.finalADI,axis=0) #combine both channels
        self.finalVarADI=np.mean(self.varADI,axis=0)         #small number approximation - 
        self.finalVarADI=(np.mean(self.finalVarADI,axis=0) / #variance on mean = mean of variances / N
                            (2*self.varDerot.shape[0]))

        if mode=='SDI':
            self.finalSADI=np.mean(self.medianSADI,axis=0)    #produces L- and R- channel images
            self.finalSADI=self.finalSADI[1]-self.finalSADI[0]#
            self.finalVarSADI=np.mean(self.varSADI,axis=0)    #small number approximation - 
            self.finalVarSADI=np.sum(self.finalVarSADI,axis=0)#variance on mean = mean of variances / N
        pass
    
    def LOCI(self,**kwargs):
        raise NotImplementedError("LOCI is not implemented yet.")#        pass

    def DPI(self,**kwargs): #test on HR 3090
        self.CI()

        #split L and R channels into O- and E-rays respectively
        self.MeanOray=np.mean(self.medianNoDerot[:,0,:,:],axis=0) #add an additional axis for DPI observations, so that it contains [+Q, -Q, +U, -U] or as many as are available
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

    def makesky(self,**kwargs):
        print self.sky
        if len(self.sky.obs_main) == 1:
            f=self.sky.datadir+self.sky.obs_main[0][:][1]+'.fits'
            data,header=self.readdata(f)
            data=self.badpixcorrect(data)
            self.skyframes=data
        else:
            for f in self.sky.obs_main[0][:][1]:
                f=self.sky.datadir+f+'.fits'
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
        pass

    def reduce(self,**kwargs):
        if self.darks is not None:
            self.masterdark,self.darkvar,self.RON=darkbias.makemasterdark(self.darks,**kwargs)
        else:
            self.masterdark=None

        if self.flats is not None:        
            self.masterflat,self.flatvar,self.badpixmap=flatfield.makemasterflat(self.flats,**kwargs)#,self.masterdark,**kwargs)
        else:
            self.masterflat=None
        if self.sky is not None:
            self.makesky()

        print np.sum(np.isnan(self.mastersky))
        self.medianNoDerot=np.array([])
        self.medianDerot=np.array([])
        self.varNoDerot=np.array([])
        self.varDerot=np.array([])
        self.headers=np.array([])
        isci=-1
        print self.science[0][1]
        #exit()
        self.centrel=[]
        self.centrer=[]
        self.shiftl=[]
        self.shiftr=[]
        for f in self.science:
            isci+=1
            #if isci == 0:
            #    continue
            f=self.datadir+f[1]+'.fits'
            #read data
            data,header=self.readdata(f)
            #print data
            #print np.sum(np.isnan(data))
            if isci==0:
                #pull important info out of header from first science file
                self.rot=header['HIERARCH ESO INS4 COMB ROT']
                self.parangInit=header['HIERARCH ESO TEL PARANG START']
                self.optics={'filt': [header['HIERARCH ESO INS1 FILT NO'],header['HIERARCH ESO INS1 FILT ID'],header['HIERARCH ESO INS1 FILT NAME']],'opti':[header['HIERARCH ESO INS1 OPTI2 NO'],header['HIERARCH ESO INS1 OPTI2 ID'],header['HIERARCH ESO INS1 OPTI2 NAME']],'stop': [header['HIERARCH ESO INS1 OPTI1 NO'],header['HIERARCH ESO INS1 OPTI1 ID'],header['HIERARCH ESO INS1 OPTI1 NAME']]}
                #coros and stops could be in IRDIS (INS1) or in CPI (INS4)
                pass
            #intermediate processing
            data=self.badpixcorrect(data)
            #print data
            #print np.sum(np.isfinite(data))
            data=self.skysub(data,self.mastersky)
            #print data
            #print np.sum(np.isfinite(data))
            mask=np.logical_not(np.isfinite(data))
            data[mask]=np.nan
            if (len(data.shape) == 3): #datacube
                for i in range(data.shape[0]):
                    data[i,:,:] = interpolate_replace_nans(data[i,:,:],Gaussian2DKernel(stddev=1))
            else:
                data[:,:] = interpolate_replace_nans(data[:,:],Gaussian2DKernel(stddev=1))
            #split channels and align
            #print data[np.logical_not(np.isfinite(data))]
            #mask=np.logical_not(np.isfinite(data))
            
            #exit()
            datal,datar=self.splitchannels(data)
            #print np.sum(np.isfinite(datal))
            #print np.sum(np.isfinite(datar))
            
#            print data.shape,datal.shape
#            exit()
            #chanshift=register_images(datal[0],datar[0],usfac=10.)
            #datar=simage.interpolation.shift(datar,[0,chanshift[0],chanshift[1]])
            #data=np.r_[datal,datar]
            #datal=None
            #datar=None

            #centreGuess=[header['HIERARCH ESO SEC CORO XC'],header['HIERARCH ESO SEC CORO XC']]
            #self.findcentre(self.medianNoDerot[isci][0],centreGuess)
            #self.centre=[self.central['xcentroid'],self.central['ycentroid']]

            #----------------------------NOW THINGS ARE DIFFERENT DEPENDING ON METHOD!!--------------------------
            try:
                self.medianNoDerot=np.r_[self.medianNoDerot,
                                         [simage.interpolation.shift(np.nanmedian(datal,axis=0),
                                                                     #[0,
                                                                     [header['HIERARCH ESO INS1 DITH POSX'],
                                                                      header['HIERARCH ESO INS1 DITH POSY']
                                                                     ]
                                         ),
                                          simage.interpolation.shift(np.nanmedian(datar,axis=0),
                                                                     #[0,
                                                                     [header['HIERARCH ESO INS1 DITH POSX'],
                                                                      header['HIERARCH ESO INS1 DITH POSY']])
                                         ]
                                        ] #make sure the sign is right here (by inspection!!)
                self.varNoDerot=np.r_[self.varNoDerot,
                                      [simage.interpolation.shift((np.pi/(2.*datal.shape[0]))*(np.std(datal,axis=0))**2,
                                                                  #[0,
                                                                  [header['HIERARCH ESO INS1 DITH POSX'],
                                                                   header['HIERARCH ESO INS1 DITH POSY']]),
                                       simage.interpolation.shift((np.pi/(2.*datar.shape[0]))*(np.std(datar,axis=0))**2,
                                                                  #[0,
                                                                  [header['HIERARCH ESO INS1 DITH POSX'],
                                                                   header['HIERARCH ESO INS1 DITH POSY']])
                                      ]
                                     ]
            except:
                self.medianNoDerot=np.array(
                                         [simage.interpolation.shift(np.nanmedian(datal,axis=0),
                                                                     #[0,
                                                                     [header['HIERARCH ESO INS1 DITH POSX'],
                                                                      header['HIERARCH ESO INS1 DITH POSY']
                                                                     ]
                                         ),
                                          simage.interpolation.shift(np.nanmedian(datar,axis=0),
                                                                     #[0,
                                                                     [header['HIERARCH ESO INS1 DITH POSX'],
                                                                      header['HIERARCH ESO INS1 DITH POSY']])
                                         ]
                                        ) #make sure the sign is right here (by inspection!!)
                self.varNoDerot=np.array(
                                      [simage.interpolation.shift((np.pi/(2.*datal.shape[0]))*(np.std(datal,axis=0))**2,
                                                                  #[0,
                                                                   [header['HIERARCH ESO INS1 DITH POSX'],
                                                                   header['HIERARCH ESO INS1 DITH POSY']]),
                                       simage.interpolation.shift((np.pi/(2.*datar.shape[0]))*(np.std(datar,axis=0))**2,
                                                                 # [0,
                                                                  [header['HIERARCH ESO INS1 DITH POSX'],
                                                                   header['HIERARCH ESO INS1 DITH POSY']])
                                      ]
                                     )

            #plt.imshow(self.medianNoDerot[0,:,:])
            print np.max(self.medianNoDerot)
            #plt.show()
            #exit()
            self.headers=np.r_[self.headers,header]
            #print header
            print self.medianNoDerot.shape
            #if isci == 0:
            centreGuess=[header['HIERARCH ESO SEQ CORO XC'],header['HIERARCH ESO SEQ CORO YC']]
            self.findcentre(self.medianNoDerot[isci*2],centreGuess,interact=True)
            self.centrel.append(self.centre)#[self.central['xcentroid'],self.central['ycentroid']]
            #centreGuess.append([header['HIERARCH ESO SEQ CORO XC'],header['HIERARCH ESO SEQ CORO YC']])
            self.findcentre(self.medianNoDerot[isci*2+1],centreGuess,interact=True)
            self.centrer.append(self.centre)#[self.central['xcentroid'],self.central['ycentroid']]
            print self.centrel[isci],self.centrer[isci]
            self.shiftl.append(np.array(511 - self.centrel[isci]))
            self.shiftr.append(np.array(511 - self.centrer[isci]))
            print self.shiftl[isci],self.shiftl[isci][0],self.shiftl[isci][1][0]
            #tform=tf.SimilarityTransform(scale=1,rotation=0,translation=(self.shiftl[isci][1],self.shiftl[isci][0]))
            self.medianNoDerot[isci*2]=simage.interpolation.shift(self.medianNoDerot[isci*2],np.array([self.shiftl[isci][1][0],self.shiftl[isci][0][0]]))#tf.warp(self.medianNoDerot[isci*2 -1],tform)# [yshift, xshift]
            self.medianNoDerot[isci*2+1]=simage.interpolation.shift(self.medianNoDerot[isci*2 + 1],np.array([self.shiftr[isci][1][0],self.shiftr[isci][0][0]]))
            self.varNoDerot[isci*2]=simage.interpolation.shift(self.varNoDerot[isci*2],np.array([self.shiftl[isci][1][0],self.shiftl[isci][0][0]]))#tf.warp(self.medianNoDerot[isci*2 -1],tform)# [yshift, xshift]
            self.varNoDerot[isci*2+1]=simage.interpolation.shift(self.varNoDerot[isci*2 + 1],np.array([self.shiftr[isci][1][0],self.shiftr[isci][0][0]]))
            #self.centre=self.centrel
            #else:
            #self.medianNoDerot[isci*2]=simage.interpolation.shift(self.medianNoDerot[isci*2],[self.shiftr[1],self.shiftr[0]])
            #print simage.interpolation.shift(self.medianNoDerot[1],self.shiftr)
            #print self.medianNoDerot[1]
            #print self.shiftr
            print isci*2,isci*2+1
            temp=self.medianNoDerot[isci*2]+self.medianNoDerot[isci*2+1]#simage.interpolation.shift(self.medianNoDerot[1],[self.shiftr[1],self.shiftr[0]])
            #temp-=np.median(temp)
            #print np.sum(np.isnan(temp))
            #print np.sum(np.logical_not(np.isfinite(temp)))
            #print np.max(temp),np.min(temp)
#            print temp.shape
            #self.findcentre(temp,[512,512],window=512,interact=True)
#            plt.imshow(temp, norm=cols.Normalize(vmin=0.1,vmax=3.),origin='lower')
#            plt.show()
#            exit()
            #now derotate cube if pupil stabilised
            if self.rot=='PUPIL':
                continue
            #    centreGuess=[header['HIERARCH ESO SEQ CORO XC'],header['HIERARCH ESO SEQ CORO XC']]
            #    self.findcentre(self.medianNoDerot[isci][0],centreGuess)
            #    self.centre=[self.central['xcentroid'],self.central['ycentroid']]
#                pass
            #first find centre of rotation
                
            #then calculate rotation as a function of time
                self.parang=[header['HIERARCH ESO TEL PARANG START'],header['HIERARCH ESO TEL PARANG START']]
                self.pdelt=(self.parang[1]-self.parang[0])/data.shape[0]
                #shift data to common reference frame
                for i in range(data.shape[0]):
                    angle=-1.*self.parang[0]-self.parangInit + (i+0.5)*self.pdelt #rotation angle at centre of exposure relative to beginning of entire sequence - add absolute rotations as well!
                    #then derotate each frame of each half of the detector
                    self.derotate(datal[i,:,:],self.centre,angle) #simage.interpolation.rotate(datal[i,:,:],angle,reshape=False) #somehow I must be able to pass in the centre of rotation...I guess I could also shift it so that it is centred correctly first.
                
                #centreGuess=[header['HIERARCH ESO SEC CORO XC'],header['HIERARCH ESO SEC CORO XC']]
                #self.findcentre(self.medianNoDerot[isci][1],centreGuess)
                #self.centre=[self.central['xcentroid'],self.central['ycentroid']]
#                pass
            #first find centre of rotation
                
            #then calculate rotation as a function of time
                self.parang=[header['HIERARCH ESO TEL PARANG START'],header['HIERARCH ESO TEL PARANG START']]
                self.pdelt=(self.parang[1]-self.parang[0])/data.shape[0]
                for i in range(data.shape[0]):
                    angle=-1.*self.parang[0]-self.parangInit + (i+0.5)*self.pdelt #rotation angle at centre of exposure relative to beginning of entire sequence - add absolute rotations as well!
                    #then derotate each frame of each half of the detector
                    self.derotate(datar[i,:,:],self.centre,angle) 
                self.medianDerot=np.r_[self.medianDerot,
                                       [simage.interpolation.shift(np.nanmedian(datal,axis=0),
                                                                   [0,header['HIERARCH ESO INS1 DITH POSX'],
                                                                    header['HIERARCH ESO INS1 DITH POSY']]),
                                        simage.interpolation.shift(np.nanmedian(datar,axis=0),
                                                                   [0,header['HIERARCH ESO INS1 DITH POSX']+self.shiftr[1],
                                                                    header['HIERARCH ESO INS1 DITH POSY']+self.shiftr[0]])
                                        ]
                                       ] #make sure the sign is right here (by inspection!!)
                self.varDerot=np.r_[self.varDerot,
                                    [simage.interpolation.shift((np.pi/(2.*datal.shape[0]))*(np.std(datal,axis=0))**2,
                                                                [0,header['HIERARCH ESO INS1 DITH POSX'],
                                                                 header['HIERARCH ESO INS1 DITH POSY']]),
                                     simage.interpolation.shift((np.pi/(2.*datar.shape[0]))*(np.std(datar,axis=0))**2,
                                                                [0,header['HIERARCH ESO INS1 DITH POSX']+self.shiftr[1],
                                                                 header['HIERARCH ESO INS1 DITH POSY']+self.shiftr[0]])
                                     ]
                                    ]
            #self.sciframes
            #isci+=1


        if self.mode=='IRD_SCI_CLI_OBJ':#'CI':
            self.CI()
            print  self.finalNoDerot.shape
#            plt.imshow(self.finalNoDerot,norm=cols.Normalize(vmin=0.1,vmax=1.))
#            plt.show()
        elif self.mode=='DPI':
            self.DPI()
        elif self.mode=='SDI':
            self.SDI()
        else:
            print 'IRDIS mode not recognised'
            
    def sciFlux(self,**kwargs):
        ''' 
        Take non-coronographic observation to compute contrast and facilitate flux calibration
        '''
        #read flux file and calibrate it
        fluxframes,fluxhdr=self.readdata(self.flux)
        #extract point-source counts and peak counts for target
        
        #scale for ND filters and exposure time
        pass

    def photCal(self,**kwargs):
        '''
        Perform photometric calibration with standard star
        '''
        pass

    def readdata(self,filename,**kwargs):
        hdu=fits.open(filename)
        cube=hdu[0].data #extract data itself
        exptime=hdu[0].header['EXPTIME']
        #then extract important header info ? (might have already done this before, not sure about architecture yet)
        cube=cube/exptime #?
        if self.masterdark is not None:
            cube=cube-self.masterdark
        if self.masterflat is not None:
            cube=cube/self.masterflat
        header=hdu[0].header
        extra='??'
        hdu.close()
        return cube,header#,extra

    def output(self,**kwargs):
        self.outfile=(self.headers[0]['HIERARCH ESO OBS NAME'] + '_' + 
                      self.optics['filt'][2] + '_' + 
                      self.headers[0]['HIERARCH ESO DET SEQ1 DIT'] + '_' + 
                      self.headers[0]['HIERARCH ESO OBS START'] +
                      '.fits')
        print 'Writing reduced data to file ',self.outfile
        if self.interProd: #write intermediate products too
            pass
        fits.writeto(self.outfile ,self.finalNoDerot,self.outheader)
        fits.update(self.outfile , self.finalvarNoDerot,'var CI')
        if self.rot=='PUPIL':
            fits.update(self.outfile , self.finalDerot, 'CI derotated')
            fits.update(self.outfile , self.finalvarDerot, 'var CI derotated')
            fits.update(self.outfile , self.finalADI, 'Classical ADI')
            fits.update(self.outfile , self.finalvarADI, 'var Classical ADI')
        if self.mode=='SDI':
            fits.update(self.outfile , self.SDINoDerot, 'SDI')
            fits.update(self.outfile , self.varSDINoDerot, 'var SDI')
            if self.rot=='PUPIL':
                fits.update(self.outfile , self.SDIDerot, 'SDI derotated')
                fits.update(self.outfile , self.varSDIDerot, 'var SDI derotated')
                fits.update(self.outfile , self.finalSADI, 'SADI')
                fits.update(self.outfile , self.finalvarSADI, 'var SADI')
        if self.mode=='DPI':
            pass
#        fits.update( , self.badpixmap, 'bad pixel map')
        pass

class IrdisLSS(Irdis):
    pass
