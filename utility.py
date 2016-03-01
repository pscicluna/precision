import numpy as np
import photutils as pu
import glob
import os
import sys
from xml.dom.minidom import parse
import xml.dom.minidom
import astropy.io.fits as fits

def psfsub(source,psfref,centres,rin,rout,**kwargs):
    sourceap=pu.CircularAnnulus(centres[0],rin,rout)
    psfap=pu.CircularAnnulus(centres[1],rin,rout)
    sourcetab=pu.aperture_photometry(source,sourceap)
    psftab=pu.aperture_photometry(psfref,psfap)
    ratio=sourcetab['aperture_sum'][0]/psftab['aperture_sum'][0]
    source=source-psfref*ratio
    return source

def findassociations(dirname,**kwargs):
    '''
    Find ESO file association tree xml files to identify which files should be reduced with which calibrations.
    '''
    files=glob.glob(dirname+'*.xml') #find all xml files in directory of interest.
    files.sort()
    observations=[]
    for f in files:
        if f.startswith('SPHER'):
            #print 'Reduction of SPHERE observations is not implemented yet'
            #print 'Skipping file ',f
            #now make it parse the file
            a=parseassociation(f)
            if a['category'] == 'IRD_SCI_CLI':
                #data can be reduced
                print 'File ',f,' will be reduced'
                observations.append(a)
            else:
                print 'Mode ',a['category'],' is not supported yet'
                print 'skipping file ',f
        elif f.startswith('NACO'):
            print 'Reduction of NACO observations is not implemented yet'
            print 'Skipping file ',f
        elif f.startswith('VISIR'):
            print 'Reduction of VISIR observations is not implemented yet'
            print 'Skipping file ',f
        else:
            print 'Skipping file ',f
    pass

def parseassociation(filename,**kwargs):
    '''
    Read association tree file to build a reduction group/object

    reads in ESO association tree xml and puts the information into a reduction object

    to do:
    figure out which files are shared between different objects, pre-process them and present higher-level calibration files to the main reduction
    '''
#    Handler=xml.sax.handler.ContentHandler
#    association=xml.sax.parse(filename,Handler)
    DOMtree=parse(filename)
    assoctree=DOMtree.firstChild
    if assoctree.attributes['category'].value.startswith('IRD_SCI'):
        mode=assoctree.attributes['category'].value
        science=assoctree.childNodes[1] #mainfiles
        calibs=assoctree.childNodes[5] #associatedfiles
        msgs=assoctree.childNodes[3]
        assocs=CalAssoc(assoctree)
        return assocs
#        pass
    else:
        print 'Currently only IRDIS science data are supported'
        print 'skipping file ',filename
        return
#    for calib in calibs.childNodes[1:-1:2]: #pull out each calibration and identify which file does which thing - some need to be generalised so they don't just refer to CLI
#        if calib.attributes['category'].endswith('DARK'): #dark frames
#            pass
#        elif calib.attributes['category'].endswith('BACK'): #sky?
#            pass
#        elif calib.attributes['category'].endswith('CLI_FLUX'): #flux
#            pass
#        elif calib.attributes['category'].endswith('FLAT'): #flatfields
#            pass
#        elif calib.attributes['category'] == 'IRD_DIST': #distortion maps
#            pass
#        elif calib.attributes['category'].endswith('CLI_PHOT'): #photometric calibration
#            pass
#        else:
#            pass
#        pass
#    scifiles=[]
#    for sci in science.childNodes[1:-1:2]:
#        scifiles.append([sci.attributes['name'].value,sci.attributes['category']])
#    pass
                        

def CalAssoc(node,**kwargs):
    a={'category': node.attributes['category'].value,'main':[],'cals':[]}#,'msg':[],'cals':[]}
    main=node.childNodes[1]
    print main.toxml()
    for m in main.childNodes[1:-1:2]:
        a['main'].append([m.attributes['category'].value,m.attributes['name'].value])
    msg=node.childNodes[3]
    if msg.hasChildNodes():
        message=[]
        for m in msg.childNodes[1:-1:2]:
            message.append(m.attributes['text'].value)
        a['msg']=message
    cals=node.childNodes[5]
    if cals.hasChildNodes():
        for c in cals.childNodes[1:-1:2]:
            a['cals'].append(CalAssoc(c))
    return a

#is this the place for a fits file finder and sorter? I guess I can always move it later
def spherefitssorter(directory,**kwargs):
    #search through a directory and find all fits files. Then search the headers for the OB properties, and produce output which groups files together for reduction
    #for ESO instruments, this can actually be replaced by just parsing the association trees (xml) to figure out which files and calibrations are associated.
    #I can worry about other observatories later
    filelist=glob.glob(directory+'/SPHER*.fits')
    filelist.sort()
    print("There are ",len(filelist)," SPHERE files in this directory")
    observations=[]
    iobs=0
    for f in filelist:
        hdu=fits.open(f)
        if iobs==0: #looks like this thing works! Might be easier/better to build instrument objects directly, though
            observations=[dict(obstemplid=[hdu[0].header['HIERARCH ESO OBS ID'],hdu[0].header['HIERARCH ESO OBS TPLNO'],hdu[0].header['HIERARCH ESO SEQ ARM']],
                               dit=hdu[0].header['HIERARCH ESO DET SEQ1 DIT'],ndit=hdu[0].header['HIERARCH ESO DET NDIT'],
                               nexp=hdu[0].header['HIERARCH ESO TPL NEXP'],filt=hdu[0].header['HIERARCH ESO INS COMB IFLT'],
                               coron=[hdu[0].header['HIERARCH ESO INS COMB ICOR'],hdu[0].header['HIERARCH ESO SEQ CORO XC'],hdu[0].header['HIERARCH ESO SEQ CORO YC'],hdu[0].header['HIERARCH ESO SEQ CORO DATE']],
                               mode=hdu[0].header['HIERARCH ESO INS1 MODE'],rot=hdu[0].header['HIERARCH ESO INS4 COMB ROT'],
                               cat=hdu[0].header['HIERARCH ESO DPR CATG'],obstypes=hdu[0].header['HIERARCH ESO OCS OBSTYPE LIST'],
                               template=hdu[0].header['HIERARCH ESO TPL ID'],
                               pola=hdu[0].header['HIERARCH ESO INS COMB POLA'],reduce=None, #bool to allow querying the user whether to reduce or not
                               filelist=[dict(filename=f,iexp=hdu[0].header['HIERARCH ESO TPL EXPNO'],obstype=hdu[0].header['HIERARCH ESO DPR TYPE'],
                                              parang=[hdu[0].header['HIERARCH ESO TEL PARANG START'],hdu[0].header['HIERARCH ESO TEL PARANG END']],
                                              dither=[hdu[0].header['HIERARCH ESO INS1 DITH POSX'],hdu[0].header['HIERARCH ESO INS1 DITH POSY']])])]  
            iobs+=1
        else:
            for i in range(iobs):
                if ([hdu[0].header['HIERARCH ESO OBS ID'],hdu[0].header['HIERARCH ESO OBS TPLNO'],hdu[0].header['HIERARCH ESO SEQ ARM']] == observations[i]['obstemplid']) & (hdu[0].header['HIERARCH ESO INS COMB IFLT'] == observations[i]['filt']) & (hdu[0].header['HIERARCH ESO INS4 COMB ROT'] == observations[i]['rot']) & (hdu[0].header['HIERARCH ESO DET SEQ1 DIT'] == observations[i]['dit']):
                    observations[i]['filelist'].append(dict(filename=f,iexp=hdu[0].header['HIERARCH ESO TPL EXPNO'],obstype=hdu[0].header['HIERARCH ESO DPR TYPE'],
                                                           parang=[hdu[0].header['HIERARCH ESO TEL PARANG START'],hdu[0].header['HIERARCH ESO TEL PARANG END']],
                                                           dither=[hdu[0].header['HIERARCH ESO INS1 DITH POSX'],hdu[0].header['HIERARCH ESO INS1 DITH POSY']]))
                    break
            else:
                observations.append(
                    dict(obstemplid=[hdu[0].header['HIERARCH ESO OBS ID'],hdu[0].header['HIERARCH ESO OBS TPLNO'],hdu[0].header['HIERARCH ESO SEQ ARM']],
                         dit=hdu[0].header['HIERARCH ESO DET SEQ1 DIT'],ndit=hdu[0].header['HIERARCH ESO DET NDIT'],
                         nexp=hdu[0].header['HIERARCH ESO TPL NEXP'],filt=hdu[0].header['HIERARCH ESO INS COMB IFLT'],
                         coron=[hdu[0].header['HIERARCH ESO INS COMB ICOR'],hdu[0].header['HIERARCH ESO SEQ CORO XC'],hdu[0].header['HIERARCH ESO SEQ CORO YC'],hdu[0].header['HIERARCH ESO SEQ CORO DATE']],
                         mode=hdu[0].header['HIERARCH ESO INS1 MODE'],rot=hdu[0].header['HIERARCH ESO INS4 COMB ROT'],
                         cat=hdu[0].header['HIERARCH ESO DPR CATG'],obstypes=hdu[0].header['HIERARCH ESO OCS OBSTYPE LIST'],
                         template=hdu[0].header['HIERARCH ESO TPL ID'],
                         pola=hdu[0].header['HIERARCH ESO INS COMB POLA'],reduce=None,
                         filelist=[dict(
                                filename=f,iexp=hdu[0].header['HIERARCH ESO TPL EXPNO'],obstype=hdu[0].header['HIERARCH ESO DPR TYPE'],
                                parang=[hdu[0].header['HIERARCH ESO TEL PARANG START'],hdu[0].header['HIERARCH ESO TEL PARANG END']],
                                dither=[hdu[0].header['HIERARCH ESO INS1 DITH POSX'],hdu[0].header['HIERARCH ESO INS1 DITH POSY']]
                                   )]
                         )
                    )
                iobs+=1
        hdu.close
    print len(observations)
    for i in range(iobs):
        print observations[i]
        print ' '
    return observations

        #what date?
#        date=hdu[0].header['DATE']

        #sci, acq, cal?

        #Mode?

        #optics?

        #exposure?

###potentially useful keywords###
#HIERARCH ESO DET SEQ1 DIT    = 0.8374640 / [s] Integration time                 
#HIERARCH ESO DET SEQ1 EXPTIME= 45.5074200 / [s] Exposure Sequence Time 
#HIERARCH ESO DPR CATG = 'SCIENCE ' / Observation category                       
#HIERARCH ESO DPR TECH = 'IMAGE,DUAL,CORONOGRAPHY' / Observation technique       
#HIERARCH ESO DPR TYPE = 'OBJECT  ' / Observation type                           
#HIERARCH ESO INS COMB ICOR = 'N_ALC_YJH_S' / Infra-Red coronograph combination n
#HIERARCH ESO INS COMB IFLT = 'DB_J23  ' / Infra-Red filter combination name.    
#HIERARCH ESO INS COMB POLA = '        ' / Polarisation combination name.    
#HIERARCH ESO INS COMB VCOR = 'V_NC_WF ' / Visible coronograph combination name.  
#HIERARCH ESO INS1 DITH POSX = 0.000 / X Position for dithering                  
#HIERARCH ESO INS1 DITH POSY = 0.000 / Y Position for dithering                  
#HIERARCH ESO INS1 FILT ID = 'FILT_BBF_J' / IRDIS filter unique id.              
#HIERARCH ESO INS1 FILT NAME = 'B_J     ' / IRDIS filter name.                   
#HIERARCH ESO INS1 FILT NO = 14 / IRDIS filter wheel position index.             
#HIERARCH ESO INS1 ID = 'IRDIS/176642' / Instrument ID.                          
#HIERARCH ESO INS1 MODE = 'DBI     ' / Instrument mode used.                     
#HIERARCH ESO INS1 OPTI1 ID = 'MASK_STOP_ALC2' / OPTIi unique ID.                
#HIERARCH ESO INS1 OPTI1 NAME = 'ST_ALC2 ' / IRDIS Lyot stop.                    
#HIERARCH ESO INS1 OPTI1 NO = 13 / OPTIi slot number.                            
#HIERARCH ESO INS1 OPTI1 TYPE = 'MASK_STOP_ALC2' / OPTIi element.                
#HIERARCH ESO INS1 OPTI2 ID = 'FILT_DBF_J23' / OPTIi unique ID.                  
#HIERARCH ESO INS1 OPTI2 NAME = 'D_J23   ' / IRDIS Dual Filter.                  
#HIERARCH ESO INS1 OPTI2 NO = 5 / OPTIi slot number.                             
#HIERARCH ESO INS1 OPTI2 TYPE = 'FILT_DBF_J23' / OPTIi element.  
#HIERARCH ESO INS4 COMB IBS = '        ' / Assembly for infrared beamsplitter    
#HIERARCH ESO INS4 COMB IND = 'ND_N_1.0' / Assembly for infrared neutral density 
#HIERARCH ESO INS4 COMB POLA_CPI = 'N_I     ' / Assembly for polarimetry at CPI l
#HIERARCH ESO INS4 COMB ROT = 'PUPIL   ' / Assembly for derotator and ADC modes  
#HIERARCH ESO INS4 MODE = 'IRDIS_DBI' / Instrument mode used.  
#HIERARCH ESO OBS NAME = 'OB-VYCMa-JH' / OB name                                 
#HIERARCH ESO OBS NTPL =      3 / Number of templates within OB   
#HIERARCH ESO OBS TARG NAME = 'VY CMa  ' / OB target name                        
#HIERARCH ESO OBS TPLNO =     2 / Template number within OB
#HIERARCH ESO OCS OBSTYPE LIST = 'S C O   ' / User defined list of observ. types 
#HIERARCH ESO OCS WAFFLE AMPL = 0.050 / User defined waffle amplitude            
#HIERARCH ESO OCS WAFFLE ORIENT = 'x       ' / User defined waffle orientation  
#HIERARCH ESO SEQ ARM = 'IRDIS   ' / Name of the sub-system.                     
#HIERARCH ESO SEQ CORO CPI_ND = 'ND_1.0  ' / ND filter used in the centring proc.
#HIERARCH ESO SEQ CORO DATE = '2014-12-04T21:06:51' / Date and time of centring. 
#HIERARCH ESO SEQ CORO XC = 476.185 / Coronagraph position at the end of centring
#HIERARCH ESO SEQ CORO YC = 528.171 / Coronagraph position at the end of centring
#HIERARCH ESO TEL AIRM END = 1.058 / Airmass at end                              
#HIERARCH ESO TEL AIRM START = 1.060 / Airmass at start 
#HIERARCH ESO TEL CHOP ST =   F / True when chopping is active   
#HIERARCH ESO TEL PARANG END = -91.075 / [deg] Parallactic angle at end          
#HIERARCH ESO TEL PARANG START = -91.207 / [deg] Parallactic angle at start  
#HIERARCH ESO TPL EXPNO =     3 / Exposure number within template                
#HIERARCH ESO TPL ID = 'SPHERE_irdis_dbi_obs' / Template signature ID            
#HIERARCH ESO TPL NAME = 'SPHERE IRDIS observation DBI submode' / Template name  
#HIERARCH ESO TPL NEXP =     18 / Number of exposures within template 


def irdisobsbuilder(hdu,filename,**kwargs):
    #build a dictionary to contain all properties of an IRDIS observation file
    obs=dict(fitsfile=filename, dit=hdu[0].header['HIERARCH ESO DET SEQ1 DIT'], ndit=hdu[0].header['HIERARCH ESO DET NDIT'],iexp=hdu[0].header['HIERARCH ESO TPL EXPNO'],nexp=hdu[0].header['HIERARCH ESO TPL NEXP'],filt=hdu[0].header['HIERARCH ESO INS COMB IFLT'],coron=[hdu[0].header['HIERARCH ESO INS COMB ICOR'],hdu[0].header['HIERARCH ESO SEQ CORO XC'],hdu[0].header['HIERARCH ESO SEQ CORO YC'],hdu[0].header['HIERARCH ESO SEQ CORO DATE']],parang=[hdu[0].header['HIERARCH ESO TEL PARANG START'],hdu[0].header['HIERARCH ESO TEL PARANG END']],mode=hdu[0].header['HIERARCH ESO INS4 MODE'],rot=hdu[0].header['HIERARCH ESO INS4 COMB ROT'],dither=[hdu[0].header['HIERARCH ESO INS1 DITH POSX'],hdu[0].header['HIERARCH ESO INS1 DITH POSY']],cat=hdu[0].header['HIERARCH ESO DPR CATG'],obstype=[hdu[0].header['HIERARCH ESO DPR TYPE'],hdu[0].header['HIERARCH ESO OCS OBSTYPE LIST']],template=hdu[0].header['HIERARCH ESO TPL ID'],arm=hdu[0].header['HIERARCH ESO SEQ ARM'],obsid=hdu[0].header['HIERARCH ESO OBS ID'],origfile=hdu[0].header['ORIGFILE'])
    return obs

def drizzle():
    pass

def jitter():
    pass

def measuresensitivity():
    pass

def measurecontrast():
    pass

def conversionfac():
    pass

