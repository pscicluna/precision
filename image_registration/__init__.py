######################################################################
#
#        This file is from Adam Ginsburg's image_registration 
#        library (originally released un the MIT license) and
#        is used to provide precise image registration. At 
#        some point, I may replace it with a more elegant 
#        solution than simply copying the files into the repo
#        but we're not there yet. Any copyright pertaining
#        to this file remains with Adam Ginsburg.
#
#    P.S. 2016 01 14
#
######################################################################

from cross_correlation_shifts import cross_correlation_shifts, cross_correlation_shifts_FITS
from chi2_shifts import chi2_shift,chi2n_map,chi2_shift_iterzoom
from register_images import register_images
import fft_tools
import tests
#import FITS_tools
from .version import __version__
