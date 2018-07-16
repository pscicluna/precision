#Placeholder
#Name stands for "Python REduction of high-Contrast Infrared Spectroscopy and (polarimetric) Imaging ObservatioNs"
"""Python Package for reduction of infrared imaging. Currently only treats SPHERE/IRDIS, as that's all I've used to so far, but will be expanded to SPHERE/IFS and NACO before long. 

API consists of classes tailored to particular instruments, along with more generic functions that can be combined for use on other instruments. The idea, however, is to extend the base class when new instruments are required.
"""

import darkbias
import irdis
import utility

version_info = (0,2,0) #version, point, revision
__version__ = "v"+'.'.join([str(i) for i in version_info])

__all__=["flatfield", "irdis","darkbias","utility"]

