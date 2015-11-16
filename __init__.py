#Placeholder
#Name stands for "Python REduction of high-Contrast Infrared Spectroscopy and (polarimetric) Imaging ObservatioNs"
"""Python Package for reduction of infrared imaging. Currently only treats SPHERE/IRDIS, as that's all I've used to so far, but will be expanded to SPHERE/IFS and NACO before long. 

API consists of classes tailored to particular instruments, along with more generic functions that can be combined for use on other instruments. The idea, however, is to extend the base class when new instruments are required.
"""

import darkbias
import irdis
import utility

__all__=["instrument", "irdis","darkbias","utility"]

