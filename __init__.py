#Placeholder
#Name stands for "Python REduction of high-Contrast Infrared Spectroscopy and (polarimetric) Imaging ObservatioNs"
"""
    Python Package for reduction of infrared imaging. Currently only treats SPHERE/IRDIS, as that's all I've used to so far, but will be expanded to SPHERE/IFS and NACO before long. 

    Copyright (C) 2018 Peter Scicluna

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

API consists of classes tailored to particular instruments, along with more generic functions that can be combined for use on other instruments. The idea, however, is to extend the base class when new instruments are required.
"""

from . import darkbias
from . import irdis
from . import utility

version_info = (0,2,1) #version, point, revision
__version__ = "v"+'.'.join([str(i) for i in version_info])

__all__=["flatfield", "irdis","darkbias","utility"]

