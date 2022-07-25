import os
import sys
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from math import atan
import numpy as np
from os.path import dirname
src_path = dirname(dirname(__file__))
sys.path.append(os.path.join(src_path, 'general_tools'))
sys.path.append(os.path.join(src_path, 'Streakdetection'))

import AstroFunctions as astro
import StreakDetection as sd


filepath = r'/media/jmwawrow/Data/DRDC Data/ORC GBO/Field 1/139/Test/2021_139_DESCENT_46927_B_BESSEL_SA26_SF1_00000356.fit.corr'
# filepath = r'/media/jmwawrow/Data/DRDC Data/2022_07_16_Kingston_Trial_Survey/Light/2022_07_16_LIGHT_05sec/2022_07_16_1x1_5.000secs_0.00C_Light_CLEAR_NoTarget_00002019.fits.corr'
hdul = fits.open(filepath)
hdr = hdul[1].header
data = hdul[1].data
print(repr(hdr))
print(data[0][13])