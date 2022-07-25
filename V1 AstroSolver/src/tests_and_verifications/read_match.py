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


filepath = r'D:\DRDC Data\ORC GBO\Field 1\139\Test\2021_139_DESCENT_46927_B_BESSEL_SA26_SF1_00000356.fit.corr'
# filepath = r'/media/jmwawrow/Data/DRDC Data/2022_07_16_Kingston_Trial_Survey/Light/2022_07_16_LIGHT_05sec/2022_07_16_1x1_5.000secs_0.00C_Light_CLEAR_NoTarget_00002019.fits.corr'

def calculate_zmag(filepath):
    hdr, data = astro.read_fits_file(filepath, extension=1)
    # print(repr(hdr))
    g_mag, bp_mag, rp_mag = read_gaia_mags(data)
    B, V, R, I = convert_gaia_to_bvri(g_mag, bp_mag, rp_mag)
    print(V)
    print(R)
    print(I)

def read_gaia_mags(imgdata):
    g_mag = np.empty(np.shape(imgdata)[0])
    bp_mag = np.empty(np.shape(imgdata)[0])
    rp_mag = np.empty(np.shape(imgdata)[0])
    for i, star in enumerate(imgdata):
        g_mag[i] = star[13]
        bp_mag[i] = star[24]
        rp_mag[i] = star[25]
    return g_mag, bp_mag, rp_mag

def convert_gaia_to_bvri(g_mag, bp_mag, rp_mag):
    B = np.empty(len(g_mag))
    V = np.empty(len(g_mag))
    R = np.empty(len(g_mag))
    I = np.empty(len(g_mag))
    # Calculation from: https://arxiv.org/pdf/1804.09368.pdf
    V = g_mag + 0.01760 + (0.01760 * (bp_mag - rp_mag)) + (0.1732 * ((bp_mag - rp_mag)**2))
    R = g_mag + 0.003226 - (0.3833 * (bp_mag - rp_mag)) + (0.1345 * ((bp_mag - rp_mag)**2))
    I = g_mag - 0.02085 - (0.7419 * (bp_mag - rp_mag)) + (0.09631 * ((bp_mag - rp_mag)**2))
    for i in range(len(g_mag)):
        g_v = g_mag[i] - V[i]
        const = -0.02907 - g_v
        coeff = [-0.001768, -0.2297, -0.02385, const]
        b_v = np.roots(coeff)
        print(b_v)
        # B[i] = b_v + V[i]
    return B, V, R, I

calculate_zmag(filepath)