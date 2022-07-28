import os
import sys
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits, ascii
from astropy.table import QTable
from astropy.wcs import WCS
from math import atan
import numpy as np
import warnings
from os.path import dirname
src_path = dirname(dirname(__file__))
sys.path.append(os.path.join(src_path, 'general_tools'))
sys.path.append(os.path.join(src_path, 'Streakdetection'))

import AstroFunctions as astro
import StreakDetection as sd


filepath = r'D:\DRDC Data\ORC GBO\Field 1\139\Test\2021_139_DESCENT_46927_B_BESSEL_SA26_SF1_00000356.fit.corr'
# filepath = r'/media/jmwawrow/Data/DRDC Data/2022_07_16_Kingston_Trial_Survey/Light/2022_07_16_LIGHT_05sec/2022_07_16_1x1_5.000secs_0.00C_Light_CLEAR_NoTarget_00002019.fits.corr'

def calculate_zmag(filepath, save_table=False):
    save_loc = f'{filepath}.csv'
    hdr, data = astro.read_fits_file(filepath, extension=1)
    print(repr(hdr))
    ra, dec, g_mag, bp_mag, rp_mag = read_gaia_mags(data)
    B, V, R, I = convert_gaia_to_bvri(g_mag, bp_mag, rp_mag)
    ref_mag_table = create_ref_table(ra, dec, g_mag, bp_mag, rp_mag, B, V, R, I, 
                                     save_table=save_table, save_loc=save_loc)
    ref_mag_table.pprint_all()

def read_gaia_mags(imgdata):
    ra = np.empty(np.shape(imgdata)[0])
    dec = np.empty(np.shape(imgdata)[0])
    g_mag = np.empty(np.shape(imgdata)[0])
    bp_mag = np.empty(np.shape(imgdata)[0])
    rp_mag = np.empty(np.shape(imgdata)[0])
    for i, star in enumerate(imgdata):
        ra[i] = star[11]
        dec[i] = star[12]
        g_mag[i] = star[13]
        bp_mag[i] = star[24]
        rp_mag[i] = star[25]
    return ra, dec, g_mag, bp_mag, rp_mag

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
        # print(b_v)
        # B[i] = b_v + V[i]
    B[:] = np.nan
    return B, V, R, I


def create_ref_table(ra, dec, g_mag, bp_mag, rp_mag, B, V, R, I, save_table=False, save_loc=None):
    ref_mag_table = QTable(
            names=['RA', 'Dec', 'G', 'BP', 'RP', 'B', 'V', 'R', 'I'],
            data=[ra, dec, g_mag, bp_mag, rp_mag, B, V, R, I],
            units=[u.degree, u.degree, u.mag, u.mag, u.mag, u.mag, u.mag, u.mag, u.mag]
        )
    if save_table:
        if save_loc is not None:
            ascii.write(ref_mag_table, save_loc, 'csv')
        else:
            warnings.warn('save_loc must be provided if save_table is True.')
    return ref_mag_table


os.remove(f'{filepath}.csv')
calculate_zmag(filepath, save_table=True)
