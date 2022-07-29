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


filepath = r'/media/jmwawrow/Data/DRDC Data/ORC GBO/Field 1/139/Test/2021_139_DESCENT_46927_B_BESSEL_SA38_00000406.fit'
# filepath = r'/media/jmwawrow/Data/DRDC Data/2022_07_16_Kingston_Trial_Survey/Light/2022_07_16_LIGHT_05sec/2022_07_16_1x1_5.000secs_0.00C_Light_CLEAR_NoTarget_00002019.fits.corr'

def calculate_zmag(filepath, save_table=False, exposure_key='EXPTIME'):
    save_loc = f'{filepath}.csv'
    hdr, _ = astro.read_fits_file(filepath)
    _, data = astro.read_fits_file(f'{filepath}.corr', extension=1)
    # print(repr(hdr))
    field_id, ra, dec, g_mag, bp_mag, rp_mag = read_gaia_mags(data)
    B, V, R, I = convert_gaia_to_bvri(g_mag, bp_mag, rp_mag)
    ref_mag_table = create_ref_table(field_id, ra, dec, g_mag, bp_mag, rp_mag, B, V, R, I, 
                                     save_table=save_table, save_loc=save_loc)
    # ref_mag_table.pprint_all()
    instr_mags = read_field_mags(filepath, hdr, exposure_key)
    # print(instr_mags)
    mag_diff = calculate_mag_diff(hdr, ref_mag_table, field_id, instr_mags)
    # print(mag_diff)
    zmag, zmag_sigma = average_mag_diff(mag_diff)
    print(f'Pinpoint zmag: {hdr["zmag"]}')
    print(f'Zmag: {zmag:.3f} +/- {zmag_sigma:.3f}')

def read_gaia_mags(imgdata):
    field_id = np.empty(np.shape(imgdata)[0], dtype=int)
    ra = np.empty(np.shape(imgdata)[0])
    dec = np.empty(np.shape(imgdata)[0])
    g_mag = np.empty(np.shape(imgdata)[0])
    bp_mag = np.empty(np.shape(imgdata)[0])
    rp_mag = np.empty(np.shape(imgdata)[0])
    for i, star in enumerate(imgdata):
        field_id[i] = star[9]
        ra[i] = star[11]
        dec[i] = star[12]
        g_mag[i] = star[13]
        bp_mag[i] = star[24]
        rp_mag[i] = star[25]
    return field_id, ra, dec, g_mag, bp_mag, rp_mag

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


def create_ref_table(field_id, ra, dec, g_mag, bp_mag, rp_mag, B, V, R, I, save_table=False, save_loc=None):
    sc = SkyCoord(ra, dec, unit='deg')
    ref_mag_table = QTable(
            names=['Field ID', 'skycoord', 'G', 'BP', 'RP', 'B', 'V', 'R', 'I'],
            data=[field_id, sc, g_mag, bp_mag, rp_mag, B, V, R, I],
            units=(None, None, u.mag, u.mag, u.mag, u.mag, u.mag, u.mag, u.mag)
        )
    if save_table:
        if save_loc is not None:
            ascii.write(ref_mag_table, save_loc, 'csv', overwrite=True)
        else:
            warnings.warn('save_loc must be provided if save_table is True.')
    return ref_mag_table


def read_field_mags(filepath, hdr, exposure_key='EXPTIME'):
    _, xyls = astro.read_fits_file(f'{filepath}.xyls', extension=1)
    fluxes = np.empty(len(xyls))
    for i, star in enumerate(xyls):
        fluxes[i] = star[2]
    exptime = hdr[exposure_key]
    instr_mags = astro.calculate_magnitudes(fluxes, exptime)
    return instr_mags


def calculate_mag_diff(hdr, ref_mag_table, field_id, instr_mags):
    corr_field_stars = instr_mags[field_id]
    # Read the filter here and use that filter's ref mag.
    indx_mags = np.array(ref_mag_table['G'])
    mag_diff = indx_mags - corr_field_stars
    return mag_diff


def average_mag_diff(mag_diff):
    zmag = np.mean(mag_diff)
    zmag_sigma = np.std(mag_diff)
    return zmag, zmag_sigma

try:
    os.remove(f'{filepath}.csv')
except FileNotFoundError:
    pass

directory = r'/media/jmwawrow/Data/DRDC Data/ORC GBO/Field 1/139/Test/'
file_suffix = '.fit'
filecount = 0
file_paths = []
file_names = []
for dirpth, _, files in os.walk(directory):
    for file in files:
        if file.endswith(file_suffix):
            file_paths.append(os.path.join(dirpth, file))
            file_names.append(file)
            filecount += 1

for file_path in file_paths:
    calculate_zmag(file_path, save_table=True)
