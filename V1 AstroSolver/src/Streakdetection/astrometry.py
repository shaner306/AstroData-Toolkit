import os
import sys
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.astrometry_net import AstrometryNet
from math import atan
from os.path import dirname
src_path = dirname(dirname(__file__))
sys.path.append(os.path.join(src_path, 'general_tools'))
sys.path.append(os.path.join(src_path, 'tests_and_verifications'))

import AstroFunctions as astro
import StreakDetectionTest as sd


def solve_plate_astrometry_net(directory, file_suffix=".fits"):
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
        hdr, _ = astro.read_fits_file(file_path)
        focal_length = hdr['FOCALLEN'] * u.mm
        xpixsz = hdr['XPIXSZ']
        ypixsz = hdr['YPIXSZ']
        if xpixsz == ypixsz:
            pixsz = xpixsz * u.um
            rad_per_pix = atan(pixsz / focal_length) * u.rad
            arcsec_per_pix = rad_per_pix.to(u.arcsec)
        low = arcsec_per_pix.value - 0.5
        high = arcsec_per_pix.value + 0.5
        img_radec = SkyCoord(ra=hdr['OBJCTRA'], dec=hdr['OBJCTDEC'], unit=(u.hourangle, u.deg))
        ra = img_radec.ra.deg
        dec = img_radec.dec.deg
        # --ra {ra} --dec {dec} --radius 360
        terminal_call = f'solve-field -p --fits-image --overwrite -D "{directory}/solved" -d 100 -u arcsecperpix \
            -L {low} -H {high} -y "{file_path}"'
        os.system(terminal_call)


def solve_plate_astrometry_net_web(directory, api_key='ldktgfflhujslcyn', file_suffix=".fits"):
    ast = AstrometryNet()
    ast.api_key = api_key
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
        with fits.open(file_path, 'update', memmap=False) as hdul:
            hdr = hdul[0].header
            wcs = WCS(hdr)
            ctypes = wcs.wcs.ctype
            ctype1 = ctypes[0]
            ctype2 = ctypes[1]
            if ctype1 == '' or ctype2 == '':
                sources = sd.StreakDetection(file_path)
                focal_length = hdr['FOCALLEN'] * u.mm
                xpixsz = hdr['XPIXSZ']
                ypixsz = hdr['YPIXSZ']
                if xpixsz == ypixsz:
                    pixsz = xpixsz * u.um
                    rad_per_pix = atan(pixsz / focal_length) * u.rad
                    arcsec_per_pix = rad_per_pix.to(u.arcsec)
                low = arcsec_per_pix.value - 0.1
                high = arcsec_per_pix.value + 0.1
                image_width = hdr['NAXIS1']
                image_height = hdr['NAXIS2']
                wcs_header = ast.solve_from_source_list(sources['xcentroid'], sources['ycentroid'], image_width,
                                                        image_height, scale_lower=low, scale_upper=high,
                                                        publicly_visible='n',
                                                        scale_units='arcsecperpix', solve_timeout=300)

                print('')
                hist_count = 0
                comment_count = 0
                if len(wcs_header) > 5:
                    for keyword in wcs_header[4:]:
                        if keyword == 'HISTORY':
                            hdr.append((keyword, wcs_header[keyword][hist_count]), end=True)
                            hist_count += 1
                        elif keyword == 'COMMENT':
                            hdr.append((keyword, wcs_header[keyword][comment_count]), end=True)
                            comment_count += 1
                        else:
                            value = wcs_header[keyword]
                            comment = wcs_header.comments[keyword]
                            hdr.append((keyword, value, comment), end=True)
                    print('Solved plate. Writing WCS to header...')
                else:
                    print('Could not plate solve.')
            else:
                print('Image is already plate solved. Skipping.')
