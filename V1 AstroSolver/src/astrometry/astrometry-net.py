import os
import sys
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import astropy.units as u
from math import atan
from os.path import dirname
src_path = dirname(dirname(__file__))
sys.path.append(os.path.join(src_path, 'general_tools'))

import AstroFunctions as astro

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
        #--ra {ra} --dec {dec} --radius 360 
        terminal_call = f'solve-field -p --fits-image --overwrite -D "{directory}/solved" -d 100 -u arcsecperpix \
            -L {low} -H {high} -y "{file_path}"'
        os.system(terminal_call)

solve_plate_astrometry_net("/media/jmwawrow/Data/DRDC Data/2021 - Suffield Sky Survey/2021 10 26 - QSI with C8/2021 10 26 - Automated Pointing Run/asatrometry-net test/")