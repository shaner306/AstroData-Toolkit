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
sys.path.append(os.path.join(src_path, 'Streakdetection'))

import AstroFunctions as astro
import StreakDetection as sd

def solve_plate_astrometry_net(directory, file_suffix=".fits", sigma=5.0, streakLength=5, TRM=True, useMask=True):
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
        _ = sd.streak_detection_single(file_path, sigma=sigma, streakLength=streakLength, TRM=TRM, useMask=useMask)
        create_xyls_call = f'text2fits -H "x y" -f dd "{file_path}.txt" "{file_path}.xyls"'
        os.system(create_xyls_call)
        image_width = hdr['NAXIS1']
        image_height = hdr['NAXIS2']
        try:
            focal_length = hdr['FOCALLEN'] * u.mm
            xpixsz = hdr['XPIXSZ']
            ypixsz = hdr['YPIXSZ']
        except KeyError:
            terminal_call = f'solve-field -p -N none --overwrite --tag-all -w {image_width} -e {image_height} -X "x" -Y "y" \
                 --objs 100 -d 100 --temp-axy -y -g -W tmp.wcs "{file_path}.xyls"'
            os.system(terminal_call)
            combine_wcs_call = f'new-wcs -i "{file_path}" -w tmp.wcs -o "{file_path}.new" -d'
            os.system(combine_wcs_call)
            if os.path.getsize(f"{file_path}.new") == 0:
                os.remove(f"{file_path}.new")
                print("Could not plate solve.")
                # Flag this and modify the source detection script to try a couple more times.
            continue
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
        terminal_call = f'solve-field -p -N none --overwrite --tag-all -w {image_width} -e {image_height} -X "x" -Y "y" \
                -u arcsecperpix -L {low} -H {high} --temp-axy --objs 100 -d 100 -y -W tmp.wcs "{file_path}.xyls"'
        os.system(terminal_call)
        combine_wcs_call = f'new-wcs -i "{file_path}" -w tmp.wcs -o "{file_path}.new" -d'
        os.system(combine_wcs_call)
        if os.path.getsize(f"{file_path}.new") == 0:
            os.remove(f"{file_path}.new")
            print("Could not plate solve.")

solve_plate_astrometry_net(r'/media/jmwawrow/Data/DRDC Data/2022 07 22/', file_suffix=".fits", TRM=False)

def solve_plate_astrometry_net_web(directory, api_key, file_suffix=".fits"):
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
                sources = sd.detect_streaks(file_path)
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
                                                        image_height, scale_lower=low, scale_upper=high, publicly_visible='n',
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
