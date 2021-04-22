# -*- coding: utf-8 -*-
"""
Test the functionality of the various functions in AstroFunctions.py.

Created on Thu Apr 22 14:58:28 2021

@author: jmwawrow
"""

import AstroFunctions as astro
from Astropy.wcs import WCS
import numpy as np
import os

ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\FITS Tutorial\Reference_stars_mod.csv'
directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-03-20 - Calibrated\Solved Images\HIP 2894'
ground_based = True

# ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Landolt Stars\2009_Landolt_Standard_Stars.txt'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Landolt Stars'
# ground_based = False

reference_stars, ref_star_positions = astro.read_ref_stars(ref_stars_file)
large_table_columns = astro.init_large_table_columns()

for dirpath, dirnames, filenames in os.walk(directory):
    for filename in filenames:
        if filename.endswith(".fits"):
        # if filename.endswith("_clean.fits"):
            filepath = os.path.join(dirpath, filename)
            hdr, imgdata = astro.read_fits_file(filepath)
            exptime = hdr['EXPTIME']
            # exptime = hdr['AEXPTIME']
            bkg, bkg_std = astro.calculate_img_bkg(imgdata)
            irafsources = astro.detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
            if not irafsources:
                continue
            fwhm, fwhm_std = astro.calculate_fwhm(irafsources)
            photometry_result = astro.perform_photometry(irafsources, fwhm, imgdata, bkg=bkg)
            fluxes = np.array(photometry_result['flux_fit'])
            instr_mags = astro.calculate_magnitudes(photometry_result, exptime)
            instr_mags_sigma = astro.calculate_magnitudes_sigma(photometry_result, exptime)
            wcs = WCS(hdr)
            skypositions = astro.convert_pixel_to_ra_dec(irafsources, wcs)
            altazpositions = None
            if ground_based:
                altazpositions = astro.convert_ra_dec_to_alt_az(skypositions, hdr)
            matched_stars = astro.find_ref_stars(reference_stars, 
                                                 ref_star_positions,
                                                 skypositions,
                                                 instr_mags,
                                                 instr_mags_sigma,
                                                 fluxes,
                                                 ground_based=ground_based,
                                                 altazpositions=altazpositions)
            if not matched_stars:
                continue
            
            large_table_columns = astro.update_large_table_columns(large_table_columns, 
                                                                   matched_stars, 
                                                                   hdr, 
                                                                   exptime, 
                                                                   ground_based=ground_based, 
                                                                   name_key='HIP')

large_stars_table = astro.create_large_stars_table(large_table_columns, ground_based=ground_based)
large_stars_table.pprint_all()
stars_table = astro.group_each_star(large_stars_table, ground_based=ground_based)
stars_table.pprint_all()
if ground_based:
    filter_fci, zprime_fci = astro.space_based_transform(stars_table, plot_results=True, instr_filter='g')
    print(f"(V-clear) = {filter_fci:.3f} * (B-V) + {zprime_fci:.3f}")
else:
    transform_index_list = ['(B-V)', '(V-R)', '(V-I)']
    for index in transform_index_list:
        filter_fci, zprime_fci = astro.space_based_transform(stars_table, plot_results=True, index=index)
        print(f"(V-clear) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")