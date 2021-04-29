# -*- coding: utf-8 -*-
"""
Test the functionality of the various functions in AstroFunctions.py.

Created on Thu Apr 22 14:58:28 2021

@author: jmwawrow
"""

import AstroFunctions as astro
from astropy.wcs import WCS
from astropy.table import unique
import numpy as np
import os

# ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\FITS Tutorial\Reference_stars_mod.csv'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-03-20 - Calibrated\Solved Images\HIP 2894'

ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\FITS Tutorial\Reference_stars_Apr29.txt'
directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-04-21\Solved Stars'
ground_based = True

# ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Landolt Stars\2009_Landolt_Standard_Stars.txt'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Landolt Stars'
# ground_based = False

reference_stars, ref_star_positions = astro.read_ref_stars(ref_stars_file)
large_table_columns = astro.init_large_table_columns()
gb_transform_table_columns = astro.init_gb_transform_table_columns()

for dirpath, dirnames, filenames in os.walk(directory):
    for filename in filenames:
        if filename.endswith(".fits"):
        # if filename.endswith("_clean.fits"):
            filepath = os.path.join(dirpath, filename)
            print(filepath)
            hdr, imgdata = astro.read_fits_file(filepath)
            exptime = hdr['EXPTIME']
            # exptime = hdr['AEXPTIME']
            bkg, bkg_std = astro.calculate_img_bkg(imgdata)
            irafsources = astro.detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
            print(len(irafsources))
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
            
            instr_filter = astro.get_instr_filter_name(hdr)
            _, _, _, colour_index = astro.get_app_mag_and_index(matched_stars.ref_star, instr_filter)
            field = astro.get_field_name(matched_stars, name_key='Name')
            if np.isnan(matched_stars.ref_star[colour_index]).any():
                no_nan_indices = np.invert(np.isnan(matched_stars.ref_star[colour_index]))
                matched_stars = matched_stars._replace(
                    ref_star_index = matched_stars.ref_star_index[no_nan_indices],
                    img_star_index = matched_stars.img_star_index[no_nan_indices],
                    ref_star = matched_stars.ref_star[no_nan_indices],
                    ref_star_loc = matched_stars.ref_star_loc[no_nan_indices],
                    img_star_loc = matched_stars.img_star_loc[no_nan_indices],
                    ang_separation = matched_stars.ang_separation[no_nan_indices],
                    img_instr_mag = matched_stars.img_instr_mag[no_nan_indices],
                    img_instr_mag_sigma = matched_stars.img_instr_mag_sigma[no_nan_indices],
                    flux = matched_stars.flux[no_nan_indices],
                    img_star_altaz = matched_stars.img_star_altaz[no_nan_indices],
                    img_star_airmass = matched_stars.img_star_airmass[no_nan_indices]
                    )
            try:
                len(matched_stars.img_instr_mag)
            except TypeError:
                print("Only 1 reference star detected in the image.")
                continue
            # print(field)
            c_fci, zprime_f = astro.ground_based_first_order_transforms(matched_stars, instr_filter)
            gb_transform_table_columns = astro.update_gb_transform_table_columns(gb_transform_table_columns,
                                                                                 field,
                                                                                 c_fci,
                                                                                 zprime_f,
                                                                                 instr_filter,
                                                                                 colour_index,
                                                                                 altazpositions)
            
            # large_table_columns = astro.update_large_table_columns(large_table_columns, 
            #                                                         matched_stars, 
            #                                                         hdr, 
            #                                                         exptime, 
            #                                                         ground_based=ground_based, 
            #                                                         name_key='Name')


gb_transform_table = astro.create_gb_transform_table(gb_transform_table_columns)
gb_transform_table.pprint_all()

mask = gb_transform_table['X'] < 2.1
# mask = gb_transform_table['Field'] != '106'
gb_transform_table = gb_transform_table[mask]
gb_transform_table.pprint_all()

gb_final_transforms = astro.ground_based_second_order_transforms(gb_transform_table, plot_results=True)
gb_final_transforms.pprint_all()

# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
# unique_filters = unique(gb_transform_table, keys='filter')
# for unique_filter_row in unique_filters:
#     unique_filter = unique_filter_row['filter']
#     current_index = unique_filter_row['CI']
#     mask = gb_transform_table['filter'] == unique_filter
#     current_filter = gb_transform_table[mask]
#     X_plot = np.arange(start=min(current_filter['X'])-0.2, stop=max(current_filter['X'])+0.2, step=0.1)
#     m, b = np.polyfit(current_filter['X'], current_filter['C_fCI'], 1)
#     plt.plot(current_filter['X'], current_filter['C_fCI'], 'o')
#     plt.plot(X_plot, m*X_plot+b)
#     plt.title(f'C_({unique_filter}{current_index}) = {m:.3f} * X + {b:.3f}')
#     plt.ylabel(f'C_({unique_filter}{current_index})')
#     plt.xlabel('X')
#     plt.show()
#     plt.close()
#     m, b = np.polyfit(current_filter['X'], current_filter['Zprime_f'], 1)
#     plt.plot(current_filter['X'], current_filter['Zprime_f'], 'o')
#     plt.plot(X_plot, m*X_plot+b)
#     plt.title(f'Z\'_({unique_filter}) = {m:.3f} * X + {b:.3f}')
#     plt.ylabel(f'Z\'_({unique_filter})')
#     plt.xlabel('X')
#     plt.show()
#     plt.close()

# large_stars_table = astro.create_large_stars_table(large_table_columns, ground_based=ground_based)
# large_stars_table.pprint_all()
# unique_fields = unique(large_stars_table, keys='Field')
# for unique_field in unique_fields['Field']:
#     if unique_field in list(gb_transform_table['Field']):
#         field_mask = large_stars_table['Field'] == unique_field
#         current_field = large_stars_table[field_mask]
#         unique_filters = unique(current_field, keys='filter')
#         for unique_filter in unique_filters['filter']:
#             filter_mask = current_field['filter'] == unique_filter
#             current_filter = current_field[filter_mask]
#             unique_stars = unique(current_filter, keys='Name')
#             colors = cm.rainbow(np.linspace(0, 1, len(unique_stars)))
#             X_plot = np.arange(start=min(current_filter['X'])-0.02, stop=max(current_filter['X'])+0.02, step=0.01)
#             for i, unique_star in enumerate(unique_stars['Name']):
#                 star_mask = current_filter['Name'] == unique_star
#                 current_star = current_filter[star_mask]
#                 m, b = np.polyfit(current_star['X'], current_star['mag_instrumental'], 1)
#                 plt.scatter(current_star['X'], current_star['mag_instrumental'], color=colors[i], label=unique_star)
#                 plt.plot(X_plot, m*X_plot+b, color=colors[i])
#             # plt.plot(current_filter['X'], current_filter['mag_instrumental'], 'o')
#             plt.xlabel('X')
#             plt.ylabel(unique_filter.lower())
#             plt.ylim([min(current_filter['mag_instrumental'])*1.05, max(current_filter['mag_instrumental'])*0.95])
#             plt.gca().invert_yaxis()
#             plt.title(unique_field)
#             # plt.legend()
#             plt.show()
#             plt.close()
# stars_table = astro.group_each_star(large_stars_table, ground_based=ground_based)
# stars_table.pprint_all()
# transform_index_list = ['(B-V)', '(V-R)', '(V-I)']
# unique_fields = unique(stars_table, keys='Field')
# for field in unique_fields['Field']:
#     mask = stars_table['Field'] == field
#     field_table = stars_table[mask]
#     if len(field_table) == 1:
#         print("Couldn't calculate the transforms as there was only 1 star in the field.")
#         continue
#     for index in transform_index_list:
#         filter_fci, zprime_fci = astro.space_based_transform(field_table, 
#                                                              plot_results=True, 
#                                                              instr_filter='g', 
#                                                              index=index, 
#                                                              field=field)
#         print(f"{field}: (V-g) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
# if ground_based:
#     filter_fci, zprime_fci = astro.space_based_transform(stars_table, plot_results=True, instr_filter='g')
#     print(f"(V-clear) = {filter_fci:.3f} * (B-V) + {zprime_fci:.3f}")
# else:
#     transform_index_list = ['(B-V)', '(V-R)', '(V-I)']
#     for index in transform_index_list:
#         filter_fci, zprime_fci = astro.space_based_transform(stars_table, plot_results=True, index=index)
#         print(f"(V-clear) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")