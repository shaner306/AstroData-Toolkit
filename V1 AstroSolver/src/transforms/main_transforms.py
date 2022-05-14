# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 11:49:14 2022

@author: mstew

Stored in this file are all the hgiher levels general_tools needed for transformation
calculations
"""
import math
import os
from math import sqrt
from random import shuffle

# from .AstroFunctions import *
# import TRMtester.py as trm
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from tqdm import tqdm

import sys
from os.path import dirname
src_path = dirname(dirname(__file__))
sys.path.append(os.path.join(src_path, 'general_tools'))
sys.path.append(os.path.join(src_path, 'transforms', 'GroundBasedTransforms','boyde_aux'))
sys.path.append(os.path.join(src_path, 'transforms', 'GroundBasedTransforms','buchheim_aux'))
sys.path.append(os.path.join(src_path, 'transforms', 'GroundBasedTransforms','warner_aux'))
sys.path.append(os.path.join(src_path, 'transforms', 'SpaceBasedTransforms'))
sys.path.append(os.path.join(src_path, 'transforms', 'GroundBasedTransforms'))
sys.path.append(os.path.join(src_path, 'photometry'))
sys.path.append(os.path.join(src_path, 'data_visualization'))

import AstroFunctions as astro
import auxilary_phot_boyde_functions as boyde_aux
import auxilary_phot_buchheim_functions as buch_aux
import auxilary_phot_warner_functions as warn_aux
import auxillary_sb_functions as auxillary_phot_sb_functions
import general_gb_functions
import perform_photometry
from Visualize import plot_match_confirmation
from astropy.io import ascii


def _main_gb_transform_calc(directory,
                            ref_stars_file,
                            plot_results=False,
                            save_plots=False,
                            remove_large_airmass_bool=False,
                            file_suffix=(".fits", ".fit", ".fts"),
                            exposure_key='EXPTIME',
                            lat_key='SITELAT',
                            lon_key='SITELONG',
                            elev_key='SITEELEV',
                            name_key='Name',
                            photometry_method='psf',
                            aperture_estimation_mode='mean',
                            **kwargs):
    """
    Taking 

    Parameters
    ----------

    directory: String of the  Image directory which includes the images to be processed
    ref_stars_file: String describing the directory of reference star files to be used in processing
    plot_results: Boolean 
    save_plots: Boolean
    remove_large_airmass_bool: Boolean
    file_suffix: tuple describing the file suffixes to be used in 
    exposure_key= FITS Keyword to be used to access the exposure time of the image


    Returns

    -------

    Nil

    """
    # TODO: Docstring.
    reference_stars, ref_star_positions = astro.read_ref_stars(ref_stars_file)
    gb_transform_table_columns = astro.init_gb_transform_table_columns()
    auxiliary_data_columns = astro.init_auxiliary_data_columns()
    large_table_columns = astro.init_large_table_columns()

    #if save_plots:
    save_loc = kwargs.get('save_loc')
    if not os.path.exists(save_loc):
        os.mkdir(save_loc)
    with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
        f.write('File')
        f.write('\t')
        f.write('Reason')
        f.write('\n')

    "Create an array of all .fits files in the directory (including subfolders)."
    excluded_files = 0
    filecount = 0
    file_paths = []
    file_names = []
    for dirpth, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_suffix):
                file_paths.append(os.path.join(dirpth, file))
                file_names.append(file)
                filecount += 1
    "Split the files into those for calculation and those for verification."
    shuffle(file_paths)
    split_decimal = 0.7
    split_filecount_location = math.ceil(split_decimal * filecount)
    calculation_files = file_paths[:split_filecount_location]
    verification_files = file_paths[split_filecount_location:]
    with open(os.path.join(save_loc, 'CalVerSplit.txt'), 'a+') as f:
        f.write('File Path'+'\t'+'Calculation/Verification')
        for calc_file in calculation_files:
            f.write('\n'+f'{calc_file}'+'\t'+'Calculation')
        for verify_file in verification_files:
            f.write('\n'+f'{verify_file}'+'\t'+'Verification')

    "Iterate over the images."
    for file_num, filepath in enumerate(tqdm(calculation_files)):
        # for dirpath, dirnames, filenames in os.walk(directory):
        #     for filename in filenames:
        #         if filename.endswith((file_suffix)):
        # filepath = os.path.join(dirpath, filename)
        # print(filepath)
        hdr, imgdata = astro.read_fits_file(filepath)
        exptime = hdr[exposure_key]
        bkg, bkg_std = astro.calculate_img_bkg(imgdata)
        irafsources = astro.detecting_stars(
            imgdata, bkg=bkg, bkg_std=bkg_std)
        if not irafsources:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No sources detected')
                f.write('\n')
            excluded_files += 1
            continue
        fwhms, fwhm, fwhm_std = astro.calculate_fwhm(irafsources)
        if photometry_method=='psf':

            photometry_result = perform_photometry.perform_PSF_photometry(
                irafsources, fwhm, imgdata, bkg=bkg,filepath=filepath,hdr=hdr)

# =============================================================================
# Experimental perform_photometry
#             photometry_result = perform_photometry.perform_PSF_photometry_2(
#                              irafsources, fwhm, imgdata, bkg=bkg,filepath=filepath,hdr=hdr)
# =============================================================================
            # Store the flux and uncertainty of the stars in a separate variable.
            fluxes = np.array(photometry_result['flux_fit'])
            fluxes_unc = np.array(photometry_result['flux_unc'])
            # Convert the flux and uncertainty to magnitude and its uncertainty.
            instr_mags = astro.calculate_magnitudes(fluxes, exptime)
            instr_mags_sigma, snr = astro.calculate_magnitudes_sigma(
                fluxes,fluxes_unc, exptime)
            
            
        elif photometry_method=='aperture':
            
            # Perform Aperture Photometry
            photometry_result=perform_photometry.perform_aperture_photometry(irafsources,fwhms,imgdata,bkg=bkg,
                                                                             bkg_std=np.ones(np.shape(imgdata))*bkg_std,
                                                                             hdr=hdr,filepath=filepath,
                                                                             aperture_estimation_mode=aperture_estimation_mode)
            
            #Re-arrange values to align with PSF Fitting standard
            fluxes_unc=np.transpose(np.array(photometry_result['flux_unc']))
            fluxes=np.array(photometry_result['flux_fit'])
            
            instr_mags=astro.calculate_magnitudes(fluxes,exptime)
            instr_mags_sigma, snr = astro.calculate_magnitudes_sigma(
                fluxes,fluxes_unc, exptime)
        
        wcs = WCS(hdr)
        skypositions = astro.convert_pixel_to_ra_dec(irafsources, wcs)
        altazpositions = None
        try:
            altazpositions = astro.convert_ra_dec_to_alt_az(skypositions,
                                                      hdr,
                                                      lat_key=lat_key,
                                                      lon_key=lon_key,
                                                      elev_key=elev_key)
        except AttributeError as e:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No plate solution found')
                f.write('\n')
                excluded_files += 1
            continue
        matched_stars = astro.find_ref_stars(reference_stars,
                                       ref_star_positions,
                                       skypositions,
                                       instr_mags,
                                       instr_mags_sigma,
                                       snr,
                                       fluxes,
                                       fluxes_unc,
                                       ground_based=True,
                                       altazpositions=altazpositions)
        if not matched_stars:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No reference star found in image')
                f.write('\n')
                excluded_files += 1
            continue

        filename = filepath.split('\\')[-1]
        unique_id = filename
        plot_match_confirmation(wcs, imgdata, matched_stars, reference_stars, unique_id, save_loc, save_plots=save_plots, name_key=name_key)

        instr_filter = astro.get_instr_filter_name(hdr)
        colour_indices = astro.get_all_colour_indices(instr_filter)
        # print("match")
        for colour_index in colour_indices:

            field = astro.get_field_name(matched_stars, name_key=name_key)

            try:
                len(matched_stars.img_instr_mag)
            except TypeError:
                with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                    f.write(f'{filepath}')
                    f.write('\t')
                    f.write('Only 1 reference star detected in the image')
                    f.write('\n')
                    excluded_files += 1
                # print("Only 1 reference star detected in the image.")
                continue

            if np.isnan(matched_stars.ref_star[colour_index]).any():
                no_nan_indices = np.invert(
                    np.isnan(matched_stars.ref_star[colour_index]))
                matched_stars = matched_stars._replace(
                    ref_star_index=matched_stars.ref_star_index[no_nan_indices],
                    img_star_index=matched_stars.img_star_index[no_nan_indices],
                    ref_star=matched_stars.ref_star[no_nan_indices],
                    ref_star_loc=matched_stars.ref_star_loc[no_nan_indices],
                    img_star_loc=matched_stars.img_star_loc[no_nan_indices],
                    ang_separation=matched_stars.ang_separation[no_nan_indices],
                    img_instr_mag=matched_stars.img_instr_mag[no_nan_indices],
                    img_instr_mag_sigma=matched_stars.img_instr_mag_sigma[no_nan_indices],
                    flux=matched_stars.flux[no_nan_indices],
                    img_star_altaz=matched_stars.img_star_altaz[no_nan_indices],
                    img_star_airmass=matched_stars.img_star_airmass[no_nan_indices]
                )
            large_table_columns = astro.update_large_table_columns(large_table_columns,
                                                             filepath,
                                                             matched_stars,
                                                             hdr,
                                                             exptime,
                                                             ground_based=True,
                                                             name_key=name_key)
            auxiliary_data_columns = astro.update_auxiliary_data_columns(auxiliary_data_columns,
                                                                   filepath,
                                                                   exptime,
                                                                   fwhm,
                                                                   fwhm_std,
                                                                   matched_stars)
            # try:
            if not save_plots:
                c_fci, c_fci_sigma, zprime_f, zprime_f_sigma = astro.ground_based_first_order_transforms(matched_stars,
                                                                                                   instr_filter,
                                                                                                   colour_index,
                                                                                                   plot_results=plot_results)
            else:
                c_fci, c_fci_sigma, zprime_f, zprime_f_sigma = astro.ground_based_first_order_transforms(matched_stars,
                                                                                                   instr_filter,
                                                                                                   colour_index,
                                                                                                   plot_results=plot_results,
                                                                                                   save_plots=save_plots,
                                                                                                   save_loc=save_loc,
                                                                                                   unique_id=unique_id)
            # except Exception:
            #     continue
            gb_transform_table_columns = astro.update_gb_transform_table_columns(gb_transform_table_columns,
                                                                           filepath,
                                                                           field,
                                                                           c_fci,
                                                                           c_fci_sigma,
                                                                           zprime_f,
                                                                           zprime_f_sigma,
                                                                           instr_filter,
                                                                           colour_index,
                                                                           altazpositions)
    with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
        f.write('Total excluded:')
        f.write('\t')
        f.write(
            f'{excluded_files} / {split_filecount_location} ({100*(excluded_files/split_filecount_location):.1f}%)')
    large_stars_table = astro.create_large_stars_table(
        large_table_columns, ground_based=True)
    gb_transform_table = astro.create_gb_transform_table(gb_transform_table_columns)
    auxiliary_data_table = astro.create_auxiliary_data_table(auxiliary_data_columns)
    if remove_large_airmass_bool:
        gb_transform_table = astro.remove_large_airmass(gb_transform_table, 3.0)
    if save_plots:
        ascii.write(gb_transform_table,
                    f"{os.path.join(save_loc, 'gb_large_transform_table')}.csv", format='csv')
        ascii.write(auxiliary_data_table,
                    f"{os.path.join(save_loc, 'auxiliary_data_table')}.csv", format='csv')
        ascii.write(large_stars_table, os.path.join(
            save_loc, 'large_stars_table.csv'), format='csv')
        gb_final_transforms = astro.ground_based_second_order_transforms(gb_transform_table,
                                                                   plot_results=plot_results,
                                                                   save_plots=save_plots,
                                                                   save_loc=save_loc)
        formats = {
            'k\'\'_fCI': '%0.3f',
            'k\'\'_fCI_sigma': '%0.3f',
            'T_fCI': '%0.3f',
            'T_fCI_sigma': '%0.3f',
            'k\'_f': '%0.3f',
            'k\'_f_sigma': '%0.3f',
            'Z_f': '%0.3f',
            'Z_f_sigma': '%0.3f'
        }

        ascii.write(gb_final_transforms,
                    f"{os.path.join(save_loc, '_gb_final_transforms')}.csv", format='csv')
        astro.write_table_to_latex(gb_final_transforms, f"{os.path.join(save_loc,'gb_final_transforms')}.txt",
                             formats=formats)
        formats = {
            'exptime': '%0.3f',
            'fwhm': '%0.3f',
            'fwhm_sigma': '%0.3f',
            'avg mag_sigma': '%0.3f',
            'std mag_sigma': '%0.3f'
        }

    else:
        gb_final_transforms = astro.ground_based_second_order_transforms(
            gb_transform_table,
            plot_results=plot_results,
            save_plots=save_plots)
    return gb_final_transforms, auxiliary_data_table

def _main_gb_transform_calc_TEST(directory,
                                 ref_stars_file,
                                 plot_results=False,
                                 save_plots=False,
                                 file_suffix=(".fits", ".fit", ".fts"),
                                 exposure_key='EXPTIME',
                                 name_key='Name',
                                 **kwargs):
    # TODO: Docstring.
    # TODO: Fix errors when save_plots = False.
    reference_stars, ref_star_positions = astro.read_ref_stars(ref_stars_file)
    large_table_columns = astro.init_large_table_columns()

    if save_plots:
        save_loc = kwargs.get('save_loc')
        unique_id = kwargs.get('unique_id')
        if not os.path.exists(save_loc):
            os.mkdir(save_loc)

    "Create the text file for logging problem files."
    with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
        f.write('File')
        f.write('\t')
        f.write('Reason')
        f.write('\n')
    "Create an array of all .fits files in the directory (including subfolders)."
    excluded_files = 0
    filecount = 0
    file_paths = []
    file_names = []
    for dirpth, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_suffix):
                file_paths.append(os.path.join(dirpth, file))
                file_names.append(file)
                filecount += 1
    "Split the files into those for calculation and those for verification."
    shuffle(file_paths)
    split_decimal = 0.7
    split_filecount_location = math.ceil(split_decimal * filecount)
    calculation_files = file_paths[:split_filecount_location]
    verification_files = file_paths[split_filecount_location:]
    with open(os.path.join(save_loc, 'CalVerSplit.txt'), 'a+') as f:
        f.write('File Path'+'\t'+'Calculation/Verification')
        for calc_file in calculation_files:
            f.write('\n'+f'{calc_file}'+'\t'+'Calculation')
        for verify_file in verification_files:
            f.write('\n'+f'{verify_file}'+'\t'+'Verification')

    "Iterate over the images."
    for file_num, filepath in enumerate(tqdm(calculation_files)):
        # filepath = os.path.join(dirpath, filename)
        hdr, imgdata = astro.read_fits_file(filepath)
        exptime = hdr[exposure_key]
        bkg, bkg_std = astro.calculate_img_bkg(imgdata)
        irafsources = astro.detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
        if not irafsources:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No sources detected')
                f.write('\n')
            excluded_files += 1
            continue
        _, fwhm, fwhm_std = astro.calculate_fwhm(irafsources)
        photometry_result = perform_photometry.perform_PSF_photometry(
            irafsources, fwhm, imgdata, bkg=bkg)
        fluxes = np.array(photometry_result['flux_fit'])
        fluxes_unc = np.array(photometry_result['flux_unc'])
        instr_mags = astro.calculate_magnitudes(fluxes, exptime)
        instr_mags_sigma, snr = astro.calculate_magnitudes_sigma(
            fluxes,fluxes_unc, exptime)
        wcs = WCS(hdr)
        skypositions = astro.convert_pixel_to_ra_dec(irafsources, wcs)
        try:
            altazpositions = astro.convert_ra_dec_to_alt_az(skypositions, hdr,
                                                      lat_key='SITELAT',
                                                      lon_key='SITELONG',
                                                      elev_key='SITEELEV')
        except AttributeError as e:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No plate solution found')
                f.write('\n')
                excluded_files += 1
            continue
        matched_stars = astro.find_ref_stars(reference_stars,
                                       ref_star_positions,
                                       skypositions,
                                       instr_mags,
                                       instr_mags_sigma,
                                       snr,
                                       fluxes,
                                       fluxes_unc,
                                       ground_based=True,
                                       altazpositions=altazpositions)
        if not matched_stars:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No reference star found in image')
                f.write('\n')
                excluded_files += 1
            continue

        large_table_columns = astro.update_large_table_columns(large_table_columns,
                                                         filepath,
                                                         matched_stars,
                                                         hdr,
                                                         exptime,
                                                         ground_based=True,
                                                         name_key=name_key)
    with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
        f.write('Total excluded:')
        f.write('\t')
        f.write(
            f'{excluded_files} / {split_filecount_location} ({100*(excluded_files/split_filecount_location):.1f}%)')
    large_stars_table = astro.create_large_stars_table(
        large_table_columns, ground_based=True)
    # large_stars_table = astro.read_fits_file(large_stars_table, max_airmass=3.0)
    stars_table, different_filter_list = astro.group_each_star_GB(large_stars_table)
    stars_table.pprint(max_lines=30, max_width=200)
    if save_plots:
        ascii.write(stars_table, os.path.join(
            save_loc, 'stars_table.csv'), format='csv')
    gb_transform_table =\
        general_gb_functions.calc_gb_first_transforms_AVG(stars_table,
                                                          different_filter_list,
                                                          save_loc,
                                                          plot_results=plot_results,
                                                          save_plots=save_plots)
    gb_transform_table.pprint(max_lines=30, max_width=-1)
    gb_final_transforms =\
        astro.ground_based_second_order_transforms(gb_transform_table,
                                             plot_results=plot_results,
                                             save_plots=save_plots,
                                             save_loc=save_loc)
    gb_final_transforms.pprint_all()
    if save_plots:
        gb_final_transforms =\
            astro.ground_based_second_order_transforms(gb_transform_table,
                                                 plot_results=plot_results,
                                                 save_plots=save_plots,
                                                 save_loc=save_loc)
        formats = {
            'k\'\'_fCI': '%0.3f',
            'k\'\'_fCI_sigma': '%0.3f',
            'T_fCI': '%0.3f',
            'T_fCI_sigma': '%0.3f',
            'k\'_f': '%0.3f',
            'k\'_f_sigma': '%0.3f',
            'Z_f': '%0.3f',
            'Z_f_sigma': '%0.3f'
        }
        ascii.write(gb_final_transforms,
                    f"{os.path.join(save_loc, '_gb_final_transforms')}.csv",
                    format='csv')
        astro.write_table_to_latex(gb_final_transforms,
                             f"{os.path.join(save_loc, 'gb_final_transforms')}.txt",
                             formats=formats)
    hidden_transform_table = None
    exoatmospheric_table = astro.exoatmospheric_mags_Warner(
        stars_table, gb_final_transforms, different_filter_list)
    hidden_transform_table =\
        astro.hidden_transform_Warner(exoatmospheric_table,
                                gb_final_transforms,
                                different_filter_list,
                                save_plots,
                                save_loc=save_loc)
    verify_save_loc = os.path.join(save_loc, 'Verification')
    app_mag_table = astro.verify_gb_transforms_auto(directory,
                                              verification_files,
                                              ref_stars_file,
                                              gb_final_transforms,
                                              hidden_transform_table,
                                              plot_results=True,
                                              save_plots=True,
                                              file_suffix=file_suffix,
                                              exposure_key=exposure_key,
                                              name_key=name_key,
                                              save_loc=verify_save_loc)
    return large_stars_table


def _main_gb_transform_calc_Warner(directory,  # Light Frames
                                   ref_stars_file,  # Reference Stars Files
                                   plot_results=False,
                                   save_plots=False,
                                   file_suffix=[],
                                   exposure_key='EXPTIME',
                                   name_key='Name',
                                   lat_key='SITELAT',
                                   lon_key='SITELONG',
                                   elev_key='SITEELEV',
                                   photometry_method='psf',
                                   aperture_estimation_error='mean',
                                   **kwargs):
    # TODO: Docstring.
    # TODO: Fix errors when save_plots = False.
    # TODO: Make Save_loc mandatory
    """
    Perform all of the beginning operations.
    Create the refrence stars table and read their positions.
    Initialize empty arrays for the star information and auxiliary info.
    """
    reference_stars, ref_star_positions = astro.read_ref_stars(ref_stars_file)
    large_table_columns = astro.init_large_table_columns()
    star_aux_table_columns = astro.init_star_aux_table_columns()

    "Create the save location if the one specified by the user doesn't exist."

    # FIXME: Add Try and Except for when sav_loc isn't an input into the kwargs in the function

    if save_plots:
        save_loc = kwargs.get('save_loc')
        unique_id = kwargs.get('unique_id')
        if not os.path.exists(save_loc):
            os.mkdir(save_loc)

    "Create the text file for logging problem files."
    with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
        f.write('File')
        f.write('\t')
        f.write('Reason')
        f.write('\n')
    "Create an array of all .fits files in the directory (including subfolders)."
    excluded_files = 0
    filecount = 0
    file_paths = []
    file_names = []
    for dirpth, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_suffix):
                file_paths.append(os.path.join(dirpth, file))
                file_names.append(file)
                filecount += 1
    "Split the files into those for calculation and those for verification."

    shuffle(file_paths)
    split_decimal = 1
    split_filecount_location = math.ceil(split_decimal * filecount)
    calculation_files = file_paths[:split_filecount_location]
    verification_files = file_paths[split_filecount_location:]
    with open(os.path.join(save_loc, 'CalVerSplit.txt'), 'a+') as f:
        f.write('File Path'+'\t'+'Calculation/Verification')
        for calc_file in calculation_files:
            f.write('\n'+f'{calc_file}'+'\t'+'Calculation')
        for verify_file in verification_files:
            f.write('\n'+f'{verify_file}'+'\t'+'Verification')
    "Iterate over the images."
    for file_num, filepath in enumerate(tqdm(calculation_files)):
        # Read the fits file. Stores the header and image to variables.
        hdr, imgdata = astro.read_fits_file(filepath)
        # Read the exposure time of the image.
        exptime = hdr[exposure_key]
        # Calculate the image background and standard deviation.
        bkg, bkg_std = astro.calculate_img_bkg(imgdata)
        # Detect point sources in the image.
        irafsources = astro.detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
        # If no stars are in the image, log it and go to the next one.
        if not irafsources:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No sources detected')
                f.write('\n')
            excluded_files += 1
            continue
        # Calculate the number of point sources detected.
        num_sources = len(irafsources)
        # Calculate the FWHM in pixels.
        fwhms, fwhm, fwhm_std = astro.calculate_fwhm(irafsources)
        # Do PSF photometry on the detected sources.
        if photometry_method=='psf':

            photometry_result = perform_photometry.perform_PSF_photometry(
                irafsources, fwhm, imgdata, bkg=bkg,filepath=filepath,hdr=hdr)

# =============================================================================
# Experimental perform_photometry
#             photometry_result = perform_photometry.perform_PSF_photometry_2(
#                              irafsources, fwhm, imgdata, bkg=bkg,filepath=filepath,hdr=hdr)
# =============================================================================
            # Store the flux and uncertainty of the stars in a separate variable.
            fluxes = np.array(photometry_result['flux_fit'])
            fluxes_unc = np.array(photometry_result['flux_unc'])
            # Convert the flux and uncertainty to magnitude and its uncertainty.
            instr_mags = astro.calculate_magnitudes(fluxes, exptime)
            instr_mags_sigma, snr = astro.calculate_magnitudes_sigma(
                fluxes,fluxes_unc, exptime)
            
            
        elif photometry_method=='aperture':
            
            # Perform Aperture Photometry
            photometry_result=perform_photometry.perform_aperture_photometry(irafsources,fwhms,imgdata,bkg=bkg,bkg_std=np.ones(np.shape(imgdata))*bkg_std,hdr=hdr,filepath=filepath,aperture_estimation_mode=aperture_estimation_error)
            
            #Re-arrange values to align with PSF Fitting standard
            fluxes_unc=np.transpose(np.array(photometry_result['flux_unc']))
            fluxes=np.array(photometry_result['flux_fit'])
            instr_mags=astro.calculate_magnitudes(fluxes,exptime)
            instr_mags_sigma, snr = astro.calculate_magnitudes_sigma(
                fluxes,fluxes_unc, exptime)
        
        
        
        
        # Read the World Coordinate System transformation added to the fits header
        # by a plate solving software (external to this program, e.g. PinPoint).

        # FIXME : Existing WCS Data must be present to plate solve.
        # PinPoint Should be called and solved within the program

        wcs = WCS(hdr)
        # Convert the stars' (x,y) location to (RA,dec).
        skypositions = astro.convert_pixel_to_ra_dec(irafsources, wcs)
        try:
            # altazpositions = astro.convert_ra_dec_to_alt_az(skypositions, hdr, lat_key='OBSGEO-B', lon_key='OBSGEO-L',
            #                                           elev_key='OBSGEO-H')

            # FIXME: SkyCoord System does not work

            # Convert the stars' (RA,dec) location to (Azimuth,Elevation).
            altazpositions = astro.convert_ra_dec_to_alt_az(skypositions, hdr, lat_key=lat_key, lon_key=lon_key,
                                                      elev_key=elev_key)
        # If there's no plate solution, it will raise and AttributeError.
        # This catches that error, logs it, and moves onto the next image.
        except AttributeError as e:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No plate solution found')
                f.write('\n')
                excluded_files += 1
            continue
        # Convert the FWHM from pixels to arcsec.
        fwhms_arcsec, fwhm_arcsec, fwhms_arcsec_std = astro.convert_fwhm_to_arcsec(
            hdr, fwhms, fwhm, fwhm_std)
        # If it can't convert from pixels to arcsec (e.g. the focal length wasn't defined in the header),
        # store it as NaN.
        if not fwhm_arcsec:
            fwhm_arcsec = np.nan
            fwhms_arcsec_std = np.nan
        # Read the time from the fits header and then convert it to Julian Date.
        t = Time(hdr['DATE-OBS'], format='fits', scale='utc')
        time = t.jd
        # Store the filter used to take the image as a variable.
        img_filter = hdr['FILTER']
        # Calculate the background sky brightness and standard deviation.
        background_sky_brightness = astro.calculate_background_sky_brightness(
            bkg, hdr, exptime)
        background_sky_brightness_sigma = astro.calculate_BSB_sigma(
            bkg, bkg_std, exptime)
        # Take the average of all stars' Az/El/airmass and store as a variable.
        azimuth = np.mean(altazpositions.az)
        elevation = np.mean(altazpositions.alt)
        airmass = np.mean(altazpositions.secz)
        # Update the table with auxiliary data on the images (FWHM, BSB, etc.)
        star_aux_table_columns = astro.update_star_aux_columns(star_aux_table_columns,
                                                         file_names[file_num],
                                                         time,
                                                         img_filter,
                                                         fwhm,
                                                         fwhm_std,
                                                         fwhm_arcsec,
                                                         fwhms_arcsec_std,
                                                         num_sources,
                                                         background_sky_brightness,
                                                         background_sky_brightness_sigma,
                                                         azimuth,
                                                         elevation,
                                                         airmass)
        # Match the detected sources with a star from the reference stars file.
        matched_stars = astro.find_ref_stars(reference_stars,
                                       ref_star_positions,
                                       skypositions,
                                       instr_mags,
                                       instr_mags_sigma,
                                       snr,
                                       fluxes,
                                       fluxes_unc,
                                       ground_based=True,
                                       altazpositions=altazpositions)
        # If no image star corresponds to a reference star, log it and go to the next image.
        if not matched_stars:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No reference star found in image')
                f.write('\n')
                excluded_files += 1
            continue
        # Update the table that contains information on each detection of a reference star.
        large_table_columns = astro.update_large_table_columns(large_table_columns,
                                                         filepath,
                                                         matched_stars,
                                                         hdr,
                                                         exptime,
                                                         ground_based=True,
                                                         name_key=name_key)
    # Complete the text file that stores information on files that were not used to calculate the transforms.
    with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
        f.write('Total excluded:')
        f.write('\t')
        f.write(
            f'{excluded_files} / {split_filecount_location} ({100*(excluded_files/split_filecount_location):.1f}%)')
    # Create an AstroPy table of the auxiliary data and write it to a .csv file.
    star_aux_table = astro.create_star_aux_table(star_aux_table_columns)
    ascii.write(star_aux_table, os.path.join(
        save_loc, 'auxiliary_table.csv'), format='csv')
    # Create an AstroPy table of each reference star detection and write it to a .csv file.
    large_stars_table = astro.create_large_stars_table(
        large_table_columns, ground_based=True)
    ascii.write(large_stars_table, os.path.join(
        save_loc, 'large_stars_table.csv'), format='csv')
    # Group each observation of a star at an airmass.
    # E.g. if there are 5 images of star X at 1.2 airmass, and 10 images of star X at 2 airmass,
    # it will produce a mean and standard deviation of the observations at both 1.2 and 2 airmass.
    # This creates that table, stores the different filters used to take the images (e.g. BVRI or BGR),
    # and writes it to a .csv file.
    stars_table, different_filter_list = astro.group_each_star_GB(large_stars_table)
    ascii.write(stars_table, os.path.join(
        save_loc, 'stars_table.csv'), format='csv')

    # TRYING SOMETHING

    # reformatted_large_stars_table = create_reformatted_large_table(large_stars_table, keys='Name')
    # ascii.write(reformatted_large_stars_table, os.path.join(save_loc, 'reformatted_large_stars_table.csv'), format='csv')

    ############# Begin the Warner Transforms #############

    # Calculate the slope of each star's instrumental magnitude vs airmass,
    # store it in a table, and write it to a .csv file.
    slopes_table = warn_aux.calculate_slopes_Warner(
        stars_table, different_filter_list, save_plots, save_loc=save_loc)
    ascii.write(slopes_table, os.path.join(
        save_loc, 'slopes_table.csv'), format='csv')
    # Calculate the first and second order extinctions.
    extinction_table_warner = warn_aux.second_order_extinction_calc_warner(slopes_table,
                                                                  different_filter_list,
                                                                  save_plots,
                                                                  save_loc=save_loc,**kwargs)
    # Calculate the exoatmospheric magnitudes (m_0).
    exoatmospheric_table = warn_aux.exoatmospheric_mags_warner(
        stars_table, extinction_table_warner, different_filter_list)
    # Finish the transform by calculating the colour transform and zero point.
    warner_final_transform_table = warn_aux.colour_transform_and_zp_calc_Warner(exoatmospheric_table,
                                                                       different_filter_list,
                                                                       extinction_table_warner, save_plots,
                                                                       save_loc=save_loc)
    # Save the transform table to a .csv file.
    ascii.write(warner_final_transform_table, os.path.join(
        save_loc, '_gb_final_transforms.csv'), format='csv')
    # Calculate the hidden transform and write it to a .csv file.
    hidden_transform_table = warn_aux.hidden_transform_Warner(exoatmospheric_table,
                                                     warner_final_transform_table,
                                                     different_filter_list,
                                                     save_plots,
                                                     save_loc=save_loc)
    ascii.write(hidden_transform_table, os.path.join(
        save_loc, 'hidden_transform_table.csv'), format='csv')
    exoatmospheric_table_verify = warn_aux.exoatmospheric_mags_verify_Warner(stars_table,
                                                                    extinction_table_warner,
                                                                    hidden_transform_table,
                                                                    different_filter_list)
    # Verify the transforms.
    verify_save_loc = os.path.join(save_loc, 'Verification')
    app_mag_table = astro.verify_gb_transforms_auto(directory,
                                              verification_files,
                                              ref_stars_file,
                                              warner_final_transform_table,
                                              hidden_transform_table,
                                              plot_results=True,
                                              save_plots=True,
                                              file_suffix=file_suffix,
                                              exposure_key=exposure_key,
                                              name_key=name_key,
                                              lat_key=lat_key,
                                              lon_key=lon_key,
                                              elev_key=elev_key,
                                              save_loc=verify_save_loc)
    return warner_final_transform_table



def _main_gb_new_boyd_method(
        directory,
        ref_stars_file,
        plot_results=False,
        save_plots=True,
        remove_large_airmass_bool=False,
        file_suffix=(".fits", ".fit", ".fts"),
        exposure_key='EXPTIME',
        lat_key='SITELAT',
        lon_key='SITELONG',
        elev_key='SITEELEV',
        name_key='Name',
        photometry_method='aperture',
        aperture_estimation_mode='mean',
        **kwargs):
    '''
    A derivative of the _main_gb_transform_calc which uses the Boyde Method for
    calculating the colour transforms to convert instrumental magntiude to standard magntiude

    $ TODO: Docstring
    Parameters
    ----------
    directory : TYPE
        describes the parent directory which holds 
    ref_stars_file : TYPE
        DESCRIPTION.
    plot_results : TYPE, optional
        DESCRIPTION. The default is False.
    save_plots : TYPE, optional
        DESCRIPTION. The default is False.
    remove_large_airmass_bool : TYPE, optional
        DESCRIPTION. The default is False.
    file_suffix : TYPE, optional
        DESCRIPTION. The default is (".fits", ".fit", ".fts").
    exposure_key : TYPE, optional
        DESCRIPTION. The default is 'EXPTIME'.
    lat_key : TYPE, optional
        DESCRIPTION. The default is 'SITELAT'.
    lon_key : TYPE, optional
        DESCRIPTION. The default is 'SITELONG'.
    elev_key : TYPE, optional
        DESCRIPTION. The default is 'SITEELEV'.
    name_key : TYPE, optional
        DESCRIPTION. The default is 'Name'.
    photometry_method : TYPE, optional
        DESCRIPTION. The default is 'psf'.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

'''

    
    
    matched_stars_dict={}
    reference_stars, ref_star_positions = astro.read_ref_stars(ref_stars_file)
    large_table_columns = astro.init_large_table_columns()
    star_aux_table_columns = astro.init_star_aux_table_columns()

    "Create the save location if the one specified by the user doesn't exist."
    #if save_plots:
    save_loc = kwargs.get('save_loc')
    unique_id = kwargs.get('unique_id')
    if not os.path.exists(save_loc):
        os.mkdir(save_loc)


    "Create the text file for logging problem files."
    with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
        f.write('File')
        f.write('\t')
        f.write('Reason')
        f.write('\n')
    "Create an array of all .fits files in the directory (including subfolders)."
    excluded_files = 0
    filecount = 0
    file_paths = []
    file_names = []
    for dirpth, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_suffix):
                file_paths.append(os.path.join(dirpth, file))
                file_names.append(file)
                filecount += 1
    "Split the files into those for calculation and those for verification."
    shuffle(file_paths)
    split_decimal = 1  # when decimal is 1 then all
    split_filecount_location = math.ceil(split_decimal * filecount)
    calculation_files = file_paths[:split_filecount_location]
    verification_files = file_paths[split_filecount_location:]
    with open(os.path.join(save_loc,
                           'CalVerSplit.txt'),
              'a+') as f:
        f.write('File Path'+'\t'+'Calculation/Verification')
        for calc_file in calculation_files:
            f.write('\n'+f'{calc_file}'+'\t'+'Calculation')
        for verify_file in verification_files:
            f.write('\n'+f'{verify_file}'+'\t'+'Verification')
    "Iterate over the images."


    for file_num, filepath in enumerate(tqdm(calculation_files)):
        # Read the fits file. Stores the header and image to variables.
        hdr, imgdata = astro.read_fits_file(filepath)
        # Read the exposure time of the image.
        exptime = hdr[exposure_key]
        # Calculate the image background and standard deviation.
        bkg, bkg_std = astro.calculate_img_bkg(imgdata)
        # Detect point sources in the image.
        try:
            irafsources = astro.detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std,fwhm=hdr['FWHM'])
        except:
            print("Could not find FWHM Tag")
            irafsources = astro.detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std, fwhm=3)
        
        
        
        # If no stars are in the image, log it and go to the next one.
        if not irafsources:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No sources detected')
                f.write('\n')
            excluded_files += 1
            continue
        # Calculate the number of point sources detected.
        num_sources = len(irafsources)
        # Calculate the FWHM in pixels.
        fwhms, fwhm, fwhm_std = astro.calculate_fwhm(irafsources)
        # Do PSF photometry on the detected sources.
        if photometry_method=='psf':

            photometry_result = perform_photometry.perform_PSF_photometry(
                irafsources, fwhm, imgdata, bkg=bkg,filepath=filepath,hdr=hdr)

# =============================================================================
# Experimental perform_photometry
#             photometry_result = perform_photometry.perform_PSF_photometry_2(
#                              irafsources, fwhm, imgdata, bkg=bkg,filepath=filepath,hdr=hdr)
# =============================================================================
            # Store the flux and uncertainty of the stars in a separate variable.
            fluxes = np.array(photometry_result['flux_fit'])
            fluxes_unc = np.array(photometry_result['flux_unc'])
            # Convert the flux and uncertainty to magnitude and its uncertainty.
            instr_mags = astro.calculate_magnitudes(fluxes, exptime)
            instr_mags_sigma, snr = astro.calculate_magnitudes_sigma(
                fluxes,fluxes_unc, exptime)
            
            
        elif photometry_method=='aperture':
            try:
                # Perform Aperture Photometry

                photometry_result=perform_photometry.perform_aperture_photometry(irafsources,fwhms,imgdata,
                                                                                 bkg=bkg,
                                                                                 bkg_std=np.ones(np.shape(imgdata))*bkg_std,
                                                                                 hdr=hdr,filepath=filepath,
                                                                                 aperture_estimation_mode=aperture_estimation_mode)
                
                #Re-arrange values to align with PSF Fitting standard
                fluxes_unc=(np.array(photometry_result['aperture_sum_err']))
                fluxes=np.array(photometry_result['aperture_sum'])
                instr_mags=astro.calculate_magnitudes(fluxes,exptime)
                instr_mags_sigma, snr = astro.calculate_magnitudes_sigma(
                    fluxes,fluxes_unc, exptime)
            except Exception as exception:
                raise RuntimeError(exception)
        
        # Read the World Coordinate System transformation added to the
        # fits header by a plate solving software (external to this program, e.g. PinPoint).
        wcs = WCS(hdr)
        # Convert the stars' (x,y) centroid locations to (RA,dec).
        skypositions = astro.convert_pixel_to_ra_dec(irafsources, wcs)
        try:
            # altazpositions = astro.convert_ra_dec_to_alt_az(skypositions, hdr,
            # lat_key='OBSGEO-B', lon_key='OBSGEO-L',
            #                                           elev_key='OBSGEO-H')

            # Convert the stars' (RA,dec) location to (Azimuth,Elevation).
            altazpositions = astro.convert_ra_dec_to_alt_az(skypositions, hdr,
                                                      lat_key=lat_key,
                                                      lon_key=lon_key,
                                                      elev_key=elev_key)
        # If there's no plate solution, it will raise and AttributeError.
        # This catches that error, logs it, and moves onto the next image.
        except AttributeError as e:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No plate solution found')
                f.write('\n')
                excluded_files += 1
            continue
        # Convert the FWHM from pixels to arcsec.
        fwhms_arcsec, fwhm_arcsec, fwhms_arcsec_std = astro.convert_fwhm_to_arcsec(
            hdr, fwhms, fwhm, fwhm_std)
        # If it can't convert from pixels to arcsec (e.g. the focal length
        # wasn't defined in the header),
        # store it as NaN.
        if not fwhm_arcsec:
            fwhm_arcsec = np.nan
            fwhms_arcsec_std = np.nan
        # Read the time from the fits header and then convert it to
        # Julian Date.
        t = Time(hdr['DATE-OBS'], format='fits', scale='utc')
        time = t.jd
        # Store the filter used to take the image as a variable.
        img_filter = hdr['FILTER']
        # Calculate the background sky brightness and standard deviation.
        background_sky_brightness = astro.calculate_background_sky_brightness(
            bkg, hdr, exptime)
        background_sky_brightness_sigma = astro.calculate_BSB_sigma(
            bkg, bkg_std, exptime)
        # Take the average of all stars' Az/El/airmass and store as a variable.
        azimuth = np.mean(altazpositions.az)
        elevation = np.mean(altazpositions.alt)
        airmass = np.mean(altazpositions.secz)

        
        try:
            predicted_airmass=float(hdr['AIRMASS'])
        except KeyError:
            # AIRMASS header not found, attempt to estimate from WCS Data
            wcs = WCS(hdr)
            altazpositions = astro.convert_ra_dec_to_alt_az((
                wcs.pixel_to_world(wcs.pixel_shape[0] / 2, wcs.pixel_shape[1] / 2)),hdr)
            predicted_airmass = altazpositions.secz.value

            # Other Method for predicting airmass
            # =============================================================================
            #         predicted_airmass = 1 / \
            #             np.cos(np.deg2rad(
            #                 90-(CCDData.read(filepath, unit='adu').header['CENTALT'])))
            # =============================================================================



        #airmass_std=np.sqrt(((altazpositions.secz-predicted_airmass)**2)/(np.count_nonzero(altazpositions.secz)-1))
        
        #if hdr['CENTALT']<30:
            # If the fits header airmass is greater than two, the airmass will start to change rapidly with the elevation angles of the stars
            
            # for now we will ignore all values above 2, but further calcualtions will need to be prodcued 
            
            # TODO: Write code that handles high variability airmasses (i.e airmass>2)
            # continue
        #else:
            # TODO: Write Code that converts the Pinpoint Solved WCS RA DEC data to AltAz instead of using the predicted data 
        airmass = predicted_airmass

            
        # Update the table with auxiliary data on the images (FWHM, BSB, etc.)
        star_aux_table_columns =\
            astro.update_star_aux_columns(star_aux_table_columns,
                                    file_names[file_num],
                                    time,
                                    img_filter,
                                    fwhm,
                                    fwhm_std,
                                    fwhm_arcsec,
                                    fwhms_arcsec_std,
                                    num_sources,
                                    background_sky_brightness,
                                    background_sky_brightness_sigma,
                                    azimuth,
                                    elevation,
                                    airmass)
        # Match the detected sources with a star from the reference stars file.
        matched_stars = astro.find_ref_stars(reference_stars,
                                       ref_star_positions,
                                       skypositions,
                                       instr_mags,
                                       instr_mags_sigma,
                                       snr,
                                       fluxes,
                                       fluxes_unc,
                                       ground_based=True,
                                       altazpositions=altazpositions)
        matched_stars_dict[filepath]=matched_stars
        # If no image star corresponds to a reference star, log it and go to
        # the next image.
        if not matched_stars:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No reference star found in image')
                f.write('\n')
                excluded_files += 1
            continue

        # Step 1

        # For each image calculate the slope and the intercept of V_ref-v_instrumental vs. Colour indices

        # Update the table that contains information on each detection of a
        # reference star.
        large_table_columns = astro.update_large_table_columns(large_table_columns,
                                                         filepath,
                                                         matched_stars,
                                                         hdr,
                                                         exptime,
                                                         ground_based=True,
                                                         name_key=name_key
                                                         )

    # Complete the text file that stores information on files that were not used to calculate the transforms.
    with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
        f.write('Total excluded:')
        f.write('\t')
        f.write(
            f'{excluded_files} / {split_filecount_location} ({100*(excluded_files/split_filecount_location):.1f}%)')
    # Create an AstroPy table of the auxiliary data and write it to a .csv file.
    
    
    
    # Group each observation of a star at an airmass.
    # E.g. if there are 5 images of star X at 1.2 airmass, and 10 images of
    # star X at 2 airmass,
    # it will produce a mean and standard deviation of the observations at
    # both 1.2 and 2 airmass.
    # This creates that table, stores the different filters used to take the
    # images (e.g. BVRI or BGR),
    # and writes it to a .csv file.
    try:
        star_aux_table = astro.create_star_aux_table(star_aux_table_columns)
        ascii.write(star_aux_table, os.path.join(
            save_loc, 'auxiliary_table.csv'), format='csv',overwrite=True)
        # Create an AstroPy table of each reference star detection and write it to a .csv file.
        large_stars_table = astro.create_large_stars_table(
            large_table_columns, ground_based=True)
        ascii.write(large_stars_table, os.path.join(
            save_loc, 'large_stars_table.csv'), format='csv',overwrite=True)
        stars_table, different_filter_list = astro.group_each_star_GB(large_stars_table)
        ascii.write(stars_table, os.path.join(
            save_loc, 'stars_table.csv'), format='csv',overwrite=True)
    except Exception as e:
        print(e)

        # raise Exception(e)

        
    ### Start Boyd Transformation ###
        # Create Boyde Table
    Boyde_Table = Table(names=['Image Name', 'C', 'Z-prime', 'Index (i.e. B-V)', 'Airmass','Air Mass Std',
                               'Colour Filter', 'Step1_Standard_Deviation', 'Number of Valid Matched Stars'],
                        dtype=('str', 'float64', 'float64', 'str', 'float64', 'float64','str', 'float64', 'int'))
    for filepath in matched_stars_dict:
        
        matched_stars=matched_stars_dict[filepath]
        
        
        
        
        #if predicted_airmass > 3:
            # If the fits header airmass is greater than two, the airmass will start to change rapidly with the elevation angles of the stars

            # for now we will ignore all values above 2, but further calcualtions will need to be prodcued

            # TODO: Write code that handles high variability airmasses (i.e airmass>2)
        #    continue
        #else:


        try:
            Boyde_Table = boyde_aux.calculate_boyde_slopes(
                matched_stars, filepath, Boyde_Table, save_plots,
                save_loc,stars_table)
        except Exception as e:
            raise KeyError(e)
            print("Could not Calculate Boyde Slopes")
            continue


    # Calculating Second Step of Boyd Method
    try:
        ascii.write(Boyde_Table, os.path.join(
            save_loc, 'Boyde_Table1.csv'), format='csv')
    except:
        print('Could Not Save Boyde Table')
    
    
    
    
    # Calculate Step 2
    match_stars_lim=4

    try:
        Boyde_Table_grouped = boyde_aux.calculate_boyde_slope_2(
            Boyde_Table, save_loc,match_stars_lim,save_plots)

        date_data=boyde_aux.create_coefficeint_output(Boyde_Table_grouped)

        ascii.write(Boyde_Table_grouped, os.path.join(
            save_loc, 'Boyde_Table2.csv'), format='csv',overwrite=True)

        ascii.write(date_data,os.path.join(
            save_loc, 'Date_data.csv'), format='csv',overwrite=True)
    except Exception:

        print(Exception)


def _main_gb_transform_calc_Buchheim(directory,

                                     ref_stars_file,
                                     plot_results=False,
                                     save_plots=False,
                                     file_suffix=(".fits", ".fit", ".fts"),
                                     exposure_key='EXPTIME',
                                     name_key='Name',
                                     lat_key='SITELAT',
                                     lon_key='SITELONG',
                                     elev_key='SITEELEV',
                                     photometry_method='psf',
                                     aperture_estimaiton_mode='mean',
                                     **kwargs):

    # TODO: Docstring.
    # TODO: Fix errors when save_plots = False.
    """
    Perform all of the beginning operations.
    Create the refrence stars table and read their positions.
    Initialize empty arrays for the star information and auxiliary info.
    """
    reference_stars, ref_star_positions = astro.read_ref_stars(ref_stars_file)
    large_table_columns = astro.init_large_table_columns()
    star_aux_table_columns = astro.init_star_aux_table_columns()

    "Create the save location if the one specified by the user doesn't exist."
    if save_plots:
        save_loc = kwargs.get('save_loc')
        unique_id = kwargs.get('unique_id')
        if not os.path.exists(save_loc):
            os.mkdir(save_loc)

    "Create the text file for logging problem files."
    with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
        f.write('File')
        f.write('\t')
        f.write('Reason')
        f.write('\n')
    "Create an array of all .fits files in the directory (including subfolders)."
    excluded_files = 0
    filecount = 0
    file_paths = []
    file_names = []
    for dirpth, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_suffix):
                file_paths.append(os.path.join(dirpth, file))
                file_names.append(file)
                filecount += 1
    "Split the files into those for calculation and those for verification."
    shuffle(file_paths)  # why shuffle?
    split_decimal = 1  # when decimal is 1 then all
    split_filecount_location = math.ceil(split_decimal * filecount)
    calculation_files = file_paths[:split_filecount_location]
    verification_files = file_paths[split_filecount_location:]
    with open(os.path.join(save_loc,
                           'CalVerSplit.txt'),
              'a+') as f:
        f.write('File Path'+'\t'+'Calculation/Verification')
        for calc_file in calculation_files:
            f.write('\n'+f'{calc_file}'+'\t'+'Calculation')
        for verify_file in verification_files:
            f.write('\n'+f'{verify_file}'+'\t'+'Verification')
    "Iterate over the images."
    for file_num, filepath in enumerate(tqdm(calculation_files)):
        # Read the fits file. Stores the header and image to variables.
        hdr, imgdata = astro.read_fits_file(filepath)
        # Read the exposure time of the image.
        exptime = hdr[exposure_key]
        # Calculate the image background and standard deviation.
        bkg, bkg_std = astro.calculate_img_bkg(imgdata)
        # Detect point sources in the image.
        irafsources = astro.detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
        # If no stars are in the image, log it and go to the next one.
        if not irafsources:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No sources detected')
                f.write('\n')
            excluded_files += 1
            continue
        # Calculate the number of point sources detected.
        num_sources = len(irafsources)
        # Calculate the FWHM in pixels.
        fwhms, fwhm, fwhm_std = astro.calculate_fwhm(irafsources)
        # Do PSF photometry on the detected sources.
        
        if photometry_method=='psf':

            photometry_result = perform_photometry.perform_PSF_photometry(
                irafsources, fwhm, imgdata, bkg=bkg,filepath=filepath,hdr=hdr)

# =============================================================================
# Experimental perform_photometry
#             photometry_result = perform_photometry.perform_PSF_photometry_2(
#                              irafsources, fwhm, imgdata, bkg=bkg,filepath=filepath,hdr=hdr)
# =============================================================================
            # Store the flux and uncertainty of the stars in a separate variable.
            fluxes = np.array(photometry_result['flux_fit'])
            fluxes_unc = np.array(photometry_result['flux_unc'])
            # Convert the flux and uncertainty to magnitude and its uncertainty.
            instr_mags = astro.calculate_magnitudes(fluxes, exptime)
            instr_mags_sigma, snr = astro.calculate_magnitudes_sigma(
                fluxes,fluxes_unc, exptime)
            
            
        elif photometry_method=='aperture':
            
            # Perform Aperture Photometry
            photometry_result=perform_photometry.perform_aperture_photometry(irafsources,fwhms,imgdata,bkg=bkg,bkg_std=np.ones(np.shape(imgdata))*bkg_std,hdr=hdr,filepath=filepath,aperture_estimation_mode=aperture_estimaiton_mode)
            
            #Re-arrange values to align with PSF Fitting standard
            fluxes_unc=np.transpose(np.array(photometry_result['flux_unc']))
            fluxes=np.array(photometry_result['flux_fit'])
            instr_mags=astro.calculate_magnitudes(fluxes,exptime)
            instr_mags_sigma, snr = astro.calculate_magnitudes_sigma(
                fluxes,fluxes_unc, exptime)
        
        
        
        
        # Read the World Coordinate System transformation added to the
        # fits header by a plate solving software (external to this program, e.g. PinPoint).
        wcs = WCS(hdr)
        # Convert the stars' (x,y) centroid locations to (RA,dec).
        skypositions = astro.convert_pixel_to_ra_dec(irafsources, wcs)
        try:
            # altazpositions = astro.convert_ra_dec_to_alt_az(skypositions, hdr,
            # lat_key='OBSGEO-B', lon_key='OBSGEO-L',
            #                                           elev_key='OBSGEO-H')

            # Convert the stars' (RA,dec) location to (Azimuth,Elevation).
            altazpositions = astro.convert_ra_dec_to_alt_az(skypositions, hdr,
                                                      lat_key=lat_key,
                                                      lon_key=lon_key,
                                                      elev_key=elev_key)
        # If there's no plate solution, it will raise and AttributeError.
        # This catches that error, logs it, and moves onto the next image.
        except AttributeError as e:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No plate solution found')
                f.write('\n')
                excluded_files += 1
            continue
        # Convert the FWHM from pixels to arcsec.
        fwhms_arcsec, fwhm_arcsec, fwhms_arcsec_std = astro.convert_fwhm_to_arcsec(
            hdr, fwhms, fwhm, fwhm_std)
        # If it can't convert from pixels to arcsec (e.g. the focal length
        # wasn't defined in the header),
        # store it as NaN.
        if not fwhm_arcsec:
            fwhm_arcsec = np.nan
            fwhms_arcsec_std = np.nan
        # Read the time from the fits header and then convert it to
        # Julian Date.
        t = Time(hdr['DATE-OBS'], format='fits', scale='utc')
        time = t.jd
        # Store the filter used to take the image as a variable.
        img_filter = hdr['FILTER']
        # Calculate the background sky brightness and standard deviation.
        background_sky_brightness = astro.calculate_background_sky_brightness(
            bkg, hdr, exptime)
        background_sky_brightness_sigma = astro.calculate_BSB_sigma(
            bkg, bkg_std, exptime)
        # Take the average of all stars' Az/El/airmass and store as a variable.
        azimuth = np.mean(altazpositions.az)
        elevation = np.mean(altazpositions.alt)
        airmass = np.mean(altazpositions.secz)
        # Update the table with auxiliary data on the images (FWHM, BSB, etc.)
        star_aux_table_columns =\
            astro.update_star_aux_columns(star_aux_table_columns,
                                    file_names[file_num],
                                    time,
                                    img_filter,
                                    fwhm,
                                    fwhm_std,
                                    fwhm_arcsec,
                                    fwhms_arcsec_std,
                                    num_sources,
                                    background_sky_brightness,
                                    background_sky_brightness_sigma,
                                    azimuth,
                                    elevation,
                                    airmass)
        # Match the detected sources with a star from the reference stars file.
        matched_stars = astro.find_ref_stars(reference_stars,
                                       ref_star_positions,
                                       skypositions,
                                       instr_mags,
                                       instr_mags_sigma,
                                       snr,
                                       fluxes,
                                       fluxes_unc,
                                       ground_based=True,
                                       altazpositions=altazpositions)
        # If no image star corresponds to a reference star, log it and go to
        # the next image.
        if not matched_stars:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No reference star found in image')
                f.write('\n')
                excluded_files += 1
            continue
        # Update the table that contains information on each detection of a
        # reference star.
        large_table_columns = astro.update_large_table_columns(large_table_columns,
                                                         filepath,
                                                         matched_stars,
                                                         hdr,
                                                         exptime,
                                                         ground_based=True,
                                                         name_key=name_key)
    # Complete the text file that stores information on files that were not used to calculate the transforms.
    with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
        f.write('Total excluded:')
        f.write('\t')
        f.write(
            f'{excluded_files} / {split_filecount_location} ({100*(excluded_files/split_filecount_location):.1f}%)')
    # Create an AstroPy table of the auxiliary data and write it to a .csv file.
    star_aux_table = astro.create_star_aux_table(star_aux_table_columns)
    ascii.write(star_aux_table, os.path.join(
        save_loc, 'auxiliary_table.csv'), format='csv')
    # Create an AstroPy table of each reference star detection and write it to a .csv file.
    large_stars_table = astro.create_large_stars_table(
        large_table_columns, ground_based=True)
    ascii.write(large_stars_table, os.path.join(
        save_loc, 'large_stars_table.csv'), format='csv')
    # Group each observation of a star at an airmass.
    # E.g. if there are 5 images of star X at 1.2 airmass, and 10 images of
    # star X at 2 airmass,
    # it will produce a mean and standard deviation of the observations at
    # both 1.2 and 2 airmass.
    # This creates that table, stores the different filters used to take the
    # images (e.g. BVRI or BGR),
    # and writes it to a .csv file.
    stars_table, different_filter_list = astro.group_each_star_GB(large_stars_table)
    ascii.write(stars_table, os.path.join(
        save_loc, 'stars_table.csv'), format='csv')

    ############# Begin the Buchheim Transforms #############

    # Calculate the slope of each star's instrumental magnitude vs airmass,
    # store it in a table,
    # and write it to a .csv file.
    slopes_table = buch_aux.calculate_slopes_Buchheim(
        stars_table, different_filter_list, save_plots, save_loc=save_loc)
    ascii.write(slopes_table, os.path.join(
        save_loc, 'slopes_table.csv'), format='csv')

    # Calculate the first and second order extinctions.
    buch_aux.second_order_extinction_calc_Buchheim(stars_table, different_filter_list,
                                          save_plots,
                                          save_loc=save_loc)
    extinction_table_Buckhheim = \
        buch_aux.extinction_calc_Buchheim_sect6(slopes_table, different_filter_list,
                                       save_plots, save_loc=save_loc)
    ascii.write(extinction_table_Buckhheim, os.path.join(
        save_loc, 'extinction_table_Buckhheim.csv'), format='csv')
    return
    # Calculate the exoatmospheric magnitudes (m_0).
    exoatmospheric_table = buch_aux.exoatmospheric_mags_Warner(
        stars_table, extinction_table_Buckhheim, different_filter_list)

    ############# Begin the Warner Transforms #############

    # Finish the transform by calculating the colour transform and zero point.
    Buckhheim_final_transform_table =\
        warn_aux.colour_transform_and_zp_calc_Warner(exoatmospheric_table,
                                            different_filter_list,
                                            extinction_table_Buckhheim,
                                            save_plots, save_loc=save_loc)
    # Save the transform table to a .csv file.
    ascii.write(Buckhheim_final_transform_table, os.path.join(
        save_loc, '_gb_final_transforms.csv'), format='csv')
    # Calculate the hidden transform and write it to a .csv file.
    hidden_transform_table =\
        warn_aux.hidden_transform_Warner(exoatmospheric_table,
                                Buckhheim_final_transform_table,
                                different_filter_list,
                                save_plots,
                                save_loc=save_loc)
    ascii.write(hidden_transform_table, os.path.join(
        save_loc, 'hidden_transform_table.csv'), format='csv')
    exoatmospheric_table_verify =\
        warn_aux.exoatmospheric_mags_verify_Warner(stars_table,
                                          extinction_table_Buckhheim,
                                          hidden_transform_table,
                                          different_filter_list)
    # Verify the transforms.
    verify_save_loc = os.path.join(save_loc, 'Verification')
    app_mag_table = astro.verify_gb_transforms_auto(directory,
                                              verification_files,
                                              ref_stars_file,
                                              Buckhheim_final_transform_table,
                                              hidden_transform_table,
                                              plot_results=True,
                                              save_plots=True,
                                              file_suffix=file_suffix,
                                              exposure_key=exposure_key,
                                              name_key=name_key,
                                              lat_key=lat_key,
                                              lon_key=lon_key,
                                              elev_key=elev_key,
                                              save_loc=verify_save_loc)
    return Buckhheim_final_transform_table

def _main_sb_transform_calc(directory,
                            ref_stars_file,
                            plot_results=False,
                            save_plots=False,
                            file_suffix=(".fits", ".fit", ".fts"),
                            exposure_key='EXPTIME',
                            name_key='Name',
                            transform_index_list=['(B-V)', '(V-R)', '(V-I)'],
                            **kwargs):
    # TODO: Docstring.
    reference_stars, ref_star_positions = astro.read_ref_stars(ref_stars_file)
    large_table_columns = astro.init_large_table_columns()

    if save_plots:
        save_loc = kwargs.get('save_loc')
        unique_id = kwargs.get('unique_id')

    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(file_suffix):
                filepath = os.path.join(dirpath, filename)
                hdr, imgdata = astro.read_fits_file(filepath)
                exptime = hdr[exposure_key]
                bkg, bkg_std = astro.calculate_img_bkg(imgdata)
                irafsources = astro.detecting_stars(
                    imgdata, bkg=bkg, bkg_std=bkg_std)
                if not irafsources:
                    continue
                _, fwhm, fwhm_std = astro.calculate_fwhm(irafsources)
                
                
                
                photometry_result = perform_photometry.perform_PSF_photometry(
                    irafsources, fwhm, imgdata, bkg=bkg)
                fluxes = np.array(photometry_result['flux_fit'])
                fluxes_unc = np.array(photometry_result['flux_unc'])
                instr_mags = astro.calculate_magnitudes(fluxes, exptime)
                instr_mags_sigma, snr = astro.calculate_magnitudes_sigma(
                    fluxes,fluxes_unc, exptime)
                wcs = WCS(hdr)
                skypositions = astro.convert_pixel_to_ra_dec(irafsources, wcs)
                matched_stars = astro.find_ref_stars(reference_stars,
                                               ref_star_positions,
                                               skypositions,
                                               instr_mags,
                                               instr_mags_sigma,
                                               snr,
                                               fluxes,
                                               fluxes_unc,
                                               ground_based=False,
                                               altazpositions=None)
                if not matched_stars:
                    continue

                large_table_columns = astro.update_large_table_columns(large_table_columns,
                                                                 filepath,
                                                                 matched_stars,
                                                                 hdr,
                                                                 exptime,
                                                                 ground_based=False,
                                                                 name_key=name_key)
    large_stars_table = astro.create_large_stars_table(
        large_table_columns, ground_based=False)
    stars_table = astro.group_each_star(large_stars_table, ground_based=False)
    sb_final_transform_columns = auxillary_phot_sb_functions.init_sb_final_transform_columns()
    if save_plots:
        astro.write_table_to_latex(stars_table, f"{os.path.join(save_loc, f'{unique_id}_stars_table')}.txt",
                             formats={'c': '%0.3f',
                                      'c_sigma': '%0.3f'})
        for index in transform_index_list:
            filter_fci, filter_fci_sigma, zprime_fci, zprime_fci_sigma = space_based_transform(stars_table,
                                                                                               plot_results=plot_results,
                                                                                               index=index,
                                                                                               save_plots=save_plots,
                                                                                               save_loc=save_loc,
                                                                                               unique_id=unique_id)
            sb_final_transform_columns = auxillary_phot_sb_functions.update_sb_final_transform_columns(sb_final_transform_columns,
                                                                                                       index,
                                                                                                       filter_fci,
                                                                                                       filter_fci_sigma,
                                                                                                       zprime_fci,
                                                                                                       zprime_fci_sigma)
            # print(f"(V-clear) = ({filter_fci:.3f} +/- {filter_fci_sigma:.3f}) * {index} + " \
            #       f"({zprime_fci:.3f} +/- {zprime_fci_sigma:.3f})")
    else:
        for index in transform_index_list:
            filter_fci, filter_fci_sigma, zprime_fci, zprime_fci_sigma = space_based_transform(stars_table,
                                                                                               plot_results=plot_results,
                                                                                               index=index,
                                                                                               save_plots=save_plots)
            sb_final_transform_columns = auxillary_phot_sb_functions.update_sb_final_transform_columns(sb_final_transform_columns,
                                                                                                       index,
                                                                                                       filter_fci,
                                                                                                       filter_fci_sigma,
                                                                                                       zprime_fci,
                                                                                                       zprime_fci_sigma)
            # print(f"(V-clear) = ({filter_fci:.3f} +/- {filter_fci_sigma:.3f}) * {index} + " \
            #       f"({zprime_fci:.3f} +/- {zprime_fci_sigma:.3f})")
    sb_final_transform_table = auxillary_phot_sb_functions.create_sb_final_transform_table(
        sb_final_transform_columns)
    if save_plots:
        formats = {
            'T_fCI': '%0.3f',
            'T_fCI_sigma': '%0.3f',
            'Z_fCI': '%0.3f',
            'Z_fCI_sigma': '%0.3f'
        }
        astro.write_table_to_latex(sb_final_transform_table, f"{os.path.join(save_loc, f'{unique_id}_transform_table')}.txt",
                             formats=formats)
    return sb_final_transform_table

def space_based_transform(stars_table,
                          plot_results=False,
                          save_plots=False,
                          index='(B-V)',
                          app_filter='V',
                          instr_filter='c',
                          field=None,
                          **kwargs):
    """
    Calculate the transforms for a space based sensor.

    Parameters
    ----------
    stars_table : astropy.table.table.Table
        Table containing the mean of the important information for each star.
        Has columns:
            Field : string
                Unique identifier of the star field that the reference star is
                in (e.g. Landolt field "108").
            Name : string
                Name/unique identifier of the reference star.
            V : numpy.float64
                Apparent V magnitude from the reference file.
            (B-V) : numpy.float64
                Apparent B-V colour index from the reference file.
            (U-B) : numpy.float64
                Apparent U-B colour index from the reference file.
            (V-R) : numpy.float64
                Apparent V-R colour index from the reference file.
            (V-I) : numpy.float64
                Apparent V-I colour index from the reference file.
            V_sigma : numpy.float64
                Standard deviation of the apparent V magnitude from the
                reference file.
            <filter> : numpy.float64
                Mean instrumental magnitude of all detections of the star in
                <filter>. There is a different column for 
                each different filter used across the images.
            <filter>_sigma : numpy.float64
                Standard deviation of the instrumental magnitudes of all
                detections of the star in <filter>. 
                There is a different column for each different filter used
                across the images.
            X_<filter> : numpy.float64
                Mean airmass of all detections of the star in <filter>.
                There is a different column for each different 
                filter used across the images. Only output if ground_based
                is True.
            X_<filter>_sigma : numpy.float64
                Standard deviation of the airmasses of all detections of
                the star in <filter>. There is a different 
                column for each different filter used across the images.
                Only output if ground_based is True.
    plot_results : bool, optional
        Controls whether or not to plot the results from the transforms.
        The default is False.
    index : string, optional
        Colour index to calculate the transform for. The default is '(B-V)'.
    app_filter : string, optional
        Apparent magnitude filter band to calculate the transform for.
        The default is 'V'.
    instr_filter : string, optional
        Instrumental filter band to calculate the transform for.
        The default is 'clear'.
    field : string, optional
        Unique identifier of the star field that the reference star is
        in (e.g. Landolt field "108"). 
        The default is None.

    Returns
    -------
    filter_fci : float
        Instrumental transform coefficient for filter f using the colour
        index CI.
    zprime_fci : float
        Zero point magnitude for filter f.

    """
    max_app_filter_sigma = max(stars_table[f'{app_filter}_sigma'])
    max_instr_filter_sigma = max(stars_table[f'{instr_filter}_sigma'])
    err_sum = np.nan_to_num(stars_table[f'{app_filter}_sigma'],
                            nan=max_app_filter_sigma) + \
        np.nan_to_num(
            stars_table[f'{instr_filter}_sigma'], nan=max_instr_filter_sigma)
    err_sum = np.array(err_sum)
    err_sum[err_sum == 0] = max(err_sum)

    x = stars_table[index]
    y = stars_table[f'{app_filter}_ref'] - stars_table[instr_filter]
    fit, or_fit, line_init = astro.init_linear_fitting(sigma=2.5)
    fitted_line, mask = or_fit(line_init, x, y, weights=1.0 / err_sum)
    filtered_data = np.ma.masked_array(y, mask=mask)
    filter_fci = fitted_line.slope.value
    zprime_fci = fitted_line.intercept.value
    cov = fit.fit_info['param_cov']

    # a_fit, cov = curve_fit(linear_func, stars_table[index],
    #                          stars_table[f'{app_filter}_ref'] -
    # stars_table[instr_filter],
    #                          sigma=err_sum)
    # filter_fci = a_fit[0]
    filter_fci_sigma = sqrt(cov[0][0])
    # zprime_fci = a_fit[1]
    zprime_fci_sigma = sqrt(cov[1][1])
    if plot_results:
        index_plot = np.arange(start=min(stars_table[index]), stop=max(
            stars_table[index]) + 0.01, step=0.01)
        plt.errorbar(x, y, yerr=err_sum, color='#1f77b4', fmt='o',
                     fillstyle='none', capsize=2, label="Clipped Data")
        plt.plot(x, filtered_data, 'o', color='#1f77b4', label="Fitted Data")
        plt.plot(index_plot, fitted_line(index_plot), '-', color='#ff7f0e',
                 label=f"({app_filter}-{instr_filter}) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")

        # plt.errorbar(stars_table[index], stars_table[f'{app_filter}_ref'] -\
        # stars_table[instr_filter],
        #              yerr=err_sum, fmt='o', capsize=2)
        # plt.plot(index_plot, filter_fci * index_plot + zprime_fci,
        #          label=f"({app_filter}-{instr_filter}) =\
        # {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        plt.ylabel(f"{app_filter}-{instr_filter}")
        plt.xlabel(f"{index}")
        plt.legend()
        plt.title(f"Space Based Transform Calculation for {index}")
        if save_plots:
            unique_id = kwargs.get('unique_id')
            save_loc = f"{os.path.join(kwargs.get('save_loc'), f'{unique_id}_TZfci_{index}')}.png"
            plt.savefig(save_loc)
        # if not field:
        #     plt.title(f"({app_filter}-{instr_filter}) =\
        # {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        # else:
        #     plt.title(f"{field}: ({app_filter}-{instr_filter}) =\
        # {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        plt.show()
        plt.close()
    return filter_fci, filter_fci_sigma, zprime_fci, zprime_fci_sigma
