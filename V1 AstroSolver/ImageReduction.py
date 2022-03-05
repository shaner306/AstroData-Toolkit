# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 16:22:33 2021

@author: shane
"""
from pathlib import Path
from astropy.nddata import CCDData
from astropy.visualization import hist
from astropy.io import fits
from astropy import units as u
from collections import Counter
import ccdproc as ccdp
import matplotlib.pyplot as plt
import numpy as np
import os
from ccdproc import ImageFileCollection
import time


# #-------------------------------------------------------------------------------
# #
# #   Initialization variables
# #
# #-------------------------------------------------------------------------------

# #  The following variables are used as switches to run the main functions.
# #
# #  0 = Don't run function
# #  1 = Run function
# #
# create_master_dir = 1
# run_master_bias = 0
# run_master_dark = 1
# run_master_flat = 1
# correct_light_frames = 1

# # The highest directoy of the .fits files to process
# # topdir = 'C:\\Astro_Data\\'

# topdir = 'D:\Image Reduction Test Images'

# if os.path.isdir(topdir) == False:
#     raise RuntimeError('WARNING -- Directory of .fits files does not exist')

# # Create output directory for master files
# if create_master_dir == True:
#     master_frame_dir = Path(topdir, 'master_frame_data')
#     master_frame_dir.mkdir(exist_ok = True)

# #  Select directory for master frames
# master_frame_directory = topdir + '\master_frame_data'

# #  If a directory already exists containing the master files, uncomment the
# #  following line and place the path as a string with double backslashes.

# #master_frame_directory = 'C:\\pineapple\\is_a_fruit'

# if os.path.isdir(master_frame_directory) == False:
#     raise RuntimeError('WARNING -- Directory of Master .fits files does not exist')


# # The extension to search for
# exten = '.fits'

# # The list of *.fits files created by the directory walk
# results = []


# -------------------------------------------------------------------------------
#
#   Functions
#
# -------------------------------------------------------------------------------

def create_master_bias(all_fits, master_dir):
    """
    Taking all Bias frames, create a Master Bias frame.

    Parameters
    ----------

    all_fits : Astropy.CCDProc.ImageFileCollection object
        This object contains all .fits files in the directory and sub-directories
        in 'topdir' given by the user in the initialization section.


    master_dir : string
        The output directory to save the Master Bias frame into.

    Returns
    -------

    Nil

    """
    unique_imagetype_list = list(set(all_fits.summary['imagetyp']))
    bias_imgtype_matches = [
        s for s in unique_imagetype_list if "bias" in s.lower()]
    bias_imgtypes_concatenateded = '|'.join(bias_imgtype_matches)
    bias_fits = all_fits.files_filtered(
        include_path=True,
        imagetyp=bias_imgtypes_concatenateded)

    biases = []
    for file in bias_fits:
        bias = CCDData.read(file, unit='adu')
        biases.append(bias)

    master_bias = ccdp.combine(biases,
                               method='average',
                               sigma_clip=True,
                               sigma_clip_low_thresh=5,
                               sigma_clip_high_thresh=5,
                               sigma_clip_func=np.ma.median,
                               sigma_clip_dev_func=np.ma.std,
                               mem_limit=1e9
                               )

    master_bias.meta['combined'] = True
    bias_file_name = '\master_bias.fits'
    master_bias.write(master_dir + bias_file_name)
    print('Saving ', bias_file_name)


def create_master_dark(all_fits, master_dir):
    """
    Taking a Master Bias frame and individual Dark frames, create Master Dark
    frames based on differing exposure times. Master Bias is subtracted from
    each Dark frame and then a Master Dark is created from a combination of
    similarly exposed Darks.

    Parameters
    ----------

    all_fits : Astropy.CCDProc.ImageFileCollection object
        This object contains all .fits files in the directory and sub-directories
        in 'topdir' given by the user in the initialization section.


    master_dir : string
        The output directory to save the Master Dark frames into.

    Returns
    -------

    Nil

    """

# master_files = ccdp.ImageFileCollection(master_dir)

    try:
        # combined_bias = CCDData.read(master_files.files_filtered(imagetyp ='bias',
        # combined=True))
        master_bias = CCDData.read(master_dir + '\\master_bias.fits', unit='adu')
    except:
        raise RuntimeError('WARNING -- Could not open Master Bias file')
    # Dynamic Dark Imagetype Reader
    unique_imagetype_list = list(set(all_fits.summary['imagetyp']))
    dark_imgtype_matches = [
        s for s in unique_imagetype_list if "dark" in s.lower()]
    dark_imgtypes_concatenateded = '|'.join(dark_imgtype_matches)

    # Enables different datasets with different dark headers

    # Create Conditonal Array with different Dark Image Data Types
    dark_mask = list(np.zeros(len(all_fits.summary), dtype=bool))
    for dark_imgtypes in dark_imgtype_matches:
        dark_mask = dark_mask +\
            (all_fits.summary['imagetyp'] == (dark_imgtypes))

    dark_times = set(all_fits.summary['exptime'][dark_mask])

    for exp_time in sorted(dark_times):

        darks_to_calibrate = all_fits.files_filtered(
            imagetyp=dark_imgtypes_concatenateded,
            exptime=exp_time,
            include_path=True)

        darks = []
        for file in darks_to_calibrate:
            dark = CCDData.read(file, unit='adu')
            dark = ccdp.subtract_bias(dark, master_bias)
            darks.append(dark)

        master_dark = ccdp.combine(darks,
                                   method='average',
                                   sigma_clip=True,
                                   sigma_clip_low_thresh=5,
                                   sigma_clip_high_thresh=5,
                                   sigma_clip_func=np.ma.median,
                                   signma_clip_dev_func=np.ma.std,
                                   mem_limit=1e9,
                                   )

        master_dark.meta['combined'] = True
        dark_file_name = '\\master_dark_{:6.3f}.fits'.format(exp_time)
        master_dark.write(master_dir + dark_file_name)
        print('Saving ', dark_file_name)


def create_master_flat(all_fits, master_dir):
    """
    Taking Master Bias, Master Darks and individual Flat frames, create
    Master Flat frames based on the filters used.  Master Bias is subtracted
    from each Flat frame and then Master Dark is scaled for differing exposure
    time, and then subtracted.  Corrected Flats are then grouped by filter type
    and then combined.

    Parameters
    ----------

    all_fits : Astropy.CCDProc.ImageFileCollection object
        This object contains all .fits files in the directory and sub-directories
        in 'topdir' given by the user in the initialization section.


    master_dir : string
        The output directory to save the Master Flat frames into.

    Returns
    -------

    Nil

    """

    unique_imagetype_list = list(set(all_fits.summary['imagetyp']))
    bias_imgtype_matches = [
        s for s in unique_imagetype_list if "bias" in s.lower()]
    bias_imgtypes_concatenateded = '|'.join(bias_imgtype_matches)
    dark_imgtype_matches = [
        s for s in unique_imagetype_list if "dark" in s.lower()]
    dark_imgtypes_concatenateded = '|'.join(dark_imgtype_matches)
    flat_imgtype_matches = [
        s for s in unique_imagetype_list if "flat" in s.lower()]
    flat_imgtypes_concatenateded = '|'.join(flat_imgtype_matches)

    master_files = ccdp.ImageFileCollection(master_dir)

    master_files_list = master_files.files_filtered(
        imagetyp=bias_imgtypes_concatenateded,
        combined=True,
        include_path=True
    )
    try:
        # combined_bias = CCDData.read(master_files_list)
        combined_bias = CCDData.read(master_dir + '\\master_bias.fits', unit='adu')
    except:
        raise RuntimeError('WARNING -- Could not open Master Bias file')

    try:
        # Get the Master Darks
        master_darks = {ccd.header['exposure']: ccd for ccd in master_files.ccds(
            imagetyp=dark_imgtypes_concatenateded, combined=True
        )}
        # Create Conditonal Array with different Dark Image Data Types
        dark_mask = list(np.zeros(len(all_fits.summary), dtype=bool))
        for dark_imgtypes in dark_imgtype_matches:
            dark_mask = dark_mask +\
                (all_fits.summary['imagetyp'] == (dark_imgtypes))

        dark_times = set(all_fits.summary['exptime'][dark_mask])

    except:
        raise RuntimeError('WARNING -- Could not open Master Dark files')

    if len(master_darks) == 0:
        raise RuntimeError('WARNING -- Could not open Master Dark files')

    flat_mask = list(np.zeros(len(all_fits.summary), dtype=bool))
    for flat_imgtypes in flat_imgtype_matches:
        flat_mask = flat_mask +\
            (all_fits.summary['imagetyp'] == (flat_imgtypes))

    flat_filters = set(all_fits.summary['filter'][flat_mask])

    for frame_filter in sorted(flat_filters):
        flats_to_combine = all_fits.files_filtered(
            imagetyp=flat_imgtypes_concatenateded,
            filter=frame_filter,
            include_path=True
        )

        flats = []
        for file in flats_to_combine:
            flat = CCDData.read(file, unit='adu')
            flat = ccdp.subtract_bias(flat, combined_bias)
            closest_dark = find_nearest_dark_exposure(flat, dark_times,
                                                      tolerance=100
                                                      )
            flat = ccdp.subtract_dark(flat, master_darks[closest_dark],
                                      exposure_time='exptime',
                                      exposure_unit=u.second,
                                      scale=True
                                      )
            flats.append(flat)

        combined_flat = ccdp.combine(flats,
                                     method='median',
                                     sigma_clip=True,
                                     sigma_clip_low_thresh=5,
                                     sigma_clip_high_thresh=5,
                                     sigma_clip_func=np.ma.median,
                                     signma_clip_dev_func=np.ma.std,
                                     mem_limit=1e9
                                     )

        combined_flat.meta['combined'] = True
        flat_file_name = '\\master_flat_filter_{}.fits'.format(
            frame_filter.replace("''", "p"))
        combined_flat.write(master_dir + flat_file_name)
        print('Saving ', flat_file_name)


def correct_lights(all_fits, master_dir, corrected_light_dir, correct_outliers_params):
    """
    Taking Master Bias, Master Darks and filter-specific Master Flats, create
    corrected individual light frames.


    Parameters
    ----------

    all_fits : Astropy.CCDProc.ImageFileCollection object
        This object contains all .fits files in the directory and sub-directories
        in 'topdir' given by the user in the initialization section.


    master_dir : string
        The output directory to get the Master frames from.

    Returns
    -------

    Nil

    """
    unique_imagetype_list = list(set(all_fits.summary['imagetyp']))
    bias_imgtype_matches = [
        s for s in unique_imagetype_list if "bias" in s.lower()]
    bias_imgtypes_concatenateded = '|'.join(bias_imgtype_matches)
    dark_imgtype_matches = [
        s for s in unique_imagetype_list if ("dark" in s.lower())]
    dark_imgtypes_concatenateded = '|'.join(dark_imgtype_matches)
    flat_imgtype_matches = [
        s for s in unique_imagetype_list if "flat" in s.lower()]
    flat_imgtypes_concatenateded = '|'.join(flat_imgtype_matches)
    light_imgtype_matches = [
        s for s in unique_imagetype_list if "light" in s.lower()]
    light_imgtypes_concatenateded = '|'.join(light_imgtype_matches)

    master_files = ccdp.ImageFileCollection(master_dir)

    # Create dictionaries of the dark and flat master frames in the master directory
    # TODO: Double Check Concatenation
    master_darks = {ccd.header['exposure']: ccd for ccd in master_files.ccds(
        imagetyp=dark_imgtypes_concatenateded, combined=True)}
    master_flats = {ccd.header['filter']: ccd for ccd in master_files.ccds(
        imagetyp=flat_imgtypes_concatenateded, combined=True)}
    # There is only one bias frame, so no need to set up a dictionary.
    master_bias = [ccd for ccd in master_files.ccds(
        imagetyp=bias_imgtypes_concatenateded, combined=True)][0]
    light_mask = list(np.zeros(len(all_fits.summary), dtype=bool))

    for light_imgtypes in light_imgtype_matches:
        light_mask = light_mask +\
            (all_fits.summary['imagetyp'] == (light_imgtypes))

    light_filter = set(all_fits.summary['filter'][light_mask])

    dark_mask = list(np.zeros(len(all_fits.summary), dtype=bool))
    for dark_imgtypes in dark_imgtype_matches:
        dark_mask = dark_mask +\
            (all_fits.summary['imagetyp'] == (dark_imgtypes))
            
    dark_times = set(all_fits.summary['exptime'][dark_mask])
    
    
    example_light = all_fits.files_filtered(imagetyp=light_imgtypes_concatenateded,
                                                filter=list(light_filter)[0],
                                                include_path=True)[0]
    
    for dark_time in dark_times:
        ### Dark Image Masking ###
        mask = np.zeros(CCDData.read(example_light, unit='adu').data.shape, dtype=bool)
        hot_pixels = np.zeros(CCDData.read(example_light, unit='adu').data.shape, dtype=bool)
        dark_threshold_mask = np.zeros(CCDData.read(example_light
            , unit='adu').data.shape, dtype=bool)
        
        
        if correct_outliers_params['Hot Pixel'] is True:
        
            # Calculate Dark Current to find hot pixels
            dark_current = master_darks[dark_time].multiply(
                float(master_darks[dark_time].header['EGAIN'])*u.electron /
                u.adu).divide(float(master_darks[dark_time].header['EXPTIME'])*u.second)
            hot_pixels = dark_current > 4*np.nanmean(dark_current)
            mask = mask | hot_pixels
        
        if correct_outliers_params['Dark Frame Threshold Bool']:
        
            dark_pix_over_max = master_darks[dark_time].data > float(
                correct_outliers_params['Dark Frame Threshold Max'])
            dark_pix_under_min = master_darks[dark_time].data < float(
                correct_outliers_params['Dark Frame Threshold Min'])
            dark_threshold_mask = dark_pix_over_max | dark_pix_under_min
            mask = mask | dark_threshold_mask
            
            correct_outlier_darks(correct_outliers_params, hot_pixels, dark_threshold_mask,
                                  master_darks, dark_time, master_dir)

    for frame_filter in sorted(light_filter):

        lights_to_correct = all_fits.files_filtered(
            imagetyp=light_imgtypes_concatenateded,
            filter=frame_filter,
            include_path=True
        )
        light_ccds = []
        try:
            if correct_outliers_params['Outlier Boolean'] is True:
                

                ### Flat Image Masking ###
                if correct_outliers_params['ccdmask']:
                    # Based on IRAF ccdmask
                    flats_to_compare = all_fits.files_filtered(
                        imagetyp=flat_imgtypes_concatenateded,
                        filter=frame_filter,
                        include_path=True
                    )
                    # TODO: Add functionality for single flat frame
                    bad_ratio = False
                    if len(flats_to_compare) > 1:
                        # Creates Dictionary Object where key is the mean of the calibrated frame data
                        calibrated_flats = {}
                        # Calibrate Flats
                        for flat in flats_to_compare:
                            flatdata = CCDData.read(flat, unit='adu')
                            flatdata = ccdp.subtract_bias(flatdata, master_bias)
                            closest_dark = find_nearest_dark_exposure(CCDData.read(lights_to_correct[0], unit='adu'), dark_times,
                                                                      tolerance=100
                                                                      )

                            ########## For Scalable darks (i.e. Bias Frame provided) ##########
                            flatdata = ccdp.subtract_dark(flatdata, master_darks[closest_dark],
                                                          exposure_time='exptime',
                                                          exposure_unit=u.second,
                                                          scale=True
                                                          )

                            calibrated_flats[flatdata.data.mean()] = flatdata

                        max_key = max(calibrated_flats.keys())
                        min_key = min(calibrated_flats.keys())
                        ratio = calibrated_flats[max_key].divide(calibrated_flats[min_key])
                        if (np.nanmean(ratio) < 1.1) and (np.nanmean(ratio) > 0.9):  # ratio thresholds
                            bad_ratio = True
                        else:
                            maskr = ccdp.ccdmask(ratio)

                            if correct_outliers_params['Replace Bool']:
                                # Replaces the values in mask with a desired value
                                if correct_outliers_params['Replace Mode'] == 'Ave':
                                    print('Calculate Average Background')

                                    print('Save New Master Flat')

                                elif correct_outliers_params['Replace Mode'] == 'Interpolate':
                                    print('Interpolation')

                    if (len(flats_to_compare) == 1) or (bad_ratio is True):
                        # Calibrate flat field
                        flatdata = CCDData.read(flats_to_compare[0], unit='adu')
                        flatdata = ccdp.subtract_bias(flatdata, master_bias)
                        closest_dark = find_nearest_dark_exposure(CCDData.read(lights_to_correct[0], unit='adu'), dark_times,
                                                                  tolerance=100
                                                                  )

                        ########## For Scalable darks (i.e. Bias Frame provided) ##########
                        flatdata = ccdp.subtract_dark(flatdata, master_darks[closest_dark],
                                                      exposure_time='exptime',
                                                      exposure_unit=u.second,
                                                      scale=True
                                                      )
                        maskr = ccdp.ccdmask(flatdata)

                        correct_outlier_flats(correct_outliers_params,
                                              maskr,
                                              flats_to_compare,
                                              frame_filter,
                                              master_dir)

                    mask = mask | maskr

                    
            for file_name in lights_to_correct:

                light = CCDData.read(file_name, unit='adu')

                # Note that the first argument in the remainder of the ccdproc calls is
                # the *reduced* image, so that the calibration steps are cumulative.

                ########## For Scalable darks (i.e. Bias Frame provided) ##########
                reduced = ccdp.subtract_bias(light, master_bias)

                closest_dark = find_nearest_dark_exposure(
                    reduced, master_darks.keys())

                reduced = ccdp.subtract_dark(reduced, master_darks[closest_dark],
                                             exposure_time='exptime', exposure_unit=u.second,
                                             scale=True)
                ########## For non-scalable darks (i.e. no Bias Fram provided) ##########
                #### You also have to set run_master_bias to 0 on line 378 of Main.py ####
                # closest_dark = find_nearest_dark_exposure(light, master_darks.keys())

                # reduced = ccdp.subtract_dark(light, master_darks[closest_dark],
                #                               exposure_time='exptime', exposure_unit=u.second,
                #                               scale=False)

                ########## End ##########

                # Flat Correction
                good_flat = master_flats[reduced.header['filter']]
                redcued = ccdp.flat_correct(reduced, good_flat)

                ### Cosmic Rays Outliers ###
                if correct_outliers_params['Outlier Boolean'] is True:

                    if correct_outliers_params['Cosmic Rays Bool']:
                        # Convert image to Electrons
                        reduced_in_e = ccdp.gain_correct(reduced, float(light.header['EGAIN'])*u.electron/u.adu)
                        reduced_in_e.mask = mask
                        new_reduced_in_e = ccdp.cosmicray_lacosmic(reduced_in_e, readnoise=10, sigclip=5, verbose=True)
                        reduced = ccdp.gain_correct(new_reduced_in_e, (u.adu/(float(light.header['EGAIN'])*u.electron)))
                        mask = reduced_in_e.mask
                        print('Removed Cosmic Rays')

                    reduced.mask = mask
                reduced.meta['correctd'] = True
                file_name = file_name.split("\\")[-1]
                try:
                    reduced.write(corrected_light_dir + '\\' + file_name)
                except OSError:
                    file_name = file_name[:-5]
                    print(file_name)
                    file_name = file_name + "1.fits"
                    reduced.write(corrected_light_dir + '\\' + file_name)

                print('Saving ', file_name)
        except Exception as e:
            print(e)
            pass


def find_nearest_dark_exposure(image, dark_exposure_times, tolerance=0.5):
    """
    Find the nearest exposure time of a dark frame to the exposure time of the image,
    raising an error if the difference in exposure time is more than tolerance.

    Parameters
    ----------

    image : astropy.nddata.CCDData
        Image for which a matching dark is needed.

    dark_exposure_times : list
        Exposure times for which there are darks.

    tolerance : float or ``None``, optional
        Maximum difference, in seconds, between the image and the closest dark. Set
        to ``None`` to skip the tolerance test.

    Returns
    -------

    float
        Closest dark exposure time to the image.
    """

    dark_exposures = np.array(list(dark_exposure_times))
    idx = np.argmin(np.abs(dark_exposures - image.header['exptime']))
    closest_dark_exposure = dark_exposures[idx]

    # if (tolerance is not None and
    #     np.abs(image.header['exptime'] - closest_dark_exposure) > tolerance):

    # raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
    # 'time {}.'.format(closest_dark_exposure, a_flat.header['exptime']))

    return closest_dark_exposure


def correct_outlier_flats(correct_outliers_params, maskr, flats_to_compare, frame_filter, master_dir):
    if correct_outliers_params['Replace Bool']:
        # Replaces the values in mask with a desired value
        if correct_outliers_params['Replace Mode'] == 'Ave':
            print('Calculate Average Background')
            # TODO: Add Script
        elif correct_outliers_params['Replace Mode'] == 'Interpolate':
            radius = 1
            coordinates = np.where(maskr == True)

            for flat in flats_to_compare:
                flatdata = CCDData.read(flat, unit='adu')
                for i in range(0, np.shape(coordinates)[1]):

                    # Create copy of flat to add mask to - To Avoid Masked Pixs
                    flatdata2 = flatdata.copy()
                    flatdata2.mask = maskr
                    try:
                        Replaceable_mean = np.nanmean(
                            flatdata2[(coordinates[0][i]-radius):(coordinates[0][i]+radius), (coordinates[1][i]-radius):(coordinates[1][i]+radius)])
                    except:  # Except faulty pix is in the corner of the image
                        Replaceable_mean = np.nanmean(flatdata2)
                    flatdata.data[coordinates[0][i]][coordinates[1][i]] = Replaceable_mean
                # Save the Data
                if len(flats_to_compare) == 1:
                    flat_file_name = '\\master_flat_outliercorrected_filter_{}.fits'.format(
                        frame_filter.replace("''", "p"))
                    flatdata.write(master_dir + flat_file_name)
                else:

                    corrected_dir = master_dir.split('\\')[0]+'//corrected_flats'
                    if os.path.isdir(corrected_dir) is False:
                        os.mkdir(corrected_dir)
                    flat_name_dir = corrected_dir + '\\' + \
                        os.path.splitext(os.path.basename(flat))[0] + 'corrected.fits'

                    flatdata.write(flat_name_dir)
            # IF ratio of flats is bad or only one flat frame we can save the flat as the master
            print('Saved New Flats')


def correct_outlier_darks(correct_outliers_params, hot_pixels, dark_threshold_mask, master_dark, closest_dark, master_dir):
    if correct_outliers_params['Replace Bool']:
        if correct_outliers_params['Replace Mode'] == 'Ave':
            print('Calculate Average Background')
            # TODO: Add Script
        elif correct_outliers_params['Replace Mode'] == 'Interpolate':
            radius = 1
            coordinates = np.where(hot_pixels == True)
            darkdata = (master_dark[closest_dark])
            darkdata2 = darkdata.copy()  # copy since we don't want new values messing up new values
            darkdata2.mask = hot_pixels.data | dark_threshold_mask
            average_dark_value=np.nanmean(darkdata2)
            for i in range(0, np.shape(coordinates)[1]):
                try:
                    Replaceable_mean = np.nanmean(
                        darkdata2[(coordinates[0][i]-radius):(coordinates[0][i]+radius), (coordinates[1][i]-radius):(coordinates[1][i]+radius)])
                except RuntimeWarning:  # Mean of empty slice
                    Replaceable_mean = average_dark_value
                except ValueError:
                    try:
                        Replaceable_mean = float(np.mean((darkdata2[(coordinates[0][i]-radius):(coordinates[0][i]+radius), (coordinates[1][i]-radius):(coordinates[1][i]+radius)])).data)
                    
                    except RuntimeWarning:
                        Replaceable_mean = average_dark_value
                        
                if str(Replaceable_mean)=='nan':
                    Replaceable_mean=average_dark_value
                    
                darkdata.data[coordinates[0][i]][coordinates[1][i]] = Replaceable_mean     
            coordinates = np.where(dark_threshold_mask == True)
            for i in range(0, np.shape(coordinates)[1]):
                try:
                    Replaceable_mean = np.nanmean(
                        darkdata2[(coordinates[0][i]-radius):(coordinates[0][i]+radius), (coordinates[1][i]-radius):(coordinates[1][i]+radius)])
                except:  # Except faulty pix is in the corner of the image
                    Replaceable_mean = np.nanmean(darkdata2)
                darkdata.data[coordinates[0][i]][coordinates[1][i]] = Replaceable_mean

        dark_name_dir = master_dir + '\\' + 'master_dark_' + str(closest_dark) + 'corrected.fits'
        darkdata.write(dark_name_dir)
