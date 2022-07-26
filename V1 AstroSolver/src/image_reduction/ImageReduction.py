# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 16:22:33 2021

@author: shane

Image Reudction Module includes all classical Astronomy methods for 
image reduction. 


Creation of Master Frames (with Dark Scaling):

┌───────────┐    ┌───────┐   ┌─────────────────┐
│Bias Frames├───►│Combine├──►│Master Bias Frame├────────────────────────────────────┐
└───────────┘    └───────┘   └─────────────────┘                                    │
                                      ▼                                             │
┌────────────┐               ┌─────────────────┐  ┌───────┐  ┌─────────────────┐    │
│Dark Frames ├──────────────►│Sub M. Bias Frame├─►│Combine├─►│Master Dark Frame│    │
└────────────┘               └─────────────────┘  └───────┘  └─────────┬───────┘    │
                                                                       │            │
                                      ┌────────────────────────────────┘            │
                                      ▼                                             │
┌────────────┐               ┌─────────────────┐  ┌─────────────────┐               │
│Flat Frames ├──────────────►│Sub M. Dark Frame├─►│Sub M. Bias Frame│◄──────────────┘
└────────────┘               └─────────────────┘  └─────────────────┘
                                                            ▼
                                                        ┌───────┐     ┌──────────────────┐
                                                        │Combine├────►│Master Flat Frame │
                                                        └───────┘     └──────────────────┘
                                                        

Creation of Master Frames (without Dark Scaling):

┌────────────┐               ┌───────┐  ┌─────────────────┐
│Dark Frames ├──────────────►│Combine├─►│Master Dark Frame│
└────────────┘               └───────┘  └────────┬────────┘
                                                 │
                                      ┌──────────┘
                                      ▼
┌────────────┐               ┌─────────────────┐    ┌───────┐   ┌──────────────────┐
│Flat Frames ├──────────────►│Sub M. Dark Frame├───►│Combine├──►│Master Flat Frame │
└────────────┘               └─────────────────┘    └───────┘   └──────────────────┘



The Master Frames are then used to calibrate the images

             *┌──────┐  ┌──────┐  ┌──────┐
              │M.    │  │M.    │  │M.    │
              │ Bias │  │ Dark │  │ Flat │
              │ Frame│  │ Frame│  │ Frame│
┌──────────┐  └──────┘  └──────┘  └──────┘  ┌─────────┐
│          │      ▼         ▼         ▼     │         │
│          │  ┌──────┐  ┌──────┐  ┌──────┐  │ Calib.  │
│  Light   │  │Sub M.│  │Sub M.│  │Div M.│  │  Light  │
│   Frames ├─►│ Bias ├─►│ Dark ├─►│ Flat ├─►│  Frames │
│          │  └──────┘  └──────┘  └──────┘  │         │
└──────────┘                                └─────────┘
             
              * M. Bias Frame Sub only with Scaling
              

See AstroPy ccdproc

"""
import astropy.nddata
import matplotlib as mpl
import numpy.ma

mpl.use('Agg')
from astropy.nddata import CCDData
from astropy.io import fits
from astropy import units as u
import ccdproc as ccdp
import numpy as np
import os
from ccdproc import ImageFileCollection


# #-------------------------------------------------------------------------------
# #
# #   Initialization variables
# #
# #-------------------------------------------------------------------------------

# #  The following variables are used as switches to run the main general_tools.
# #
# #  0 = Don't run function
# #  1 = Run function
# #
# create_master_dir = 1
# run_master_bias = 0
# run_master_dark = 1
# run_master_flat = 1
# correct_light_frames = 1

# # The highest directory of the .fits files to process
# # topdir = 'C:\\Astro_Data\\'


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
# extension = '.fits'

# # The list of *.fits files created by the directory walk
# results = []


# -------------------------------------------------------------------------------
#
#   Functions
#
# -------------------------------------------------------------------------------

def create_master_bias(all_fits: ccdp.ImageFileCollection ,
                       master_dir: str):
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
    
    Outputs
    -------
    master_bias files being written to master_dir
    """
    unique_imagetype_list = list(set(all_fits.summary['imagetyp'])) # defines the imagetypes extracted from the headers
    bias_imgtype_matches = [
        s for s in unique_imagetype_list if "bias" in s.lower()]
    bias_imgtypes_concatenated = '|'.join(bias_imgtype_matches)
    bias_fits = ImageFileCollection(
        filenames=(all_fits.files_filtered(include_path=True, imagetyp=bias_imgtypes_concatenated)))

    try:
        unique_bin_list = list(set(bias_fits.summary['ybinning'])) # defines the number of unqiue bins
    except KeyError as e:
        raise e
        print("Could not find binning keyword")

    for binning in unique_bin_list: #iterate through bins
        try:
            bin_bias_fits = bias_fits.files_filtered(include_path=True, ybinning=binning)
            biases = []
            for file in bin_bias_fits:
                # bias = CCDData.read(file, unit='adu')

                # Work around for applying BZERO and BSCALE
                hdul = fits.open(file)
                bzero = 32768
                bscale = 1
                if 'BZERO' in hdul[0].header:
                    bzero=hdul[0].header['BZERO']
                if 'BSCALE' in hdul[0].header:
                    bscale=hdul[0].header['BSCALE']
                if str(hdul[0].data.dtype)=='uint16' and hdul[0].header['BITPIX']==16:
                    hdul[0].data = ((hdul[0].data).astype('float64'))
                    hdul[0].data = hdul[0].data*bscale + bzero

                bias = CCDData(hdul[0].data, unit='adu')
                bias.header = hdul[0].header

                hdul.close()

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
            try:
                master_bias.meta['combined'] = True
                bias_file_name = str(binning) + 'x' + str(binning) + '_master_bias.fits'
                master_bias.write(os.path.join(master_dir, bias_file_name),overwrite=True)
                print('Saving ', bias_file_name)
            except OSError as e:
                raise RuntimeError('WARNING -- Could not Save Master Bias File--'+e)
        except e:
            raise RuntimeError(e)
            continue


def create_master_dark(all_fits: ccdp.ImageFileCollection ,
                       master_dir: str,
                       scalable_dark_bool: bool):
    """
    Taking a Master Bias frame and individual Dark frames, create Master Dark
    frames based on differing exposure times. Master Bias is subtracted from
    each Dark frame and then a Master Dark is created from a combination of
    similarly exposed Darks.

    Parameters
    ----------

    scalable_dark_bool:
        Boolean to determine if the Dark frames are scalable. If True, the
        Dark frames are scaled to the same exposure time as the Master Bias.
        If False, the Dark frames are not scaled.


    all_fits : Astropy.CCDProc.ImageFileCollection object
        This object contains all .fits files in the directory and sub-directories
        in 'topdir' given by the user in the initialization section.


    master_dir : string
        The output directory to save the Master Dark frames into.

    Returns
    -------
    #TODO: Return master_dark as an object

    Output
    -------
    master_dark saved in the master_dir

    """
    bscale=1
    bzero=32768
    unique_imagetype_list = list(set(all_fits.summary['imagetyp']))
    dark_imgtype_matches = [
        s for s in unique_imagetype_list if "dark" in s.lower()]
    dark_imgtypes_concatenated = '|'.join(dark_imgtype_matches)

    # Create Conditional Array with different Dark Image Data Types
    dark_mask = list(np.zeros(len(all_fits.summary), dtype=bool))
    for dark_imgtypes in dark_imgtype_matches:
        dark_mask = dark_mask + \
                    (all_fits.summary['imagetyp'] == (dark_imgtypes))

    dark_times = set(all_fits.summary['exptime'][dark_mask])

    dark_fits = ImageFileCollection(
        filenames=(all_fits.files_filtered(include_path=True, imagetyp=dark_imgtypes_concatenated)))

    try:
        unique_bin_list = list(set(dark_fits.summary['ybinning']))
    except KeyError as e:
        raise e
        print("WARNING -- Could not find binning keyword")

    # master_files = ccdp.ImageFileCollection(master_dir)
    for binning in unique_bin_list:

        if scalable_dark_bool is True:

            try:

                # combined_bias = CCDData.read(master_files.files_filtered(imagetyp ='bias',
                # combined=True))
                bin_master_string = str(binning) + 'x' + str(binning) + '_master_bias.fits'
                master_bias = CCDData.read(os.path.join(master_dir, bin_master_string), unit='adu')
            except RuntimeError:
                raise RuntimeError('WARNING -- Could not find/open Master Bias file')
        # Dynamic Dark imageType Reader

        else:
            master_bias = None

        for exp_time in sorted(dark_times):

            darks_images_to_calibrate = all_fits.files_filtered(
                imagetyp=dark_imgtypes_concatenated,
                exptime=exp_time,
                ybinning=binning,
                include_path=True)

            darks = []
            for file in darks_images_to_calibrate:

                # dark = CCDData.read(file, unit='adu')

                # Work around for applying BZERO and BSCALE
                hdul = fits.open(file)
                bzero = 32768
                bscale = 1
                if 'BZERO' in hdul[0].header:
                    bzero=hdul[0].header['BZERO']
                if 'BSCALE' in hdul[0].header:
                    bscale=hdul[0].header['BSCALE']
                if str(hdul[0].data.dtype) == 'uint16' and hdul[0].header['BITPIX'] == 16:
                    hdul[0].data = ((hdul[0].data).astype('float64'))
                    hdul[0].data = hdul[0].data*bscale + bzero
                dark = CCDData(hdul[0].data, unit='adu')
                dark.header = hdul[0].header
                hdul.close()

                if scalable_dark_bool is True:
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
            try:
                master_dark.meta['combined'] = True
                dark_file_name = str(binning) + 'x' + str(binning) + '_master_dark_{:6.3f}.fits'.format(exp_time)
                master_dark.write(os.path.join(master_dir, dark_file_name),overwrite=True)
            except:
                raise RuntimeError('WARNING -- Could not Save Master Dark File')


def create_master_flat(all_fits: ccdp.ImageFileCollection,
                       master_dir: str,
                       scalable_dark_bool: str):
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
    #TODO: Return combined flat as an object

    Nil

    Outputs
    -----
    combined_flat saved in master_dir
    """

    unique_imagetype_list = list(set(all_fits.summary['imagetyp']))

    dark_imgtype_matches = [
        s for s in unique_imagetype_list if "dark" in s.lower()]
    dark_imgtypes_concatenateded = '|'.join(dark_imgtype_matches)
    flat_imgtype_matches = [
        s for s in unique_imagetype_list if "flat" in s.lower()]
    flat_imgtypes_concatenateded = '|'.join(flat_imgtype_matches)

    master_files = ccdp.ImageFileCollection(master_dir)

    flat_fits = ImageFileCollection(
        filenames=(all_fits.files_filtered(include_path=True, imagetyp=flat_imgtypes_concatenateded)))

    try:
        unique_bin_list = list(set(flat_fits.summary['ybinning']))
    except KeyError as e:
        raise e
        print("Could not find binning keyword")

    for binning in unique_bin_list:
        if scalable_dark_bool is True:
            try:
                # TODO: Standardize the Master Frame Gathering Process
                # combined_bias = CCDData.read(master_files_list)
                master_bias_string = str(binning) + 'x' + str(binning) + '_master_bias.fits'
                combined_bias = CCDData.read(os.path.join(master_dir, master_bias_string), unit='adu')
            except:
                raise RuntimeError('WARNING -- Could not open Master Bias file')

        try:
            # Get the Master Dark(s)
            master_darks = {ccd.header['exptime']: ccd for ccd in master_files.ccds(
                imagetyp=dark_imgtypes_concatenateded, ybinning=binning, combined=True
            )}
            # Create Conditonal Array with different Dark Image Data Types
            # TODO: Use np.where for this expression
            dark_mask = list(np.zeros(len(all_fits.summary), dtype=bool))
            for dark_imgtypes in dark_imgtype_matches:
                dark_mask = dark_mask + \
                            (all_fits.summary['imagetyp'] == (dark_imgtypes))

            dark_times = set(all_fits.summary['exptime'][dark_mask])

        except:
            print('WARNING -- Could not open Master Dark files')
            continue

        if len(master_darks) == 0:
            print('WARNING -- Could not open Master Dark files')
            continue

        flat_mask = list(np.zeros(len(all_fits.summary), dtype=bool))
        for flat_imgtypes in flat_imgtype_matches:
            flat_mask = flat_mask + \
                        (all_fits.summary['imagetyp'] == (flat_imgtypes))

        flat_filters = set(all_fits.summary['filter'][flat_mask])

        for frame_filter in sorted(flat_filters):
            flats_to_combine = all_fits.files_filtered(
                imagetyp=flat_imgtypes_concatenateded,
                filter=frame_filter,
                ybinning=binning,
                include_path=True
            )

            flats = []
            for file in flats_to_combine:

                # flat = CCDData.read(file, unit='adu')

                # Work around for applying BSCALE and BZERO
                hdul = fits.open(file)
                bzero = 32768
                bscale = 1
                if 'BZERO' in hdul[0].header:
                    bzero=hdul[0].header['BZERO']
                if 'BSCALE' in hdul[0].header:
                    bscale=hdul[0].header['BSCALE']

                if str(hdul[0].data.dtype)=='uint16' and hdul[0].header['BITPIX']==16:
                    hdul[0].data = ((hdul[0].data).astype('float64'))
                    hdul[0].data = hdul[0].data*bscale + bzero
                flat = CCDData(hdul[0].data, unit='adu')
                flat.header = hdul[0].header
                hdul.close()

                closest_dark = find_nearest_dark_exposure(flat, dark_times,
                                                          tolerance=100
                                                          )
                if scalable_dark_bool:
                    flat = ccdp.subtract_bias(flat, combined_bias)
                    flat = ccdp.subtract_dark(flat, master_darks[closest_dark],
                                              exposure_time='exptime',
                                              exposure_unit=u.second,
                                              scale=True
                                              )
                else:
                    flat = ccdp.subtract_dark(flat, master_darks[closest_dark],
                                              exposure_time='exptime',
                                              exposure_unit=u.second,
                                              scale=False
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
            try:
                combined_flat.meta['combined'] = True
                flat_file_name = '\\' + str(binning) + 'x' + str(binning) + '_master_flat_filter_{}.fits'.format(
                    frame_filter.replace("''", "p"))
                combined_flat.write((str(master_dir) + flat_file_name),overwrite=True)

                print('Saving ', flat_file_name)
                #TODO: Return combined_flat as an object
            except:
                raise RuntimeError('WARNING -- Could not Save Master Flat File')


def correct_lights(all_fits: ccdp.ImageFileCollection,
                   master_dir: str,
                   corrected_light_dir: str,
                   correct_outliers_params: str,
                   use_existing_masters: bool,
                   scalable_dark_bool: bool,
                   del_uncertainty_data=True,
                   convert_back_to_int=True):
    """
    Taking Master Bias, Master Darks and filter-specific Master Flats, create
    corrected individual light frames.


    Parameters
    ----------

    all_fits : ccdproc.ImageFileCollection object
        This object contains all .fits files in the directory and sub-directories
        in 'topdir' given by the user in the initialization section.


    master_dir : str
        The output directory to get the Master frames from.

    corrected_light_dir: str
       All corrected lights will be saved into this directory

    use_existing_masters: bool
        Boolean defining whether or not to use existing masters -- currently not in use

    scalable_dark_bool: bool
        Scale the dark frame with the exposure time.

     del_uncertainty_data: bool
        Delete the uncertainty data tied to the image reduction

    convert_back_to_int: bool
        Default is True. Converts the data back to wither signed 32 bit or 16 bit interger to save storage space


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

    light_imgs = ImageFileCollection(
        filenames=(all_fits.files_filtered(include_path=True, imagetyp=light_imgtypes_concatenateded)))
    try:

        unique_bin_list_y = list(set(light_imgs.summary['ybinning']))

    except KeyError as e:
        raise e
        print("Could not find binning keyword")
    # Create dictionaries of the dark and flat master frames in the master directory

    for binning in unique_bin_list_y:
        try:
            master_darks = {ccd.header['exptime']: ccd for ccd in master_files.ccds(
                imagetyp=dark_imgtypes_concatenateded, ybinning=binning, combined=True)}
            master_flats = {ccd.header['filter']: ccd for ccd in master_files.ccds(
                imagetyp=flat_imgtypes_concatenateded, ybinning=binning, combined=True)}

            # There is only one bias frame, so no need to set up a dictionary.
            if scalable_dark_bool:
                master_bias = [ccd for ccd in master_files.ccds(
                    imagetyp=bias_imgtypes_concatenateded, ybinning=binning, combined=True)][0]
            else:
                master_bias = None
            light_mask = list(np.zeros(len(all_fits.summary), dtype=bool))

            for light_imgtypes in light_imgtype_matches:
                light_mask = light_mask + \
                             (all_fits.summary['imagetyp'] == (light_imgtypes))

            light_filter = set(all_fits.summary['filter'][light_mask])

            dark_mask = list(np.zeros(len(all_fits.summary), dtype=bool))
            for dark_imgtypes in dark_imgtype_matches:
                dark_mask = dark_mask + \
                            (all_fits.summary['imagetyp'] == (dark_imgtypes))

            dark_times = set(all_fits.summary['exptime'][dark_mask])



            example_light = all_fits.files_filtered(imagetyp=light_imgtypes_concatenateded,
                                                    ybinning=binning,
                                                    filter=list(light_filter)[0],
                                                    include_path=True)[0]


            # Set the default data type for the light images, this will be used to reconvert the image after calibration
            bzero=32768 # Default value for BZERO
            bscale=1 # Default value for BScale



            corrected_master_darks = {}
            for dark_time in dark_times:

                ### Dark Image Masking ###
                '''
                
                '''

                mask = np.zeros(CCDData.read(example_light, unit='adu').data.shape, dtype=bool)

                if correct_outliers_params['Hot Pixel'] and correct_outliers_params['Outlier Boolean']:
                    # Find Hot Pixels using the master_dark
                    mask = find_hot_pixels(master_darks[dark_time], mask)

                if correct_outliers_params['Dark Frame Threshold Bool'] and correct_outliers_params['Outlier Boolean']:
                    # Apply A dark frame threshold to only allow some pixels to be passed through
                    mask = dark_frame_threshold(master_darks[dark_time], correct_outliers_params, mask)

                if correct_outliers_params['Outlier Boolean'] and (
                        correct_outliers_params['Hot Pixel'] or correct_outliers_params['Dark Frame Threshold Bool']):

                    # TODO : Add overwrite option
                    # Check to see if corrected file already exists
                    corrected_dark_dir = (str(master_dir) + '\\master_dark_' + str(dark_time) + 'corrected.fits')

                    ### Master Dark Image Correction ###

                    if os.path.isfile(corrected_dark_dir) is False:
                        #Correct the outliers in the master darks
                        corrected_master_dark_data = correct_outlier_darks(correct_outliers_params,
                                                                           master_darks[dark_time],dark_time,
                                                                           master_dir, mask)
                    else:
                        corrected_master_dark_data = CCDData.read(corrected_dark_dir, unit='adu')
                        print('WARNING -- Corrected Master Dark Already Exists, Using Existing Corrected Master Dark')
                    corrected_master_darks[dark_time] = corrected_master_dark_data
                ### Done Dark Image Masking and Correction ###

            # Iterate over the filter in the image
            for frame_filter in sorted(light_filter):
                lights_to_correct = all_fits.files_filtered(
                    imagetyp=light_imgtypes_concatenateded, ybinning=binning,
                    filter=frame_filter,
                    include_path=True)
                light_ccds = []
                try:
                    corrected_master_flat = []
                    if correct_outliers_params['Outlier Boolean'] is True:

                        ### Flat Image Masking ###
                        if correct_outliers_params['ccdmask']:
                            example_light = lights_to_correct[0]

                            if (correct_outliers_params['Hot Pixel'] or correct_outliers_params[
                                'Dark Frame Threshold Bool']):
                                # Hot Pixel or Dark Frame Threshold Bool is true then we must use the
                                # corrected_master_darks instead of the regular non-corrected darks
                                maskr, corrected_master_flat = flat_image_masking(flat_imgtypes_concatenateded,
                                                                                  lights_to_correct,
                                                                                  dark_times,
                                                                                  example_light,
                                                                                  master_bias,
                                                                                  corrected_master_darks,
                                                                                  scalable_dark_bool,
                                                                                  correct_outliers_params,
                                                                                  frame_filter,
                                                                                  master_dir,
                                                                                  corrected_light_dir,
                                                                                  all_fits, binning)
                            else:
                                # If Dark Masters have not been Outlier Corrected, use regular master_darks
                                maskr, corrected_master_flat = flat_image_masking(flat_imgtypes_concatenateded,
                                                                                  lights_to_correct,
                                                                                  dark_times,
                                                                                  example_light,
                                                                                  master_bias,
                                                                                  master_darks,
                                                                                  scalable_dark_bool,
                                                                                  correct_outliers_params,
                                                                                  frame_filter,
                                                                                  master_dir,
                                                                                  corrected_light_dir,
                                                                                  all_fits, binning)

                            # Combine masks
                            mask = mask | maskr

                    for file_name in lights_to_correct:
                        # Work around for applying BZERO and BSCALE
                        hdul = fits.open(file_name)
                        if 'BZERO' in hdul[0].header:
                            bzero = hdul[0].header['BZERO']
                        if 'BSCALE' in hdul[0].header:
                            bscale = hdul[0].header['BSCALE']
                        #TODO: Find a better way to work BZERO and BSCALE
                        if str(hdul[0].data.dtype) == 'uint16' and hdul[0].header['BITPIX'] == 16:
                            hdul[0].data = ((hdul[0].data).astype('float64'))
                            hdul[0].data = hdul[0].data*bscale + bzero
                        # hdul[0].data=((hdul[0].data).astype('float64'))
                        # hdul[0].data=hdul[0].data+hdul[0].header['BZERO']

                        light = CCDData(hdul[0].data, unit='adu')
                        light.header = hdul[0].header



                        # Note that the first argument in the remainder of the ccdproc calls is
                        # the *reduced* image, so that the calibration steps are cumulative.

                        ########## For Scalable darks (i.e. Bias Frame provided) ##########
                        if scalable_dark_bool == True:
                            reduced = ccdp.subtract_bias(light, master_bias)

                            if correct_outliers_params['Outlier Boolean'] and (
                                    correct_outliers_params['Dark Frame Threshold Bool'] or correct_outliers_params[
                                'Hot Pixel']):

                                closest_dark = find_nearest_dark_exposure(reduced, corrected_master_darks.keys())
                                reduced = ccdp.subtract_dark(reduced, corrected_master_darks[closest_dark],
                                                             exposure_time='exptime', exposure_unit=u.second,
                                                             scale=True)
                            else:

                                closest_dark = find_nearest_dark_exposure(
                                    reduced, master_darks.keys())

                                reduced = ccdp.subtract_dark(reduced,master_darks[closest_dark],
                                                             exposure_time='exptime', exposure_unit=u.second,
                                                             scale=True)

                        ########## For non-scalable darks (i.e. no Bias Fram provided) ##########
                        #### You also have to set run_master_bias to 0 on line 378 of Main.py ####
                        if scalable_dark_bool == False:

                            if correct_outliers_params['Outlier Boolean'] and (
                                    correct_outliers_params['Dark Frame Threshold Bool'] or correct_outliers_params[
                                'Hot Pixel']):
                                # if master_dark data should be outlier corrected
                                closest_dark = find_nearest_dark_exposure(light, corrected_master_darks.keys())
                                reduced = ccdp.subtract_dark(light, corrected_master_darks[closest_dark],
                                                             exposure_time='exptime', exposure_unit=u.second,
                                                             scale=False)
                            else:

                                closest_dark = find_nearest_dark_exposure(
                                    light, master_darks.keys())

                                reduced = ccdp.subtract_dark(light, master_darks[closest_dark],
                                                             exposure_time='exptime', exposure_unit=u.second,
                                                             scale=False)

                        ########## End Non ##########

                        ##### Flat Correction #####

                        if correct_outliers_params['Outlier Boolean'] and correct_outliers_params['ccdmask']:
                            #If Flat should be outlier corrected
                            try:
                                reduced = ccdp.flat_correct(reduced, corrected_master_flat)
                            except:
                                print('Corrected Master Flat Not Found, reverting to non-outlier_corrected flat')
                                good_flat = master_flats[reduced.header['filter']]
                                reduced = ccdp.flat_correct(reduced, good_flat)
                        else:
                            good_flat = master_flats[reduced.header['filter']]
                            reduced = ccdp.flat_correct(reduced, good_flat)

                        #### Cosmic Rays Outliers Correction ####
                        if correct_outliers_params['Outlier Boolean']:
                            if correct_outliers_params['Cosmic Rays Bool']:
                                if 'EGAIN' in light.header or 'GAINADU' in light.header:
                                    if 'EGAIN' in light.header:
                                        gain=light.header['EGAIN']
                                    elif 'GAINADU' in light.header:
                                        gain=light.header['GAINADU']
                                    try:
                                        reduced_in_e = ccdp.gain_correct(reduced, float(gain)*u.electron/u.adu)
                                        reduced_in_e.mask = mask
                                        new_reduced_in_e = ccdp.cosmicray_lacosmic(reduced_in_e, readnoise=0, sigclip=4, verbose=True)
                                        reduced = ccdp.gain_correct(new_reduced_in_e, (u.adu/(float(gain)*u.electron)))
                                        mask = reduced_in_e.mask
                                        print('Removed Cosmic Rays')
                                    except KeyError:
                                        print('Could Not Remove Cosmic Rays')
                            reduced.mask = mask

                        if correct_outliers_params['Force Offset']:
                            # Force the data to be shifted by the minimum value
                            if np.min(reduced.data) < 0:
                                reduced.data = reduced.data + abs(np.min(reduced.data))
                            else:
                                reduced.data = reduced.data - abs(np.min(reduced.data))


                        ### Final Steps in Processing ###
                        reduced.meta['correctd'] = True
                        reduced_hdul = reduced.to_hdu()
                        # Set all values below 0 equal to zero
                        reduced_hdul[0].data = np.where(reduced_hdul[0].data < 0, 0, reduced_hdul[0].data)

                        if convert_back_to_int:
                            if np.max(reduced_hdul[0].data) > 32767 and np.max(reduced_hdul[0].data) < 2147483647:
                                reduced_hdul[0].scale('int32')
                            elif np.max(reduced_hdul[0].data) < 32767:
                                reduced_hdul[0].scale('int16')




                        # If mask layer is empty, delete it to free up space
                        if np.sum(reduced_hdul[1].data) == 0:
                            del reduced_hdul[1]
                        if del_uncertainty_data: # Save Data by deleting the uncertainty layer
                            del reduced_hdul[1]
                        file_name = file_name.split("\\")[-1]
                        try:
                            reduced_hdul.writeto(str(corrected_light_dir) + '\\' + file_name,overwrite=True)
                        except OSError:
                            #File Name Already Taken
                            file_name = file_name[:-5]
                            file_name = file_name + "1.fits"
                            reduced_hdul.writeto(str(corrected_light_dir) + '\\' + file_name,overwrite=True)
                        print('Saved ', file_name)
                except Exception as e:
                    raise(e)
                    pass
        except Exception as e:
            raise(e)
            continue


def find_nearest_dark_exposure(image: astropy.nddata.CCDData,
                               dark_exposure_times: list,
                               tolerance: float =0.5) -> float :
    """
    Find the nearest exposure time of a dark frame to the exposure time of the image,
    raising an error if the difference in exposure time is more than tolerance.

    Parameters
    ----------

    image : astropy.nddata.CCDData
        Image for which a matching dark is needed.

    dark_exposure_times : list
        Exposure times for which there are darks.

    tolerance : numpy.float64 or ``None``, optional
        Maximum difference, in seconds, between the image and the closest dark. Set
        to ``None`` to skip the tolerance test.

    Returns
    -------

    float: numpy.float64
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


# %% Outlier Correction
def flat_image_masking(flat_imgtypes_concatenateded: str,
                       lights_to_correct: list,
                       dark_times: list,
                       example_light:  str,
                       master_bias: astropy.nddata.CCDData,
                       master_darks: ccdp.ImageFileCollection,
                       scalable_dark_bool: bool,
                       correct_outliers_params: dict,
                       frame_filter: str,
                       master_dir: str,
                       corrected_light_dir: str,
                       all_fits , binning: str):
    '''
    
    Produces a Mask on top of the flat frame image (or images) to identify hot pixels
    
    Parameters
    ----------
    flat_imgtypes_concatenateded : str
        A concatenated string describing the matching keywords for Flat fields images.
        Note: Dark flats will occassionally be passed through this. This will result in faulty calibration.
    lights_to_correct : list
        List of directory names that define the lights to be corrected.
    dark_times : list
        A list of the dark times in the master darks files
    example_light : str
        a string defining a sample frame
    master_bias : CCDATA object
        CCDATA of the master bias frame that has previously been generated
    master_darks : CCDATA object(s)
        CCDATA describing the master dark frame at various exposure times
    scalable_dark_bool : bool
        when true the dark frames are scaled to match the exposure time of 
        the ligth (science) frames
    correct_outliers_params : dict
        correct_outlier_params is a dictionary object which tells the program whether
        or not to calcualte 
    frame_filter : str
        string describing the filter for the image
    master_dir : str
        string defining the path for the directory contianing master frames
    corrected_light_dir : str
        string defining the path for the corrected light frames
    all_fits : ccdproc.ImageFileCollection object
        This object contains all .fits files in the directory and sub-directories
        in 'topdir' given by the user in the initialization section.
        
    binning : string
        defines the binnning of the images. Binning must be the same in both dimensions
        i.e. 3x3 

    Returns
    -------
    maskr : numpy.ma.masked_array
        Array which identifies hot/bad pixels. bad pixel=1, good pixel=0

    corrected_master_flat : astro.nddata.CCDData
        A new Master flat that has had it's outlier corrected for

    '''


    # Based on IRAF ccdmask
    flats_to_compare = all_fits.files_filtered(
        imagetyp=flat_imgtypes_concatenateded,
        filter=frame_filter,
        ybinning=binning,
        include_path=True
    )

    bad_ratio = False

    if len(flats_to_compare) > 1:
        # Creates Dictionary Object where key is the mean of the calibrated frame data
        calibrated_flats = {}
        # Calibrate Flats
        for flat in flats_to_compare:

            # flatdata = CCDData.read(flat, unit='adu')

            # Work around to apply  BZERO and BSCALE
            hdul = fits.open(flat)
            bzero = 32768
            bscale = 1
            if 'BZERO' in hdul[0].header:
                bzero = hdul[0].header['BZERO']
            if 'BSCALE' in hdul[0].header:
                bscale = hdul[0].header['BSCALE']

            if str(hdul[0].data.dtype) == 'uint16' and hdul[0].header['BITPIX'] == 16:
                hdul[0].data = ((hdul[0].data).astype('float64'))
                hdul[0].data = hdul[0].data*bscale + bzero

            flatdata = CCDData(hdul[0].data, unit='adu')
            flatdata.header = hdul[0].header
            hdul.close()

            #FIXME: Assumes all example_lights are the same exposure time
            closest_dark = find_nearest_dark_exposure(CCDData.read(example_light, unit='adu'), dark_times,
                                                      tolerance=100
                                                      )
            if scalable_dark_bool:
                ########## For Scalable darks (i.e. Bias Frame provided) ##########
                flatdata = ccdp.subtract_bias(flatdata, master_bias)

                flatdata = ccdp.subtract_dark(flatdata, master_darks[closest_dark],
                                              exposure_time='exptime',
                                              exposure_unit=u.second,
                                              scale=True
                                              )
            else:
                ########## For Non Scalable darks (i.e. Bias Frame not provided) ##########
                flatdata = ccdp.subtract_dark(flatdata, master_darks[closest_dark],
                                              exposure_time='exptime',
                                              exposure_unit=u.second,
                                              scale=False
                                              )

            calibrated_flats[flatdata.data.mean()] = flatdata

        max_key = max(calibrated_flats.keys())
        min_key = min(calibrated_flats.keys())
        ratio = calibrated_flats[max_key].divide(calibrated_flats[min_key])
        if (np.nanmean(ratio) < 1.1) and (np.nanmean(ratio) > 0.9):  # ratio thresholds
            bad_ratio = True
        else:
            maskr = ccdp.ccdmask(ratio)
            # TODO: Run this scenario

            corrected_master_flat = correct_outlier_flats(correct_outliers_params,
                                                          maskr,
                                                          flats_to_compare,
                                                          frame_filter,
                                                          master_dir,
                                                          corrected_light_dir)

    if (len(flats_to_compare) == 1) or (bad_ratio is True):
        # Calibrate flat field
        # flatdata = CCDData.read(flats_to_compare[0], unit='adu')

        # Work around to deal with BSCALE and BZERO
        hdul = fits.open(flats_to_compare[0])
        if 'BZERO' in hdul[0].header:
            bzero = hdul[0].header['BZERO']
        if 'BSCALE' in hdul[0].header:
            bscale = hdul[0].header['BSCALE']


        if str(hdul[0].data.dtype) == 'uint16' and hdul[0].header['BITPIX'] == 16:
            hdul[0].data = ((hdul[0].data).astype('float64'))
            hdul[0].data = hdul[0].data*bscale + bzero

        flatdata = CCDData(hdul[0].data, unit='adu')
        flatdata.header = hdul[0].header
        hdul.close()

        closest_dark = find_nearest_dark_exposure(flatdata, dark_times,
                                                  tolerance=100
                                                  )

        ########## For Scalable darks (i.e. Bias Frame provided) ##########
        if scalable_dark_bool:
            flatdata = ccdp.subtract_bias(flatdata, master_bias)
            flatdata = ccdp.subtract_dark(flatdata, master_darks[closest_dark],
                                          exposure_time='exptime',
                                          exposure_unit=u.second,
                                          scale=True
                                          )
        else:
            flatdata = ccdp.subtract_dark(flatdata, master_darks[closest_dark],
                                          exposure_time='exptime',
                                          exposure_unit=u.second,
                                          scale=False
                                          )

        maskr = ccdp.ccdmask(flatdata)

        flat_file_name = '\\master_flat_filter_corrected_{}.fits'.format(
            frame_filter.replace("''", "p"))
        if os.path.isfile((str(master_dir) + flat_file_name)) is False:

            corrected_master_flat = correct_outlier_flats(correct_outliers_params,
                                                          maskr,
                                                          flats_to_compare,
                                                          frame_filter,
                                                          master_dir,
                                                          corrected_light_dir)
        else:
            corrected_master_flat = CCDData.read((str(master_dir) + flat_file_name), unit='adu')

    return maskr, corrected_master_flat


def correct_outlier_flats(correct_outliers_params, maskr: numpy.ma , flats_to_compare, frame_filter, master_dir,
                          corrected_light_dir: str):
    '''
    This function will correct each flat passed into it by using the method ccdp.ccdmask. Multiple flats can be used
    for calibration if they have not been combined yet or existing master flats can be used.

    Parameters
    ----------
    correct_outliers_params : dict
        Describes the Correct outlier Paramers defined in the GUI
    maskr : numpy masked array 
        A mask for the bad pixels found in the flat images
    flats_to_compare : list
        list of strings defining the flats to be processed 
    frame_filter : string
        
    master_dir : string
        Directory of Master Directory 
    corrected_light_dir : string
        String deigning directory where the corrected lights are to be saved

    Returns
    -------
    result : astropy.nddata.CCData
        The result of the correct_outlier_flats is a combined flat that has had their outlier pixels corrected for.
    '''

    if correct_outliers_params['Replace Bool']:
        # Replaces the values in mask with a desired value
        flats = []
        if correct_outliers_params['Replace Mode'] == 'Ave':
            print('Replacing Outliers with Average')

            coordinates = np.where(maskr.data == True)

            if (correct_outliers_params['Multiple Flat Combination'] is False):
                if (flats_to_compare != ()):
                    flats_to_compare = flats_to_compare[0]
                flat = flats_to_compare
                # Work around to deal with BZERO and BSCALE
                hdul = fits.open(flat)
                bzero = 32768
                bscale = 1
                if 'BZERO' in hdul[0].header:
                    bzero=hdul[0].header['BZERO']
                if 'BSCALE' in hdul[0].header:
                    bscale=hdul[0].header['BSCALE']

                if str(hdul[0].data.dtype)=='uint16' and hdul[0].header['BITPIX']==16:

                    hdul[0].data = ((hdul[0].data).astype('float64'))
                    hdul[0].data = hdul[0].data*bscale + bzero
                flatdata = CCDData(hdul[0].data, unit='adu')
                flatdata.header = hdul[0].header
                hdul.close()
                flatdata.mask = maskr
                Replaceable_mean = np.nanmean(flatdata)
                for i in range(0, np.shape(coordinates)[1]):
                    flatdata.data[coordinates[0][i]][coordinates[1][i]] = Replaceable_mean
                flats.append(flatdata)
            else:
                #TODO: Check to see if this works correctly
                for flat in flats_to_compare:
                    # Work around to apply BZERO and BSCALE
                    hdul = fits.open(flat)
                    if 'BZERO' in hdul[0].header:
                        bzero = hdul[0].header['BZERO']
                    if 'BSCALE' in hdul[0].header:
                        bscale = hdul[0].header['BSCALE']
                    if str(hdul[0].data.dtype) == 'uint16' and hdul[0].header['BITPIX'] == 16:
                        hdul[0].data = ((hdul[0].data).astype('float64'))
                        hdul[0].data = hdul[0].data*bscale + bzero
                    flatdata = CCDData(hdul[0].data, unit='adu')
                    flatdata.header = hdul[0].header
                    hdul.close()
                    Replaceable_mean = np.nanmean(flatdata)
                    for i in range(0, np.shape(coordinates)[1]):
                        flatdata.data[coordinates[0][i]][coordinates[1][i]] = Replaceable_mean
                    flats.append(flatdata)

        elif correct_outliers_params['Replace Mode'] == 'Interpolate':
            radius = int(correct_outliers_params['Radius of local Averaging'])
            coordinates = np.where(maskr == True)
            flats = []
            if (correct_outliers_params['Multiple Flat Combination'] is False):
                if (flats_to_compare != ()):
                    flats_to_compare = flats_to_compare[0]
                flat = flats_to_compare

                # Work around to apply BZERO or BSCALE
                hdul = fits.open(flat)
                bzero = 32768
                bscale = 1
                if 'BZERO' in hdul[0].header:
                    bzero=hdul[0].header['BZERO']
                if 'BSCALE' in hdul[0].header:
                    bscale=hdul[0].header['BSCALE']
                if str(hdul[0].data.dtype)=='uint16' and hdul[0].header['BITPIX']==16:
                    hdul[0].data = ((hdul[0].data).astype('float64'))
                    hdul[0].data = hdul[0].data*bscale + bzero
                flatdata = CCDData(hdul[0].data, unit='adu')
                flatdata.header = hdul[0].header
                hdul.close()
                for i in range(0, np.shape(coordinates)[1]):
                    # Create copy of flat to add mask to - To Avoid Masked Pixs
                    flatdata2 = flatdata.copy()
                    flatdata2.mask = maskr
                    flatdata.mask = maskr
                    if (((coordinates[0][i] - radius) < 0) or ((coordinates[1][i] - radius) < 0) or (
                            (coordinates[0][i] + radius) > (np.shape(maskr))[0]) or (
                            (coordinates[1][0] + radius) > (np.shape(maskr))[1])):
                        # Condition for left most coordinates that have a radius beyond the left of the image
                        Replaceable_mean = np.nanmean(flatdata2)
                    else:
                        try:

                            Replaceable_mean = np.nanmean(
                                flatdata2[(coordinates[0][i] - radius):(coordinates[0][i] + radius),
                                (coordinates[1][i] - radius):(coordinates[1][i] + radius)])
                        except:  # Except faulty pix is in the corner of the image on the right side of the image
                            Replaceable_mean = np.nanmean(flatdata)
                    flatdata.data[coordinates[0][i]][coordinates[1][i]] = Replaceable_mean
                flats.append(flatdata)

            else:

                # Correct Flats
                for flat in flats_to_compare:
                    # flatdata = CCDData.read(flat, unit='adu')

                    # Work around to apply BZERO and BSCALE

                    hdul = fits.open(flat)
                    if 'BZERO' in hdul[0].header:
                        bzero = hdul[0].header['BZERO']
                    if 'BSCALE' in hdul[0].header:
                        bscale = hdul[0].header['BSCALE']
                    if str(hdul[0].data.dtype) == 'uint16' and hdul[0].header['BITPIX'] == 16:
                        hdul[0].data = ((hdul[0].data).astype('float64'))
                        hdul[0].data = hdul[0].data*bscale + bzero
                    flatdata = CCDData(hdul[0].data, unit='adu')
                    flatdata.header = hdul[0].header
                    hdul.close()

                    for i in range(0, np.shape(coordinates)[1]):

                        # Create copy of flat to add mask to - To Avoid Masked Pixs
                        flatdata2 = flatdata.copy()
                        flatdata2.mask = maskr
                        flatdata.mask = maskr

                        if (((coordinates[0][i] - radius) < 0) or ((coordinates[1][i] - radius) < 0)):
                            # Condition for left most coordinates that have a radius beyond the left of the image
                            Replaceable_mean = np.nanmean(flatdata2)
                        else:
                            try:
                                Replaceable_mean = np.nanmean(
                                    flatdata2[(coordinates[0][i] - radius):(coordinates[0][i] + radius),
                                    (coordinates[1][i] - radius):(coordinates[1][i] + radius)])
                            except:  # Except faulty pix is in the corner of the image on the right side of the image
                                Replaceable_mean = np.nanmean(flatdata2)
                        flatdata.data[coordinates[0][i]][coordinates[1][i]] = Replaceable_mean
                    flats.append(flatdata)

        # Save the Data
        if len(flats) == 1:
            flat_file_name = '\\master_flat_filter_corrected_{}.fits'.format(
                frame_filter.replace("''", "p"))
            flatdata.write(str(master_dir) + flat_file_name,overwrite=True)
            result = flatdata
        else:
            # TODO: Test This out
            if correct_outliers_params['Save Corrected Flats'] == True:
                corrected_dir = master_dir.split('\\')[0] + '\\corrected_flats'
                if os.path.isdir(corrected_dir) is False:
                    os.mkdir(corrected_dir)
                flat_name_dir = os.path.join(corrected_dir,
                                             os.path.join(os.path.splitext(os.path.basename(flat))[0],
                                                          'corrected.fits'))
                flatdata.write(flat_name_dir,overwrite=True)
            # IF ratio of flats is bad or only one flat frame we can save the flat as the master
            print('Saved New Flats')

            # Combine Flats 

            if correct_outliers_params['Multiple Flat Combination'] == True:
                corrected_dir = master_dir.split('\\')[0] + '\\corrected_flats'
                combined_flat = ccdp.combine(flats,
                                             method='median',
                                             sigma_clip=True,
                                             sigma_clip_low_thresh=5,
                                             sigma_clip_high_thresh=5,
                                             sigma_clip_func=np.ma.median,
                                             signma_clip_dev_func=np.ma.std,
                                             mem_limit=1e9
                                             )
                flat_file_name = '\\master_flat_filter_corrected_{}.fits'.format(
                    frame_filter.replace("''", "p"))
                combined_flat.write((str(master_dir) + flat_file_name),overwrite=True)
                result = combined_flat
        return result


def correct_outlier_darks(correct_outliers_params: dict, master_dark: astropy.nddata.CCDData,closest_dark: float,
                          master_dir:str, mask: np.ma.masked_array):
    '''
    This functions corrects the master_dark supplied in the function using either 'Ave' or 'interpolate' methods.
    'Ave' uses the average image data  to replace all the pixels that have been flagged as bad in the mask.
    'Interpolate' uses the local average of the image data to replace all the pixels that have been flagged as bad in the mask


    Parameters
    ----------
    correct_outliers_params dict: dict
        A dictionary containing the parameters used for correcting outliers,

    master_dark: astropy.nddata.CCDData
        A CCDData object containing image data that is used to to idenitify the outlier pixels in the image
    master_dir
    mask: np.ma.masked_array
        A numpy mask describing the bad pixels in the image

    Returns
    -------

    '''
    if correct_outliers_params['Replace Bool']:
        if correct_outliers_params['Replace Mode'] == 'Ave':
            # Correct the outliers using the average values
            print('Calculate Average Background')
            coordinates = np.where(mask.data == True)
            darkdata = (master_dark)
            Replaceable_mean = np.nanmean(darkdata)
            for i in range(0, np.shape(coordinates)[1]):
                darkdata.data[coordinates[0][i]][coordinates[1][i]] = Replaceable_mean


        elif correct_outliers_params['Replace Mode'] == 'Interpolate':
            radius = int(correct_outliers_params['Radius of local Averaging'])
            coordinates = np.where(mask.data == True)
            darkdata = (master_dark)
            darkdata2 = darkdata.copy()  # copy since we don't want new values messing up new values
            darkdata2.mask = mask.data
            average_dark_value = np.nanmean(darkdata2)
            for i in range(0, np.shape(coordinates)[1]):
                if (((coordinates[0][i] - radius) < 0) or ((coordinates[1][i] - radius) < 0) or (
                        (coordinates[0][i] + radius) > np.shape(mask)[0]) or (
                        (coordinates[1][i] + radius) > np.shape(mask)[1])):
                    Replaceable_mean = average_dark_value
                else:
                    try:
                        Replaceable_mean = np.nanmean(
                            darkdata2[(coordinates[0][i] - radius):(coordinates[0][i] + radius),
                            (coordinates[1][i] - radius):(coordinates[1][i] + radius)])
                    except:  # Mean of empty slice
                        Replaceable_mean = average_dark_value

                if str(Replaceable_mean) == 'nan':
                    Replaceable_mean = average_dark_value

                darkdata.data[coordinates[0][i]][coordinates[1][i]] = Replaceable_mean

        dark_name_dir = (str(master_dir) + '\\master_dark_' + str(closest_dark) + 'corrected.fits')
        darkdata.mask = mask
        darkdata.write(dark_name_dir,overwrite=True)

    return darkdata


def find_hot_pixels(master_dark: astropy.nddata.CCDData , mask: np.ndarray):

    '''

    Identifies hot pixels in the master_dark file passed into the function. Identified hot pixels are passed into a
    mask array and corrected for later.

    Parameters
    ----------
    master_dark: astropy.nddata.ccddata.CCDData
        The CCDData corresponding to the closest matched master dark frame

    mask: numpy.array
        Represents the masked pixel already in the data

    Returns
    -------
    mask: numpy.ma.core.MaskedArray
        Represents the updated masked pixel array

    '''


    # Not the best solution. Implement better handling of fits header inconsistencies.

    if 'EGAIN' in master_dark.header:
        gain=master_dark.header['EGAIN']
    elif 'GAINADU' in master_dark.header:
        gain=master_dark.header['GAINADU']

    try:
        dark_current = master_dark.multiply(
            float(gain)*u.electron /
            u.adu).divide(float(master_dark.header['EXPTIME'])*u.second)
    except KeyError:
        dark_current = master_dark.multiply(
            float(gain)*u.electron /
            u.adu).divide(float(master_dark.header['EXPTIME'])*u.second)

    if (dark_current):

        # If a pixel is 4 times the mean value of the average of the dark current frame, consider it a hot pixel

        hot_pixels = abs(dark_current.data) > abs(4 * np.nanmean(dark_current))
        mask = mask | hot_pixels
        mask = np.ma.masked_array(mask)


    else:
        mask = mask
    return mask


def dark_frame_threshold(master_dark: astropy.nddata.CCDData,
                         correct_outliers_params:dict,
                         mask: np.ma.masked_array):
    '''
    Performs a threshold filter on the dark frames. Only necessary if the dark image data deviates from expected
    significantly

    Parameters
    ----------
    master_dark: CCDData
        The master_dark CCData.
    correct_outliers_params: dict
        dictionary list containing the correct for outlier parameters
        'Dark Frame Threshold Max': maximum value for the filter
        'Dark Frame Threshold Min': minimum value for the filter

    mask: numpy.ma.masked_array
        An array describing the previously identified bad pixels.
    Returns
    -------
    mask: numpy.ma.masked_array
        A union mask between the inputted mask and the mask created from the threshold filter
    '''

    try:
        dark_pix_over_max = master_dark.data > float(
            correct_outliers_params['Dark Frame Threshold Max'])
        dark_pix_under_min = master_dark.data < float(
            correct_outliers_params['Dark Frame Threshold Min'])
        dark_threshold_mask = dark_pix_over_max | dark_pix_under_min
        mask = mask | dark_threshold_mask

    except Exception as e:
        raise RuntimeError(e)
    return mask
