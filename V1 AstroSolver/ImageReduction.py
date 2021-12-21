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



#-------------------------------------------------------------------------------
#
#   Functions
#
#-------------------------------------------------------------------------------

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
    ########## For St. John's ##########
    # bias_fits = all_fits.files_filtered(include_path=True, imagetyp = 'bias')
    ########## For ORC GBO files ##########
    bias_fits = all_fits.files_filtered(include_path=True, imagetyp = 'Bias Frame')
    ########## End ##########

    biases = []
    for file in bias_fits:
        bias = CCDData.read(file, unit = 'adu')
        biases.append(bias)


    master_bias = ccdp.combine(biases,
                            method='average',
                            sigma_clip=True,
                            sigma_clip_low_thresh=5,
                            sigma_clip_high_thresh=5,
                            sigma_clip_func=np.ma.median,
                            sigma_clip_dev_func = np.ma.std,
                            mem_limit = 1e9
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
##    master_files = ccdp.ImageFileCollection(master_dir)



    try:
##        combined_bias = CCDData.read(master_files.files_filtered(imagetyp ='bias',
##                                                                combined=True))
        master_bias = CCDData.read(master_dir + '\master_bias.fits')
    except:
        raise RuntimeError('WARNING -- Could not open Master Bias file')
    
    ########## For St. John's ##########
    # dark_mask = all_fits.summary['imagetyp'] == 'DARK'
    ########## For ORC GBO files ##########
    dark_mask = all_fits.summary['imagetyp'] == 'Dark Frame'
    ########## End ##########
    
    dark_times = set(all_fits.summary['exptime'][dark_mask])

    for exp_time in sorted(dark_times):
        
        ########## For St. John's ##########
        # darks_to_calibrate = all_fits.files_filtered(imagetyp='dark',
        #                                             exptime=exp_time,
        #                                             include_path=True)
        ########## For ORC GBO files ##########
        darks_to_calibrate = all_fits.files_filtered(imagetyp='Dark Frame',
                                                    exptime=exp_time,
                                                    include_path=True)
        ########## End ##########
        
        darks = []
        for file in darks_to_calibrate:
            dark = CCDData.read(file, unit = 'adu')
            dark = ccdp.subtract_bias(dark, master_bias)
            darks.append(dark)

        master_dark = ccdp.combine(darks,
                                     method='average',
                                     sigma_clip=True,
                                     sigma_clip_low_thresh=5,
                                     sigma_clip_high_thresh=5,
                                     sigma_clip_func=np.ma.median,
                                     signma_clip_dev_func=np.ma.std,
                                     mem_limit=1e9
                                    )

        master_dark.meta['combined'] = True

        dark_file_name = '\master_dark_{:6.3f}.fits'.format(exp_time)
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

    master_files = ccdp.ImageFileCollection(master_dir)
    
    ########## For St. John's ##########
    # master_files_list = master_files.files_filtered(imagetyp ='bias',
    #                                                 combined=True,
    #                                                 include_path = True
    #                                                 )
    ########## For ORC GBO files ##########
    master_files_list = master_files.files_filtered(imagetyp ='Bias Frame',
                                                    combined=True,
                                                    include_path = True
                                                    )
    ########## End ##########


    try:
        #combined_bias = CCDData.read(master_files_list)
        combined_bias = CCDData.read(master_dir + '\master_bias.fits')
    except:
        raise RuntimeError('WARNING -- Could not open Master Bias file')


    try:
        # Get the Master Darks
        
        ########## For St. John's ##########
        # master_darks = {ccd.header['exptime']: ccd for ccd in master_files.ccds(
        #                                         imagetyp='dark', combined=True
        #                                         )}
        ########## For ORC GBO files ##########
        master_darks = {ccd.header['exptime']: ccd for ccd in master_files.ccds(
                                                imagetyp='Dark Frame', combined=True
                                                )}
        ########## End ##########
        
        # Get a list of dark frame exposures from all dark frames
        
        ########## For St. John's ##########
        # dark_mask = all_fits.summary['imagetyp'] == 'DARK'
        ########## For ORC GBO files ##########
        dark_mask = all_fits.summary['imagetyp'] == 'Dark Frame'
        ########## End ##########
        
        dark_times = set(all_fits.summary['exptime'][dark_mask])
    except:
        raise RuntimeError('WARNING -- Could not open Master Dark files')

    if len(master_darks) == 0:
        raise RuntimeError('WARNING -- Could not open Master Dark files')
    
    ########## For St. John's ##########
    # flat_mask = all_fits.summary['imagetyp'] == 'FLAT'
    ########## For ORC GBO files ##########
    flat_mask = all_fits.summary['imagetyp'] == 'Flat Field'
    ########## End ##########
    
    flat_filters = set(all_fits.summary['filter'][flat_mask])

    for frame_filter in sorted(flat_filters):
        
        ########## For St. John's ##########
        # flats_to_combine = all_fits.files_filtered(imagetyp = 'flat',
        #                                             filter = frame_filter,
        #                                             include_path = True
        #                                             )
        ########## For ORC GBO files ##########
        flats_to_combine = all_fits.files_filtered(imagetyp = 'Flat Field',
                                                    filter = frame_filter,
                                                    include_path = True
                                                    )
        ########## End ##########
        
        flats = []
        for file in flats_to_combine:
            flat = CCDData.read(file, unit = 'adu')
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
                                     method='average',
                                     sigma_clip=True,
                                     sigma_clip_low_thresh=5,
                                     sigma_clip_high_thresh=5,
                                     sigma_clip_func=np.ma.median,
                                     signma_clip_dev_func=np.ma.std,
                                     mem_limit=1e9
                                     )

        combined_flat.meta['combined'] = True
        flat_file_name = '\master_flat_filter_{}.fits'.format(frame_filter.replace("''", "p"))
        combined_flat.write(master_dir + flat_file_name)
        print('Saving ',flat_file_name)





def correct_lights(all_fits, master_dir, corrected_light_dir):
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

    master_files = ccdp.ImageFileCollection(master_dir)

    # Create dictionaries of the dark and flat master frames in the master directory
    
    ########## For St. John's ##########
    # master_darks = {ccd.header['exposure']: ccd for ccd in master_files.ccds(imagetyp='dark', combined=True)}
    # master_flats = {ccd.header['filter']: ccd for ccd in master_files.ccds(imagetyp='flat', combined=True)}
    ########## For ORC GBO files ##########
    master_darks = {ccd.header['EXPTIME']: ccd for ccd in master_files.ccds(imagetyp='Dark Frame', combined=True)}
    master_flats = {ccd.header['filter']: ccd for ccd in master_files.ccds(imagetyp='Flat Frame', combined=True)}
    ########## For 8" Celestron ##########
    # master_flats = {ccd.header['filter']: ccd for ccd in master_files.ccds(imagetyp='Flat Field', combined=True)}
    ########## End ##########

    # There is only one bias frame, so no need to set up a dictionary.
    
    ########## For St. John's ##########
    # master_bias = [ccd for ccd in master_files.ccds(imagetyp='bias', combined=True)][0]
    ########## For ORC GBO files ##########
    master_bias = [ccd for ccd in master_files.ccds(imagetyp='Bias Frame', combined=True)][0]
    ########## End ##########
    
    ########## For St. John's ##########
    # light_mask = all_fits.summary['imagetyp'] == 'LIGHT'
    ########## For ORC GBO files ##########
    light_mask = all_fits.summary['imagetyp'] == 'Light Frame'
    ########## End ##########
    
    light_filter = set(all_fits.summary['filter'][light_mask])


    for frame_filter in sorted(light_filter):
        
        ########## For St. John's ##########
        # lights_to_correct = all_fits.files_filtered(imagetyp = 'light',
        #                                             filter = frame_filter,
        #                                             include_path = True
        #                                             )
        ########## For ORC GBO files ##########
        lights_to_correct = all_fits.files_filtered(imagetyp = 'Light Frame',
                                                    filter = frame_filter,
                                                    include_path = True
                                                    )
        ########## End ##########
        

        light_ccds = []

        try:
            for file_name in lights_to_correct:
                light = CCDData.read(file_name, unit = 'adu')
        
                # Note that the first argument in the remainder of the ccdproc calls is
                # the *reduced* image, so that the calibration steps are cumulative.
                
                ########## For Scalable darks (i.e. Bias Frame provided) ##########
                reduced = ccdp.subtract_bias(light, master_bias)
                
                closest_dark = find_nearest_dark_exposure(reduced, master_darks.keys())
        
                reduced = ccdp.subtract_dark(reduced, master_darks[closest_dark],
                                              exposure_time='exptime', exposure_unit=u.second,
                                              scale=True)
                ########## For non-scalable darks (i.e. no Bias Fram provided) ##########
                # closest_dark = find_nearest_dark_exposure(light, master_darks.keys())
        
                # reduced = ccdp.subtract_dark(light, master_darks[closest_dark],
                #                               exposure_time='exptime', exposure_unit=u.second,
                #                               scale=False)
                
                ########## End ##########
        
                good_flat = master_flats[reduced.header['filter']]
                reduced = ccdp.flat_correct(reduced, good_flat)
        
                reduced.meta['correctd'] = True
                file_name = file_name.split("\\")[-1]
                try:
                    reduced.write(corrected_light_dir + '\\' + file_name)
                except OSError:
                    file_name = file_name[:-5]
                    print(file_name)
                    file_name = file_name + "1.fits"
                    reduced.write(corrected_light_dir + '\\' + file_name)
        
        
                print('Saving ',file_name)
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



# def main():
#     pass



# #-------------------------------------------------------------------------------
# #
# #   Main Program
# #
# #-------------------------------------------------------------------------------



# if __name__ == '__main__':
#     main()

# #  Start timer

# start_time = time.time()

# # The first step is to recursively search a directory and subdirs to find all .fits files

# for dirpath, dirnames, files in os.walk(topdir):
#     for name in files:
#         if name.lower().endswith(exten):
#             results.append('%s' % os.path.join(dirpath, name))
# print('Have list of all .fits files')


# # %%

# # Using ImageFileCollection, gather all fits files

# all_fits = ImageFileCollection(filenames = results)
# print('Files sorted into ImageFileCollection object')


# # %%
# if run_master_bias == True:
#     print('\n')
#     print('Calling run_master_bias')
#     create_master_bias(all_fits, master_frame_directory)



# if run_master_dark == True:
#     print('\n')
#     print('Calling run_master_dark')
#     create_master_dark(all_fits, master_frame_directory)



# if run_master_flat == True:
#     print('\n')
#     print('Calling run_master_flat')
#     create_master_flat(all_fits, master_frame_directory)



# if correct_light_frames == True:
#     print('\n')
#     print('Creating output directory:', topdir + '\corrected_lights')
#     print('Calling correct_light_frames')

#     # Make output directory
#     correct_light_dir = Path(topdir, 'corrected_lights')
#     correct_light_dir.mkdir(exist_ok = True)
#     correct_light_directory = topdir + '\corrected_lights'

#     #  If a specific directory is desired for the corrected light frames, uncomment the
#     #  following line and place the path as a string with double backslashes.

#     #correct_light_directory = 'C:\\apple\\is_also_a_fruit'

#     #  Call function
#     correct_lights(all_fits, master_frame_directory, correct_light_directory)


# stop_time = time.time()
# elapsed_time = stop_time - start_time

# print('Elapsed time:', elapsed_time, 'seconds.')