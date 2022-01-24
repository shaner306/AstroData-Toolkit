# import pywin32_system32
import math
import multiprocessing as mp
import os
import time
from pathlib import Path
# from .AstroFunctions import *
# import TRMtester.py as trm
import PySimpleGUI as sg
import numpy as np
import pandas as pd
import tqdm
import win32com
import win32com.client
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.wcs import WCS
from ccdproc import ImageFileCollection
from numpy import mean
from photutils.background import Background2D
from photutils.background import MeanBackground
from photutils.background import MedianBackground
from photutils.background import ModeEstimatorBackground
from photutils.background import SExtractorBackground
import AstroFunctions as astro
import ImageReduction as IR
import SBDarkSub
import utils

"""
---------------------
AstroSolver
---------------------
    Functions:
        1. GUI Creation 
        2. Pinpoint Solver 
        3. Ground-Based Transform Calculation
        4. Space-Based Transform Calculation
        5. Track Rate Mode Photometry and Lightcurve Generator
        6. Track Rate Mode Astrometry (Experimental) - NOT WORKING CURRENTLY
        7. Ground-Based Image Reduction and Master Files Creation 
        8. NEOSSAT Dark Subtraction
        
        
   Function 1: GUI
   Function 2: PINPOINT SOLVER
   Function 3. Ground-Based Transform Calculation
   Function 4. Space-Based Transform Calculation
   Function 5. Track Rate Mode Photometry
   Function 6. Image Reduction
   Function 7. NEOSSAT Dark Subtraction
            
            
1. Space Based - Airmass not a factor in determining transforms
2. Ground Based - Multiple Order Transforms with Airmass and Extinctions

        ------------------------    
        STAR STARE MODE (SSM)
        ------------------------
            1.Solve Using pinpoint and IRAF
                - iraf_Sources= detecting_stars(fitsdata, bkg, bkg_std)
            
            2.Convert Pixel Location to RA and DEC
                - skypositions= convert_pixel_to_ra_dec(iraf_Sources, wcs)
            
            3.Convert RA Dec to Altitude and Azimuth
                - altazpositions = convert_ra_dec_to_alt_az(skypositions, header)
            
            4.Calculate FWHM
                - fwhm, fwhm_stdev= calculate_fwhm(iraf_Sources)
            
            5.Produce Photometry data Fluxes and Centroids
                - photometry_result= perform_photometry(iraf_Sources, fwhm, fitsdata, bkg)
            
            6.Calculate Instrumental Magnitudes + Sigmas of Matched Stars
                - calculate_magnitudes(photometry_result, exposure_Time)
                - calculate_magnitudes_sigma(photometry_result, exposure_Time)          
                
            7. Calculate Transforms

                a. Space Based Sensor
                   
                    Construct Tables to Store Data
                        large_table_columns= update_large_table_columns(large_table_columns, iraf_Sources, header, exposure_Time, ground_based=False, name_key='Name')
                        large_stars_table = create_large_stars_table(large_table_columns, ground_based=False)
                        stars_table= group_each_star(large_stars_table, ground_based=False, keys='Name')
                
                    Calculate Space Based Transforms    
                        filter_fci, zprime_fci = space_based_transform(stars_table, plot_results=False,index='(B-V)', app_filter='V', instr_filter='clear', field=None)
        
                    Calculate Standard Magnitude
                
                
                b. Ground Based Sensor
                    avg_Airmass= get_avg_airmass(altazpositions)
                
            8. Output Table to Excel
            
            9. Produce Plots and Save to PNG
                - With error bars
    2. Ground Based - Multiple Order Transforms with Airmass and Extinctions
        avg_Airmass= get_avg_airmass(altazpositions)
"""


def Gui():
    sg.theme("Default1")

    layout = [[sg.T("AstroSolver Processor V0.1")], [sg.T("   ")],
              [sg.Text("Image Folder: "),
               sg.Input(key="-IN2-", change_submits=True),
               sg.FolderBrowse(key="-IN1-")],
              [sg.Text("Catalog Folder: "),
               sg.Input(key="-IN3-", change_submits=True),
               sg.FolderBrowse(key="-IN4-")],
              [sg.Text("Reference Stars: "),
               sg.Input(key="-IN5-", change_submits=True),
               sg.FileBrowse(key="-IN6-")],
              [sg.Button("Submit")]]

    # layout = [[sg.T(" ")],
    #           [sg.Text("Image Folder: "),
    #            sg.Input(key="-IN2-" ,change_submits=True),
    #            sg.FolderBrowse(key="-IN1-")],
    #           [sg.Text("Catalog Folder: "),
    #            sg.Input(key="-IN3-" ,change_submits=True),
    #            sg.FolderBrowse(key="-IN4-")],
    #           [sg.Text("Reference Stars: "),
    #            sg.Input(key="-IN5-" ,change_submits=True),
    #            sg.FileBrowse(key="-IN6-")],
    #          [sg.T("                   "), sg.Checkbox('Print On:', default=True, key="-IN7-")],
    #           [sg.T("         "), sg.Radio('Permission Granted', "RADIO1", default=False, key="-IN8-")],
    #           [sg.T("         "), sg.Radio('Permission not Granted', "RADIO1", default=True)],
    #           [sg.Radio('Track Rate Mode', "RADIO3", default=False, key="-IN9-"),
    #            sg.Radio('Star Stare Mode', "RADIO3", default=True)],
    #           [sg.Button("Submit")]]

    # Building Window
    window = sg.Window('My File Browser', layout, size=(600, 250))

    while True:
        event, values = window.read()
        # print(values["-IN2-"])
        if event == sg.WIN_CLOSED or event == "Exit":
            break
        elif event == "Submit":

            imagefolder = values["-IN2-"]
            catalogfolder = values["-IN3-"]
            refdoc = values["-IN5-"]
            # print(values["-IN2-"])
            window.close()
            return imagefolder, catalogfolder, refdoc

# Calculate Various Image Background Values


def BackgroundEstimationMulti(fitsdata, sigma_clip, bkgmethod, printval):
    # Specify Sigma Clipping Value and Calculate the Background using "SExtractor Algorithm"
    sigma_clip = SigmaClip(sigma=2.5)
    bkg = SExtractorBackground(sigma_clip)

    # Calculate Mean Background
    bkg = MeanBackground(sigma_clip)
    bkg_value2 = bkg.calc_background(fitsdata)

    # Calculate Median Background
    bkg = MedianBackground(sigma_clip)
    bkg_value3 = bkg.calc_background(fitsdata)

    # Calculate Mode Background
    bkg = ModeEstimatorBackground(sigma_clip)
    bkg_value4 = bkg.calc_background(fitsdata)

    # Calculate SExtractor Background
    bkg_estimator2 = SExtractorBackground()

    # Remove Background Value
    # bkg = Background2D(fitsdata, (2, 2), filter_size=(3,3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2) Closest Approximate to Matlab Result
    bkg = Background2D(fitsdata, (50, 50), filter_size=(
        3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2)
    bg_rem = fitsdata - bkg.background

    if printval == 1:
        print("Background Solutions")
        print("Sigma Clip: " + str(sigma_clip))
        print("Using: " + bkgmethod)
        print("---------------------------")
        print("SExtractor Background: " + str(mean(bkg.background)))
        print("SExtractor Background(Filtered): " + str(mean(
            bkg.background)) + "\n " + "     " + "Box Size: " + "50x50" + "\n " + "     " + "Filter Size: " + "3x3")
        print("Mean Background: " + str(bkg_value2))
        print("Median Background: " + str(bkg_value3))
        print("Mode Estimator Background: " + str(bkg_value4))
        print("Remaining Background (subtracted): " + str(bg_rem))
        print("Polyfit Background: Not Implemented Yet")

    else:
        return

# Read Reference Star File and Extract Data


def ref_star_folder_read(refstars_doc):
    refstars = pd.read_excel(refstars_doc)
    refstars.head()
    HIP = refstars["HIP"]
    erad = refstars["erad"]
    edec = refstars["edec"]
    vref = refstars["V"]
    bvindex = refstars["(B-V)"]
    vrindex = refstars["(V-R)"]
    refstarsfin = np.column_stack((HIP, erad, edec, vref))
    return HIP, erad, edec, vref, bvindex, vrindex, refstarsfin

# Use Pinpoint to Locate Reference Stars on Image


def ref_star_search(s, f, erad, edec, HIP, vref, bvindex, vrindex, refstarsfin):
    refx = []
    refy = []
    vref2 = []
    HIP2 = []
    vrindexdet = []
    bvindexdet = []

    for s in range(89):

        try:
            f.SkyToXy(erad[s], edec[s])
            refSTAR_X = f.ScratchX  # the x result of skytoxy
            refSTAR_Y = f.ScratchY
            if refSTAR_X > 0 and refSTAR_X < 1900:
                refx.append(refSTAR_X)
                refy.append(refSTAR_Y)
                vref2.append(vref[s])
                HIP2.append(HIP[s])
                vrindexdet.append(vrindex[s])
                bvindexdet.append(bvindex[s])
            # else:
            #   print("Star Found outside bounds of image")
            # s=s+1
        except:
            if s > 89:
                print("Stop")
            # else:
            #     print('Pinpoint can''t process coords')

    nmstars = f.MatchedStars.Count
    mstars = f.MatchedStars
    print("Matched Stars:" + str(nmstars))
    print("Reference Stars Located:")
    print("")
    for i in range(1, nmstars):
        # print(i)
        mstar = mstars.Item(i)
        X_min = mstar.X - 0.5 * mstar.Width
        X_max = mstar.X + 0.5 * mstar.Width
        Y_min = mstar.Y - 0.5 * mstar.Height
        Y_max = mstar.Y + 0.5 * mstar.Height
        # print(i)
        length = len(refx)

        exptime = f.ExposureInterval
        rawflux = mstar.RawFlux
        Zp = f.MagZeroPoint

        vmag = Zp - 2.5 * (math.log10(rawflux / exptime))
        # starx= mstar.X
        # stary=mstar.Y
        # print(vmag)

        for j in range(length):
            # print("ref" +str(refx[j]))
            # print(X_max)
            if (refx[j] > X_min) and (refx[j] < X_max):
                if (refy[j] > Y_min) and (refy[j] < Y_max):
                    if abs(vmag - vref2[j]) < 0.5:
                        print("HIP: " + str(HIP2[j]))
                        print("Located at: X: " + str(mstar.X) +
                              " Y: " + str(mstar.Y))
                        # print("matched X:" + str(X_max))
                        # print(str(vref2[j]))
                        # print(mstar.ColorMagnitude)
                        print(
                            "Reference Mag: " + str(vref2[j]) + " vs " + "Detected Mag: " + str(vmag))
                        print("")
                        Bvtransform = (vref2[j] - vmag) / bvindexdet[j]
                        print("B-V Transform: " + str(Bvtransform))
                        Vrtransform = (vref2[j] - vmag) / vrindexdet[j]
                        print("V-R Transform: " + str(Vrtransform))


# Initialize Pinpoint Software
def pinpoint_init():
    # f = win32com.client.Dispatch("Pinpoint.plate")
    return


# Get List of Files/ Images for Pinpoint Solving
def getFileList(inbox):
    filepathall = []

    list1 = os.listdir(inbox)  # List of Files
    listSize = len(list1)  # Number of Files
    c = list1

    for i in range(0, listSize):
        filepath2 = inbox + "\\" + c[i]
        filepathall.append(filepath2)
    return filepathall


# Get Header Data from Individual .fits Image file
def fits_header_import(filepath, sb, filter_key='FILTER'):

    # Open HDUList in Astropy
    imagehdularray = fits.open(filepath)
    header = imagehdularray[0].header
    date = imagehdularray[0].header['DATE-OBS']

    # Determine if Image is Space-Based or Ground-Based
    if sb == 1:
        exposuretime = imagehdularray[0].header['AEXPTIME']
        XPIXSZ = 0
        YPIXSZ = 0

    else:
        exposuretime = imagehdularray[0].header['EXPTIME']
        XPIXSZ = imagehdularray[0].header['XPIXSZ']
        YPIXSZ = imagehdularray[0].header['YPIXSZ']
        focal_Length = imagehdularray[0].header['FOCALLEN']

    # Import Image Size in pixels on X and Y axis
    imagesizeX = imagehdularray[0].header['NAXIS1']
    imagesizeY = imagehdularray[0].header['NAXIS2']

    # Import Photometry Data and CCD Filter
    fitsdata = imagehdularray[0].data
    filt = imagehdularray[0].header['FILTER']

    # Convert header data to World Coordinate System (WCS) format
    wcs = WCS(header)
    return imagehdularray, date, exposuretime, imagesizeX, imagesizeY, fitsdata, filt, header, XPIXSZ, YPIXSZ, wcs

# Convert Image Size from Pixel to Arc Seconds and Ratio ArcSec/pix


def calc_ArcsecPerPixel(header):
    focal_Length = header['FOCALLEN']
    xpix_size = header['XPIXSZ']
    ypix_size = header['XPIXSZ']
    xbin = header['XPIXSZ']
    ybin = header['XPIXSZ']
    print(str(focal_Length) + " " + str(xpix_size))
    x_arcsecperpixel = math.atan(xpix_size / (focal_Length)) * 206.265
    y_arcsecperpixel = math.atan(ypix_size / (focal_Length)) * 206.265

    return x_arcsecperpixel, y_arcsecperpixel

# Clip Edges on Image for Source detection


def edge_Protect(bg_rem, edge_prot, imagesizeX, imagesizeY, fitsdata):
    bg_rem[1:edge_protect, 1:edge_protect] = 0
    bg_rem[imagesizeX - edge_protect:imagesizeX, :] = 0
    bg_rem[:, 1:edge_protect] = 0
    bg_rem[:, imagesizeY - edge_protect:imagesizeY] = 0
    im_mean = mean(bg_rem)
    im_rms = np.std(fitsdata)
    return im_mean, bg_rem, im_rms


# inbox, catloc, refstars_doc = Gui()
# print(imagefolder, catalogfolder, refdoc)
# Image Location of .fits Format
inbox = 'D:\\Wawrow\\2. Observational Data\\2021-03-10 - Calibrated\\HIP 46066\\LIGHT\\B'
inbox1 = r'D:\NEOSSat-SA-111\test'
ref_stars_file = r'D:\Astro2\Reference Star Files\Reference_stars_Apr29.txt'

# refstars_doc = 'D:\\Reference_stars.xlsx'
# refstars_csv='D:\\Reference_stars.csv' #Reference Star List
catloc1 = "D:\\squid\\USNOA20-All"
catloc2 = 'D:\\squid\\UCAC4'

save_loc = os.path.join(inbox, 'Outputs')  # Output Folder for Files

# Switches for Master Frame creation: 0 is Do Not Run, 1 is Run
create_master_dir = 1
run_master_bias = 1
run_master_dark = 1
run_master_flat = 1
correct_light_frames = 1
OutputsaveLoc = 0  # 0 Default will save outputs in image folder
reduce_dir = 'D:\Image Reduction Test Images'

# Start Pinpoint Software in Python
f = pinpoint_init()  # Start Pinpoint

# Set Image Processing Variables
streak_array = []  # Streak Detection for Track Rate Mode
edge_protect = 10  # Img Edge Clipping
min_obj_pixels = 5  # Min Pixels to qualify as a Point Source
SNRLimit = 0  # Signal-To-Noise Ratio

# Set Ground-Based or Space-Based, Chose only one
ground_based = False
space_based = False  # TODO Add Ground Based Selection Input

# Whether to Solve in Pinpoint before conducting Photometry
pinpoint = True  # TODO Add Pinpoint Solve or Not to GUI
plot_results = True  # Plot results on a light curve
TRM = False  # Conduct Astrometry or not
save_plots = True  # Save PNG files in directory after Photometry processing
remove_large_airmass = False
image_reduce = False  # Reduce Images before Solving

# Function #1: Pinpoint Solving


def pinpoint_solve(inbox, catloc, max_mag, sigma, catexp, match_residual, max_solve_time, cat, sb, use_sextractor, all_sky_solve):
    f = pinpoint_init()  # FIXME: Doesn't do anything
    filepathall = getFileList(inbox)  # FIXME: Doesn't Do Anything
    file_suffices = [".fits", ".fit", ".fts"]

    for dirpath, dirnames, filenames in os.walk(inbox):
        for filename in filenames:
            if (filename.endswith(any(file_suffices))):
                filepath = os.path.join(dirpath, filename)
                f = win32com.client.Dispatch("Pinpoint.plate")
                print("Processing Image: " + filepath)

                "Import Data from FITS Image"
                header = 0
                imagehdularray, date, exposure_Time, imagesizeX, imagesizeY, fitsdata, filt, header, XPIXSZ, YPIXSZ, wcs = fits_header_import(
                    filepath, sb)

                """Pinpoint Solve"""
                # FIXME: Why is pinpont variable here when always TRUE?
                if pinpoint:
                    try:
                        f.AttachFITS(filepath)
                        f.Declination = f.targetDeclination
                        f.RightAscension = f.targetRightAscension

                        x_arcsecperpixel, y_arcsecperpixel = calc_ArcsecPerPixel(
                            header)
                        # yBin = 4.33562092816E-004*3600;
                        # xBin =  4.33131246330E-004*3600;
                        # f.ArcsecperPixelHoriz  = 4.556
                        # f.ArcsecperPixelVert =  4.556
                        # if sb==1:
                        #     yBin = 4.33562092816E-004*3600*2 #%Image Specific Pixel Size in arcsec / Obtained from FITS Header and converted from deg
                        #     xBin =  4.33131246330E-004*3600*2#%Image Specific Pixel Size in arcsec / Obtained from FITS Header and converted from deg
                        # else:
                        #     yBin = 4.33562092816E-004*3600 #%Image Specific Pixel Size in arcsec / Obtained from FITS Header and converted from deg
                        #     xBin =  4.33131246330E-004*3600 #%Image Specific Pixel Size in arcsec / Obtained from FITS Header and converted from deg
                        # f.ArcsecperPixelHoriz  = xBin    #%CCD Pixel scale on CD
                        # f.ArcsecperPixelVert = yBin
                        f.ArcsecperPixelHoriz = x_arcsecperpixel  # %CCD Pixel scale on CD
                        f.ArcsecperPixelVert = y_arcsecperpixel

                        "Pinpoint Solve Inputs"
                        # TODO Add Inputs for pinpoint solving to GUI
                        f.Catalog = 11
                        f.CatalogPath = catloc
                        f.UseSExtractor = use_sextractor
                        f.CatalogMaximumMagnitude = max_mag
                        # print(f.CatalogMaximumMagnitude)
                        f.CatalogExpansion = catexp
                        f.MaxSolveTime = max_solve_time
                        f.MaxMatchResidual = match_residual
                        f.SigmaAboveMean = sigma
                        f.RemoveHotPixels()

                        "Pinpoint Solving"
                        f.FindCatalogStars()
                        # print(f.CatalogStars.Count)
                        f.FindImageStars()
                        # print(f.ImageStars.Count)
                        f.Solve()
                        f.UpdateFITS()
                        # print( f.MatchedStars.count)
                        # f.FindImageStars()
                        # print(f.ImageStars)

                        f.DetachFITS()
                        f = None
                        print("Pinpoint Solved")
                    except:
                        print("Could Not Solve")
                        continue

# Calculate Ground-Based Transforms


def Ground_based_transforms(directory, ref_stars_file):
    plot_results = True
    save_plots = True
    remove_large_airmass = False
    # file_suffix=".fit"
    # exposure_key='EXPTIME'
    lat_key = 'OBSGEO-B'
    lon_key = 'OBSGEO-L'
    elev_key = 'OBSGEO-H'
    file_suffix = ".fits"
    exposure_key = 'EXPTIME'
    # lat_key='SITELAT'
    # lon_key='SITELONG'
    # elev_key='SITEELEV'
    name_key = 'Name'
    # ref_stars_file = r'D:\Astro2\Reference Star Files\Reference_stars_Apr29.txt'
    save_loc = os.path.join(directory, 'Outputs')
    if not os.path.exists(save_loc):
        os.makedirs(save_loc)
    gb_final_transforms, \
        auxiliary_data_table = \
        astro._main_gb_transform_calc(directory,
                                      ref_stars_file,
                                      plot_results=plot_results,
                                      save_plots=save_plots,
                                      remove_large_airmass_bool=remove_large_airmass,
                                      file_suffix=file_suffix,
                                      exposure_key=exposure_key,
                                      lat_key=lat_key,
                                      lon_key=lon_key,
                                      elev_key=elev_key,
                                      name_key=name_key,
                                      save_loc=save_loc)

    gb_final_transforms.pprint_all()
    auxiliary_data_table.pprint_all()
    return

# Calculate Space-Based Transforms


def space_based_transform(directory, ref_stars_file):
    plot_results = True
    save_plots = True
    file_suffix = "_clean.fits"
    exposure_key = 'AEXPTIME'
    name_key = 'Name'
    transform_index_list = ['(B-V)', '(V-R)', '(V-I)']

    # for subfolder in subfolder_list:
    unique_id = ''
    directory = r'D:\NEOSSat-SA-111\test'
    save_loc = os.path.join(directory, 'Outputs')
    if not os.path.exists(save_loc):
        os.makedirs(save_loc)

    sb_final_transform_table =\
        astro._main_sb_transform_calc(directory,
                                      ref_stars_file,
                                      plot_results=plot_results,
                                      save_plots=save_plots,
                                      file_suffix=file_suffix,
                                      exposure_key=exposure_key,
                                      name_key=name_key,
                                      transform_index_list=transform_index_list,
                                      save_loc=save_loc,
                                      unique_id=unique_id)
    sb_final_transform_table.pprint_all()

    return

# Conduct Track Rate Mode (TRM) Photometry on satellites


def trm_photometry(directory):

    # Set
    temp_dir = 'tmp'
    max_distance_from_sat = 20
    size = 20
    max_num_nan = 5
    plot_results = 0

    # Produce Data Tables and Conduct Photometry on Images
    sats_table, uncertainty_table, sat_fwhm_table = astro._main_sc_lightcurve(directory,
                                                                              temp_dir=temp_dir,
                                                                              max_distance_from_sat=max_distance_from_sat,
                                                                              size=size,
                                                                              max_num_nan=max_num_nan,
                                                                              plot_results=plot_results)

    # Print Data
    sat_fwhm_table.pprint_all()
    uncertainty_table.pprint_all()
    sats_table.pprint_all()

# Function #6: GB Image Reduction


def Image_reduce(reduce_dir, create_master_dark, create_master_flat, create_master_bias, create_master_dir=True):

    if os.path.isdir(reduce_dir) == False:
        raise RuntimeError(
            'WARNING -- Directory of .fits files does not exist')

    # Create output directory for master files
    if create_master_dir == True:
        master_frame_dir = Path(reduce_dir, 'master_frame_data')
        master_frame_dir.mkdir(exist_ok=True)

    #  Select directory for master frames
    master_frame_directory = reduce_dir + '\master_frame_data'

    #  If a directory already exists containing the master files, uncomment the
    #  following line and place the path as a string with double backslashes.

    # master_frame_directory = 'C:\\pineapple\\is_a_fruit'
    if os.path.isdir(master_frame_directory) == False:
        raise RuntimeError(
            'WARNING -- Directory of Master .fits files does not exist')

    # The extension to search for
    exten = '.fits'

    # Fits files found through directory search
    results = []

    # Track the Runtime
    start_time = time.time()

    # Find all fits files in subdirectories
    for dirpath, dirnames, files in os.walk(reduce_dir):
        for name in files:
            if name.lower().endswith(exten):
                results.append('%s' % os.path.join(dirpath, name))
    print('Have list of all .fits files')

    # Using ImageFileCollection, gather all fits files

    all_fits = ImageFileCollection(filenames=results)
    print('Files sorted into ImageFileCollection object')

    # %% Image Reduction
    # Create Master Bias
    if run_master_bias == True:
        print('\n')
        print('Calling run_master_bias')
        IR.create_master_bias(all_fits, master_frame_directory)

    # Create Master Dark
    if run_master_dark == True:
        print('\n')
        print('Calling run_master_dark')
        IR.create_master_dark(all_fits, master_frame_directory)

    # Create Master Flat
    if run_master_flat == True:
        print('\n')
        print('Calling run_master_flat')
        IR.create_master_flat(all_fits, master_frame_directory)

    # Correct Light Frames with Master Files
    if correct_light_frames == True:
        print('\n')
        print('Creating output directory:', reduce_dir + '\corrected_lights')
        print('Calling correct_light_frames')

        # Make output directory
        correct_light_dir = Path(reduce_dir, 'corrected_lights')
        correct_light_dir.mkdir(exist_ok=True)
        correct_light_directory = reduce_dir + '\corrected_lights'

        #  If a specific directory is desired for the corrected light frames, uncomment the
        #  following line and place the path as a string with double backslashes.

        # correct_light_directory = 'C:\\apple\\is_also_a_fruit'

        #  Call function
        IR.correct_lights(all_fits, master_frame_directory,
                          correct_light_directory)

    stop_time = time.time()
    elapsed_time = stop_time - start_time

    print('Elapsed time:', elapsed_time, 'seconds.')

    # %% NEOSSAT Dark Subtraction


def DarkSub(target, obspath, **kwargs):
    """Process all observations of a specific target in a specific directory."""

    # Make sure output directory exists.
    save_loc = os.path.join(obspath, 'Outputs')
    if not os.path.exists(save_loc):
        os.makedirs(save_loc)

    # Unpack parameters.
    T = kwargs.pop('T', 8)
    bpix = kwargs.pop('bpix', -1.0e10)
    snrcut = kwargs.pop('snrcut', 10.0)
    fmax = kwargs.pop('fmax', 2)
    xoff = kwargs.pop('xoff', 0)
    yoff = kwargs.pop('yoff', 0)
    nproc = kwargs.pop('nproc', 4)

    print('Processing observations in directory {}'.format(obspath))

    # Create a table of all fits files in the specified diretcory.
    obs_table = utils.observation_table(obspath)

    # Parse the observation to select the desired target and appropriate darks.
    light_table, dark_table = utils.parse_observation_table(obs_table, target)
    nlight, ndark = len(light_table), len(dark_table)

    # Process the dark images.
    print('Processing {} dark images to create the masterdark.'.format(ndark))

    # Use multiprocessing to process all dark images.
    results = []
    pbar = tqdm.tqdm(total=ndark)
    with mp.Pool(nproc) as p:

        for i in range(ndark):
            darkfile = dark_table['FILENAME'][i]
            xsc, ysc = dark_table['xsc'][i], dark_table['ysc'][i]
            xov, yov = dark_table['xov'][i], dark_table['yov'][i]

            args = (obspath, darkfile, xsc, ysc, xov, yov,
                    snrcut, fmax, xoff, yoff, T, bpix)

            results.append(p.apply_async(SBDarkSub.darkprocess,
                           args=args, callback=lambda x: pbar.update()))

        p.close()
        p.join()

        alldarkdata = [result.get() for result in results]

    pbar.close()

    # Combine the processed darks to obtain a master dark.
    masterdark = SBDarkSub.combinedarks(alldarkdata)

    # Display the master dark.
    #    imstat = utils.imagestat(masterdark, bpix)
    #    visualize.plot_image(masterdark, imstat, 0.3, 10.0)

    # Save the masterdark.
    head, tail = os.path.split(obspath)
    darkname = 'masterdark_{}_{}.fits'.format(tail, target)
    darkname = os.path.join(save_loc, darkname)
    hdu = fits.PrimaryHDU(masterdark)
    hdu.writeto(darkname, overwrite=True)

    # Clear memory.
    del alldarkdata

    # Process the light images.
    print('Processing {} light images.'.format(nlight))

    # Use multiproessing to process all light images.
    pbar = tqdm.tqdm(total=nlight)
    with mp.Pool(nproc) as p:

        for i in range(nlight):
            filename = os.path.join(obspath, light_table['FILENAME'][i])
            xsc, ysc = light_table['xsc'][i], light_table['ysc'][i]
            xov, yov = light_table['xov'][i], light_table['yov'][i]

            args = (filename, save_loc, masterdark, xsc, ysc,
                    xov, yov, snrcut, fmax, xoff, yoff, T, bpix)

            p.apply_async(SBDarkSub.lightprocess_save, args=args,
                          callback=lambda x: pbar.update())

        p.close()
        p.join()

    pbar.close()

    return


class AstroReducer:

    def __init__(self, targets, start_time, stop_time, *args):
        return
