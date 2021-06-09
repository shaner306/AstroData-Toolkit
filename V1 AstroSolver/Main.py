import pandas as pd
import win32com.client as win32
import win32com
import os
#import pywin32_system32
import math
import numpy as np
from numpy import mean
from astropy.stats import SigmaClip
from photutils.background import SExtractorBackground
from photutils.background import MeanBackground
from photutils.background import Background2D
from photutils.background import ModeEstimatorBackground
from photutils.background import MedianBackground
from astropy.wcs import WCS
import AstroFunctions as astro
import ImageReduction as IR
#from .AstroFunctions import *
#import TRMtester.py as trm
import numpy
import PySimpleGUI as sg
from pathlib import Path
from astropy.io import fits
from ccdproc import ImageFileCollection
import time
import utils
import Visualize
import SBDarkSub
import multiprocessing as mp

import tqdm
import numpy as np
from scipy import optimize
from scipy import fftpack
import sys

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval
from photutils import DAOStarFinder, CircularAperture, aperture_photometry
"""
Classes:
---------------------
AstroSolver
---------------------
    Stage 1
    - Open Image folder and determine the number of files
    - Open first file and store FITS data
    - Generate output Logs
    - Import ref star file
    Stage 2
    - Solve image with pinpoint and IRAFstarfinder (data from both)
    - Scan image for ref star
    - Report pinpoint and IRAF data to report log and table
    Stage 3
    - Calculate transforms if possible (GBO only?)
    - Background subtraction and estimation
    - Fit Moffat, gaussian, and lorentzian profiles for PSF
    - Output all data to spread sheet
    Stage 4
    - Produce report log with errors, random data, image quality report
    - Produce light curve charts, graphs, and images and save them to an output folder
    (ME ONLY) Comment and document entire usage and process
    (ME ONLY) PACKAGE V1 properly
    ADD additional photometry data
    
    STAR STARE MODE (SSM)
        Input Images Folder
        Images are Indexed
        First Image Selected
        Store Header and Data
        
        Estimate Background and Background Standard Deviation
            - bkg = BackgroundEstimationMulti(fitsdata, 2.5, 1, 0)
            - backg,bkg_std = calculate_img_bkg(fitsdata)
            
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
    
    TRACK RATE MODE (TRM)
        
-------------------    
AstroReducer    
-------------------
"""

def Gui ():

    sg.theme("Default1")
    
    layout = [[sg.T("AstroSolver Processor V0.1")], [sg.T("   ")],
              [sg.Text("Image Folder: "), 
               sg.Input(key="-IN2-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN1-")],
              [sg.Text("Catalog Folder: "), 
               sg.Input(key="-IN3-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN4-")],
              [sg.Text("Reference Stars: "), 
               sg.Input(key="-IN5-" ,change_submits=True), 
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
    
    ###Building Window
    window = sg.Window('My File Browser', layout, size=(600,250))
        
    while True:
        event, values = window.read()
        #print(values["-IN2-"])
        if event == sg.WIN_CLOSED or event=="Exit":
            break
        elif event == "Submit":
            
            imagefolder =values["-IN2-"]
            catalogfolder = values["-IN3-"]
            refdoc = values["-IN5-"]
            #print(values["-IN2-"])
            window.close()
            return imagefolder, catalogfolder, refdoc

def BackgroundEstimationMulti(fitsdata, sigma_clip, bkgmethod, printval):
   
    sigma_clip = SigmaClip(sigma=2.5)
    bkg = SExtractorBackground(sigma_clip)
    #bkg_value1 = bkg.calc_background(fitsdata)
    
    #print(bkg_value)
    bkg = MeanBackground(sigma_clip)
    bkg_value2 = bkg.calc_background(fitsdata)
    
    bkg = MedianBackground(sigma_clip)
    bkg_value3 = bkg.calc_background(fitsdata)
    
    bkg = ModeEstimatorBackground(sigma_clip)
    bkg_value4 = bkg.calc_background(fitsdata)
    
    

    bkg_estimator2 = SExtractorBackground()
    #bkg = Background2D(fitsdata, (2, 2), filter_size=(3,3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2) Closest Approximate to Matlab Result
    bkg = Background2D(fitsdata, (50,50), filter_size=(3,3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2)
    bg_rem = fitsdata - bkg.background

    if printval == 1:
        print("Background Solutions")
        print("Sigma Clip: " + str(sigma_clip))
        print("Using: "+ bkgmethod)
        print("---------------------------")
        print("SExtractor Background: " + str(mean(bkg.background)))
        print("SExtractor Background(Filtered): " + str(mean(bkg.background))+"\n "+ "     " + "Box Size: " +"50x50" +"\n "+ "     " + "Filter Size: " +"3x3")
        print("Mean Background: " + str(bkg_value2))
        print("Median Background: " + str(bkg_value3))
        print("Mode Estimator Background: " + str(bkg_value4))
        print("Remaining Background (subtracted): " + str(bg_rem))
        print("Polyfit Background: Not Implemented Yet")
        
    else:
        return

    
def ref_star_folder_read(refstars_doc):
    
    refstars = pd.read_excel(refstars_doc)
    refstars.head()
    HIP= refstars["HIP"]
    erad = refstars["erad"]
    edec= refstars["edec"]
    vref= refstars["V"]
    bvindex=refstars["(B-V)"]
    vrindex=refstars["(V-R)"]
    refstarsfin= np.column_stack((HIP, erad,edec,vref))
    return HIP, erad,edec,vref,bvindex,vrindex,refstarsfin



def ref_star_search(s,f,erad,edec, HIP, vref,bvindex,vrindex,refstarsfin):
    refx=[]
    refy=[]
    vref2=[]
    HIP2=[]
    vrindexdet=[]
    bvindexdet=[]
    
    
    for s in range(89):
        
        try:
            f.SkyToXy(erad[s],edec[s]);
            refSTAR_X  = f.ScratchX   #the x result of skytoxy
            refSTAR_Y = f.ScratchY
            if refSTAR_X>0 and refSTAR_X<1900:
                
                refx.append(refSTAR_X)
                refy.append(refSTAR_Y)
                vref2.append(vref[s])
                HIP2.append(HIP[s])
                vrindexdet.append(vrindex[s])
                bvindexdet.append(bvindex[s])
            # else:
            #   print("Star Found outside bounds of image")
            #s=s+1
        except:
            if s > 89:
                print("Stop")
            # else:
            #     print('Pinpoint can''t process coords')
                                                      
    nmstars = f.MatchedStars.Count
    mstars = f.MatchedStars;
    print("Matched Stars:"+ str(nmstars))
    print("Reference Stars Located:")
    print("")
    for i in range(1,nmstars):
            #print(i)
            mstar = mstars.Item(i)
            X_min = mstar.X - 0.5*mstar.Width
            X_max = mstar.X + 0.5*mstar.Width
            Y_min = mstar.Y - 0.5*mstar.Height
            Y_max = mstar.Y + 0.5*mstar.Height
            #print(i)
            length = len(refx)
            
            exptime = f.ExposureInterval  
            rawflux = mstar.RawFlux;            
            Zp = f.MagZeroPoint;
            
            vmag= Zp - 2.5*(math.log10(rawflux/exptime))
            #starx= mstar.X
            #stary=mstar.Y
            #print(vmag)

            for j in range(length):
                #print("ref" +str(refx[j]))
                #print(X_max)
                if (refx[j] > X_min) and (refx[j] < X_max):
                    if (refy[j] > Y_min) and (refy[j] < Y_max):
                          if abs(vmag - vref2[j]) < 0.5:
                                  
                                  print("HIP: " +str(HIP2[j]))
                                  print("Located at: X: " +str(mstar.X) + " Y: " + str(mstar.Y))
                                  #print("matched X:" + str(X_max))
                                  #print(str(vref2[j]))
                                  #print(mstar.ColorMagnitude)
                                  print("Reference Mag: "+str(vref2[j]) + " vs " + "Detected Mag: " + str(vmag))
                                  print("")
                                  Bvtransform=(vref2[j]-vmag)/bvindexdet[j]
                                  print("B-V Transform: " + str(Bvtransform))
                                  Vrtransform=(vref2[j]-vmag)/vrindexdet[j]
                                  print("V-R Transform: " + str(Vrtransform))
    
    


def pinpoint_init():
    #f = win32com.client.Dispatch("Pinpoint.plate")
    return

def getFileList(inbox):
    filepathall = []
    
    list1 = os.listdir(inbox) #List of Files
    listSize = len(list1) #Number of Files
    print(listSize)
    c=list1
    print(c)
    #o=0;
    for i in range(1,listSize):
        print(c[i])
        filepath2 = inbox+"\\"+c[i]
        filepathall.append(filepath2)
        #o=o+1;
    return filepathall
    #o=0;
def fits_header_import(filepath, filter_key='FILTER'):
        imagehdularray = fits.open(filepath)
        header = imagehdularray[0].header
        date=imagehdularray[0].header['DATE-OBS']
        exposuretime=imagehdularray[0].header['EXPTIME']
        imagesizeX=imagehdularray[0].header['NAXIS1']
        imagesizeY=imagehdularray[0].header['NAXIS2']
        fitsdata =  imagehdularray[0].data
        #focal_Length=imagehdularray[0].header['FOCALLEN']
        XPIXSZ=imagehdularray[0].header['XPIXSZ']
        YPIXSZ=imagehdularray[0].header['YPIXSZ']
        filt=imagehdularray[0].header['FILTER']
        wcs = WCS(header)
        return imagehdularray,date,exposuretime,imagesizeX,imagesizeY, fitsdata, filt,header, XPIXSZ, YPIXSZ, wcs

def calc_ArcsecPerPixel(header):
   focal_Length= header['FOCALLEN']
   xpix_size=header['XPIXSZ']
   ypix_size=header['XPIXSZ']
   xbin=header['XPIXSZ']
   ybin=header['XPIXSZ']
   x_arcsecperpixel = math.atan(xpix_size/focal_Length)*3600*xbin
   y_arcsecperpixel = math.atan(ypix_size/focal_Length)*3600*ybin
   
   return x_arcsecperpixel, y_arcsecperpixel

def edge_Protect (bg_rem, edge_prot, imagesizeX, imagesizeY, fitsdata):
                    
        bg_rem[1:edge_protect,1:edge_protect] = 0;
        bg_rem[imagesizeX - edge_protect:imagesizeX, :] = 0
        bg_rem[:, 1:edge_protect] = 0
        bg_rem[:, imagesizeY - edge_protect:imagesizeY] = 0
        im_mean = mean(bg_rem)
        im_rms=np.std(fitsdata)
        return im_mean, bg_rem, im_rms

"#TODO GUI SETUP"

# inbox, catloc, refstars_doc = Gui()
# print(imagefolder, catalogfolder, refdoc)

inbox = 'D:\\Wawrow\\2. Observational Data\\2021-03-10 - Calibrated\\HIP 46066\\LIGHT' #Image Location of .fits Format
refstars_doc = 'D:\\Reference_stars.xlsx'
refstars_csv='D:\\Reference_stars.csv' #Reference Star List
catloc = 'D:\squid\\USNOA20-All'; #Catalog #TODO Change Catalog to UCAC3 
save_loc = os.path.join(inbox, 'Outputs') # Output Folder for Files

"#TODO Initialize these once astroreduction is included"
create_master_dir = 1
run_master_bias = 1
run_master_dark = 1
run_master_flat = 1
correct_light_frames = 1
OutputsaveLoc = 0 ; #0 Default will save outputs in image folder
reduce_dir= 'D:\Image Reduction Test Images'






"Read Ref Doc"

HIP, erad, edec, vref, bvindex, vrindex, refstarsfin = ref_star_folder_read(refstars_doc)
reference_stars, ref_star_positions = astro.read_ref_stars(refstars_csv) # Reading the Reference Star Doc
f = pinpoint_init() #Start Pinpoint 

"Set Variables"
streak_array= [] #Streak Detection for Track Rate Mode
#sigma_clip = 3.5; #Clipping Factor
edge_protect = 10; #Img Edge Clipping
min_obj_pixels = 5 #Min Pixels to qualify as a Point Source
SNRLimit = 0; #Signal-To-Noise Ratio
ground_based = False
space_based = False # TODO Add Ground Based Selection Input
pinpoint = False #TODO Add Pinpoint Solve or Not to GUI
plot_results = True
TRM = False
save_plots = True
remove_large_airmass = False
image_reduce=False
"Opening Image Folder and Determing the number of files"
filepathall = getFileList(inbox); #Get List of Images

# -*- coding: utf-8 -*-
"""
Run the Space-based transform calculation.
Created on Mon May 31 12:35:34 2021

@author: jmwawrow
"""

def pinpoint_solve(inbox, catloc, max_mag, sigma, catexp, match_residual, max_solve_time, cat):
        f = pinpoint_init()
        filepathall = getFileList(inbox)
        file_suffix=".fits"
        
        for dirpath, dirnames, filenames in os.walk(inbox):
            for filename in filenames:
                if filename.endswith(file_suffix):
                    filepath = os.path.join(dirpath, filename)
                    f = win32com.client.Dispatch("Pinpoint.plate")
                    print("Processing Image: " + filepath)
                    
                    "Import Data from FITS Image"
                    imagehdularray, date, exposure_Time,imagesizeX, imagesizeY, fitsdata, filt, header,XPIXSZ, YPIXSZ,wcs = fits_header_import(filepath)
                
                    """Pinpoint Solve"""
                    if pinpoint:
                        try:
                            f.AttachFITS(filepath)
                            f.Declination = f.targetDeclination
                            f.RightAscension = f.targetRightAscension 
                            x_arcsecperpixel, y_arcsecperpixel = calc_ArcsecPerPixel(header)
                            # yBin = 4.33562092816E-004*3600;
                            # xBin =  4.33131246330E-004*3600; 
                            f.ArcsecperPixelHoriz  =  x_arcsecperpixel
                            f.ArcsecperPixelVert =  y_arcsecperpixel
                            
                            
                            "Pinpoint Solve Inputs"
                            #TODO Add Inputs for pinpoint solving to GUI
                            f.Catalog = 5
                            f.CatalogPath = catloc
                            f.CatalogMaximumMagnitude = max_mag
                            f.CatalogExpansion = catexp
                            f.SigmaAboveMean = sigma
                            f.FindImageStars
                            f.FindCatalogStars
                            f.MaxSolveTime = max_solve_time 
                            f.MaxMatchResidual = match_residual
                    
                            
                            "Pinpoint Solving"
                            f.FindCatalogStars()
                            print(f.CatalogStars.Count)
                            f.FindImageStars()
                            print(f.ImageStars.Count)
                            f.Solve()
                            #f.MatchedStars.count
                            #f.FindImageStars()
                            #print(f.ImageStars)
                            
                            
                            
                            f.DetachFITS()
                            f=None
                            print ("Pinpoint Solved")
                        except:
                            print("Could Not Solve")
                            continue

            
"Import Data from FITS Image"


"""
1. Space Based - Airmass not a factor in determining transforms
2. Ground Based - Multiple Order Transforms with Airmass and Extinctions

"""

"""Running Functions"""
       
directory = r'D:\Solved Stars\Tycho 3023_1724'
ref_stars_file = r'D:\Astro2\Reference Star Files\Reference_stars_Apr29.txt'

def Ground_based_transforms(directory, ref_stars_file):
    plot_results = True
    save_plots = True
    remove_large_airmass = False
    # file_suffix=".fit"
    # exposure_key='EXPTIME'
    # lat_key='OBSGEO-B'
    # lon_key='OBSGEO-L'
    # elev_key='OBSGEO-H'
    file_suffix=".fits"
    exposure_key='EXPTIME'
    lat_key='SITELAT'
    lon_key='SITELONG'
    elev_key='SITEELEV'
    name_key='Name'
    reference_stars, ref_star_positions = astro.read_ref_stars(ref_stars_file)
    save_loc = os.path.join(directory, 'Outputs')
    gb_final_transforms = astro._main_gb_transform_calc(directory, 
                                                        reference_stars, 
                                                        ref_star_positions,
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
    return
    
def space_based_transform(directory, ref_stars_file):
    plot_results = True
    save_plots = True
    file_suffix = "_clean_cord.fits"
    exposure_key = 'AEXPTIME'
    name_key = 'Name'
    transform_index_list = ['(B-V)', '(V-R)', '(V-I)']
    
    #for subfolder in subfolder_list:
    unique_id = ''
    directory = r'D:\NEOSSat-SA-111\clean'
    save_loc = os.path.join(directory, 'Outputs')
    
    sb_final_transform_table = astro._main_sb_transform_calc(directory, 
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

def trm_photometry():
    directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-03-20 - Calibrated\Intelsat 10-02\LIGHT'
    temp_dir = 'tmp'
    max_distance_from_sat = 20
    size = 20
    max_num_nan = 5
    plot_results = 0
    
    sats_table, uncertainty_table, sat_fwhm_table = astro._main_sc_lightcurve(directory, 
                                                                              temp_dir=temp_dir, 
                                                                              max_distance_from_sat=max_distance_from_sat, 
                                                                              size=size, 
                                                                              max_num_nan=max_num_nan, 
                                                                              plot_results=plot_results)
    
    sat_fwhm_table.pprint_all()
    uncertainty_table.pprint_all()
    sats_table.pprint_all()



def Image_reduce(reduce_dir, create_master_dark, create_master_flat,create_master_bias, create_master_dir=True):


    if os.path.isdir(reduce_dir) == False:
        raise RuntimeError('WARNING -- Directory of .fits files does not exist')
    
    # Create output directory for master files
    if create_master_dir == True:
        master_frame_dir = Path(reduce_dir, 'master_frame_data')
        master_frame_dir.mkdir(exist_ok = True)
    
    #  Select directory for master frames
    master_frame_directory = reduce_dir + '\master_frame_data'
    
    #  If a directory already exists containing the master files, uncomment the
    #  following line and place the path as a string with double backslashes.
    
    #master_frame_directory = 'C:\\pineapple\\is_a_fruit'
    
    if os.path.isdir(master_frame_directory) == False:
        raise RuntimeError('WARNING -- Directory of Master .fits files does not exist')
    
    
    # The extension to search for
    exten = '.fits'
    
    # The list of *.fits files created by the directory walk
    results = []




    
    #  Start timer
    
    start_time = time.time()
    
    # The first step is to recursively search a directory and subdirs to find all .fits files
    
    for dirpath, dirnames, files in os.walk(reduce_dir):
        for name in files:
            if name.lower().endswith(exten):
                results.append('%s' % os.path.join(dirpath, name))
    print('Have list of all .fits files')
    
    
    # %% 
    
    # Using ImageFileCollection, gather all fits files
    
    all_fits = ImageFileCollection(filenames = results)
    print('Files sorted into ImageFileCollection object')
    
    
    # %% Image Reduction 
    if run_master_bias == True:
        print('\n')
        print('Calling run_master_bias')
        IR.create_master_bias(all_fits, master_frame_directory)
    
    
    
    if run_master_dark == True:
        print('\n')
        print('Calling run_master_dark')
        IR.create_master_dark(all_fits, master_frame_directory)
    
    
    
    if run_master_flat == True:
        print('\n')
        print('Calling run_master_flat')
        IR.create_master_flat(all_fits, master_frame_directory)
    
    
    
    if correct_light_frames == True:
        print('\n')
        print('Creating output directory:', reduce_dir + '\corrected_lights')
        print('Calling correct_light_frames')
    
        # Make output directory
        correct_light_dir = Path(reduce_dir, 'corrected_lights')
        correct_light_dir.mkdir(exist_ok = True)
        correct_light_directory = reduce_dir + '\corrected_lights'
    
        #  If a specific directory is desired for the corrected light frames, uncomment the
        #  following line and place the path as a string with double backslashes.
    
        #correct_light_directory = 'C:\\apple\\is_also_a_fruit'
    
        #  Call function
        IR.correct_lights(all_fits, master_frame_directory, correct_light_directory)
    
    
    stop_time = time.time()
    elapsed_time = stop_time - start_time
    
    print('Elapsed time:', elapsed_time, 'seconds.')

    #%% NEOSSAT Dark Subtraction
    
def DarkSub(target, obspath, savedir, **kwargs):
    """Process all observations of a specific target in a specific directory."""

    # Make sure output directory exists.
    utils.ensure_dir(savedir)

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

            args = (obspath, darkfile, xsc, ysc, xov, yov, snrcut, fmax, xoff, yoff, T, bpix)

            results.append(p.apply_async(SBDarkSub.darkprocess, args=args, callback=lambda x: pbar.update()))

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
    darkname = os.path.join(savedir, darkname)
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

            args = (filename, savedir, masterdark, xsc, ysc, xov, yov, snrcut, fmax, xoff, yoff, T, bpix)

            p.apply_async(SBDarkSub.lightprocess_save, args=args, callback=lambda x: pbar.update())

        p.close()
        p.join()

    pbar.close()

    return
    
    
    
    
    
    
    
   #%% 
    #Track Rate Mode Astrometry







    #TODO Add Moffat and Gaussian Fit, SimpleSquid data collection


#         im_rms=np.std(fitsdata)
#         im_mean, im_rms = trm.BackgroundIteration(bg_rem, 0.1);
#         low_clip = im_mean + 2.5 * im_rms;
#         high_clip = 60000
#         binary_image = np.zeros((imagesizeX,imagesizeY))
#         bg_rem[bg_rem<= low_clip]

# #binary_image = (binary_image * bg_rem[bg_rem<= low_clip]) + (1 * bg_rem[bg_rem> low_clip])
#         th, im_th = cv2.threshold(bg_rem, low_clip, 1, cv2.THRESH_BINARY)
# #print(im_mean)
#         connected_image = measure.label(im_th, background=0)
# # plt.subplot(133)
# # plt.imshow(connected_image, cmap='nipy_spectral')
# # plt.axis('off')
# # plt.tight_layout()
# # plt.show()
# #im = cv2.imread(bg_rem)
# # th, im_th = cv2.threshold(im, 128, 255, cv2.THRESH_BINARY)

# #num_labels, labels_im = cv2.connectedComponents(im_th)
# #num_sourcepix = cv2.connectedComponentsWithStats(binary_image, np.array(connected_image), np.array(stats), np.array(centroids), 4, np.int(CV_32S))
#         num_sourcepix =numpy.zeros(shape=(100000,1))
#         [size_x, size_y] = imagesizeX,imagesizeY
            
#         for x in range(0,size_x):
#             for y in range(0,size_y):
#                 pixval = connected_image[x,y]
                
#                 if (pixval != 0):
#                     num_sourcepix[pixval, 0] = num_sourcepix[pixval, 0] + 1;
         
#         [valid_sources, temp] = numpy.nonzero(num_sourcepix > min_obj_pixels)
#         num_valid_sources = valid_sources.size
        
#         centroid_x = np.zeros((num_valid_sources,1));
#         centroid_y = np.zeros((num_valid_sources,1));
#         rms_x_pos = np.zeros((num_valid_sources,1));
#         rms_y_pos = np.zeros((num_valid_sources,1));
#         m11        = np.zeros((num_valid_sources,1));
#         m02        = np.zeros((num_valid_sources,1));
#         m20        = np.zeros((num_valid_sources,1));
#         ecct       = np.zeros((num_valid_sources,1));
#         compact    = np.zeros((num_valid_sources,1));
#         obj_flux   = np.zeros((num_valid_sources,1));
#         obj_max1   = np.zeros((num_valid_sources,1));
#         length     = np.zeros((num_valid_sources,1));
        
        
#         for j in range(num_valid_sources):
        
#             vsj = valid_sources[j];  
        
#             [mask_x, mask_y] = numpy.nonzero(connected_image == vsj);
#             obj_flux[j], obj_max1[j] = trm.PointSourceFluxExtraction(mask_x, mask_y, bg_rem);
            
          
#             centroid_x[j] = mean(mask_x);
#             rms_x_pos[j] = numpy.std(mask_x);
#             centroid_y[j] = mean(mask_y);
#             rms_y_pos[j] = numpy.std(mask_y);
        
#             m11[j] = MomentCalculation(mask_x, mask_y, centroid_x[j], centroid_y[j], 1,1);
#             m02[j] = MomentCalculation(mask_x, mask_y, centroid_x[j], centroid_y[j], 0,2);
#             m20[j] = MomentCalculation(mask_x, mask_y, centroid_x[j], centroid_y[j], 2,0);
#             compact[j] = Compact(num_sourcepix[vsj], m02[j], m20[j]);
#             ecct[j] = EccentricityCalculation(m11[j], m02[j], m20[j]);
        
#             x_length = (max(mask_x) - min(mask_x));
#             y_length = (max(mask_y) - min(mask_y));
#             length[j] = numpy.sqrt(x_length**2 + y_length**2);
        
                    
        
#         [compact_mean, compact_rms] = BackgroundIteration(compact,0.1);
#         [ecct_mean, ecct_rms] = BackgroundIteration(ecct,0.1)
#         compact_cut = compact_mean  + 1 * compact_rms;  
#         ecct_cut = 0.7; 
#         stars = numpy.nonzero(ecct < ecct_cut);
#         streaks = numpy.nonzero(ecct > ecct_cut);
#         stars= np.delete(stars, 1,0)
#         streaks= np.delete(streaks, 1,0)
#         sda = valid_sources[stars]
#         num_pix_in_stars = num_sourcepix[sda]
#         [mean_starpix, rms_starpix] = BackgroundIteration(num_pix_in_stars, 0.1);
        
#         pix_cutoff = mean_starpix + 10 * rms_starpix;
        
#         num_stars = stars.size;
        
#         stellar_flux_SNR = np.zeros((num_valid_sources,1));
                    
#         xmin = edge_protect;
#         xmax = imagesizeX - edge_protect;
#         ymin = edge_protect;
#         ymax = imagesizeY - edge_protect;
#         streaksize = streaks.size
#     #         st_index = connected_image(Y,X);
#     #         [mask_x, mask_y] = find(connected_image == st_index);
#     #         [mobj_flux, mobj_max] = find_point_source_flux(mask_x, mask_y, bg_rem);
#         [mask_x, mask_y] = numpy.nonzero(connected_image == vsj);
#         obj_flux[j], obj_max1[j] = trm.PointSourceFluxExtraction(mask_x, mask_y, bg_rem);
               
#         if obj_max1[j] < 60000: 
#             if (X > 10) && (X < (image_size_y-10)) && (Y > 10) && (Y < (image_size_x-10))
#     #                     %Find middle pixel value
#                          [cen_x, rms_x, cen_y, rms_y] = trm.find_weighted_centroid(mask_x, mask_y, 0*bg_rem+1);

#                            if (cen_x > 10) && (cen_x < (image_size_x-10)) && (cen_y > 10) && (cen_y < (image_size_y-10))
#                             mid_pix_val = bg_rem(round(cen_x),round(cen_y));
#                             mid_pix_valPP = bg_rem(Y,X);
#                              if vmag < mag
#                                  r = zeros(1,size(mask_x,1));
#                                  S = zeros(1,size(mask_x,1));
#                                  for q in range(1,size(mask_x,1)):
    #                                r(q) = sqrt((mask_x(q)+0.5-(mstar.Y+1))^2 + (mask_y(q)+0.5-(mstar.X+1))^2);
    #                                S(q) = bg_rem(mask_x(q),mask_y(q));
    #                             end
   
    #                             C_index = find(r==min(r),1);
    #                             r(C_index) = 0; %centroid radial value
    #                             C = S(C_index);
    
    #                             %a holds [alpha Beta] moffat parameters
    #                             %Fix a(2) Beta parameter to 1.5
    
    
    #                             fun = @(a) sum((S - (C./((1+(r.^2)/(a(1)^2)).^1.5))).^2);
    #                             aguess = 1;
    #                             [a,fminres] = fminsearch(fun,aguess);
    
    
    #                             %b holds [alpha Beta] moffat parameters
    #                             fung = @(b) sum((S - (C*exp(-(r.^2)/(2*(b^2))))).^2);
    #                             bguess = 2;
    #                             [b,fminresg] = fminsearch(fung,bguess);
    #                             %Optional plot the fits:
    #                             scatter(r,S);
    #                             E = @(a,r) (C./((1+(r.^2)/(a(1)^2)).^1.5));
    #                             hold on
    #                             ezplot(@(r)E(a,r),[0,max(r)]);  
    #                             h=get(gca,'children');
    #                             set(h(1),'color','red')
    #                             F = @(b,r) (C*exp(-(r.^2)/(2*(b^2))));
    #                             hold on
    #                             ezplot(@(r)F(b,r),[0,max(r)]);
    #                             axis([0,max(r),0,60000]);  %C+10 changed to 60000
    #                             h=get(gca,'children');
    #                             set(h(1),'color','green')
    #                 %           Output results
    #                             fprintf(starlog, [num2str(StarXY(1)) ',' num2str(StarXY(2)) ',' num2str(vmag) ',' num2str(rawflux) ',' num2str(mobj_flux) ',' num2str(mobj_max) ',' num2str(mid_pix_val) ',' ...
    #                                 num2str(mid_pix_valPP) ',' num2str(Zp) ',' num2str(ppbgmean) ',' num2str(ppbgsigma) ',' num2str(p.ImageStatisticalMode) ','  num2str(SQmean) ',' num2str(SQsigma) ','  ...
    #                                 num2str(p.ExposureInterval) ',' num2str(p.FullWidthHalfMax) ',' num2str(a(1)) ',' num2str(1.5) ',' num2str(b) ',' num2str(InstrumentalMag) ',' fpath1(i).name '\r\n']);
    #                             else %don't fit the star profile
    #                                 a = [0,0];
    #                                 b = 0;
    #                             end
    #                         end
                     
    #                     star_count = star_count +1;
    #                     pix_frac = pix_frac + mid_pix_valPP/rawflux;
    #                     if vmag < mag
    #                         mstar_count = mstar_count +1;
    #                         moffat_avg = moffat_avg + a(1);
    #                         gauss_avg = gauss_avg + b;
    #                     end
                       
    #                 end
    #             end
    #         end
    #         avg_pix_frac = pix_frac/star_count;
    #         moffat_avg = moffat_avg/mstar_count;
    #         gauss_avg = gauss_avg/mstar_count;
    #         FWHM = 2*moffat_avg*0.7664;  %derived by hand ;)
    #         DOY = str2num(fpath1(i).name(14:16))+str2num(fpath1(i).name(17:18))/24 + str2num(fpath1(i).name(19:20))/1440 + str2num(fpath1(i).name(21:22))/86400;
    #         fprintf(statslog, [num2str(DOY) ',' num2str(avg_pix_frac) ',' num2str(moffat_avg) ',' num2str(gauss_avg) ',' num2str(FWHM) ',' p.ExposureStartTime ',' num2str(p.ExposureInterval) ',' num2str(p.FullWidthHalfMax) ',' fpath1(i).name '\r\n']);
    #     end
    #        %if flag = 0

# large_stars_table = create_large_stars_table(large_table_columns, ground_based=False)
# stars_table= group_each_star(large_stars_table, ground_based=False, keys='Name')

# transform_index_list = ['(B-V)', '(V-R)', '(V-I)']
# for index in transform_index_list:
#     filter_fci, zprime_fci = space_based_transform(stars_table, plot_results=True, index=index)
#     print(f"(V-clear) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")

"""Background Subtraction"""
# TODO

"""Generating Images and Plots"""
# TODO

"""Creating PSF Profiles: Moffat and Gaussian"""
# TODO

"""Storing Data to Docs"""

        
class AstroReducer:
    
    def __init__(self, targets, start_time, stop_time, *args):
        return