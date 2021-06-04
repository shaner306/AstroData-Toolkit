import pandas as pd
import win32com.client as win32
import win32com
import os
import pywin32_system32
import math
import numpy as np
from numpy import mean
import pandas as pd
from pandas import DataFrame
import scipy
from scipy import ndimage
import skimage
from skimage import measure, filters
import matplotlib
import matplotlib.pyplot as plt
import datetime
import astropy
from astropy.stats import SigmaClip
from astropy.io import fits
import PIL
import cv2
from astropy.stats import SigmaClip
from photutils.background import SExtractorBackground
from photutils.background import MeanBackground
from photutils.background import Background2D
from photutils.background import ModeEstimatorBackground
from photutils.background import MedianBackground
from astropy import table
from astropy.coordinates import EarthLocation, AltAz, SkyCoord, match_coordinates_sky
from astropy.io import fits, ascii
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.table import Table, QTable
import astropy.units as u
from astropy.time import Time
from collections import namedtuple
from matplotlib import pyplot as plt
from photutils.detection import IRAFStarFinder
from photutils.psf import DAOGroup, BasicPSFPhotometry, IntegratedGaussianPRF
from astropy.wcs import WCS
import astropy.wcs
import re
import AstroFunctions as astro

#import TRMtester.py as trm
import numpy
import PySimpleGUI as sg

from astropy import table
from astropy.coordinates import EarthLocation, AltAz, SkyCoord, match_coordinates_sky
from astropy.io import fits, ascii
from astropy.modeling.fitting import LevMarLSQFitter, LinearLSQFitter, FittingWithOutlierRemoval
from astropy.modeling.models import Linear1D
from astropy.stats import sigma_clip, sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.table import Table, QTable, hstack
import astropy.units as u
from astropy.time import Time
from astropy.wcs import WCS
from collections import namedtuple
import ctypes
import cv2 as cv
from math import sqrt, atan
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from photutils.aperture import RectangularAperture
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from photutils.detection import IRAFStarFinder
from photutils.psf import DAOGroup, BasicPSFPhotometry, IntegratedGaussianPRF
import numpy as np
import tkinter as tk
import re
import os
from shutil import copy2, rmtree
from scipy.optimize import curve_fit

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

def linear_func(x, m, b):
    y = (m * x) + b
    return y


def init_linear_fitting(niter=3, sigma=3.0):
    fit = LevMarLSQFitter(calc_uncertainties=True)
    or_fit = FittingWithOutlierRemoval(fit, sigma_clip, niter=niter, sigma=sigma)
    line_init = Linear1D()
    return fit, or_fit, line_init

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
def calculate_img_bkg(imgdata, sigma=3.0):
    """
    Calculate the median and standard deviation of the background of the sigma clipped image.

    Parameters
    ----------
    imgdata : numpy.ndarray
        Data from the fits file.
    sigma : float
        Number of standard deviations to use for both the lower and upper clipping limit. The default is 3.0.

    Returns
    -------
    bkg : float
        Median background value of the image in ADU.
    bkg_std : float
        Standard deviation of the image background in ADU.3

    """
    _, bkg, bkg_std = sigma_clipped_stats(imgdata, sigma=3.0)
    return bkg, bkg_std

def detecting_stars(imgdata, bkg, bkg_std, fwhm=2.0):
    """
    Detect stars using IRAFStarFinder.

    Parameters
    ----------
    imgdata : numpy.ndarray
        Data from the fits file.
    fwhm : float, optional
        Initial guess of FWHM of the image stars in pixels. The default is 2.0.
    bkg : float
        Background value of the image in ADU.
    bkg_std : float
        Standard deviation of the image background in ADU.

    Returns
    -------
    irafsources : astropy.table.Table
        Table containing information of all stars detected in the image.
        Has columns:
            id,
            xcentroid,
            ycentroid,
            fwhm,
            sharpness,
            roundness,
            pa,
            npix,
            sky,
            peak,
            flux,
            mag

    """
    iraffind = IRAFStarFinder(threshold=bkg+3*bkg_std, fwhm=fwhm)
    irafsources = iraffind(imgdata - bkg)
    return irafsources

def get_instr_filter_name(hdr, filter_key='FILTER'):
    """
    Get the name of the filter used when taking the image.

    Parameters
    ----------
    hdr : astropy.io.fits.header.Header
        Header from the fits file.
    filter_key : string, optional
        Keyword of the entry in the FITS header that contains the filter information. The default is 'FILTER'.

    Returns
    -------
    instr_filter : string
        Instrumental filter band of the image.

    """
    instr_filter = hdr[filter_key].lower()
    return instr_filter

def BackgroundEstimationMulti(fitsdata, sigma_clip, bkgmethod, printval):
   
    sigma_clip = SigmaClip(sigma=2.5)
    bkg = SExtractorBackground(sigma_clip)
    bkg_value1 = bkg.calc_background(fitsdata)
    
    #print(bkg_value)
    bkg = MeanBackground(sigma_clip)
    bkg_value2 = bkg.calc_background(fitsdata)
    
    bkg = MedianBackground(sigma_clip)
    bkg_value3 = bkg.calc_background(fitsdata)
    
    bkg = ModeEstimatorBackground(sigma_clip)
    bkg_value4 = bkg.calc_background(fitsdata)
    
    
    bkg_estimator1 = SExtractorBackground()
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
def convert_pixel_to_ra_dec(irafsources,wcs):
    """
    Convert the locations of the sources from pixels to RA/dec.

    Parameters
    ----------
    irafsources : astropy.table.Table
        Table containing information of all stars detected in the image.
        Has columns:
            id,
            xcentroid,
            ycentroid,
            fwhm,
            sharpness,
            roundness,
            pa,
            npix,
            sky,
            peak,
            flux,
            mag
    wcs : astropy.wcs.wcs.WCS
        AstroPy World Coordinate System object for the image.

    Returns
    -------
    skypositions : astropy.coordinates.sky_coordinate.SkyCoord
        AstroPy SkyCoord object containing the RA/dec positions of all sources in the image.

    """
    skypositions = wcs.pixel_to_world(irafsources['xcentroid'], irafsources['ycentroid'])
    return skypositions


def convert_ra_dec_to_alt_az(skypositions, hdr, lat_key='SITELAT', lon_key='SITELONG', elev_key='SITEELEV'):
    """
    Convert RA/dec locations to Altitude/Azimuth (Azimuth/Elevation).

    Parameters
    ----------
    skypositions : astropy.coordinates.sky_coordinate.SkyCoord
        AstroPy SkyCoord object containing the RA/dec position(s) to be converted to Alt/Az.
    hdr : astropy.io.fits.header.Header
        Header from the fits file.

    Returns
    -------
    altazpositions : astropy.coordinates.sky_coordinate.SkyCoord
        Alt/Az position(s) of the skyposition(s)

    """
    obstime = Time(hdr['DATE-OBS'], format='fits')
    lat = hdr[lat_key]
    lon = hdr[lon_key]
    height = hdr[elev_key]
    location = EarthLocation.from_geodetic(lon=lon, lat=lat, height=height)
    current_aa = AltAz(location=location, obstime=obstime)
    altazpositions = skypositions.transform_to(current_aa)
    return altazpositions


def calculate_fwhm(irafsources):
    """
    Calculate the mean and standard deviation FWHM of all sources in the image in pixels.

    Parameters
    ----------
    irafsources : astropy.table.Table
        Table containing information of all stars detected in the image.
        Has columns:
            id,
            xcentroid,
            ycentroid,
            fwhm,
            sharpness,
            roundness,
            pa,
            npix,
            sky,
            peak,
            flux,
            mag

    Returns
    -------
    fwhm : float
        Mean FWHM of all sources in the image.
    fwhm_std : float
        Standard deviation of the FWHM of all sources in the image.

    """
    fwhms = np.array(irafsources['fwhm'])
    fwhm = fwhms.mean()
    fwhm_std = fwhms.std()
    return fwhm, fwhm_std


def perform_photometry(irafsources, fwhm, imgdata, bkg, fitter=LevMarLSQFitter(), fitshape=25):
    """
    Perform PSF photometry on all sources in a selected image.

    Parameters
    ----------
    irafsources : astropy.table.Table
        Table containing information of all stars detected in the image.
        Has columns:
            id,
            xcentroid,
            ycentroid,
            fwhm,
            sharpness,
            roundness,
            pa,
            npix,
            sky,
            peak,
            flux,
            mag
    fwhm : float
        Mean FWHM of all sources in the image.
    imgdata : numpy.ndarray
        Data from the fits file.
    bkg : float, optional
        Background value of the image in ADU. The default is None.
    bkg_estimator : callable, instance of any photutils.background.BackgroundBase subclass, optional
        bkg_estimator should be able to compute either a scalar background or a 2D background of a given 2D image. 
        If None, no background subtraction is performed. The default is None.
    fitter : astropy.modeling.fitting.Fitter instance, optional
        Fitter object used to compute the optimized centroid positions and/or flux of the identified sources. 
        The default is LevMarLSQFitter().
    fitshape : int or length-2 array-like, optional
        Rectangular shape around the center of a star which will be used to collect the data to do the fitting. 
        Can be an integer to be the same along both axes. For example, 5 is the same as (5, 5), 
        which means to fit only at the following relative pixel positions: [-2, -1, 0, 1, 2]. 
        Each element of fitshape must be an odd number. The default is 25.

    Returns
    -------
    photometry_result : astropy.table.Table or None
        Table with the photometry results, i.e., centroids and fluxes estimations and the initial estimates used to 
        start the fitting process. Uncertainties on the fitted parameters are reported as columns called 
        <paramname>_unc provided that the fitter object contains a dictionary called fit_info with the key param_cov, 
        which contains the covariance matrix. If param_cov is not present, uncertanties are not reported.

    """
    daogroup = DAOGroup(2 * fwhm)
    psf_model = IntegratedGaussianPRF(sigma=fwhm * gaussian_fwhm_to_sigma)
    psf_model.x_0.fixed = True
    psf_model.y_0.fixed = True
    pos = Table(names=['x_0', 'y_0', 'flux_0'],
                data=[irafsources['xcentroid'], irafsources['ycentroid'], irafsources['flux']])
    bkg_estimator = MedianBackground()
    photometry = BasicPSFPhotometry(group_maker=daogroup, 
                                    bkg_estimator=None, 
                                    psf_model=psf_model, 
                                    fitter=fitter,
                                    fitshape=fitshape)
    photometry_result = photometry(image=imgdata - bkg, init_guesses=pos)
    return photometry_result



def normalize_flux_by_time(fluxes_tab, exptime):
    """
    Calculate flux per 1 second of time.

    Parameters
    ----------
    fluxes_tab : array-like
        1D array-like containing the fluxes for all sources in the image in counts.
    exptime : float
        Expsure time of the image in seconds.

    Returns
    -------
    fluxes : numpy array
        Array containing the time normalized fluxes with units counts / second.

    """
    fluxes_units = np.array(fluxes_tab) * u.ct
    exptime_units = exptime * u.s
    fluxes = fluxes_units / exptime_units
    return fluxes


def calculate_magnitudes(photometry_result, exptime):
    """
    Convert the flux of all sources in the image to instrumental magnitudes.

    Parameters
    ----------
    photometry_result : astropy.table.Table
        Table with the photometry results, i.e., centroids and fluxes estimations and the initial estimates used to 
        start the fitting process.
    exptime : float
        Exposure time of the image in seconds.

    Returns
    -------
    instr_mags : numpy array
        Array containing the instrumental magnitudes of the sources in the image.

    """
    fluxes = normalize_flux_by_time(photometry_result['flux_fit'], exptime)
    instr_mags_units = u.Magnitude(fluxes)
    instr_mags = instr_mags_units.value
    return instr_mags


def calculate_magnitudes_sigma(photometry_result, exptime):
    """
    Calculate the standard deviation of the instrumental magnitudes of the sources in the image.

    Parameters
    ----------
    photometry_result : astropy.table.Table
        Table with the photometry results, i.e., centroids, fluxes, amd flux uncertainty estimations and the initial 
        estimates used to start the fitting process.
    exptime : float
        Expsure time of the image in seconds.

    Returns
    -------
    instr_mags_sigma : numpy array
        Array containing the standard deviation of the instrumental magnitudes of the sources in the image.

    """
    fluxes = normalize_flux_by_time(photometry_result['flux_fit'], exptime)
    flux_uncs = normalize_flux_by_time(photometry_result['flux_unc'], exptime)
    snr = (fluxes / flux_uncs).value
    instr_mags_sigma = 1.0857 / np.sqrt(snr)
    return instr_mags_sigma
def find_ref_stars(reference_stars, 
                   ref_star_positions, 
                   skypositions, 
                   instr_mags, 
                   instr_mags_sigma, 
                   fluxes, 
                   ground_based=False, 
                   altazpositions=None, 
                   max_ref_sep=10.0):
    """
    Match the stars detected in the image to those provided in the reference star file.
    
    Parameters
    ----------
    reference_stars : astropy.table.Table
        Table with the data extracted from ref_stars_file.
    ref_star_positions : astropy.coordinates.sky_coordinate.SkyCoord
        AstroPy SkyCoord object containing the RA/dec positions of all reference stars in the file.
    skypositions : astropy.coordinates.sky_coordinate.SkyCoord
        AstroPy SkyCoord object containing the RA/dec positions of all sources in the image.
    instr_mags : numpy array
        Array containing the instrumental magnitudes of the sources in the image.
    instr_mags_sigma : numpy array
        Array containing the standard deviation of the instrumental magnitudes of the sources in the image.
    fluxes : numpy array
        Array containing the non-normalized fluxes of all sources.
    ground_based : bool, optional
        Whether or not the image is from a ground-based sensor. The default is False.
    altazpositions : astropy.coordinates.sky_coordinate.SkyCoord, optional
        Alt/Az position(s) of the skyposition(s). Must be provided if ground_based is true. The default is None.
    max_ref_sep : float, optional
        Maximum angular separtation to consider an image star to be a reference star in arcsec. The default is 10.0.

    Returns
    -------
    matched_stars : namedtuple
        Attributes:
            ref_star_index : array-like or int
                Index of the star(s) in ref_stars_file that correspond to a star in the image.
            img_star_index : array-like or int
                Index of the star(s) detected in the image that correspond to a star in ref_stars_file.
            ref_star : astropy.table.Table
                Rows from ref_stars_file that correspond to a matched reference star.
            ref_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched reference star(s) from ref_stars_file.
            img_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched star(s) detected in the image.
            ang_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the image and ref_stars_file.
            img_instr_mag : numpy array
                Array containing the instrumental magnitudes of the matched star(s) detected in the image.
            img_instr_mag_sigma : numpy array
                Array containing the standard deviations of the instrumental magnitudes of the matched star(s) 
                detected in the image.
            flux : numpy array
                Array containing the non-normalized fluxes of the matched star(s) detected in the image.
            img_star_altaz : astropy.coordinates.sky_coordinate.SkyCoord
                Alt/Az of the matched star(s) detected in the image. None if ground_based is False.
            img_star_airmass : float
                sec(z) of img_star_altaz. None if ground_based is False.
    None : if no stars are matched.

    """
    max_ref_sep=10.0
    # reference_stars, ref_star_positions = read_ref_stars(ref_stars_file)
    max_ref_sep = max_ref_sep * u.arcsec
    
    idx, sep2d, _ = match_coordinates_sky(ref_star_positions, skypositions)
    ref_star_index = np.where(sep2d < max_ref_sep)
    
    if len(ref_star_index[0]) == 0:
        print("No reference star detected in the image.")
        return
    elif len(ref_star_index[0]) == 1:
        ref_star_index = int(ref_star_index[0])
    else:
        ref_star_index = ref_star_index[0]
    img_star_index = idx[ref_star_index]
    ref_star = reference_stars[ref_star_index]
    ref_star_loc = ref_star_positions[ref_star_index]
    img_star_loc = skypositions[img_star_index]
    ang_separation = sep2d[ref_star_index]
    img_instr_mag = instr_mags[img_star_index]
    img_instr_mag_sigma = instr_mags_sigma[img_star_index]
    img_fluxes = fluxes[img_star_index]
    img_star_altaz = None
    img_star_airmass = None
    if ground_based:
        if not altazpositions:
            print('altazpositions should be provided if ground_based is True.')
            return
        img_star_altaz = altazpositions[img_star_index]
        img_star_airmass = img_star_altaz.secz.value
    matched_stars = namedtuple('matched_stars', 
                               ['ref_star_index',
                                'img_star_index',
                                'ref_star',
                                'ref_star_loc',
                                'img_star_loc',
                                'ang_separation',
                                'img_instr_mag',
                                'img_instr_mag_sigma',
                                'flux',
                                'img_star_altaz',
                                'img_star_airmass'])
    return matched_stars(
        ref_star_index, 
        img_star_index, 
        ref_star, 
        ref_star_loc, 
        img_star_loc, 
        ang_separation, 
        img_instr_mag, 
        img_instr_mag_sigma,
        img_fluxes,
        img_star_altaz, 
        img_star_airmass)



def get_field_name(matched_stars, name_key='Name'):
    """
    Get the name of the field of the matched stars in the image.

    Parameters
    ----------
    matched_stars : namedtuple
        Attributes:
            ref_star_index : array-like or int
                Index of the star(s) in ref_stars_file that correspond to a star in the image.
            img_star_index : array-like or int
                Index of the star(s) detected in the image that correspond to a star in ref_stars_file.
            ref_star : astropy.table.Table
                Rows from ref_stars_file that correspond to a matched reference star.
            ref_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched reference star(s) from ref_stars_file.
            img_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched star(s) detected in the image.
            ang_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the image and ref_stars_file.
            img_instr_mag : numpy array
                Array containing the instrumental magnitudes of the matched star(s) detected in the image.
            img_instr_mag_sigma : numpy array
                Array containing the standard deviations of the instrumental magnitudes of the matched star(s) 
                detected in the image.
            flux : numpy array
                Array containing the non-normalized fluxes of the matched star(s) detected in the image.
            img_star_altaz : astropy.coordinates.sky_coordinate.SkyCoord
                Alt/Az of the matched star(s) detected in the image. None if ground_based is False.
            img_star_airmass : float
                sec(z) of img_star_altaz. None if ground_based is False.
    name_key : string, optional
        The column name of the unique identifier/star name in matched_stars.ref_star. The default is 'Name'.

    Returns
    -------
    field : string
        Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").

    """
    try:
        num_stars = len(matched_stars.img_instr_mag)
    except TypeError:
        num_stars = 1
    if num_stars > 1:
        for row in matched_stars.ref_star:
            split_string = re.split('[^a-zA-Z0-9]', str(row[name_key]))
            if len(split_string) > 1:
                return ' '.join(split_string[:-1])
            elif len(split_string) == 1:
                return split_string[0]
            else:
                return
    elif num_stars == 1:
        split_string = re.split('[^a-zA-Z0-9]', str(matched_stars.ref_star[name_key]))
        if len(split_string) > 1:
            return ' '.join(split_string[:-1])
        elif len(split_string) == 1:
            return split_string[0]
        else:
            return
    else:
        return
def init_large_table_columns():
    """
    Create all of the columns that will be turned into the large stars table.

    Returns
    -------
    large_table_columns : namedtuple
        Attributes:
            field : empty list
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            ref_star_name : empty list
                Name/unique identifier of the reference star.
            times : empty list
                Time of the observation.
            flux_table : empty list
                Flux of the source.
            exposure : empty list
                Exposure time of the image.
            ref_star_RA : empty list
                Right ascension of the reference star(s) from the file.
            ref_star_dec : empty list
                Declination of the reference star(s) from the file.
            img_star_RA : empty list
                Right ascension of the reference star(s) found in the image.
            img_star_dec : empty list
                Declination of the reference star(s) found in the image.
            angular_separation : empty list
                Angular distance between the matched star(s) detected in the image and ref_stars_file.
            img_star_mag : empty list
                Instrumental magnitude of the reference star(s) found in the image.
            img_star_mag_sigma : empty list
                Standard deviation of the instrumental magnitude of the reference star(s) found in the image.
            filters : empty list
                Filter used during for the image.
            V_apparents : empty list
                Apparent V magnitude from the reference file.
            B_V_apparents : empty list
                Apparent B-V colour index from the reference file.
            U_B_apparents : empty list
                Apparent U-B colour index from the reference file.
            V_R_apparents : empty list
                Apparent V-R colour index from the reference file.
            V_I_apparents : empty list
                Apparent V-I colour index from the reference file.
            V_sigma_apparents : empty list
                Standard deviation of the apparent V magnitude from the reference file.
            img_star_airmass : empty list
                Sec(z) of the reference star(s) found in the image.
            
    """
    field = []
    ref_star_name = []
    times = []
    flux_table = []
    exposure = []
    ref_star_RA = []
    ref_star_dec = []
    img_star_RA = []
    img_star_dec = []
    angular_separation = []
    V_apparents = []
    img_star_mag = []
    img_star_mag_sigma = []
    filters = []
    B_V_apparents = []
    U_B_apparents = []
    V_R_apparents = []
    V_I_apparents = []
    V_sigma_apparents = []
    img_star_airmass = []
    # X_rounded = []
    large_table_columns = namedtuple('large_table_columns',
                                     ['field',
                                      'ref_star_name',
                                      'times',
                                      'flux_table',
                                      'exposure',
                                      'ref_star_RA',
                                      'ref_star_dec',
                                      'img_star_RA',
                                      'img_star_dec',
                                      'angular_separation',
                                      'img_star_mag',
                                      'img_star_mag_sigma',
                                      'filters',
                                      'V_apparents',
                                      'B_V_apparents',
                                      'U_B_apparents',
                                      'V_R_apparents',
                                      'V_I_apparents',
                                      'V_sigma_apparents',
                                      'img_star_airmass',
                                      # 'X_rounded'
                                      ])
    return large_table_columns(field,
                               ref_star_name, 
                               times, 
                               flux_table, 
                               exposure, 
                               ref_star_RA, 
                               ref_star_dec, 
                               img_star_RA, 
                               img_star_dec, 
                               angular_separation, 
                               img_star_mag, 
                               img_star_mag_sigma, 
                               filters, 
                               V_apparents, 
                               B_V_apparents, 
                               U_B_apparents, 
                               V_R_apparents, 
                               V_I_apparents, 
                               V_sigma_apparents, 
                               img_star_airmass, 
                               # X_rounded
                               )


def update_large_table_columns(large_table_columns, matched_stars, header, exptime,ref_star, ground_based=False,  name_key='Name'):
    """
    Update columns to be used for the large stars table based on information from the current image.

    Parameters
    ----------
    large_table_columns : namedtuple
        Attributes:
            field : string list
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            ref_star_name : string list
                Name/unique identifier of the reference star.
            times : numpy.float64
                Time of the observation.
            flux_table : numpy.float64
                Flux of the source.
            exposure : numpy.float64
                Exposure time of the image.
            ref_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) from the file in unit hourangle.
            ref_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) from the file in unit degree.
            img_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) found in the image in unit hourangle.
            img_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) found in the image in unit degree.
            angular_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the image and ref_stars_file.
            img_star_mag : numpy.float64
                Instrumental magnitude of the reference star(s) found in the image.
            img_star_mag_sigma : numpy.float64
                Standard deviation of the instrumental magnitude of the reference star(s) found in the image.
            filters : string list
                Filter used during for the image.
            V_apparents : numpy.float64
                Apparent V magnitude from the reference file.
            B_V_apparents : numpy.float64
                Apparent B-V colour index from the reference file.
            U_B_apparents : numpy.float64
                Apparent U-B colour index from the reference file.
            V_R_apparents : numpy.float64
                Apparent V-R colour index from the reference file.
            V_I_apparents : numpy.float64
                Apparent V-I colour index from the reference file.
            V_sigma_apparents : numpy.float64
                Standard deviation of the apparent V magnitude from the reference file.
            img_star_airmass : numpy.float64
                Sec(z) of the reference star(s) found in the image.
    matched_stars : namedtuple
        Attributes:
            ref_star_index : array-like or int
                Index of the star(s) in ref_stars_file that correspond to a star in the image.
            img_star_index : array-like or int
                Index of the star(s) detected in the image that correspond to a star in ref_stars_file.
            ref_star : astropy.table.Table
                Rows from ref_stars_file that correspond to a matched reference star.
            ref_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched reference star(s) from ref_stars_file.
            img_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched star(s) detected in the image.
            ang_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the image and ref_stars_file.
            img_instr_mag : numpy array
                Array containing the instrumental magnitudes of the matched star(s) detected in the image.
            img_instr_mag_sigma : numpy array
                Array containing the standard deviations of the instrumental magnitudes of the matched star(s) 
                detected in the image.
            flux : numpy array
                Array containing the non-normalized fluxes of the matched star(s) detected in the image.
            img_star_altaz : astropy.coordinates.sky_coordinate.SkyCoord
                Alt/Az of the matched star(s) detected in the image. None if ground_based is False.
            img_star_airmass : float
                sec(z) of img_star_altaz. None if ground_based is False.
    hdr : astropy.io.fits.header.Header
        Header from the fits file.
    exptime : float
        Expsure time of the image in seconds.
    ground_based : bool, optional
        Whether or not the image is from a ground-based sensor. The default is False.
    name_key : string, optional
        The column name of the unique identifier/star name in matched_stars.ref_star. The default is 'Name'.

    Returns
    -------
    updated_large_table_columns : namedtuple
        Attributes:
            field : string list
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            ref_star_name : string list
                Name/unique identifier of the reference star.
            times : numpy.float64
                Time of the observation.
            flux_table : numpy.float64
                Flux of the source.
            exposure : numpy.float64
                Exposure time of the image.
            ref_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) from the file in unit hourangle.
            ref_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) from the file in unit degree.
            img_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) found in the image in unit hourangle.
            img_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) found in the image in unit degree.
            angular_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the image and ref_stars_file.
            img_star_mag : numpy.float64
                Instrumental magnitude of the reference star(s) found in the image.
            img_star_mag_sigma : numpy.float64
                Standard deviation of the instrumental magnitude of the reference star(s) found in the image.
            filters : string list
                Filter used during for the image.
            V_apparents : numpy.float64
                Apparent V magnitude from the reference file.
            B_V_apparents : numpy.float64
                Apparent B-V colour index from the reference file.
            U_B_apparents : numpy.float64
                Apparent U-B colour index from the reference file.
            V_R_apparents : numpy.float64
                Apparent V-R colour index from the reference file.
            V_I_apparents : numpy.float64
                Apparent V-I colour index from the reference file.
            V_sigma_apparents : numpy.float64
                Standard deviation of the apparent V magnitude from the reference file.
            img_star_airmass : numpy.float64
                Sec(z) of the reference star(s) found in the image. Only output if ground_based is True.

    """
    #exptime= 5.0
    updated_large_table_columns = large_table_columns
    try:
        num_stars = len(matched_stars.img_instr_mag)
    except TypeError:
        num_stars = 1
    if num_stars > 1:
        for row in matched_stars.ref_star:
            split_string = re.split('[^a-zA-Z0-9]', str(row['HIP']))
            if len(split_string) > 1:
                updated_large_table_columns.field.append(' '.join(split_string[:-1]))
            elif len(split_string) == 1:
                updated_large_table_columns.field.append(split_string[0])
            else:
                print('Could not find the name of the field.')
                
                
                
        updated_large_table_columns.ref_star_name.extend(matched_stars.ref_star['HIP'])
        updated_large_table_columns.flux_table.extend(matched_stars.flux)
        time = Time(header['DATE-OBS'], format='fits')
        time_repeat = np.full(num_stars, time.jd)
        updated_large_table_columns.times.extend(time_repeat)
        exposure_repeat = np.full(num_stars, exptime)
        updated_large_table_columns.exposure.extend(exposure_repeat)
        updated_large_table_columns.ref_star_RA.extend(matched_stars.ref_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.ref_star_dec.extend(matched_stars.ref_star_loc.dec)
        updated_large_table_columns.img_star_RA.extend(matched_stars.img_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.img_star_dec.extend(matched_stars.img_star_loc.dec)
        updated_large_table_columns.angular_separation.extend(matched_stars.ang_separation.to(u.arcsec))
        updated_large_table_columns.img_star_mag.extend(matched_stars.img_instr_mag)
        updated_large_table_columns.img_star_mag_sigma.extend(matched_stars.img_instr_mag_sigma)
        filter_name_repeat = np.full(num_stars, header['FILTER'])
        updated_large_table_columns.filters.extend(filter_name_repeat)
        updated_large_table_columns.V_apparents.extend(ref_star['V'])
        try:
            updated_large_table_columns.B_V_apparents.extend(ref_star['(B-V)'])
            updated_large_table_columns.U_B_apparents.extend(ref_star['(U-B)'])
            updated_large_table_columns.V_R_apparents.extend(ref_star['(V-R)'])
            updated_large_table_columns.V_I_apparents.extend(ref_star['(V-I)'])
            updated_large_table_columns.V_sigma_apparents.extend(ref_star['V_sigma'])
        except KeyError:
            updated_large_table_columns.B_V_apparents.extend(ref_star['B-V'])
            updated_large_table_columns.U_B_apparents.extend(ref_star['U-B'])
            updated_large_table_columns.V_R_apparents.extend(ref_star['V-R'])
            updated_large_table_columns.V_I_apparents.extend(ref_star['V-I'])
            updated_large_table_columns.V_sigma_apparents.extend(ref_star['e_V'])
        if not ground_based:
            return updated_large_table_columns
        updated_large_table_columns.img_star_airmass.extend(matched_stars.img_star_airmass)
        # updated_large_table_columns.X_rounded.extend(round(matched_stars.img_star_airmass, 1))
    elif num_stars == 1:
        
        star= matched_stars.ref_star['HIP']
        try:
            split_string = re.split('[^a-zA-Z0-9]', str(star[name_key]))
        except:
            split_string = []
            
        if len(split_string) > 1:
            updated_large_table_columns.field.append(' '.join(split_string[:-1]))
        elif len(split_string) == 1:
            updated_large_table_columns.field.append(split_string[0])
        else:
            print('Could not find the name of the field.')
            updated_large_table_columns.field.append("Not Found")
        updated_large_table_columns.ref_star_name.append(matched_stars.ref_star['HIP'])   
        updated_large_table_columns.flux_table.append(matched_stars.flux)
        time = Time(header['DATE-OBS'], format='fits')
        
        updated_large_table_columns.times.append(time.jd)
        updated_large_table_columns.exposure.append(exptime)
        updated_large_table_columns.ref_star_RA.append(matched_stars.ref_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.ref_star_dec.append(matched_stars.ref_star_loc.dec)
        updated_large_table_columns.img_star_RA.append(matched_stars.img_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.img_star_dec.append(matched_stars.img_star_loc.dec)
        updated_large_table_columns.angular_separation.append(matched_stars.ang_separation.to(u.arcsec))
        updated_large_table_columns.img_star_mag.append(matched_stars.img_instr_mag)
        updated_large_table_columns.img_star_mag_sigma.append(matched_stars.img_instr_mag_sigma)
        updated_large_table_columns.filters.append(header['FILTER'])
        updated_large_table_columns.V_apparents.append(ref_star['V'])
        try:
            
            updated_large_table_columns.B_V_apparents.append(ref_star['(B-V)'])
            updated_large_table_columns.U_B_apparents.append(ref_star['(U-B)'])
            updated_large_table_columns.V_R_apparents.append(ref_star['(V-R)'])
            updated_large_table_columns.V_I_apparents.append(ref_star['(V-I)'])
            updated_large_table_columns.V_sigma_apparents.append(ref_star['V_sigma'])
        except KeyError:
            updated_large_table_columns.B_V_apparents.append(ref_star['B-V'])
            updated_large_table_columns.U_B_apparents.append(ref_star['U-B'])
            updated_large_table_columns.V_R_apparents.append(ref_star['V-R'])
            updated_large_table_columns.V_I_apparents.append(ref_star['V-I'])
            updated_large_table_columns.V_sigma_apparents.append(ref_star['e_V'])
        if not ground_based:
            return updated_large_table_columns
        updated_large_table_columns.img_star_airmass.append(matched_stars.img_star_airmass)
        # updated_large_table_columns.X_rounded.append(round(matched_stars.img_star_airmass, 1))
    else:
        return
    return updated_large_table_columns


def create_large_stars_table(large_table_columns, ground_based=False):
    """
    Convert the large table columns into an AstroPy.table QTable.

    Parameters
    ----------
    large_table_columns : namedtuple
        Attributes:
            field : string list
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            ref_star_name : string list
                Name/unique identifier of the reference star.
            times : numpy.float64
                Time of the observation.
            flux_table : numpy.float64
                Flux of the source.
            exposure : numpy.float64
                Exposure time of the image.
            ref_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) from the file in unit hourangle.
            ref_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) from the file in unit degree.
            img_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) found in the image in unit hourangle.
            img_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) found in the image in unit degree.
            angular_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the image and ref_stars_file.
            img_star_mag : numpy.float64
                Instrumental magnitude of the reference star(s) found in the image.
            img_star_mag_sigma : numpy.float64
                Standard deviation of the instrumental magnitude of the reference star(s) found in the image.
            filters : string list
                Filter used during for the image.
            V_apparents : numpy.float64
                Apparent V magnitude from the reference file.
            B_V_apparents : numpy.float64
                Apparent B-V colour index from the reference file.
            U_B_apparents : numpy.float64
                Apparent U-B colour index from the reference file.
            V_R_apparents : numpy.float64
                Apparent V-R colour index from the reference file.
            V_I_apparents : numpy.float64
                Apparent V-I colour index from the reference file.
            V_sigma_apparents : numpy.float64
                Standard deviation of the apparent V magnitude from the reference file.
            img_star_airmass : numpy.float64
                Sec(z) of the reference star(s) found in the image.
    ground_based : bool, optional
        Whether or not the image is from a ground-based sensor. The default is False.

    Returns
    -------
    large_stars_table : astropy.table.table.QTable
        Table containing all information on the detected stars. Has columns:
            Field : string
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            Name : string
                Name/unique identifier of the reference star.
            Time (JD) : numpy.float64
                Time of the observation.
            RA_ref : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) from the file in unit hourangle.
            dec_ref : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) from the file in unit degree.
            RA_img : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) found in the image in unit hourangle.
            dec_img : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) found in the image in unit degree.
            angular_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the image and ref_stars_file.
            flux : numpy.float64
                Flux of the source.
            exposure : numpy.float64
                Exposure time of the image.
            mag_instrumental : numpy.float64
                Instrumental magnitude of the reference star(s) found in the image.
            mag_instrumental_sigma : numpy.float64
                Standard deviation of the instrumental magnitude of the reference star(s) found in the image.
            filter : string
                Filter used during for the image.
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
                Standard deviation of the apparent V magnitude from the reference file.
            X : numpy.float64
                Sec(z) of the reference star(s) found in the image. Only output if ground_based is True.

    """
    if ground_based:
        large_stars_table = QTable(
            data=[
                large_table_columns.field,
                large_table_columns.ref_star_name, 
                large_table_columns.times, 
                large_table_columns.ref_star_RA, 
                large_table_columns.ref_star_dec, 
                large_table_columns.img_star_RA, 
                large_table_columns.img_star_dec, 
                large_table_columns.angular_separation, 
                large_table_columns.flux_table, 
                large_table_columns.exposure, 
                large_table_columns.img_star_mag, 
                large_table_columns.img_star_mag_sigma, 
                large_table_columns.filters, 
                large_table_columns.V_apparents, 
                large_table_columns.B_V_apparents, 
                large_table_columns.U_B_apparents, 
                large_table_columns.V_R_apparents, 
                large_table_columns.V_I_apparents, 
                large_table_columns.V_sigma_apparents, 
                #large_table_columns.img_star_airmass, 
                # large_table_columns.X_rounded
                ],
            names=[
                'Field',
                'Name',
                'Time (JD)',
                'RA_ref',
                'dec_ref',
                'RA_img',
                'dec_img',
                'angular_separation',
                'flux',
                'exposure',
                'mag_instrumental',
                'mag_instrumental_sigma',
                'filter',
                'V',
                '(B-V)',
                '(U-B)',
                '(V-R)',
                '(V-I)',
                'V_sigma',
                #'X',
                # 'X_rounded'
                ]
            )
    else:
        large_stars_table = QTable(
            data=[
                large_table_columns.field,
                large_table_columns.ref_star_name, 
                large_table_columns.times, 
                large_table_columns.ref_star_RA, 
                large_table_columns.ref_star_dec, 
                large_table_columns.img_star_RA, 
                large_table_columns.img_star_dec, 
                large_table_columns.angular_separation, 
                large_table_columns.flux_table, 
                large_table_columns.exposure, 
                large_table_columns.img_star_mag, 
                large_table_columns.img_star_mag_sigma, 
                large_table_columns.filters, 
                large_table_columns.V_apparents, 
                large_table_columns.B_V_apparents, 
                large_table_columns.U_B_apparents, 
                large_table_columns.V_R_apparents, 
                large_table_columns.V_I_apparents, 
                large_table_columns.V_sigma_apparents,
                ],
            names=[
                'Field',
                'Name',
                'Time (JD)',
                'RA_ref',
                'dec_ref',
                'RA_img',
                'dec_img',
                'angular_separation',
                'flux',
                'exposure',
                'mag_instrumental',
                'mag_instrumental_sigma',
                'filter',
                'V',
                '(B-V)',
                '(U-B)',
                '(V-R)',
                '(V-I)',
                'V_sigma',
                ]
            )
    return large_stars_table


def group_each_star(large_stars_table, vref, ground_based=False, keys='Name'):
    """
    Group all detections of unique reference stars together.

    Parameters
    ----------
    large_stars_table : astropy.table.table.QTable
        Table containing all information on the detected stars. Has columns:
            Field : string
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            Name : string
                Name/unique identifier of the reference star.
            Time (JD) : numpy.float64
                Time of the observation.
            RA_ref : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) from the file in unit hourangle.
            dec_ref : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) from the file in unit degree.
            RA_img : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) found in the image in unit hourangle.
            dec_img : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) found in the image in unit degree.
            angular_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the image and ref_stars_file.
            flux : numpy.float64
                Flux of the source.
            exposure : numpy.float64
                Exposure time of the image.
            mag_instrumental : numpy.float64
                Instrumental magnitude of the reference star(s) found in the image.
            mag_instrumental_sigma : numpy.float64
                Standard deviation of the instrumental magnitude of the reference star(s) found in the image.
            filter : string
                Filter used during for the image.
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
                Standard deviation of the apparent V magnitude from the reference file.
            X : numpy.float64
                Sec(z) of the reference star(s) found in the image.
    ground_based : bool, optional
        Whether or not the image is from a ground-based sensor. The default is False.
    keys : string, optional
        Table column(s) to use to group the different reference stars. The default is 'Name'.

    Returns
    -------
    stars_table : astropy.table.table.Table
        Table containing the mean of the important information for each star. Has columns:
            Field : string
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
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
                Standard deviation of the apparent V magnitude from the reference file.
            <filter> : numpy.float64
                Mean instrumental magnitude of all detections of the star in <filter>. There is a different column for 
                each different filter used across the images.
            <filter>_sigma : numpy.float64
                Standard deviation of the instrumental magnitudes of all detections of the star in <filter>. 
                There is a different column for each different filter used across the images.
            X_<filter> : numpy.float64
                Mean airmass of all detections of the star in <filter>. There is a different column for each different 
                filter used across the images. Only output if ground_based is True.
            X_<filter>_sigma : numpy.float64
                Standard deviation of the airmasses of all detections of the star in <filter>. There is a different 
                column for each different filter used across the images. Only output if ground_based is True.

    """
    
    unique_stars = table.unique(large_stars_table, keys=keys)
    N = len(unique_stars)
    nan_array = np.empty(N)
    nan_array.fill(np.nan)
    apparent_mags_table = Table(
        names=[
            'Field',
            'Name',
            'V',
            '(B-V)',
            '(U-B)',
            '(V-R)',
            '(V-I)',
            'V_sigma'
            ],
        data=[
            np.empty(N, dtype=object),
            np.empty(N, dtype=object),
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array
            ]
        )
    different_filters = table.unique(large_stars_table, keys='filter')
    different_filter_list = list(different_filters['filter'])
    different_filter_list = [different_filter.lower() for different_filter in different_filter_list]
    different_filter_data = np.empty((N, len(different_filter_list)))
    different_filter_data.fill(np.nan)
    filter_sigma_list = []
    for different_filter in different_filter_list:
        filter_sigma_list.append(f"{different_filter}_sigma")
    different_filter_table = Table(data=different_filter_data, names=different_filter_list)
    different_filter_sigma_table = Table(data=different_filter_data, names=filter_sigma_list)
    if ground_based:
        filter_X_list = []
        for different_filter in different_filter_list:
            filter_X_list.append(f"X_{different_filter}")
        filter_X_sigma_list = []
        for different_filter in different_filter_list:
            filter_X_sigma_list.append(f"X_{different_filter}_sigma")
        filter_X_table = Table(data=different_filter_data, names=filter_X_list)
        filter_X_sigma_table = Table(data=different_filter_data, names=filter_X_sigma_list)
        stars_table = table.hstack([apparent_mags_table,
                                    different_filter_table,
                                    different_filter_sigma_table,
                                    filter_X_table,
                                    filter_X_sigma_table],
                                   join_type='exact')
    else:
        stars_table = table.hstack([apparent_mags_table,
                                    different_filter_table,
                                    different_filter_sigma_table],
                                   join_type='exact')
    stars_table['Field'] = unique_stars['Field']
    stars_table['V'] = unique_stars['V']
    stars_table['(B-V)'] = unique_stars['(B-V)']
    stars_table['(U-B)'] = unique_stars['(U-B)']
    stars_table['(V-R)'] = unique_stars['(V-R)']
    stars_table['(V-I)'] = unique_stars['(V-I)']
    stars_table['V_sigma'] = unique_stars['V_sigma']
    i = 0
    for star in unique_stars['Name']:
        stars_table['Name'][i] = star
        mask = large_stars_table['Name'] == star
        current_star_table = large_stars_table[mask]
        unique_filters = table.unique(current_star_table, keys='filter')
        for unique_filter in unique_filters['filter']:
            mask = ((current_star_table['filter'] == unique_filter))
            current_star_filter_table = current_star_table[mask]
            mags_numpy = np.array(current_star_filter_table['mag_instrumental'])
            mean_mag = mags_numpy.mean()
            std_mag = mags_numpy.std()
            filter_column = unique_filter.lower()
            sigma_column = f'{filter_column}_sigma'
            stars_table[filter_column][i] = mean_mag
            stars_table[sigma_column][i] = std_mag
            if ground_based:
                X_numpy = np.array(current_star_filter_table['X'])
                mean_X = X_numpy.mean()
                std_X = X_numpy.std()
                X_column = f'X_{filter_column}'
                X_std_column = f'X_{filter_column}_sigma'
                stars_table[X_column][i] = mean_X
                stars_table[X_std_column][i] = std_X
        i += 1 
    return stars_table


def write_table_to_latex(table, output_file, formats=None):
    if not formats:
        ascii.write(table, output=output_file, format='latex')
    else:
        ascii.write(table, output=output_file, format='latex', formats=formats)


def space_based_transform(stars_table, 
                          plot_results=False, 
                          save_plots=False,
                          index='(B-V)', 
                          app_filter='V', 
                          instr_filter='b', 
                          field=None,
                          **kwargs):
    """
    Calculate the transforms for a space based sensor.

    Parameters
    ----------
    stars_table : astropy.table.table.Table
        Table containing the mean of the important information for each star. Has columns:
            Field : string
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
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
                Standard deviation of the apparent V magnitude from the reference file.
            <filter> : numpy.float64
                Mean instrumental magnitude of all detections of the star in <filter>. There is a different column for 
                each different filter used across the images.
            <filter>_sigma : numpy.float64
                Standard deviation of the instrumental magnitudes of all detections of the star in <filter>. 
                There is a different column for each different filter used across the images.
            X_<filter> : numpy.float64
                Mean airmass of all detections of the star in <filter>. There is a different column for each different 
                filter used across the images. Only output if ground_based is True.
            X_<filter>_sigma : numpy.float64
                Standard deviation of the airmasses of all detections of the star in <filter>. There is a different 
                column for each different filter used across the images. Only output if ground_based is True.
    plot_results : bool, optional
        Controls whether or not to plot the results from the transforms. The default is False.
    index : string, optional
        Colour index to calculate the transform for. The default is '(B-V)'.
    app_filter : string, optional
        Apparent magnitude filter band to calculate the transform for. The default is 'V'.
    instr_filter : string, optional
        Instrumental filter band to calculate the transform for. The default is 'clear'.
    field : string, optional
        Unique identifier of the star field that the reference star is in (e.g. Landolt field "108"). 
        The default is None.

    Returns
    -------
    filter_fci : float
        Instrumental transform coefficient for filter f using the colour index CI.
    zprime_fci : float
        Zero point magnitude for filter f.

    """
    max_app_filter_sigma = max(stars_table[f'{app_filter}_sigma'])
    max_instr_filter_sigma = max(stars_table[f'{instr_filter}_sigma'])
    err_sum = np.nan_to_num(stars_table[f'{app_filter}_sigma'], nan=max_app_filter_sigma) + \
        np.nan_to_num(stars_table[f'{instr_filter}_sigma'], nan=max_instr_filter_sigma)
    err_sum = np.array(err_sum)
    err_sum[err_sum == 0] = max(err_sum)
    a_fit, cov = curve_fit(linear_func, stars_table[index], 
                             stars_table[f'{app_filter}'] - stars_table[instr_filter], p0=None,
                             sigma=err_sum)
    filter_fci = a_fit[0]
    filter_fci_sigma = sqrt(cov[0][0])
    zprime_fci = a_fit[1]
    zprime_fci_sigma = sqrt(cov[1][1])
    if plot_results:
        index_plot = np.arange(start=min(stars_table[index])-0.2, stop=max(stars_table[index])+0.2, step=0.1)
        plt.errorbar(stars_table[index], stars_table[f'{app_filter}'] - stars_table[instr_filter], 
                     yerr=err_sum, fmt='o', capsize=2)
        plt.plot(index_plot, filter_fci * index_plot + zprime_fci, 
                 label=f"({app_filter}-{instr_filter}) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        plt.ylabel(f"{app_filter}-{instr_filter}")
        plt.xlabel(f"{index}")
        plt.legend()
        plt.title(f"Space Based Transform Calculation for {index}")
        if save_plots:
            save_loc = f"{os.path.join(kwargs.get('save_loc'), f'TZfci_{index}')}.png"
            plt.savefig(save_loc)
        # if not field:
        #     plt.title(f"({app_filter}-{instr_filter}) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        # else:
        #     plt.title(f"{field}: ({app_filter}-{instr_filter}) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        plt.show()
        plt.close()
    return filter_fci, filter_fci_sigma, zprime_fci, zprime_fci_sigma


def get_avg_airmass(altazpositions):
    """
    Get the average airmass for all sources detected in the image.

    Parameters
    ----------
    altazpositions : astropy.coordinates.sky_coordinate.SkyCoord
        Alt/Az positions of all of the sources in the image.

    Returns
    -------
    avg_airmass : float
        Mean airmass of all of the sources in the image.

    """
    airmasses = np.array(altazpositions.secz.value)
    avg_airmass = airmasses.mean()
    return avg_airmass





def BackgroundEstimationMulti(fitsdata, sigma_clip, bkgmethod, printval):
   
    sigma_clip = SigmaClip(sigma=2.5)
    bkg = SExtractorBackground(sigma_clip)
    bkg_value1 = bkg.calc_background(fitsdata)
    
    #print(bkg_value)
    bkg = MeanBackground(sigma_clip)
    bkg_value2 = bkg.calc_background(fitsdata)
    
    bkg = MedianBackground(sigma_clip)
    bkg_value3 = bkg.calc_background(fitsdata)
    
    #bkg = ModeEstimatorBackground(sigma_clip)
   #bkg_value4 = bkg.calc_background(fitsdata)
    
    
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
        #print("Mode Estimator Background: " + str(bkg_value4))
        print("Remaining Background (subtracted): " + str(bg_rem))
        print("Polyfit Background: Not Implemented Yet")
        
    else:
        return bkg
    
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
            starx= mstar.X
            stary=mstar.Y
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

def pinpoint_init():
    f = win32com.client.Dispatch("Pinpoint.plate")
    return f

def getFileList(inbox):
    filepathall = []
    directory = os.path.dirname(inbox)
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
        focal_Length=imagehdularray[0].header['FOCALLEN']
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
   y_arcsecperpixel = math.atan(xpix_size/focal_Length)*3600*ybin
   
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

inbox = 'D:\\Wawrow\\2. Observational Data\\2021-03-10 - Calibrated\\HIP 46066\\LIGHT\\B\\' #Image Location of .fits Format
refstars_doc = 'D:\\Reference_stars.xlsx'
refstars_csv='D:\\Reference_stars.csv' #Reference Star List
catloc = 'D:\squid\\USNOA20-All'; #Catalog #TODO Change Catalog to UCAC3 
save_loc = os.path.join(inbox, 'Outputs') # Output Folder for Files

"#TODO Initialize these once astroreduction is included"
bias = 0 
darks = 0 
flats = 0
OutputsaveLoc = 0 ; #0 Default will save outputs in image folder

"Read Ref Doc"

HIP, erad, edec, vref, bvindex, vrindex, refstarsfin = ref_star_folder_read(refstars_doc)
reference_stars, ref_star_positions = astro.read_ref_stars(refstars_csv) # Reading the Reference Star Doc
f = pinpoint_init() #Start Pinpoint 

"Set Variables"
streak_array= [] #Streak Detection for Track Rate Mode
sigma_clip = 3.5; #Clipping Factor
edge_protect = 10; #Img Edge Clipping
min_obj_pixels = 5 #Min Pixels to qualify as a Point Source
SNRLimit = 0; #Signal-To-Noise Ratio
ground_based = False # TODO Add Ground Based Selection Input
pinpoint = False #TODO Add Pinpoint Solve or Not to GUI

"Opening Image Folder and Determing the number of files"
filepathall = getFileList(inbox); #Get List of Images
large_table_columns= init_large_table_columns() #Create Table Template for Space-Based Detected Star Transform storage
gb_transform_table_columns = astro.init_gb_transform_table_columns() #Create Table Template for Detected Star Transform storage
large_stars_table = create_large_stars_table(large_table_columns, ground_based=False) #Create Table Template for Space-Based Detected Star Transform storage



for i in range(1,len(filepathall)):
    if filepathall[i].endswith(".fits"):
        #f = win32com.client.Dispatch("Pinpoint.plate")
        print("Processing Image: " + filepathall[i])
        
        "Import Data from FITS Image"
        imagehdularray, date, exposure_Time,imagesizeX, imagesizeY, fitsdata, filt, header,XPIXSZ, YPIXSZ,wcs = fits_header_import(filepathall[i])
    
        """Pinpoint Solve"""
        if pinpoint:
            try:
                f.AttachFITS(filepathall[i])
                f.Declination = f.targetDeclination;
                f.RightAscension = f.targetRightAscension; 
                x_arcsecperpixel, y_arcsecperpixel = calc_ArcsecPerPixel(header)
                # yBin = 4.33562092816E-004*3600;
                # xBin =  4.33131246330E-004*3600; 
                f.ArcsecperPixelHoriz  =  x_arcsecperpixel;
                f.ArcsecperPixelVert =  y_arcsecperpixel;
                
                
                "Pinpoint Solve Inputs"
                #TODO Add Inputs for pinpoint solving to GUI
                f.Catalog = 5;
                f.CatalogPath = catloc;
                f.CatalogMaximumMagnitude = 13;
                f.CatalogExpansion = 0.8;
                f.SigmaAboveMean = 2.5; 
                f.FindImageStars; 
                f.FindCatalogStars; 
                f.MaxSolveTime = 60; 
                f.MaxMatchResidual = 1.5; 
        
                
                "Pinpoint Solving"
                f.FindCatalogStars()
                f.Solve()
                #f.MatchedStars.count
                #f.FindImageStars()
                #print(f.ImageStars)
                
                "Searching for Ref Stars"
                #TODO Fix ref star search
                #ref_star_search(s,f,erad,edec, HIP, vref,bvindex,vrindex,refstarsfin)
                f.DetachFITS()
                f=None
                print ("Pinpoint Solved")
            except:
                    continue

            
        "Import Data from FITS Image"
        
    
        """
        1. Space Based - Airmass not a factor in determining transforms
        2. Ground Based - Multiple Order Transforms with Airmass and Extinctions
    
        """
        
        """Running Functions"""
        bkg = BackgroundEstimationMulti(fitsdata, 2.5, 1, 0) #Background Estimation
        backg, bkg_std = calculate_img_bkg(fitsdata) # Find Median Bkg and bkg_std
        iraf_Sources= detecting_stars(fitsdata, backg, bkg_std) # Solve Image using IRAF Starfinder
        if not iraf_Sources:
                continue
            
        fwhm, fwhm_stdev= calculate_fwhm(iraf_Sources) #Calculate Full-Width Half-Max
        photometry_result= perform_photometry(iraf_Sources, fwhm, fitsdata, backg) #Determine Photometry
        calculate_magnitudes(photometry_result, exposure_Time) #Determine Magnitudes
        calculate_magnitudes_sigma(photometry_result, exposure_Time) 
        fluxes = np.array(photometry_result['flux_fit']) #
        
        skypositions= convert_pixel_to_ra_dec(iraf_Sources,wcs) #Calculate Sky Position in RA and DEC
        altazpositions = None
        
        if ground_based:
                try:
                    altazpositions = astro.convert_ra_dec_to_alt_az(skypositions, header, lat_key='OBSGEO-B', 
                                                                    lon_key='OBSGEO-L', elev_key='OBSGEO-H')
                   
                except AttributeError as e:
                    print(e)
                    continue
            
        
        instr_mags = calculate_magnitudes(photometry_result, exposure_Time)
        instr_mags_sigma = calculate_magnitudes_sigma(photometry_result, exposure_Time)
        matched_stars = astro.find_ref_stars(reference_stars, 
                                                 ref_star_positions,
                                                 skypositions,
                                                 instr_mags,
                                                 instr_mags_sigma,
                                                 fluxes,
                                                 ground_based=ground_based,
                                                 altazpositions=False)
        
        ref_star = reference_stars[matched_stars.ref_star_index]
        vref=ref_star['V']
        
        #print (matched_stars)
        
        if not matched_stars:
                continue
        
        instr_filter = astro.get_instr_filter_name(header)
        colour_indices = astro.get_all_colour_indices(instr_filter)
        
        
        if ground_based:
            for colour_index in colour_indices:
                    # _, _, _, colour_index = astro.get_app_mag_and_index(matched_stars.ref_star, instr_filter)
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
                        # print(len(matched_stars.img_instr_mag))
                    except TypeError:
                        print("Only 1 reference star detected in the image.")
                        continue
                    avg_airmass = astro.get_avg_airmass(altazpositions)
                            # if avg_airmass > 2.0:
                            #     continue
                    c_fci, c_fci_sigma, zprime_f, zprime_f_sigma = astro.ground_based_first_order_transforms(matched_stars, 
                                                                                                                      instr_filter, 
                                                                                                                      colour_index, 
                                                                                                                      plot_results=True,
                                                                                                                      save_plots=True,
                                                                                                                      save_loc=save_loc,
                                                                                                                      unique_id=filepathall[i])
                    gb_transform_table_columns = astro.update_gb_transform_table_columns(gb_transform_table_columns,
                                                                                                  field,
                                                                                                  c_fci,
                                                                                                  c_fci_sigma,
                                                                                                  zprime_f,
                                                                                                  zprime_f_sigma,
                                                                                                  instr_filter,
                                                                                                  colour_index,
                                                                                                  altazpositions)
        else:   
             large_table_columns= update_large_table_columns(large_table_columns, matched_stars, header, exposure_Time, ref_star,  ground_based=False, name_key='Name')
                #large_table_columns = astro.update_large_table_columns(large_table_columns, 
                                                                              # matched_stars, 
                                                                              # header, 
                                                                              # exposure_Time, 
                                                                              # ground_based=ground_based, 
                                                                              # name_key='Name')

if ground_based:
    gb_transform_table = astro.create_gb_transform_table(gb_transform_table_columns)
    gb_transform_table.pprint_all()
    gb_transform_table = astro.remove_large_airmass(gb_transform_table)
    gb_transform_table.pprint_all()
    formats = {
        'C_fCI': '%0.3f',
        'C_fCI_sigma': '%0.3f',
        'Zprime_f': '%0.3f',
        'Zprime_f_sigma': '%0.3f',
        'X': '%0.3f'
        }
    astro.write_table_to_latex(gb_transform_table, f"{os.path.join(save_loc, 'gb_transform_table')}.txt", formats=formats)
    gb_final_transforms = astro.ground_based_second_order_transforms(gb_transform_table, 
                                                                      plot_results=True, save_plots=True, save_loc=save_loc)
    gb_final_transforms.pprint_all()
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
    astro.write_table_to_latex(gb_final_transforms, f"{os.path.join(save_loc, 'gb_final_transforms')}.txt", formats=formats)
            
            
        
else:
    large_stars_table = create_large_stars_table(large_table_columns, ground_based=False)
    stars_table= group_each_star(large_stars_table, vref, ground_based=False, keys='Name')
    #print(stars_table)
    transform_index_list = ['(B-V)', '(V-R)', '(V-I)']
    astro.write_table_to_latex(gb_transform_table, f"{os.path.join(save_loc, 'gb_transform_table')}.txt", formats=formats)
    
    for index in transform_index_list:
        filter_fci, zprime_fci = space_based_transform(stars_table, plot_results=True, index=index)
        print(f"(V-clear) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
            

    

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



n=0;
flag = 1;

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