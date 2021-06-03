# -*- coding: utf-8 -*-
"""
AstroFunctions.py.

This file holds all of the functions that will likely be implemented in a final version of the image processor.

Created on Thu Apr 15 10:14:43 2021

@author: Jack Wawrow
"""
from astropy import table
from astropy.coordinates import EarthLocation, AltAz, SkyCoord, match_coordinates_sky
from astropy.io import fits, ascii
from astropy.modeling.fitting import LevMarLSQFitter, LinearLSQFitter, FittingWithOutlierRemoval
from astropy.modeling.models import Linear1D
from astropy.stats import sigma_clip, sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.table import Table, QTable, hstack
from astropy.time import Time
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


def linear_func(x, m, b):
    y = (m * x) + b
    return y


def init_linear_fitting(niter=3, sigma=3.0):
    fit = LevMarLSQFitter(calc_uncertainties=True)
    or_fit = FittingWithOutlierRemoval(fit, sigma_clip, niter=niter, sigma=sigma)
    line_init = Linear1D()
    return fit, or_fit, line_init


def read_ref_stars(ref_stars_file):
    """
    Read a file containing information regarding the reference stars to be used to calculate the transforms.
    
    Parameters
    ----------
    ref_stars_file : string
        Location of the reference stars file.
            
    Returns
    -------
    reference_stars : astropy.table.Table
        Table with the data extracted from ref_stars_file.
    ref_star_positions : astropy.coordinates.sky_coordinate.SkyCoord
        AstroPy SkyCoord object containing the RA/dec positions of all reference stars in the file.
    
    """
    try:
        reference_stars = ascii.read(ref_stars_file, format='basic', delimiter='\t', guess=False, encoding='UTF-8')
    except Exception:
        reference_stars = ascii.read(ref_stars_file, encoding='UTF-8')
    reference_stars = reference_stars.filled(np.nan)
    ref_star_positions = SkyCoord(ra=reference_stars['RA'], dec=reference_stars['Dec'], unit=(u.hourangle, u.deg))
    return reference_stars, ref_star_positions


def read_fits_file(filepath):
    """
    Read a fits file and return its header and data.

    Parameters
    ----------
    filepath : string
        Location of the fits file.

    Returns
    -------
    hdr : astropy.io.fits.header.Header
        Header from the fits file.
    imgdata : numpy.ndarray
        Data from the fits file.

    """
    with fits.open(filepath) as hdul:
        hdr = hdul[0].header
        imgdata = hdul[0].data
    return hdr, imgdata

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
    instr_filter = hdr[filter_key][0].lower()
    return instr_filter


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


def detecting_stars(imgdata, bkg, bkg_std, fwhm=2.0, sigma=4.0):
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
    sigma : float
        Number of standard deviations above the background to use as the threshold.

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
    # iraffind = IRAFStarFinder(threshold=bkg+3*bkg_std, fwhm=fwhm)
    iraffind = IRAFStarFinder(threshold=sigma*bkg_std, fwhm=fwhm)
    irafsources = iraffind(imgdata - bkg)
    return irafsources


def convert_pixel_to_ra_dec(irafsources, wcs):
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
    TODO
    fwhm : float
        Mean FWHM of all sources in the image.
    fwhm_std : float
        Standard deviation of the FWHM of all sources in the image.

    """
    fwhms = np.array(irafsources['fwhm'])
    fwhm = fwhms.mean()
    fwhm_std = fwhms.std()
    return fwhms, fwhm, fwhm_std


def convert_fwhm_to_arcsec(hdr, fwhms, fwhm, fwhm_std, 
                           focal_length_key='FOCALLEN', xpixsz_key='XPIXSZ', ypixsz_key='YPIXSZ'):
    try:
        focal_length = hdr[focal_length_key] * u.mm
        xpixsz = hdr[xpixsz_key]
        ypixsz = hdr[ypixsz_key]
        if xpixsz == ypixsz:
            pixsz = xpixsz * u.um
            rad_per_pix = atan(pixsz / focal_length) * u.rad
            arcsec_per_pix = rad_per_pix.to(u.arcsec)
            fwhm_arcsec = fwhm * arcsec_per_pix.value
            fwhm_std_arcsec = fwhm_std * arcsec_per_pix.value
            fwhms_arcsec = fwhms * arcsec_per_pix.value
    except KeyError:
        fwhms_arcsec = fwhms
    return fwhms_arcsec, fwhm_arcsec, fwhm_std_arcsec


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
        Expsure time of the image in seconds.

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


def update_large_table_columns(large_table_columns, matched_stars, hdr, exptime, ground_based=False, name_key='Name'):
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
    updated_large_table_columns = large_table_columns
    try:
        num_stars = len(matched_stars.img_instr_mag)
    except TypeError:
        num_stars = 1
    if num_stars > 1:
        for row in matched_stars.ref_star:
            split_string = re.split('[^a-zA-Z0-9]', str(row[name_key]))
            if len(split_string) > 1:
                updated_large_table_columns.field.append(' '.join(split_string[:-1]))
            elif len(split_string) == 1:
                updated_large_table_columns.field.append(split_string[0])
            else:
                print('Could not find the name of the field.')
        updated_large_table_columns.ref_star_name.extend(matched_stars.ref_star[name_key])
        updated_large_table_columns.flux_table.extend(matched_stars.flux)
        time = Time(hdr['DATE-OBS'], format='fits')
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
        filter_name_repeat = np.full(num_stars, hdr['FILTER'][0])
        updated_large_table_columns.filters.extend(filter_name_repeat)
        updated_large_table_columns.V_apparents.extend(matched_stars.ref_star['V_ref'])
        try:
            updated_large_table_columns.B_V_apparents.extend(matched_stars.ref_star['(B-V)'])
            updated_large_table_columns.U_B_apparents.extend(matched_stars.ref_star['(U-B)'])
            updated_large_table_columns.V_R_apparents.extend(matched_stars.ref_star['(V-R)'])
            updated_large_table_columns.V_I_apparents.extend(matched_stars.ref_star['(V-I)'])
            updated_large_table_columns.V_sigma_apparents.extend(matched_stars.ref_star['V_sigma'])
        except KeyError:
            updated_large_table_columns.B_V_apparents.extend(matched_stars.ref_star['B-V'])
            updated_large_table_columns.U_B_apparents.extend(matched_stars.ref_star['U-B'])
            updated_large_table_columns.V_R_apparents.extend(matched_stars.ref_star['V-R'])
            updated_large_table_columns.V_I_apparents.extend(matched_stars.ref_star['V-I'])
            updated_large_table_columns.V_sigma_apparents.extend(matched_stars.ref_star['e_V'])
        if not ground_based:
            return updated_large_table_columns
        updated_large_table_columns.img_star_airmass.extend(matched_stars.img_star_airmass)
        # updated_large_table_columns.X_rounded.extend(round(matched_stars.img_star_airmass, 1))
    elif num_stars == 1:
        split_string = re.split('[^a-zA-Z0-9]', str(matched_stars.ref_star[name_key]))
        if len(split_string) > 1:
            updated_large_table_columns.field.append(' '.join(split_string[:-1]))
        elif len(split_string) == 1:
            updated_large_table_columns.field.append(split_string[0])
        else:
            print('Could not find the name of the field.')
        updated_large_table_columns.ref_star_name.append(matched_stars.ref_star[name_key])
        updated_large_table_columns.flux_table.append(matched_stars.flux)
        time = Time(hdr['DATE-OBS'], format='fits')
        updated_large_table_columns.times.append(time.jd)
        updated_large_table_columns.exposure.append(exptime)
        updated_large_table_columns.ref_star_RA.append(matched_stars.ref_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.ref_star_dec.append(matched_stars.ref_star_loc.dec)
        updated_large_table_columns.img_star_RA.append(matched_stars.img_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.img_star_dec.append(matched_stars.img_star_loc.dec)
        updated_large_table_columns.angular_separation.append(matched_stars.ang_separation.to(u.arcsec))
        updated_large_table_columns.img_star_mag.append(matched_stars.img_instr_mag)
        updated_large_table_columns.img_star_mag_sigma.append(matched_stars.img_instr_mag_sigma)
        updated_large_table_columns.filters.append(hdr['FILTER'][0])
        updated_large_table_columns.V_apparents.append(matched_stars.ref_star['V_ref'])
        try:
            updated_large_table_columns.B_V_apparents.append(matched_stars.ref_star['(B-V)'])
            updated_large_table_columns.U_B_apparents.append(matched_stars.ref_star['(U-B)'])
            updated_large_table_columns.V_R_apparents.append(matched_stars.ref_star['(V-R)'])
            updated_large_table_columns.V_I_apparents.append(matched_stars.ref_star['(V-I)'])
            updated_large_table_columns.V_sigma_apparents.append(matched_stars.ref_star['V_sigma'])
        except KeyError:
            updated_large_table_columns.B_V_apparents.append(matched_stars.ref_star['B-V'])
            updated_large_table_columns.U_B_apparents.append(matched_stars.ref_star['U-B'])
            updated_large_table_columns.V_R_apparents.append(matched_stars.ref_star['V-R'])
            updated_large_table_columns.V_I_apparents.append(matched_stars.ref_star['V-I'])
            updated_large_table_columns.V_sigma_apparents.append(matched_stars.ref_star['e_V'])
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
                large_table_columns.img_star_airmass, 
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
                'V_ref',
                '(B-V)',
                '(U-B)',
                '(V-R)',
                '(V-I)',
                'V_sigma',
                'X',
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
                'V_ref',
                '(B-V)',
                '(U-B)',
                '(V-R)',
                '(V-I)',
                'V_sigma',
                ]
            )
    return large_stars_table


def group_each_star(large_stars_table, ground_based=False, keys='Name'):
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
            'V_ref',
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
    stars_table['V_ref'] = unique_stars['V_ref']
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
                          instr_filter='c', 
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
    
    
    x = stars_table[index]
    y = stars_table[f'{app_filter}_ref'] - stars_table[instr_filter]
    fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
    fitted_line, mask = or_fit(line_init, x, y, weights=1.0/err_sum)
    filtered_data = np.ma.masked_array(y, mask=mask)
    filter_fci = fitted_line.slope.value
    zprime_fci = fitted_line.intercept.value
    cov = fit.fit_info['param_cov']
    
    # a_fit, cov = curve_fit(linear_func, stars_table[index], 
    #                          stars_table[f'{app_filter}_ref'] - stars_table[instr_filter], 
    #                          sigma=err_sum)
    # filter_fci = a_fit[0]
    filter_fci_sigma = sqrt(cov[0][0])
    # zprime_fci = a_fit[1]
    zprime_fci_sigma = sqrt(cov[1][1])
    if plot_results:
        index_plot = np.arange(start=min(stars_table[index]), stop=max(stars_table[index])+0.01, step=0.01)
        plt.errorbar(x, y, yerr=err_sum, color='#1f77b4', fmt='o', fillstyle='none', capsize=2, label="Clipped Data")
        plt.plot(x, filtered_data, 'o', color='#1f77b4', label="Fitted Data")
        plt.plot(index_plot, fitted_line(index_plot), '-', color='#ff7f0e', 
                 label=f"({app_filter}-{instr_filter}) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        
        # plt.errorbar(stars_table[index], stars_table[f'{app_filter}_ref'] - stars_table[instr_filter], 
        #              yerr=err_sum, fmt='o', capsize=2)
        # plt.plot(index_plot, filter_fci * index_plot + zprime_fci, 
        #          label=f"({app_filter}-{instr_filter}) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        plt.ylabel(f"{app_filter}-{instr_filter}")
        plt.xlabel(f"{index}")
        plt.legend()
        plt.title(f"Space Based Transform Calculation for {index}")
        if save_plots:
            unique_id = kwargs.get('unique_id')
            save_loc = f"{os.path.join(kwargs.get('save_loc'), f'{unique_id}_TZfci_{index}')}.png"
            plt.savefig(save_loc)
        # if not field:
        #     plt.title(f"({app_filter}-{instr_filter}) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        # else:
        #     plt.title(f"{field}: ({app_filter}-{instr_filter}) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        plt.show()
        plt.close()
    return filter_fci, filter_fci_sigma, zprime_fci, zprime_fci_sigma


def init_sb_final_transform_columns():
    index = []
    filter_fci = []
    filter_fci_sigma = []
    zprime_fci = [] 
    zprime_fci_sigma = []
    sb_final_transform_columns = namedtuple('sb_final_transform_columns', 
                                            ['index',
                                             'filter_fci',
                                             'filter_fci_sigma',
                                             'zprime_fci',
                                             'zprime_fci_sigma'])
    return sb_final_transform_columns(index, 
                                      filter_fci, 
                                      filter_fci_sigma, 
                                      zprime_fci, 
                                      zprime_fci_sigma)


def update_sb_final_transform_columns(sb_final_transform_columns,
                                      index,
                                      filter_fci, 
                                      filter_fci_sigma, 
                                      zprime_fci, 
                                      zprime_fci_sigma):
    updated_sb_final_transform_columns = sb_final_transform_columns
    updated_sb_final_transform_columns.index.append(index)
    updated_sb_final_transform_columns.filter_fci.append(filter_fci)
    updated_sb_final_transform_columns.filter_fci_sigma.append(filter_fci_sigma)
    updated_sb_final_transform_columns.zprime_fci.append(zprime_fci)
    updated_sb_final_transform_columns.zprime_fci_sigma.append(zprime_fci_sigma)
    return updated_sb_final_transform_columns


def create_sb_final_transform_table(sb_final_transform_columns):
    sb_final_transform_table = Table(
        names=[
            'CI',
            'T_fCI',
            'T_fCI_sigma',
            'Z_fCI',
            'Z_fCI_sigma'
            ],
        data=[
            sb_final_transform_columns.index,
            sb_final_transform_columns.filter_fci,
            sb_final_transform_columns.filter_fci_sigma,
            sb_final_transform_columns.zprime_fci,
            sb_final_transform_columns.zprime_fci_sigma
            ]
        )
    return sb_final_transform_table

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


def init_gb_transform_table_columns():
    """
    Initialize the columns that will create the table used for the ground-based transforms.

    Returns
    -------
    gb_transform_table_columns : namedtuple
        Attributes:
            field : empty list
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            c_fci : empty list
                C coefficient for filter f with colour index ci.
                TODO
            zprime_f : empty list
                Z' coefficient for filter f.
            instr_filter : empty list
                Instrumental filter band to calculate the transform for.
            colour_index : empty list
                Name of the colour index used to calculate c_fci and zprime_f.
            airmass : empty list
                The mean airmass for all sources in the image.

    """
    field = []
    c_fci = []
    # TODO: Calculate c_fci_simga and zprime_f_sigma.
    c_fci_simga = []
    zprime_f = []
    zprime_f_sigma = []
    instr_filter = []
    colour_index = []
    airmass = []
    gb_transform_table_columns = namedtuple('gb_transform_table_columns', 
                                            ['field',
                                             'c_fci',
                                             'c_fci_sigma',
                                             'zprime_f',
                                             'zprime_f_sigma',
                                             'instr_filter',
                                             'colour_index',
                                             'airmass'])
    return gb_transform_table_columns(field, 
                                      c_fci, 
                                      c_fci_simga,
                                      zprime_f, 
                                      zprime_f_sigma, 
                                      instr_filter, 
                                      colour_index, 
                                      airmass)


def update_gb_transform_table_columns(gb_transform_table_columns,
                                      field,
                                      c_fci,
                                      c_fci_sigma,
                                      zprime_f,
                                      zprime_f_sigma,
                                      instr_filter,
                                      colour_index,
                                      altazpositions):
    """
    Update columns to be used for the transform table based on information from the current image.

    Parameters
    ----------
    gb_transform_table_columns : namedtuple
        Attributes:
            field : string list
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            c_fci : np.float64
                C coefficient for filter f with colour index ci.
                TODO
            zprime_f : numpy.float64
                Z' coefficient for filter f.
            instr_filter : string list
                Instrumental filter band to calculate the transform for.
            colour_index : string list
                Name of the colour index used to calculate c_fci and zprime_f.
            airmass : numpy.float64
                The mean airmass for all sources in the image.
            
    field : string
        Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
    c_fci : float
        C coefficient for filter f with colour index ci.
    zprime_f : float
        Z' coefficient for filter f.
        TODO
    instr_filter : string
        Instrumental filter band to calculate the transform for.
    colour_index : string
        Name of the colour index used to calculate c_fci and zprime_f.
    altazpositions : astropy.coordinates.sky_coordinate.SkyCoord
        Alt/Az positions of all of the sources in the image.

    Returns
    -------
    updated_gb_transform_table_columns : namedtuple
        Attributes:
            field : string list
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            c_fci : np.float64
                C coefficient for filter f with colour index ci.
            zprime_f : numpy.float64
                Z' coefficient for filter f.
            instr_filter : string list
                Instrumental filter band to calculate the transform for.
            colour_index : string list
                Name of the colour index used to calculate c_fci and zprime_f.
            airmass : numpy.float64
                The mean airmass for all sources in the image.

    """
    updated_gb_transform_table_columns = gb_transform_table_columns
    updated_gb_transform_table_columns.field.append(field)
    updated_gb_transform_table_columns.c_fci.append(c_fci)
    updated_gb_transform_table_columns.c_fci_sigma.append(c_fci_sigma)
    updated_gb_transform_table_columns.zprime_f.append(zprime_f)
    updated_gb_transform_table_columns.zprime_f_sigma.append(zprime_f_sigma)
    updated_gb_transform_table_columns.instr_filter.append(instr_filter)
    updated_gb_transform_table_columns.colour_index.append(colour_index)
    avg_airmass = get_avg_airmass(altazpositions)
    updated_gb_transform_table_columns.airmass.append(avg_airmass)
    return updated_gb_transform_table_columns


def create_gb_transform_table(gb_transform_table_columns):
    """
    Convert the columns of the ground-based transform table into an AstroPy table.

    Parameters
    ----------
    gb_transform_table_columns : namedtuple
        Attributes:
            field : string list
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            c_fci : np.float64
                C coefficient for filter f with colour index ci.
                TODO
            zprime_f : numpy.float64
                Z' coefficient for filter f.
            instr_filter : string list
                Instrumental filter band to calculate the transform for.
            colour_index : string list
                Name of the colour index used to calculate c_fci and zprime_f.
            airmass : numpy.float64
                The mean airmass for all sources in the image.

    Returns
    -------
    gb_transform_table : astropy.table.Table
        Table containing the results from the intermediary step in calculating the transforms. Has columns:
            field : string list
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            C_fCI : np.float64
                C coefficient for filter f with colour index ci.
                TODO
            Zprime_f : numpy.float64
                Z' coefficient for filter f.
            filter : string list
                Instrumental filter band to calculate the transform for.
            CI : string list
                Name of the colour index used to calculate c_fci and zprime_f.
            X : numpy.float64
                The mean airmass for all sources in the image.

    """
    gb_transform_table = Table(
        names=[
            'Field',
            'C_fCI',
            'C_fCI_sigma',
            'Zprime_f',
            'Zprime_f_sigma',
            'filter',
            'CI',
            'X'
            ],
        data=[
            gb_transform_table_columns.field,
            gb_transform_table_columns.c_fci,
            gb_transform_table_columns.c_fci_sigma,
            gb_transform_table_columns.zprime_f,
            gb_transform_table_columns.zprime_f_sigma,
            gb_transform_table_columns.instr_filter,
            gb_transform_table_columns.colour_index,
            gb_transform_table_columns.airmass
            ]
        )
    return gb_transform_table


def get_app_mag_and_index(ref_star, instr_filter):
    """
    Get the desired apparent magnitude and colour filter for the particular instr_filter.
    
    Parameters
    ----------
    ref_star : astropy.table.Table
        Rows from ref_stars_file that correspond to a matched reference star.
    instr_filter : string
        Instrumental filter band to calculate the transform for.

    Returns
    -------
    app_mag : numpy.float64
        Apparent magnitude of the reference star(s) in the desired filter (e.g. B for b, V for v, etc.).
    app_mag_sigma : numpy.float64
        Standard deviation of the apparent magnitude of the reference star(s) in the same filter as app_mag.
    colour_index : string
        Name of the colour index to use when calculating the transform (e.g. B-V for b, V-R for r).

    """
    if instr_filter == 'b':
        colour_index = 'B-V'
        app_filter = 'B'
        app_mag = np.array(ref_star['V_ref'] + ref_star[colour_index])
        app_mag_sigma = np.nan_to_num(ref_star['e_V'], nan=max(ref_star['e_V'])) + \
            np.nan_to_num(ref_star[f'e_{colour_index}'], nan=max(ref_star[f'e_{colour_index}']))
        # app_mag_sigma = np.array(ref_star['e_V'] + ref_star[f'e_{colour_index}'])
    elif instr_filter == 'v' or instr_filter == 'g':
        colour_index = 'B-V'
        app_filter = 'V'
        app_mag = np.array(ref_star['V_ref'])
        app_mag_sigma = np.nan_to_num(ref_star['e_V'], nan=max(ref_star['e_V']))
        # app_mag_sigma = np.array(ref_star['e_V'])
    elif instr_filter == 'r':
        colour_index = 'V-R'
        app_filter = 'R'
        app_mag = np.array(ref_star['V_ref'] - ref_star[colour_index])
        app_mag_sigma = np.nan_to_num(ref_star['e_V'], nan=max(ref_star['e_V'])) + \
            np.nan_to_num(ref_star[f'e_{colour_index}'], nan=max(ref_star[f'e_{colour_index}']))
        # app_mag_sigma = np.array(ref_star['e_V'] + ref_star[f'e_{colour_index}'])
    elif instr_filter == 'i':
        colour_index = 'V-I'
        app_filter = 'I'
        app_mag = np.array(ref_star['V_ref'] - ref_star[colour_index])
        app_mag_sigma = np.nan_to_num(ref_star['e_V'], nan=max(ref_star['e_V'])) + \
            np.nan_to_num(ref_star[f'e_{colour_index}'], nan=max(ref_star[f'e_{colour_index}']))
        # app_mag_sigma = np.array(ref_star['e_V'] + ref_star[f'e_{colour_index}'])
    else:
        colour_index = None
        app_filter = None
        app_mag = None
        app_mag_sigma = None
    return app_mag, app_mag_sigma, app_filter, colour_index


def ground_based_first_order_transforms(matched_stars, instr_filter, colour_index, field=None, 
                                        plot_results=False, save_plots=False, **kwargs):
    """
    Perform the intermediary step to calculating the ground based transforms.

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
                Alt/Az of the matched star(s) detected in the image.
            img_star_airmass : float
                sec(z) of img_star_altaz.
    instr_filter : string
        Instrumental filter band to calculate the transform for.
    colour_index : string
        TODO
    plot_results : bool, optional
        Controls whether or not to plot the results from the transforms. The default is False.
    field : string, optional
        Unique identifier of the star field that the reference star is in (e.g. Landolt field "108"). 
        The default is None.

    Returns
    -------
    c_fci : float
        C coefficient for filter f with colour index ci.
    zprime_fci : TYPE
        Z' coefficient for filter f.

    """
    try:
        len(matched_stars.img_instr_mag)
    except TypeError:
        return
    app_mag, app_mag_sigma, app_filter, _ = get_app_mag_and_index(matched_stars.ref_star, instr_filter)
    max_instr_filter_sigma = max(matched_stars.img_instr_mag_sigma)
    err_sum = app_mag_sigma + np.nan_to_num(matched_stars.img_instr_mag_sigma, nan=max_instr_filter_sigma)
    err_sum = np.array(err_sum)
    err_sum[err_sum == 0] = max(err_sum)
    x = matched_stars.ref_star[colour_index]
    y = app_mag - matched_stars.img_instr_mag
    fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
    fitted_line, mask = or_fit(line_init, x, y, weights=1.0/err_sum)
    filtered_data = np.ma.masked_array(y, mask=mask)
    c_fci = fitted_line.slope.value
    zprime_f = fitted_line.intercept.value
    cov = fit.fit_info['param_cov']
    if cov is None:
        c_fci_sigma = 0.0
        zprime_f_sigma = 0.0
    else:
        c_fci_sigma = sqrt(cov[0][0])
        zprime_f_sigma = sqrt(cov[1][1])
    if plot_results:
        index_plot = np.arange(start=min(matched_stars.ref_star[colour_index]), 
                               stop=max(matched_stars.ref_star[colour_index])+0.01, 
                               step=0.01)
        plt.errorbar(x, y, yerr=err_sum, color='#1f77b4', fmt='o', fillstyle='none', capsize=2, label="Clipped Data")
        plt.plot(x, filtered_data, 'o', color='#1f77b4', label="Fitted Data")
        plt.plot(index_plot, fitted_line(index_plot), '-', color='#ff7f0e', 
                 label=f"({app_filter}-{instr_filter}) = {c_fci:.3f} * {colour_index} + {zprime_f:.3f}")
        # plt.plot(index_plot, c_fci * index_plot + zprime_f, 
        #          label=f"({app_filter}-{instr_filter}) = {c_fci:.3f} * {colour_index} + {zprime_f:.3f}")
        plt.ylabel(f"{app_filter}-{instr_filter}")
        plt.xlabel(f"{colour_index}")
        plt.legend()
        plt.title(f"C and Z' Coefficient Calculations for {colour_index}")
        # if not field:
        #     plt.title(f"({app_filter}-{instr_filter}) = {c_fci:.3f} * {colour_index} + {zprime_f:.3f}")
        # else:
        #     plt.title(f"{field}: ({app_filter}-{instr_filter}) = {c_fci:.3f} * {colour_index} + {zprime_f:.3f}")
        if save_plots:
            unique_id = kwargs.get('unique_id')
            save_loc = f"{os.path.join(kwargs.get('save_loc'), f'CZprime{colour_index}_{unique_id}')}.png"
            plt.savefig(save_loc)
        plt.show()
        plt.close()
    return c_fci, c_fci_sigma, zprime_f, zprime_f_sigma


def get_all_colour_indices(instr_filter):
    if instr_filter == 'b':
        colour_indices = ['B-V']
    elif instr_filter == 'v' or instr_filter == 'g':
        colour_indices = ['B-V', 'V-R', 'V-I']
    elif instr_filter == 'r':
        colour_indices = ['V-R']
    elif instr_filter == 'i':
        colour_indices = ['V-I']
    else:
        colour_indices = None
    return colour_indices


def remove_large_airmass(gb_transform_table, max_airmass=2.0):
    gb_transform_table = gb_transform_table[gb_transform_table['X'] <= 2.0]
    return gb_transform_table


def ground_based_second_order_transforms(gb_transform_table, plot_results=False, save_plots=False, **kwargs):
    """
    Perform the final step in calculating the transforms for a ground-based observatory.

    Parameters
    ----------
    gb_transform_table : astropy.table.Table
        Table containing the results from the intermediary step in calculating the transforms. Has columns:
            field : string list
                Unique identifier of the star field that the reference star is in (e.g. Landolt field "108").
            C_fCI : np.float64
                C coefficient for filter f with colour index ci.
                TODO
            Zprime_f : numpy.float64
                Z' coefficient for filter f.
            filter : string list
                Instrumental filter band to calculate the transform for.
            CI : string list
                Name of the colour index used to calculate c_fci and zprime_f.
            X : numpy.float64
                The mean airmass for all sources in the image.
    plot_results : bool, optional
        Controls whether or not to plot the results from the transforms. The default is False.

    Returns
    -------
    gb_final_transforms : astropy.table.Table
        Table containing the results for the final transforms. Has columns:
            filter : string
                Instrumental filter band used to calculate the transform.
            CI : string
                Name of the colour index used to calculate the transform (e.g. B-V for b, V-R for r).
            k''_fCI : float
                The second order atmospheric extinction coefficient for filter f using the colour index CI.
                TODO
            T_fCI : float
                The instrumental transform coefficient for filter f using the colour index CI.
            k'_f : float
                The first order atmospheric extinction coefficient for filter f.
            Z_f : float
                The zero point magnitude for filter f.

    """
    gb_final_transforms = Table()
    unique_filters = table.unique(gb_transform_table, keys=['filter', 'CI'])
    num_filters = len(unique_filters)
    nan_array = np.empty(num_filters)
    nan_array.fill(np.nan)
    gb_final_transforms = Table(
        names=[
            'filter',
            'CI',
            'k\'\'_fCI',
            'k\'\'_fCI_sigma',
            'T_fCI',
            'T_fCI_sigma',
            'k\'_f',
            'k\'_f_sigma',
            'Z_f',
            'Z_f_sigma'
            ],
        data=[
            np.empty(num_filters, dtype=object),
            np.empty(num_filters, dtype=object),
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array
            ]
        )
    for unique_filter_index, unique_filter_row in enumerate(unique_filters):
        unique_filter = unique_filter_row['filter']
        current_index = unique_filter_row['CI']
        gb_final_transforms['filter'][unique_filter_index] = unique_filter
        gb_final_transforms['CI'][unique_filter_index] = current_index
        mask = ((gb_transform_table['filter'] == unique_filter) & (gb_transform_table['CI'] == current_index))
        current_filter = gb_transform_table[mask]
        x = current_filter['X']
        y = current_filter['C_fCI']
        
        max_c_fci_sigma = max(current_filter['C_fCI_sigma'])
        sigma = np.nan_to_num(current_filter['C_fCI_sigma'], nan=max_c_fci_sigma)
        sigma = np.array(sigma)
        sigma[sigma == 0] = max(sigma)
        
        
        # sigma = current_filter['C_fCI_sigma']
        fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
        fitted_line_c, mask = or_fit(line_init, x, y, weights=1.0/sigma)
        filtered_data_c = np.ma.masked_array(y, mask=mask)
        kprimeprime_fci = fitted_line_c.slope.value
        t_fci = fitted_line_c.intercept.value
        cov_c = fit.fit_info['param_cov']
        # c_fci_sigma = cov[0][0]
        # zprime_f_sigma = cov[1][1]
        # a_fit_c, cov_c = curve_fit(linear_func, current_filter['X'], current_filter['C_fCI'], 
        #                            sigma=current_filter['C_fCI_sigma'])
        # kprimeprime_fci = a_fit_c[0]
        # t_fci = a_fit_c[1]
        if cov_c is None:
            kprimeprime_fci_sigma = np.nan
            t_fci_sigma = np.nan
        else:
            kprimeprime_fci_sigma = sqrt(cov_c[0][0])
            t_fci_sigma = sqrt(cov_c[1][1])
        # kprimeprime_fci, t_fci = np.polyfit(current_filter['X'], current_filter['C_fCI'], 1)
        y = current_filter['Zprime_f']
        
        max_zprime_f_sigma = max(current_filter['Zprime_f_sigma'])
        sigma = np.nan_to_num(current_filter['Zprime_f_sigma'], nan=max_zprime_f_sigma)
        sigma = np.array(sigma)
        sigma[sigma == 0] = max(sigma)
        
        
        # sigma = current_filter['Zprime_f_sigma']
        fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
        fitted_line_z, mask = or_fit(line_init, x, y, weights=1.0/sigma)
        filtered_data_z = np.ma.masked_array(y, mask=mask)
        kprime_f = fitted_line_z.slope.value
        zprime_f = fitted_line_z.intercept.value
        cov_z = fit.fit_info['param_cov']
        # a_fit_z, cov_z = curve_fit(linear_func, current_filter['X'], current_filter['Zprime_f'], 
        #                            sigma=current_filter['Zprime_f_sigma'])
        # kprime_f = a_fit_z[0]
        # zprime_f = a_fit_z[1]
        if cov_z is None:
            kprime_f_sigma = np.nan
            zprime_f_sigma = np.nan
        else:
            kprime_f_sigma = sqrt(cov_z[0][0])
            zprime_f_sigma = sqrt(cov_z[1][1])
        # kprime_f, zprime_f = np.polyfit(current_filter['X'], current_filter['Zprime_f'], 1)
        gb_final_transforms['k\'\'_fCI'][unique_filter_index] = kprimeprime_fci
        gb_final_transforms['k\'\'_fCI_sigma'][unique_filter_index] = kprimeprime_fci_sigma
        gb_final_transforms['T_fCI'][unique_filter_index] = t_fci
        gb_final_transforms['T_fCI_sigma'][unique_filter_index] = t_fci_sigma
        gb_final_transforms['k\'_f'][unique_filter_index] = kprime_f
        gb_final_transforms['k\'_f_sigma'][unique_filter_index] = kprime_f_sigma
        gb_final_transforms['Z_f'][unique_filter_index] = zprime_f
        gb_final_transforms['Z_f_sigma'][unique_filter_index] = zprime_f_sigma
        if plot_results:
            X_plot = np.arange(start=min(current_filter['X']), stop=max(current_filter['X'])+0.01, step=0.01)
            ci_plot = re.sub('[^a-zA-Z]+', '', current_index)
            ci_plot = ci_plot.lower()
            plt.errorbar(x, current_filter['C_fCI'], yerr=current_filter['C_fCI_sigma'], 
                         color='#1f77b4', fmt='o', fillstyle='none', capsize=2, label="Clipped Data")
            plt.plot(x, filtered_data_c, 'o', color='#1f77b4', label="Fitted Data")
            plt.plot(X_plot, fitted_line_c(X_plot), '-', color='#ff7f0e',
                     label=f'C_{unique_filter}{ci_plot} = {kprimeprime_fci:.3f} * X + {t_fci:.3f}')
            plt.legend()
            plt.title(f"k''_{unique_filter}{ci_plot} and T_{unique_filter}{ci_plot} Coefficient calculations")
            # plt.title(f'C_{unique_filter}{ci_plot} = {kprimeprime_fci:.3f} * X + {t_fci:.3f}')
            plt.ylabel(f'C_{unique_filter}{ci_plot}')
            plt.xlabel('X')
            if save_plots:
                save_loc = f"{os.path.join(kwargs.get('save_loc'), f'KprimeprimeT{unique_filter}{ci_plot}')}.png"
                plt.savefig(save_loc)
            plt.show()
            plt.close()
            plt.errorbar(x, current_filter['Zprime_f'], yerr=current_filter['Zprime_f_sigma'], 
                         color='#1f77b4', fmt='o', fillstyle='none', capsize=2, label="Clipped Data")
            plt.plot(x, filtered_data_z, 'o', color='#1f77b4', label="Fitted Data")
            plt.plot(X_plot, fitted_line_z(X_plot), '-', color='#ff7f0e',
                     label=f'Z\'_{unique_filter} = {kprime_f:.3f} * X + {zprime_f:.3f}')
            plt.legend()
            plt.title(f"Z_{unique_filter} and k'_{unique_filter} Coefficient Calculations")
            # plt.title(f'Z\'_{unique_filter} = {kprime_f:.3f} * X + {zprime_f:.3f}')
            plt.ylabel(f'Z\'_{unique_filter}')
            plt.xlabel('X')
            if save_plots:
                save_loc = f"{os.path.join(kwargs.get('save_loc'), f'KprimeZ{unique_filter}{ci_plot}')}.png"
                plt.savefig(save_loc)
            plt.show()
            plt.close()
    
    return gb_final_transforms


# TODO
# This will be harder than I thought. Namely, keeping track of which star was which when it's not known (unless maybe 
# I do assume to know it because it is just a check? Maybe by using matched_stars or large_stars_table instead of 
# creating a whole new table?). Either way, I'll start to convert the light curve code first and then come back to this 
# later.


def get_colour_index_lower(instr_filter):
    """
    Get the colour index for the desired instrumental filter band.

    Parameters
    ----------
    instr_filter : TYPE
        DESCRIPTION.

    Returns
    -------
    colour_index : TYPE
        DESCRIPTION.
    ci : TYPE
        DESCRIPTION.

    """
    if instr_filter == 'b':
        colour_index = 'B-V'
    elif instr_filter == 'v' or instr_filter == 'g':
        colour_index = 'B-V'
    elif instr_filter == 'r':
        colour_index = 'V-R'
    elif instr_filter == 'i':
        colour_index = 'V-I'
    else:
        colour_index = None
    ci = re.sub('[^a-zA-Z]+', '', colour_index)
    ci = ci.lower()
    return colour_index, ci


def init_unknown_object_table_columns():
    """
    Initialize the columns that will create the unknown object table.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    instr_filter = []
    name = []
    instr_mag = []
    
    unknown_object_table_columns = namedtuple('unknown_object_table_columns', 
                                              ['instr_filter',
                                              'name',
                                              'instr_mag'])
    return unknown_object_table_columns(instr_filter, 
                                        name, 
                                        instr_mag)


def update_unknown_table_columns(unknown_object_table_columns, hdr, instr_mag, filter_key='FILTER'):
    """
    Update the columns tha will create the unknown object table.

    Parameters
    ----------
    unknown_object_table_columns : TYPE
        DESCRIPTION.
    hdr : TYPE
        DESCRIPTION.
    instr_mag : TYPE
        DESCRIPTION.
    filter_key : TYPE, optional
        DESCRIPTION. The default is 'FILTER'.

    Returns
    -------
    updated_unknown_object_table_columns : TYPE
        DESCRIPTION.

    """
    updated_unknown_object_table_columns = unknown_object_table_columns
    instr_filter = get_instr_filter_name(hdr=hdr, filter_key=filter_key)
    name = hdr['OBJECT']
    updated_unknown_object_table_columns.instr_filter.append(instr_filter)
    updated_unknown_object_table_columns.name.append(name)
    updated_unknown_object_table_columns.instr_mag.append(instr_mag)
    return updated_unknown_object_table_columns

def calculate_c_fci(gb_final_transforms, instr_filter, airmass, colour_index):
    mask = ((gb_final_transforms['filter'] == instr_filter) & (gb_final_transforms['CI'] == colour_index))
    row_of_transforms = gb_final_transforms[mask]
    if len(row_of_transforms) == 0 and instr_filter == 'v':
        instr_filter = 'g'
        mask = ((gb_final_transforms['filter'] == instr_filter) & (gb_final_transforms['CI'] == colour_index))
        row_of_transforms = gb_final_transforms[mask]
    c_fci = float(row_of_transforms['T_fCI']) - (float(row_of_transforms["k''_fCI"]) * airmass)
    return c_fci


def calculate_c_prime(gb_final_transforms, instr_filter, airmass):
    """
    Calculate the C' coefficient to use when applying the transforms.

    Parameters
    ----------
    gb_final_transforms : TYPE
        DESCRIPTION.
    instr_filter : TYPE
        DESCRIPTION.

    Returns
    -------
    c_prime_fci : TYPE
        DESCRIPTION.
        TODO

    """
    colour_index, ci = get_colour_index_lower(instr_filter)
    c_numerator = calculate_c_fci(gb_final_transforms, instr_filter, airmass, colour_index)
    c_negative = calculate_c_fci(gb_final_transforms, ci[0], airmass, colour_index)
    c_positive = calculate_c_fci(gb_final_transforms, ci[1], airmass, colour_index)
    c_prime_fci = c_numerator / (1 - c_negative + c_positive)
    return c_prime_fci


def calculate_z_prime_f(gb_final_transforms, instr_filter, airmass, colour_index):
    mask = ((gb_final_transforms['filter'] == instr_filter) & (gb_final_transforms['CI'] == colour_index))
    row_of_transforms = gb_final_transforms[mask]
    if len(row_of_transforms) == 0 and instr_filter == 'v':
        instr_filter = 'g'
        mask = ((gb_final_transforms['filter'] == instr_filter) & (gb_final_transforms['CI'] == colour_index))
        row_of_transforms = gb_final_transforms[mask]
    z_prime_f = float(row_of_transforms['Z_f']) - (float(row_of_transforms["k'_f"]) * airmass)
    return z_prime_f


def calculate_lower_z_f(gb_final_transforms, c_prime_fci, instr_filter, airmass):
    colour_index, ci = get_colour_index_lower(instr_filter)
    c_prime_fci = calculate_c_prime(gb_final_transforms, instr_filter, airmass)
    z_prime_positive = calculate_z_prime_f(gb_final_transforms, ci[0], airmass, colour_index)
    z_prime_negative = calculate_z_prime_f(gb_final_transforms, ci[1], airmass, colour_index)
    z_prime_no_brackets = calculate_z_prime_f(gb_final_transforms, instr_filter, airmass, colour_index)
    lower_z_f = c_prime_fci * (z_prime_positive - z_prime_negative) + z_prime_no_brackets
    return lower_z_f


def apply_gb_transforms_VERIFICATION(gb_final_transforms, stars_table, instr_filter):
    """
    Apply the transforms to a source with an unknown standard magnitude.

    Parameters
    ----------
    gb_final_transforms : astropy.table.Table
        Table containing the results for the final transforms. Has columns:
            filter : string
                Instrumental filter band used to calculate the transform.
            CI : string
                Name of the colour index used to calculate the transform (e.g. B-V for b, V-R for r).
            k''_fCI : float
                The second order atmospheric extinction coefficient for filter f using the colour index CI.
            T_fCI : float
                The instrumental transform coefficient for filter f using the colour index CI.
            k'_f : float
                The first order atmospheric extinction coefficient for filter f.
            Z_f : float
                The zero point magnitude for filter f.
    unknown_object_table : astropy.table.Table
        Table containing the stars to calculate the apparent magnitudes for. Has columns:
            filter : string
                Instrumental filter band used to take the image.
            name : string
                Unique identifier of the source to apply the transform to.
            instrumental mag : float
                Instrumental magnitude of the source.
                TODO: Change this so that it accepts multiple sources somehow.

    Returns
    -------
    app_mag_table : astropy.table.Table
        Table containing the apparent magnitudes of the object after applying the transforms. Has columns:
            time : float
                Julian date that the image was taken.
            filter : string
                Apparent filter that the image was transformed to.
            apparent mag : float
                Apparent magnitude of the source after applying the transforms.

    """
    # If there is only 1 entry per filter per source, assume that they are averages and don't need any interpolation.
    app_mag_first_columns = Table(stars_table['Field', 'Name', 'V_ref', '(B-V)', '(U-B)', '(V-R)', '(V-I)', 'V_sigma'])
    colour_index, ci = get_colour_index_lower(instr_filter)
    instr_mag = stars_table[instr_filter]
    airmass = stars_table[f'X_{instr_filter}']
    c_prime_fci = calculate_c_prime(gb_final_transforms, instr_filter, airmass)
    try:
        positive_instr_mag = stars_table[ci[0]]
        negative_instr_mag = stars_table[ci[1]]
    except KeyError:
        table_ci = ci.replace('v', 'g')
        positive_instr_mag = stars_table[table_ci[0]]
        negative_instr_mag = stars_table[table_ci[1]]
    lower_z_f = calculate_lower_z_f(gb_final_transforms, c_prime_fci, instr_filter, airmass)
    app_mag_list = instr_mag + c_prime_fci * (positive_instr_mag - negative_instr_mag) + lower_z_f
    app_mag_filter = instr_filter.upper()
    app_mag_column = Table(names=[app_mag_filter], data=[app_mag_list])
    app_mag_table = table.hstack([app_mag_first_columns, app_mag_column])
    return app_mag_table


def apply_gb_timeseries_transforms(gb_final_transforms, large_sats_table):
    """
    Apply the transforms to a timeseries of observations.

    Parameters
    ----------
    gb_final_transforms : astropy.table.Table
        Table containing the results for the final transforms. Has columns:
            filter : string
                Instrumental filter band used to calculate the transform.
            CI : string
                Name of the colour index used to calculate the transform (e.g. B-V for b, V-R for r).
            k''_fCI : float
                The second order atmospheric extinction coefficient for filter f using the colour index CI.
            T_fCI : float
                The instrumental transform coefficient for filter f using the colour index CI.
            k'_f : float
                The first order atmospheric extinction coefficient for filter f.
            Z_f : float
                The zero point magnitude for filter f.
    large_sats_table : astropy.table.Table
        DESCRIPTION.

    Returns
    -------
    app_large_sats_table : TYPE
        DESCRIPTION.

    """
    app_large_sats_table = large_sats_table
    return app_large_sats_table


def copy_and_rename(directory, 
                    file_suffix=".fits", 
                    time_key='DATE-OBS', 
                    filter_key='FILTER', 
                    temp_dir='tmp', 
                    debugging=False):
    if not debugging:
        try:
            os.mkdir(temp_dir)
        except FileExistsError as e:
            print(e)
    else:
        if os.path.exists(temp_dir):
            rmtree(temp_dir)
        os.mkdir(temp_dir)
    
    filecount = 0
    for dirpth, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_suffix):
                with fits.open(os.path.join(dirpth, file)) as image:
                    hdr = image[0].header
                t = Time(hdr[time_key], format='fits', scale='utc')
                filter_name = hdr[filter_key]
                t_datetime = t.to_datetime()
                new_filename = f'{t_datetime.strftime("%Y%m%d%H%M%S")}{filter_name}.fits'
                copy2(os.path.join(dirpth, file), f'{temp_dir}/{new_filename}')
                filecount += 1
    filenames = sorted(os.listdir(temp_dir))
    return filecount, filenames


def remove_temp_dir(temp_dir='tmp'):
    rmtree(temp_dir)


def set_sat_positions(imgdata, filecount, set_sat_positions_bool, max_distance_from_sat=25, norm=LogNorm(), cmap_set='Set1'):
    def mbox(title, text, style):
        return ctypes.windll.user32.MessageBoxW(0, text, title, style)
    def set_sat_position(event, x, y, flags, params):
        global sat_locs
        if event == cv.EVENT_LBUTTONDOWN:
            sat_locs.append([x, y])
    def return_entry(event=None):
        """Get and print the content of the entry."""
        # global entry
        global content
        content = entry.get()
        root.destroy()
    # global set_sat_positions_bool
    while set_sat_positions_bool:
        global sat_locs
        sat_locs = []
        mbox('Information',
             'Please select the positions of the satellites on the following image. Press any key when finished.',
             0)
        cv.namedWindow('TestImage')
        cv.setMouseCallback('TestImage', set_sat_position)
        logdata = cv.normalize(imgdata, None, alpha=0, beta=10, norm_type=cv.NORM_MINMAX, dtype=cv.CV_32F)
        cv.imshow('TestImage', logdata)
        cv.waitKey(0)
        cv.destroyAllWindows()
        sat_locs = np.array(sat_locs)
        print(sat_locs)
        num_sats = len(sat_locs)
        num_nans = np.zeros(num_sats, dtype=int)
        names = np.empty(num_sats + 2, dtype=object)
        names[0] = 'Time (JD)'
        names[1] = 'Filter'
        date_col = np.empty((filecount, 1))
        date_col.fill(np.nan)
        filter_col = np.empty((filecount, 1), dtype=object)
        data = np.empty((filecount, num_sats))
        data.fill(np.nan)

        for i, name in enumerate(names[2:]):
            fig = Figure()
            ax = fig.add_subplot()
            ax.imshow(imgdata, cmap='gray', norm=norm, interpolation='nearest')
            sat_aperture = RectangularAperture(sat_locs[i], w=max_distance_from_sat * 2,
                                               h=max_distance_from_sat * 2)
            sat_aperture.plot(axes=ax, color='r', lw=1.5, alpha=0.5)
            root = tk.Tk()
            root.title("Set Satellite Position")
            img_frame = tk.Frame(root)
            img_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
            input_frame = tk.Frame(root)
            input_frame.pack(side=tk.RIGHT)
            label = tk.Label(input_frame, text='Enter Satellite')
            label.pack()
            entry = tk.Entry(input_frame)
            entry.bind("<Return>", return_entry)
            entry.pack(padx=5)
            button = tk.Button(input_frame, text="OK", command=return_entry)
            button.pack()
            canvas = FigureCanvasTkAgg(fig, master=img_frame)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, img_frame)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
            root.update()
            root.focus_force()
            entry.focus_set()
            root.mainloop()
            names[i + 2] = content
            print(f"Satellite {names[i + 2]} at location ({sat_locs[i, 0]}, {sat_locs[i, 1]})")
        print(names)
        
        cmap = plt.get_cmap(cmap_set) # TODO: move this into a separate function?
        colours = [cmap(i) for i in range(0, num_sats)]
        legend_elements = []
        window = tk.Tk()
        window.title('Plotting in Tkinter Test')
        img_frame = tk.Frame(window)
        img_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        input_frame = tk.Frame(window)
        input_frame.pack(side=tk.RIGHT)
        label = tk.Label(input_frame, text='Are the satellite positions correct?')
        label.pack()
        yes_no = tk.IntVar()
        yes_btn = tk.Radiobutton(input_frame, text='Yes', variable=yes_no, value=1)
        yes_btn.pack(anchor=tk.W, padx=5)
        no_btn = tk.Radiobutton(input_frame, text='No', variable=yes_no, value=2)
        no_btn.pack(anchor=tk.W, padx=5)
        closebutton = tk.Button(input_frame, text='OK', command=window.destroy)
        closebutton.pack()
        fig = Figure()
        ax = fig.add_subplot()
        ax.imshow(imgdata, cmap='gray', norm=LogNorm(), interpolation='nearest')
        for i in range(0, num_sats):
            sat_aperture = RectangularAperture(sat_locs[i], w=max_distance_from_sat * 2,
                                               h=max_distance_from_sat * 2)
            sat_aperture.plot(axes=ax, color=colours[i], lw=1.5, alpha=0.5)
            legend_elements.append(Line2D([0], [0], color='w', marker='s', markerfacecolor=colours[i], markersize=7,
                                          label=names[i + 2]))
        fig.legend(handles=legend_elements, framealpha=1)
        canvas = FigureCanvasTkAgg(fig, master=img_frame)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, img_frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        window.mainloop()
        # Have a way for the user to confirm the satellite locations. If it is wrong, then decide whether to change
        # set_sat_positions to True/False or change_sat_positions
        if yes_no.get() == 1:
            set_sat_positions_bool = False
        else:
            continue
        # elif yes_no.get() == 2:
        #     sat_locs = []
        #     names = []
        # print(names[0])
        # print(date_col)
        sat_names = names[2:]
        date_table = Table(names=[names[0]], data=date_col)
        filter_table = Table(names=[names[1]], data=filter_col)
        data_table = Table(names=names[2:], data=data)
        sats_table = hstack([date_table, filter_table, data_table], join_type='exact')
        uncertainty_table = hstack([date_table, filter_table, data_table], join_type='exact')
        sat_fwhm_table = hstack([date_table, filter_table, data_table], join_type='exact')
        sats_table.pprint_all()
    # TODO: Figure out what needs to be returned from this function and return it (maybe as a nametuple?).
    sat_information = namedtuple('sat_information', 
                                 ['sats_table',
                                  'uncertainty_table',
                                  'sat_fwhm_table',
                                  'sat_locs',
                                  'num_sats',
                                  'num_nans',
                                  'sat_names'])
    return set_sat_positions_bool, sat_information(sats_table, 
                                                   uncertainty_table, 
                                                   sat_fwhm_table, 
                                                   sat_locs, 
                                                   num_sats, 
                                                   num_nans, 
                                                   sat_names)
    # return sats_table, uncertainty_table, sat_fwhm_table, sat_locs, num_sats, num_nans, sat_names


def plot_detected_sats(filename,
                       plot_results, 
                       imgdata, 
                       irafsources, 
                       sat_information, 
                       max_distance_from_sat=20, 
                       norm=LogNorm()):
    if plot_results != 0:
        fig, ax = plt.subplots()
        ax.imshow(imgdata, cmap='gray', norm=norm, interpolation='nearest')
        ax.scatter(irafsources['xcentroid'], irafsources['ycentroid'],
                    s=100, edgecolor='red', facecolor='none')
        for i in range(0, sat_information.num_sats):
            rect = patches.Rectangle((sat_information.sat_locs[i,0] - max_distance_from_sat, 
                                      sat_information.sat_locs[i,1] - max_distance_from_sat),
                                     width=max_distance_from_sat * 2, height=max_distance_from_sat * 2,
                                     edgecolor='green', facecolor='none')
            ax.add_patch(rect)
        plt.title(filename)
        if plot_results == 1:
            plt.show(block=False)
            plt.pause(2)
        elif plot_results == 2:
            plt.show()
        plt.close()


def check_if_sat(sat_information, index, irafsources, instr_mags, instr_mags_sigma, fwhms_arcsec, max_distance_from_sat=20):
    for obj_index, obj in enumerate(irafsources):
        obj_x = obj['xcentroid']
        obj_y = obj['ycentroid']
        for sat_num, sat in enumerate(sat_information.sat_locs, start=2):
            sat_x = sat[0]
            sat_y = sat[1]
            if abs(sat_x - obj_x) < max_distance_from_sat and abs(sat_y - obj_y) < max_distance_from_sat:
                sat_information.sats_table[index][sat_num] = instr_mags[obj_index]
                sat_information.uncertainty_table[index][sat_num] = instr_mags_sigma[obj_index]
                sat_information.sat_fwhm_table[index][sat_num] = fwhms_arcsec[obj_index]
                sat[0] = obj_x
                sat[1] = obj_y
    print(sat_information.sats_table[index])
    return sat_information


def determine_if_change_sat_positions(sat_information, filenum, change_sat_positions_bool, max_num_nan=5):
    sat_mags = np.array(list(sat_information.sats_table[filenum]))
    mask = np.isnan(sat_mags[2:].astype(float))
    sat_information.num_nans[mask] += 1
    sat_information.num_nans[~mask] = 0
    print(sat_information.num_nans)
    num_nan = max(sat_information.num_nans)
    # print(num_nan)
    if num_nan >= max_num_nan:
        change_sat_positions_bool = True
    return change_sat_positions_bool, num_nan


def change_sat_positions(filenames, 
                         filenum, 
                         num_nan, 
                         sat_information,
                         change_sat_positions_bool,
                         max_distance_from_sat=20, 
                         size=25, 
                         temp_dir='tmp', 
                         cmap_set='Set1', 
                         plot_results=0):
    # Change the position
    # Display filenames[filenum - num_nan]
    def mbox(title, text, style):
        return ctypes.windll.user32.MessageBoxW(0, text, title, style)
    def change_sat_position(event, x, y, flags, params):
        global sat_locs
        if event == cv.EVENT_LBUTTONDOWN:
            sat_locs[params[0]] = [x, y]
            cv.destroyAllWindows()
    print(filenames[filenum - num_nan])
    # sat_checked = np.zeros(num_sats, dtype=IntVar())
    with fits.open(f"{temp_dir}/{filenames[filenum - num_nan]}") as image:
        hdr = image[0].header
        imgdata = image[0].data
    root = tk.Tk()
    root.title("Current satellite positions")
    input_frame = tk.Frame(root)
    input_frame.pack(side=tk.RIGHT)
    img_frame = tk.Frame(root)
    img_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    label = tk.Label(input_frame, text='Select the satellite(s) whose position you would like to change.')
    label.grid(row=0)
    sat_checked = []
    for sat_num, sat in enumerate(sat_information.sat_names):
        sat_checked.append(tk.IntVar())
        checkbutton = tk.Checkbutton(input_frame, text=sat, variable=sat_checked[sat_num])
        checkbutton.grid(row=sat_num+1, sticky=tk.W, padx=5)
    none_select = tk.IntVar()
    checkbutton = tk.Checkbutton(input_frame, text="None", variable=none_select)
    checkbutton.grid(row=sat_information.num_sats+2, sticky=tk.W, padx=5)
    closebutton = tk.Button(input_frame, text='OK', command=root.destroy)
    closebutton.grid(row=sat_information.num_sats+3)
    legend_elements = []
    fig = Figure()
    ax = fig.add_subplot()
    ax.imshow(imgdata, cmap='gray', norm=LogNorm(), interpolation='nearest')
    cmap = plt.get_cmap(cmap_set) # TODO: move this into a separate function?
    colours = [cmap(i) for i in range(0, sat_information.num_sats)]
    for i in range(0, sat_information.num_sats):
        sat_aperture = RectangularAperture(sat_locs[i], w=max_distance_from_sat * 2,
                                           h=max_distance_from_sat * 2)
        sat_aperture.plot(axes=ax, color=colours[i], lw=1.5, alpha=0.5)
        legend_elements.append(Line2D([0], [0], color='w', marker='s', markerfacecolor=colours[i], markersize=7,
                                      label=sat_information.sat_names[i]))
    fig.legend(handles=legend_elements, framealpha=1)
    canvas = FigureCanvasTkAgg(fig, master=img_frame)
    canvas.draw()
    toolbar = NavigationToolbar2Tk(canvas, img_frame)
    toolbar.update()
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    root.mainloop()
    if none_select.get() == 1:
        none_sats = True
    elif none_select.get() == 0:
        none_sats = False
    if not none_sats:
        sat_checked_int = np.empty(len(sat_checked), dtype=int)
        for sat_num, sat in enumerate(sat_checked):
            sat_checked_int[sat_num] = sat.get()
        print(sat_checked_int)
        sat_checked_mask = sat_checked_int == 1
        print(sat_checked_mask)
        for sat in sat_information.sat_names[sat_checked_mask]:
            print(sat)
            index = np.where(sat_information.sat_names == sat)
            mbox('Information',
                 f'Please select the new position of {sat} on the following image.',
                 0)
            cv.namedWindow('TestImage')
            cv.setMouseCallback('TestImage', change_sat_position, index)
            logdata = cv.normalize(imgdata, None, alpha=0, beta=10, norm_type=cv.NORM_MINMAX, dtype=cv.CV_32F)
            cv.imshow('TestImage', logdata)
            cv.waitKey(0)
        print(sat_locs)
        # If some selection is none
        # none_sats = True
        for reversing_index in range(1, num_nan+1):
            filepath = f"{temp_dir}/{filenames[filenum - reversing_index]}"
            hdr, imgdata = read_fits_file(filepath)
            print(filenames[filenum - reversing_index])
            filename = filenames[filenum - reversing_index]
            # Do the photometry on the images. Should change this to a bunch of function calls, rather than copy
            # and pasting the code.
            exptime = hdr['EXPTIME'] * u.s  # Store the exposure time with unit seconds.
            bkg, bkg_std = calculate_img_bkg(imgdata)
            # mean_val, median_val, std_val = sigma_clipped_stats(imgdata)  # Calculate background stats.
            irafsources = detecting_stars(imgdata, bkg, bkg_std)
            # iraffind = IRAFStarFinder(threshold=4*std_val, fwhm=2)  # Find stars using IRAF.
            # irafsources = iraffind(imgdata - median_val)  # Subtract background median value.
            if not irafsources:
                sat_information.num_nans[:] = 0
                continue
            # try:
            #     irafsources.sort('flux', reverse=True)  # Sort the stars by flux, from greatest to least.
            # except Exception as e:
            #     # Reset the NaN boolean as it doesn't count.
            #     num_nans[:] = 0
            #     continue
            # irafpositions = np.transpose((irafsources['xcentroid'],
            #                               irafsources['ycentroid']))  # Store source positions as a numpy array.
            plot_detected_sats(filename,
                               plot_results, 
                               imgdata, 
                               irafsources, 
                               sat_information, 
                               max_distance_from_sat=max_distance_from_sat, 
                               norm=LogNorm())
            # if plot_results != 0:
            #     fig, ax = plt.subplots()
            #     ax.imshow(imgdata, cmap='gray', norm=LogNorm(), interpolation='nearest')
            #     ax.scatter(irafsources['xcentroid'], irafsources['ycentroid'],
            #                s=100, edgecolor='red', facecolor='none')
            #     for i in range(0, sat_information.num_sats):
            #         rect = patches.Rectangle(
            #             (sat_locs[i, 0] - max_distance_from_sat, sat_locs[i, 1] - max_distance_from_sat),
            #             width=max_distance_from_sat * 2, height=max_distance_from_sat * 2,
            #             edgecolor='green', facecolor='none')
            #         ax.add_patch(rect)
            #     plt.title(filenames[filenum - reversing_index])
            #     if plot_results == 1:
            #         plt.show(block=False)
            #         plt.pause(2)
            #     elif plot_results == 2:
            #         plt.show()
            #     plt.close()
            # TODO: Calculate array of FWHM to be used later.
            fwhms, fwhm, fwhm_std = calculate_fwhm(irafsources)
            fwhms_arcsec, fwhm_arcsec, fwhm_std_arcsec = convert_fwhm_to_arcsec(hdr, fwhms, fwhm, fwhm_std)
            # iraf_fwhms = irafsources['fwhm']  # Save FWHM in a list.
            # iraf_fwhms = np.array(iraf_fwhms)  # Convert FWHM list to numpy array.
            # # Print information about the file                                                                              #####
            # # Calculate statistics of the FWHMs given by the loop over each source.                                         #####
            # iraf_sdom = iraf_fwhms.std() / sqrt(len(iraf_fwhms))  # Standard deviation of the mean. Not used.
            # num_IRAF_sources = len(irafsources)  # Number of stars IRAF found.
            # print(f"No. of IRAF sources: {num_IRAF_sources}")  # Print number of IRAF stars found.
            # iraf_fwhm = iraf_fwhms.mean()  # Calculate IRAF FWHM mean.
            # iraf_std = iraf_fwhms.std()  # Calculate IRAF standard deviation.
            # print(f"IRAF Calculated FWHM (pixels): {iraf_fwhm:.3f} +/- {iraf_std:.3f}")  # Print IRAF FWHM.
            # Calculate the flux of each star using PSF photometry. This uses the star positions calculated by
            # IRAFStarFinder earlier.
            photometry_result = perform_photometry(irafsources, fwhm, imgdata, bkg)
            # daogroup = DAOGroup(2 * iraf_fwhm)  # Groups overlapping stars together.
            # # sigma_clip = SigmaClip()
            # # mmm_bkg = MMMBackground(sigma_clip=sigma_clip)
            # # bkg_value = mmm_bkg(imgdata)
            # psf_model = IntegratedGaussianPRF(
            #     sigma=iraf_fwhm * gaussian_fwhm_to_sigma)  # Defime the PSF model to be used for photometry.
            # psf_model.x_0.fixed = True  # Don't change the initial 'guess' of the star x positions to be provided.
            # psf_model.y_0.fixed = True  # Don't change the initial 'guess' of the star y positions to be provided.
            # # Provide the initial guesses for the x-y positions and the flux. The flux will be fit using the psf_model,
            # # so that will change, but the star positions will remain the same.
            # pos = Table(names=['x_0', 'y_0', 'flux_0'],
            #             data=[irafsources['xcentroid'], irafsources['ycentroid'], irafsources['flux']])
            # # Initialize the photometry to be performed. Do not estimate the background, as it will be subtracted from
            # # the image when fitting the PSF.
            # photometry = BasicPSFPhotometry(group_maker=daogroup,
            #                                 bkg_estimator=None,
            #                                 psf_model=psf_model,
            #                                 fitter=LevMarLSQFitter(),
            #                                 fitshape=size)
            # Perform the photometry on the background subtracted image. Also pass the fixed x-y positions and the
            # initial guess for the flux.
            # result_tab = photometry(image=imgdata - median_val, init_guesses=pos)
            # fluxes = photometry_result['flux_fit']  # Store the fluxes as a list.
            instr_mags = calculate_magnitudes(photometry_result, exptime)
            instr_mags_sigma = calculate_magnitudes_sigma(photometry_result, exptime)
            # fluxes = np.array(
            #     fluxes) * u.ct  # Convert the fluxes to a numpy array and add the unit of count to it.
            # fluxes = fluxes / exptime  # Normalize the fluxes by exposure time (unit is now counts / second)
            # flux_uncs = result_tab['flux_unc']
            # flux_uncs = np.array(flux_uncs) * u.ct
            # flux_uncs = flux_uncs / exptime
            # snr = (fluxes / flux_uncs).value
            # instr_mags_sigma = 1.0857 / np.sqrt(snr)
            # instr_mags_units = u.Magnitude(fluxes)  # Convert the fluxes to an instrumental magnitude.
            # instr_mags = instr_mags_units.value  # Store the magnitudes without the unit attached.
            # Calculate the FWHM in units of arcseconds as opposed to pixels.
            # TODO: Calculate FWHM in arcsec (create function if one doesn't already exist).
            # try:
            #     focal_length = hdr['FOCALLEN'] * u.mm  # Store the telescope's focal length with unit millimetres.
            #     xpixsz = hdr['XPIXSZ']  # Store the size of the x pixels.
            #     ypixsz = hdr['XPIXSZ']  # Store the size of the y pixels.
            #     if xpixsz == ypixsz:  # If the pixels are square.
            #         pixsz = xpixsz * u.um  # Store the pixel size with unit micrometre.
            #         # Can find FOV by finding deg/pix and then multiplying by the x and y number of pix (NAXIS).
            #         rad_per_pix = atan(
            #             pixsz / focal_length) * u.rad  # Calculate the angular resolution of each pixel. Store with unit radians.
            #         arcsec_per_pix = rad_per_pix.to(
            #             u.arcsec)  # Convert the per pixel angular resultion to arcseconds.
            #         iraf_FWHM_arcsec = iraf_fwhm * arcsec_per_pix.value  # Convert the IRAFStarFinder FWHM from pixels to arcsec.
            #         iraf_std_arcsec = iraf_std * arcsec_per_pix  # Convert the IRAFStarFinder FWHM standard deviation from pixels to arcsec.
            #         iraf_FWHMs_arcsec = iraf_fwhms * arcsec_per_pix.value
            #         print(
            #             f"IRAF Calculated FWHM (arcsec): {iraf_FWHM_arcsec:.3f} +/- {iraf_std_arcsec:.3f}")  # Print the IRAFStarFinder FWHM in arcsec.
            # except KeyError:
            #     iraf_FWHMs_arcsec = iraf_fwhms
            # print(irafsources['peak'] + median_val)                                                                         # Akin to 'max_pixel' from Shane's spreadsheet.
            # print(result_tab['x_0', 'y_0', 'flux_fit', 'flux_unc'])                                                         # Print the fluxes and their uncertainty for the current image.
            sat_information = check_if_sat(sat_information, 
                                           filenum - reversing_index, 
                                           irafsources, 
                                           instr_mags, 
                                           instr_mags_sigma,
                                           fwhms_arcsec,
                                           max_distance_from_sat=max_distance_from_sat)
            # for obj_index, obj in enumerate(irafsources):
            #     obj_x = obj['xcentroid']
            #     obj_y = obj['ycentroid']
            #     for sat_num, sat in enumerate(sat_locs, start=2):
            #         sat_x = sat[0]
            #         sat_y = sat[1]
            #         if abs(sat_x - obj_x) < max_distance_from_sat and abs(
            #                 sat_y - obj_y) < max_distance_from_sat:
            #             sat_information.sats_table[filenum - reversing_index][sat_num] = instr_mags[obj_index]
            #             sat_information.uncertainty_table[filenum - reversing_index][sat_num] = instr_mags_sigma[obj_index]
            #             # TODO: array FWHMs
            #             # sat_information.sat_fwhm_table[filenum - reversing_index][sat_num] = iraf_FWHMs_arcsec[obj_index]
            print(sat_information.sats_table[filenum - reversing_index])
        sat_information.num_nans[sat_checked_mask] = 0
    change_sat_positions_bool = False
    return change_sat_positions_bool, sat_information


def add_new_time_and_filter(hdr, sat_information, filenum):
    t = Time(hdr['DATE-OBS'], format='fits', scale='utc')
    sat_information.sats_table['Time (JD)'][filenum] = t.jd
    try:
        sat_information.sats_table['Filter'][filenum] = hdr['FILTER']
    except KeyError:
        sat_information.sats_table['Filter'][filenum] = 'C'
    sat_information.uncertainty_table['Time (JD)'][filenum] = t.jd
    try:
        sat_information.uncertainty_table['Filter'][filenum] = hdr['FILTER']
    except KeyError:
        sat_information.uncertainty_table['Filter'][filenum] = 'C'
    sat_information.sat_fwhm_table['Time (JD)'][filenum] = t.jd
    try:
        sat_information.sat_fwhm_table['Filter'][filenum] = hdr['FILTER']
    except KeyError:
        sat_information.sat_fwhm_table['Filter'][filenum] = 'C'
    return sat_information


def _main_gb_transform_calc(directory, 
                            ref_stars_file, 
                            plot_results=False, 
                            save_plots=False, 
                            remove_large_airmass_bool=True, 
                            file_suffix=".fits", 
                            exposure_key='EXPTIME', 
                            lat_key='SITELAT', 
                            lon_key='SITELONG', 
                            elev_key='SITEELEV', 
                            name_key='Name',
                            **kwargs):
    # TODO: Docstring.
    reference_stars, ref_star_positions = read_ref_stars(ref_stars_file)
    gb_transform_table_columns = init_gb_transform_table_columns()
    
    if save_plots:
        save_loc = kwargs.get('save_loc')
    
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(file_suffix):
                filepath = os.path.join(dirpath, filename)
                hdr, imgdata = read_fits_file(filepath)
                exptime = hdr[exposure_key]
                bkg, bkg_std = calculate_img_bkg(imgdata)
                irafsources = detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
                if not irafsources:
                    continue
                _, fwhm, fwhm_std = calculate_fwhm(irafsources)
                photometry_result = perform_photometry(irafsources, fwhm, imgdata, bkg=bkg)
                fluxes = np.array(photometry_result['flux_fit'])
                instr_mags = calculate_magnitudes(photometry_result, exptime)
                instr_mags_sigma = calculate_magnitudes_sigma(photometry_result, exptime)
                wcs = WCS(hdr)
                skypositions = convert_pixel_to_ra_dec(irafsources, wcs)
                altazpositions = None
                try:
                    altazpositions = convert_ra_dec_to_alt_az(skypositions, hdr, lat_key=lat_key, 
                                                              lon_key=lon_key, elev_key=elev_key)
                except AttributeError as e:
                    print(e)
                    continue
                matched_stars = find_ref_stars(reference_stars, 
                                               ref_star_positions,
                                               skypositions,
                                               instr_mags,
                                               instr_mags_sigma,
                                               fluxes,
                                               ground_based=True,
                                               altazpositions=altazpositions)
                if not matched_stars:
                    continue
                
                instr_filter = get_instr_filter_name(hdr)
                colour_indices = get_all_colour_indices(instr_filter)
                for colour_index in colour_indices:
                    field = get_field_name(matched_stars, name_key=name_key)
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
                    if not save_plots:
                        c_fci, c_fci_sigma, zprime_f, zprime_f_sigma = ground_based_first_order_transforms(matched_stars, 
                                                                                                           instr_filter, 
                                                                                                           colour_index, 
                                                                                                           plot_results=plot_results)
                    else:
                        unique_id = filename
                        c_fci, c_fci_sigma, zprime_f, zprime_f_sigma = ground_based_first_order_transforms(matched_stars, 
                                                                                                           instr_filter, 
                                                                                                           colour_index, 
                                                                                                           plot_results=plot_results,
                                                                                                           save_plots=save_plots,
                                                                                                           save_loc=save_loc,
                                                                                                           unique_id=unique_id)
                    gb_transform_table_columns = update_gb_transform_table_columns(gb_transform_table_columns,
                                                                                   field,
                                                                                   c_fci,
                                                                                   c_fci_sigma,
                                                                                   zprime_f,
                                                                                   zprime_f_sigma,
                                                                                   instr_filter,
                                                                                   colour_index,
                                                                                   altazpositions)
    gb_transform_table = create_gb_transform_table(gb_transform_table_columns)
    if remove_large_airmass_bool:
        gb_transform_table = remove_large_airmass(gb_transform_table)
    if save_plots:
        gb_final_transforms = ground_based_second_order_transforms(gb_transform_table, 
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
        write_table_to_latex(gb_final_transforms, f"{os.path.join(save_loc, 'gb_final_transforms')}.txt", formats=formats)
    else:
        gb_final_transforms = ground_based_second_order_transforms(gb_transform_table, 
                                                                   plot_results=plot_results, 
                                                                   save_plots=save_plots)
    # Test the Transforms.
    # directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021_J132_46927_DESCENT\May 18 2021\Landolt Fields\Solved TEST'
    # large_table_columns = init_large_table_columns()
    # for dirpath, dirnames, filenames in os.walk(directory):
    #     for filename in filenames:
    #         if filename.endswith(".fit"):
    #             filepath = os.path.join(dirpath, filename)
    #             hdr, imgdata = read_fits_file(filepath)
    #             exptime = hdr['EXPTIME']
    #             bkg, bkg_std = calculate_img_bkg(imgdata)
    #             irafsources = detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
    #             if not irafsources:
    #                 continue
    #             _, fwhm, fwhm_std = calculate_fwhm(irafsources)
    #             photometry_result = perform_photometry(irafsources, fwhm, imgdata, bkg=bkg)
    #             fluxes = np.array(photometry_result['flux_fit'])
    #             instr_mags = calculate_magnitudes(photometry_result, exptime)
    #             instr_mags_sigma = calculate_magnitudes_sigma(photometry_result, exptime)
    #             wcs = WCS(hdr)
    #             skypositions = convert_pixel_to_ra_dec(irafsources, wcs)
    #             altazpositions = None
    #             try:
    #                 altazpositions = convert_ra_dec_to_alt_az(skypositions, hdr, lat_key='OBSGEO-B', 
    #                                                                 lon_key= 'OBSGEO-L', elev_key='OBSGEO-H')
    #             except AttributeError as e:
    #                 print(e)
    #                 continue
    #             matched_stars = find_ref_stars(reference_stars, 
    #                                                   ref_star_positions,
    #                                                   skypositions,
    #                                                   instr_mags,
    #                                                   instr_mags_sigma,
    #                                                   fluxes,
    #                                                   ground_based=True,
    #                                                   altazpositions=altazpositions)
    #             if not matched_stars:
    #                 continue
    #             large_table_columns = update_large_table_columns(large_table_columns, 
    #                                                                     matched_stars, 
    #                                                                     hdr, 
    #                                                                     exptime, 
    #                                                                     ground_based=True, 
    #                                                                     name_key='Name')
    
    # large_stars_table = create_large_stars_table(large_table_columns, ground_based=True)
    # large_stars_table.pprint_all()
    # # large_stars_table = remove_large_airmass(large_stars_table)
    # stars_table = group_each_star(large_stars_table, ground_based=True)
    # stars_table.pprint_all()
    
    # instr_filters = ['b', 'v', 'r', 'i']
    # app_mag_table = Table(stars_table['Field', 'Name', 'V_ref', '(B-V)', '(U-B)', '(V-R)', '(V-I)', 'V_sigma'])
    # for instr_filter in instr_filters:
    #     app_mag_table_filter = apply_gb_transforms_VERIFICATION(gb_final_transforms, stars_table, instr_filter)
    #     app_mag_table = hstack([app_mag_table, app_mag_table_filter[instr_filter.upper()]])
    
    # app_mag_table.pprint_all()
    
    # import matplotlib.pyplot as plt
    # # import matplotlib
    # # matplotlib.use('TkAgg')
    # plt.plot(app_mag_table['V_ref'] + app_mag_table['(B-V)'], app_mag_table['B'], 'o')
    # m, b = np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['B'])] + app_mag_table['(B-V)'][~np.isnan(app_mag_table['B'])], app_mag_table['B'][~np.isnan(app_mag_table['B'])], 1)
    # plt.plot(app_mag_table['V_ref'] + app_mag_table['(B-V)'], m*(app_mag_table['V_ref'] + app_mag_table['(B-V)'])+b, '-', label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'] + app_mag_table['(B-V)'], app_mag_table['V_ref'] + app_mag_table['(B-V)'], '-', label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('B (calculated)')
    # plt.xlabel('B (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()
    
    # plt.plot(app_mag_table['V_ref'], app_mag_table['V'], 'o')
    # m, b = np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['V'])], app_mag_table['V'][~np.isnan(app_mag_table['V'])], 1)
    # plt.plot(app_mag_table['V_ref'], m*(app_mag_table['V_ref'])+b, '-', label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'], app_mag_table['V_ref'], '-', label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('V (calculated)')
    # plt.xlabel('V (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()
    
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-R)'], app_mag_table['R'], 'o')
    # m, b = np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['R'])] - app_mag_table['(V-R)'][~np.isnan(app_mag_table['R'])], 
    #                   app_mag_table['R'][~np.isnan(app_mag_table['R'])], 1)
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-R)'], 
    #           m*(app_mag_table['V_ref'] - app_mag_table['(V-R)'])+b, 
    #           '-', label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-R)'], app_mag_table['V_ref'] - app_mag_table['(V-R)'], '-', label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('R (calculated)')
    # plt.xlabel('R (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()
    
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-I)'], app_mag_table['I'], 'o')
    # m, b = np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['I'])] - app_mag_table['(V-I)'][~np.isnan(app_mag_table['I'])], 
    #                   app_mag_table['I'][~np.isnan(app_mag_table['I'])], 1)
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-I)'], 
    #           m*(app_mag_table['V_ref'] - app_mag_table['(V-I)'])+b, 
    #           '-', label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-I)'], app_mag_table['V_ref'] - app_mag_table['(V-I)'], '-', label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('I (calculated)')
    # plt.xlabel('I (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()
    return gb_final_transforms


def _main_sb_transform_calc(directory, 
                            ref_stars_file, 
                            plot_results=False, 
                            save_plots=False,
                            file_suffix=".fits", 
                            exposure_key='EXPTIME',  
                            name_key='Name',
                            transform_index_list=['(B-V)', '(V-R)', '(V-I)'],
                            **kwargs):
    # TODO: Docstring.
    reference_stars, ref_star_positions = read_ref_stars(ref_stars_file)
    large_table_columns = init_large_table_columns()
    
    if save_plots:
        save_loc = kwargs.get('save_loc')
        unique_id = kwargs.get('unique_id')
    
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(file_suffix):
                filepath = os.path.join(dirpath, filename)
                hdr, imgdata = read_fits_file(filepath)
                exptime = hdr[exposure_key]
                bkg, bkg_std = calculate_img_bkg(imgdata)
                irafsources = detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
                if not irafsources:
                    continue
                _, fwhm, fwhm_std = calculate_fwhm(irafsources)
                photometry_result = perform_photometry(irafsources, fwhm, imgdata, bkg=bkg)
                fluxes = np.array(photometry_result['flux_fit'])
                instr_mags = calculate_magnitudes(photometry_result, exptime)
                instr_mags_sigma = calculate_magnitudes_sigma(photometry_result, exptime)
                wcs = WCS(hdr)
                skypositions = convert_pixel_to_ra_dec(irafsources, wcs)
                matched_stars = find_ref_stars(reference_stars, 
                                                     ref_star_positions,
                                                     skypositions,
                                                     instr_mags,
                                                     instr_mags_sigma,
                                                     fluxes,
                                                     ground_based=False,
                                                     altazpositions=None)
                if not matched_stars:
                    continue
                                
                large_table_columns = update_large_table_columns(large_table_columns, 
                                                                 matched_stars, 
                                                                 hdr, 
                                                                 exptime, 
                                                                 ground_based=False, 
                                                                 name_key=name_key)
    large_stars_table = create_large_stars_table(large_table_columns, ground_based=False)
    stars_table = group_each_star(large_stars_table, ground_based=False)
    sb_final_transform_columns = init_sb_final_transform_columns()
    if save_plots:
        write_table_to_latex(stars_table, f"{os.path.join(save_loc, f'{unique_id}_stars_table')}.txt", 
                             formats={'c': '%0.3f',
                                      'c_sigma': '%0.3f'})
        for index in transform_index_list:
            filter_fci, filter_fci_sigma, zprime_fci, zprime_fci_sigma = space_based_transform(stars_table, 
                                                                                               plot_results=plot_results, 
                                                                                               index=index,
                                                                                               save_plots=save_plots, 
                                                                                               save_loc=save_loc,
                                                                                               unique_id=unique_id)
            sb_final_transform_columns = update_sb_final_transform_columns(sb_final_transform_columns, 
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
            sb_final_transform_columns = update_sb_final_transform_columns(sb_final_transform_columns, 
                                                                           index, 
                                                                           filter_fci, 
                                                                           filter_fci_sigma, 
                                                                           zprime_fci, 
                                                                           zprime_fci_sigma)
            # print(f"(V-clear) = ({filter_fci:.3f} +/- {filter_fci_sigma:.3f}) * {index} + " \
            #       f"({zprime_fci:.3f} +/- {zprime_fci_sigma:.3f})")
    sb_final_transform_table = create_sb_final_transform_table(sb_final_transform_columns)
    if save_plots:
        formats = {
        'T_fCI': '%0.3f',
        'T_fCI_sigma': '%0.3f',
        'Z_fCI': '%0.3f',
        'Z_fCI_sigma': '%0.3f'
        }
        write_table_to_latex(sb_final_transform_table, f"{os.path.join(save_loc, f'{unique_id}_transform_table')}.txt",
                             formats=formats)
    return sb_final_transform_table

def _main_sc_lightcurve(directory, temp_dir='tmp', max_distance_from_sat=20, size=25, max_num_nan=5, plot_results=0):
    filecount, filenames = copy_and_rename(directory=directory, temp_dir=temp_dir, debugging=True)
    set_sat_positions_bool = True
    change_sat_positions_bool = False
    num_nan = 0
    for filenum, file in enumerate(filenames):
        filepath = f"{temp_dir}/{file}"
        hdr, imgdata = read_fits_file(filepath)
        if set_sat_positions_bool:
            set_sat_positions_bool, sat_information = set_sat_positions(imgdata, filecount, set_sat_positions_bool)
        sat_information = add_new_time_and_filter(hdr, sat_information, filenum)
        if change_sat_positions_bool:
            change_sat_positions_bool, sat_information = change_sat_positions(filenames, 
                                                                              filenum, 
                                                                              num_nan, 
                                                                              sat_information,
                                                                              change_sat_positions_bool,
                                                                              max_distance_from_sat=max_distance_from_sat, 
                                                                              size=size, 
                                                                              temp_dir=temp_dir, 
                                                                              cmap_set='Set1', 
                                                                              plot_results=plot_results)
        exptime = hdr['EXPTIME'] * u.s
        bkg, bkg_std = calculate_img_bkg(imgdata)
        irafsources = detecting_stars(imgdata, bkg, bkg_std)
        if not irafsources:
            sat_information.num_nans[:] = 0
            continue
        plot_detected_sats(filenames[filenum],
                           plot_results, 
                           imgdata, 
                           irafsources, 
                           sat_information, 
                           max_distance_from_sat=max_distance_from_sat, 
                           norm=LogNorm())
        fwhms, fwhm, fwhm_std = calculate_fwhm(irafsources)
        fwhms_arcsec, fwhm_arcsec, fwhm_std_arcsec = convert_fwhm_to_arcsec(hdr, fwhms, fwhm, fwhm_std)
        photometry_result = perform_photometry(irafsources, fwhm, imgdata, bkg)
        instr_mags = calculate_magnitudes(photometry_result, exptime)
        instr_mags_sigma = calculate_magnitudes_sigma(photometry_result, exptime)
        check_if_sat(sat_information, 
                     filenum, 
                     irafsources, 
                     instr_mags, 
                     instr_mags_sigma, 
                     fwhms_arcsec,
                     max_distance_from_sat=max_distance_from_sat)
        change_sat_positions_bool, num_nan = determine_if_change_sat_positions(sat_information, 
                                                                               filenum, 
                                                                               change_sat_positions_bool, 
                                                                               max_num_nan=max_num_nan)
    remove_temp_dir(temp_dir=temp_dir)
    sats_table = sat_information.sats_table
    uncertainty_table = sat_information.uncertainty_table
    sat_fwhm_table = sat_information.sat_fwhm_table
    return sats_table, uncertainty_table, sat_fwhm_table
        

def __debugging__(gb_final_transforms, save_loc):
    """
    Debug the main functions. This should simplify git commits by only needing to edit this file.

    Returns
    -------
    None.

    """
    # Maybe put the write table option in here?
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
    write_table_to_latex(gb_final_transforms, f"{os.path.join(save_loc, 'gb_final_transforms')}.txt", formats=formats)
    return


# TODO: change the SCInstrMagLightCurve.py code into functions in this file.