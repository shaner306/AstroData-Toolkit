# -*- coding: utf-8 -*-
"""
AstroFunctions.py.

This file holds all of the functions that will likely be implemented in a 
final version of the image processor.

Created on Thu Apr 15 10:14:43 2021

@author: Jack Wawrow
"""
import ctypes
import math
import os
import re
import tkinter as tk
from collections import namedtuple, Counter
from itertools import permutations, groupby, combinations
from math import sqrt, atan, pi
from shutil import copy2, rmtree
from warnings import warn

import astropy.units as u
import cv2 as cv
import matplotlib.cm as cm
import matplotlib.dates as mdates
import numpy
import numpy as np
import pandas as pd
import scipy
import win32com
from astropy import table
from astropy.coordinates import EarthLocation, AltAz, SkyCoord,\
    match_coordinates_sky
from astropy.io import fits, ascii
from astropy.modeling.fitting import LevMarLSQFitter, FittingWithOutlierRemoval
from astropy.modeling.models import Linear1D
from astropy.stats import SigmaClip
from astropy.stats import sigma_clip,\
    sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.table import Table, QTable, hstack
from astropy.time import Time
from astropy.utils import iers
from astropy.wcs import WCS
import matplotlib
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,\
    NavigationToolbar2Tk
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
from photutils.aperture import RectangularAperture
from photutils.background import Background2D
from photutils.background import MeanBackground
from photutils.background import MedianBackground
from photutils.background import SExtractorBackground
from photutils.detection import IRAFStarFinder
from photutils.psf import DAOGroup, BasicPSFPhotometry, IntegratedGaussianPRF
from skimage import measure
from tqdm import tqdm
from random import shuffle, choice

# from scipy.optimize import curve_fit

# If you are working offline, set the following to False.
# In that case, it will use the International Earth Rotation and Reference
# Systems (IERS) data
# that was included with the astropy release installed on the system.
# That may not be as accurate
# for times that are more recent thatn the astropy release date.
iers.conf.auto_download = True

# Testing Github Desktop functionality.


plt.rcParams.update({"text.usetex": False})


# def linear_func(x, m, b):
#     y = (m * x) + b
#     return y


def init_linear_fitting(niter=3, sigma=3.0, slope=1.0, intercept=0.0):
    """
    Initilize the parameters needed to make a linear fit.

    Parameters
    ----------
    niter : int, optional
        Maximum number of iterations. The default is 3.
    sigma : float, optional
        The number of standard deviations to use for both the lower and upper 
        clipping limit. The default is 3.0.

    Returns
    -------
    fit : An Astropy fitter
        An instance of an Astropy LevMarLSQFitter.
    or_fit : astropy.modeling.fitting.FittingWithOutlierRemoval
        An instance of an FittingWithOutlierRemoval using sigma_clip as the 
        outlier remover.
    line_init : astropy.modeling.Fittable1DModel
        An instance of an Astropy Linear1D.

    """
    fit = LevMarLSQFitter(calc_uncertainties=True)
    or_fit = FittingWithOutlierRemoval(
        fit, sigma_clip, niter=niter, sigma=sigma)
    line_init = Linear1D(slope=slope, intercept=intercept)
    return fit, or_fit, line_init


def BackgroundIteration(image, tolerance):
    old_mean = 1e9
    old_rms = 1e9

    new_mean = 2e9
    new_rms = 2e9

    while abs(new_rms - old_rms) > (tolerance * old_rms):
        # old_mean = float(new_mean)
        old_rms = float(new_rms)
        # image = myclip(image, (old_mean - 2 * old_rms),
        # (old_mean + 2 * old_rms))

        if (np.size(image) == 0):
            new_mean = 0
            new_rms = 2e9
            break

        new_mean = np.mean(image)
        new_rms = numpy.std(image)
        # retval = [new_mean, new_rms]
        return new_mean, new_rms


def myclip(x1, lo, hi):
    vector = np.vectorize(np.float)
    x = vector(x1)

    float(hi)
    float(lo)

    y = (x * np.any(x <= hi)) + ((hi) * np.any(x > hi))
    y = (y * np.any(y > lo)) + ((lo) * np.any(y <= lo))
    return y


def PointSourceFluxExtraction(mask_x, mask_y, flux_image):
    num_elem_x = mask_x.size
    # num_elem_y = mask_y.size
    sum1 = 0
    pix_flux = np.zeros((num_elem_x))
    for i in range(num_elem_x):
        x_pix = mask_x[i]
        y_pix = mask_y[i]

        sum1 = sum1 + flux_image[x_pix, y_pix]
        pix_flux[i] = flux_image[x_pix, y_pix]

    object_flux = sum1
    max_pixel_flux = max(pix_flux)
    return object_flux, max_pixel_flux


def MomentCalculation(xmask, ymask, xc, yc, p, q):
    num_pix = xmask.size
    mom = sum((xmask - xc) ** p * (ymask - yc) ** q) / num_pix
    moment = mom
    return moment


def EccentricityCalculation(m11, m02, m20):
    eccent = numpy.sqrt((m20 - m02) ** 2 + (4 * m11 ** 2)) / (m20 + m02)
    return eccent


def Compact(num_pix, m02, m20):
    compact = (num_pix / (m02 + m20))
    return compact


def WeightedCentroid(mask_x, mask_y, flux_image):
    num_elem_x = mask_x.size
    num_elem_y = mask_y.size
    x_wt_sum = 0
    y_wt_sum = 0
    flux_sum = 0
    # print("2")
    if num_elem_x != num_elem_y:
        # object_flux = -999
        # print("3")
        return
    else:
        for i in range(num_elem_x):
            x_pix = mask_x[i]
            y_pix = mask_y[i]

            x_wt_sum = x_wt_sum + (x_pix * flux_image[x_pix, y_pix])
            y_wt_sum = y_wt_sum + (y_pix * flux_image[x_pix, y_pix])
            flux_sum = flux_sum + flux_image[x_pix, y_pix]

    x_centroid = x_wt_sum / flux_sum
    y_centroid = y_wt_sum / flux_sum

    x_var_sum = 0
    y_var_sum = 0
    flux_sum = 0
    # print("2")
    for i in range(num_elem_x):
        x_pix = mask_x[i]
        y_pix = mask_y[i]

        x_var_sum = x_var_sum + ((x_pix - x_centroid)
                                 ** 2 * flux_image[x_pix, y_pix])
        y_var_sum = y_var_sum + ((y_pix - y_centroid)
                                 ** 2 * flux_image[x_pix, y_pix])
        flux_sum = flux_sum + flux_image[x_pix, y_pix]

    x_rms = numpy.sqrt(x_var_sum / flux_sum)
    y_rms = numpy.sqrt(y_var_sum / flux_sum)
    return x_centroid, x_rms, y_centroid, y_rms


def read_ref_stars(ref_stars_file):
    """
    Read a file containing information regarding the reference stars to be
    used to calculate the transforms.

    Parameters
    ----------
    ref_stars_file : string
        Location of the reference stars file.

    Returns
    -------
    reference_stars : astropy.table.Table
        Table with the data extracted from ref_stars_file.
    ref_star_positions : astropy.coordinates.sky_coordinate.SkyCoord
        AstroPy SkyCoord object containing the RA/dec positions of all
        reference stars in the file.

    """
    try:
        reference_stars = ascii.read(
            ref_stars_file, format='basic',
            delimiter='\t', guess=False, encoding='UTF-8')
    except Exception:
        reference_stars = ascii.read(ref_stars_file, encoding='UTF-8')
    reference_stars = reference_stars.filled(np.nan)
    ref_star_positions = SkyCoord(
        ra=reference_stars['RA'],
        dec=reference_stars['Dec'], unit=(u.hourangle, u.deg))
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
    with fits.open(filepath, memmap=False) as hdul:
        hdr = hdul[0].header
        imgdata = hdul[0].data
    # hdul = fits.open(filepath)
    # hdr = hdul[0].header
    # imgdata = hdul[0].data
    # hdul.close()
    return hdr, imgdata


def get_instr_filter_name(hdr, filter_key='FILTER'):
    """
    Get the name of the filter used when taking the image.

    Parameters
    ----------
    hdr : astropy.io.fits.header.Header
        Header from the fits file.
    filter_key : string, optional
        Keyword of the entry in the FITS header that 
        contains the filter information. The default is 'FILTER'.

    Returns
    -------
    instr_filter : string
        Instrumental filter band of the image.

    """
    instr_filter = hdr[filter_key][0].lower()
    return instr_filter


def calculate_img_bkg(imgdata, sigma=3.0):
    """
    Calculate the median and standard deviation of the background of the sigma
    clipped image.

    Parameters
    ----------
    imgdata : numpy.ndarray
        Data from the fits file.
    sigma : float
        Number of standard deviations to use for both the lower and
        upper clipping limit. The default is 3.0.

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
        Number of standard deviations above the background to use as the
        threshold.

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
    iraffind = IRAFStarFinder(
        threshold=sigma * bkg_std, fwhm=fwhm, brightest=150)
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
        AstroPy SkyCoord object containing the RA/dec positions of all 
        sources in the image.

    """
    skypositions = wcs.pixel_to_world(
        irafsources['xcentroid'], irafsources['ycentroid'])
    return skypositions


def convert_ra_dec_to_alt_az(skypositions, hdr, lat_key='SITELAT',
                             lon_key='SITELONG', elev_key='SITEELEV'):
    # TODO: Reduce Calculations if Headers have ALTAZ
    """
    Convert RA/dec locations to Altitude/Azimuth (Azimuth/Elevation).

    Parameters
    ----------
    skypositions : astropy.coordinates.sky_coordinate.SkyCoord
        AstroPy SkyCoord object containing the RA/dec position(s) to be 
        converted to Alt/Az.
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
    Calculate the mean and standard deviation FWHM of all sources in the image
    in pixels.

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
    fwhms : numpy array
        Array containing the FWHM of all sources in the image (in pixels).
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
                           focal_length_key='FOCALLEN', xpixsz_key='XPIXSZ',
                           ypixsz_key='YPIXSZ'):
    """
    Convert the FWHM from pixels to arcsec.

    Parameters
    ----------
    hdr : astropy.io.fits.header.Header
        Header from the fits file.
    fwhms : numpy array
        Array containing the FWHM of all sources in the image (in pixels).
    fwhm : float
        Mean FWHM of all sources in the image.
    fwhm_std : float
        Standard deviation of the FWHM of all sources in the image.
    focal_length_key : string, optional
        Key in the FITS header that contains the focal length. The default is 
        'FOCALLEN'.
    xpixsz_key : string, optional
        Key in the FITS header that contains the x pixel size. The default is
        'XPIXSZ'.
    ypixsz_key : string, optional
        Key in the FITS header that contains the y pixel size. The default is 
        'YPIXSZ'.

    Returns
    -------
    fwhms_arcsec : numpy array
        Array containing the FWHM of all sources in the image in arcsec.
    fwhm_arcsec : float
        Mean FWHM of all sources in the image in arcsec.
    fwhm_std_arcsec : float
        Standard deviation of the FWHM of all sources in the image in arcsec.

    """
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
        print("FWHM in pixels.")
        fwhms_arcsec = fwhms
        fwhm_arcsec = None
        fwhm_std_arcsec = None
    return fwhms_arcsec, fwhm_arcsec, fwhm_std_arcsec


def convert_fwhm_to_arcsec_trm(hdr, fwhm,
                               focal_length_key='FOCALLEN',
                               xpixsz_key='XPIXSZ', ypixsz_key='YPIXSZ'):
    """
    Convert the FWHM from pixels to arcsec.

    Parameters
    ----------
    hdr : astropy.io.fits.header.Header
        Header from the fits file.
    fwhm : float
        Mean FWHM of all sources in the image.
    focal_length_key : string, optional
        Key in the FITS header that contains the focal length.
        The default is 'FOCALLEN'.
    xpixsz_key : string, optional
        Key in the FITS header that contains the x pixel size.
        The default is 'XPIXSZ'.
    ypixsz_key : string, optional
        Key in the FITS header that contains the y pixel size.
        The default is 'YPIXSZ'.

    Returns
    -------
    fwhm_arcsec : float
        Mean FWHM of all sources in the image in arcsec.

    """
    try:
        focal_length = hdr[focal_length_key] * u.mm
        xpixsz = hdr[xpixsz_key]
        ypixsz = hdr[ypixsz_key]
        if xpixsz == ypixsz:
            pixsz = xpixsz * u.um
            rad_per_pix = atan(pixsz / focal_length) * u.rad
            arcsec_per_pix = rad_per_pix.to(u.arcsec)
            fwhm_arcsec = fwhm * arcsec_per_pix.value
    except KeyError:
        print("FWHM in pixels.")
        fwhm_arcsec = fwhm
    return fwhm_arcsec


def perform_photometry(irafsources, fwhm, imgdata, bkg,
                       fitter=LevMarLSQFitter(), fitshape=5):
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
    bkg_estimator : callable, instance of any
    photutils.background.BackgroundBase subclass, optional
        bkg_estimator should be able to compute either a scalar background or
        a 2D background of a given 2D image.
        If None, no background subtraction is performed. The default is None.
    fitter : astropy.modeling.fitting.Fitter instance, optional
        Fitter object used to compute the optimized centroid positions and/or
        flux of the identified sources.
        The default is LevMarLSQFitter().
    fitshape : int or length-2 array-like, optional
        Rectangular shape around the center of a star which will be used to
        collect the data to do the fitting.
        Can be an integer to be the same along both axes. For example, 5 is
        the same as (5, 5),
        which means to fit only at the following relative pixel positions:
            [-2, -1, 0, 1, 2].
        Each element of fitshape must be an odd number. The default is 25.

    Returns
    -------
    photometry_result : astropy.table.Table or None
        Table with the photometry results, i.e., centroids and fluxes
        estimations and the initial estimates used to
        start the fitting process. Uncertainties on the fitted parameters
        are reported as columns called
        <paramname>_unc provided that the fitter object contains a dictionary
        called fit_info with the key param_cov,
        which contains the covariance matrix. If param_cov is not present,
        uncertanties are not reported.

    """
    daogroup = DAOGroup(2 * fwhm)
    psf_model = IntegratedGaussianPRF(sigma=fwhm * gaussian_fwhm_to_sigma)
    psf_model.x_0.fixed = True
    psf_model.y_0.fixed = True
    pos = Table(names=['x_0', 'y_0', 'flux_0'],
                data=[irafsources['xcentroid'],
                      irafsources['ycentroid'],
                      irafsources['flux']])

    photometry = BasicPSFPhotometry(group_maker=daogroup,
                                    bkg_estimator=None,
                                    psf_model=psf_model,
                                    fitter=fitter,
                                    fitshape=fitshape)
    photometry_result = photometry(image=imgdata - bkg, init_guesses=pos)
    return photometry_result


def perform_photometry_sat(sat_x, sat_y, fwhm, imgdata, bkg,
                           fitter=LevMarLSQFitter(), fitshape=5):
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
    bkg_estimator : callable, instance of any
    photutils.background.BackgroundBase subclass, optional
        bkg_estimator should be able to compute either a scalar background or a
        2D background of a given 2D image.
        If None, no background subtraction is performed. The default is None.
    fitter : astropy.modeling.fitting.Fitter instance, optional
        Fitter object used to compute the optimized centroid positions and/or
        flux of the identified sources.
        The default is LevMarLSQFitter().
    fitshape : int or length-2 array-like, optional
        Rectangular shape around the center of a star which will be used to
        collect the data to do the fitting.
        Can be an integer to be the same along both axes. For example, 5 is
        the same as (5, 5),
        which means to fit only at the following relative pixel positions:
            [-2, -1, 0, 1, 2].
        Each element of fitshape must be an odd number. The default is 25.

    Returns
    -------
    photometry_result : astropy.table.Table or None
        Table with the photometry results, i.e., centroids and fluxes
        estimations and the initial estimates used to
        start the fitting process. Uncertainties on the fitted parameters are
        reported as columns called
        <paramname>_unc provided that the fitter object contains a dictionary
        called fit_info with the key param_cov,
        which contains the covariance matrix. If param_cov is not present,
        uncertanties are not reported.

    """
    daogroup = DAOGroup(2 * fwhm)  # The 2 is critical seperation
    psf_model = IntegratedGaussianPRF(sigma=fwhm * gaussian_fwhm_to_sigma)
    psf_model.x_0.fixed = True
    psf_model.y_0.fixed = True
    pos = Table(names=['x_0', 'y_0'],
                data=[sat_x, sat_y])

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
        1D array-like containing the fluxes for all sources in the image
        in counts.
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
        Table with the photometry results, i.e.,
        centroids and fluxes estimations and the initial estimates used to
        start the fitting process.
    exptime : float
        Expsure time of the image in seconds.

    Returns
    -------
    instr_mags : numpy array
        Array containing the instrumental magnitudes of the sources in the
        image.

    """
    fluxes = normalize_flux_by_time(photometry_result['flux_fit'], exptime)
    instr_mags_units = u.Magnitude(fluxes)
    instr_mags = instr_mags_units.value
    return instr_mags


def calculate_magnitudes_sigma(photometry_result, exptime):
    """
    Calculate the standard deviation of the instrumental magnitudes of the
    sources in the image.

    Parameters
    ----------
    photometry_result : astropy.table.Table
        Table with the photometry results, i.e., centroids, fluxes, amd flux
        uncertainty estimations and the initial
        estimates used to start the fitting process.
    exptime : float
        Expsure time of the image in seconds.

    Returns
    -------
    instr_mags_sigma : numpy array
        Array containing the standard deviation of the instrumental magnitudes
        of the sources in the image.

    """
    fluxes = normalize_flux_by_time(photometry_result['flux_fit'], exptime)
    flux_uncs = normalize_flux_by_time(photometry_result['flux_unc'], exptime)
    snr = (fluxes / flux_uncs).value  # FIXME : SNR
    instr_mags_sigma = 1.0857 / np.sqrt(snr)  # FIXME : Instrumental magntiude
    return instr_mags_sigma


def calculate_background_sky_brightness(bkg, hdr, exptime,
                                        gb_final_transforms=None,
                                        focal_length_key='FOCALLEN',
                                        xpixsz_key='XPIXSZ',
                                        ypixsz_key='YPIXSZ',
                                        filter_key='Filter'):
    bkg_flux = normalize_flux_by_time(bkg, exptime)
    # bkg_std_flux = normalize_flux_by_time(bkg_std, exptime)
    try:
        focal_length = hdr[focal_length_key] * u.mm
        xpixsz = hdr[xpixsz_key] * u.um
        ypixsz = hdr[ypixsz_key] * u.um
        x_rad_per_pix = atan(xpixsz / focal_length) * u.rad
        x_arcsec_per_pix = x_rad_per_pix.to(u.arcsec)
        y_rad_per_pix = atan(ypixsz / focal_length) * u.rad
        y_arcsec_per_pix = y_rad_per_pix.to(u.arcsec)
        square_arcsec_per_pix = x_arcsec_per_pix * y_arcsec_per_pix
    except KeyError:
        print("Could not determine arcsec^2 / pix.")
        bsb = np.nan
        return bsb
    bkg_flux_per_area = bkg_flux / square_arcsec_per_pix
    # bkg_flux_std_per_area = bkg_std_flux / square_arcsec_per_pix
    bsb_units = u.Magnitude(bkg_flux_per_area)
    instr_bsb = bsb_units.value
    if not gb_final_transforms:
        if 'ZMAG' in hdr.keys():
            bsb = instr_bsb + hdr['ZMAG']
        else:
            bsb = instr_bsb
    else:
        instr_filter = get_instr_filter_name(hdr, filter_key)
        if instr_filter == 'v' or instr_filter == 'g':
            mask = (gb_final_transforms['filter'] == instr_filter) & (
                gb_final_transforms['CI'] == 'B-V')
        else:
            mask = gb_final_transforms['filter'] == instr_filter
        zpoint = float(gb_final_transforms['Z_f'][mask])
        bsb = instr_bsb + zpoint
    return bsb


def calculate_BSB_sigma(bkg, bkg_std, exptime):
    fluxes = normalize_flux_by_time(bkg, exptime)
    flux_uncs = normalize_flux_by_time(bkg_std, exptime)
    snr = (fluxes / flux_uncs).value  # FIXME: SNR
    instr_mags_sigma = 1.0857 / np.sqrt(snr)  # FIXME: Instrumental magntitude
    return instr_mags_sigma


def calculate_limiting_magnitude(bkg, hdr, exptime,
                                 fwhm, gb_final_transforms=None,
                                 focal_length_key='FOCALLEN',
                                 xpixsz_key='XPIXSZ',
                                 ypixsz_key='YPIXSZ',
                                 filter_key='Filter'):
    bkg_flux = normalize_flux_by_time(bkg, exptime)
    square_fwhm = pi * ((fwhm / 2.0) ** 2)
    bkg_flux_per_square_fwhm = bkg_flux * square_fwhm
    limiting_mag_units = u.Magnitude(bkg_flux_per_square_fwhm)
    limiting_mag = limiting_mag_units.value
    # if not gb_final_transforms:
    #     limiting_mag = limiting_mag_units.value
    # else:
    #     instr_filter = get_instr_filter_name(hdr, filter_key)
    #     if instr_filter == 'v' or instr_filter =='g':
    #         mask = (gb_final_transforms['filter'] == instr_filter) &
    # (gb_final_transforms['CI'] == 'B-V')
    #     else:
    #         mask = gb_final_transforms['filter'] == instr_filter
    #     zpoint = float(gb_final_transforms['Z_f'][mask])
    #     limiting_mag = zpoint + limiting_mag_units.value
    return limiting_mag


# def add_zpoint(sat_auxiliary_table, unique_filters,
#  gb_final_transforms=None):
#     if not gb_final_transforms:
#         return


def find_ref_stars(reference_stars,
                   ref_star_positions,
                   skypositions,
                   instr_mags,
                   instr_mags_sigma,
                   fluxes,
                   ground_based=False,
                   altazpositions=None,
                   max_ref_sep=1.0):
    """
    Match the stars detected in the image to those provided in the
    reference star file.

    Parameters
    ----------
    reference_stars : astropy.table.Table
        Table with the data extracted from ref_stars_file.
    ref_star_positions : astropy.coordinates.sky_coordinate.SkyCoord
        AstroPy SkyCoord object containing the RA/dec positions of all
        reference stars in the file.
    skypositions : astropy.coordinates.sky_coordinate.SkyCoord
        AstroPy SkyCoord object containing the RA/dec positions of all
        sources in the image.
    instr_mags : numpy array
        Array containing the instrumental magnitudes of the sources
        in the image.
    instr_mags_sigma : numpy array
        Array containing the standard deviation of the instrumental
        magnitudes of the sources in the image.
    fluxes : numpy array
        Array containing the non-normalized fluxes of all sources.
    ground_based : bool, optional
        Whether or not the image is from a ground-based sensor.
        The default is False.
    altazpositions : astropy.coordinates.sky_coordinate.SkyCoord, optional
        Alt/Az position(s) of the skyposition(s).
        Must be provided if ground_based is true. The default is None.
    max_ref_sep : float, optional
        Maximum angular separtation to consider an image star to be a
        reference star in arcsec. The default is 10.0.

    Returns
    -------
    matched_stars : namedtuple
        Attributes:
            ref_star_index : array-like or int
                Index of the star(s) in ref_stars_file that correspond to a
                star in the image.
            img_star_index : array-like or int
                Index of the star(s) detected in the image that correspond to
                a star in ref_stars_file.
            ref_star : astropy.table.Table
                Rows from ref_stars_file that correspond to a matched
                reference star.
            ref_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched reference star(s) from ref_stars_file.
            img_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched star(s) detected in the image.
            ang_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the
                image and ref_stars_file.
            img_instr_mag : numpy array
                Array containing the instrumental magnitudes of the matched
                star(s) detected in the image.
            img_instr_mag_sigma : numpy array
                Array containing the standard deviations of the instrumental
                magnitudes of the matched star(s) 
                detected in the image.
            flux : numpy array
                Array containing the non-normalized fluxes of the matched
                star(s) detected in the image.
            img_star_altaz : astropy.coordinates.sky_coordinate.SkyCoord
                Alt/Az of the matched star(s) detected in the image.
                None if ground_based is False.
            img_star_airmass : float
                sec(z) of img_star_altaz. None if ground_based is False.
    None : if no stars are matched.

    """
    # reference_stars, ref_star_positions = read_ref_stars(ref_stars_file)
    max_ref_sep = max_ref_sep * u.arcsec
    idx, sep2d, _ = match_coordinates_sky(skypositions, ref_star_positions)
    img_star_index = np.where(sep2d < max_ref_sep)
    if len(img_star_index[0]) == 0:
        warn("No reference star detected in the image.")
        return
    elif len(img_star_index[0]) == 1:
        img_star_index = int(img_star_index[0])
    else:
        img_star_index = img_star_index[0]
    # ref_star_index = np.where(sep2d < max_ref_sep)
    # if len(ref_star_index[0]) == 0:
    #     print("No reference star detected in the image.")
    #     return
    # elif len(ref_star_index[0]) == 1:
    #     ref_star_index = int(ref_star_index[0])
    # else:
    #     ref_star_index = ref_star_index[0]
    ref_star_index = idx[img_star_index]
    # img_star_index = idx[ref_star_index]
    ref_star = reference_stars[ref_star_index]
    ref_star_loc = ref_star_positions[ref_star_index]
    img_star_loc = skypositions[img_star_index]
    ang_separation = sep2d[img_star_index]
    # ang_separation = sep2d[ref_star_index]
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
                Index of the star(s) in ref_stars_file that correspond to a
                star in the image.
            img_star_index : array-like or int
                Index of the star(s) detected in the image that correspond to
                a star in ref_stars_file.
            ref_star : astropy.table.Table
                Rows from ref_stars_file that correspond to a matched
                reference star.
            ref_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched reference star(s) from ref_stars_file.
            img_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched star(s) detected in the image.
            ang_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the
                image and ref_stars_file.
            img_instr_mag : numpy array
                Array containing the instrumental magnitudes of the matched
                star(s) detected in the image.
            img_instr_mag_sigma : numpy array
                Array containing the standard deviations of the instrumental
                magnitudes of the matched star(s)
                detected in the image.
            flux : numpy array
                Array containing the non-normalized fluxes of the matched
                star(s) detected in the image.
            img_star_altaz : astropy.coordinates.sky_coordinate.SkyCoord
                Alt/Az of the matched star(s) detected in the image.
                None if ground_based is False.
            img_star_airmass : float
                sec(z) of img_star_altaz. None if ground_based is False.
    name_key : string, optional
        The column name of the unique identifier/star name in
        matched_stars.ref_star. The default is 'Name'.

    Returns
    -------
    field : string
        Unique identifier of the star field that the reference star is
        in (e.g. Landolt field "108").

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
        split_string = re.split(
            '[^a-zA-Z0-9]', str(matched_stars.ref_star[name_key]))
        if len(split_string) > 1:
            return ' '.join(split_string[:-1])
        elif len(split_string) == 1:
            return split_string[0]
        else:
            return
    else:
        return


def init_auxiliary_data_columns():
    """
    Initialize the columns that will create the auxiliary data table.

    Returns
    -------
    auxiliary_data_columns : namedtuple
        Attributes:
            filename : empty list
                Name of the file to calculate the transforms for.
            exposure_time : empty list
                Exposure time of the image.
            fwhm : empty list
                Mean FWHM of all sources in the image.
            fwhm_std : empty list
                Standard deviation of the FWHM of all sources in the image.
            avg_mag_sigma : empty list
                Average uncertainty of the magnitudes of all of the reference
                stars in the image.
            std_mag_sigma : empty list
                Standard deviation of the uncertainties of the magnitudes of
                all of the reference stars in the image.

    """
    filename = []
    exposure_time = []
    fwhm = []
    fwhm_std = []
    avg_mag_sigma = []
    std_mag_sigma = []
    auxiliary_data_columns = namedtuple('auxiliary_data_columns',
                                        ['filename',
                                         'exposure_time',
                                         'fwhm',
                                         'fwhm_std',
                                         'avg_mag_sigma',
                                         'std_mag_sigma'])
    return auxiliary_data_columns(filename, exposure_time, fwhm,
                                  fwhm_std, avg_mag_sigma, std_mag_sigma)


def update_auxiliary_data_columns(auxiliary_data_columns, filename,
                                  exptime, fwhm, fwhm_std, matched_stars):
    """
    Update the columns that will create the auxiliary data table.

    Parameters
    ----------
    auxiliary_data_columns : namedtuple
        Attributes:
            filename : string list
                Name of the file to calculate the transforms for.
            exposure_time : float or int list
                Exposure time of the image.
            fwhm : float list
                Mean FWHM of all sources in the image.
            fwhm_std : float list
                Standard deviation of the FWHM of all sources in the image.
            avg_mag_sigma : float list
                Average uncertainty of the magnitudes of all of the reference
                stars in the image.
            std_mag_sigma : float list
                Standard deviation of the uncertainties of the magnitudes of
                all of the reference stars in the image.
    filename : string
        Name of the file to calculate the transforms for.
    exptime : float or int
        Exposure time of the image.
    fwhm : float
        Mean FWHM of all sources in the image.
    fwhm_std : float
        Standard deviation of the FWHM of all sources in the image.
    matched_stars : namedtuple
        Attributes:
            ref_star_index : array-like or int
                Index of the star(s) in ref_stars_file that correspond to a
                star in the image.
            img_star_index : array-like or int
                Index of the star(s) detected in the image that correspond to
                a star in ref_stars_file.
            ref_star : astropy.table.Table
                Rows from ref_stars_file that correspond to a matched
                reference star.
            ref_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched reference star(s) from ref_stars_file.
            img_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched star(s) detected in the image.
            ang_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the
                image and ref_stars_file.
            img_instr_mag : numpy array
                Array containing the instrumental magnitudes of the matched
                star(s) detected in the image.
            img_instr_mag_sigma : numpy array
                Array containing the standard deviations of the instrumental
                magnitudes of the matched star(s) 
                detected in the image.
            flux : numpy array
                Array containing the non-normalized fluxes of the matched
                star(s) detected in the image.
            img_star_altaz : astropy.coordinates.sky_coordinate.SkyCoord
                Alt/Az of the matched star(s) detected in the image.
                None if ground_based is False.
            img_star_airmass : float
                sec(z) of img_star_altaz. None if ground_based is False.

    Returns
    -------
    updated_auxiliary_data_columns : namedtuple
        Attributes:
            filename : string list
                Name of the file to calculate the transforms for.
            exposure_time : float or int list
                Exposure time of the image.
            fwhm : float list
                Mean FWHM of all sources in the image.
            fwhm_std : float list
                Standard deviation of the FWHM of all sources in the image.
            avg_mag_sigma : float list
                Average uncertainty of the magnitudes of all of the reference
                stars in the image.
            std_mag_sigma : float list
                Standard deviation of the uncertainties of the magnitudes of
                all of the reference stars in the image.

    """
    updated_auxiliary_data_columns = auxiliary_data_columns
    updated_auxiliary_data_columns.filename.append(filename)
    updated_auxiliary_data_columns.exposure_time.append(exptime)
    updated_auxiliary_data_columns.fwhm.append(fwhm)
    updated_auxiliary_data_columns.fwhm_std.append(fwhm_std)
    avg_mag_sigma = np.mean(matched_stars.img_instr_mag_sigma)
    std_mag_sigma = np.std(matched_stars.img_instr_mag_sigma)
    updated_auxiliary_data_columns.avg_mag_sigma.append(avg_mag_sigma)
    updated_auxiliary_data_columns.std_mag_sigma.append(std_mag_sigma)
    return updated_auxiliary_data_columns


def create_auxiliary_data_table(auxiliary_data_columns):
    """
    Create the auxiliary data table from the auxiliary data columns.

    Parameters
    ----------
    auxiliary_data_columns : namedtuple
        Attributes:
            filename : string list
                Name of the file to calculate the transforms for.
            exposure_time : float or int list
                Exposure time of the image.
            fwhm : float list
                Mean FWHM of all sources in the image.
            fwhm_std : float list
                Standard deviation of the FWHM of all sources in the image.
            avg_mag_sigma : float list
                Average uncertainty of the magnitudes of all of the reference
                stars in the image.
            std_mag_sigma : float list
                Standard deviation of the uncertainties of the magnitudes of
                all of the reference stars in the image.

    Returns
    -------
    auxiliary_data_table : astropy.table.Table
        Table containing additional information about the images used to
        create the transforms. Has columns:
            filename : string
                Name of the file the transforms were calculated for.
            exptime : float or int
                Exposure time of the image.
            fwhm : float
                Mean FWHM of all sources in the image.
            fwhm_sigma : float
                Standard deviation of the FWHM of all sources in the image.
            avg mag_sigma : float
                Average uncertainty of the magnitudes of all of the reference
                stars in the image.
            std mag_sigma : float
                Standard deviation of the uncertainties of the magnitudes of
                all of the reference stars in the image.

    """
    auxiliary_data_table = Table(
        names=[
            'filename',
            'exptime',
            'fwhm',
            'fwhm_sigma',
            'avg mag_sigma',
            'std mag_sigma'
        ],
        data=[
            auxiliary_data_columns.filename,
            auxiliary_data_columns.exposure_time,
            auxiliary_data_columns.fwhm,
            auxiliary_data_columns.fwhm_std,
            auxiliary_data_columns.avg_mag_sigma,
            auxiliary_data_columns.std_mag_sigma
        ]
    )
    return auxiliary_data_table


def init_large_table_columns():
    """
    Create all of the columns that will be turned into the large stars table.

    Returns
    -------
    large_table_columns : namedtuple
        Attributes:
            field : empty list
                Unique identifier of the star field that the reference star
                is in (e.g. Landolt field "108").
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
                Angular distance between the matched star(s) detected in the
                image and ref_stars_file.
            img_star_mag : empty list
                Instrumental magnitude of the reference star(s) found in
                the image.
            img_star_mag_sigma : empty list
                Standard deviation of the instrumental magnitude of the
                reference star(s) found in the image.
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
                Standard deviation of the apparent V magnitude from the
                reference file.
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
    B_V_sigma_apparents = []
    U_B_sigma_apparents = []
    V_R_sigma_apparents = []
    V_I_sigma_apparents = []
    img_star_airmass = []
    img_names = []
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
                                      'B_V_sigma_apparents',
                                      'U_B_sigma_apparents',
                                      'V_R_sigma_apparents',
                                      'V_I_sigma_apparents',
                                      'img_star_airmass',
                                      'img_name'
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
                               B_V_sigma_apparents,
                               U_B_sigma_apparents,
                               V_R_sigma_apparents,
                               V_I_sigma_apparents,
                               img_star_airmass,
                               img_names
                               # X_rounded
                               )


def update_large_table_columns(large_table_columns, filename,
                               matched_stars, hdr, exptime, ground_based=False,
                               name_key='Name'):
    """
    Update columns to be used for the large stars table based on information
    from the current image.

    Parameters
    ----------
    large_table_columns : namedtuple
        Attributes:
            field : string list
                Unique identifier of the star field that the reference
                star is in (e.g. Landolt field "108").
            ref_star_name : string list
                Name/unique identifier of the reference star.
            times : numpy.float64
                Time of the observation.
            flux_table : numpy.float64
                Flux of the source.
            exposure : numpy.float64
                Exposure time of the image.
            ref_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) from the file in
                unit hourangle.
            ref_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) from the file in
                unit degree.
            img_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) found in the image in
                unit hourangle.
            img_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) found in the image in
                unit degree.
            angular_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in
                the image and ref_stars_file.
            img_star_mag : numpy.float64
                Instrumental magnitude of the reference star(s) found in
                the image.
            img_star_mag_sigma : numpy.float64
                Standard deviation of the instrumental magnitude of the
                reference star(s) found in the image.
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
                Standard deviation of the apparent V magnitude from the
                reference file.
            img_star_airmass : numpy.float64
                Sec(z) of the reference star(s) found in the image.
    matched_stars : namedtuple
        Attributes:
            ref_star_index : array-like or int
                Index of the star(s) in ref_stars_file that correspond 
                to a star in the image.
            img_star_index : array-like or int
                Index of the star(s) detected in the image that correspond 
                to a star in ref_stars_file.
            ref_star : astropy.table.Table
                Rows from ref_stars_file that correspond to a matched 
                reference star.
            ref_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched reference star(s) from ref_stars_file.
            img_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched star(s) detected in the image.
            ang_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in
                the image and ref_stars_file.
            img_instr_mag : numpy array
                Array containing the instrumental magnitudes of the
                matched star(s) detected in the image.
            img_instr_mag_sigma : numpy array
                Array containing the standard deviations of the
                instrumental magnitudes of the matched star(s) 
                detected in the image.
            flux : numpy array
                Array containing the non-normalized fluxes of the
                matched star(s) detected in the image.
            img_star_altaz : astropy.coordinates.sky_coordinate.SkyCoord
                Alt/Az of the matched star(s) detected in the image. 
                None if ground_based is False.
            img_star_airmass : float
                sec(z) of img_star_altaz. None if ground_based is False.
    hdr : astropy.io.fits.header.Header
        Header from the fits file.
    exptime : float
        Expsure time of the image in seconds.
    ground_based : bool, optional
        Whether or not the image is from a ground-based sensor.
        The default is False.
    name_key : string, optional
        The column name of the unique identifier/star name in
        matched_stars.ref_star. The default is 'Name'.

    Returns
    -------
    updated_large_table_columns : namedtuple
        Attributes:
            field : string list
                Unique identifier of the star field that the
                reference star is in (e.g. Landolt field "108").
            ref_star_name : string list
                Name/unique identifier of the reference star.
            times : numpy.float64
                Time of the observation.
            flux_table : numpy.float64
                Flux of the source.
            exposure : numpy.float64
                Exposure time of the image.
            ref_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) from the file in
                unit hourangle.
            ref_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) from the file in
                unit degree.
            img_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) found in
                the image in unit hourangle.
            img_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) found in the image in
                unit degree.
            angular_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in
                the image and ref_stars_file.
            img_star_mag : numpy.float64
                Instrumental magnitude of the reference star(s) found in
                the image.
            img_star_mag_sigma : numpy.float64
                Standard deviation of the instrumental magnitude of the
                reference star(s) found in the image.
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
                Standard deviation of the apparent V magnitude from the
                reference file.
            img_star_airmass : numpy.float64
                Sec(z) of the reference star(s) found in the image.
                Only output if ground_based is True.

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
                updated_large_table_columns.field.append(
                    ' '.join(split_string[:-1]))
            elif len(split_string) == 1:
                updated_large_table_columns.field.append(split_string[0])
            else:
                print('Could not find the name of the field.')
        updated_large_table_columns.ref_star_name.extend(
            matched_stars.ref_star[name_key])
        updated_large_table_columns.flux_table.extend(matched_stars.flux)
        time = Time(hdr['DATE-OBS'], format='fits')
        time_repeat = np.full(num_stars, time.jd)
        updated_large_table_columns.times.extend(time_repeat)
        exposure_repeat = np.full(num_stars, exptime)
        updated_large_table_columns.exposure.extend(exposure_repeat)
        updated_large_table_columns.ref_star_RA.extend(
            matched_stars.ref_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.ref_star_dec.extend(
            matched_stars.ref_star_loc.dec)
        updated_large_table_columns.img_star_RA.extend(
            matched_stars.img_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.img_star_dec.extend(
            matched_stars.img_star_loc.dec)
        updated_large_table_columns.angular_separation.extend(
            matched_stars.ang_separation.to(u.arcsec))
        updated_large_table_columns.img_star_mag.extend(
            matched_stars.img_instr_mag)
        updated_large_table_columns.img_star_mag_sigma.extend(
            matched_stars.img_instr_mag_sigma)
        filter_name_repeat = np.full(num_stars, hdr['FILTER'][0])
        updated_large_table_columns.filters.extend(filter_name_repeat)
        filename_repeat = np.full(num_stars, filename)
        updated_large_table_columns.img_name.extend(filename_repeat)
        updated_large_table_columns.V_apparents.extend(
            matched_stars.ref_star['V_ref'])
        try:
            updated_large_table_columns.B_V_apparents.extend(
                matched_stars.ref_star['(B-V)'])
            updated_large_table_columns.U_B_apparents.extend(
                matched_stars.ref_star['(U-B)'])
            updated_large_table_columns.V_R_apparents.extend(
                matched_stars.ref_star['(V-R)'])
            updated_large_table_columns.V_I_apparents.extend(
                matched_stars.ref_star['(V-I)'])
            updated_large_table_columns.V_sigma_apparents.extend(
                matched_stars.ref_star['V_sigma'])
            updated_large_table_columns.B_V_sigma_apparents.extend(
                matched_stars.ref_star['(B-V)_sigma'])
            updated_large_table_columns.U_B_sigma_apparents.extend(
                matched_stars.ref_star['(U-B)_sigma'])
            updated_large_table_columns.V_R_sigma_apparents.extend(
                matched_stars.ref_star['(V-R)_sigma'])
            updated_large_table_columns.V_I_sigma_apparents.extend(
                matched_stars.ref_star['(V-I)_sigma'])
        except KeyError:
            updated_large_table_columns.B_V_apparents.extend(
                matched_stars.ref_star['B-V'])
            updated_large_table_columns.U_B_apparents.extend(
                matched_stars.ref_star['U-B'])
            updated_large_table_columns.V_R_apparents.extend(
                matched_stars.ref_star['V-R'])
            updated_large_table_columns.V_I_apparents.extend(
                matched_stars.ref_star['V-I'])
            updated_large_table_columns.V_sigma_apparents.extend(
                matched_stars.ref_star['e_V'])
            updated_large_table_columns.B_V_sigma_apparents.extend(
                matched_stars.ref_star['e_B-V'])
            updated_large_table_columns.U_B_sigma_apparents.extend(
                matched_stars.ref_star['e_U-B'])
            updated_large_table_columns.V_R_sigma_apparents.extend(
                matched_stars.ref_star['e_V-R'])
            updated_large_table_columns.V_I_sigma_apparents.extend(
                matched_stars.ref_star['e_V-I'])
        if not ground_based:
            return updated_large_table_columns
        updated_large_table_columns.img_star_airmass.extend(
            matched_stars.img_star_airmass)
        # updated_large_table_columns.X_rounded.extend
        # (round(matched_stars.img_star_airmass, 1))
    elif num_stars == 1:
        split_string = re.split(
            '[^a-zA-Z0-9]', str(matched_stars.ref_star[name_key]))
        if len(split_string) > 1:
            updated_large_table_columns.field.append(
                ' '.join(split_string[:-1]))
        elif len(split_string) == 1:
            updated_large_table_columns.field.append(split_string[0])
        else:
            print('Could not find the name of the field.')
        updated_large_table_columns.ref_star_name.append(
            matched_stars.ref_star[name_key])
        updated_large_table_columns.flux_table.append(matched_stars.flux)
        time = Time(hdr['DATE-OBS'], format='fits')
        updated_large_table_columns.times.append(time.jd)
        updated_large_table_columns.exposure.append(exptime)
        updated_large_table_columns.ref_star_RA.append(
            matched_stars.ref_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.ref_star_dec.append(
            matched_stars.ref_star_loc.dec)
        updated_large_table_columns.img_star_RA.append(
            matched_stars.img_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.img_star_dec.append(
            matched_stars.img_star_loc.dec)
        updated_large_table_columns.angular_separation.append(
            matched_stars.ang_separation.to(u.arcsec))
        updated_large_table_columns.img_star_mag.append(
            matched_stars.img_instr_mag)
        updated_large_table_columns.img_star_mag_sigma.append(
            matched_stars.img_instr_mag_sigma)
        updated_large_table_columns.filters.append(hdr['FILTER'][0])
        updated_large_table_columns.img_name.append(filename)
        updated_large_table_columns.V_apparents.append(
            matched_stars.ref_star['V_ref'])
        try:
            updated_large_table_columns.B_V_apparents.append(
                matched_stars.ref_star['(B-V)'])
            updated_large_table_columns.U_B_apparents.append(
                matched_stars.ref_star['(U-B)'])
            updated_large_table_columns.V_R_apparents.append(
                matched_stars.ref_star['(V-R)'])
            updated_large_table_columns.V_I_apparents.append(
                matched_stars.ref_star['(V-I)'])
            updated_large_table_columns.V_sigma_apparents.append(
                matched_stars.ref_star['V_sigma'])
            updated_large_table_columns.B_V_sigma_apparents.append(
                matched_stars.ref_star['(B-V)_sigma'])
            updated_large_table_columns.U_B_sigma_apparents.append(
                matched_stars.ref_star['(U-B)_sigma'])
            updated_large_table_columns.V_R_sigma_apparents.append(
                matched_stars.ref_star['(V-R)_sigma'])
            updated_large_table_columns.V_I_sigma_apparents.append(
                matched_stars.ref_star['(V-I)_sigma'])
        except KeyError:
            updated_large_table_columns.B_V_apparents.append(
                matched_stars.ref_star['B-V'])
            updated_large_table_columns.U_B_apparents.append(
                matched_stars.ref_star['U-B'])
            updated_large_table_columns.V_R_apparents.append(
                matched_stars.ref_star['V-R'])
            updated_large_table_columns.V_I_apparents.append(
                matched_stars.ref_star['V-I'])
            updated_large_table_columns.V_sigma_apparents.append(
                matched_stars.ref_star['e_V'])
            updated_large_table_columns.B_V_sigma_apparents.append(
                matched_stars.ref_star['e_B-V'])
            updated_large_table_columns.U_B_sigma_apparents.append(
                matched_stars.ref_star['e_U-B'])
            updated_large_table_columns.V_R_sigma_apparents.append(
                matched_stars.ref_star['e_V-R'])
            updated_large_table_columns.V_I_sigma_apparents.append(
                matched_stars.ref_star['e_V-I'])
        if not ground_based:
            return updated_large_table_columns
        updated_large_table_columns.img_star_airmass.append(
            matched_stars.img_star_airmass)
        # updated_large_table_columns.X_rounded.append
        # (round(matched_stars.img_star_airmass, 1))
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
                Unique identifier of the star field that the reference star
                is in (e.g. Landolt field "108").
            ref_star_name : string list
                Name/unique identifier of the reference star.
            times : numpy.float64
                Time of the observation.
            flux_table : numpy.float64
                Flux of the source.
            exposure : numpy.float64
                Exposure time of the image.
            ref_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) from the file in
                unit hourangle.
            ref_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) from the file in
                unit degree.
            img_star_RA : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) found in the image in
                unit hourangle.
            img_star_dec : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) found in the image in
                unit degree.
            angular_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in
                the image and ref_stars_file.
            img_star_mag : numpy.float64
                Instrumental magnitude of the reference star(s) found in
                the image.
            img_star_mag_sigma : numpy.float64
                Standard deviation of the instrumental magnitude of the
                reference star(s) found in the image.
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
                Standard deviation of the apparent V magnitude from the
                reference file.
            img_star_airmass : numpy.float64
                Sec(z) of the reference star(s) found in the image.
    ground_based : bool, optional
        Whether or not the image is from a ground-based sensor.
        The default is False.

    Returns
    -------
    large_stars_table : astropy.table.table.QTable
        Table containing all information on the detected stars. Has columns:
            Field : string
                Unique identifier of the star field that the reference star is
                in (e.g. Landolt field "108").
            Name : string
                Name/unique identifier of the reference star.
            Time (JD) : numpy.float64
                Time of the observation.
            RA_ref : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) from the file in
                unit hourangle.
            dec_ref : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) from the file in
                unit degree.
            RA_img : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) found in the image in
                unit hourangle.
            dec_img : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) found in the image in
                unit degree.
            angular_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in
                the image and ref_stars_file.
            flux : numpy.float64
                Flux of the source.
            exposure : numpy.float64
                Exposure time of the image.
            mag_instrumental : numpy.float64
                Instrumental magnitude of the reference star(s) found in
                the image.
            mag_instrumental_sigma : numpy.float64
                Standard deviation of the instrumental magnitude of the
                reference star(s) found in the image.
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
                Standard deviation of the apparent V magnitude from the
                reference file.
            X : numpy.float64
                Sec(z) of the reference star(s) found in the image.
                Only output if ground_based is True.

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
                large_table_columns.B_V_sigma_apparents,
                large_table_columns.U_B_sigma_apparents,
                large_table_columns.V_R_sigma_apparents,
                large_table_columns.V_I_sigma_apparents,
                large_table_columns.img_star_airmass,
                large_table_columns.img_name
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
                'B-V',
                'U-B',
                'V-R',
                'V-I',
                'V_sigma',
                'e_B-V',
                'e_U-B',
                'e_V-R',
                'e_V-I',
                'X',
                'img_name'
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
                large_table_columns.B_V_sigma_apparents,
                large_table_columns.U_B_sigma_apparents,
                large_table_columns.V_R_sigma_apparents,
                large_table_columns.V_I_sigma_apparents
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
                '(B-V)_sigma',
                '(U-B)_sigma',
                '(V-R)_sigma',
                '(V-I)_sigma',
            ]
        )
    return large_stars_table


# class Delta:
#     def __init__(self, delta):
#         self.last = None
#         self.delta = delta
#         self.key = 1
#     def __call__(self, value):
#         if self.last is not None and abs(self.last - value[1]) > self.delta:
#             # Compare with the last value (`self.last`)
#             # If difference is larger than 20, advance to next project
#             self.key += 1
#         self.last = value[1]  # Remeber the last value.
#         return self.key


def Delta(delta):
    last = None
    key = 1

    def keyfunc(value):
        nonlocal last, key
        if last is not None and abs(last - value) > delta:
            key += 1
        last = value
        return key

    return keyfunc


def group_each_star_GB(large_stars_table, keys='Name'):
    """
    Group all detections of unique reference stars together.

    Parameters
    ----------
    large_stars_table : astropy.table.table.QTable
        Table containing all information on the detected stars. Has columns:
            Field : string
                Unique identifier of the star field that the reference star
                is in (e.g. Landolt field "108").
            Name : string
                Name/unique identifier of the reference star.
            Time (JD) : numpy.float64
                Time of the observation.
            RA_ref : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) from the file in
                unit hourangle.
            dec_ref : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) from the file in
                unit degree.
            RA_img : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) found in
                the image in unit hourangle.
            dec_img : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) found in the image in
                unit degree.
            angular_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in
                the image and ref_stars_file.
            flux : numpy.float64
                Flux of the source.
            exposure : numpy.float64
                Exposure time of the image.
            mag_instrumental : numpy.float64
                Instrumental magnitude of the reference star(s) found in
                the image.
            mag_instrumental_sigma : numpy.float64
                Standard deviation of the instrumental magnitude of the
                reference star(s) found in the image.
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
                Standard deviation of the apparent V magnitude from the
                reference file.
            X : numpy.float64
                Sec(z) of the reference star(s) found in the image.
    ground_based : bool, optional
        Whether or not the image is from a ground-based sensor.
        The default is False.
    keys : string, optional
        Table column(s) to use to group the different reference stars.
        The default is 'Name'.

    Returns
    -------
    stars_table : astropy.table.table.Table
        Table containing the mean of the important information for each star.
        Has columns:
            Field : string
                Unique identifier of the star field that the reference star
                is in (e.g. Landolt field "108").
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
                Standard deviation of the apparent V magnitude from
                the reference file.
            <filter> : numpy.float64
                Mean instrumental magnitude of all detections of the star in
                <filter>. There is a different column for
                each different filter used across the images.
            <filter>_sigma : numpy.float64
                Standard deviation of the instrumental magnitudes of all
                detections of the star in <filter>.
                There is a different column for each different filter
                used across the images.
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

    """

    unique_stars = table.unique(large_stars_table, keys=keys)
    N = 1
    nan_array = np.empty(N)
    nan_array.fill(np.nan)
    apparent_mags_table = Table(
        names=[
            'Field',
            'Name',
            'V_ref',
            'B-V',
            'U-B',
            'V-R',
            'V-I',
            'V_sigma',
            'e_B-V',
            'e_U-B',
            'e_V-R',
            'e_V-I',
        ],
        data=[
            np.empty(N, dtype=object),
            np.empty(N, dtype=object),
            nan_array,
            nan_array,
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
    different_filters = table.unique(large_stars_table, keys='filter')
    different_filter_list = list(different_filters['filter'])
    different_filter_list = [different_filter.lower()
                             for different_filter in different_filter_list]
    different_filter_data = np.empty((N, len(different_filter_list)))
    different_filter_data.fill(np.nan)
    filter_sigma_list = []
    for different_filter in different_filter_list:
        filter_sigma_list.append(f"{different_filter}_sigma")
    different_filter_table = Table(
        data=different_filter_data, names=different_filter_list)
    different_filter_sigma_table = Table(
        data=different_filter_data, names=filter_sigma_list)
    filter_X_list = []
    for different_filter in different_filter_list:
        filter_X_list.append(f"X_{different_filter}")
    filter_X_sigma_list = []
    for different_filter in different_filter_list:
        filter_X_sigma_list.append(f"X_{different_filter}_sigma")
    filter_X_table = Table(data=different_filter_data, names=filter_X_list)
    filter_X_sigma_table = Table(
        data=different_filter_data, names=filter_X_sigma_list)
    all_indices, all_indices_formatted = get_all_indicies_combinations(
        different_filter_list, len(different_filter_list),
        multiple_filters=True)
    data_instr_index_table = [
        nan_array for different_index in all_indices_formatted]
    instr_index_table = Table(
        names=all_indices_formatted, data=data_instr_index_table)
    stars_table = table.hstack([apparent_mags_table,
                                different_filter_table,
                                different_filter_sigma_table,
                                filter_X_table,
                                filter_X_sigma_table,
                                instr_index_table],
                               join_type='exact')
    num_columns = len(stars_table.columns)
    stars_table_nan_array = np.empty(num_columns)
    stars_table_nan_array.fill(np.nan)
    # else:
    #     stars_table = table.hstack([apparent_mags_table,
    #                                 different_filter_table,
    #                                 different_filter_sigma_table],
    #                                join_type='exact')
    # stars_table['Field'] = unique_stars['Field']
    # stars_table['V_ref'] = unique_stars['V_ref']
    # stars_table['(B-V)'] = unique_stars['(B-V)']
    # stars_table['(U-B)'] = unique_stars['(U-B)']
    # stars_table['(V-R)'] = unique_stars['(V-R)']
    # stars_table['(V-I)'] = unique_stars['(V-I)']
    # stars_table['V_sigma'] = unique_stars['V_sigma']
    i = 0
    for star in unique_stars['Name']:
        mask = large_stars_table['Name'] == star
        current_star_table = large_stars_table[mask]
        current_star_table.sort('X')
        # current_star_table.pprint(max_lines=-1, max_width=300)
        # data = np.array(current_star_table['X'])
        # gaps = [y - x for x, y in zip(data[:-1], data[1:])]
        # # have python calculate the standard deviation for the gaps
        # sd = np.std(gaps)

        # # create a list of lists, put the first value of the source data in
        # the first
        # lists = [[data[0]]]
        # for x in data[1:]:
        #     # if the gap from the current item to the previous is more than
        # 1 SD
        #     # Note: the previous item is the last item in the last list
        #     # Note: the '> 1' is the part you'd modify to make it stricter
        # or more relaxed
        #     if (x - lists[-1][-1]) / sd > 3:
        #         # then start a new list
        #         lists.append([])
        #     # add the current item to the last list in the list
        #     lists[-1].append(x)

        # print(lists)
        # if len(current_star_table) > 1:
        for k, g in groupby(np.array(current_star_table['X']), key=Delta(0.1)):
            stars_table.add_row(stars_table_nan_array)
            stars_table['Field'][i] = current_star_table['Field'][0]
            stars_table['V_ref'][i] = current_star_table['V_ref'][0]
            stars_table['B-V'][i] = current_star_table['B-V'][0]
            stars_table['U-B'][i] = current_star_table['U-B'][0]
            stars_table['V-R'][i] = current_star_table['V-R'][0]
            stars_table['V-I'][i] = current_star_table['V-I'][0]
            stars_table['V_sigma'][i] = current_star_table['V_sigma'][0]
            stars_table['e_B-V'][i] = current_star_table['e_B-V'][0]
            stars_table['e_U-B'][i] = current_star_table['e_U-B'][0]
            stars_table['e_V-R'][i] = current_star_table['e_V-R'][0]
            stars_table['e_V-I'][i] = current_star_table['e_V-I'][0]
            stars_table['Name'][i] = star
            list_of_airmasses = list(g)
            mask = np.in1d(current_star_table['X'], list_of_airmasses)
            # print(mask)
            current_star_airmass_table = current_star_table[mask]
            # print(current_star_airmass_table)
            # print(np.array(current_star_table['X']))
            # print(list_of_airmasses)
            unique_filters = table.unique(
                current_star_airmass_table, keys='filter')
            for unique_filter in unique_filters['filter']:
                mask = (
                    (current_star_airmass_table['filter'] == unique_filter))
                current_star_filter_table = current_star_airmass_table[mask]
                mags_numpy = np.array(
                    current_star_filter_table['mag_instrumental'])
                mags_std_numpy = np.array(
                    current_star_filter_table['mag_instrumental_sigma'])
                mags_weights = 1 / (mags_std_numpy ** 2)
                # mean_mag = mags_numpy.mean()
                masked_mags = sigma_clip(mags_numpy)
                sig_clip_mask = np.ma.getmask(masked_mags)
                # print(sig_clip_mask)
                # if len(mags_numpy) != len(masked_mags):
                #     print(mags_numpy)
                #     print(masked_mags)
                if np.all(sig_clip_mask):
                    continue
                mean_mag = np.average(
                    mags_numpy[~sig_clip_mask],
                    weights=mags_weights[~sig_clip_mask])
                std_mag = mags_numpy[~sig_clip_mask].std()
                # mean_mag = np.average(mags_numpy, weights=mags_weights)
                # std_mag = mags_numpy.std()
                filter_column = unique_filter.lower()
                sigma_column = f'{filter_column}_sigma'
                stars_table[filter_column][i] = mean_mag
                stars_table[sigma_column][i] = std_mag
                X_numpy = np.array(current_star_filter_table['X'])
                mean_X = X_numpy[~sig_clip_mask].mean()
                std_X = X_numpy[~sig_clip_mask].std()
                # mean_X = X_numpy.mean()
                # std_X = X_numpy.std()
                X_column = f'X_{filter_column}'
                X_std_column = f'X_{filter_column}_sigma'
                stars_table[X_column][i] = mean_X
                stars_table[X_std_column][i] = std_X

            # for i, unique_star in enumerate(multiple_stars):
            # star_mask = stars_table['Name'] == unique_star
            # current_star = stars_table[star_mask]
            for instr_index_name in all_indices_formatted:
                first_mag = stars_table[instr_index_name[0]][i]
                second_mag = stars_table[instr_index_name[-1]][i]
                instr_index = first_mag - second_mag
                stars_table[instr_index_name][i] = instr_index
            i += 1
    stars_table.remove_row(-1)
    # print(different_filter_list)
    return stars_table, different_filter_list


def create_reformatted_large_table(large_stars_table, keys='Name'):
    unique_stars = table.unique(large_stars_table, keys=keys)
    N = len(large_stars_table)
    nan_array = np.empty(N)
    nan_array.fill(np.nan)
    apparent_mags_table = Table(
        names=[
            'Field',
            'Name',
            'V_ref',
            'B-V',
            'U-B',
            'V-R',
            'V-I',
            'V_sigma',
            'e_B-V',
            'e_U-B',
            'e_V-R',
            'e_V-I',
        ],
        data=[
            np.empty(N, dtype=object),
            np.empty(N, dtype=object),
            nan_array,
            nan_array,
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
    different_filters = table.unique(large_stars_table, keys='filter')
    different_filter_list = list(different_filters['filter'])
    different_filter_list = [different_filter.lower()
                             for different_filter in different_filter_list]
    different_filter_data = np.empty((N, len(different_filter_list)))
    different_filter_data.fill(np.nan)
    filter_sigma_list = []
    for different_filter in different_filter_list:
        filter_sigma_list.append(f"{different_filter}_sigma")
    different_filter_table = Table(
        data=different_filter_data, names=different_filter_list)
    different_filter_sigma_table = Table(
        data=different_filter_data, names=filter_sigma_list)
    filter_X_list = []
    for different_filter in different_filter_list:
        filter_X_list.append(f"X_{different_filter}")
    filter_X_sigma_list = []
    for different_filter in different_filter_list:
        filter_X_sigma_list.append(f"X_{different_filter}_sigma")
    filter_X_table = Table(data=different_filter_data, names=filter_X_list)
    filter_X_sigma_table = Table(
        data=different_filter_data, names=filter_X_sigma_list)
    all_indices, all_indices_formatted = get_all_indicies_combinations(
        different_filter_list, len(different_filter_list),
        multiple_filters=True)
    data_instr_index_table = [
        nan_array for different_index in all_indices_formatted]
    instr_index_table = Table(
        names=all_indices_formatted, data=data_instr_index_table)
    reformatted_large_stars_table = table.hstack([apparent_mags_table,
                                                  different_filter_table,
                                                  different_filter_sigma_table,
                                                  filter_X_table,
                                                  filter_X_sigma_table,
                                                  instr_index_table],
                                                 join_type='exact')
    num_columns = len(reformatted_large_stars_table.columns)
    stars_table_nan_array = np.empty(num_columns)
    stars_table_nan_array.fill(np.nan)
    # else:
    #     stars_table = table.hstack([apparent_mags_table,
    #                                 different_filter_table,
    #                                 different_filter_sigma_table],
    #                                join_type='exact')
    # stars_table['Field'] = unique_stars['Field']
    # stars_table['V_ref'] = unique_stars['V_ref']
    # stars_table['(B-V)'] = unique_stars['(B-V)']
    # stars_table['(U-B)'] = unique_stars['(U-B)']
    # stars_table['(V-R)'] = unique_stars['(V-R)']
    # stars_table['(V-I)'] = unique_stars['(V-I)']
    # stars_table['V_sigma'] = unique_stars['V_sigma']
    i = 0
    for current_star_table in large_stars_table:
        # mask = large_stars_table['Name'] == star
        # current_star_table = large_stars_table[mask]
        reformatted_large_stars_table['Field'][i] = current_star_table['Field']
        reformatted_large_stars_table['V_ref'][i] = current_star_table['V_ref']
        reformatted_large_stars_table['B-V'][i] = current_star_table['B-V']
        reformatted_large_stars_table['U-B'][i] = current_star_table['U-B']
        reformatted_large_stars_table['V-R'][i] = current_star_table['V-R']
        reformatted_large_stars_table['V-I'][i] = current_star_table['V-I']
        reformatted_large_stars_table['V_sigma'][i] =\
            current_star_table['V_sigma']
        reformatted_large_stars_table['e_B-V'][i] = current_star_table['e_B-V']
        reformatted_large_stars_table['e_U-B'][i] = current_star_table['e_U-B']
        reformatted_large_stars_table['e_V-R'][i] = current_star_table['e_V-R']
        reformatted_large_stars_table['e_V-I'][i] = current_star_table['e_V-I']
        reformatted_large_stars_table['Name'][i] = current_star_table['Name']
        unique_filter = current_star_table['filter']
        # mags_numpy = np.array(current_star_table['mag_instrumental'])
        # mags_std_numpy = np.array(current_star_table\
        # ['mag_instrumental_sigma'])
        # mags_weights = 1 / (mags_std_numpy ** 2)
        # # mean_mag = mags_numpy.mean()
        # masked_mags = sigma_clip(mags_numpy)
        # sig_clip_mask = np.ma.getmask(masked_mags)
        # # print(sig_clip_mask)
        # # if len(mags_numpy) != len(masked_mags):
        # #     print(mags_numpy)
        # #     print(masked_mags)
        # if np.all(sig_clip_mask):
        #     continue
        # mean_mag = np.average(mags_numpy[~sig_clip_mask],
        # weights=mags_weights[~sig_clip_mask])
        # std_mag = mags_numpy[~sig_clip_mask].std()
        mean_mag = current_star_table['mag_instrumental']
        std_mag = current_star_table['mag_instrumental_sigma']
        # mean_mag = np.average(mags_numpy, weights=mags_weights)
        # std_mag = mags_numpy.std()
        filter_column = unique_filter.lower()
        sigma_column = f'{filter_column}_sigma'
        reformatted_large_stars_table[filter_column][i] = mean_mag
        reformatted_large_stars_table[sigma_column][i] = std_mag
        # X_numpy = np.array(current_star_table['X'])
        # mean_X = X_numpy[~sig_clip_mask].mean()
        # std_X = X_numpy[~sig_clip_mask].std()
        # # mean_X = X_numpy.mean()
        # # std_X = X_numpy.std()
        mean_X = current_star_table['X']
        std_X = np.nan
        X_column = f'X_{filter_column}'
        X_std_column = f'X_{filter_column}_sigma'
        reformatted_large_stars_table[X_column][i] = mean_X
        reformatted_large_stars_table[X_std_column][i] = std_X

        # for i, unique_star in enumerate(multiple_stars):
        # star_mask = stars_table['Name'] == unique_star
        # current_star = stars_table[star_mask]
        for instr_index_name in all_indices_formatted:
            first_mag = reformatted_large_stars_table[instr_index_name[0]][i]
            second_mag = reformatted_large_stars_table[instr_index_name[-1]][i]
            instr_index = first_mag - second_mag
            reformatted_large_stars_table[instr_index_name][i] = instr_index
        i += 1
    # print(different_filter_list)
    return reformatted_large_stars_table


def calc_gb_first_transforms_AVG(stars_table, different_filter_list, save_loc,
                                 plot_results=False, save_plots=False):
    gb_transform_table_columns = init_gb_transform_table_columns_AVG()
    for different_filter in different_filter_list:
        # sigma_column = f'{different_filter}_sigma'
        X_column = f'X_{different_filter}'
        X_std_column = f'X_{different_filter}_sigma'
        stars_table.sort(X_column)
        for k, g in groupby(np.array(stars_table[X_column]
                                     [~np.isnan(stars_table[X_column])]),
                            key=Delta(0.05)):
            list_of_airmasses = list(g)
            mask = np.in1d(stars_table[X_column], list_of_airmasses)
            current_airmass_table = stars_table[mask]
            try:
                avg_airmass = np.average(
                    current_airmass_table[X_column],
                    weights=current_airmass_table[X_std_column])
            except ZeroDivisionError:
                avg_airmass = np.mean(current_airmass_table[X_column])
            unique_id = f"airmsass_{avg_airmass:0.3f}"
            # current_airmass_table.pprint(max_width=-1)
            # app_mag, app_mag_sigma, app_filter, colour_index =
            # get_app_mag_and_index_AVG(current_airmass_table,
            # different_filter)
            colour_indices = get_all_colour_indices(different_filter)
            for colour_index in colour_indices:
                try:
                    c_fci, c_fci_sigma, zprime_f, zprime_f_sigma =\
                        ground_based_first_order_transforms_AVG(
                            current_airmass_table,
                            different_filter,
                            colour_index,
                            plot_results=plot_results,
                            save_plots=save_plots,
                            unique_id=unique_id,
                            save_loc=save_loc)
                except TypeError:
                    print(f"Only 1 star at airmass: {avg_airmass:0.3f}.")
                    continue
                gb_transform_table_columns =\
                    update_gb_transform_table_columns_AVG(
                        gb_transform_table_columns,
                        c_fci,
                        c_fci_sigma,
                        zprime_f,
                        zprime_f_sigma,
                        different_filter,
                        colour_index,
                        avg_airmass)
    gb_transform_table = create_gb_transform_table_AVG(
        gb_transform_table_columns)
    return gb_transform_table


def init_gb_transform_table_columns_AVG():
    c_fci = []
    c_fci_simga = []
    zprime_f = []
    zprime_f_sigma = []
    instr_filter = []
    colour_index = []
    airmass = []
    gb_transform_table_columns = namedtuple('gb_transform_table_columns',
                                            ['c_fci',
                                             'c_fci_sigma',
                                             'zprime_f',
                                             'zprime_f_sigma',
                                             'instr_filter',
                                             'colour_index',
                                             'airmass'])
    return gb_transform_table_columns(c_fci,
                                      c_fci_simga,
                                      zprime_f,
                                      zprime_f_sigma,
                                      instr_filter,
                                      colour_index,
                                      airmass)


def get_app_mag_and_index_AVG(stars_table, instr_filter):
    if instr_filter == 'b' or instr_filter == 'u':
        colour_index = 'B-V'
        app_filter = 'B'
        app_mag = np.array(stars_table['V_ref'] + stars_table[colour_index])
        app_mag_sigma = np.nan_to_num(stars_table['V_sigma'],
                                      nan=max(stars_table['V_sigma'])) + \
            np.nan_to_num(stars_table[f'e_{colour_index}'], nan=max(
                stars_table[f'e_{colour_index}']))
        # app_mag_sigma =\
        # np.array(ref_star['e_V'] + ref_star[f'e_{colour_index}'])
    elif instr_filter == 'v' or instr_filter == 'g':
        colour_index = 'B-V'
        app_filter = 'V'
        app_mag = np.array(stars_table['V_ref'])
        app_mag_sigma = np.nan_to_num(
            stars_table['V_sigma'], nan=max(stars_table['V_sigma']))
        # app_mag_sigma = np.array(ref_star['e_V'])
    elif instr_filter == 'r':
        colour_index = 'V-R'
        app_filter = 'R'
        app_mag = np.array(stars_table['V_ref'] - stars_table[colour_index])
        app_mag_sigma = np.nan_to_num(stars_table['V_sigma'],
                                      nan=max(stars_table['V_sigma'])) + \
            np.nan_to_num(stars_table[f'e_{colour_index}'], nan=max(
                stars_table[f'e_{colour_index}']))
        # app_mag_sigma = np.array(ref_star['e_V'] + ref_star[f'e_{colour_index}'])
    elif instr_filter == 'i':
        colour_index = 'V-I'
        app_filter = 'I'
        app_mag = np.array(stars_table['V_ref'] - stars_table[colour_index])
        app_mag_sigma = np.nan_to_num(stars_table['V_sigma'],
                                      nan=max(stars_table['V_sigma'])) + \
            np.nan_to_num(stars_table[f'e_{colour_index}'], nan=max(
                stars_table[f'e_{colour_index}']))
        # app_mag_sigma = np.array(ref_star['e_V'] + \
        # ref_star[f'e_{colour_index}'])
    else:
        colour_index = None
        app_filter = None
        app_mag = None
        app_mag_sigma = None
    return app_mag, app_mag_sigma, app_filter, colour_index


def update_gb_transform_table_columns_AVG(gb_transform_table_columns,
                                          c_fci,
                                          c_fci_sigma,
                                          zprime_f,
                                          zprime_f_sigma,
                                          instr_filter,
                                          colour_index,
                                          avg_airmass):
    updated_gb_transform_table_columns = gb_transform_table_columns
    updated_gb_transform_table_columns.c_fci.append(c_fci)
    updated_gb_transform_table_columns.c_fci_sigma.append(c_fci_sigma)
    updated_gb_transform_table_columns.zprime_f.append(zprime_f)
    updated_gb_transform_table_columns.zprime_f_sigma.append(zprime_f_sigma)
    updated_gb_transform_table_columns.instr_filter.append(instr_filter)
    updated_gb_transform_table_columns.colour_index.append(colour_index)
    # avg_airmass = get_avg_airmass(altazpositions)
    updated_gb_transform_table_columns.airmass.append(avg_airmass)
    return updated_gb_transform_table_columns


def create_gb_transform_table_AVG(gb_transform_table_columns):
    gb_transform_table = Table(
        names=[
            'C_fCI',
            'C_fCI_sigma',
            'Zprime_f',
            'Zprime_f_sigma',
            'filter',
            'CI',
            'X'
        ],
        data=[
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


def ground_based_first_order_transforms_AVG(stars_table,
                                            instr_filter,
                                            colour_index,
                                            plot_results=False,
                                            save_plots=False, **kwargs):
    try:
        len(stars_table)
    except TypeError:
        return
    app_mag, app_mag_sigma, app_filter, _ = get_app_mag_and_index_AVG(
        stars_table, instr_filter)
    sigma_column = f'{instr_filter}_sigma'
    max_instr_filter_sigma = max(stars_table[sigma_column])
    err_sum = app_mag_sigma + \
        np.nan_to_num(stars_table[sigma_column], nan=max_instr_filter_sigma)
    err_sum = np.array(err_sum)
    err_sum[err_sum == 0] = max(err_sum)
    x = stars_table[colour_index][~np.isnan(stars_table[colour_index])]
    y = app_mag[~np.isnan(stars_table[colour_index])] - \
        stars_table[instr_filter][~np.isnan(stars_table[colour_index])]
    fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
    # print(len(x))
    # print(len(y))
    # print(stars_table[colour_index])
    if len(x) > 1 and len(y) > 1:
        # , weights=1.0 / (err_sum[~np.isnan(stars_table[colour_index])]))
        fitted_line, mask = or_fit(line_init, x, y)
    else:
        return
    filtered_data = np.ma.masked_array(y, mask=mask)
    c_fci = fitted_line.slope.value
    zprime_f = fitted_line.intercept.value
    if c_fci == 1 and zprime_f == 0:
        return
    cov = fit.fit_info['param_cov']
    if cov is None:
        c_fci_sigma = 0.0
        zprime_f_sigma = 0.0
    else:
        c_fci_sigma = sqrt(cov[0][0])
        zprime_f_sigma = sqrt(cov[1][1])
    if plot_results:
        # print(min(stars_table[colour_index]))
        # print(max(stars_table[colour_index]))
        index_plot\
            = np.arange(start=min(stars_table[colour_index]
                                  [~np.isnan(stars_table[colour_index])]),
                        stop=max(stars_table[colour_index][~np.isnan(
                            stars_table[colour_index])]) + 0.01,
                        step=0.01)
        plt.errorbar(x, y, yerr=err_sum[~np.isnan(stars_table[colour_index])],
                     color='#1f77b4', fmt='o',
                     fillstyle='none', capsize=2, label="Clipped Data")
        plt.plot(x, filtered_data, 'o', color='#1f77b4', label="Fitted Data")
        plt.plot(index_plot, fitted_line(index_plot), '-', color='#ff7f0e',
                 label=f"({app_filter}-{instr_filter}) = {c_fci:.3f} * {colour_index} + {zprime_f:.3f}")
        # plt.plot(index_plot, c_fci * index_plot + zprime_f,
        #          label=f"({app_filter}-{instr_filter}) = {c_fci:.3f} * \
        # {colour_index} + {zprime_f:.3f}")
        plt.ylabel(f"{app_filter}-{instr_filter}")
        plt.xlabel(f"{colour_index}")
        plt.legend()
        plt.title(f"C and Z' Coefficient Calculations for {colour_index}")
        # if not field:
        #     plt.title(f"({app_filter}-{instr_filter}) = {c_fci:.3f} * \
        # {colour_index} + {zprime_f:.3f}")
        # else:
        #     plt.title(f"{field}: ({app_filter}-{instr_filter})\
        # = {c_fci:.3f} * {colour_index} + {zprime_f:.3f}")
        if save_plots:
            unique_id = kwargs.get('unique_id')
            save_loc\
                = f"{os.path.join(kwargs.get('save_loc'),f'CZprime{app_filter}-{colour_index}_{unique_id}')}.png"
            plt.savefig(save_loc)
        plt.show()
        plt.close()
    return c_fci, c_fci_sigma, zprime_f, zprime_f_sigma


def group_each_star(large_stars_table, ground_based=False, keys='Name'):
    """
    Group all detections of unique reference stars together.

    Parameters
    ----------
    large_stars_table : astropy.table.table.QTable
        Table containing all information on the detected stars. Has columns:
            Field : string
                Unique identifier of the star field that the reference star is
                in (e.g. Landolt field "108").
            Name : string
                Name/unique identifier of the reference star.
            Time (JD) : numpy.float64
                Time of the observation.
            RA_ref : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) from the file in unit
                hourangle.
            dec_ref : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) from the file in
                unit degree.
            RA_img : astropy.coordinates.angles.Longitude
                Right ascension of the reference star(s) found in the 
                image in unit hourangle.
            dec_img : astropy.coordinates.angles.Latitude
                Declination of the reference star(s) found in the image
                in unit degree.
            angular_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in
                the image and ref_stars_file.
            flux : numpy.float64
                Flux of the source.
            exposure : numpy.float64
                Exposure time of the image.
            mag_instrumental : numpy.float64
                Instrumental magnitude of the reference star(s) found
                in the image.
            mag_instrumental_sigma : numpy.float64
                Standard deviation of the instrumental magnitude of the
                reference star(s) found in the image.
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
                Standard deviation of the apparent V magnitude from the
                reference file.
            X : numpy.float64
                Sec(z) of the reference star(s) found in the image.
    ground_based : bool, optional
        Whether or not the image is from a ground-based sensor.
        The default is False.
    keys : string, optional
        Table column(s) to use to group the different reference stars.
        The default is 'Name'.

    Returns
    -------
    stars_table : astropy.table.table.Table
        Table containing the mean of the important information for each star.
        Has columns:
            Field : string
                Unique identifier of the star field that the reference star
                is in (e.g. Landolt field "108").
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
                Mean instrumental magnitude of all detections of the star
                in <filter>. There is a different column for 
                each different filter used across the images.
            <filter>_sigma : numpy.float64
                Standard deviation of the instrumental magnitudes of all
                detections of the star in <filter>. 
                There is a different column for each different filter
                used across the images.
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
    different_filter_list = [different_filter.lower()
                             for different_filter in different_filter_list]
    different_filter_data = np.empty((N, len(different_filter_list)))
    different_filter_data.fill(np.nan)
    filter_sigma_list = []
    for different_filter in different_filter_list:
        filter_sigma_list.append(f"{different_filter}_sigma")
    different_filter_table = Table(
        data=different_filter_data, names=different_filter_list)
    different_filter_sigma_table = Table(
        data=different_filter_data, names=filter_sigma_list)
    if ground_based:
        filter_X_list = []
        for different_filter in different_filter_list:
            filter_X_list.append(f"X_{different_filter}")
        filter_X_sigma_list = []
        for different_filter in different_filter_list:
            filter_X_sigma_list.append(f"X_{different_filter}_sigma")
        filter_X_table = Table(data=different_filter_data, names=filter_X_list)
        filter_X_sigma_table = Table(
            data=different_filter_data, names=filter_X_sigma_list)
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
    if not ground_based:
        for star in unique_stars['Name']:
            stars_table['Name'][i] = star
            mask = large_stars_table['Name'] == star
            current_star_table = large_stars_table[mask]
            unique_filters = table.unique(current_star_table, keys='filter')
            for unique_filter in unique_filters['filter']:
                mask = ((current_star_table['filter'] == unique_filter))
                current_star_filter_table = current_star_table[mask]
                mags_numpy = np.array(
                    current_star_filter_table['mag_instrumental'])
                mags_std_numpy = np.array(
                    current_star_filter_table['mag_instrumental_sigma'])
                mags_weights = 1 / (mags_std_numpy ** 2)
                # mean_mag = mags_numpy.mean()
                mean_mag = np.average(mags_numpy, weights=mags_weights)
                std_mag = mags_numpy.std()
                filter_column = unique_filter.lower()
                sigma_column = f'{filter_column}_sigma'
                stars_table[filter_column][i] = mean_mag
                stars_table[sigma_column][i] = std_mag
            i += 1
    else:
        for star in unique_stars['Name']:
            mask = large_stars_table['Name'] == star
            current_star_table = large_stars_table[mask]
            current_star_table.sort('X')
            # current_star_table.pprint(max_lines=-1, max_width=300)
            # data = np.array(current_star_table['X'])
            # gaps = [y - x for x, y in zip(data[:-1], data[1:])]
            # # have python calculate the standard deviation for the gaps
            # sd = np.std(gaps)

            # # create a list of lists, put the first value of the source data
            # in the first
            # lists = [[data[0]]]
            # for x in data[1:]:
            #     # if the gap from the current item to the previous is more
            # than 1 SD
            #     # Note: the previous item is the last item in the last list
            #     # Note: the '> 1' is the part you'd modify to make it
            # stricter or more relaxed
            #     if (x - lists[-1][-1]) / sd > 3:
            #         # then start a new list
            #         lists.append([])
            #     # add the current item to the last list in the list
            #     lists[-1].append(x)

            # print(lists)
            # if len(current_star_table) > 1:
            for k, g in groupby(np.array(current_star_table['X']),
                                key=Delta(0.1)):
                stars_table['Name'][i] = star
                list_of_airmasses = list(g)
                mask = np.in1d(current_star_table['X'], list_of_airmasses)
                # print(mask)
                current_star_airmass_table = current_star_table[mask]
                print(current_star_airmass_table)
                # print(np.array(current_star_table['X']))
                print(list_of_airmasses)
                unique_filters = table.unique(
                    current_star_table, keys='filter')
                for unique_filter in unique_filters['filter']:
                    mask = ((current_star_table['filter'] == unique_filter))
                    current_star_filter_table = current_star_table[mask]
                    mags_numpy = np.array(
                        current_star_filter_table['mag_instrumental'])
                    mags_std_numpy = np.array(
                        current_star_filter_table['mag_instrumental_sigma'])
                    mags_weights = 1 / (mags_std_numpy ** 2)
                    # mean_mag = mags_numpy.mean()
                    mean_mag = np.average(mags_numpy, weights=mags_weights)
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
    """
    Write an Astropy table formatted to be inserted into a LaTeX document.

    Parameters
    ----------
    table : astropy.table.Table
        Table for which to write to LaTeX.
    output_file : string
        File location to write the table to.
    formats : dict, optional
        Dictionary of format specifiers or formatting functions.
        The default is None.

    Returns
    -------
    None.

    """
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
    fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
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


def init_sb_final_transform_columns():
    """
    Initialize the columns that will create the table used for the space-based
    transforms.

    Returns
    -------
    sb_final_transform_columns : namedtuple
        Attributes:
            index : empty list
                Name of the colour index used to calculate filter_fci and
                zprime_fci.
            filter_fci : empty list
                T coefficient for index.
            filter_fci_sigma : empty list
                Standard deviation of the T coefficient for index.
            zprime_fci : empty list
                Zero point for index.
            zprime_fci_sigma : empty list
                Standard deviation of the zero point for index.
    """
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
    """
    Update columns to be used for the transform table based on information
    from the current image.

    Parameters
    ----------
    sb_final_transform_columns : namedtuple
        Attributes:
            index : array-like
                Name of the colour index used to calculate filter_fci and
                zprime_fci.
            filter_fci : np.float64
                T coefficient for index.
            filter_fci_sigma : np.float64
                Standard deviation of the T coefficient for index.
            zprime_fci : np.float64
                Zero point for index.
            zprime_fci_sigma : np.float64
                Standard deviation of the zero point for index.
    index : string
        Name of the colour index used to calculate filter_fci and zprime_fci.
    filter_fci : float
        T coefficient for index.
    filter_fci_sigma : float
        Standard deviation of the T coefficient for index.
    zprime_fci : float
        Zero point for index.
    zprime_fci_sigma : float
        Standard deviation of the zero point for index.

    Returns
    -------
    updated_sb_final_transform_columns : namedtuple
        Attributes:
            index : array-like
                Name of the colour index used to calculate filter_fci and
                zprime_fci.
            filter_fci : np.float64
                T coefficient for index.
            filter_fci_sigma : np.float64
                Standard deviation of the T coefficient for index.
            zprime_fci : np.float64
                Zero point for index.
            zprime_fci_sigma : np.float64
                Standard deviation of the zero point for index.

    """
    updated_sb_final_transform_columns = sb_final_transform_columns
    updated_sb_final_transform_columns.index.append(index)
    updated_sb_final_transform_columns.filter_fci.append(filter_fci)
    updated_sb_final_transform_columns.filter_fci_sigma.append(
        filter_fci_sigma)
    updated_sb_final_transform_columns.zprime_fci.append(zprime_fci)
    updated_sb_final_transform_columns.zprime_fci_sigma.append(
        zprime_fci_sigma)
    return updated_sb_final_transform_columns


def create_sb_final_transform_table(sb_final_transform_columns):
    """
    Convert the columns of the space-based transform table into an AstroPy
    table.

    Parameters
    ----------
    sb_final_transform_columns : namedtuple
        Attributes:
            index : array-like
                Name of the colour index used to calculate filter_fci and
                zprime_fci.
            filter_fci : np.float64
                T coefficient for index.
            filter_fci_sigma : np.float64
                Standard deviation of the T coefficient for index.
            zprime_fci : np.float64
                Zero point for index.
            zprime_fci_sigma : np.float64
                Standard deviation of the zero point for index.

    Returns
    -------
    sb_final_transform_table : astropy.table.Table
        Table containing results of the final space-based transforms.
        Has columns:
            CI : string
                Name of the colour index used to calculate the corresponding
                T_fCI and Z_fCI.
            T_fCI : float
                T coefficient for the corresponding CI.
            T_fCI_sigma : float
                Standard deviation of the T coefficient for the corresponding
                CI.
            Z_fCI : float
                Zero point for the corresponding CI.
            Z_fCI_sigma : float
                Standard deviation of the Zero point for the corresponding
                CI.

    """
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
    Initialize the columns that will create the table used for the
    ground-based transforms.

    Returns
    -------
    gb_transform_table_columns : namedtuple
        Attributes:
            field : empty list
                Unique identifier of the star field that the reference
                star is in (e.g. Landolt field "108").
            c_fci : empty list
                C coefficient for filter f with colour index ci.
            c_fci_sigma : empty list
                Standard deviation of the C coefficient for filter f
                with colour index ci.
            zprime_f : empty list
                Z' coefficient for filter f.
            zprime_f_sigma : empty list
                Standard deviation of the Z' coefficient for filter f.
            instr_filter : empty list
                Instrumental filter band to calculate the transform for.
            colour_index : empty list
                Name of the colour index used to calculate c_fci and zprime_f.
            airmass : empty list
                The mean airmass for all sources in the image.

    """
    field = []
    c_fci = []
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
    Update columns to be used for the transform table based on information
    from the current image.

    Parameters
    ----------
    gb_transform_table_columns : namedtuple
        Attributes:
            field : string list
                Unique identifier of the star field that the reference star
                is in (e.g. Landolt field "108").
            c_fci : np.float64
                C coefficient for filter f with colour index ci.
            c_fci_sigma : np.float64
                Standard deviation of the C coefficient for filter f with
                colour index ci.
            zprime_f : numpy.float64
                Z' coefficient for filter f.
            zprime_f_sigma : numpy.float64
                Standard deviation of the Z' coefficient for filter f.
            instr_filter : string list
                Instrumental filter band to calculate the transform for.
            colour_index : string list
                Name of the colour index used to calculate c_fci and zprime_f.
            airmass : numpy.float64
                The mean airmass for all sources in the image.

    field : string
        Unique identifier of the star field that the reference star is in
        (e.g. Landolt field "108").
    c_fci : float
        C coefficient for filter f with colour index ci.
    c_fci_sigma : float
        Standard deviation of the C coefficient for filter f with colour
        index ci.
    zprime_f : float
        Z' coefficient for filter f.
    zprime_f_sigma : float
        Standard deviation of the Z' coefficient for filter f.
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
                Unique identifier of the star field that the reference
                star is in (e.g. Landolt field "108").
            c_fci : np.float64
                C coefficient for filter f with colour index ci.
            c_fci_sigma : np.float64
                Standard deviation of the C coefficient for filter f
                with colour index ci.
            zprime_f : numpy.float64
                Z' coefficient for filter f.
            zprime_f_sigma : numpy.float64
                Standard deviation of the Z' coefficient for filter f.
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
                Unique identifier of the star field that the reference star
                is in (e.g. Landolt field "108").
            c_fci : np.float64
                C coefficient for filter f with colour index ci.
            c_fci_sigma : np.float64
                Standard deviation of the C coefficient for filter f with
                colour index ci.
            zprime_f : numpy.float64
                Z' coefficient for filter f.
            zprime_f_sigma : numpy.float64
                Standard deviation of the Z' coefficient for filter f.
            instr_filter : string list
                Instrumental filter band to calculate the transform for.
            colour_index : string list
                Name of the colour index used to calculate c_fci and zprime_f.
            airmass : numpy.float64
                The mean airmass for all sources in the image.

    Returns
    -------
    gb_transform_table : astropy.table.Table
        Table containing the results from the intermediary step in
        calculating the transforms. Has columns:
            field : string list
                Unique identifier of the star field that the reference star
                is in (e.g. Landolt field "108").
            C_fCI : np.float64
                C coefficient for filter f with colour index ci.
            C_fCI_sigma : np.float64
                Standard deviation of the C coefficient for filter f with
                colour index ci.
            Zprime_f : numpy.float64
                Z' coefficient for filter f.
            Zprime_f_sigma : numpy.float64
                Standard deviation of the Z' coefficient for filter f.
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
    Get the desired apparent magnitude and colour filter for the
    particular instr_filter.

    Parameters
    ----------
    ref_star : astropy.table.Table
        Rows from ref_stars_file that correspond to a matched reference star.
    instr_filter : string
        Instrumental filter band to calculate the transform for.

    Returns
    -------
    app_mag : numpy.float64
        Apparent magnitude of the reference star(s) in the desired filter
        (e.g. B for b, V for v, etc.).
    app_mag_sigma : numpy.float64
        Standard deviation of the apparent magnitude of the reference star(s)
        in the same filter as app_mag.
    colour_index : string
        Name of the colour index to use when calculating the transform
        (e.g. B-V for b, V-R for r).

    """
    if instr_filter == 'b' or instr_filter == 'u':
        colour_index = 'B-V'
        app_filter = 'B'
        app_mag = np.array(ref_star['V_ref'] + ref_star[colour_index])
        app_mag_sigma = np.nan_to_num(ref_star['e_V'],
                                      nan=max(ref_star['e_V'])) + \
            np.nan_to_num(ref_star[f'e_{colour_index}'], nan=max(
                ref_star[f'e_{colour_index}']))
        # app_mag_sigma = np.array(ref_star['e_V'] + \
        # ref_star[f'e_{colour_index}'])
    elif instr_filter == 'v' or instr_filter == 'g':
        colour_index = 'B-V'
        app_filter = 'V'
        app_mag = np.array(ref_star['V_ref'])
        app_mag_sigma = np.nan_to_num(
            ref_star['e_V'], nan=max(ref_star['e_V']))
        # app_mag_sigma = np.array(ref_star['e_V'])
    elif instr_filter == 'r':
        colour_index = 'V-R'
        app_filter = 'R'
        app_mag = np.array(ref_star['V_ref'] - ref_star[colour_index])
        app_mag_sigma = np.nan_to_num(ref_star['e_V'],
                                      nan=max(ref_star['e_V'])) + \
            np.nan_to_num(ref_star[f'e_{colour_index}'], nan=max(
                ref_star[f'e_{colour_index}']))
        # app_mag_sigma = np.array(ref_star['e_V'] +\
        # ref_star[f'e_{colour_index}'])
    elif instr_filter == 'i':
        colour_index = 'V-I'
        app_filter = 'I'
        app_mag = np.array(ref_star['V_ref'] - ref_star[colour_index])
        app_mag_sigma = np.nan_to_num(ref_star['e_V'],
                                      nan=max(ref_star['e_V'])) + \
            np.nan_to_num(ref_star[f'e_{colour_index}'], nan=max(
                ref_star[f'e_{colour_index}']))
        # app_mag_sigma = np.array(ref_star['e_V'] + \
        # ref_star[f'e_{colour_index}'])
    else:
        colour_index = None
        app_filter = None
        app_mag = None
        app_mag_sigma = None
    return app_mag, app_mag_sigma, app_filter, colour_index


def ground_based_first_order_transforms(matched_stars, instr_filter,
                                        colour_index, field=None,
                                        plot_results=False,
                                        save_plots=False, **kwargs):
    """
    Perform the intermediary step to calculating the ground based transforms.

    Parameters
    ----------
    matched_stars : namedtuple
        Attributes:
            ref_star_index : array-like or int
                Index of the star(s) in ref_stars_file that correspond to a
                star in the image.
            img_star_index : array-like or int
                Index of the star(s) detected in the image that correspond to
                a star in ref_stars_file.
            ref_star : astropy.table.Table
                Rows from ref_stars_file that correspond to a matched
                reference star.
            ref_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched reference star(s) from ref_stars_file.
            img_star_loc : astropy.coordinates.sky_coordinate.SkyCoord
                RA/dec of the matched star(s) detected in the image.
            ang_separation : astropy.coordinates.angles.Angle
                Angular distance between the matched star(s) detected in the
                image and ref_stars_file.
            img_instr_mag : numpy array
                Array containing the instrumental magnitudes of the matched
                star(s) detected in the image.
            img_instr_mag_sigma : numpy array
                Array containing the standard deviations of the instrumental
                magnitudes of the matched star(s) 
                detected in the image.
            flux : numpy array
                Array containing the non-normalized fluxes of the matched
                star(s) detected in the image.
            img_star_altaz : astropy.coordinates.sky_coordinate.SkyCoord
                Alt/Az of the matched star(s) detected in the image.
            img_star_airmass : float
                sec(z) of img_star_altaz.
    instr_filter : string
        Instrumental filter band to calculate the transform for.
    colour_index : string
        The colour index to calculate the transforms for.
    plot_results : bool, optional
        Controls whether or not to plot the results from the transforms.
        The default is False.
    field : string, optional
        Unique identifier of the star field that the reference star is in
        (e.g. Landolt field "108"). 
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
    app_mag, app_mag_sigma, app_filter, _ = get_app_mag_and_index(
        matched_stars.ref_star, instr_filter)
    max_instr_filter_sigma = max(matched_stars.img_instr_mag_sigma)
    err_sum = app_mag_sigma + \
        np.nan_to_num(matched_stars.img_instr_mag_sigma,
                      nan=max_instr_filter_sigma)
    err_sum = np.array(err_sum)
    err_sum[err_sum == 0] = max(err_sum)
    x = matched_stars.ref_star[colour_index]
    y = app_mag - matched_stars.img_instr_mag
    fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
    fitted_line, mask = or_fit(line_init, x, y, weights=1.0 / err_sum)
    filtered_data = np.ma.masked_array(y, mask=mask)
    c_fci = fitted_line.slope.value
    zprime_f = fitted_line.intercept.value
    if c_fci == 1 and zprime_f == 0:
        return
    cov = fit.fit_info['param_cov']
    if cov is None:
        c_fci_sigma = 0.0
        zprime_f_sigma = 0.0
    else:
        c_fci_sigma = sqrt(cov[0][0])
        zprime_f_sigma = sqrt(cov[1][1])
    if plot_results:
        index_plot = np.arange(start=min(matched_stars.ref_star[colour_index]),
                               stop=max(
                                   matched_stars.ref_star[colour_index]) + 0.01,
                               step=0.01)
        plt.errorbar(x, y, yerr=err_sum, color='#1f77b4', fmt='o',
                     fillstyle='none', capsize=2, label="Clipped Data")
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
    """
    Get all of the colour indices to calculate the transforms with the
    instrumental mag filter.

    Parameters
    ----------
    instr_filter : string
        Instrumental filter band to calculate the transform for.

    Returns
    -------
    colour_indices : string list
        List of all of the colour indices to use when calculating the
        transforms.

    """
    if instr_filter == 'b' or instr_filter == 'u':
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
    """
    Remove all observations with airmass above a set point.

    Parameters
    ----------
    gb_transform_table : TYPE
        DESCRIPTION.
    max_airmass : float, optional
        Airmass value for which to remove all observations greater than.
        The default is 2.0.

    Returns
    -------
    gb_transform_table : astropy.table.Table
        Table containing the results from the intermediary step in
        calculating the transforms. Has columns:
            field : string list
                Unique identifier of the star field that the reference star
                is in (e.g. Landolt field "108").
            C_fCI : np.float64
                C coefficient for filter f with colour index ci.
            C_fCI_sigma : np.float64
                Standard deviation of the C coefficient for filter f with
                colour index ci.
            Zprime_f : numpy.float64
                Z' coefficient for filter f.
            Zprime_f_sigma : numpy.float64
                Standard deviation of the Z' coefficient for filter f.
            filter : string list
                Instrumental filter band to calculate the transform for.
            CI : string list
                Name of the colour index used to calculate c_fci and zprime_f.
            X : numpy.float64
                The mean airmass for all sources in the image.

    """
    gb_transform_table\
        = gb_transform_table[gb_transform_table['X'] <= max_airmass]
    return gb_transform_table


def ground_based_second_order_transforms(gb_transform_table,
                                         plot_results=False,
                                         save_plots=False,
                                         **kwargs):
    """
    Perform the final step in calculating the transforms for a ground-based
    observatory.

    Parameters
    ----------
    gb_transform_table : astropy.table.Table
        Table containing the results from the intermediary step in calculating
        the transforms. Has columns:
            field : string list
                Unique identifier of the star field that the reference star
                is in (e.g. Landolt field "108").
            C_fCI : np.float64
                C coefficient for filter f with colour index ci.
            C_fCI_sigma : np.float64
                Standard deviation of the C coefficient for filter f with
                colour index ci.
            Zprime_f : numpy.float64
                Z' coefficient for filter f.
            Zprime_f_sigma : numpy.float64
                Standard deviation of the Z' coefficient for filter f.
            filter : string list
                Instrumental filter band to calculate the transform for.
            CI : string list
                Name of the colour index used to calculate c_fci and zprime_f.
            X : numpy.float64
                The mean airmass for all sources in the image.
    plot_results : bool, optional
        Controls whether or not to plot the results from the transforms.
        The default is False.

    Returns
    -------
    gb_final_transforms : astropy.table.Table
        Table containing the results for the final transforms. Has columns:
            filter : string
                Instrumental filter band used to calculate the transform.
            CI : string
                Name of the colour index used to calculate the transform
                (e.g. B-V for b, V-R for r).
            k''_fCI : float
                The second order atmospheric extinction coefficient for
                filter f using the colour index CI.
            k''_fCI_sigma : float
                The standard deviation of the second order atmospheric
                extinction coefficient for filter f using the 
                colour index CI.
            T_fCI : float
                The instrumental transform coefficient for filter f
                using the colour index CI.
            T_fCI_sigma : float
                The standard deviation of the instrumental transform
                coefficient for filter f using the colour index CI.
            k'_f : float
                The first order atmospheric extinction coefficient
                for filter f.
            k'_f_sigma : float
                The standard deviation of the first order atmospheric
                extinction coefficient for filter f.
            Z_f : float
                The zero point magnitude for filter f.
            Z_f_sigma : float
                The standard deviation of the zero point magnitude for
                filter f.

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
        mask = ((gb_transform_table['filter'] == unique_filter) & (
            gb_transform_table['CI'] == current_index))
        current_filter = gb_transform_table[mask]
        x = current_filter['X']
        y = current_filter['C_fCI']
        max_c_fci_sigma = max(current_filter['C_fCI_sigma'])
        sigma = np.nan_to_num(
            current_filter['C_fCI_sigma'], nan=max_c_fci_sigma)
        sigma = np.array(sigma)
        sigma[sigma == 0] = max(sigma)
        # sigma = current_filter['C_fCI_sigma']
        fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
        if np.sum(np.abs(sigma)) == 0:
            fitted_line_c, mask = or_fit(line_init, x, y)
        elif len(x) > 1 and len(y) > 1:
            fitted_line_c, mask = or_fit(line_init, x, y, weights=1.0 / sigma)
        else:
            continue
        filtered_data_c = np.ma.masked_array(y, mask=mask)
        kprimeprime_fci = fitted_line_c.slope.value
        t_fci = fitted_line_c.intercept.value
        cov_c = fit.fit_info['param_cov']
        # c_fci_sigma = cov[0][0]
        # zprime_f_sigma = cov[1][1]
        # a_fit_c, cov_c = curve_fit(linear_func, current_filter['X'],
        # current_filter['C_fCI'],
        #                            sigma=current_filter['C_fCI_sigma'])
        # kprimeprime_fci = a_fit_c[0]
        # t_fci = a_fit_c[1]
        if cov_c is None:
            kprimeprime_fci_sigma = np.nan
            t_fci_sigma = np.nan
        else:
            kprimeprime_fci_sigma = sqrt(cov_c[0][0])
            t_fci_sigma = sqrt(cov_c[1][1])
        # kprimeprime_fci, t_fci = np.polyfit(current_filter['X'],
        # current_filter['C_fCI'], 1)
        y = current_filter['Zprime_f']
        max_zprime_f_sigma = max(current_filter['Zprime_f_sigma'])
        sigma = np.nan_to_num(
            current_filter['Zprime_f_sigma'], nan=max_zprime_f_sigma)
        sigma = np.array(sigma)
        sigma[sigma == 0] = max(sigma)
        # sigma = current_filter['Zprime_f_sigma']
        fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
        if np.sum(np.abs(sigma)) == 0:
            fitted_line_z, mask = or_fit(line_init, x, y)
        else:
            fitted_line_z, mask = or_fit(line_init, x, y, weights=1.0 / sigma)
        filtered_data_z = np.ma.masked_array(y, mask=mask)
        kprime_f = fitted_line_z.slope.value
        zprime_f = fitted_line_z.intercept.value
        cov_z = fit.fit_info['param_cov']
        # a_fit_z, cov_z = curve_fit(linear_func, current_filter['X'],
        # current_filter['Zprime_f'],
        #                            sigma=current_filter['Zprime_f_sigma'])
        # kprime_f = a_fit_z[0]
        # zprime_f = a_fit_z[1]
        if cov_z is None:
            kprime_f_sigma = np.nan
            zprime_f_sigma = np.nan
        else:
            kprime_f_sigma = sqrt(cov_z[0][0])
            zprime_f_sigma = sqrt(cov_z[1][1])
        # kprime_f, zprime_f = np.polyfit(current_filter['X'],
        # current_filter['Zprime_f'], 1)
        gb_final_transforms['k\'\'_fCI'][unique_filter_index] =\
            kprimeprime_fci
        gb_final_transforms['k\'\'_fCI_sigma'][unique_filter_index] =\
            kprimeprime_fci_sigma
        gb_final_transforms['T_fCI'][unique_filter_index] = t_fci
        gb_final_transforms['T_fCI_sigma'][unique_filter_index] =\
            t_fci_sigma
        gb_final_transforms['k\'_f'][unique_filter_index] =\
            kprime_f
        gb_final_transforms['k\'_f_sigma'][unique_filter_index] =\
            kprime_f_sigma
        gb_final_transforms['Z_f'][unique_filter_index] =\
            zprime_f
        gb_final_transforms['Z_f_sigma'][unique_filter_index] =\
            zprime_f_sigma
        if plot_results:
            X_plot = np.arange(start=min(current_filter['X']), stop=max(
                current_filter['X']) + 0.01, step=0.01)
            ci_plot = re.sub('[^a-zA-Z]+', '', current_index)
            ci_plot = ci_plot.lower()
            plt.errorbar(x, current_filter['C_fCI'],
                         yerr=current_filter['C_fCI_sigma'],
                         color='#1f77b4',
                         fmt='o',
                         fillstyle='none',
                         capsize=2, label="Clipped Data")
            plt.plot(x, filtered_data_c,
                     'o',
                     color='#1f77b4',
                     label="Fitted Data")
            plt.plot(X_plot, fitted_line_c(X_plot),
                     '-',
                     color='#ff7f0e',
                     label=f'C_{unique_filter}{ci_plot} = {kprimeprime_fci:.3f} * X + {t_fci:.3f}')
            plt.legend()
            plt.title(
                f"k''_{unique_filter}{ci_plot} and T_{unique_filter}{ci_plot} Coefficient calculations")
            # plt.title(f'C_{unique_filter}{ci_plot} = {kprimeprime_fci:.3f} * X + {t_fci:.3f}')
            plt.ylabel(f'C_{unique_filter}{ci_plot}')
            plt.xlabel('X')
            if save_plots:
                save_loc = f"{os.path.join(kwargs.get('save_loc'), f'KprimeprimeT{unique_filter}{ci_plot}')}.png"
                plt.savefig(save_loc)
            plt.show()
            plt.close()
            plt.errorbar(x, current_filter['Zprime_f'],
                         yerr=current_filter['Zprime_f_sigma'],
                         color='#1f77b4',
                         fmt='o',
                         fillstyle='none',
                         capsize=2,
                         label="Clipped Data")
            plt.plot(x, filtered_data_z, 'o',
                     color='#1f77b4',
                     label="Fitted Data")
            plt.plot(X_plot, fitted_line_z(X_plot),
                     '-',
                     color='#ff7f0e',
                     label=f'Z\'_{unique_filter} = {kprime_f:.3f} * X + {zprime_f:.3f}')
            plt.legend()
            plt.title(
                f"Z_{unique_filter} and k'_{unique_filter} Coefficient Calculations")
            # plt.title(f'Z\'_{unique_filter} = {kprime_f:.3f} * X + {zprime_f:.3f}')
            plt.ylabel(f'Z\'_{unique_filter}')
            plt.xlabel('X')
            if save_plots:
                save_loc = f"{os.path.join(kwargs.get('save_loc'), f'KprimeZ{unique_filter}{ci_plot}')}.png"
                plt.savefig(save_loc)
            plt.show()
            plt.close()

    return gb_final_transforms


def get_colour_index_lower(instr_filter):
    """
    Get the colour index for the desired instrumental filter band.

    Parameters
    ----------
    instr_filter : string
        Instrumental filter band used to calculate the transform.

    Returns
    -------
    colour_index : string
        Colour index of the form X-Y'.
    ci : string
        Colour index of the form 'xy'.

    """
    if instr_filter == 'b' or instr_filter == 'u':
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


# def init_unknown_object_table_columns():
#     """
#     Initialize the columns that will create the unknown object table.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     """
#     instr_filter = []
#     name = []
#     instr_mag = []

#     unknown_object_table_columns = namedtuple('unknown_object_table_columns',
#                                               ['instr_filter',
#                                               'name',
#                                               'instr_mag'])
#     return unknown_object_table_columns(instr_filter,
#                                         name,
#                                         instr_mag)


# def update_unknown_table_columns(unknown_object_table_columns, hdr, instr_mag, filter_key='FILTER'):
#     """
#     Update the columns tha will create the unknown object table.

#     Parameters
#     ----------
#     unknown_object_table_columns : TYPE
#         DESCRIPTION.
#     hdr : TYPE
#         DESCRIPTION.
#     instr_mag : TYPE
#         DESCRIPTION.
#     filter_key : TYPE, optional
#         DESCRIPTION. The default is 'FILTER'.

#     Returns
#     -------
#     updated_unknown_object_table_columns : TYPE
#         DESCRIPTION.

#     """
#     updated_unknown_object_table_columns = unknown_object_table_columns
#     instr_filter = get_instr_filter_name(hdr=hdr, filter_key=filter_key)
#     name = hdr['OBJECT']
#     updated_unknown_object_table_columns.instr_filter.append(instr_filter)
#     updated_unknown_object_table_columns.name.append(name)
#     updated_unknown_object_table_columns.instr_mag.append(instr_mag)
#     return updated_unknown_object_table_columns


def calculate_c_fci(gb_final_transforms, instr_filter, airmass, colour_index):
    """
    Calculate the C coefficient to use when calculating the C' coefficient when applying the transforms.

    Parameters
    ----------
    gb_final_transforms : astropy.table.Table
        Table containing the results for the final transforms. Has columns:
            filter : string
                Instrumental filter band used to calculate the transform.
            CI : string
                Name of the colour index used to calculate the transform 
                (e.g. B-V for b, V-R for r).
            k''_fCI : float
                The second order atmospheric extinction coefficient for 
                filter f using the colour index CI.
            k''_fCI_sigma : float
                The standard deviation of the second order atmospheric 
                extinction coefficient for filter f using the 
                colour index CI.
            T_fCI : float
                The instrumental transform coefficient for filter f using
                the colour index CI.
            T_fCI_sigma : float
                The standard deviation of the instrumental transform
                coefficient for filter f using the colour index CI.
            k'_f : float
                The first order atmospheric extinction coefficient for
                filter f.
            k'_f_sigma : float
                The standard deviation of the first order atmospheric
                extinction coefficient for filter f.
            Z_f : float
                The zero point magnitude for filter f.
            Z_f_sigma : float
                The standard deviation of the zero point magnitude for
                filter f.
    instr_filter : string
        The instrumental filter to use to calculate C_fCI.
    airmass : float
        Average airmass of the image.
    colour_index : string
        The colour index to use to calculate C_fCI.

    Returns
    -------
    c_fci : float
        C coefficient to use when calculating C'.

    """
    mask = ((gb_final_transforms['filter'] == instr_filter) & (
        gb_final_transforms['CI'] == colour_index))
    row_of_transforms = gb_final_transforms[mask]
    if len(row_of_transforms) == 0 and instr_filter == 'v':
        instr_filter = 'g'
        mask = ((gb_final_transforms['filter'] == instr_filter) & (
            gb_final_transforms['CI'] == colour_index))
        row_of_transforms = gb_final_transforms[mask]
    if len(row_of_transforms) == 0 and instr_filter == 'b':
        instr_filter = 'u'
        mask = ((gb_final_transforms['filter'] == instr_filter) & (
            gb_final_transforms['CI'] == colour_index))
        row_of_transforms = gb_final_transforms[mask]
    c_fci = float(row_of_transforms['T_fCI']) - \
        (float(row_of_transforms["k''_fCI"]) * airmass)
    return c_fci


def calculate_c_prime(gb_final_transforms, instr_filter, airmass):
    """
    Calculate the C' coefficient to use when applying the transforms.

    Parameters
    ----------
    gb_final_transforms : astropy.table.Table
        Table containing the results for the final transforms. Has columns:
            filter : string
                Instrumental filter band used to calculate the transform.
            CI : string
                Name of the colour index used to calculate the transform
                (e.g. B-V for b, V-R for r).
            k''_fCI : float
                The second order atmospheric extinction coefficient for
                filter f using the colour index CI.
            k''_fCI_sigma : float
                The standard deviation of the second order atmospheric
                extinction coefficient for filter f using the 
                colour index CI.
            T_fCI : float
                The instrumental transform coefficient for filter f using
                the colour index CI.
            T_fCI_sigma : float
                The standard deviation of the instrumental transform
                coefficient for filter f using the colour index CI.
            k'_f : float
                The first order atmospheric extinction coefficient for
                filter f.
            k'_f_sigma : float
                The standard deviation of the first order atmospheric
                extinction coefficient for filter f.
            Z_f : float
                The zero point magnitude for filter f.
            Z_f_sigma : float
                The standard deviation of the zero point magnitude
                for filter f.
    instr_filter : string
        The instrumental filter used to calculate C'.

    Returns
    -------
    c_prime_fci : float
        The C' coefficient for filter instr_filter to use when
        applying the transforms.

    """
    colour_index, ci = get_colour_index_lower(instr_filter)
    c_numerator = calculate_c_fci(
        gb_final_transforms, instr_filter, airmass, colour_index)
    c_negative = calculate_c_fci(
        gb_final_transforms, ci[0], airmass, colour_index)
    c_positive = calculate_c_fci(
        gb_final_transforms, ci[1], airmass, colour_index)
    c_prime_fci = c_numerator / (1 - c_negative + c_positive)
    return c_prime_fci


def calculate_z_prime_f(gb_final_transforms,
                        instr_filter,
                        airmass,
                        colour_index):
    """
    Calculate the Z' coefficient to use when applying the transforms.

    Parameters
    ----------
    gb_final_transforms : astropy.table.Table
        Table containing the results for the final transforms. Has columns:
            filter : string
                Instrumental filter band used to calculate the transform.
            CI : string
                Name of the colour index used to calculate the transform
                (e.g. B-V for b, V-R for r).
            k''_fCI : float
                The second order atmospheric extinction coefficient for
                filter f using the colour index CI.
            k''_fCI_sigma : float
                The standard deviation of the second order atmospheric
                extinction coefficient for filter f using the 
                colour index CI.
            T_fCI : float
                The instrumental transform coefficient for filter f using
                the colour index CI.
            T_fCI_sigma : float
                The standard deviation of the instrumental transform
                coefficient for filter f using the colour index CI.
            k'_f : float
                The first order atmospheric extinction coefficient
                for filter f.
            k'_f_sigma : float
                The standard deviation of the first order atmospheric
                extinction coefficient for filter f.
            Z_f : float
                The zero point magnitude for filter f.
            Z_f_sigma : float
                The standard deviation of the zero point magnitude
                for filter f.
    instr_filter : string
        The instrumental filter to use to calculate Z'_f.
    airmass : float
        Average airmass of the image.
    colour_index : string
        The colour index to use to calculate Z'_f.

    Returns
    -------
    z_prime_f : float
        Z' coefficient to use when calculating z.

    """
    mask = ((gb_final_transforms['filter'] == instr_filter) & (
        gb_final_transforms['CI'] == colour_index))
    row_of_transforms = gb_final_transforms[mask]
    if len(row_of_transforms) == 0 and instr_filter == 'v':
        instr_filter = 'g'
        mask = ((gb_final_transforms['filter'] == instr_filter) & (
            gb_final_transforms['CI'] == colour_index))
        row_of_transforms = gb_final_transforms[mask]
    if len(row_of_transforms) == 0 and instr_filter == 'b':
        instr_filter = 'u'
        mask = ((gb_final_transforms['filter'] == instr_filter) & (
            gb_final_transforms['CI'] == colour_index))
        row_of_transforms = gb_final_transforms[mask]
    z_prime_f = float(row_of_transforms['Z_f']) - \
        (float(row_of_transforms["k'_f"]) * airmass)
    return z_prime_f


def calculate_lower_z_f(gb_final_transforms,
                        c_prime_fci,
                        instr_filter,
                        airmass):
    """
    Calculate the z_f coeffient to use when applying the transforms.

    Parameters
    ----------
    gb_final_transforms : astropy.table.Table
        Table containing the results for the final transforms. Has columns:
            filter : string
                Instrumental filter band used to calculate the transform.
            CI : string
                Name of the colour index used to calculate the transform
                (e.g. B-V for b, V-R for r).
            k''_fCI : float
                The second order atmospheric extinction coefficient for
                filter f using the colour index CI.
            k''_fCI_sigma : float
                The standard deviation of the second order atmospheric
                extinction coefficient for filter f using the 
                colour index CI.
            T_fCI : float
                The instrumental transform coefficient for filter f
                using the colour index CI.
            T_fCI_sigma : float
                The standard deviation of the instrumental transform
                coefficient for filter f using the colour index CI.
            k'_f : float
                The first order atmospheric extinction coefficient for
                filter f.
            k'_f_sigma : float
                The standard deviation of the first order atmospheric
                extinction coefficient for filter f.
            Z_f : float
                The zero point magnitude for filter f.
            Z_f_sigma : float
                The standard deviation of the zero point magnitude for
                filter f.
    c_prime_fci : float
        The C' coefficient for filter instr_filter to use when applying
        the transforms.
    instr_filter : string
        The instrumental filter to use to calculate z_f.
    airmass : float
        Average airmass of the image.

    Returns
    -------
    lower_z_f : float
        The z coefficient for filter instr_filter to use when applying
        the transforms.

    """
    colour_index, ci = get_colour_index_lower(instr_filter)
    c_prime_fci = calculate_c_prime(gb_final_transforms, instr_filter, airmass)
    z_prime_positive = calculate_z_prime_f(
        gb_final_transforms, ci[0], airmass, colour_index)
    z_prime_negative = calculate_z_prime_f(
        gb_final_transforms, ci[1], airmass, colour_index)
    z_prime_no_brackets = calculate_z_prime_f(
        gb_final_transforms, instr_filter, airmass, colour_index)
    lower_z_f = c_prime_fci * \
        (z_prime_positive - z_prime_negative) + z_prime_no_brackets
    return lower_z_f


def apply_gb_transforms_VERIFICATION(gb_final_transforms,
                                     stars_table,
                                     instr_filter):
    """
    Apply the transforms to a source with an unknown standard magnitude.

    Parameters
    ----------
    gb_final_transforms : astropy.table.Table
        Table containing the results for the final transforms. Has columns:
            filter : string
                Instrumental filter band used to calculate the transform.
            CI : string
                Name of the colour index used to calculate the transform
                (e.g. B-V for b, V-R for r).
            k''_fCI : float
                The second order atmospheric extinction coefficient for
                filter f using the colour index CI.
            T_fCI : float
                The instrumental transform coefficient for filter f using
                the colour index CI.
            k'_f : float
                The first order atmospheric extinction coefficient for
                filter f.
            Z_f : float
                The zero point magnitude for filter f.
    unknown_object_table : astropy.table.Table
        Table containing the stars to calculate the apparent magnitudes for.
        Has columns:
            filter : string
                Instrumental filter band used to take the image.
            name : string
                Unique identifier of the source to apply the transform to.
            instrumental mag : float
                Instrumental magnitude of the source.
                : Change this so that it accepts multiple sources somehow.

    Returns
    -------
    app_mag_table : astropy.table.Table
        Table containing the apparent magnitudes of the object after applying
        the transforms. Has columns:
            time : float
                Julian date that the image was taken.
            filter : string
                Apparent filter that the image was transformed to.
            apparent mag : float
                Apparent magnitude of the source after applying the transforms.

    """
    # If there is only 1 entry per filter per source, assume that they are averages and don't need any interpolation.
    app_mag_first_columns = Table(
        stars_table['Field',
                    'Name',
                    'V_ref',
                    'B-V',
                    'U-B',
                    'V-R',
                    'V-I',
                    'V_sigma'])
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
    lower_z_f = calculate_lower_z_f(
        gb_final_transforms, c_prime_fci, instr_filter, airmass)
    app_mag_list = instr_mag + c_prime_fci * \
        (positive_instr_mag - negative_instr_mag) + lower_z_f
    app_mag_filter = instr_filter.upper()
    instr_filter_sigma = f"{instr_filter}_sigma"
    app_filter_sigma = f"{app_mag_filter}_sigma"
    app_mag_column = Table(names=[app_mag_filter], data=[app_mag_list])
    app_mag_sigma_column = Table(names=[app_filter_sigma], data=[
                                 stars_table[instr_filter_sigma]])
    app_mag_table = table.hstack(
        [app_mag_first_columns, app_mag_column, app_mag_sigma_column])
    return app_mag_table


def apply_gb_timeseries_transforms(gb_final_transforms,
                                   sat_dict,
                                   sat_auxiliary_table,
                                   unique_filters):
    """
    Apply the transforms to a timeseries of observations.

    TODO.

    Parameters
    ----------
    gb_final_transforms : astropy.table.Table
        Table containing the results for the final transforms. Has columns:
            filter : string
                Instrumental filter band used to calculate the transform.
            CI : string
                Name of the colour index used to calculate the transform 
                (e.g. B-V for b, V-R for r).
            k''_fCI : float
                The second order atmospheric extinction coefficient for filter
                f using the colour index CI.
            T_fCI : float
                The instrumental transform coefficient for filter f using the
                colour index CI.
            k'_f : float
                The first order atmospheric extinction coefficient for
                filter f.
            Z_f : float
                The zero point magnitude for filter f.
    large_sats_table : astropy.table.Table
        DESCRIPTION.

    Returns
    -------
    app_large_sats_table : TYPE
        DESCRIPTION.

    """
    app_sat_dict = dict.fromkeys(sat_dict.keys())
    for sat, sat_table in sat_dict.items():
        sat_table_columns = sat_table.columns
        num_columns = len(sat_table_columns)
        app_table_data = np.empty((len(sat_table), num_columns))
        app_sat_table = Table(names=sat_table_columns, data=app_table_data)
        app_sat_table['Time (JD)'] = sat_table['Time (JD)']
        for unique_filter in unique_filters:
            instr_filter = unique_filter.lower()
            try:
                app_sat_table[f'{unique_filter}_sigma'] =\
                    sat_table[f'{unique_filter}_sigma']
            except KeyError:
                pass
            colour_index, ci = get_colour_index_lower(instr_filter)
            CI = ci.upper()
            instr_mag = sat_table[unique_filter]
            airmass = sat_auxiliary_table['Airmass']
            c_prime_fci = calculate_c_prime(
                gb_final_transforms, instr_filter, airmass)
            try:
                positive_instr_mag = sat_table[CI[0]]
                negative_instr_mag = sat_table[CI[1]]
            except KeyError:
                # table_ci = CI.replace('V', 'G')
                if 'v' in ci:
                    table_ci = ci.replace('v', 'g')
                else:
                    table_ci = ci
                if 'b' in ci:
                    table_ci = table_ci.replace('b', 'u')
                positive_instr_mag = sat_table[table_ci[0]]
                negative_instr_mag = sat_table[table_ci[1]]
            lower_z_f = calculate_lower_z_f(
                gb_final_transforms, c_prime_fci, instr_filter, airmass)
            app_mag_list = instr_mag + c_prime_fci * \
                (positive_instr_mag - negative_instr_mag) + lower_z_f
            app_sat_table[unique_filter] = app_mag_list
        app_sat_dict[sat] = app_sat_table
    return app_sat_dict


def copy_and_rename(directory,
                    file_suffix=(".fits", ".fit", ".fts"),
                    time_key='DATE-OBS',
                    filter_key='FILTER',
                    temp_dir='tmp',
                    debugging=False):
    """
    Copy and rename all files in a directory to a new location.

    Parameters
    ----------
    directory : string
        Location of the files to be copied and modified.
    file_suffix : string, optional
        Considers all files that end with this suffix. The default is ".fits".
    time_key : string, optional
        Key in the FITS header that contains the time of the observation.
        The default is 'DATE-OBS'.
    filter_key : string, optional
        Key in the FITS header that contains the filter used for the
        observation. The default is 'FILTER'.
    temp_dir : string, optional
        Location to copy the files to. The default is 'tmp'.
    debugging : bool, optional
        Keyword specifying if the code is being run in a debugging mode.
        If True, it removes temp_dir if it exists. If 
        False, it tries to make temp_dir and uses the current one
        if it exists. The default is False.

    Returns
    -------
    filecount : int
        Number of files that were modified.
    filenames : string list
        List of all of the modified files sorted alphabetically.

    """
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
            if file.endswith((file_suffix)):
                image = fits.open(os.path.join(dirpth, file))
                hdr = image[0].header
                image.close()
                # with fits.open(os.path.join(dirpth, file)) as image:
                #     hdr = image[0].header
                t = Time(hdr[time_key], format='fits', scale='utc')
                try:
                    filter_name = hdr[filter_key]
                except KeyError:
                    filter_name = 'C'
                t_datetime = t.to_datetime()
                new_filename = f'{t_datetime.strftime("%Y%m%d%H%M%S")}{filter_name}.fits'
                copy2(os.path.join(dirpth, file), f'{temp_dir}/{new_filename}')
                filecount += 1
    filenames = sorted(os.listdir(temp_dir))
    return filecount, filenames


def remove_temp_dir(temp_dir='tmp'):
    """
    Remove the temp directory.

    Parameters
    ----------
    temp_dir : string, optional
        Directory to remove. The default is 'tmp'.

    Returns
    -------
    None.

    """
    rmtree(temp_dir)


def TRM_sat_detection(filepath,
                      ecct_cut=0.5,
                      sigma_clip=2.5,
                      edge_protect=10,
                      SNRLimit=0,
                      min_obj_pixels=5,
                      pix_frac=0,
                      moffat_avg=0,
                      gauss_avg=0,
                      star_count=0,
                      mstar_count=0,
                      count=0):
    def BackgroundIteration(image, tolerance):
        # old_mean = 1e9
        old_rms = 1e9

        new_mean = 2e9
        new_rms = 2e9

        while abs(new_rms - old_rms) > (tolerance * old_rms):
            # old_mean = float(new_mean)
            old_rms = float(new_rms)
            # image = myclip(image, (old_mean - 2 * old_rms),
            # (old_mean + 2 * old_rms))

            if (np.size(image) == 0):
                new_mean = 0
                new_rms = 2e9
                break

            new_mean = np.mean(image)
            new_rms = np.std(image)
            # retval = [new_mean, new_rms]
            return new_mean, new_rms

    def myclip(x1, lo, hi):
        vector = np.vectorize(np.float)
        x = vector(x1)

        float(hi)
        float(lo)
        # print(x)
        # print(hi)
        # print(lo)

        y = (x * np.any(x <= hi)) + ((hi) * np.any(x > hi))
        y = (y * np.any(y > lo)) + ((lo) * np.any(y <= lo))
        return y

    def PointSourceFluxExtraction(mask_x, mask_y, flux_image):

        num_elem_x = mask_x.size
        # num_elem_y = mask_y.size
        sum1 = 0
        pix_flux = np.zeros((num_elem_x))
        for i in range(num_elem_x):
            x_pix = mask_x[i]
            y_pix = mask_y[i]

            sum1 = sum1 + flux_image[x_pix, y_pix]
            pix_flux[i] = flux_image[x_pix, y_pix]

        object_flux = sum1
        max_pixel_flux = max(pix_flux)
        return object_flux, max_pixel_flux

    def MomentCalculation(xmask, ymask, xc, yc, p, q):
        num_pix = xmask.size
        mom = sum((xmask - xc) ** p * (ymask - yc) ** q) / num_pix
        moment = mom
        return moment

    def EccentricityCalculation(m11, m02, m20):
        eccent = np.sqrt((m20 - m02) ** 2 + (4 * m11 ** 2)) / (m20 + m02)
        return eccent

    def Compact(num_pix, m02, m20):
        compact = (num_pix / (m02 + m20))
        return compact

    def WeightedCentroid(mask_x, mask_y, flux_image):

        num_elem_x = mask_x.size
        num_elem_y = mask_y.size
        x_wt_sum = 0
        y_wt_sum = 0
        flux_sum = 0
        # print("2")
        if num_elem_x != num_elem_y:
            # object_flux = -999
            # print("3")
            return
        else:
            for i in range(num_elem_x):
                x_pix = mask_x[i]
                y_pix = mask_y[i]

                x_wt_sum = x_wt_sum + (x_pix * flux_image[x_pix, y_pix])
                y_wt_sum = y_wt_sum + (y_pix * flux_image[x_pix, y_pix])
                flux_sum = flux_sum + flux_image[x_pix, y_pix]

        x_centroid = x_wt_sum / flux_sum
        y_centroid = y_wt_sum / flux_sum

        x_var_sum = 0
        y_var_sum = 0
        flux_sum = 0
        # print("2")
        for i in range(num_elem_x):
            x_pix = mask_x[i]
            y_pix = mask_y[i]

            x_var_sum = x_var_sum + \
                ((x_pix - x_centroid) ** 2 * flux_image[x_pix, y_pix])
            y_var_sum = y_var_sum + \
                ((y_pix - y_centroid) ** 2 * flux_image[x_pix, y_pix])
            flux_sum = flux_sum + flux_image[x_pix, y_pix]

        x_rms = np.sqrt(x_var_sum / flux_sum)
        y_rms = np.sqrt(y_var_sum / flux_sum)
        return x_centroid, x_rms, y_centroid, y_rms

    # streak1 = r'D:\Transfer to mac\2021-03-10 - Calibrated\Intelsat 10-02
    # Post Eclipse\LIGHT\B_lim\0066_3x3_-10.00_5.00_B_21-23-04.fits'
    # streak = \
    # 'D:\\Breeze-M_R_B_38746U\\CAN_OTT.00018674.BREEZE-M_R_B_#38746U.FIT'
    # streak12 =\
    #  r'D:\Transfer to mac\trm-stars-images\NEOS_SCI_2021099173229frame.fits'
    hdr, fitsdata = read_fits_file(filepath)
    # STARS = open(filepath, "w")
    # imagehdularray = fits.open(filepath)

    streak_array = []
    # sigma_clip = 2.5
    # edge_protect = 10
    # min_obj_pixels = 5
    # SNRLimit = 0
    # pix_frac = 0;
    # moffat_avg = 0;
    # gauss_avg = 0;
    # star_count = 0;
    # mstar_count = 0;
    # count =0

    # date=imagehdularray[0].header['DATE-OBS']
    # exposuretime=imagehdularray[0].header['EXPTIME']
    # imagesizeX=imagehdularray[0].header['NAXIS1']
    # imagesizeY=imagehdularray[0].header['NAXIS2']
    imagesizeX = hdr['NAXIS1']
    imagesizeY = hdr['NAXIS2']
    exposuretime = hdr['EXPTIME']
    # fitsdata =  imagehdularray[0].data
    sigma_clip = SigmaClip(sigma=2.5)
    # bkg = SExtractorBackground(sigma_clip)
    # bkg_value = bkg.calc_background(fitsdata)

    # print(bkg_value)
    # bkg = MeanBackground(sigma_clip)
    # bkg_value = bkg.calc_background(fitsdata)
    # bkg_estimator1 = SExtractorBackground()
    bkg_estimator2 = SExtractorBackground()
    # bkg = Background2D(fitsdata, (2, 2),
    # filter_size=(3,3),sigma_clip=sigma_clip,
    # bkg_estimator=bkg_estimator2) Closest Approximate to Matlab Result
    bkg = Background2D(fitsdata, (50, 50), filter_size=(
        3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2)
    bg_rem = fitsdata - bkg.background

    # print(np.mean(bkg.background))

    bg_rem[1:edge_protect, 1:edge_protect] = 0
    bg_rem[imagesizeX - edge_protect:imagesizeX, :] = 0
    bg_rem[:, 1:edge_protect] = 0
    bg_rem[:, imagesizeY - edge_protect:imagesizeY] = 0
    im_mean = np.mean(bg_rem)

    im_rms = np.std(fitsdata)
    im_mean, im_rms = BackgroundIteration(bg_rem, 0.1)
    low_clip = im_mean + 2.5 * im_rms
    # high_clip = 161

    # binary_image = np.zeros((imagesizeX,imagesizeY))

    bg_rem[bg_rem <= low_clip]
    # binary_image = (binary_image * bg_rem[bg_rem<= low_clip]) +
    # (1 * bg_rem[bg_rem> low_clip])

    th, im_th = cv.threshold(bg_rem, low_clip, 1, cv.THRESH_BINARY)
    # print(im_mean)
    connected_image = measure.label(im_th, background=0)
    # plt.subplot(133)
    # plt.imshow(connected_image, cmap='nipy_spectral')
    # plt.axis('off')
    # plt.tight_layout()
    # plt.show()
    # im = cv2.imread(bg_rem)
    # th, im_th = cv2.threshold(im, 128, 255, cv2.THRESH_BINARY)

    # num_labels, labels_im = cv2.connectedComponents(im_th)
    # num_sourcepix = cv2.connectedComponentsWithStats(binary_image,
    # np.array(connected_image), np.array(stats), np.array(centroids),
    # 4, np.int(CV_32S))
    num_sourcepix = np.zeros(shape=(100000, 1))
    [size_x, size_y] = imagesizeX, imagesizeY

    for x in range(0, size_y):
        for y in range(0, size_x):
            pixval = connected_image[x, y]

            if (pixval != 0):
                num_sourcepix[pixval, 0] = num_sourcepix[pixval, 0] + 1

    [valid_sources, temp] = np.nonzero(num_sourcepix > min_obj_pixels)
    num_valid_sources = valid_sources.size
    if num_valid_sources == 0:
        return

    centroid_x = np.zeros((num_valid_sources, 1))
    centroid_y = np.zeros((num_valid_sources, 1))
    rms_x_pos = np.zeros((num_valid_sources, 1))
    rms_y_pos = np.zeros((num_valid_sources, 1))
    m11 = np.zeros((num_valid_sources, 1))
    m02 = np.zeros((num_valid_sources, 1))
    m20 = np.zeros((num_valid_sources, 1))
    ecct = np.zeros((num_valid_sources, 1))
    compact = np.zeros((num_valid_sources, 1))
    obj_flux = np.zeros((num_valid_sources, 1))
    obj_max1 = np.zeros((num_valid_sources, 1))
    length = np.zeros((num_valid_sources, 1))

    for j in range(num_valid_sources):

        vsj = valid_sources[j]
        a = [0, 0]

        [mask_x, mask_y] = np.nonzero(connected_image == vsj)
        obj_flux[j], obj_max1[j] = PointSourceFluxExtraction(
            mask_x, mask_y, bg_rem)

        centroid_x[j] = np.mean(mask_x)
        rms_x_pos[j] = np.std(mask_x)
        centroid_y[j] = np.mean(mask_y)
        rms_y_pos[j] = np.std(mask_y)

        m11[j] = MomentCalculation(
            mask_x, mask_y, centroid_x[j], centroid_y[j], 1, 1)
        m02[j] = MomentCalculation(
            mask_x, mask_y, centroid_x[j], centroid_y[j], 0, 2)
        m20[j] = MomentCalculation(
            mask_x, mask_y, centroid_x[j], centroid_y[j], 2, 0)
        compact[j] = Compact(num_sourcepix[vsj], m02[j], m20[j])
        ecct[j] = EccentricityCalculation(m11[j], m02[j], m20[j])

        x_length = (max(mask_x) - min(mask_x))
        y_length = (max(mask_y) - min(mask_y))
        length[j] = np.sqrt(x_length ** 2 + y_length ** 2)

        Zp = 21
        vmag = Zp - 2.5 * np.log10(obj_flux[j] / exposuretime)

        if obj_max1[j] < 60000:  # %fit unsaturated stars only
            if (centroid_x[j] > 10) and (centroid_x[j] < (imagesizeY - 10)) \
                and (centroid_y[j] > 10) and (
                    centroid_y[j] < (imagesizeX - 10)):
                # %Find middle pixel value

                [cen_x, rms_x, cen_y, rms_y] = WeightedCentroid(
                    mask_x, mask_y, 0 * bg_rem + 1)

                if (centroid_x[j] > 10) and\
                    (centroid_x[j] < (imagesizeX - 10)) and\
                        (centroid_y[j] > 10) and\
                        (centroid_y[j] < (imagesizeY - 10)):
                    # mid_pix_val = bg_rem(round(cen_x),round(cen_y))
                    cenx = int(centroid_y[j])
                    ceny = int(centroid_x[j])
                    mid_pix_valPP = bg_rem[ceny, cenx]

                    if vmag < 13:
                        # Fit a moffat profile
                        # holds radial distance from centroid
                        r = np.zeros(len(mask_x))
                        S = np.zeros(len(mask_x))  # holds intensity
                        np.delete(S, -1)
                        for q in range(0, len(mask_x)):
                            r[q] = np.sqrt(
                                (mask_x[q] + 0.5 - (ceny + 1)) ** 2 +
                                (mask_y[q] + 0.5 - (cenx + 1)) ** 2)
                            S[q] = bg_rem[mask_x[q], mask_y[q]]

                        C_index = np.argmin(r)
                        r[C_index] = 0  # %centroid radial value
                        C = S[C_index]
                        # %a holds [alpha Beta] moffat parameters
                        # %Fix a(2) Beta parameter to 1.5
                        # a = [0,0]
                        # print(C)

                        def fun(a): return sum(
                            (S - (C / ((1 + (r ** 2) / (a[0] ** 2)) ** 1.5)))
                            ** 2)
                        aguess = 1
                        a = scipy.optimize.fmin(func=fun, x0=aguess, disp=0)
                        # print(a)

                        # %b holds [alpha Beta] moffat parameters

                        def fung(b): return sum(
                            (S - (C * np.exp(-(r ** 2) / (2 * (b ** 2))))) **
                            2)
                        bguess = 2
                        b = scipy.optimize.fmin(func=fung, x0=bguess, disp=0)
                        # print(b)
                        # %Optional plot the fits:

                        # plt.scatter(r,S);
                        # E = lambda a,r: (C/((1+(r**2)/(a[0]^2))**1.5))
                        # F = lambda b,r:(C*np.exp(-(r**2)/(2*(b**2))))
                        # plot=plt(E,[0,max(r)])

                        # h = plt.gca().get_children()

                        # plot.set(h(1),'color','red')

                        # plot= plt(F,[0,max(r)])
                        # plt.axis([0,max(r),0,60000])

                        # h = plt.gca().get_children()
                        # plot.set(h(1),'color','green')
                        # Output results

                    else:
                        a = [0, 0]
                        b = 0

                    pix_frac = pix_frac + mid_pix_valPP / obj_flux[j]

                    if vmag < 13 and a[0] < 4:
                        # mstar_count = mstar_count +1;
                        # print(a[0])
                        count = count + 1
                        moffat_avg = moffat_avg + a[0]
                        gauss_avg = gauss_avg + b

    # [compact_mean, compact_rms] = BackgroundIteration(compact,0.1)
    [ecct_mean, ecct_rms] = BackgroundIteration(ecct, 0.1)
    # compact_cut = compact_mean  + 1 * compact_rms
    # ecct_cut = 0.5

    stars = np.nonzero(ecct < ecct_cut)
    streaks = np.nonzero(ecct > ecct_cut)
    stars = np.delete(stars, 1, 0)
    streaks = np.delete(streaks, 1, 0)
    sda = valid_sources[stars]
    if len(sda) == 0 or count == 0:
        return
    num_pix_in_stars = num_sourcepix[sda]
    [mean_starpix, rms_starpix] = BackgroundIteration(num_pix_in_stars, 0.1)

    sat_x = centroid_y[stars].flatten()
    sat_y = centroid_x[stars].flatten()
    # sat_array = np.empty((len(stars[0]), 2))
    # for sat_index in range(len(stars[0])):
    #     sat_array[sat_index][0] = float(sat_x[0][sat_index])
    #     sat_array[sat_index][1] = float(sat_y[0][sat_index])

    pix_cutoff = mean_starpix + 10 * rms_starpix

    num_stars = stars.size

    stellar_flux_SNR = np.zeros((num_valid_sources, 1))

    xmin = edge_protect
    xmax = imagesizeX - edge_protect
    ymin = edge_protect
    ymax = imagesizeY - edge_protect
    streaksize = streaks.size

    for k in range(streaksize):

        real_star_num = streaks[0, k]
        vsj = valid_sources[real_star_num]
        [mask_x, mask_y] = np.nonzero(connected_image == vsj)

        [cen_x, rms_x, cen_y, rms_y] = WeightedCentroid(mask_x, mask_y, bg_rem)

        temp_SNR = obj_max1[real_star_num] / im_rms
        stellar_flux_SNR[k] = temp_SNR
        # print(temp_SNR)
        if temp_SNR > SNRLimit:
            if ((cen_x > xmin) & (cen_x < xmax) & (cen_y > ymin) &
                    (cen_y < ymax)):
                stellar_flux_SNR[k] = temp_SNR

                if streak_array == []:
                    streak_arrayelement = [cen_x, rms_x, cen_y, rms_y,
                                           obj_flux[real_star_num],
                                           stellar_flux_SNR[k],
                                           exposuretime]
                    streak_array.append(streak_arrayelement)
                    # flux=float(obj_flux[real_star_num,0])
                    # streak_line='{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(cen_y), float(cen_x),  flux)

                    # STARS.write(streak_line+"\n")
                    # print(Streaks_Detected,
                    # [num2str(cen_x) ','
                    # num2str(rms_x) ','
                    # num2str(cen_y) ','
                    # num2str(rms_y) ','
                    # num2str(obj_max1(1,rsn))
                    # ',' num2str(temp_SNR)
                    # ',' fpath1(i).name '\r\n'])
                else:
                    new_element = [cen_x,
                                   rms_x,
                                   cen_y,
                                   rms_y,
                                   obj_flux[real_star_num],
                                   stellar_flux_SNR[k],
                                   exposuretime]
                    # flux=float(obj_flux[real_star_num,0])
                    # streak_line='{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(cen_y), float(cen_x),  flux)
                    # STARS.write(streak_line+"\n")
                    streak_array.append(new_element)

    # avg_pix_frac = pix_frac/star_count
    moffat_avg = moffat_avg / count
    gauss_avg = gauss_avg / count
    FWHM = float(2 * gauss_avg * 0.7664)

    # STARS.close()
    return sat_x, sat_y, bkg.background, FWHM


def set_sat_positions(imgdata,
                      filecount,
                      set_sat_positions_bool,
                      max_distance_from_sat=25,
                      norm=LogNorm(),
                      cmap_set='Set1'):
    """
    Initialize the names and locations of the satellites to create a light
    curve for.

    Parameters
    ----------
    imgdata : numpy.ndarray
        Data from the fits file.
    filecount : float
        Number of files that will be used to create the light curve.
    set_sat_positions_bool : bool
        Decides whether or not to initialize the satellite positions.
    max_distance_from_sat : int, optional
        Maximimum number of pixels that a source can be away from the defined
        sat position to be considered the sat. 
        The default is 25.
    norm : TYPE, optional
        DESCRIPTION. The default is LogNorm(). #TODO
    cmap_set : string, optional
        CMAP set to use for plotting the satellite positions. The default is
        'Set1'.

    Returns
    -------
    TYPE
        DESCRIPTION.
        TODO

    """

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

    def onclick(event, sat_locs):
        # global num_sat_event
        sat_locs.append([event.xdata, event.ydata])

    # global set_sat_positions_bool
    while set_sat_positions_bool:
        # global sat_locs
        sat_locs = []
        # global num_sat_event
        # num_sat_event = 0
        # return satloc
        # mbox('Information',
        #       'Please select the positions of the satellites on the
        # following image. Press any key when finished.',
        #       0)
        # cv.namedWindow('TestImage')
        # cv.setMouseCallback('TestImage', set_sat_position)
        # logdata = cv.normalize(imgdata, None, alpha=0, beta=5,
        # norm_type=cv.NORM_MINMAX, dtype=cv.CV_32F)
        # cv.imshow('TestImage', logdata)
        # cv.waitKey(0)
        # cv.destroyAllWindows()
        try:
            fig = Figure()
            ax = fig.add_subplot()
            root = tk.Tk()
            ax.imshow(imgdata, cmap='gray', norm=norm, interpolation='nearest')
            canvas = FigureCanvasTkAgg(fig, master=root)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
            canvas.mpl_connect('button_press_event',
                               lambda event: onclick(event, sat_locs))
            root.update()
            root.mainloop()
        except ValueError:
            fig = Figure()
            ax = fig.add_subplot()
            root = tk.Tk()
            ax.imshow(abs(imgdata), cmap='gray', norm=LogNorm(
                vmin=1), interpolation='nearest')
            canvas = FigureCanvasTkAgg(fig, master=root)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
            canvas.mpl_connect('button_press_event',
                               lambda event: onclick(event, sat_locs))
            root.update()
            root.mainloop()
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
        auxiliary_data = np.empty((filecount, 3))
        auxiliary_data.fill(np.nan)

        for i, name in enumerate(names[2:]):
            try:
                fig = Figure()
                ax = fig.add_subplot()
                ax.imshow(imgdata, cmap='gray', norm=norm,
                          interpolation='nearest')
                sat_aperture = RectangularAperture(sat_locs[i],
                                                   w=max_distance_from_sat * 2,
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
                button = tk.Button(input_frame, text="OK",
                                   command=return_entry)
                button.pack()
                canvas = FigureCanvasTkAgg(fig, master=img_frame)
                canvas.draw()
                toolbar = NavigationToolbar2Tk(canvas, img_frame)
                toolbar.update()
                canvas.get_tk_widget().pack(side=tk.LEFT,
                                            fill=tk.BOTH,
                                            expand=1)
                root.update()
                root.focus_force()
                entry.focus_set()
                root.mainloop()
            except ValueError:
                fig = Figure()
                ax = fig.add_subplot()
                ax.imshow(abs(imgdata), cmap='gray', norm=LogNorm(
                    vmin=1), interpolation='nearest')
                sat_aperture =\
                    RectangularAperture(sat_locs[i],
                                        w=max_distance_from_sat * 2,
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
                button = tk.Button(input_frame, text="OK",
                                   command=return_entry)
                button.pack()
                canvas = FigureCanvasTkAgg(fig, master=img_frame)
                canvas.draw()
                toolbar = NavigationToolbar2Tk(canvas, img_frame)
                toolbar.update()
                canvas.get_tk_widget().pack(side=tk.LEFT,
                                            fill=tk.BOTH,
                                            expand=1)
                root.update()
                root.focus_force()
                entry.focus_set()
                root.mainloop()
            names[i + 2] = content
            print(
                f"Satellite {names[i + 2]} at location ({sat_locs[i, 0]}, {sat_locs[i, 1]})")
        print(names)

        cmap = plt.get_cmap(cmap_set)
        colours = [cmap(i) for i in range(0, num_sats)]
        legend_elements = []
        window = tk.Tk()
        window.title('Plotting in Tkinter Test')
        img_frame = tk.Frame(window)
        img_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        input_frame = tk.Frame(window)
        input_frame.pack(side=tk.RIGHT)
        label = tk.Label(
            input_frame, text='Are the satellite positions correct?')
        label.pack()
        yes_no = tk.IntVar()
        yes_btn = tk.Radiobutton(
            input_frame, text='Yes', variable=yes_no, value=1)
        yes_btn.pack(anchor=tk.W, padx=5)
        no_btn = tk.Radiobutton(input_frame, text='No',
                                variable=yes_no, value=2)
        no_btn.pack(anchor=tk.W, padx=5)
        closebutton = tk.Button(input_frame, text='OK', command=window.destroy)
        closebutton.pack()
        try:
            fig = Figure()
            ax = fig.add_subplot()
            ax.imshow(imgdata, cmap='gray', norm=LogNorm(),
                      interpolation='nearest')
            for i in range(0, num_sats):
                sat_aperture = RectangularAperture(sat_locs[i],
                                                   w=max_distance_from_sat * 2,
                                                   h=max_distance_from_sat * 2)
                sat_aperture.plot(axes=ax, color=colours[i], lw=1.5, alpha=0.5)
                legend_elements.append(Line2D([0],
                                              [0],
                                              color='w',
                                              marker='s',
                                              markerfacecolor=colours[i],
                                              markersize=7,
                                              label=names[i + 2]))
            fig.legend(handles=legend_elements, framealpha=1)
            canvas = FigureCanvasTkAgg(fig, master=img_frame)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, img_frame)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            window.mainloop()
        except ValueError:
            fig = Figure()
            ax = fig.add_subplot()
            ax.imshow(abs(imgdata), cmap='gray', norm=LogNorm(
                vmin=1), interpolation='nearest')
            for i in range(0, num_sats):
                sat_aperture = RectangularAperture(sat_locs[i],
                                                   w=max_distance_from_sat * 2,
                                                   h=max_distance_from_sat * 2)
                sat_aperture.plot(axes=ax, color=colours[i], lw=1.5, alpha=0.5)
                legend_elements.append(Line2D([0],
                                              [0],
                                              color='w',
                                              marker='s',
                                              markerfacecolor=colours[i],
                                              markersize=7,
                                              label=names[i + 2]))
            fig.legend(handles=legend_elements, framealpha=1)
            canvas = FigureCanvasTkAgg(fig, master=img_frame)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, img_frame)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            window.mainloop()
        # Have a way for the user to confirm the satellite locations.
        # If it is wrong, then decide whether to change
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
        auxiliary_column_names = ["FWHM", "Airmass", "BSB"]
        ausiliary_columns = Table(
            names=auxiliary_column_names, data=auxiliary_data)
        sats_table = hstack(
            [date_table, filter_table, data_table], join_type='exact')
        uncertainty_table = hstack(
            [date_table, filter_table, data_table], join_type='exact')
        sat_auxiliary_table = hstack(
            [date_table, filter_table, ausiliary_columns], join_type='exact')
        sats_table.pprint_all()
    sat_information = namedtuple('sat_information',
                                 ['sats_table',
                                  'uncertainty_table',
                                  'sat_auxiliary_table',
                                  'sat_locs',
                                  'num_sats',
                                  'num_nans',
                                  'sat_names'])
    return set_sat_positions_bool, sat_information(sats_table,
                                                   uncertainty_table,
                                                   sat_auxiliary_table,
                                                   sat_locs,
                                                   num_sats,
                                                   num_nans,
                                                   sat_names)
    # return sats_table, uncertainty_table, sat_fwhm_table, sat_locs,
    # num_sats, num_nans, sat_names


def plot_detected_sats(filename,
                       plot_results,
                       imgdata,
                       irafsources,
                       sat_information,
                       max_distance_from_sat=20,
                       norm=LogNorm()):
    """
    TODO.

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    plot_results : TYPE
        DESCRIPTION.
    imgdata : TYPE
        DESCRIPTION.
    irafsources : TYPE
        DESCRIPTION.
    sat_information : TYPE
        DESCRIPTION.
    max_distance_from_sat : TYPE, optional
        DESCRIPTION. The default is 20.
    norm : TYPE, optional
        DESCRIPTION. The default is LogNorm().

    Returns
    -------
    None.

    """
    if plot_results != 0:
        fig, ax = plt.subplots()
        ax.imshow(imgdata, cmap='gray', norm=norm, interpolation='nearest')
        ax.scatter(irafsources['xcentroid'], irafsources['ycentroid'],
                   s=100, edgecolor='red', facecolor='none')
        for i in range(0, sat_information.num_sats):
            rect = patches.Rectangle((sat_information.sat_locs[i, 0] -
                                      max_distance_from_sat,
                                      sat_information.sat_locs[i, 1] -
                                      max_distance_from_sat),
                                     width=max_distance_from_sat * 2,
                                     height=max_distance_from_sat * 2,
                                     edgecolor='green', facecolor='none')
            ax.add_patch(rect)
        plt.title(filename)
        if plot_results == 1:
            plt.show(block=False)
            plt.pause(2)
        elif plot_results == 2:
            plt.show()
        plt.close()


def get_image_airmass(hdr,
                      airmass_key='AIRMASS',
                      az_key='CENAZ',
                      alt_key='CENALT',
                      ra_key='RA',
                      dec_key='DEC',
                      unit=(u.deg, u.deg),
                      lat_key='SITELAT',
                      lon_key='SITELONG',
                      elev_key='SITEELEV',
                      time_key='DATE-OBS'):
    if airmass_key in hdr.keys():
        airmass = hdr[airmass_key]
    elif az_key and alt_key in hdr.keys():
        altaz = AltAz(az=hdr[az_key], alt=hdr[alt_key])
        airmass = altaz.secz
    elif ra_key and dec_key in hdr.keys():
        skypositions = SkyCoord(ra=hdr[ra_key], dec=hdr[dec_key], unit=unit)
        altaz = convert_ra_dec_to_alt_az(
            skypositions,
            hdr,
            lat_key=lat_key,
            lon_key=lon_key,
            elev_key=elev_key)
        airmass = altaz.secz
    else:
        return
    return airmass


def check_if_sat(sat_information,
                 index,
                 sat_x_array,
                 sat_y_array,
                 instr_mags,
                 instr_mags_sigma,
                 fwhm_arcsec,
                 airmass,
                 bsb,
                 max_distance_from_sat=20):
    for obj_index in range(len(sat_x_array)):
        obj_x = sat_x_array[obj_index]
        obj_y = sat_y_array[obj_index]
        for sat_num, sat in enumerate(sat_information.sat_locs, start=2):
            sat_x = sat[0]
            sat_y = sat[1]
            if abs(sat_x - obj_x) < max_distance_from_sat and\
                    abs(sat_y - obj_y) < max_distance_from_sat:
                sat_information.sats_table[index][sat_num] =\
                    instr_mags[obj_index]
                sat_information.uncertainty_table[index][sat_num] =\
                    instr_mags_sigma[obj_index]
                sat[0] = obj_x
                sat[1] = obj_y
    sat_information.sat_auxiliary_table[index][2] = fwhm_arcsec
    sat_information.sat_auxiliary_table[index][3] = airmass
    sat_information.sat_auxiliary_table[index][4] = bsb
    print(sat_information.sats_table[index])
    return sat_information


def determine_if_change_sat_positions(sat_information,
                                      filenum,
                                      change_sat_positions_bool,
                                      max_num_nan=5):
    sat_mags = np.array(list(sat_information.sats_table[filenum]))
    mask = np.isnan(sat_mags[2:].astype(float))
    sat_information.num_nans[mask] += 1
    sat_information.num_nans[~mask] = 0
    print(sat_information.num_nans)
    num_nan = max(sat_information.num_nans)
    # print(num_nan)
    if num_nan != 0 and (num_nan % max_num_nan) == 0:
        change_sat_positions_bool = True
    return change_sat_positions_bool, num_nan


def change_sat_positions(filenames,
                         filenum,
                         num_nan,
                         sat_information,
                         change_sat_positions_bool,
                         gb_final_transforms=None,
                         max_distance_from_sat=20,
                         size=25,
                         temp_dir='tmp',
                         cmap_set='Set1',
                         plot_results=0):
    # Change the position
    # Display filenames[filenum - num_nan]
    def mbox(title, text, style):
        return ctypes.windll.user32.MessageBoxW(0, text, title, style)

    def onclick(event, sat_locs, index):
        # global num_sat_event
        sat_locs[index] = [event.xdata, event.ydata]

    def change_sat_position(event, x, y, flags, params):
        global sat_locs
        if event == cv.EVENT_LBUTTONDOWN:
            sat_locs[params[0]] = [x, y]
            cv.destroyAllWindows()

    print(filenames[filenum - num_nan])
    # sat_checked = np.zeros(num_sats, dtype=IntVar())
    with fits.open(f"{temp_dir}/{filenames[filenum - num_nan]}",
                   memmap=False) as image:
        hdr = image[0].header
        imgdata = image[0].data
    # image = fits.open(f"{temp_dir}/{filenames[filenum - num_nan]}")
    # hdr = image[0].header
    # imgdata = image[0].data
    # image.close()
    root = tk.Tk()
    root.title("Current satellite positions")
    input_frame = tk.Frame(root)
    input_frame.pack(side=tk.RIGHT)
    img_frame = tk.Frame(root)
    img_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
    label = tk.Label(
        input_frame,
        text='Select the satellite(s) whose position you would like to change.')
    label.grid(row=0)
    sat_checked = []
    for sat_num, sat in enumerate(sat_information.sat_names):
        sat_checked.append(tk.IntVar())
        checkbutton = tk.Checkbutton(
            input_frame, text=sat, variable=sat_checked[sat_num])
        checkbutton.grid(row=sat_num + 1, sticky=tk.W, padx=5)
    none_select = tk.IntVar()
    checkbutton = tk.Checkbutton(
        input_frame, text="None", variable=none_select)
    checkbutton.grid(row=sat_information.num_sats + 2, sticky=tk.W, padx=5)
    closebutton = tk.Button(input_frame, text='OK', command=root.destroy)
    closebutton.grid(row=sat_information.num_sats + 3)
    try:
        legend_elements = []
        fig = Figure()
        ax = fig.add_subplot()
        ax.imshow(imgdata, cmap='gray', norm=LogNorm(
            vmin=1), interpolation='nearest')
        # TODO: move this into a separate function?
        cmap = plt.get_cmap(cmap_set)
        colours = [cmap(i) for i in range(0, sat_information.num_sats)]
        for i in range(0, sat_information.num_sats):
            sat_aperture = RectangularAperture(sat_information.sat_locs[i],
                                               w=max_distance_from_sat * 2,
                                               h=max_distance_from_sat * 2)
            sat_aperture.plot(axes=ax, color=colours[i], lw=1.5, alpha=0.5)
            legend_elements.append(Line2D([0], [0], color='w',
                                          marker='s',
                                          markerfacecolor=colours[i],
                                          markersize=7,
                                          label=sat_information.sat_names[i]))
        fig.legend(handles=legend_elements, framealpha=1)
        canvas = FigureCanvasTkAgg(fig, master=img_frame)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, img_frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        root.mainloop()
    except ValueError:
        legend_elements = []
        fig = Figure()
        ax = fig.add_subplot()
        ax.imshow(abs(imgdata), cmap='gray', norm=LogNorm(
            vmin=1), interpolation='nearest')
        # TODO: move this into a separate function?
        cmap = plt.get_cmap(cmap_set)
        colours = [cmap(i) for i in range(0, sat_information.num_sats)]
        for i in range(0, sat_information.num_sats):
            sat_aperture = RectangularAperture(sat_information.sat_locs[i],
                                               w=max_distance_from_sat * 2,
                                               h=max_distance_from_sat * 2)
            sat_aperture.plot(axes=ax, color=colours[i], lw=1.5, alpha=0.5)
            legend_elements.append(Line2D([0], [0], color='w', marker='s',
                                          markerfacecolor=colours[i],
                                          markersize=7,
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
        sat_checked_mask = sat_checked_int == 1
        for sat in sat_information.sat_names[sat_checked_mask]:
            print(sat)
            index = np.where(sat_information.sat_names == sat)
            try:
                fig = Figure()
                ax = fig.add_subplot()
                root = tk.Tk()
                ax.imshow(imgdata, cmap='gray', norm=LogNorm(),
                          interpolation='nearest')
                canvas = FigureCanvasTkAgg(fig, master=root)
                canvas.draw()
                toolbar = NavigationToolbar2Tk(canvas, root)
                toolbar.update()
                canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH,
                                            expand=1)
                canvas.mpl_connect('button_press_event', lambda event: onclick(
                    event, sat_information.sat_locs, index))
                root.update()
                # root.focus_force()
                root.mainloop()
            except ValueError:
                fig = Figure()
                ax = fig.add_subplot()
                root = tk.Tk()
                ax.imshow(abs(imgdata), cmap='gray', norm=LogNorm(
                    vmin=1), interpolation='nearest')
                canvas = FigureCanvasTkAgg(fig, master=root)
                canvas.draw()
                toolbar = NavigationToolbar2Tk(canvas, root)
                toolbar.update()
                canvas.get_tk_widget().pack(side=tk.LEFT,
                                            fill=tk.BOTH,
                                            expand=1)
                canvas.mpl_connect('button_press_event', lambda event: onclick(
                    event, sat_information.sat_locs, index))
                root.update()
                # root.focus_force()
                root.mainloop()
            # mbox('Information',
            #       f'Please select the new position of {sat} on the following image.',
            #       0)
            # cv.namedWindow('TestImage')
            # cv.setMouseCallback('TestImage', change_sat_position, index)
            # logdata = cv.normalize(imgdata, None, alpha=0, beta=10, norm_type=cv.NORM_MINMAX, dtype=cv.CV_32F)
            # cv.imshow('TestImage', logdata)
            # cv.waitKey(0)
        print(sat_information.sat_locs)
        # If some selection is none
        # none_sats = True
        for reversing_index in range(1, num_nan + 1):
            filepath = f"{temp_dir}/{filenames[filenum - reversing_index]}"
            hdr, imgdata = read_fits_file(filepath)
            print(filenames[filenum - reversing_index])
            filename = filenames[filenum - reversing_index]
            exptime = hdr['EXPTIME'] * u.s
            try:
                sat_x, sat_y, bkg_trm, fwhm = TRM_sat_detection(filepath)
            except TypeError:
                print("No satellites detected.")
                continue
            bkg = np.median(bkg_trm)
            bsb = calculate_background_sky_brightness(
                bkg, hdr, exptime, gb_final_transforms)
            fwhm_arcsec = convert_fwhm_to_arcsec_trm(hdr, fwhm)
            airmass = get_image_airmass(hdr)
            photometry_result = perform_photometry_sat(
                sat_x, sat_y, fwhm, imgdata, bkg_trm)
            instr_mags = calculate_magnitudes(photometry_result, exptime)
            instr_mags_sigma = calculate_magnitudes_sigma(
                photometry_result, exptime)
            sat_information = check_if_sat(sat_information,
                                           filenum - reversing_index,
                                           sat_x,
                                           sat_y,
                                           instr_mags,
                                           instr_mags_sigma,
                                           fwhm_arcsec,
                                           airmass,
                                           bsb,
                                           max_distance_from_sat=max_distance_from_sat)
        sat_information.num_nans[sat_checked_mask] = 0
    change_sat_positions_bool = False
    return change_sat_positions_bool, sat_information


def add_new_time_and_filter(hdr, sat_information, filenum):
    t = Time(hdr['DATE-OBS'], format='fits', scale='utc')
    sat_information.sats_table['Time (JD)'][filenum] = t.jd
    try:
        sat_information.sats_table['Filter'][filenum] = hdr['FILTER'][0]
    except KeyError:
        sat_information.sats_table['Filter'][filenum] = 'C'
    sat_information.uncertainty_table['Time (JD)'][filenum] = t.jd
    try:
        sat_information.uncertainty_table['Filter'][filenum] = hdr['FILTER'][0]
    except KeyError:
        sat_information.uncertainty_table['Filter'][filenum] = 'C'
    sat_information.sat_auxiliary_table['Time (JD)'][filenum] = t.jd
    try:
        sat_information.sat_auxiliary_table['Filter'][filenum] =\
            hdr['FILTER'][0]
    except KeyError:
        sat_information.sat_auxiliary_table['Filter'][filenum] = 'C'
    return sat_information


def interpolate_sats(sats_table, uncertainty_table, unique_filters):
    dict_keys = sats_table.columns[2:]
    # unique_filter_table = table.unique(sats_table, keys='Filter')
    # unique_filters = list(unique_filter_table['Filter'])
    unique_filters_sigma = [
        f"{unique_filter}_sigma" for unique_filter in unique_filters]
    filter_table_columns = unique_filters + unique_filters_sigma
    num_columns = len(filter_table_columns)
    time_table = sats_table['Time (JD)']
    times_list = np.array(time_table)
    num_obs = len(sats_table)
    data = np.empty((num_obs, num_columns))
    data.fill(np.nan)
    filter_table = Table(names=filter_table_columns, data=data)
    sat_dict = dict.fromkeys(dict_keys)
    for sat, sat_table in sat_dict.items():
        sat_dict[sat] = hstack([time_table, filter_table])
        for unique_filter in unique_filters:
            mask = sats_table['Filter'] == unique_filter
            filter_all_sats_table = sats_table[mask]
            try:
                filter_interpolated =\
                    np.interp(times_list,
                              filter_all_sats_table['Time (JD)'][
                                  ~np.isnan(filter_all_sats_table[sat])],
                              filter_all_sats_table[sat]
                              [~np.isnan(filter_all_sats_table[sat])])
            except ValueError:
                continue
            filter_interpolated[np.isnan(sats_table[sat])] = np.nan
            sat_dict[sat][unique_filter] = filter_interpolated
            # Uncertainty interpolation
            mask = uncertainty_table['Filter'] == unique_filter
            filter_all_uncertainty_table = uncertainty_table[mask]
            filter_uncertainty_interpolated =\
                np.interp(times_list,
                          filter_all_uncertainty_table['Time (JD)'][
                              ~np.isnan(filter_all_uncertainty_table[sat])],
                          filter_all_uncertainty_table[sat][
                              ~np.isnan(filter_all_uncertainty_table[sat])])
            filter_uncertainty_interpolated[np.isnan(
                uncertainty_table[sat])] = np.nan
            sat_dict[sat][f"{unique_filter}_sigma"] =\
                filter_uncertainty_interpolated
    return sat_dict


def interpolate_aux_data(sat_auxiliary_table, unique_filters, aux_key):
    time_table = sat_auxiliary_table['Time (JD)']
    times_list = np.array(time_table)
    interpolated_table_columns = np.empty(
        (len(sat_auxiliary_table), len(unique_filters)))
    interpolated_table_columns.fill(np.nan)
    interpolated_table_data = Table(
        data=interpolated_table_columns, names=unique_filters)
    interpolated_table = hstack([time_table, interpolated_table_data])
    for unique_filter in unique_filters:
        mask = sat_auxiliary_table['Filter'] == unique_filter
        filter_all_aux_table = sat_auxiliary_table[mask]
        filter_interpolated =\
            np.interp(times_list,
                      filter_all_aux_table['Time (JD)'][~np.isnan(
                          filter_all_aux_table[aux_key])],
                      filter_all_aux_table[aux_key]
                      [~np.isnan(filter_all_aux_table[aux_key])])
        filter_interpolated[np.isnan(sat_auxiliary_table[aux_key])] = np.nan
        interpolated_table[unique_filter] = filter_interpolated
    return interpolated_table


def create_limiting_mag_dict(limiting_mag_table):
    limiting_mag_dict = {"Limiting Mag": limiting_mag_table}
    return limiting_mag_dict


def remove_dim_sats(app_sat_dict, limiting_mag_table, unique_filters):
    for sat, sat_table in app_sat_dict.items():
        for unique_filter in unique_filters:
            mask = sat_table[unique_filter] > limiting_mag_table[unique_filter]
            sat_table[unique_filter][mask] = np.nan
            sat_table[f"{unique_filter}_sigma"][mask] = np.nan


def interpolate_bsb(sat_auxiliary_table, filters_to_plot, times_list):
    bsb_key = "Background Sky Brightness"
    bsb_table_data = np.empty((len(times_list), len(filters_to_plot)))
    bsb_table_data.fill(np.nan)
    bsb_table = Table(names=filters_to_plot, data=bsb_table_data)
    for unique_filter in filters_to_plot:
        mask = sat_auxiliary_table['Filter'] == unique_filter
        filter_all_aux_table = sat_auxiliary_table[mask]
        filter_interpolated =\
            np.interp(times_list,
                      filter_all_aux_table['Time (JD)'][~np.isnan(
                          filter_all_aux_table[bsb_key])],
                      filter_all_aux_table[bsb_key]
                      [~np.isnan(filter_all_aux_table[bsb_key])])
        filter_interpolated[np.isnan(sat_auxiliary_table[bsb_key])] = np.nan
        bsb_table[unique_filter] = filter_interpolated
    return bsb_table


def determine_num_filters(sats_table):
    unique_filter_table = table.unique(sats_table, keys='Filter')
    unique_filters = list(unique_filter_table['Filter'])
    num_filters = len(unique_filters)
    multiple_filters = True
    if num_filters == 1:
        multiple_filters = False
    return unique_filters, num_filters, multiple_filters


def get_all_indicies_combinations(unique_filters,
                                  num_filters,
                                  multiple_filters):
    if multiple_filters:
        all_indices = list(permutations(unique_filters, 2))
        all_indices_formatted = [
            f'{index[0]}-{index[1]}' for index in all_indices]
        return all_indices, all_indices_formatted
    else:
        return


def calculate_timeseries_colour_indices(sat_dict, all_indices):
    index_column_names = [f"{index[0]}-{index[1]}" for index in all_indices]
    index_uncertainty_column_names = [
        f"{index}_sigma" for index in index_column_names]
    index_table_columns = index_column_names + index_uncertainty_column_names
    num_columns = len(index_table_columns)
    colour_indices_dict = dict.fromkeys(sat_dict.keys())
    for sat, sat_table in colour_indices_dict.items():
        time_table = sat_dict[sat]['Time (JD)']
        data = np.empty((len(time_table), num_columns))
        data.fill(np.nan)
        index_table = Table(names=index_table_columns, data=data)
        colour_indices_dict[sat] = hstack([time_table, index_table])
        for index_number, index in enumerate(all_indices):
            colour_indices_dict[sat][index_column_names[index_number]] =\
                sat_dict[sat][index[0]] \
                - sat_dict[sat][index[1]]
            colour_indices_dict[sat][index_uncertainty_column_names[index_number]] =\
                sat_dict[sat][f"{index[0]}_sigma"] \
                + sat_dict[sat][
                f"{index[1]}_sigma"]
    return colour_indices_dict


def choose_indices_to_plot(unique_filters,
                           num_filters,
                           all_indices_formatted,
                           sat_auxiliary_table):
    rootx = tk.Tk()
    rootx.title("Plot Options")
    label = tk.Label(
        rootx, text='Please select the items that you wish to plot.')
    label.grid(row=0, column=0, columnspan=3)
    label = tk.Label(rootx, text='Filters:')
    label.grid(row=1, column=0)
    filter_checked = []
    for filter_num, unique_filter in enumerate(unique_filters):
        filter_checked.append(tk.IntVar())
        checkbutton = tk.Checkbutton(
            rootx, text=unique_filter, variable=filter_checked[filter_num])
        checkbutton.grid(row=filter_num + 2, column=0, sticky=tk.W)
    # none_filter = tk.IntVar()
    # checkbutton = tk.Checkbutton(rootx, text="None", variable=none_filter)
    # checkbutton.grid(row=num_filters+2, column=0, sticky=tk.W)
    label = tk.Label(rootx, text='Colour Indices:')
    label.grid(row=1, column=1)
    index_checked = []
    for index_num, index in enumerate(all_indices_formatted):
        index_checked.append(tk.IntVar())
        checkbutton = tk.Checkbutton(
            rootx, text=index, variable=index_checked[index_num])
        checkbutton.grid(row=index_num + 2, column=1, sticky=tk.W)
    label = tk.Label(rootx, text='Auxiliary Data:')
    label.grid(row=1, column=2)
    aux_checked = []
    for aux_num, aux_type in enumerate(sat_auxiliary_table.columns[2:]):
        aux_checked.append(tk.IntVar())
        checkbutton = tk.Checkbutton(
            rootx, text=aux_type, variable=aux_checked[aux_num])
        checkbutton.grid(row=aux_num + 2, column=2, sticky=tk.W)
    closebutton = tk.Button(rootx, text='OK', command=rootx.destroy)
    closebutton.grid(row=len(all_indices_formatted) + 3, column=1)
    rootx.mainloop()
    filter_checked_int = np.empty(len(filter_checked), dtype=int)
    for filter_num, unique_filter in enumerate(filter_checked):
        filter_checked_int[filter_num] = unique_filter.get()
    filter_checked_mask = filter_checked_int == 1
    unique_filters = np.array(unique_filters)
    filters_to_plot = unique_filters[filter_checked_mask]
    index_checked_int = np.empty(len(index_checked), dtype=int)
    for index_num, index in enumerate(index_checked):
        index_checked_int[index_num] = index.get()
    index_checked_mask = index_checked_int == 1
    all_indices_formatted = np.array(all_indices_formatted)
    indices_to_plot = all_indices_formatted[index_checked_mask]
    aux_checked_int = np.empty(len(aux_checked), dtype=int)
    for aux_num, aux_type in enumerate(aux_checked):
        aux_checked_int[aux_num] = aux_type.get()
    aux_checked_mask = aux_checked_int == 1
    sat_aux_columns_formatted = np.array(list(sat_auxiliary_table.columns[2:]))
    aux_data_to_plot = sat_aux_columns_formatted[aux_checked_mask]
    return filters_to_plot, indices_to_plot, aux_data_to_plot


# def remove_light_curve_outliers(app_sat_dict, unique_filters):
#     unique_filters_sigma = [f"{unique_filter}_sigma" for unique_filter in\
# unique_filters]
#     for sat, sat_table in app_sat_dict.items():
#         sat_table_pandas = sat_table[unique_filters].to_pandas()
#         # sat_table_pandas.set_index('Time (JD)')
#         sat_table_rolling_median = sat_table_pandas.rolling(45).median()
#         sat_table_rolling_std = sat_table_pandas.rolling(45).std()
#         points_to_remove = (sat_table_pandas >= sat_table_rolling_median +\
# (2.5 * sat_table_rolling_std)) | (sat_table_pandas <=\
# sat_table_rolling_median - (2.5 * sat_table_rolling_std))
#         points_to_remove_numpy = points_to_remove.to_numpy()
#         # print(points_to_remove_numpy)
#         sat_table_numpy = sat_table[unique_filters].to_pandas().to_numpy()
#         # print(sat_table_numpy)
#         sat_table_numpy[points_to_remove_numpy] = np.nan
#         # print(sat_table_numpy)
#         uncertainty_table_numpy =\
# sat_table[unique_filters_sigma].to_pandas().to_numpy()
#         uncertainty_table_numpy[points_to_remove_numpy] = np.nan
#         # time_numpy = np.array(sat_table['Time (JD)'])
#         data_table = Table(names=unique_filters, data=sat_table_numpy)
#         uncertainty_table = Table(names=unique_filters_sigma, data=\
#         uncertainty_table_numpy)
#         new_sat_table = hstack([sat_table['Time (JD)'], data_table,
#         uncertainty_table])
#         app_sat_dict[sat] = new_sat_table


# def remove_outliers(app_sat_dict, unique_filters):
#     # unique_filters_sigma = [f"{unique_filter}_sigma" for unique_filter in\
# unique_filters]
#     for sat, sat_table in app_sat_dict.items():
#         for unique_filter in unique_filters:
#             q75, q25 = np.percentile(sat_table[unique_filter], [75, 25])
#             intr_qtr = q75 - q25
#             upper_bound = q75 + (1.5 * intr_qtr)
#             lower_bound = q25 - (1.5 * intr_qtr)
#             # mask = (sat_table[unique_filter] > upper_bound) | \
# (sat_table[unique_filter] < lower_bound)
#             sat_table[unique_filter][sat_table[unique_filter] > upper_bound]\
# = np.nan
#             sat_table[f"{unique_filter}_sigma"]\
# [sat_table[unique_filter] > upper_bound] = np.nan
#             sat_table[unique_filter][sat_table[unique_filter] < lower_bound]\
# = np.nan
#             sat_table[f"{unique_filter}_sigma"]\
# [sat_table[unique_filter] < lower_bound] = np.nan


def plot_light_curve_multiband(sat_table,
                               sat,
                               colour_indices_dict,
                               sat_auxiliary_table,
                               filters_to_plot,
                               indices_to_plot,
                               aux_data_to_plot):
    nrows = 2 + len(aux_data_to_plot)
    height_ratios = [1] * nrows
    height_ratios[0] = 3
    gs_kw = dict(height_ratios=height_ratios)
    times_list = np.array(sat_table['Time (JD)'])
    times_obj = Time(times_list, format='jd', scale='utc')
    times_datetime = times_obj.to_value('datetime')
    fig, axs = plt.subplots(nrows=nrows, sharex=True,
                            gridspec_kw=gs_kw, figsize=(7.5, 9.5))
    for filter_num, unique_filter in enumerate(filters_to_plot):
        _, _, bars = axs[0].errorbar(times_datetime, sat_table[unique_filter],
                                     yerr=sat_table[f"{unique_filter}_sigma"],
                                     fmt='o',
                                     markersize=2,
                                     capsize=0,
                                     elinewidth=0.75,
                                     label=unique_filter)
        [bar.set_alpha(0.3) for bar in bars]
    for index_num, colour_index in enumerate(indices_to_plot):
        _, _, bars = axs[1].errorbar(times_datetime,
                                     colour_indices_dict[sat][colour_index],
                                     yerr=colour_indices_dict[sat]
                                     [f"{colour_index}_sigma"],
                                     fmt='o',
                                     markersize=2,
                                     capsize=0,
                                     elinewidth=0.75,
                                     label=colour_index)
        [bar.set_alpha(0.3) for bar in bars]
    for aux_num, aux_data in enumerate(aux_data_to_plot, start=2):
        axs[aux_num].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        axs[aux_num].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        axs[aux_num].set_ylabel(aux_data)
        for filter_num, unique_filter in enumerate(filters_to_plot):
            mask = sat_auxiliary_table['Filter'] == unique_filter
            filter_all_aux_table = sat_auxiliary_table[mask]
            aux_times_list = filter_all_aux_table['Time (JD)']
            aux_times_obj = Time(aux_times_list, format='jd', scale='utc')
            aux_times_datetime = aux_times_obj.to_value('datetime')
            axs[aux_num].plot(
                aux_times_datetime, filter_all_aux_table[aux_data], 'o', ms=2, label=unique_filter)
        if len(filters_to_plot) > 1:
            axs[aux_num].legend()
        if aux_data == "BSB":
            axs[aux_num].invert_yaxis()
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    axs[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axs[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axs[0].legend()
    axs[1].legend()
    axs[0].invert_yaxis()
    axs[1].invert_yaxis()
    axs[nrows - 1].set_xlabel("Time (UTC)")
    axs[0].set_ylabel("Magnitude")
    axs[1].set_ylabel("Colour Index")
    plt.tight_layout()
    return fig, axs


def choose_aux_data_to_plot(sat_auxiliary_table):
    rootx = tk.Tk()
    rootx.title("Plot Options")
    label = tk.Label(
        rootx, text='Please select the items that you wish to plot.')
    label.grid(row=0, column=0)
    label = tk.Label(rootx, text='Auxiliary Data:')
    label.grid(row=1, column=0)
    aux_checked = []
    for aux_num, aux_type in enumerate(sat_auxiliary_table.columns[2:]):
        aux_checked.append(tk.IntVar())
        checkbutton = tk.Checkbutton(
            rootx, text=aux_type, variable=aux_checked[aux_num])
        checkbutton.grid(row=aux_num + 2, column=0, sticky=tk.W)
    closebutton = tk.Button(rootx, text='OK', command=rootx.destroy)
    closebutton.grid(row=len(sat_auxiliary_table.columns[2:]) + 2, column=0)
    rootx.mainloop()
    aux_checked_int = np.empty(len(aux_checked), dtype=int)
    for aux_num, aux_type in enumerate(aux_checked):
        aux_checked_int[aux_num] = aux_type.get()
    aux_checked_mask = aux_checked_int == 1
    sat_aux_columns_formatted = np.array(list(sat_auxiliary_table.columns[2:]))
    aux_data_to_plot = sat_aux_columns_formatted[aux_checked_mask]
    return aux_data_to_plot


def plot_light_curve_singleband(sats_table,
                                sat,
                                uncertainty_table,
                                sat_auxiliary_table,
                                aux_data_to_plot):
    nrows = 1 + len(aux_data_to_plot)
    height_ratios = [1] * nrows
    height_ratios[0] = 3
    times_list = np.array(sats_table['Time (JD)'])
    times_obj = Time(times_list, format='jd', scale='utc')
    times_datetime = times_obj.to_value('datetime')
    gs_kw = dict(height_ratios=height_ratios)
    fig, axs = plt.subplots(nrows=nrows, sharex=True,
                            gridspec_kw=gs_kw, figsize=(7.5, 9.5))
    # fig = plt.figure(figsize=(7.5, 9.5))
    # spec = gridspec.GridSpec(nrows=nrows, ncols=1, height_ratios=height_ratios)
    if nrows > 1:
        _, _, bars = axs[0].errorbar(times_datetime,
                                     sats_table[sat],
                                     yerr=uncertainty_table[sat],
                                     fmt='o',
                                     markersize=2,
                                     capsize=0,
                                     elinewidth=0.75)
        [bar.set_alpha(0.3) for bar in bars]
        for aux_num, aux_data in enumerate(aux_data_to_plot, start=1):
            axs[aux_num].plot(
                times_datetime, sat_auxiliary_table[aux_data],
                'o',
                ms=2,
                label=aux_data)
            axs[aux_num].xaxis.set_major_formatter(
                mdates.DateFormatter('%H:%M'))
            axs[aux_num].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            axs[aux_num].set_ylabel(aux_data)
            if aux_data == "BSB":
                axs[aux_num].invert_yaxis()
        axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        axs[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        axs[0].invert_yaxis()
        axs[nrows - 1].set_xlabel("Time (UTC)")
        axs[0].set_ylabel("Instrumental Magnitude")
        plt.tight_layout()
    else:
        _, _, bars = axs.errorbar(times_datetime, sats_table[sat],
                                  yerr=uncertainty_table[sat],
                                  fmt='o',
                                  markersize=2,
                                  capsize=0,
                                  elinewidth=0.75)
        [bar.set_alpha(0.3) for bar in bars]
        axs.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        axs.invert_yaxis()
        axs.set_xlabel("Time (UTC)")
        axs.set_ylabel("Instrumental Magnitude")
        plt.tight_layout()
    return fig, axs


def axis_limits_singleband_gui(sats_table,
                               uncertainty_table,
                               sat_auxiliary_table,
                               aux_data_to_plot,
                               save_loc):
    def refresh_plots(canvas):
        fig_list[0], axs = plot_light_curve_singleband(sats_table,
                                                       sat,
                                                       uncertainty_table,
                                                       sat_auxiliary_table,
                                                       aux_data_to_plot)
        for entry_num, entry in enumerate(entries):
            if entry_num % 2 == 0:
                axs[entry_num // 2].set_ylim(bottom=float(entry.get()))
            elif entry_num % 2 == 1:
                axs[entry_num // 2].set_ylim(top=float(entry.get()))
        canvas.get_tk_widget().delete()
        canvas = FigureCanvasTkAgg(fig_list[0], master=root)
        canvas.draw()
        canvas.get_tk_widget().grid(column=0, row=0, rowspan=nrows,
                                    sticky=tk.NSEW)
        plt.close()

    def reset_plots(canvas):
        fig_list[0], axs = plot_light_curve_singleband(sats_table,
                                                       sat,
                                                       uncertainty_table,
                                                       sat_auxiliary_table,
                                                       aux_data_to_plot)
        row_num = 1
        for ax in axs:
            entries[row_num - 1].delete(0, tk.END)
            entries[row_num - 1].insert(tk.END, f"{ax.get_ylim()[0]:0.3f}")
            entries[row_num].delete(0, tk.END)
            entries[row_num].insert(tk.END, f"{ax.get_ylim()[1]:0.3f}")
            row_num += 2
        canvas.get_tk_widget().delete()
        canvas = FigureCanvasTkAgg(fig_list[0], master=root)
        canvas.draw()
        canvas.get_tk_widget().grid(column=0, row=0, rowspan=nrows, sticky=tk.NSEW)
        plt.close()

    def save_plots(fig):
        fig.set_size_inches(7.5, 9.5)
        fig.savefig(f"{save_loc}/{sat} Light Curve.pdf",
                    format='pdf', bbox_inches='tight')
        root.destroy()

    try:
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman"]})
    except Exception:
        pass
    nrows = 10 * (2 + len(aux_data_to_plot))
    for sat in sats_table.columns[2:]:
        fig_list = [None]
        fig_list[0], axs = plot_light_curve_singleband(sats_table,
                                                       sat,
                                                       uncertainty_table,
                                                       sat_auxiliary_table,
                                                       aux_data_to_plot)

        root = tk.Tk()
        for x in range(nrows + 1):
            tk.Grid.rowconfigure(root, x, weight=1)
        root.title(sat)
        canvas = FigureCanvasTkAgg(fig_list[0], master=root)
        canvas.draw_idle()
        toolbarFrame = tk.Frame(master=root)
        tk.Grid.columnconfigure(root, 0, weight=1)
        toolbarFrame.grid(row=nrows + 1, column=0, sticky=tk.NSEW)
        toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
        # toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas.get_tk_widget().grid(column=0, row=0, rowspan=nrows,
                                    sticky=tk.NSEW)
        root.update()
        tk.Grid.columnconfigure(root, 1, weight=1)
        label = tk.Label(text="Please set the y axis limits.")
        label.grid(column=1, row=0, columnspan=2)
        row_num = 1
        entries = [tk.Entry(root, width=10) for _ in range(2 * len(axs))]
        for ax in axs:
            label = tk.Label(root, text=ax.get_ylabel())
            label.grid(column=1, row=row_num, columnspan=2, sticky=tk.N)
            # entries[row_num-1] = tk.Entry(root, width=10)
            entries[row_num - 1].insert(tk.END, f"{ax.get_ylim()[0]:0.3f}")
            entries[row_num - 1].grid(column=1, row=row_num + 1,
                                      padx=(10, 5), pady=0, sticky=tk.N)
            # entries[row_num] = tk.Entry(root, width=10)
            entries[row_num].insert(tk.END, f"{ax.get_ylim()[1]:0.3f}")
            entries[row_num].grid(column=2, row=row_num + 1,
                                  padx=(5, 10), pady=0, sticky=tk.N)
            # print(ax.get_ylim())
            row_num += 2
        button = tk.Button(root, text="Refresh",
                           command=lambda: refresh_plots(canvas))
        button.grid(column=1, row=row_num)
        button = tk.Button(root, text="Reset",
                           command=lambda: reset_plots(canvas))
        button.grid(column=2, row=row_num)
        button = tk.Button(root, text="Save and Close",
                           command=lambda: save_plots(fig_list[0]))
        button.grid(column=1, row=row_num + 1, columnspan=2)
        root.mainloop()
        # fig.show()
        plt.close()
    return


def axis_limits_multiband_gui(app_sat_dict,
                              colour_indices_dict,
                              sat_auxiliary_table,
                              filters_to_plot,
                              indices_to_plot,
                              aux_data_to_plot,
                              save_loc):
    def refresh_plots(canvas):
        # fig.clear()
        fig_list[0], axs = plot_light_curve_multiband(sat_table,
                                                      sat,
                                                      colour_indices_dict,
                                                      sat_auxiliary_table,
                                                      filters_to_plot,
                                                      indices_to_plot,
                                                      aux_data_to_plot)
        for entry_num, entry in enumerate(entries):
            if entry_num % 2 == 0:
                axs[entry_num // 2].set_ylim(bottom=float(entry.get()))
            elif entry_num % 2 == 1:
                axs[entry_num // 2].set_ylim(top=float(entry.get()))
        canvas.get_tk_widget().delete()
        canvas = FigureCanvasTkAgg(fig_list[0], master=root)
        canvas.draw()
        canvas.get_tk_widget().grid(column=0, row=0, rowspan=nrows,
                                    sticky=tk.NSEW)
        plt.close()

    def reset_plots(canvas):
        fig_list[0], axs = plot_light_curve_multiband(sat_table,
                                                      sat,
                                                      colour_indices_dict,
                                                      sat_auxiliary_table,
                                                      filters_to_plot,
                                                      indices_to_plot,
                                                      aux_data_to_plot)
        row_num = 1
        for ax in axs:
            entries[row_num - 1].delete(0, tk.END)
            entries[row_num - 1].insert(tk.END, f"{ax.get_ylim()[0]:0.3f}")
            entries[row_num].delete(0, tk.END)
            entries[row_num].insert(tk.END, f"{ax.get_ylim()[1]:0.3f}")
            row_num += 2
        canvas.get_tk_widget().delete()
        canvas = FigureCanvasTkAgg(fig_list[0], master=root)
        canvas.draw()
        canvas.get_tk_widget().grid(column=0, row=0, rowspan=nrows,
                                    sticky=tk.NSEW)
        plt.close()

    def save_plots(fig):
        fig.set_size_inches(7.5, 9.5)
        fig.savefig(f"{save_loc}/{sat} Light Curve.pdf",
                    format='pdf', bbox_inches='tight')
        root.destroy()

    try:
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman"]})
    except Exception:
        pass
    nrows = 10 * (2 + len(aux_data_to_plot))
    for sat, sat_table in app_sat_dict.items():
        fig_list = [None]
        fig_list[0], axs = plot_light_curve_multiband(sat_table,
                                                      sat,
                                                      colour_indices_dict,
                                                      sat_auxiliary_table,
                                                      filters_to_plot,
                                                      indices_to_plot,
                                                      aux_data_to_plot)

        root = tk.Tk()
        for x in range(nrows + 1):
            tk.Grid.rowconfigure(root, x, weight=1)
        root.title(sat)
        canvas = FigureCanvasTkAgg(fig_list[0], master=root)
        canvas.draw_idle()
        toolbarFrame = tk.Frame(master=root)
        tk.Grid.columnconfigure(root, 0, weight=1)
        toolbarFrame.grid(row=nrows + 1, column=0, sticky=tk.NSEW)
        toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
        # toolbar = NavigationToolbar2Tk(canvas, root)
        toolbar.update()
        canvas.get_tk_widget().grid(column=0, row=0, rowspan=nrows,
                                    sticky=tk.NSEW)
        root.update()
        tk.Grid.columnconfigure(root, 1, weight=1)
        label = tk.Label(text="Please set the y axis limits.")
        label.grid(column=1, row=0, columnspan=2)
        row_num = 1
        entries = [tk.Entry(root, width=10) for _ in range(2 * len(axs))]
        for ax in axs:
            label = tk.Label(root, text=ax.get_ylabel())
            label.grid(column=1, row=row_num, columnspan=2, sticky=tk.N)
            # entries[row_num-1] = tk.Entry(root, width=10)
            entries[row_num - 1].insert(tk.END, f"{ax.get_ylim()[0]:0.3f}")
            entries[row_num - 1].grid(column=1, row=row_num + 1,
                                      padx=(10, 5), pady=0, sticky=tk.N)
            # entries[row_num] = tk.Entry(root, width=10)
            entries[row_num].insert(tk.END, f"{ax.get_ylim()[1]:0.3f}")
            entries[row_num].grid(column=2, row=row_num + 1,
                                  padx=(5, 10), pady=0, sticky=tk.N)
            # print(ax.get_ylim())
            row_num += 2
        button = tk.Button(root, text="Refresh",
                           command=lambda: refresh_plots(canvas))
        button.grid(column=1, row=row_num)
        button = tk.Button(root, text="Reset",
                           command=lambda: reset_plots(canvas))
        button.grid(column=2, row=row_num)
        button = tk.Button(root, text="Save and Close",
                           command=lambda: save_plots(fig_list[0]))
        button.grid(column=1, row=row_num + 1, columnspan=2)
        root.mainloop()
        # fig.show()
        plt.close()
    return


def save_interpolated_light_curve(sat_dict, save_loc, suffix=None):
    if not suffix:
        for sat, sat_table in sat_dict.items():
            ascii.write(
                sat_table, output=f"{save_loc}/{sat}_{suffix}.csv",
                format='csv')
    else:
        for sat, sat_table in sat_dict.items():
            ascii.write(
                sat_table, output=f"{save_loc}/{sat}.csv",
                format='csv')


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
                            **kwargs):
    # TODO: Docstring.
    reference_stars, ref_star_positions = read_ref_stars(ref_stars_file)
    gb_transform_table_columns = init_gb_transform_table_columns()
    auxiliary_data_columns = init_auxiliary_data_columns()

    if save_plots:
        save_loc = kwargs.get('save_loc')
        if not os.path.exists(save_loc):
            os.mkdir(save_loc)

    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith((file_suffix)):
                filepath = os.path.join(dirpath, filename)
                print(filepath)
                hdr, imgdata = read_fits_file(filepath)
                exptime = hdr[exposure_key]
                bkg, bkg_std = calculate_img_bkg(imgdata)
                irafsources = detecting_stars(
                    imgdata, bkg=bkg, bkg_std=bkg_std)
                if not irafsources:
                    continue
                _, fwhm, fwhm_std = calculate_fwhm(irafsources)
                photometry_result = perform_photometry(
                    irafsources, fwhm, imgdata, bkg=bkg)
                fluxes = np.array(photometry_result['flux_fit'])
                instr_mags = calculate_magnitudes(photometry_result, exptime)
                instr_mags_sigma = calculate_magnitudes_sigma(
                    photometry_result, exptime)
                wcs = WCS(hdr)
                skypositions = convert_pixel_to_ra_dec(irafsources, wcs)
                altazpositions = None
                try:
                    altazpositions = convert_ra_dec_to_alt_az(skypositions,
                                                              hdr,
                                                              lat_key=lat_key,
                                                              lon_key=lon_key,
                                                              elev_key=elev_key)
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
                print("match")
                for colour_index in colour_indices:

                    field = get_field_name(matched_stars, name_key=name_key)

                    try:
                        len(matched_stars.img_instr_mag)
                    except TypeError:
                        print("Only 1 reference star detected in the image.")
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

                    auxiliary_data_columns =\
                        update_auxiliary_data_columns(auxiliary_data_columns,
                                                      filename,
                                                      exptime,
                                                      fwhm,
                                                      fwhm_std,
                                                      matched_stars)
                    try:
                        if not save_plots:
                            c_fci,\
                                c_fci_sigma,\
                                zprime_f,\
                                zprime_f_sigma\
                                = ground_based_first_order_transforms(
                                    matched_stars,
                                    instr_filter,
                                    colour_index,
                                    plot_results=plot_results)
                        else:
                            unique_id = filename
                            c_fci,\
                                c_fci_sigma,\
                                zprime_f,\
                                zprime_f_sigma\
                                = ground_based_first_order_transforms(
                                    matched_stars,
                                    instr_filter,
                                    colour_index,
                                    plot_results=plot_results,
                                    save_plots=save_plots,
                                    save_loc=save_loc,
                                    unique_id=unique_id)
                    except Exception:
                        continue
                    gb_transform_table_columns =\
                        update_gb_transform_table_columns(
                            gb_transform_table_columns,
                            field,
                            c_fci,
                            c_fci_sigma,
                            zprime_f,
                            zprime_f_sigma,
                            instr_filter,
                            colour_index,
                            altazpositions)
    gb_transform_table = create_gb_transform_table(gb_transform_table_columns)
    auxiliary_data_table = create_auxiliary_data_table(auxiliary_data_columns)
    if remove_large_airmass_bool:
        gb_transform_table = remove_large_airmass(gb_transform_table, 3.0)
    if save_plots:
        gb_final_transforms =\
            ground_based_second_order_transforms(
                gb_transform_table,
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
        write_table_to_latex(gb_final_transforms,
                             f"{os.path.join(save_loc,'gb_final_transforms')}.txt",
                             formats=formats)
        formats = {
            'exptime': '%0.3f',
            'fwhm': '%0.3f',
            'fwhm_sigma': '%0.3f',
            'avg mag_sigma': '%0.3f',
            'std mag_sigma': '%0.3f'
        }
        ascii.write(auxiliary_data_table,
                    f"{os.path.join(save_loc, 'auxiliary_data_table')}.csv",
                    format='csv',
                    formats=formats)
    else:
        gb_final_transforms = ground_based_second_order_transforms(
            gb_transform_table,
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
    #             large_table_columns = update_large_table_columns(
    #             large_table_columns,
    #             matched_stars,
    #             hdr,
    #             exptime,
    #             ground_based=True,
    #             name_key='Name')

    # large_stars_table = create_large_stars_table(large_table_columns,
    #                                               ground_based=True)
    # large_stars_table.pprint_all()
    # # large_stars_table = remove_large_airmass(large_stars_table)
    # stars_table = group_each_star(large_stars_table, ground_based=True)
    # stars_table.pprint_all()

    # instr_filters = ['b', 'v', 'r', 'i']
    # app_mag_table = Table(stars_table['Field', 'Name', 'V_ref', '(B-V)',
    # '(U-B)', '(V-R)', '(V-I)', 'V_sigma'])
    # for instr_filter in instr_filters:
    #     app_mag_table_filter =\
    #    apply_gb_transforms_VERIFICATION(gb_final_transforms,
    #                                     stars_table, instr_filter)
    #     app_mag_table = hstack([app_mag_table,
    # app_mag_table_filter[instr_filter.upper()]])

    # app_mag_table.pprint_all()

    # import matplotlib.pyplot as plt
    # # import matplotlib
    # # matplotlib.use('TkAgg')
    # plt.plot(app_mag_table['V_ref'] +\
    #    app_mag_table['(B-V)'], app_mag_table['B'], 'o')
    # m, b = np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['B'])] +
    # app_mag_table['(B-V)'][~np.isnan(app_mag_table['B'])],
    # app_mag_table['B'][~np.isnan(app_mag_table['B'])], 1)
    # plt.plot(app_mag_table['V_ref'] + app_mag_table['(B-V)'],
    # m*(app_mag_table['V_ref'] + app_mag_table['(B-V)'])+b, '-',
    # label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'] +
    # app_mag_table['(B-V)'], app_mag_table['V_ref'] +
    # app_mag_table['(B-V)'], '-', label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('B (calculated)')
    # plt.xlabel('B (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()

    # plt.plot(app_mag_table['V_ref'], app_mag_table['V'], 'o')
    # m, b = np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['V'])],
    # app_mag_table['V'][~np.isnan(app_mag_table['V'])], 1)
    # plt.plot(app_mag_table['V_ref'], m*(app_mag_table['V_ref'])+b, '-',
    # label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'], app_mag_table['V_ref'], '-', label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('V (calculated)')
    # plt.xlabel('V (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()

    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-R)'],
    # app_mag_table['R'], 'o')
    # m, b = np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['R'])] -
    # app_mag_table['(V-R)'][~np.isnan(app_mag_table['R'])],
    #                   app_mag_table['R'][~np.isnan(app_mag_table['R'])], 1)
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-R)'],
    #           m*(app_mag_table['V_ref'] - app_mag_table['(V-R)'])+b,
    #           '-', label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-R)'],
    # app_mag_table['V_ref'] - app_mag_table['(V-R)'], '-', label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('R (calculated)')
    # plt.xlabel('R (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()

    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-I)'],
    # app_mag_table['I'], 'o')
    # m, b = np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['I'])]
    # - app_mag_table['(V-I)'][~np.isnan(app_mag_table['I'])],
    #                   app_mag_table['I'][~np.isnan(app_mag_table['I'])], 1)
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-I)'],
    #           m*(app_mag_table['V_ref'] - app_mag_table['(V-I)'])+b,
    #           '-', label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-I)'],
    # app_mag_table['V_ref'] - app_mag_table['(V-I)'], '-', label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('I (calculated)')
    # plt.xlabel('I (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()
    return gb_final_transforms, auxiliary_data_table


def verify_gb_transforms(directory,
                         ref_stars_file,
                         gb_final_transforms,
                         plot_results=False,
                         save_plots=False,
                         file_suffix=(".fits", ".fit", ".fts"),
                         exposure_key='EXPTIME',
                         name_key='Name',
                         **kwargs):
    # TODO: Docstring.
    reference_stars, ref_star_positions = read_ref_stars(ref_stars_file)
    large_table_columns = init_large_table_columns()

    if save_plots:
        save_loc = kwargs.get('save_loc')
        unique_id = kwargs.get('unique_id')
        if not os.path.exists(save_loc):
            os.mkdir(save_loc)

    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(file_suffix):
                filepath = os.path.join(dirpath, filename)
                hdr, imgdata = read_fits_file(filepath)
                exptime = hdr[exposure_key]
                bkg, bkg_std = calculate_img_bkg(imgdata)
                irafsources = detecting_stars(
                    imgdata, bkg=bkg, bkg_std=bkg_std)
                if not irafsources:
                    continue
                _, fwhm, fwhm_std = calculate_fwhm(irafsources)
                photometry_result = perform_photometry(
                    irafsources, fwhm, imgdata, bkg=bkg)
                fluxes = np.array(photometry_result['flux_fit'])
                instr_mags = calculate_magnitudes(photometry_result, exptime)
                instr_mags_sigma = calculate_magnitudes_sigma(
                    photometry_result, exptime)
                wcs = WCS(hdr)
                skypositions = convert_pixel_to_ra_dec(irafsources, wcs)
                try:
                    altazpositions =\
                        convert_ra_dec_to_alt_az(skypositions,
                                                 hdr,
                                                 lat_key='SITELAT',
                                                 lon_key='SITELONG',
                                                 elev_key='SITEELEV')
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

                large_table_columns =\
                    update_large_table_columns(large_table_columns,
                                               filepath,
                                               matched_stars,
                                               hdr,
                                               exptime,
                                               ground_based=True,
                                               name_key=name_key)
    large_stars_table = create_large_stars_table(
        large_table_columns, ground_based=True)
    # large_stars_table = remove_large_airmass(large_stars_table, max_airmass=2.0)
    stars_table, different_filter_list = group_each_star_GB(large_stars_table)
    stars_table.pprint(max_lines=30, max_width=200)

    # instr_filters = ['b', 'v', 'r', 'i']
    app_mag_table = Table(stars_table[
        'Field', 'Name', 'V_ref', 'B-V', 'U-B', 'V-R', 'V-I', 'V_sigma',
        'e_B-V', 'e_U-B', 'e_V-R', 'e_V-I'])
    for instr_filter in different_filter_list:
        app_mag_filter = instr_filter.upper()
        app_filter_sigma = f"{app_mag_filter}_sigma"
        app_mag_table_filter = apply_gb_transforms_VERIFICATION(
            gb_final_transforms, stars_table, instr_filter)
        app_mag_table = hstack(
            [app_mag_table,
             app_mag_table_filter[app_mag_filter, app_filter_sigma]])

    app_mag_table.pprint(max_lines=30, max_width=200)

    import matplotlib.pyplot as plt
    # import matplotlib
    # matplotlib.use('TkAgg')
    for instr_filter in different_filter_list:
        app_filter = instr_filter.upper()
        try:
            ref_app_mag, ref_app_mag_sigma, _, _ = get_app_mag_and_index_AVG(
                app_mag_table, instr_filter)
        except KeyError:
            ref_app_mag, ref_app_mag_sigma, _, _ = get_app_mag_and_index(
                app_mag_table, instr_filter)
        app_filter_sigma = f"{app_filter}_sigma"
        x = ref_app_mag
        x_err = ref_app_mag_sigma
        y = app_mag_table[app_filter]
        # y_err = app_mag_table[app_filter_sigma]
        max_y_err = max(app_mag_table[app_filter_sigma])
        print(max_y_err)
        y_err = np.nan_to_num(app_mag_table[app_filter_sigma], nan=max_y_err)
        y_err = np.array(y_err)
        y_err[y_err == 0] = max_y_err
        fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
        try:
            fitted_line, mask = or_fit(line_init, x, y)
            filtered_data = np.ma.masked_array(y, mask=mask)
            m = fitted_line.slope.value
            b = fitted_line.intercept.value
            # fitted_line, mask = or_fit(line_init, x, y, weights=1.0/y_err)
            # filtered_data = np.ma.masked_array(y, mask=mask)
            # m = fitted_line.slope.value
            # b = fitted_line.intercept.value
            # if m == 1 and b == 0:
            #     fitted_line, mask = or_fit(line_init, x, y)
            #     filtered_data = np.ma.masked_array(y, mask=mask)
            #     m = fitted_line.slope.value
            #     b = fitted_line.intercept.value
        except TypeError:
            continue
        plt.errorbar(x, y, yerr=y_err, fmt='o', fillstyle='none',
                     capsize=0, label="Clipped Data")
        plt.plot(x, filtered_data, 'o', color='#1f77b4', label="Fitted Data")
        plt.plot(x, m * x + b, '-', label=f'y={m:.3f}x+{b:.3f}')
        plt.plot(x, x, '-', label='y=x')
        plt.title('Calculated Magnitude vs. Reference Magnitude')
        plt.ylabel(f"{app_filter} (Calculated)")
        plt.xlabel(f"{app_filter} (Reference)")
        plt.legend()
        if save_plots:
            save_loc = f"{os.path.join(kwargs.get('save_loc'), f'{app_filter}_CalcVsRefMag')}.png"
            plt.savefig(save_loc)
        if plot_results:
            plt.show(block=True)
        plt.close()

    # plt.plot(app_mag_table['V_ref'] + app_mag_table['(B-V)'],
    # app_mag_table['B'], 'o')
    # m, b = np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['B'])]\
    # + app_mag_table['(B-V)'][~np.isnan(app_mag_table['B'])],
    # app_mag_table['B'][~np.isnan(app_mag_table['B'])], 1)
    # plt.plot(app_mag_table['V_ref'] + app_mag_table['(B-V)'],
    # m*(app_mag_table['V_ref'] + app_mag_table['(B-V)'])+b, '-',
    # label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'] + app_mag_table['(B-V)'],
    # app_mag_table['V_ref'] + app_mag_table['(B-V)'], '-', label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('B (calculated)')
    # plt.xlabel('B (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()

    # plt.plot(app_mag_table['V_ref'], app_mag_table['V'], 'o')
    # m, b = np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['V'])],
    # app_mag_table['V'][~np.isnan(app_mag_table['V'])], 1)
    # plt.plot(app_mag_table['V_ref'], m*(app_mag_table['V_ref'])+b, '-',
    # label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'], app_mag_table['V_ref'], '-',
    # label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('V (calculated)')
    # plt.xlabel('V (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()

    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-R)'],
    # app_mag_table['R'], 'o')
    # m, b =\
    # np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['R'])] -\
    # app_mag_table['(V-R)'][~np.isnan(app_mag_table['R'])],
    #                   app_mag_table['R'][~np.isnan(app_mag_table['R'])], 1)
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-R)'],
    #           m*(app_mag_table['V_ref'] - app_mag_table['(V-R)'])+b,
    #           '-', label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-R)'],
    # app_mag_table['V_ref'] - app_mag_table['(V-R)'], '-', label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('R (calculated)')
    # plt.xlabel('R (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()

    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-I)'],
    # app_mag_table['I'], 'o')
    # m, b = np.polyfit(app_mag_table['V_ref'][~np.isnan(app_mag_table['I'])]\
    # - app_mag_table['(V-I)'][~np.isnan(app_mag_table['I'])],
    #                   app_mag_table['I'][~np.isnan(app_mag_table['I'])], 1)
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-I)'],
    #           m*(app_mag_table['V_ref'] - app_mag_table['(V-I)'])+b,
    #           '-', label=f'y={m:.3f}x+{b:.3f}')
    # plt.plot(app_mag_table['V_ref'] - app_mag_table['(V-I)'],
    # app_mag_table['V_ref'] - app_mag_table['(V-I)'], '-', label='y=x')
    # plt.title('Calculated Magnitude vs. Reference Magnitude')
    # plt.ylabel('I (calculated)')
    # plt.xlabel('I (Reference)')
    # plt.legend()
    # plt.show(block=True)
    # plt.close()
    return app_mag_table


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
    reference_stars, ref_star_positions = read_ref_stars(ref_stars_file)
    large_table_columns = init_large_table_columns()

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
    with open(os.path.join(save_loc, 'CalculationVerificationSplit.txt'), 'a') as f:
        f.write('File Path'+'\t'+'Calculation/Verification')
        for calc_file in calculation_files:
            f.write('\n'+f'{calc_file}'+'\t'+'Calculation')
        for verify_file in verification_files:
            f.write('\n'+f'{verify_file}'+'\t'+'Verification')
    "Iterate over the images."
    for file_num, filepath in enumerate(tqdm(calculation_files)):
        # filepath = os.path.join(dirpath, filename)
        hdr, imgdata = read_fits_file(filepath)
        exptime = hdr[exposure_key]
        bkg, bkg_std = calculate_img_bkg(imgdata)
        irafsources = detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
        if not irafsources:
            with open(os.path.join(save_loc, 'ExcludedFiles.txt'), 'a') as f:
                f.write(f'{filepath}')
                f.write('\t')
                f.write('No sources detected')
                f.write('\n')
            excluded_files += 1
            continue
        _, fwhm, fwhm_std = calculate_fwhm(irafsources)
        photometry_result = perform_photometry(
            irafsources, fwhm, imgdata, bkg=bkg)
        fluxes = np.array(photometry_result['flux_fit'])
        instr_mags = calculate_magnitudes(photometry_result, exptime)
        instr_mags_sigma = calculate_magnitudes_sigma(
            photometry_result, exptime)
        wcs = WCS(hdr)
        skypositions = convert_pixel_to_ra_dec(irafsources, wcs)
        try:
            altazpositions = convert_ra_dec_to_alt_az(skypositions, hdr,
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
        matched_stars = find_ref_stars(reference_stars,
                                       ref_star_positions,
                                       skypositions,
                                       instr_mags,
                                       instr_mags_sigma,
                                       fluxes,
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

        large_table_columns = update_large_table_columns(large_table_columns,
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
    large_stars_table = create_large_stars_table(
        large_table_columns, ground_based=True)
    # large_stars_table = remove_large_airmass(large_stars_table, max_airmass=3.0)
    stars_table, different_filter_list = group_each_star_GB(large_stars_table)
    stars_table.pprint(max_lines=30, max_width=200)
    if save_plots:
        ascii.write(stars_table, os.path.join(
            save_loc, 'stars_table.csv'), format='csv')
    gb_transform_table =\
        calc_gb_first_transforms_AVG(stars_table,
                                     different_filter_list,
                                     save_loc,
                                     plot_results=plot_results,
                                     save_plots=save_plots)
    gb_transform_table.pprint(max_lines=30, max_width=-1)
    gb_final_transforms =\
        ground_based_second_order_transforms(gb_transform_table,
                                             plot_results=plot_results,
                                             save_plots=save_plots,
                                             save_loc=save_loc)
    gb_final_transforms.pprint_all()
    if save_plots:
        gb_final_transforms =\
            ground_based_second_order_transforms(gb_transform_table,
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
        write_table_to_latex(gb_final_transforms,
                             f"{os.path.join(save_loc, 'gb_final_transforms')}.txt",
                             formats=formats)
    hidden_transform_table = None
    exoatmospheric_table = exoatmospheric_mags_Warner(
        stars_table, gb_final_transforms, different_filter_list)
    hidden_transform_table =\
        hidden_transform_Warner(exoatmospheric_table,
                                gb_final_transforms,
                                different_filter_list,
                                save_plots,
                                save_loc=save_loc)
    verify_save_loc = os.path.join(save_loc, 'Verification')
    app_mag_table = verify_gb_transforms_auto(directory,
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


def init_star_aux_table_columns():
    filenames = []
    times = []
    filters = []
    fwhms_pixels = []
    fwhm_pixels_sigma = []
    fwhms_arcsec = []
    fwhm_arcsec_sigma = []
    num_sources = []
    background_sky_brightness = []
    background_sky_brightness_sigma = []
    azimuth = []
    elevation = []
    airmass = []
    star_aux_table_columns = namedtuple('sat_aux_table_columns',
                                        ['filenames',
                                         'times',
                                         'filters',
                                         'fwhms_pixels',
                                         'fwhm_pixels_sigma',
                                         'fwhms_arcsec',
                                         'fwhm_arcsec_sigma',
                                         'num_sources',
                                         'background_sky_brightness',
                                         'background_sky_brightness_sigma',
                                         'azimuth',
                                         'elevation',
                                         'airmass'])
    return star_aux_table_columns(filenames,
                                  times,
                                  filters,
                                  fwhms_pixels,
                                  fwhm_pixels_sigma,
                                  fwhms_arcsec,
                                  fwhm_arcsec_sigma,
                                  num_sources,
                                  background_sky_brightness,
                                  background_sky_brightness_sigma,
                                  azimuth,
                                  elevation,
                                  airmass)


def update_star_aux_columns(star_aux_table_columns,
                            filename,
                            time,
                            img_filter,
                            fwhm,
                            fwhm_sigma,
                            fwhm_arcsec,
                            fwhm_arcsec_sigma,
                            num_sources,
                            background_sky_brightness,
                            background_sky_brightness_sigma,
                            azimuth,
                            elevation,
                            airmass):
    updated_star_aux_table_columns = star_aux_table_columns
    updated_star_aux_table_columns.filenames.append(filename)
    updated_star_aux_table_columns.times.append(time)
    updated_star_aux_table_columns.filters.append(img_filter)
    updated_star_aux_table_columns.fwhms_pixels.append(fwhm)
    updated_star_aux_table_columns.fwhm_pixels_sigma.append(fwhm_sigma)
    updated_star_aux_table_columns.fwhms_arcsec.append(fwhm_arcsec)
    updated_star_aux_table_columns.fwhm_arcsec_sigma.append(fwhm_arcsec_sigma)
    updated_star_aux_table_columns.num_sources.append(num_sources)
    updated_star_aux_table_columns.background_sky_brightness.append(
        background_sky_brightness)
    updated_star_aux_table_columns.background_sky_brightness_sigma.append(
        background_sky_brightness_sigma)
    updated_star_aux_table_columns.azimuth.append(azimuth)
    updated_star_aux_table_columns.elevation.append(elevation)
    updated_star_aux_table_columns.airmass.append(airmass)
    return updated_star_aux_table_columns


def create_star_aux_table(star_aux_table_columns):
    star_aux_table = Table(
        names=[
            'filename',
            'Time (JD)',
            'filter',
            'FWHM_pixel',
            'FWHM_pixel_sigma',
            'FWHM_arcsec',
            'FWHM_arcsec_sigma',
            'Number of Sources',
            'BSB',
            'BSB_sigma',
            'Azimuth',
            'Elevation',
            'X'
        ],
        data=[
            star_aux_table_columns.filenames,
            star_aux_table_columns.times,
            star_aux_table_columns.filters,
            star_aux_table_columns.fwhms_pixels,
            star_aux_table_columns.fwhm_pixels_sigma,
            star_aux_table_columns.fwhms_arcsec,
            star_aux_table_columns.fwhm_arcsec_sigma,
            star_aux_table_columns.num_sources,
            star_aux_table_columns.background_sky_brightness,
            star_aux_table_columns.background_sky_brightness_sigma,
            star_aux_table_columns.azimuth,
            star_aux_table_columns.elevation,
            star_aux_table_columns.airmass
        ]
    )
    return star_aux_table


def calculate_slopes_Warner(stars_table, different_filter_list, save_plots, **kwargs):
    stars_for_second_order_extinction, multiple_stars = get_stars_with_multiple_observations(
        stars_table)
    slope_filters = [
        f"slope_{different_filter}" for different_filter in different_filter_list]
    intercept_filters = [
        f"intercept_{different_filter}" for different_filter in different_filter_list]
    slope_filters_sigma = [
        f"slope_{different_filter}_sigma" for different_filter in different_filter_list]
    intercept_filters_sigma = [
        f"intercept_{different_filter}_sigma" for different_filter in different_filter_list]
    nan_array = np.empty(len(multiple_stars))
    nan_array.fill(np.nan)
    data_filter_table = [
        nan_array for different_filter in different_filter_list]
    star_index_columns = [
        'Field',
        'Name',
        'V_ref',
        'B-V',
        'U-B',
        'V-R',
        'V-I',
        'V_sigma',
        'e_B-V',
        'e_U-B',
        'e_V-R',
        'e_V-I'
    ]
    star_index_table = Table(
        names=star_index_columns,
        data=[
            np.empty(len(multiple_stars), dtype=object),
            np.empty(len(multiple_stars), dtype=object),
            nan_array,
            nan_array,
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
    slope_table = Table(names=slope_filters, data=data_filter_table)
    intercept_table = Table(names=intercept_filters, data=data_filter_table)
    slope_sigma_table = Table(
        names=slope_filters_sigma, data=data_filter_table)
    intercept_sigma_table = Table(
        names=intercept_filters_sigma, data=data_filter_table)
    slopes_table = hstack([star_index_table, slope_table,
                          intercept_table, slope_sigma_table,
                          intercept_sigma_table])
    colors = cm.rainbow(np.linspace(0, 1, len(multiple_stars)))
    for unique_filter in different_filter_list:
        # current_filter = stars_table[unique_filter]
        x_current_filter = f"X_{unique_filter}"
        # unique_stars = table.unique(current_filter, keys='Name')
        # colors = cm.rainbow(np.linspace(0, 1, len(multiple_stars)))
        for i, unique_star in enumerate(multiple_stars):
            # print(current_filter)
            star_mask = stars_table['Name'] == unique_star
            current_star = stars_table[star_mask]
            # print(current_star[star_index_columns][0])
            # print(current_star[star_index_columns])
            for field in star_index_columns:
                slopes_table[field][i] = current_star[field][0]
            # slopes_table[star_index_columns][i] = list(current_star[star_index_columns][0])
            # print(unique_filter)
            # print(current_star)
            x_nan_indices = np.isnan(current_star[x_current_filter])
            y_nan_indices = np.isnan(current_star[unique_filter])
            nan_indices = (x_nan_indices | y_nan_indices)
            try:
                X_plot = np.arange(
                    start=min(current_star[x_current_filter][~np.isnan(
                        current_star[x_current_filter])]) - 0.02,
                    stop=max(current_star[x_current_filter][~np.isnan(
                        current_star[x_current_filter])]) + 0.02,
                    step=0.01)
            except ValueError:
                slopes_table[f"slope_{unique_filter}"][i] = np.nan
                slopes_table[f"intercept_{unique_filter}"][i] = np.nan
                slopes_table[f"slope_{unique_filter}_sigma"][i] = np.nan
                slopes_table[f"intercept_{unique_filter}_sigma"][i] = np.nan
                continue
            # m, b = np.polyfit(current_star[x_current_filter], current_star[unique_filter], 1)
            fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
            try:
                fitted_line, mask = or_fit(line_init, current_star[x_current_filter][~nan_indices],
                                           current_star[unique_filter][~nan_indices])
            except TypeError:
                # print(current_star[x_current_filter])
                # print(current_star[unique_filter])
                slopes_table[f"slope_{unique_filter}"][i] = np.nan
                slopes_table[f"intercept_{unique_filter}"][i] = np.nan
                slopes_table[f"slope_{unique_filter}_sigma"][i] = np.nan
                slopes_table[f"intercept_{unique_filter}_sigma"][i] = np.nan
                continue
            filtered_data = np.ma.masked_array(
                current_star[unique_filter][~nan_indices], mask=mask)
            m = fitted_line.slope.value
            b = fitted_line.intercept.value
            cov = fit.fit_info['param_cov']
            try:
                m_sigma = sqrt(cov[0][0])
                b_sigma = sqrt(cov[1][1])
            except TypeError:
                m_sigma = np.nan
                b_sigma = np.nan
            plt.plot(current_star[x_current_filter],
                     current_star[unique_filter], 'o', fillstyle='none',
                     color=colors[i], label="Clipped Data")
            plt.plot(current_star[x_current_filter][~nan_indices], filtered_data, 'o', color=colors[i],
                     label="Fitted Data")
            # plt.scatter(current_star[x_current_filter], current_star[unique_filter], color=colors[i], label=unique_star)
            plt.plot(X_plot, m * X_plot + b, color=colors[i])
            slopes_table[f"slope_{unique_filter}"][i] = m
            slopes_table[f"intercept_{unique_filter}"][i] = b
            slopes_table[f"slope_{unique_filter}_sigma"][i] = m_sigma
            slopes_table[f"intercept_{unique_filter}_sigma"][i] = b_sigma
        # plt.plot(current_filter['X'], current_filter['mag_instrumental'], 'o')
        plt.xlabel('X')
        plt.ylabel(unique_filter.lower())
        plt.title("Slope of a Star's magnitude vs. airmass")
        plt.ylim([min(stars_table[unique_filter][~np.isnan(stars_table[unique_filter])]) * 1.05,
                  max(stars_table[unique_filter][~np.isnan(stars_table[unique_filter])]) * 0.95])
        plt.gca().invert_yaxis()
        # plt.title(unique_field)
        # plt.legend()
        if save_plots:
            save_loc = f"{os.path.join(kwargs.get('save_loc'), f'Slopes{unique_filter}')}.png"
            plt.savefig(save_loc)
        plt.show()
        plt.close()
    # Column names to keep.
    # 'Field','Name','V_ref','B-V','U-B','V-R','V-I','V_sigma','e_B-V','e_U-B','e_V-R','e_V-I'
    # print(stars_for_second_order_extinction.columns)
    # slopes_table.pprint(max_lines=-1, max_width=250)
    return slopes_table


def get_stars_with_multiple_observations(stars_table):
    list_of_stars = stars_table['Name']
    count_of_stars = Counter(list_of_stars)
    multiple_stars = [
        star for star in count_of_stars if count_of_stars[star] > 1]
    mask = np.in1d(np.array(stars_table['Name']), np.array(multiple_stars))
    stars_for_second_order_extinction = stars_table[list(mask)]
    # stars_for_second_order_extinction.pprint(max_lines=-1, max_width=-1)
    return stars_for_second_order_extinction, multiple_stars


def second_order_extinction_calc_Warner(slopes_table,
                                        different_filter_list,
                                        save_plots, **kwargs):
    filter_column = []
    CI_column = []
    k_primeprime_column = []
    k_prime_column = []
    k_primeprime_sigma_column = []
    k_prime_sigma_column = []
    for different_filter in different_filter_list:
        colour_indices = get_all_colour_indices(different_filter)
        for colour_index in colour_indices:
            # print(different_filter)
            # print(colour_index)
            # print(slopes_table.pprint_all())
            filter_column.append(different_filter)
            CI_column.append(colour_index)
            # print(slopes_table[colour_index])
            ci_plot = np.arange(min(slopes_table[colour_index][~np.isnan(slopes_table[colour_index])]) - 0.1,
                                max(slopes_table[colour_index][~np.isnan(
                                    slopes_table[colour_index])]) + 0.1,
                                step=0.01)
            # fit, or_fit, line_init = init_linear_fitting(sigma=1.5)
            fit, or_fit, line_init = init_linear_fitting(
                niter=100, sigma=2.0, slope=0.0, intercept=0.5)
            # try:
            x_nan_indices = np.isnan(slopes_table[colour_index])
            y_nan_indices = np.isnan(slopes_table[f"slope_{different_filter}"])
            nan_indices = (x_nan_indices | y_nan_indices)
            fitted_line, mask =\
                or_fit(line_init,
                       slopes_table[colour_index][~nan_indices],
                       slopes_table[f"slope_{different_filter}"][~nan_indices])
            # except TypeError as e:
            #     # print(current_star[x_current_filter])
            #     # print(current_star[unique_filter])
            #     # slopes_table[f"slope_{unique_filter}"][i] = np.nan
            #     # slopes_table[f"intercept_{unique_filter}"][i] = np.nan
            #     # slopes_table[f"slope_{unique_filter}_sigma"][i] = np.nan
            #     # slopes_table[f"intercept_{unique_filter}_sigma"][i] = np.nan
            #     k_primeprime_column.append(np.nan)
            #     k_prime_column.append(np.nan)
            #     k_primeprime_sigma_column.append(np.nan)
            #     k_prime_sigma_column.append(np.nan)
            #     print(e)
            #     continue
            filtered_data = np.ma.masked_array(
                slopes_table[f"slope_{different_filter}"][~nan_indices], mask=mask)
            m = fitted_line.slope.value
            b = fitted_line.intercept.value
            cov = fit.fit_info['param_cov']
            try:
                m_sigma = sqrt(cov[0][0])
                b_sigma = sqrt(cov[1][1])
            except TypeError:
                m_sigma = np.nan
                b_sigma = np.nan
            k_primeprime_column.append(m)
            k_prime_column.append(b)
            k_primeprime_sigma_column.append(m_sigma)
            k_prime_sigma_column.append(b_sigma)
            # plt.plot(slopes_table[colour_index], slopes_table[f"slope_{different_filter}"], 'o')
            plt.plot(slopes_table[colour_index], slopes_table[f"slope_{different_filter}"], 'o',
                     fillstyle='none', label="Clipped Data")
            plt.plot(slopes_table[colour_index][~nan_indices],
                     filtered_data, 'o', color='#1f77b4', label="Fitted Data")
            plt.plot(ci_plot, m * ci_plot + b, '-',
                     label=f"k''={m:0.3f}, k'={b:0.3f}")
            plt.ylabel(f'slope$_{{{different_filter}}}$')
            plt.xlabel(colour_index)
            plt.title(
                f"Second Order extinction (slope$_{{{different_filter}}}$ v. {colour_index})")
            plt.legend()
            if save_plots:
                save_loc = f"{os.path.join(kwargs.get('save_loc'), f'SecondOrderExtinction{different_filter}{colour_index}')}.png"
                plt.savefig(save_loc)
            plt.show()
            plt.close()
    # print(filter_column)
    # print(CI_column)
    # print(k_primeprime_column)
    # print(k_prime_column)
    # print(k_primeprime_sigma_column)
    # print(k_prime_sigma_column)
    extinction_table_Warner = Table(
        names=[
            'filter',
            'CI',
            'k\'\'_fCI',
            'k\'\'_fCI_sigma',
            'k\'_f',
            'k\'_f_sigma'
        ],
        data=[
            filter_column,
            CI_column,
            k_primeprime_column,
            k_primeprime_sigma_column,
            k_prime_column,
            k_prime_sigma_column
        ])
    return extinction_table_Warner


def calculate_slopes_Buchheim(stars_table,
                              different_filter_list,
                              save_plots,
                              **kwargs):
    stars_for_second_order_extinction,\
        multiple_stars = get_stars_with_multiple_observations(
            stars_table)
    slope_filters = [
        f"slope_{different_filter}" for different_filter in different_filter_list]
    intercept_filters = [
        f"intercept_{different_filter}" for different_filter in different_filter_list]
    slope_filters_sigma = [
        f"slope_{different_filter}_sigma" for different_filter in different_filter_list]
    intercept_filters_sigma = [
        f"intercept_{different_filter}_sigma" for different_filter in different_filter_list]
    nan_array = np.empty(len(multiple_stars))
    nan_array.fill(np.nan)
    data_filter_table = [
        nan_array for different_filter in different_filter_list]
    num_filters = len(different_filter_list)
    all_indices, all_indices_formatted = get_all_indicies_combinations(
        different_filter_list, num_filters, multiple_filters=True)
    data_instr_index_table = [
        nan_array for different_index in all_indices_formatted]

    star_index_columns = [
        'Field',
        'Name',
        'V_ref',
        'B-V',
        'U-B',
        'V-R',
        'V-I',
        'V_sigma',
        'e_B-V',
        'e_U-B',
        'e_V-R',
        'e_V-I'
    ]
    star_index_table = Table(
        names=star_index_columns,
        data=[
            np.empty(len(multiple_stars), dtype=object),
            np.empty(len(multiple_stars), dtype=object),
            nan_array,
            nan_array,
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
    slope_table = Table(names=slope_filters, data=data_filter_table)
    intercept_table = Table(names=intercept_filters, data=data_filter_table)
    slope_sigma_table = Table(
        names=slope_filters_sigma, data=data_filter_table)
    intercept_sigma_table = Table(
        names=intercept_filters_sigma, data=data_filter_table)
    instr_index_table = Table(
        names=all_indices_formatted, data=data_instr_index_table)
    for i, unique_star in enumerate(multiple_stars):
        star_mask = stars_table['Name'] == unique_star
        current_star = stars_table[star_mask]
        for instr_index_name in all_indices_formatted:
            # first_mag = current_star[instr_index_name[0]]
            # second_mag = current_star[instr_index_name[-1]]
            instr_index = np.mean(current_star[instr_index_name])
            instr_index_table[instr_index_name][i] = instr_index

    instr_index_table.pprint(max_lines=-1, max_width=150)
    slopes_table = hstack([star_index_table,
                           slope_table,
                           intercept_table,
                           slope_sigma_table,
                           intercept_sigma_table,
                           instr_index_table])
    colors = cm.rainbow(np.linspace(0, 1, len(multiple_stars)))
    for unique_filter in different_filter_list:
        # current_filter = stars_table[unique_filter]
        x_current_filter = f"X_{unique_filter}"
        # unique_stars = table.unique(current_filter, keys='Name')
        # colors = cm.rainbow(np.linspace(0, 1, len(multiple_stars)))
        for i, unique_star in enumerate(multiple_stars):
            # print(current_filter)
            star_mask = stars_table['Name'] == unique_star
            current_star = stars_table[star_mask]
            # print(current_star[star_index_columns][0])
            # print(current_star[star_index_columns])
            for field in star_index_columns:
                slopes_table[field][i] = current_star[field][0]
            # slopes_table[star_index_columns][i] = list(current_star[star_index_columns][0])
            # print(unique_filter)
            # print(current_star)
            x_nan_indices = np.isnan(current_star[x_current_filter])
            y_nan_indices = np.isnan(current_star[unique_filter])
            nan_indices = (x_nan_indices | y_nan_indices)
            try:
                X_plot = np.arange(
                    start=min(current_star[x_current_filter][~np.isnan(
                        current_star[x_current_filter])]) - 0.02,
                    stop=max(current_star[x_current_filter][~np.isnan(
                        current_star[x_current_filter])]) + 0.02,
                    step=0.01)
            except ValueError:
                slopes_table[f"slope_{unique_filter}"][i] = np.nan
                slopes_table[f"intercept_{unique_filter}"][i] = np.nan
                slopes_table[f"slope_{unique_filter}_sigma"][i] = np.nan
                slopes_table[f"intercept_{unique_filter}_sigma"][i] = np.nan
                continue
            # m, b = np.polyfit(current_star[x_current_filter], current_star[unique_filter], 1)
            fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
            try:
                fitted_line, mask = or_fit(line_init, current_star[x_current_filter][~nan_indices],
                                           current_star[unique_filter][~nan_indices])
            except TypeError:
                # print(current_star[x_current_filter])
                # print(current_star[unique_filter])
                slopes_table[f"slope_{unique_filter}"][i] = np.nan
                slopes_table[f"intercept_{unique_filter}"][i] = np.nan
                slopes_table[f"slope_{unique_filter}_sigma"][i] = np.nan
                slopes_table[f"intercept_{unique_filter}_sigma"][i] = np.nan
                continue
            filtered_data = np.ma.masked_array(
                current_star[unique_filter][~nan_indices], mask=mask)
            m = fitted_line.slope.value
            b = fitted_line.intercept.value
            cov = fit.fit_info['param_cov']
            try:
                m_sigma = sqrt(cov[0][0])
                b_sigma = sqrt(cov[1][1])
            except TypeError:
                m_sigma = np.nan
                b_sigma = np.nan
            plt.plot(current_star[x_current_filter],
                     current_star[unique_filter],
                     'o', fillstyle='none',
                     color=colors[i], label="Clipped Data")
            plt.plot(current_star[x_current_filter][~nan_indices],
                     filtered_data, 'o', color=colors[i],
                     label="Fitted Data")
            # plt.scatter(current_star[x_current_filter],
            # current_star[unique_filter], color=colors[i], label=unique_star)
            plt.plot(X_plot, m * X_plot + b, color=colors[i])
            slopes_table[f"slope_{unique_filter}"][i] = m
            slopes_table[f"intercept_{unique_filter}"][i] = b
            slopes_table[f"slope_{unique_filter}_sigma"][i] = m_sigma
            slopes_table[f"intercept_{unique_filter}_sigma"][i] = b_sigma
        # plt.plot(current_filter['X'], current_filter['mag_instrumental'], 'o')
        plt.xlabel('X')
        plt.ylabel(unique_filter.lower())
        plt.title("Slope of a Star's magnitude vs. airmass")
        plt.ylim([min(stars_table[unique_filter][~np.isnan(stars_table[unique_filter])]) * 1.05,
                  max(stars_table[unique_filter][~np.isnan(stars_table[unique_filter])]) * 0.95])
        plt.gca().invert_yaxis()
        # plt.title(unique_field)
        # plt.legend()
        if save_plots:
            save_loc = f"{os.path.join(kwargs.get('save_loc'), f'Slopes{unique_filter}')}.png"
            plt.savefig(save_loc)
        plt.show()
        plt.close()
    # Column names to keep.
    # 'Field','Name','V_ref','B-V','U-B','V-R','V-I','V_sigma','e_B-V','e_U-B','e_V-R','e_V-I'
    # print(stars_for_second_order_extinction.columns)
    # slopes_table.pprint(max_lines=-1, max_width=250)
    return slopes_table


def second_order_extinction_calc_Buchheim(stars_table,
                                          different_filter_list,
                                          save_plots, **kwargs):
    x_list = [
        f'X_{different_filter}' for different_filter in different_filter_list]
    stars_table.sort(x_list)
    unique_fields = table.unique(stars_table, keys=['Field'])
    unique_stars = table.unique(stars_table, keys=['Name'])
    mags_list = np.empty(2)
    ci_list = np.empty(2)
    x_list = np.empty(2)
    ci = 'b-v'
    # print(unique_fields)
    for different_filter in different_filter_list:
        for field in unique_fields['Field']:
            current_field_index = stars_table['Field'] == field
            current_field = stars_table[current_field_index]
            names_list = list(current_field['Name'])
            star_combinations = list(combinations(names_list, 2))
            if len(star_combinations) > 1:
                delta_mags_list = np.empty(len(star_combinations))
                delta_ci_X_list = np.empty(len(star_combinations))
                delta_mags_list.fill(np.nan)
                delta_ci_X_list.fill(np.nan)
                for k, combination in enumerate(star_combinations):
                    mags_list.fill(np.nan)
                    ci_list.fill(np.nan)
                    x_list.fill(np.nan)
                    for i, star in enumerate(combination):
                        current_star_index = stars_table['Name'] == star
                        if sum(current_star_index) > 1:
                            rand_index = choice(
                                np.where(current_star_index)[0])
                            current_star_index[np.where(current_star_index)[0][~np.where(
                                np.where(current_star_index)[0] == rand_index)[0]]] = False
                        current_star = stars_table[current_star_index]
                        mags_list[i] = current_star[different_filter]
                        x_list[i] = current_star[f"X_{different_filter}"]
                        try:
                            ci_list[i] = current_star[ci]
                            table_ci = ci
                        except KeyError:
                            if 'v' in ci:
                                table_ci = ci.replace('v', 'g')
                            else:
                                table_ci = ci
                            try:
                                ci_list[i] = current_star[table_ci]
                            except KeyError:
                                if 'b' in ci:
                                    table_ci = table_ci.replace('b', 'u')
                                ci_list[i] = current_star[table_ci]
                    delta_mag = np.abs(mags_list[1] - mags_list[0])
                    delta_ci = np.abs(ci_list[1] - ci_list[0])
                    x_diff = np.abs(x_list[1] - x_list[0])
                    if delta_mag == 0 and delta_ci == 0 and x_diff == 0:
                        continue
                    if x_diff > 0.5:
                        continue
                    avg_x = np.mean(x_list)
                    delta_mags_list[k] = delta_mag
                    ########## FIXME: Fix this. #############
                    delta_ci_X_list[k] = delta_ci * avg_x
                fit, or_fit, line_init = init_linear_fitting(
                    niter=100, sigma=2.5, slope=0.0)
                # try:
                x_nan_indices = np.isnan(delta_ci_X_list)
                y_nan_indices = np.isnan(delta_mags_list)
                nan_indices = (x_nan_indices | y_nan_indices)
                try:
                    fitted_line, mask = or_fit(line_init, delta_ci_X_list[~nan_indices],
                                               delta_mags_list[~nan_indices])
                except TypeError:
                    continue
                filtered_data = np.ma.masked_array(
                    delta_mags_list[~nan_indices], mask=mask)
                m = fitted_line.slope.value
                b = fitted_line.intercept.value
                cov = fit.fit_info['param_cov']
                try:
                    m_sigma = sqrt(cov[0][0])
                    b_sigma = sqrt(cov[1][1])
                except TypeError:
                    m_sigma = np.nan
                    b_sigma = np.nan
                delta_ci_X_plot = np.arange(min(delta_ci_X_list[~np.isnan(delta_ci_X_list)]) - 0.1,
                                            max(delta_ci_X_list[~np.isnan(delta_ci_X_list)]) + 0.1, step=0.01)
                plt.plot(delta_ci_X_list, delta_mags_list, 'o',
                         fillstyle='none', label="Clipped Data")
                plt.plot(delta_ci_X_list[~nan_indices], filtered_data,
                         'o', color='#1f77b4', label="Fitted Data")
                plt.plot(delta_ci_X_list, m * delta_ci_X_list + b, '-',
                         label=f"k''={m:0.3f}, $\Delta{{{different_filter}}}_0$={b:0.3f}")
                plt.ylabel(f'$\Delta{{{different_filter}}}$')
                plt.xlabel(f"$\Delta({{{table_ci}}})$ $\cdot X$")
                plt.title(
                    f"Second Order extinction ($\Delta{{{different_filter}}}$ v. $\Delta$ $({{{table_ci}}})$ $\cdot X$)")
                plt.legend()
                if save_plots:
                    save_loc = f"{os.path.join(kwargs.get('save_loc'), f'SecondOrderExtinction-delta_{different_filter}_{table_ci}_{field}')}.png"
                    plt.savefig(save_loc)
                plt.show()
                plt.close()
                # plt.plot(delta_ci_X_list, delta_mags_list, 'o')
                # plt.show()
                # plt.close()
        # print(star_combinations)
    return


def extinction_calc_Buchheim_sect6(slopes_table, different_filter_list, save_plots, **kwargs):
    filter_column = []
    CI_column = []
    k_primeprime_column = []
    k_prime_column = []
    k_primeprime_sigma_column = []
    k_prime_sigma_column = []
    for different_filter in different_filter_list:
        colour_indices = get_all_colour_indices(different_filter)
        colour_index, _ = get_colour_index_lower(different_filter)
        # colour_index = 'B-V'
        # for colour_index in colour_indices:
        filter_column.append(different_filter)
        # ci = re.sub('[^a-zA-Z]+', '', colour_index)
        ci = colour_index.lower()
        # ci = ci.lower()
        try:
            ci_plot = np.arange(min(slopes_table[ci][~np.isnan(slopes_table[ci])]) - 0.1,
                                max(slopes_table[ci][~np.isnan(slopes_table[ci])]) + 0.1, step=0.01)
            table_ci = ci
            print(table_ci)
        except KeyError:
            if 'v' in ci:
                table_ci = ci.replace('v', 'g')
            else:
                table_ci = ci
            try:
                ci_plot = np.arange(min(slopes_table[table_ci][~np.isnan(slopes_table[table_ci])]) - 0.1,
                                    max(slopes_table[table_ci][~np.isnan(slopes_table[table_ci])]) + 0.1, step=0.01)
            except KeyError:
                if 'b' in ci:
                    table_ci = table_ci.replace('b', 'u')
                ci_plot = np.arange(min(slopes_table[table_ci][~np.isnan(slopes_table[table_ci])]) - 0.1,
                                    max(slopes_table[table_ci][~np.isnan(slopes_table[table_ci])]) + 0.1, step=0.01)
        CI_column.append(table_ci)
        # ci_plot = np.arange(min(slopes_table[ci]) - 0.1, max(slopes_table[ci]) + 0.1, step=0.01)
        fit, or_fit, line_init = init_linear_fitting(
            niter=100, sigma=2.0, slope=0.0, intercept=0.5)
        # try:
        x_nan_indices = np.isnan(slopes_table[table_ci])
        y_nan_indices = np.isnan(slopes_table[f"slope_{different_filter}"])
        nan_indices = (x_nan_indices | y_nan_indices)
        fitted_line, mask =\
            or_fit(line_init,
                   slopes_table[table_ci][~nan_indices],
                   slopes_table[f"slope_{different_filter}"][~nan_indices])
        # except TypeError as e:
        #     # print(current_star[x_current_filter])
        #     # print(current_star[unique_filter])
        #     # slopes_table[f"slope_{unique_filter}"][i] = np.nan
        #     # slopes_table[f"intercept_{unique_filter}"][i] = np.nan
        #     # slopes_table[f"slope_{unique_filter}_sigma"][i] = np.nan
        #     # slopes_table[f"intercept_{unique_filter}_sigma"][i] = np.nan
        #     k_primeprime_column.append(np.nan)
        #     k_prime_column.append(np.nan)
        #     k_primeprime_sigma_column.append(np.nan)
        #     k_prime_sigma_column.append(np.nan)
        #     print(e)
        #     continue
        filtered_data = np.ma.masked_array(
            slopes_table[f"slope_{different_filter}"][~nan_indices], mask=mask)
        m = fitted_line.slope.value
        b = fitted_line.intercept.value
        cov = fit.fit_info['param_cov']
        try:
            m_sigma = sqrt(cov[0][0])
            b_sigma = sqrt(cov[1][1])
        except TypeError:
            m_sigma = np.nan
            b_sigma = np.nan
        k_primeprime_column.append(m)
        k_prime_column.append(b)
        k_primeprime_sigma_column.append(m_sigma)
        k_prime_sigma_column.append(b_sigma)
        # plt.plot(slopes_table[colour_index], slopes_table[f"slope_{different_filter}"], 'o')
        plt.plot(slopes_table[table_ci],
                 slopes_table[f"slope_{different_filter}"],
                 'o',
                 fillstyle='none', label="Clipped Data")
        plt.plot(slopes_table[table_ci][~nan_indices],
                 filtered_data,
                 'o', color='#1f77b4', label="Fitted Data")
        plt.plot(ci_plot, m * ci_plot + b, '-',
                 label=f"k''={m:0.3f}, k'={b:0.3f}")
        plt.ylabel(f'slope$_{{{different_filter}}}$')
        plt.xlabel(table_ci)
        plt.title(
            f"Second Order extinction (slope$_{{{different_filter}}}$ v. {table_ci})")
        plt.legend()
        if save_plots:
            save_loc = f"{os.path.join(kwargs.get('save_loc'), f'SecondOrderExtinction-{different_filter}_{table_ci}')}.png"
            plt.savefig(save_loc)
        plt.show()
        plt.close()
    # print(filter_column)
    # print(CI_column)
    # print(k_primeprime_column)
    # print(k_prime_column)
    # print(k_primeprime_sigma_column)
    # print(k_prime_sigma_column)
    extinction_table_Buchheim = Table(
        names=[
            'filter',
            'CI',
            'k\'\'_fCI',
            'k\'\'_fCI_sigma',
            'k\'_f',
            'k\'_f_sigma'
        ],
        data=[
            filter_column,
            CI_column,
            k_primeprime_column,
            k_primeprime_sigma_column,
            k_prime_column,
            k_prime_sigma_column
        ])
    return extinction_table_Buchheim


def exoatmospheric_mags_Buchheim(stars_table,
                                 extinction_table_Warner,
                                 different_filter_list):
    # m_0 = m - k'_f * X - k''_fCI * X * CI
    nan_array = np.empty(len(stars_table))
    nan_array.fill(np.nan)
    star_index_columns = [
        'Field',
        'Name',
        'V_ref',
        'B-V',
        'U-B',
        'V-R',
        'V-I',
        'V_sigma',
        'e_B-V',
        'e_U-B',
        'e_V-R',
        'e_V-I'
    ]
    # try:
    exoatmospheric_table_begin = Table(
        names=[
            'Field',
            'Name',
            'V_ref',
            'B-V',
            'U-B',
            'V-R',
            'V-I',
            'V_sigma',
            'e_B-V',
            'e_U-B',
            'e_V-R',
            'e_V-I',
            # 'b B-V',
            # 'g B-V',
            # 'g V-R',
            # 'g V-I',
            # 'r V-R'
        ],
        data=[
            np.empty(len(stars_table), dtype=object),
            np.empty(len(stars_table), dtype=object),
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            # nan_array,
            # nan_array,
            # nan_array,
            # nan_array,
            # nan_array
        ])
    exoatmospheric_table_filter_ci_columns = []
    for instr_filter in different_filter_list:
        colour_indices_list = get_all_colour_indices(instr_filter)
        for colour_index in colour_indices_list:
            exoatmospheric_table_filter_ci_columns.append(
                f"{instr_filter} {colour_index}")
    exoatmospheric_table_filter_ci_data = np.empty(
        (len(nan_array), len(exoatmospheric_table_filter_ci_columns)))
    exoatmospheric_table_filter_ci_data.fill(np.nan)
    exoatmospheric_table_filter_ci =\
        Table(names=exoatmospheric_table_filter_ci_columns,
              data=exoatmospheric_table_filter_ci_data)
    exoatmospheric_table = hstack(
        (exoatmospheric_table_begin, exoatmospheric_table_filter_ci))
    i = 0
    for star in stars_table:
        for field in star_index_columns:
            exoatmospheric_table[field][i] = star[field]
        # exoatmospheric_table['Name'][i] = star['Name']
        for different_filter in different_filter_list:
            colour_indices = get_all_colour_indices(different_filter)
            x_column_name = f"X_{different_filter}"
            for colour_index in colour_indices:
                mask = ((extinction_table_Warner['filter'] == different_filter) & (
                    extinction_table_Warner['CI'] == colour_index))
                row_of_extinctions = extinction_table_Warner[mask]
                # if len(row_of_transforms) == 0 and instr_filter == 'v':
                #     instr_filter = 'g'
                #     mask = ((gb_final_transforms['filter'] == instr_filter) & (gb_final_transforms['CI'] == colour_index))
                #     row_of_transforms = gb_final_transforms[mask]
                instr_mag = star[different_filter]
                X = star[x_column_name]
                CI = star[colour_index]
                k_primeprime = row_of_extinctions['k\'\'_fCI']
                k_prime = row_of_extinctions['k\'_f']
                exoatmospheric_mag = float(
                    instr_mag - (k_prime * X) - (k_primeprime * X * CI))
                exoatmospheric_table[f"{different_filter} {colour_index}"][i]\
                    = exoatmospheric_mag
                # if different_filter == 'g':
                #     print(f"Instrumental mag for {different_filter} and\
                #    {colour_index}:")
                #     print(f"{instr_mag:0.3f}")
                #     print(f"Exoatmospheric mag for {different_filter} and\
                #    {colour_index}:")
                #     print(f"{exoatmospheric_mag:0.3f}")
        i += 1
    # except KeyError:
    #     exoatmospheric_table = Table(
    #         names=[
    #             'Field',
    #             'Name',
    #             'V_ref',
    #             'B-V',
    #             'U-B',
    #             'V-R',
    #             'V-I',
    #             'V_sigma',
    #             'e_B-V',
    #             'e_U-B',
    #             'e_V-R',
    #             'e_V-I',
    #             'u B-V',
    #             'g B-V',
    #             'g V-R',
    #             'g V-I',
    #             'r V-R',
    #             'i V-I'
    #         ],
    #         data=[
    #             np.empty(len(stars_table), dtype=object),
    #             np.empty(len(stars_table), dtype=object),
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array
    #         ])
    #     i = 0
    #     for star in stars_table:
    #         for field in star_index_columns:
    #             exoatmospheric_table[field][i] = star[field]
    #         # exoatmospheric_table['Name'][i] = star['Name']
    #         for different_filter in different_filter_list:
    #             colour_indices = get_all_colour_indices(different_filter)
    #             x_column_name = f"X_{different_filter}"
    #             for colour_index in colour_indices:
    #                 mask =\
    #   ((extinction_table_Warner['filter'] == different_filter) &\
    #   (extinction_table_Warner['CI'] == colour_index))
    #                 row_of_extinctions = extinction_table_Warner[mask]
    #                 # if len(row_of_transforms) == 0 and instr_filter == 'v':
    #                 #     instr_filter = 'g'
    #                 #     mask =\
    #    ((gb_final_transforms['filter'] == instr_filter) &\
    #   (gb_final_transforms['CI'] == colour_index))
    #                 #     row_of_transforms = gb_final_transforms[mask]
    #                 instr_mag = star[different_filter]
    #                 X = star[x_column_name]
    #                 CI = star[colour_index]
    #                 k_primeprime = row_of_extinctions['k\'\'_fCI']
    #                 k_prime = row_of_extinctions['k\'_f']
    #                 exoatmospheric_mag = float(instr_mag - (k_prime * X) -\
    #   (k_primeprime * X * CI))
    #                 exoatmospheric_table[f"{different_filter} \
    # {colour_index}"][i] = exoatmospheric_mag
    #                 # if different_filter == 'g':
    #                 #     print(f"Instrumental mag for {different_filter} and {colour_index}:")
    #                 #     print(f"{instr_mag:0.3f}")
    #                 #     print(f"Exoatmospheric mag for {different_filter} and {colour_index}:")
    #                 #     print(f"{exoatmospheric_mag:0.3f}")
    #         i += 1
    # exoatmospheric_table.pprint(max_lines=-1, max_width=250)
    return exoatmospheric_table


def exoatmospheric_mags_Warner(stars_table,
                               extinction_table_Warner,
                               different_filter_list):
    # m_0 = m - k'_f * X - k''_fCI * X * CI
    nan_array = np.empty(len(stars_table))
    nan_array.fill(np.nan)
    star_index_columns = [
        'Field',
        'Name',
        'V_ref',
        'B-V',
        'U-B',
        'V-R',
        'V-I',
        'V_sigma',
        'e_B-V',
        'e_U-B',
        'e_V-R',
        'e_V-I'
    ]
    # try:
    exoatmospheric_table_begin = Table(
        names=[
            'Field',
            'Name',
            'V_ref',
            'B-V',
            'U-B',
            'V-R',
            'V-I',
            'V_sigma',
            'e_B-V',
            'e_U-B',
            'e_V-R',
            'e_V-I',
            # 'b B-V',
            # 'g B-V',
            # 'g V-R',
            # 'g V-I',
            # 'r V-R'
        ],
        data=[
            np.empty(len(stars_table), dtype=object),
            np.empty(len(stars_table), dtype=object),
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            # nan_array,
            # nan_array,
            # nan_array,
            # nan_array,
            # nan_array
        ])
    exoatmospheric_table_filter_ci_columns = []
    for instr_filter in different_filter_list:
        colour_indices_list = get_all_colour_indices(instr_filter)
        for colour_index in colour_indices_list:
            exoatmospheric_table_filter_ci_columns.append(
                f"{instr_filter} {colour_index}")
    exoatmospheric_table_filter_ci_data = np.empty(
        (len(nan_array), len(exoatmospheric_table_filter_ci_columns)))
    exoatmospheric_table_filter_ci_data.fill(np.nan)
    exoatmospheric_table_filter_ci =\
        Table(names=exoatmospheric_table_filter_ci_columns,
              data=exoatmospheric_table_filter_ci_data)
    exoatmospheric_table = hstack(
        (exoatmospheric_table_begin, exoatmospheric_table_filter_ci))
    i = 0
    for star in stars_table:
        for field in star_index_columns:
            exoatmospheric_table[field][i] = star[field]
        # exoatmospheric_table['Name'][i] = star['Name']
        for different_filter in different_filter_list:
            colour_indices = get_all_colour_indices(different_filter)
            x_column_name = f"X_{different_filter}"
            for colour_index in colour_indices:
                mask = ((extinction_table_Warner['filter'] == different_filter) & (
                    extinction_table_Warner['CI'] == colour_index))
                row_of_extinctions = extinction_table_Warner[mask]
                # if len(row_of_transforms) == 0 and instr_filter == 'v':
                #     instr_filter = 'g'
                #     mask = ((gb_final_transforms['filter'] == instr_filter) & (gb_final_transforms['CI'] == colour_index))
                #     row_of_transforms = gb_final_transforms[mask]
                instr_mag = star[different_filter]
                X = star[x_column_name]
                CI = star[colour_index]
                k_primeprime = row_of_extinctions['k\'\'_fCI']
                k_prime = row_of_extinctions['k\'_f']
                exoatmospheric_mag = float(
                    instr_mag - (k_prime * X) - (k_primeprime * X * CI))
                exoatmospheric_table[f"{different_filter} {colour_index}"][i] = exoatmospheric_mag
                # if different_filter == 'g':
                #     print(f"Instrumental mag for {different_filter} and {colour_index}:")
                #     print(f"{instr_mag:0.3f}")
                #     print(f"Exoatmospheric mag for {different_filter} and {colour_index}:")
                #     print(f"{exoatmospheric_mag:0.3f}")
        i += 1
    # except KeyError:
    #     exoatmospheric_table = Table(
    #         names=[
    #             'Field',
    #             'Name',
    #             'V_ref',
    #             'B-V',
    #             'U-B',
    #             'V-R',
    #             'V-I',
    #             'V_sigma',
    #             'e_B-V',
    #             'e_U-B',
    #             'e_V-R',
    #             'e_V-I',
    #             'u B-V',
    #             'g B-V',
    #             'g V-R',
    #             'g V-I',
    #             'r V-R',
    #             'i V-I'
    #         ],
    #         data=[
    #             np.empty(len(stars_table), dtype=object),
    #             np.empty(len(stars_table), dtype=object),
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array,
    #             nan_array
    #         ])
    #     i = 0
    #     for star in stars_table:
    #         for field in star_index_columns:
    #             exoatmospheric_table[field][i] = star[field]
    #         # exoatmospheric_table['Name'][i] = star['Name']
    #         for different_filter in different_filter_list:
    #             colour_indices = get_all_colour_indices(different_filter)
    #             x_column_name = f"X_{different_filter}"
    #             for colour_index in colour_indices:
    #                 mask =\
    #     ((extinction_table_Warner['filter'] == different_filter) & (
    #                             extinction_table_Warner['CI'] == colour_index))
    #                 row_of_extinctions = extinction_table_Warner[mask]
    #                 # if len(row_of_transforms) == 0 and instr_filter == 'v':
    #                 #     instr_filter = 'g'
    #                 #     mask = ((gb_final_transforms['filter'] == instr_filter) & (gb_final_transforms['CI'] == colour_index))
    #                 #     row_of_transforms = gb_final_transforms[mask]
    #                 instr_mag = star[different_filter]
    #                 X = star[x_column_name]
    #                 CI = star[colour_index]
    #                 k_primeprime = row_of_extinctions['k\'\'_fCI']
    #                 k_prime = row_of_extinctions['k\'_f']
    #                 exoatmospheric_mag = float(instr_mag - (k_prime * X) - (k_primeprime * X * CI))
    #                 exoatmospheric_table[f"{different_filter} {colour_index}"][i] = exoatmospheric_mag
    #                 # if different_filter == 'g':
    #                 #     print(f"Instrumental mag for {different_filter} and {colour_index}:")
    #                 #     print(f"{instr_mag:0.3f}")
    #                 #     print(f"Exoatmospheric mag for {different_filter} and {colour_index}:")
    #                 #     print(f"{exoatmospheric_mag:0.3f}")
    #         i += 1
    # exoatmospheric_table.pprint(max_lines=-1, max_width=250)
    return exoatmospheric_table


def colour_transform_and_zp_calc_Warner(exoatmospheric_table,
                                        different_filter_list,
                                        extinction_table_Warner,
                                        save_plots, **kwargs):
    nan_array = np.empty(len(extinction_table_Warner))
    nan_array.fill(np.nan)
    transform_zp_table = Table(
        names=[
            'T_fCI',
            'T_fCI_sigma',
            'Z_f',
            'Z_f_sigma'
        ],
        data=[
            nan_array,
            nan_array,
            nan_array,
            nan_array
        ])
    Warner_final_transform_table = hstack(
        [extinction_table_Warner, transform_zp_table])
    # print(Warner_final_transform_table)
    for different_filter in different_filter_list:
        try:
            app_mag, app_mag_sigma, app_mag_filter, _ = get_app_mag_and_index(
                exoatmospheric_table, different_filter)
        except KeyError:
            exoatmospheric_table.rename_column('V_sigma', 'e_V')
            app_mag, app_mag_sigma, app_mag_filter, _ = get_app_mag_and_index(
                exoatmospheric_table, different_filter)
        colour_indices = get_all_colour_indices(different_filter)
        for colour_index in colour_indices:
            ci_plot = np.arange(min(exoatmospheric_table[colour_index]) - 0.1,
                                max(exoatmospheric_table[colour_index]) + 0.1,
                                step=0.01)
            exoatmospheric_mag = np.array(
                exoatmospheric_table[f"{different_filter} {colour_index}"])
            y = app_mag - exoatmospheric_mag
            x = exoatmospheric_table[colour_index]
            fit, or_fit, line_init = init_linear_fitting(niter=100, sigma=2.5)
            # print(x)
            # print(y)
            # try:
            x_nan_indices = np.isnan(x)
            y_nan_indices = np.isnan(y)
            nan_indices = (x_nan_indices | y_nan_indices)
            fitted_line, mask = or_fit(
                line_init, x[~nan_indices], y[~nan_indices])
            # except TypeError as e:
            #     # print(current_star[x_current_filter])
            #     # print(current_star[unique_filter])
            #     # slopes_table[f"slope_{unique_filter}"][i] = np.nan
            #     # slopes_table[f"intercept_{unique_filter}"][i] = np.nan
            #     # slopes_table[f"slope_{unique_filter}_sigma"][i] = np.nan
            #     # slopes_table[f"intercept_{unique_filter}_sigma"][i] = np.nan
            #     # k_primeprime_column.append(np.nan)
            #     # k_prime_column.append(np.nan)
            #     # k_primeprime_sigma_column.append(np.nan)
            #     # k_prime_sigma_column.append(np.nan)
            #     print(e)
            #     continue
            filtered_data = np.ma.masked_array(y[~nan_indices], mask=mask)
            t_fci = fitted_line.slope.value
            z_f = fitted_line.intercept.value
            cov = fit.fit_info['param_cov']
            try:
                t_fci_sigma = sqrt(cov[0][0])
                z_f_sigma = sqrt(cov[1][1])
            except TypeError:
                t_fci_sigma = np.nan
                z_f_sigma = np.nan

            mask = ((Warner_final_transform_table['filter'] == different_filter) & (
                Warner_final_transform_table['CI'] == colour_index))
            # current_filter_index = Warner_final_transform_table[mask]
            Warner_final_transform_table["T_fCI"][mask] = t_fci
            Warner_final_transform_table["T_fCI_sigma"][mask] = t_fci_sigma
            Warner_final_transform_table["Z_f"][mask] = z_f
            Warner_final_transform_table["Z_f_sigma"][mask] = z_f_sigma
            plt.plot(x, y, 'o', fillstyle='none', label="Clipped Data")
            plt.plot(x[~nan_indices], filtered_data, 'o',
                     color='#1f77b4', label="Fitted Data")
            plt.plot(ci_plot, t_fci * ci_plot + z_f, '-',
                     label=f"t_fci={t_fci:0.3f}, ZP_f={z_f:0.3f}")
            plt.ylabel(f"{app_mag_filter} - {different_filter}$_0$")
            plt.xlabel(colour_index)
            plt.title(
                f"Colour Transform and Zero Point ({different_filter}$_0$ v. {colour_index})")
            plt.legend()
            if save_plots:
                save_loc = f"{os.path.join(kwargs.get('save_loc'), f'ColourTransformZeroPoint{different_filter}{colour_index}')}.png"
                plt.savefig(save_loc)
            plt.show()
            plt.close()
    return Warner_final_transform_table


def hidden_transform_Warner(exoatmospheric_table,
                            Warner_final_transform_table,
                            different_filter_list, save_plots,
                            **kwargs):
    # exoatmospheric_table = exoatmospheric_mags_Warner(stars_table, Warner_final_transform_table, different_filter_list)
    # exoatmospheric_table.pprint(max_lines=-1, max_width=-1)
    instr_ci_list = []
    app_ci_list = []
    t_ci = []
    t_ci_sigma = []
    zp_ci = []
    zp_ci_sigma = []
    for different_filter in different_filter_list:
        colour_index, ci = get_colour_index_lower(different_filter)
        try:
            ci0_index, _ = get_colour_index_lower(ci[0])
            ci1_index, _ = get_colour_index_lower(ci[1])
            positive_instr_mag = exoatmospheric_table[f"{ci[0]} {ci0_index}"]
            negative_instr_mag = exoatmospheric_table[f"{ci[1]} {ci1_index}"]
            table_ci = ci
        except KeyError:
            if 'v' in ci:
                table_ci = ci.replace('v', 'g')
            else:
                table_ci = ci
            try:
                ci0_index, _ = get_colour_index_lower(ci[0])
                ci1_index, _ = get_colour_index_lower(ci[1])
                positive_instr_mag =\
                    exoatmospheric_table[f"{table_ci[0]} {ci0_index}"]
                negative_instr_mag =\
                    exoatmospheric_table[f"{table_ci[1]} {ci1_index}"]
            except KeyError:
                if 'b' in ci:
                    table_ci = table_ci.replace('b', 'u')
                # table_ci = ci.replace('v', 'g')
                ci0_index, _ = get_colour_index_lower(ci[0])
                ci1_index, _ = get_colour_index_lower(ci[1])
                positive_instr_mag =\
                    exoatmospheric_table[f"{table_ci[0]} {ci0_index}"]
                negative_instr_mag =\
                    exoatmospheric_table[f"{table_ci[1]} {ci1_index}"]
        instr_colour_index_mags = positive_instr_mag - negative_instr_mag
        standard_colour_index_mags = exoatmospheric_table[colour_index]
        # plt.plot(instr_colour_index_mags, standard_colour_index_mags, 'o')
        plt.plot(instr_colour_index_mags, standard_colour_index_mags,
                 'o', fillstyle='none', label="Clipped Data")
        # m, b = np.polyfit(instr_colour_index_mags, standard_colour_index_mags, 1)
        fit, or_fit, line_init = init_linear_fitting(
            niter=100, sigma=2.5, slope=1)
        # try:
        x_nan_indices = np.isnan(instr_colour_index_mags)
        y_nan_indices = np.isnan(standard_colour_index_mags)
        nan_indices = (x_nan_indices | y_nan_indices)
        fitted_line, mask = or_fit(line_init, instr_colour_index_mags[~nan_indices],
                                   standard_colour_index_mags[~nan_indices])
        # except TypeError as e:
        #     print(e)
        #     continue
        filtered_data = np.ma.masked_array(
            standard_colour_index_mags[~nan_indices], mask=mask)
        plt.plot(instr_colour_index_mags[~nan_indices],
                 filtered_data, 'o', color='#1f77b4', label="Fitted Data")
        m = fitted_line.slope.value
        b = fitted_line.intercept.value
        cov = fit.fit_info['param_cov']
        try:
            m_sigma = sqrt(cov[0][0])
            b_sigma = sqrt(cov[1][1])
        except TypeError:
            m_sigma = np.nan
            b_sigma = np.nan
        if not (colour_index in app_ci_list and table_ci in instr_ci_list):
            app_ci_list.append(colour_index)
            instr_ci_list.append(table_ci)
            t_ci.append(m)
            t_ci_sigma.append(m_sigma)
            zp_ci.append(b)
            zp_ci_sigma.append(b_sigma)
        plt.plot(instr_colour_index_mags, m * instr_colour_index_mags +
                 b, '-', label=f"y={m:0.3f}x+{b:0.3f}")
        plt.ylabel(colour_index)
        plt.xlabel(f"{table_ci[0]}-{table_ci[1]}")
        plt.title("Hidden Transform")
        plt.legend()
        if save_plots:
            save_loc = f"{os.path.join(kwargs.get('save_loc'), f'HiddenTransform{colour_index}')}.png"
            plt.savefig(save_loc)
        plt.show()
        plt.close()
    hidden_transform_table = Table(
        names=[
            "Apparent CI",
            "Instrumental CI",
            "T_CI",
            "T_CI_sigma",
            "ZP_CI",
            "ZP_CI_sigma"
        ],
        data=[
            app_ci_list,
            instr_ci_list,
            t_ci,
            t_ci_sigma,
            zp_ci,
            zp_ci_sigma
        ]
    )
    return hidden_transform_table


def calculate_standard_CI_Warner(stars_table,
                                 hidden_transform_table,
                                 different_filter):
    colour_index, ci = get_colour_index_lower(different_filter)
    try:
        ci0_index, _ = get_colour_index_lower(ci[0])
        ci1_index, _ = get_colour_index_lower(ci[1])
        positive_instr_mag = stars_table[ci[0]]
        negative_instr_mag = stars_table[ci[1]]
        table_ci = ci
    except KeyError:
        if 'v' in ci:
            table_ci = ci.replace('v', 'g')
        else:
            table_ci = ci
        try:
            ci0_index, _ = get_colour_index_lower(ci[0])
            ci1_index, _ = get_colour_index_lower(ci[1])
            positive_instr_mag = stars_table[table_ci[0]]
            negative_instr_mag = stars_table[table_ci[1]]
        except KeyError:
            if 'b' in ci:
                table_ci = table_ci.replace('b', 'u')
            # table_ci = ci.replace('v', 'g')
            ci0_index, _ = get_colour_index_lower(ci[0])
            ci1_index, _ = get_colour_index_lower(ci[1])
            positive_instr_mag = stars_table[table_ci[0]]
            negative_instr_mag = stars_table[table_ci[1]]
    instr_colour_index_mags = positive_instr_mag - negative_instr_mag
    mask = ((hidden_transform_table['Apparent CI'] == colour_index) & (
        hidden_transform_table['Instrumental CI'] == table_ci))
    row_of_hidden_transforms = hidden_transform_table[mask]
    t_ci = row_of_hidden_transforms["T_CI"]
    zp_ci = row_of_hidden_transforms["ZP_CI"]
    CI = t_ci * instr_colour_index_mags + zp_ci
    return CI


def exoatmospheric_mags_verify_Warner(stars_table,
                                      extinction_table_Warner,
                                      hidden_transform_table,
                                      different_filter_list):
    # m_0 = m - k'_f * X - k''_fCI * X * CI
    nan_array = np.empty(len(stars_table))
    nan_array.fill(np.nan)
    star_index_columns = [
        'Field',
        'Name',
        'V_ref',
        'B-V',
        'U-B',
        'V-R',
        'V-I',
        'V_sigma',
        'e_B-V',
        'e_U-B',
        'e_V-R',
        'e_V-I'
    ]
    # exoatmospheric_table_verify = Table(
    #     names=[
    #         'Field',
    #         'Name',
    #         'V_ref',
    #         'B-V',
    #         'U-B',
    #         'V-R',
    #         'V-I',
    #         'V_sigma',
    #         'e_B-V',
    #         'e_U-B',
    #         'e_V-R',
    #         'e_V-I',
    #         'b',
    #         'g',
    #         'r'
    #     ],
    #     data=[
    #         np.empty(len(stars_table), dtype=object),
    #         np.empty(len(stars_table), dtype=object),
    #         nan_array,
    #         nan_array,
    #         nan_array,
    #         nan_array,
    #         nan_array,
    #         nan_array,
    #         nan_array,
    #         nan_array,
    #         nan_array,
    #         nan_array,
    #         nan_array,
    #         nan_array,
    #         nan_array
    #     ])
    exoatmospheric_table_begin = Table(
        names=[
            'Field',
            'Name',
            'V_ref',
            'B-V',
            'U-B',
            'V-R',
            'V-I',
            'V_sigma',
            'e_B-V',
            'e_U-B',
            'e_V-R',
            'e_V-I',
            # 'b B-V',
            # 'g B-V',
            # 'g V-R',
            # 'g V-I',
            # 'r V-R'
        ],
        data=[
            np.empty(len(stars_table), dtype=object),
            np.empty(len(stars_table), dtype=object),
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            # nan_array,
            # nan_array,
            # nan_array,
            # nan_array,
            # nan_array
        ])
    exoatmospheric_table_filter_columns = []
    for instr_filter in different_filter_list:
        exoatmospheric_table_filter_columns.append(f"{instr_filter}")
    exoatmospheric_table_filter_data = np.empty(
        (len(nan_array), len(exoatmospheric_table_filter_columns)))
    exoatmospheric_table_filter_data.fill(np.nan)
    exoatmospheric_table_filter =\
        Table(names=exoatmospheric_table_filter_columns,
              data=exoatmospheric_table_filter_data)
    exoatmospheric_table_verify = hstack(
        (exoatmospheric_table_begin, exoatmospheric_table_filter))
    i = 0
    for star in stars_table:
        for field in star_index_columns:
            exoatmospheric_table_verify[field][i] = star[field]
        # exoatmospheric_table['Name'][i] = star['Name']
        for different_filter in different_filter_list:
            colour_indices = get_all_colour_indices(different_filter)
            colour_index, _ = get_colour_index_lower(different_filter)
            x_column_name = f"X_{different_filter}"
            # for colour_index in colour_indices:
            mask = ((extinction_table_Warner['filter'] == different_filter) & (
                extinction_table_Warner['CI'] == colour_index))
            row_of_extinctions = extinction_table_Warner[mask]
            # if len(row_of_transforms) == 0 and instr_filter == 'v':
            #     instr_filter = 'g'
            #     mask = ((gb_final_transforms['filter'] == instr_filter) & (gb_final_transforms['CI'] == colour_index))
            #     row_of_transforms = gb_final_transforms[mask]
            instr_mag = star[different_filter]
            X = star[x_column_name]
            CI = calculate_standard_CI_Warner(
                star, hidden_transform_table, different_filter)
            # print(CI)
            # CI = np.array(CI)
            # print(CI)
            # CI = star[colour_index]
            k_primeprime = row_of_extinctions['k\'\'_fCI']
            k_prime = row_of_extinctions['k\'_f']
            try:
                exoatmospheric_mag = float(
                    instr_mag - (k_prime * X) - (k_primeprime * X * CI))
            except TypeError:
                continue
            exoatmospheric_table_verify[different_filter][i] =\
                exoatmospheric_mag
            # if different_filter == 'g':
            #     print(f"Instrumental mag for {different_filter} and {colour_index}:")
            #     print(f"{instr_mag:0.3f}")
            #     print(f"Exoatmospheric mag for {different_filter} and {colour_index}:")
            #     print(f"{exoatmospheric_mag:0.3f}")
        i += 1
    # exoatmospheric_table.pprint(max_lines=-1, max_width=250)
    return exoatmospheric_table_verify


def apply_transforms_Warner(stars_table,
                            exoatmospheric_table_verify,
                            Warner_final_transform_table,
                            hidden_transform_table,
                            different_filter_list,
                            save_plots,
                            **kwargs):
    for different_filter in different_filter_list:
        app_mag, app_mag_sigma, app_filter, colour_index = get_app_mag_and_index_AVG(
            stars_table, different_filter)
        app_mag_filter = different_filter.upper()
        CI = calculate_standard_CI_Warner(
            stars_table, hidden_transform_table, different_filter)
        exoatmospheric_mag = exoatmospheric_table_verify[different_filter]
        mask = ((Warner_final_transform_table['filter'] == different_filter) & (
            Warner_final_transform_table['CI'] == colour_index))
        row_of_extinctions = Warner_final_transform_table[mask]
        # print(row_of_extinctions)
        t_fci = float(row_of_extinctions["T_fCI"])
        z_f = float(row_of_extinctions["Z_f"])
        # print(exoatmospheric_mag)
        # print(CI)
        # print(z_f)
        app_mag_calculated = exoatmospheric_mag + (t_fci * CI) + z_f
        fit, or_fit, line_init = init_linear_fitting(
            niter=100, sigma=2.5, slope=1, intercept=0)
        # try:
        x_nan_indices = np.isnan(app_mag)
        y_nan_indices = np.isnan(app_mag_calculated)
        nan_indices = (x_nan_indices | y_nan_indices)
        fitted_line, mask = or_fit(
            line_init, app_mag[~nan_indices], app_mag_calculated[~nan_indices])
        # except TypeError as e:
        #     print(e)
        #     continue
        filtered_data = np.ma.masked_array(
            app_mag_calculated[~nan_indices], mask=mask)
        plt.plot(app_mag, app_mag_calculated, 'o',
                 fillstyle='none', label="Clipped Data")
        plt.plot(app_mag[~nan_indices], filtered_data,
                 'o', color='#1f77b4', label="Fitted Data")
        m = fitted_line.slope.value
        b = fitted_line.intercept.value
        cov = fit.fit_info['param_cov']
        try:
            m_sigma = sqrt(cov[0][0])
            b_sigma = sqrt(cov[1][1])
        except TypeError:
            m_sigma = np.nan
            b_sigma = np.nan
        # if not (colour_index in app_ci_list and table_ci in instr_ci_list):
        #     app_ci_list.append(colour_index)
        #     instr_ci_list.append(table_ci)
        #     t_ci.append(m)
        #     t_ci_sigma.append(m_sigma)
        #     zp_ci.append(b)
        #     zp_ci_sigma.append(b_sigma)
        plt.plot(app_mag, m * app_mag + b, '-', label=f"y={m:0.3f}x+{b:0.3f}")
        plt.plot(app_mag, app_mag, label=f"y=x")
        plt.ylabel(f"{app_mag_filter} (Calculated)")
        plt.xlabel(f"{app_filter} (Reference)")
        plt.title("Calculated Magnitude vs. Reference Magnitude")
        plt.legend()
        if save_plots:
            if not os.path.exists(os.path.join(kwargs.get('save_loc'),
                                               'VERIFICATION')):
                os.mkdir(os.path.join(kwargs.get('save_loc'),
                                      'VERIFICATION'))
            save_loc = f"{os.path.join(kwargs.get('save_loc'), 'VERIFICATION', f'{app_mag_filter}_CalcVsRefMag')}.png"
            plt.savefig(save_loc)
        plt.show()
        plt.close()


def verify_gb_transforms_auto(directory,
                              stars_list,
                              ref_stars_file,
                              gb_final_transforms,
                              hidden_transform_table,
                              plot_results=False,
                              save_plots=False,
                              file_suffix=(".fits", ".fit", ".fts"),
                              exposure_key='EXPTIME',
                              name_key='Name',
                              lat_key='SITELAT',
                              lon_key='SITELONG',
                              elev_key='SITEELEV',
                              **kwargs):

    # TODO: Docstring.
    reference_stars, ref_star_positions = read_ref_stars(ref_stars_file)
    large_table_columns = init_large_table_columns()

    if save_plots:
        save_loc = kwargs.get('save_loc')
        # unique_id = kwargs.get('unique_id')
        if not os.path.exists(save_loc):
            os.mkdir(save_loc)

    for file_num, filepath in enumerate(tqdm(stars_list)):
        # filepath = os.path.join(dirpath, filename)
        hdr, imgdata = read_fits_file(filepath)
        exptime = hdr[exposure_key]
        bkg, bkg_std = calculate_img_bkg(imgdata)
        irafsources = detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
        if not irafsources:
            continue
        _, fwhm, fwhm_std = calculate_fwhm(irafsources)
        photometry_result = perform_photometry(
            irafsources, fwhm, imgdata, bkg=bkg)
        fluxes = np.array(photometry_result['flux_fit'])
        instr_mags = calculate_magnitudes(photometry_result, exptime)
        instr_mags_sigma = calculate_magnitudes_sigma(
            photometry_result, exptime)
        wcs = WCS(hdr)
        skypositions = convert_pixel_to_ra_dec(irafsources, wcs)
        try:
            altazpositions = convert_ra_dec_to_alt_az(skypositions,
                                                      hdr,
                                                      lat_key=lat_key,
                                                      lon_key=lon_key,
                                                      elev_key=elev_key)
            # altazpositions = convert_ra_dec_to_alt_az(skypositions, hdr,
            # lat_key='SITELAT', lon_key='SITELONG',
            #                                           elev_key='SITEELEV')
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

        large_table_columns = update_large_table_columns(large_table_columns,
                                                         filepath,
                                                         matched_stars,
                                                         hdr,
                                                         exptime,
                                                         ground_based=True,
                                                         name_key=name_key)
    large_stars_table = create_large_stars_table(
        large_table_columns, ground_based=True)
    # large_stars_table = remove_large_airmass(large_stars_table,
    # max_airmass=2.0)
    stars_table, different_filter_list = group_each_star_GB(large_stars_table)
    # print(different_filter_list)
    if not hidden_transform_table:
        pass
    else:
        exoatmospheric_table_verify = exoatmospheric_mags_verify_Warner(
            stars_table,
            gb_final_transforms,
            hidden_transform_table,
            different_filter_list)
        warner_verfication = apply_transforms_Warner(
            stars_table,
            exoatmospheric_table_verify,
            gb_final_transforms,
            hidden_transform_table,
            different_filter_list,
            save_plots=save_plots,
            save_loc=save_loc)
    stars_table.pprint(max_lines=30, max_width=200)

    # instr_filters = ['b', 'v', 'r', 'i']
    app_mag_table = Table(stars_table[
        'Field', 'Name', 'V_ref', 'B-V', 'U-B', 'V-R',
        'V-I', 'V_sigma', 'e_B-V', 'e_U-B', 'e_V-R', 'e_V-I'])
    for instr_filter in different_filter_list:
        app_mag_filter = instr_filter.upper()
        app_filter_sigma = f"{app_mag_filter}_sigma"
        app_mag_table_filter = apply_gb_transforms_VERIFICATION(
            gb_final_transforms, stars_table, instr_filter)
        app_mag_table = hstack(
            [app_mag_table, app_mag_table_filter[app_mag_filter,
                                                 app_filter_sigma]])

    app_mag_table.pprint(max_lines=30, max_width=200)
    import matplotlib.pyplot as plt
    # import matplotlib
    # matplotlib.use('TkAgg')
    for instr_filter in different_filter_list:
        app_filter = instr_filter.upper()
        try:
            ref_app_mag, ref_app_mag_sigma, _, _ = get_app_mag_and_index_AVG(
                app_mag_table, instr_filter)
        except KeyError:
            ref_app_mag, ref_app_mag_sigma, _, _ = get_app_mag_and_index(
                app_mag_table, instr_filter)
        app_filter_sigma = f"{app_filter}_sigma"
        x = ref_app_mag
        x_err = ref_app_mag_sigma
        y = app_mag_table[app_filter]
        # y_err = app_mag_table[app_filter_sigma]
        max_y_err = max(app_mag_table[app_filter_sigma])
        # print(max_y_err)
        y_err = np.nan_to_num(app_mag_table[app_filter_sigma], nan=max_y_err)
        y_err = np.array(y_err)
        y_err[y_err == 0] = max_y_err
        fit, or_fit, line_init = init_linear_fitting(sigma=2.5)
        x_nan_indices = np.isnan(x)
        y_nan_indices = np.isnan(y)
        nan_indices = (x_nan_indices | y_nan_indices)
        # try:
        # , weights=1.0/y_err[~nan_indices])
        fitted_line, mask = or_fit(line_init, x[~nan_indices], y[~nan_indices])
        filtered_data = np.ma.masked_array(y[~nan_indices], mask=mask)
        m = fitted_line.slope.value
        b = fitted_line.intercept.value
        # fitted_line, mask = or_fit(line_init, x, y, weights=1.0/y_err)
        # filtered_data = np.ma.masked_array(y, mask=mask)
        # m = fitted_line.slope.value
        # b = fitted_line.intercept.value
        # if m == 1 and b == 0:
        #     fitted_line, mask = or_fit(line_init, x, y)
        #     filtered_data = np.ma.masked_array(y, mask=mask)
        #     m = fitted_line.slope.value
        #     b = fitted_line.intercept.value
        # except TypeError as e:
        #     print(e)
        #     continue
        plt.errorbar(x[~nan_indices],
                     y[~nan_indices],
                     yerr=y_err[~nan_indices],
                     fmt='o',
                     fillstyle='none',
                     capsize=0,
                     label="Clipped Data")
        plt.plot(x[~nan_indices], filtered_data, 'o',
                 color='#1f77b4', label="Fitted Data")
        plt.plot(x[~nan_indices], m * x[~nan_indices] +
                 b, '-', label=f'y={m:.3f}x+{b:.3f}')
        plt.plot(x[~nan_indices], x[~nan_indices], '-', label='y=x')
        plt.title('Calculated Magnitude vs. Reference Magnitude')
        plt.ylabel(f"{app_filter} (Calculated)")
        plt.xlabel(f"{app_filter} (Reference)")
        plt.legend()
        if save_plots:
            save_filename = f"{os.path.join(save_loc, f'{app_filter}_CalcVsRefMag')}.png"
            ascii.write(app_mag_table, os.path.join(
                save_loc, 'app_mag_table.csv'), format='csv')
            # print(save_filename)
        plt.savefig(save_filename)
        if plot_results:
            plt.show(block=True)
        plt.close()
    return app_mag_table


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
                                   **kwargs):
    # TODO: Docstring.
    # TODO: Fix errors when save_plots = False.
    """
    Perform all of the beginning operations.
    Create the refrence stars table and read their positions.
    Initialize empty arrays for the star information and auxiliary info.
    """
    reference_stars, ref_star_positions = read_ref_stars(ref_stars_file)
    large_table_columns = init_large_table_columns()
    star_aux_table_columns = init_star_aux_table_columns()

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

    # What is this?
    shuffle(file_paths)
    split_decimal = 1
    split_filecount_location = math.ceil(split_decimal * filecount)
    calculation_files = file_paths[:split_filecount_location]
    verification_files = file_paths[split_filecount_location:]
    with open(os.path.join(save_loc, 'CalculationVerificationSplit.txt'), 'a') as f:
        f.write('File Path'+'\t'+'Calculation/Verification')
        for calc_file in calculation_files:
            f.write('\n'+f'{calc_file}'+'\t'+'Calculation')
        for verify_file in verification_files:
            f.write('\n'+f'{verify_file}'+'\t'+'Verification')
    "Iterate over the images."
    for file_num, filepath in enumerate(tqdm(calculation_files)):
        # Read the fits file. Stores the header and image to variables.
        hdr, imgdata = read_fits_file(filepath)
        # Read the exposure time of the image.
        exptime = hdr[exposure_key]
        # Calculate the image background and standard deviation.
        bkg, bkg_std = calculate_img_bkg(imgdata)
        # Detect point sources in the image.
        irafsources = detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
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
        fwhms, fwhm, fwhm_std = calculate_fwhm(irafsources)
        # Do PSF photometry on the detected sources.
        photometry_result = perform_photometry(
            irafsources, fwhm, imgdata, bkg=bkg)
        # Store the flux of the stars in a separate variable.
        fluxes = np.array(photometry_result['flux_fit'])
        # Convert the flux and uncertainty to magnitude and its uncertainty.
        instr_mags = calculate_magnitudes(photometry_result, exptime)
        instr_mags_sigma = calculate_magnitudes_sigma(
            photometry_result, exptime)
        # Read the World Coordinate System transformation added to the fits header
        # by a plate solving software (external to this program, e.g. PinPoint).

        # FIXME : Existing WCS Data must be present to plate solve.
        # PinPoint Should be called and solved within the program

        wcs = WCS(hdr)
        # Convert the stars' (x,y) location to (RA,dec).
        skypositions = convert_pixel_to_ra_dec(irafsources, wcs)
        try:
            # altazpositions = convert_ra_dec_to_alt_az(skypositions, hdr, lat_key='OBSGEO-B', lon_key='OBSGEO-L',
            #                                           elev_key='OBSGEO-H')

            # FIXME: SkyCoord System does not work

            # Convert the stars' (RA,dec) location to (Azimuth,Elevation).
            altazpositions = convert_ra_dec_to_alt_az(skypositions, hdr, lat_key=lat_key, lon_key=lon_key,
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
        fwhms_arcsec, fwhm_arcsec, fwhms_arcsec_std = convert_fwhm_to_arcsec(
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
        background_sky_brightness = calculate_background_sky_brightness(
            bkg, hdr, exptime)
        background_sky_brightness_sigma = calculate_BSB_sigma(
            bkg, bkg_std, exptime)
        # Take the average of all stars' Az/El/airmass and store as a variable.
        azimuth = np.mean(altazpositions.az)
        elevation = np.mean(altazpositions.alt)
        airmass = np.mean(altazpositions.secz)
        # Update the table with auxiliary data on the images (FWHM, BSB, etc.)
        star_aux_table_columns =\
            update_star_aux_columns(star_aux_table_columns,
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
        matched_stars = find_ref_stars(reference_stars,
                                       ref_star_positions,
                                       skypositions,
                                       instr_mags,
                                       instr_mags_sigma,
                                       fluxes,
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
        large_table_columns =\
            update_large_table_columns(large_table_columns,
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
    star_aux_table = create_star_aux_table(star_aux_table_columns)
    ascii.write(star_aux_table, os.path.join(
        save_loc, 'auxiliary_table.csv'), format='csv')
    # Create an AstroPy table of each reference star detection and write it to a .csv file.
    large_stars_table = create_large_stars_table(
        large_table_columns, ground_based=True)
    ascii.write(large_stars_table, os.path.join(
        save_loc, 'large_stars_table.csv'), format='csv')
    # Group each observation of a star at an airmass.
    # E.g. if there are 5 images of star X at 1.2 airmass, and 10 images of star X at 2 airmass,
    # it will produce a mean and standard deviation of the observations at both 1.2 and 2 airmass.
    # This creates that table, stores the different filters used to take the images (e.g. BVRI or BGR),
    # and writes it to a .csv file.
    stars_table, different_filter_list = group_each_star_GB(large_stars_table)
    ascii.write(stars_table, os.path.join(
        save_loc, 'stars_table.csv'), format='csv')

    # TRYING SOMETHING

    # reformatted_large_stars_table = create_reformatted_large_table(large_stars_table, keys='Name')
    # ascii.write(reformatted_large_stars_table, os.path.join(save_loc, 'reformatted_large_stars_table.csv'), format='csv')

    ############# Begin the Warner Transforms #############

    # Calculate the slope of each star's instrumental magnitude vs airmass,
    # store it in a table, and write it to a .csv file.
    slopes_table = calculate_slopes_Warner(
        stars_table, different_filter_list, save_plots, save_loc=save_loc)
    ascii.write(slopes_table, os.path.join(
        save_loc, 'slopes_table.csv'), format='csv')
    # Calculate the first and second order extinctions.
    extinction_table_Warner =\
        second_order_extinction_calc_Warner(
            slopes_table,
            different_filter_list,
            save_plots,
            save_loc=save_loc)
    # Calculate the exoatmospheric magnitudes (m_0).
    exoatmospheric_table = exoatmospheric_mags_Warner(
        stars_table, extinction_table_Warner, different_filter_list)
    # Finish the transform by calculating the colour transform and zero point.
    Warner_final_transform_table =\
        colour_transform_and_zp_calc_Warner(
            exoatmospheric_table,
            different_filter_list,
            extinction_table_Warner, save_plots,
            save_loc=save_loc)
    # Save the transform table to a .csv file.
    ascii.write(Warner_final_transform_table, os.path.join(
        save_loc, '_gb_final_transforms.csv'), format='csv')
    # Calculate the hidden transform and write it to a .csv file.
    hidden_transform_table =\
        hidden_transform_Warner(
            exoatmospheric_table,
            Warner_final_transform_table,
            different_filter_list,
            save_plots,
            save_loc=save_loc)
    ascii.write(hidden_transform_table, os.path.join(
        save_loc, 'hidden_transform_table.csv'), format='csv')
    exoatmospheric_table_verify =\
        exoatmospheric_mags_verify_Warner(
            stars_table,
            extinction_table_Warner,
            hidden_transform_table, different_filter_list)
    # Verify the transforms.
    verify_save_loc = os.path.join(save_loc, 'Verification')
    app_mag_table = verify_gb_transforms_auto(directory,
                                              verification_files,
                                              ref_stars_file,
                                              Warner_final_transform_table,
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
    return Warner_final_transform_table


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
                                     **kwargs):

    # TODO: Docstring.
    # TODO: Fix errors when save_plots = False.
    """
    Perform all of the beginning operations.
    Create the refrence stars table and read their positions.
    Initialize empty arrays for the star information and auxiliary info.
    """
    reference_stars, ref_star_positions = read_ref_stars(ref_stars_file)
    large_table_columns = init_large_table_columns()
    star_aux_table_columns = init_star_aux_table_columns()

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
                           'CalculationVerificationSplit.txt'),
              'a') as f:
        f.write('File Path'+'\t'+'Calculation/Verification')
        for calc_file in calculation_files:
            f.write('\n'+f'{calc_file}'+'\t'+'Calculation')
        for verify_file in verification_files:
            f.write('\n'+f'{verify_file}'+'\t'+'Verification')
    "Iterate over the images."
    for file_num, filepath in enumerate(tqdm(calculation_files)):
        # Read the fits file. Stores the header and image to variables.
        hdr, imgdata = read_fits_file(filepath)
        # Read the exposure time of the image.
        exptime = hdr[exposure_key]
        # Calculate the image background and standard deviation.
        bkg, bkg_std = calculate_img_bkg(imgdata)
        # Detect point sources in the image.
        irafsources = detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
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
        fwhms, fwhm, fwhm_std = calculate_fwhm(irafsources)
        # Do PSF photometry on the detected sources.
        photometry_result = perform_photometry(
            irafsources, fwhm, imgdata, bkg=bkg)
        # Store the flux of the stars in a separate variable.
        fluxes = np.array(photometry_result['flux_fit'])
        # Convert the flux and uncertainty to magnitude and its uncertainty.
        instr_mags = calculate_magnitudes(photometry_result, exptime)
        instr_mags_sigma = calculate_magnitudes_sigma(
            photometry_result, exptime)
        # Read the World Coordinate System transformation added to the
        # fits header
        # by a plate solving software (external to this program, e.g. PinPoint).
        wcs = WCS(hdr)
        # Convert the stars' (x,y) centroid locations to (RA,dec).
        skypositions = convert_pixel_to_ra_dec(irafsources, wcs)
        try:
            # altazpositions = convert_ra_dec_to_alt_az(skypositions, hdr,
            # lat_key='OBSGEO-B', lon_key='OBSGEO-L',
            #                                           elev_key='OBSGEO-H')

            # Convert the stars' (RA,dec) location to (Azimuth,Elevation).
            altazpositions = convert_ra_dec_to_alt_az(skypositions, hdr,
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
        fwhms_arcsec, fwhm_arcsec, fwhms_arcsec_std = convert_fwhm_to_arcsec(
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
        background_sky_brightness = calculate_background_sky_brightness(
            bkg, hdr, exptime)
        background_sky_brightness_sigma = calculate_BSB_sigma(
            bkg, bkg_std, exptime)
        # Take the average of all stars' Az/El/airmass and store as a variable.
        azimuth = np.mean(altazpositions.az)
        elevation = np.mean(altazpositions.alt)
        airmass = np.mean(altazpositions.secz)
        # Update the table with auxiliary data on the images (FWHM, BSB, etc.)
        star_aux_table_columns =\
            update_star_aux_columns(star_aux_table_columns,
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
        matched_stars = find_ref_stars(reference_stars,
                                       ref_star_positions,
                                       skypositions,
                                       instr_mags,
                                       instr_mags_sigma,
                                       fluxes,
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
        large_table_columns = update_large_table_columns(large_table_columns,
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
    star_aux_table = create_star_aux_table(star_aux_table_columns)
    ascii.write(star_aux_table, os.path.join(
        save_loc, 'auxiliary_table.csv'), format='csv')
    # Create an AstroPy table of each reference star detection and write it to a .csv file.
    large_stars_table = create_large_stars_table(
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
    stars_table, different_filter_list = group_each_star_GB(large_stars_table)
    ascii.write(stars_table, os.path.join(
        save_loc, 'stars_table.csv'), format='csv')

    ############# Begin the Buchheim Transforms #############

    # Calculate the slope of each star's instrumental magnitude vs airmass,
    # store it in a table,
    # and write it to a .csv file.
    slopes_table = calculate_slopes_Buchheim(
        stars_table, different_filter_list, save_plots, save_loc=save_loc)
    ascii.write(slopes_table, os.path.join(
        save_loc, 'slopes_table.csv'), format='csv')
    # Calculate the first and second order extinctions.
    second_order_extinction_calc_Buchheim(stars_table, different_filter_list,
                                          save_plots,
                                          save_loc=save_loc)
    extinction_table_Buckhheim = \
        extinction_calc_Buchheim_sect6(slopes_table, different_filter_list,
                                       save_plots, save_loc=save_loc)
    ascii.write(extinction_table_Buckhheim, os.path.join(
        save_loc, 'extinction_table_Buckhheim.csv'), format='csv')
    return
    # Calculate the exoatmospheric magnitudes (m_0).
    exoatmospheric_table = exoatmospheric_mags_Warner(
        stars_table, extinction_table_Buckhheim, different_filter_list)

    ############# Begin the Warner Transforms #############

    # Finish the transform by calculating the colour transform and zero point.
    Buckhheim_final_transform_table =\
        colour_transform_and_zp_calc_Warner(exoatmospheric_table,
                                            different_filter_list,
                                            extinction_table_Buckhheim,
                                            save_plots, save_loc=save_loc)
    # Save the transform table to a .csv file.
    ascii.write(Buckhheim_final_transform_table, os.path.join(
        save_loc, '_gb_final_transforms.csv'), format='csv')
    # Calculate the hidden transform and write it to a .csv file.
    hidden_transform_table =\
        hidden_transform_Warner(exoatmospheric_table,
                                Buckhheim_final_transform_table,
                                different_filter_list,
                                save_plots,
                                save_loc=save_loc)
    ascii.write(hidden_transform_table, os.path.join(
        save_loc, 'hidden_transform_table.csv'), format='csv')
    exoatmospheric_table_verify =\
        exoatmospheric_mags_verify_Warner(stars_table,
                                          extinction_table_Buckhheim,
                                          hidden_transform_table,
                                          different_filter_list)
    # Verify the transforms.
    verify_save_loc = os.path.join(save_loc, 'Verification')
    app_mag_table = verify_gb_transforms_auto(directory,
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


def _sky_survey_calc(directory,
                     plot_results=False,
                     save_plots=False,
                     file_suffix=(".fits", ".fit", ".fts"),
                     exposure_key='EXPTIME',
                     gb_final_transforms=None,
                     lat_key='SITELAT',
                     lon_key='SITELONG',
                     elev_key='SITEELEV',
                     **kwargs):
    # TODO: Docstring.
    # TODO: Fix errors when save_plots = False.
    # large_table_columns = init_large_table_columns()
    star_aux_table_columns = init_star_aux_table_columns()

    filecount = 0
    file_paths = []
    file_names = []
    for dirpth, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_suffix):

                file_paths.append(os.path.join(dirpth, file))
                file_names.append(file)
    #             filecount += 1
    # filenames = sorted(os.listdir(temp_dir))

    save_loc = kwargs.get('save_loc')
    if not os.path.exists(save_loc):
        os.mkdir(save_loc)

    # for dirpath, dirnames, filenames in os.walk(directory):
    #     for filename in tqdm(filenames):
    #         if filename.endswith(file_suffix):
    # TODO: Edit this line to properly read the checkpoint
    # [100:]), start=100):
    for file_num, filepath in enumerate(tqdm(file_paths)):
        # filepath = os.path.join(dirpath, filename)
        hdr, imgdata = read_fits_file(filepath)
        exptime = hdr[exposure_key]
        bkg, bkg_std = calculate_img_bkg(imgdata)
        # print(bkg)
        irafsources = detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
        if not irafsources:
            warn("Could not detect sources.")
            fwhm = np.nan
            fwhm_std = np.nan
            fwhm_arcsec = np.nan
            fwhm_arcsec_std = np.nan
            num_sources = 0
            t = Time(hdr['DATE-OBS'], format='fits', scale='utc')
            time = t.jd
            # time = hdr['DATE-OBS']
            img_filter = hdr['FILTER']
            background_sky_brightness = calculate_background_sky_brightness(
                bkg, hdr, exptime, gb_final_transforms)
            background_sky_brightness_sigma = calculate_BSB_sigma(
                bkg, bkg_std, exptime)
            try:
                azimuth = hdr['CENTAZ']
                elevation = hdr['CENTALT']
                airmass = hdr['AIRMASS']
            except KeyError:
                continue
            star_aux_table_columns =\
                update_star_aux_columns(star_aux_table_columns,
                                        file_names[file_num],
                                        time,
                                        img_filter,
                                        fwhm,
                                        fwhm_std,
                                        fwhm_arcsec,
                                        fwhm_arcsec_std,
                                        num_sources,
                                        background_sky_brightness,
                                        background_sky_brightness_sigma,
                                        azimuth,
                                        elevation,
                                        airmass)
            continue
        num_sources = len(irafsources)
        fwhms, fwhm, fwhm_std = calculate_fwhm(irafsources)
        # photometry_result = perform_photometry(irafsources, fwhm, imgdata, bkg=bkg)
        # fluxes = np.array(photometry_result['flux_fit'])
        # instr_mags = calculate_magnitudes(photometry_result, exptime)
        # instr_mags_sigma = calculate_magnitudes_sigma(photometry_result, exptime)
        try:
            wcs = WCS(hdr)
            skypositions = convert_pixel_to_ra_dec(irafsources, wcs)
            try:
                altazpositions = convert_ra_dec_to_alt_az(skypositions,
                                                          hdr,
                                                          lat_key=lat_key,
                                                          lon_key=lon_key,
                                                          elev_key=elev_key)
                azimuth = np.mean(altazpositions.az)
                elevation = np.mean(altazpositions.alt)
                airmass = np.mean(altazpositions.secz)
            except AttributeError as e:
                azimuth = hdr['CENTAZ']
                elevation = hdr['CENTALT']
                airmass = hdr['AIRMASS']
                # continue
        except Exception:
            # TODO: change to an if/else statement.
            try:
                azimuth = hdr['CENTAZ']
                elevation = hdr['CENTALT']
                airmass = hdr['AIRMASS']
            except KeyError:
                continue
        fwhms_arcsec, fwhm_arcsec, fwhm_arcsec_std = convert_fwhm_to_arcsec(
            hdr, fwhms, fwhm, fwhm_std)
        t = Time(hdr['DATE-OBS'], format='fits', scale='utc')
        time = t.jd
        # time = hdr['DATE-OBS']
        img_filter = hdr['FILTER']
        background_sky_brightness = calculate_background_sky_brightness(
            bkg, hdr, exptime, gb_final_transforms)
        background_sky_brightness_sigma = calculate_BSB_sigma(
            bkg, bkg_std, exptime)
        # azimuth = hdr['CENTAZ']
        # elevation = hdr['CENTALT']
        # airmass = hdr['AIRMASS']
        star_aux_table_columns =\
            update_star_aux_columns(star_aux_table_columns,
                                    file_names[file_num],
                                    time,
                                    img_filter,
                                    fwhm,
                                    fwhm_std,
                                    fwhm_arcsec,
                                    fwhm_arcsec_std,
                                    num_sources,
                                    background_sky_brightness,
                                    background_sky_brightness_sigma,
                                    azimuth,
                                    elevation,
                                    airmass)
        # if (file_num % 10) == 0:
        #     star_aux_table = create_star_aux_table(star_aux_table_columns)
        #     ascii.write(star_aux_table, os.path.join(save_loc, 'auxiliary_table.csv'), format='csv')
        #     with open(os.path.join(save_loc, 'checkpoint.txt'), 'a') as f:
        #         f.write(str(file_num))
        #         f.write('\n')
        #         f.write(filepath)
        #         f.write('\n')
    star_aux_table = create_star_aux_table(star_aux_table_columns)
    ascii.write(star_aux_table, os.path.join(
        save_loc, 'auxiliary_table.csv'), format='csv')

    theta = star_aux_table['Azimuth'][star_aux_table['BSB'] > 5]
    r = star_aux_table['Elevation'][star_aux_table['BSB'] > 5]
    z = star_aux_table['BSB'][star_aux_table['BSB'] > 5]
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(7, 7))
    ax.set_theta_zero_location("N")
    norm = matplotlib.colors.Normalize(vmin=np.percentile(
        z[~np.isnan(z)], 7), vmax=max(z[~np.isnan(z)]))
    m = cm.ScalarMappable(cmap=plt.get_cmap('viridis_r'), norm=norm)
    m.set_array([])
    plt.colorbar(m)
    # Change contourf in the line below to scatter if you have only 1D theta,
    # r and brightness values
    ax.scatter(theta[~np.isnan(z)], r[~np.isnan(z)], c=z[~np.isnan(
        z)], cmap=plt.get_cmap('viridis_r'), norm=norm)
    rlabels = ax.get_ymajorticklabels()
    ax.set_rlim(bottom=90, top=15)
    for label in rlabels:
        label.set_color('black')
    plt.savefig(os.path.join(save_loc, 'BSB_plot.png'))
    plt.show()
    plt.close()

    fwhm_arcsec = star_aux_table['FWHM_arcsec']
    fwhm_arcsec_sigma = star_aux_table['FWHM_arcsec_sigma']
    times_list = np.array(star_aux_table['Time (JD)'])
    times_obj = Time(times_list, format='jd', scale='utc')
    times_datetime = times_obj.to_value('datetime')
    fig, ax = plt.subplots()
    _, _, bars = ax.errorbar(times_datetime,
                             fwhm_arcsec,
                             yerr=fwhm_arcsec_sigma,
                             fmt='o',
                             markersize=2,
                             capsize=0,
                             elinewidth=0.75)
    [bar.set_alpha(0.3) for bar in bars]
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.set_ylabel('FWHM (arcsec)')
    ax.set_xlabel('Time (UTC)')
    plt.title('FWHM (arcsec) v. Time')
    plt.savefig(os.path.join(save_loc, 'FWHM_arcsec.png'))
    plt.show()
    plt.close()

    fwhm = star_aux_table['FWHM_pixel']
    fwhm_sigma = star_aux_table['FWHM_pixel_sigma']
    fig, ax = plt.subplots()
    _, _, bars = ax.errorbar(times_datetime,
                             fwhm,
                             yerr=fwhm_sigma,
                             fmt='o',
                             markersize=2,
                             capsize=0,
                             elinewidth=0.75)
    [bar.set_alpha(0.3) for bar in bars]
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.set_ylabel('FWHM (pixels)')
    ax.set_xlabel('Time (UTC)')
    plt.title('FWHM (pixels) v. Time')
    plt.savefig(os.path.join(save_loc, 'FWHM.png'))
    plt.show()
    plt.close()

    return star_aux_table


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
                irafsources = detecting_stars(
                    imgdata, bkg=bkg, bkg_std=bkg_std)
                if not irafsources:
                    continue
                _, fwhm, fwhm_std = calculate_fwhm(irafsources)
                photometry_result = perform_photometry(
                    irafsources, fwhm, imgdata, bkg=bkg)
                fluxes = np.array(photometry_result['flux_fit'])
                instr_mags = calculate_magnitudes(photometry_result, exptime)
                instr_mags_sigma = calculate_magnitudes_sigma(
                    photometry_result, exptime)
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

                large_table_columns =\
                    update_large_table_columns(large_table_columns,
                                               filepath,
                                               matched_stars,
                                               hdr,
                                               exptime,
                                               ground_based=False,
                                               name_key=name_key)
    large_stars_table = create_large_stars_table(
        large_table_columns, ground_based=False)
    stars_table = group_each_star(large_stars_table, ground_based=False)
    sb_final_transform_columns = init_sb_final_transform_columns()
    if save_plots:
        write_table_to_latex(stars_table, f"{os.path.join(save_loc, f'{unique_id}_stars_table')}.txt",
                             formats={'c': '%0.3f',
                                      'c_sigma': '%0.3f'})
        for index in transform_index_list:
            filter_fci,\
                filter_fci_sigma,\
                zprime_fci,\
                zprime_fci_sigma\
                = space_based_transform(stars_table,
                                        plot_results=plot_results,
                                        index=index,
                                        save_plots=save_plots,
                                        save_loc=save_loc,
                                        unique_id=unique_id)
            sb_final_transform_columns =\
                update_sb_final_transform_columns(sb_final_transform_columns,
                                                  index,
                                                  filter_fci,
                                                  filter_fci_sigma,
                                                  zprime_fci,
                                                  zprime_fci_sigma)
            # print(f"(V-clear) = ({filter_fci:.3f} +/- {filter_fci_sigma:.3f}) * {index} + " \
            #       f"({zprime_fci:.3f} +/- {zprime_fci_sigma:.3f})")
    else:
        for index in transform_index_list:
            filter_fci,\
                filter_fci_sigma,\
                zprime_fci,\
                zprime_fci_sigma\
                = space_based_transform(stars_table,
                                        plot_results=plot_results,
                                        index=index,
                                        save_plots=save_plots)
            sb_final_transform_columns\
                = update_sb_final_transform_columns(sb_final_transform_columns,
                                                    index,
                                                    filter_fci,
                                                    filter_fci_sigma,
                                                    zprime_fci,
                                                    zprime_fci_sigma)
            # print(f"(V-clear) = ({filter_fci:.3f} +/- {filter_fci_sigma:.3f}) * {index} + " \
            #       f"({zprime_fci:.3f} +/- {zprime_fci_sigma:.3f})")
    sb_final_transform_table = create_sb_final_transform_table(
        sb_final_transform_columns)
    if save_plots:
        formats = {
            'T_fCI': '%0.3f',
            'T_fCI_sigma': '%0.3f',
            'Z_fCI': '%0.3f',
            'Z_fCI_sigma': '%0.3f'
        }
        write_table_to_latex(sb_final_transform_table,
                             f"{os.path.join(save_loc, f'{unique_id}_transform_table')}.txt",
                             formats=formats)
    return sb_final_transform_table


def _main_sc_lightcurve(directory,
                        gb_final_transforms=None,
                        temp_dir='tmp',
                        save_loc='Outputs',
                        file_suffix=(".fits", ".fit", ".fts"),
                        ecct_cut=0.5,
                        max_distance_from_sat=20,
                        size=25,
                        max_num_nan=5,
                        plot_results=0):
    filecount, filenames = copy_and_rename(directory=directory,
                                           file_suffix=file_suffix,
                                           temp_dir=temp_dir,
                                           debugging=True)
    # remove_temp_dir(temp_dir=temp_dir)
    # return
    set_sat_positions_bool = True
    change_sat_positions_bool = False
    num_nan = 0
    for filenum, file in enumerate(filenames):
        filepath = f"{temp_dir}/{file}"
        hdr, imgdata = read_fits_file(filepath)
        if set_sat_positions_bool:
            set_sat_positions_bool, sat_information = set_sat_positions(
                imgdata, filecount, set_sat_positions_bool)
        sat_information = add_new_time_and_filter(
            hdr, sat_information, filenum)
        if change_sat_positions_bool:
            change_sat_positions_bool,\
                sat_information\
                = change_sat_positions(filenames,
                                       filenum,
                                       num_nan,
                                       sat_information,
                                       change_sat_positions_bool,
                                       gb_final_transforms=gb_final_transforms,
                                       max_distance_from_sat=max_distance_from_sat,
                                       size=size,
                                       temp_dir=temp_dir,
                                       cmap_set='Set1',
                                       plot_results=plot_results)
        exptime = hdr['EXPTIME'] * u.s
        # bkg, bkg_std = calculate_img_bkg(imgdata)

        try:
            sat_x, sat_y, bkg_trm, fwhm = TRM_sat_detection(
                filepath, ecct_cut=ecct_cut)
        except TypeError:
            print("No satellites detected.")
            continue
        if fwhm < 0:
            continue
        # bkg, bkg_std = calculate_img_bkg(imgdata)
        # irafsources = detecting_stars(imgdata, bkg, bkg_std)
        # if not irafsources:
        #     sat_information.num_nans[:] = 0
        #     continue
        # plot_detected_sats(filenames[filenum],
        #                    plot_results,
        #                    imgdata,
        #                    irafsources,
        #                    sat_information,
        #                    max_distance_from_sat=max_distance_from_sat,
        #                    norm=LogNorm())
        # fwhms, fwhm, fwhm_std = calculate_fwhm(irafsources)
        bkg = np.median(bkg_trm)
        bsb = calculate_background_sky_brightness(bkg,
                                                  hdr,
                                                  exptime,
                                                  gb_final_transforms,
                                                  focal_length_key='FOCALLEN',
                                                  xpixsz_key='XPIXSZ',
                                                  ypixsz_key='YPIXSZ')
        fwhm_arcsec = convert_fwhm_to_arcsec_trm(hdr, fwhm)
        airmass = get_image_airmass(hdr)
        photometry_result = perform_photometry_sat(
            sat_x, sat_y, fwhm, imgdata, bkg_trm)
        instr_mags = calculate_magnitudes(photometry_result, exptime)
        instr_mags_sigma = calculate_magnitudes_sigma(
            photometry_result, exptime)
        sat_information = check_if_sat(sat_information,
                                       filenum,
                                       sat_x,
                                       sat_y,
                                       instr_mags,
                                       instr_mags_sigma,
                                       fwhm_arcsec,
                                       airmass,
                                       bsb,
                                       max_distance_from_sat=max_distance_from_sat)
        change_sat_positions_bool, num_nan\
            = determine_if_change_sat_positions(sat_information,
                                                filenum,
                                                change_sat_positions_bool,
                                                max_num_nan=max_num_nan)
        # del hdr
        # del imgdata
    # try:
    remove_temp_dir(temp_dir=temp_dir)
    # except PermissionError as e:
    #     print(e)
    if not os.path.exists(save_loc):
        os.mkdir(save_loc)
    sats_table = sat_information.sats_table
    ascii.write(
        sats_table,
        output=f"{save_loc}/Measured_Magnitudes.csv",
        format='csv')
    uncertainty_table = sat_information.uncertainty_table
    ascii.write(
        uncertainty_table,
        output=f"{save_loc}/Measured_Magnitude_Uncertainties.csv",
        format='csv')
    sat_auxiliary_table = sat_information.sat_auxiliary_table
    for aux_data in sat_auxiliary_table.columns[2:]:
        if all(np.isnan(sat_auxiliary_table[aux_data])):
            sat_auxiliary_table.remove_column(aux_data)
    ascii.write(sat_auxiliary_table,
                output=f"{save_loc}/Auxiliary_Information.csv", format='csv')
    unique_filters, num_filters, multiple_filters = determine_num_filters(
        sats_table)
    if multiple_filters:
        sat_dict = interpolate_sats(
            sats_table, uncertainty_table, unique_filters)
        if not gb_final_transforms:
            app_sat_dict = None
            save_interpolated_light_curve(sat_dict, save_loc)
            all_indices, all_indices_formatted\
                = get_all_indicies_combinations(unique_filters,
                                                num_filters,
                                                multiple_filters)
            colour_indices_dict = calculate_timeseries_colour_indices(
                sat_dict, all_indices)
            save_interpolated_light_curve(
                colour_indices_dict, save_loc, suffix="Colour Indices")
            filters_to_plot, indices_to_plot, aux_data_to_plot\
                = choose_indices_to_plot(unique_filters,
                                         num_filters,
                                         all_indices_formatted,
                                         sat_auxiliary_table)
            fig = axis_limits_multiband_gui(sat_dict,
                                            colour_indices_dict,
                                            sat_auxiliary_table,
                                            filters_to_plot,
                                            indices_to_plot,
                                            aux_data_to_plot,
                                            save_loc)
        else:
            app_sat_dict = apply_gb_timeseries_transforms(gb_final_transforms,
                                                          sat_dict,
                                                          sat_auxiliary_table,
                                                          unique_filters)
            save_interpolated_light_curve(app_sat_dict, save_loc)
            all_indices, all_indices_formatted\
                = get_all_indicies_combinations(unique_filters,
                                                num_filters,
                                                multiple_filters)
            colour_indices_dict = calculate_timeseries_colour_indices(
                app_sat_dict, all_indices)
            save_interpolated_light_curve(
                colour_indices_dict, save_loc, suffix="Colour Indices")
            filters_to_plot,\
                indices_to_plot,\
                aux_data_to_plot\
                = choose_indices_to_plot(unique_filters,
                                         num_filters,
                                         all_indices_formatted,
                                         sat_auxiliary_table)
            fig = axis_limits_multiband_gui(app_sat_dict,
                                            colour_indices_dict,
                                            sat_auxiliary_table,
                                            filters_to_plot,
                                            indices_to_plot,
                                            aux_data_to_plot,
                                            save_loc)
    else:
        sat_dict = None
        app_sat_dict = None
        aux_data_to_plot = choose_aux_data_to_plot(sat_auxiliary_table)
        axis_limits_singleband_gui(sats_table,
                                   uncertainty_table,
                                   sat_auxiliary_table,
                                   aux_data_to_plot,
                                   save_loc)
        # plot_light_curve_singleband(sats_table, uncertainty_table,
        # sat_auxiliary_table, aux_data_to_plot, save_loc)
    return sat_dict,\
        app_sat_dict,\
        sats_table,\
        uncertainty_table,\
        sat_auxiliary_table


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
    write_table_to_latex(
        gb_final_transforms,
        f"{os.path.join(save_loc, 'gb_final_transforms')}.txt",
        formats=formats)
    return


def BackgroundEstimationMulti(fitsdata, sigma_clip, bkgmethod, printval):
    sigma_clip = SigmaClip(sigma=2.5)
    bkg = SExtractorBackground(sigma_clip)
    # bkg_value1 = bkg.calc_background(fitsdata)

    # print(bkg_value)
    bkg = MeanBackground(sigma_clip)
    bkg_value2 = bkg.calc_background(fitsdata)

    bkg = MedianBackground(sigma_clip)
    bkg_value3 = bkg.calc_background(fitsdata)

    # bkg = ModeEstimatorBackground(sigma_clip)
    # bkg_value4 = bkg.calc_background(fitsdata)

    bkg_estimator2 = SExtractorBackground()
    # bkg = Background2D(fitsdata, (2, 2),
    # filter_size=(3,3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2)
    # Closest Approximate to Matlab Result
    bkg = Background2D(fitsdata, (50, 50), filter_size=(
        3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2)
    bg_rem = fitsdata - bkg.background

    if printval == 1:
        print("Background Solutions")
        print("Sigma Clip: " + str(sigma_clip))
        print("Using: " + bkgmethod)
        print("---------------------------")
        print("SExtractor Background: " + str(np.mean(bkg.background)))
        print("SExtractor Background(Filtered): " + str(np.mean(
            bkg.background)) + "\n " + "     " + "Box Size: " + "50x50" +
            "\n " + "     " + "Filter Size: " + "3x3")
        print("Mean Background: " + str(bkg_value2))
        print("Median Background: " + str(bkg_value3))
        # print("Mode Estimator Background: " + str(bkg_value4))
        print("Remaining Background (subtracted): " + str(bg_rem))
        print("Polyfit Background: Not Implemented Yet")

    else:
        return bkg


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


def ref_star_search(s,
                    f,
                    erad,
                    edec,
                    HIP,
                    vref,
                    bvindex,
                    vrindex,
                    refstarsfin):
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
        # starx = mstar.X
        # stary = mstar.Y
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
                            "Reference Mag: " + str(vref2[j]) + " vs " +
                            "Detected Mag: " + str(vmag))
                        print("")
                        Bvtransform = (vref2[j] - vmag) / bvindexdet[j]
                        print("B-V Transform: " + str(Bvtransform))
                        Vrtransform = (vref2[j] - vmag) / vrindexdet[j]
                        print("V-R Transform: " + str(Vrtransform))


def ref_read(refstars_doc):
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


def pinpoint_init():
    f = win32com.client.Dispatch("Pinpoint.plate")
    return f


def getFileList(inbox):
    filepathall = []
    # directory = os.path.dirname(inbox)
    list1 = os.listdir(inbox)  # List of Files
    listSize = len(list1)  # Number of Files
    print(listSize)
    c = list1
    print(c)
    # o=0;
    for i in range(1, listSize):
        print(c[i])
        filepath2 = inbox + "\\" + c[i]
        filepathall.append(filepath2)
        # o=o+1;
    return filepathall
    # o=0;


def fits_header_import(filepath, filter_key='FILTER'):
    imagehdularray = fits.open(filepath)
    header = imagehdularray[0].header
    date = imagehdularray[0].header['DATE-OBS']
    exposuretime = imagehdularray[0].header['EXPTIME']
    imagesizeX = imagehdularray[0].header['NAXIS1']
    imagesizeY = imagehdularray[0].header['NAXIS2']
    fitsdata = imagehdularray[0].data
    # focal_Length=imagehdularray[0].header['FOCALLEN']
    XPIXSZ = imagehdularray[0].header['XPIXSZ']
    YPIXSZ = imagehdularray[0].header['YPIXSZ']
    filt = imagehdularray[0].header['FILTER']
    wcs = WCS(header)
    return imagehdularray,\
        date,\
        exposuretime,\
        imagesizeX,\
        imagesizeY,\
        fitsdata,\
        filt,\
        header,\
        XPIXSZ,\
        YPIXSZ,\
        wcs


def calc_ArcsecPerPixel(header):
    focal_Length = header['FOCALLEN']
    xpix_size = header['XPIXSZ']
    ypix_size = header['XPIXSZ']
    # xbin = header['XPIXSZ']
    # ybin = header['XPIXSZ']
    x_arcsecperpixel = 206.2648 * ((xpix_size) / (focal_Length))
    y_arcsecperpixel = 206.2648 * ((ypix_size) / (focal_Length))
    # x_arcsecperpixel = math.atan(xpix_size/focal_Length)*3600*xbin
    # y_arcsecperpixel = math.atan(ypix_size/focal_Length)*3600*ybin

    return x_arcsecperpixel, y_arcsecperpixel


def edge_Protect(bg_rem, edge_protect, imagesizeX, imagesizeY, fitsdata):
    bg_rem[1:edge_protect, 1:edge_protect] = 0
    bg_rem[imagesizeX - edge_protect:imagesizeX, :] = 0
    bg_rem[:, 1:edge_protect] = 0
    bg_rem[:, imagesizeY - edge_protect:imagesizeY] = 0
    im_mean = np.mean(bg_rem)
    im_rms = np.std(fitsdata)

    return im_mean, bg_rem, im_rms

# TODO: change the SCInstrMagLightCurve.py code into functions in this file.
