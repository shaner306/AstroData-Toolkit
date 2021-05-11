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
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.table import Table, QTable
import astropy.units as u
from astropy.time import Time
from collections import namedtuple
from matplotlib import pyplot as plt
from photutils.detection import IRAFStarFinder
from photutils.psf import DAOGroup, BasicPSFPhotometry, IntegratedGaussianPRF
import numpy as np
import re


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
    instr_filter = hdr[filter_key].lower()
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
    # iraffind = IRAFStarFinder(threshold=bkg+3*bkg_std, fwhm=fwhm)
    iraffind = IRAFStarFinder(threshold=4*bkg_std, fwhm=fwhm)
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


def convert_ra_dec_to_alt_az(skypositions, hdr):
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
    lat = hdr['SITELAT']
    lon = hdr['SITELONG']
    height = hdr['SITEELEV']
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
        filter_name_repeat = np.full(num_stars, hdr['FILTER'])
        updated_large_table_columns.filters.extend(filter_name_repeat)
        updated_large_table_columns.V_apparents.extend(matched_stars.ref_star['V'])
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
        updated_large_table_columns.filters.append(hdr['FILTER'])
        updated_large_table_columns.V_apparents.append(matched_stars.ref_star['V'])
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
                'V',
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
                'V',
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


def space_based_transform(stars_table, 
                          plot_results=False, 
                          index='(B-V)', 
                          app_filter='V', 
                          instr_filter='clear', 
                          field=None):
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
    filter_fci, zprime_fci = np.polyfit(stars_table[index], stars_table[app_filter] - stars_table[instr_filter], 1, 
                                        full=False, w=1/err_sum)
    if plot_results:
        index_plot = np.arange(start=min(stars_table[index])-0.2, stop=max(stars_table[index])+0.2, step=0.1)
        plt.errorbar(stars_table[index], stars_table[app_filter] - stars_table[instr_filter], 
                     yerr=err_sum, fmt='o', capsize=2)
        plt.plot(index_plot, filter_fci * index_plot + zprime_fci)
        plt.ylabel(f"{app_filter}-{instr_filter}")
        plt.xlabel(f"{index}")
        if not field:
            plt.title(f"({app_filter}-{instr_filter}) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        else:
            plt.title(f"{field}: ({app_filter}-{instr_filter}) = {filter_fci:.3f} * {index} + {zprime_fci:.3f}")
        plt.show()
        plt.close()
    return filter_fci, zprime_fci


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
    # c_fci_simga = []
    zprime_f = []
    # zprime_f_simga = []
    instr_filter = []
    colour_index = []
    airmass = []
    gb_transform_table_columns = namedtuple('gb_transform_table_columns', 
                                            ['field',
                                             'c_fci',
                                             'zprime_f',
                                             'instr_filter',
                                             'colour_index',
                                             'airmass'])
    return gb_transform_table_columns(field, 
                                      c_fci, 
                                      zprime_f, 
                                      instr_filter, 
                                      colour_index, 
                                      airmass)


def update_gb_transform_table_columns(gb_transform_table_columns,
                                      field,
                                      c_fci,
                                      zprime_f,
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
    updated_gb_transform_table_columns.zprime_f.append(zprime_f)
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
            'Zprime_f',
            'filter',
            'CI',
            'X'
            ],
        data=[
            gb_transform_table_columns.field,
            gb_transform_table_columns.c_fci,
            gb_transform_table_columns.zprime_f,
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
        app_mag = np.array(ref_star['V'] + ref_star[colour_index])
        app_mag_sigma = np.nan_to_num(ref_star['e_V'], nan=max(ref_star['e_V'])) + \
            np.nan_to_num(ref_star[f'e_{colour_index}'], nan=max(ref_star[f'e_{colour_index}']))
        # app_mag_sigma = np.array(ref_star['e_V'] + ref_star[f'e_{colour_index}'])
    elif instr_filter == 'v' or instr_filter == 'g':
        colour_index = 'B-V'
        app_filter = 'V'
        app_mag = np.array(ref_star['V'])
        app_mag_sigma = np.nan_to_num(ref_star['e_V'], nan=max(ref_star['e_V']))
        # app_mag_sigma = np.array(ref_star['e_V'])
    elif instr_filter == 'r':
        colour_index = 'V-R'
        app_filter = 'R'
        app_mag = np.array(ref_star['V'] - ref_star[colour_index])
        app_mag_sigma = np.nan_to_num(ref_star['e_V'], nan=max(ref_star['e_V'])) + \
            np.nan_to_num(ref_star[f'e_{colour_index}'], nan=max(ref_star[f'e_{colour_index}']))
        # app_mag_sigma = np.array(ref_star['e_V'] + ref_star[f'e_{colour_index}'])
    elif instr_filter == 'i':
        colour_index = 'V-I'
        app_filter = 'I'
        app_mag = np.array(ref_star['V'] - ref_star[colour_index])
        app_mag_sigma = np.nan_to_num(ref_star['e_V'], nan=max(ref_star['e_V'])) + \
            np.nan_to_num(ref_star[f'e_{colour_index}'], nan=max(ref_star[f'e_{colour_index}']))
        # app_mag_sigma = np.array(ref_star['e_V'] + ref_star[f'e_{colour_index}'])
    else:
        colour_index = None
        app_filter = None
        app_mag = None
        app_mag_sigma = None
    return app_mag, app_mag_sigma, app_filter, colour_index


def ground_based_first_order_transforms(matched_stars, instr_filter, plot_results=False, field=None):
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
    app_mag, app_mag_sigma, app_filter, colour_index = get_app_mag_and_index(matched_stars.ref_star, instr_filter)
    max_instr_filter_sigma = max(matched_stars.img_instr_mag_sigma)
    err_sum = app_mag_sigma + np.nan_to_num(matched_stars.img_instr_mag_sigma, nan=max_instr_filter_sigma)
    err_sum = np.array(err_sum)
    err_sum[err_sum == 0] = max(err_sum)
    c_fci, zprime_fci = np.polyfit(matched_stars.ref_star[colour_index], app_mag - matched_stars.img_instr_mag, 1, 
                                        full=False, w=1/err_sum)
    if plot_results:
        index_plot = np.arange(start=min(matched_stars.ref_star[colour_index])-0.2, 
                               stop=max(matched_stars.ref_star[colour_index])+0.2, 
                               step=0.1)
        plt.errorbar(matched_stars.ref_star[colour_index], app_mag - matched_stars.img_instr_mag, 
                     yerr=err_sum, fmt='o', capsize=2)
        plt.plot(index_plot, c_fci * index_plot + zprime_fci)
        plt.ylabel(f"{app_filter}-{instr_filter}")
        plt.xlabel(f"{colour_index}")
        if not field:
            plt.title(f"({app_filter}-{instr_filter}) = {c_fci:.3f} * {colour_index} + {zprime_fci:.3f}")
        else:
            plt.title(f"{field}: ({app_filter}-{instr_filter}) = {c_fci:.3f} * {colour_index} + {zprime_fci:.3f}")
        plt.show()
        plt.close()
    return c_fci, zprime_fci


def ground_based_second_order_transforms(gb_transform_table, plot_results=False):
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
            T_fCI : float
                The instrumental transform coefficient for filter f using the colour index CI.
            k'_f : float
                The first order atmospheric extinction coefficient for filter f.
            Z_f : float
                The zero point magnitude for filter f.

    """
    gb_final_transforms = Table()
    unique_filters = table.unique(gb_transform_table, keys='filter')
    num_filters = len(unique_filters)
    nan_array = np.empty(num_filters)
    nan_array.fill(np.nan)
    gb_final_transforms = Table(
        names=[
            'filter',
            'CI',
            'k\'\'_fCI',
            'T_fCI',
            'k\'_f',
            'Z_f'
            ],
        data=[
            np.empty(num_filters, dtype=object),
            np.empty(num_filters, dtype=object),
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
        mask = gb_transform_table['filter'] == unique_filter
        current_filter = gb_transform_table[mask]
        kprimeprime_fci, t_fci = np.polyfit(current_filter['X'], current_filter['C_fCI'], 1)
        kprime_f, zprime_f = np.polyfit(current_filter['X'], current_filter['Zprime_f'], 1)
        gb_final_transforms['k\'\'_fCI'][unique_filter_index] = kprimeprime_fci
        gb_final_transforms['T_fCI'][unique_filter_index] = t_fci
        gb_final_transforms['k\'_f'][unique_filter_index] = kprime_f
        gb_final_transforms['Z_f'][unique_filter_index] = zprime_f
        if plot_results:
            X_plot = np.arange(start=min(current_filter['X'])-0.2, stop=max(current_filter['X'])+0.2, step=0.1)
            ci_plot = re.sub('[^a-zA-Z]+', '', current_index)
            ci_plot = ci_plot.lower()
            plt.plot(current_filter['X'], current_filter['C_fCI'], 'o')
            plt.plot(X_plot, kprimeprime_fci*X_plot+t_fci)
            plt.title(f'C_{unique_filter}{ci_plot} = {kprimeprime_fci:.3f} * X + {t_fci:.3f}')
            plt.ylabel(f'C_{unique_filter}{ci_plot}')
            plt.xlabel('X')
            plt.show()
            plt.close()
            plt.plot(current_filter['X'], current_filter['Zprime_f'], 'o')
            plt.plot(X_plot, kprime_f*X_plot+zprime_f)
            plt.title(f'Z\'_{unique_filter} = {kprime_f:.3f} * X + {zprime_f:.3f}')
            plt.ylabel(f'Z\'_{unique_filter}')
            plt.xlabel('X')
            plt.show()
            plt.close()
    
    return gb_final_transforms


# TODO
# This will be harder than I thought. Namely, keeping track of which star was which when it's not known (unless maybe 
# I do assume to know it because it is just a check? Maybe by using matched_stars or large_stars_table instead of 
# creating a whole new table?). Either way, I'll start to convert the light curve code first and then come back to this 
# later.


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


def apply_gb_transforms(gb_final_transforms, unknown_object_table):
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
    app_mag_table = Table(names=['time', 'filter', 'apparent mag'])
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


# TODO: change the SCInstrMagLightCurve.py code into functions in this file.