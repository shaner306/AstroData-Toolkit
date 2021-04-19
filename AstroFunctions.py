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
from astropy.wcs import WCS
from collections import namedtuple
from photutils.detection import IRAFStarFinder
from photutils.psf import DAOGroup, BasicPSFPhotometry, IntegratedGaussianPRF
import numpy as np
import os


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
                   ground_based=False, 
                   altazpositions=None, 
                   max_ref_sep=10.0):
    """
    Match the stars detected in the image to those provided in the reference star file.

    Parameters
    ----------
    **************************************************************************
    ref_stars_file : string
        Location of the reference stars file.
    **************************************************************************
    skypositions : astropy.coordinates.sky_coordinate.SkyCoord
        AstroPy SkyCoord object containing the RA/dec positions of all sources in the image.
    instr_mags : numpy array
        Array containing the instrumental magnitudes of the sources in the image.
    instr_mags_sigma : numpy array
        Array containing the standard deviation of the instrumental magnitudes of the sources in the image.
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
                Angular distance between the matched star(s) detected in the image and ref_stars_file
            img_instr_mag : numpy array
                Array containing the instrumental magnitudes of the matched star(s) detected in the image.
            img_instr_mag_sigma : numpy array
                Array containing the standard deviations of the instrumental magnitudes of the matched star(s) 
                detected in the image.
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
        img_star_altaz, 
        img_star_airmass)


def init_large_table_columns():
    ref_star_name = []
    times = []
    flux_table = []
    exposure = []
    ref_star_RA = []
    ref_star_dec = []
    img_star_RA = []
    img_star_dec = []
    angular_separation = []
    ref_star_x = []
    ref_star_y = []
    img_star_x = []
    img_star_y = []
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
    X_rounded = []
    large_table_columns = namedtuple('large_table_columns',
                                     ['ref_star_name',
                                      'times',
                                      'flux_table',
                                      'exposure',
                                      'ref_star_RA',
                                      'ref_star_dec',
                                      'img_star_RA',
                                      'img_star_dec',
                                      'angular_separation',
                                      'ref_star_x',
                                      'ref_star_y',
                                      'img_star_x',
                                      'img_star_y',
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
                                      'X_rounded'])
    return large_table_columns(ref_star_name, 
                               times, 
                               flux_table, 
                               exposure, 
                               ref_star_RA, 
                               ref_star_dec, 
                               img_star_RA, 
                               img_star_dec, 
                               angular_separation, 
                               ref_star_x, 
                               ref_star_y, 
                               img_star_x, 
                               img_star_y, 
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
                               X_rounded)


def update_large_table_columns(large_table_columns, matched_stars, hdr, exptime, name_key='Name'):
    updated_large_table_columns = large_table_columns
    try:
        num_stars = len(matched_stars.img_instr_mag)
    except TypeError:
        num_stars = 1
    if num_stars > 1:
        updated_large_table_columns.ref_star_name.extend(matched_stars.ref_star[name_key])
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
        # X and Y to be updated
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
        if not matched_stars.img_star_airmass:
            return updated_large_table_columns
        updated_large_table_columns.img_star_airmass.extend(matched_stars.img_star_airmass)
        updated_large_table_columns.X_rounded.extend(round(matched_stars.img_star_airmass, 1))
    elif num_stars == 1:
        updated_large_table_columns.ref_star_name.append(matched_stars.ref_star[name_key])
        time = Time(hdr['DATE-OBS'], format='fits')
        updated_large_table_columns.times.append(time.jd)
        updated_large_table_columns.exposure.append(exptime)
        updated_large_table_columns.ref_star_RA.append(matched_stars.ref_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.ref_star_dec.append(matched_stars.ref_star_loc.dec)
        updated_large_table_columns.img_star_RA.append(matched_stars.img_star_loc.ra.to(u.hourangle))
        updated_large_table_columns.img_star_dec.append(matched_stars.img_star_loc.dec)
        updated_large_table_columns.angular_separation.append(matched_stars.ang_separation.to(u.arcsec))
        # X and Y to be updated
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
        if not matched_stars.img_star_airmass:
            return updated_large_table_columns
        updated_large_table_columns.img_star_airmass.append(matched_stars.img_star_airmass)
        updated_large_table_columns.X_rounded.append(round(matched_stars.img_star_airmass, 1))
    else:
        return
    
    # num_stars = len(matched_stars.img_instr_mag)
    # large_table_columns.ref_star_name.extend(list(matched_stars.ref_star[name_key]))
    # time = Time(hdr['DATE-OBS'], format='fits')
    # time_repeat = np.full(num_stars, time.jd)
    # large_table_columns.times.extend(list(time_repeat))
    # exposure_repeat = np.full(num_stars, exptime)
    # large_table_columns.exposure.extend(list(exposure_repeat))
    # large_table_columns.ref_star_RA.extend(list(matched_stars.ref_star_loc.ra.to(u.hourangle)))
    # large_table_columns.ref_star_dec.extend(list(matched_stars.ref_star_loc.dec))
    # large_table_columns.img_star_RA.extend(list(matched_stars.img_star_loc.ra.to(u.hourangle)))
    # large_table_columns.img_star_dec.extend(list(matched_stars.img_star_loc.dec))
    # large_table_columns.angular_separation.extend(list(matched_stars.ang_separation.to(u.arcsec)))
    # # X and Y to be updated
    # large_table_columns.img_star_mag.extend(list(matched_stars.img_instr_mag))
    # large_table_columns.img_star_mag_sigma.extend(list(matched_stars.img_instr_mag_sigma))
    # filter_name_repeat = np.full(num_stars, hdr['FILTER'])
    # large_table_columns.filters.extend(list(filter_name_repeat))
    # large_table_columns.V_apparents.extend(list(matched_stars.ref_star['V']))
    # try:
    #     large_table_columns.B_V_apparents.extend(list(matched_stars.ref_star['(B-V)']))
    #     large_table_columns.U_B_apparents.extend(list(matched_stars.ref_star['(U-B)']))
    #     large_table_columns.V_R_apparents.extend(list(matched_stars.ref_star['(V-R)']))
    #     large_table_columns.V_I_apparents.extend(list(matched_stars.ref_star['(V-I)']))
    #     large_table_columns.V_sigma_apparents.extend(list(matched_stars.ref_star['V_sigma']))
    # except KeyError:
    #     large_table_columns.B_V_apparents.extend(list(matched_stars.ref_star['B-V']))
    #     large_table_columns.U_B_apparents.extend(list(matched_stars.ref_star['U-B']))
    #     large_table_columns.V_R_apparents.extend(list(matched_stars.ref_star['V-R']))
    #     large_table_columns.V_I_apparents.extend(list(matched_stars.ref_star['V-I']))
    #     large_table_columns.V_sigma_apparents.extend(list(matched_stars.ref_star['e_V']))
    # if not matched_stars.img_star_airmass:
    #     pass
    #     # return updated_large_table_columns
    # else:
    #     large_table_columns.img_star_airmass.extend(list(matched_stars.img_star_airmass))
    #     large_table_columns.X_rounded.extend(list(round(matched_stars.img_star_airmass.value, 1)))
    return updated_large_table_columns


def create_large_stars_table(large_table_columns, ground_based=False):
    if ground_based:
        large_stars_table = QTable(
            data=[
                large_table_columns.ref_star_name, 
                large_table_columns.times, 
                large_table_columns.ref_star_RA, 
                large_table_columns.ref_star_dec, 
                large_table_columns.img_star_RA, 
                large_table_columns.img_star_dec, 
                large_table_columns.angular_separation, 
                # large_table_columns.ref_star_x, 
                # large_table_columns.ref_star_y, 
                # large_table_columns.img_star_x, 
                # large_table_columns.img_star_y,  
                # large_table_columns.flux_table, 
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
                large_table_columns.X_rounded
                ],
            names=[
                'Name',
                'Time (JD)',
                'RA_ref',
                'dec_ref',
                'RA_img',
                'dec_img',
                'angular_separation',
                # 'x_ref',
                # 'y_ref',
                # 'x_img',
                # 'y_img',
                # 'flux',
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
                'X_rounded'
                ]
            )
    else:
        large_stars_table = QTable(
            data=[
                large_table_columns.ref_star_name, 
                large_table_columns.times, 
                large_table_columns.ref_star_RA, 
                large_table_columns.ref_star_dec, 
                large_table_columns.img_star_RA, 
                large_table_columns.img_star_dec, 
                large_table_columns.angular_separation, 
                # large_table_columns.ref_star_x, 
                # large_table_columns.ref_star_y, 
                # large_table_columns.img_star_x, 
                # large_table_columns.img_star_y,  
                # large_table_columns.flux_table, 
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
                'Name',
                'Time (JD)',
                'RA_ref',
                'dec_ref',
                'RA_img',
                'dec_img',
                'angular_separation',
                # 'x_ref',
                # 'y_ref',
                # 'x_img',
                # 'y_img',
                # 'flux',
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
    unique_stars = table.unique(large_stars_table, keys=keys)
    N = len(unique_stars)
    nan_array = np.empty(N)
    nan_array.fill(np.nan)
    apparent_mags_table = Table(
        names=[
            'Name',
            'V',
            '(B-V)',
            '(U-B)',
            '(V-R)',
            '(V-I)',
            ],
        data=[
            np.empty(N, dtype=object),
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            ]
        )
    different_filters = table.unique(large_stars_table, keys='filter')
    different_filter_list = list(different_filters['filter'])
    different_filter_list = different_filter_list.lower()
    different_filter_data = np.empty((N, len(different_filter_list)))
    different_filter_data.fill(np.nan)
    filter_sigma_list = []
    for different_filter in different_filter_list:
        filter_sigma_list.append(f"{different_filter}_sigma")
    different_filter_table = Table(data=different_filter_data, 
                                   names=different_filter_list)
    different_filter_sigma_table = Table(data=different_filter_data, 
                                         names=filter_sigma_list)
    if ground_based:
        filter_X_list = []
        for different_filter in different_filter_list:
            filter_X_list.append(f"X_{different_filter}")
        filter_X_sigma_list = []
        for different_filter in different_filter_list:
            filter_X_sigma_list.append(f"X_{different_filter}_sigma")
        filter_X_table = Table(data=different_filter_data, 
                               names=filter_X_list)
        filter_X_sigma_table = Table(data=different_filter_data, 
                                     names=filter_X_sigma_list)
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
    return stars_table

ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\FITS Tutorial\Reference_stars.csv'
# filepath = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-03-20 - Calibrated\Solved Images\HIP 2894\LIGHT\B\0001_3x3_-10.00_5.00_B_21-20-52.fits'
directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-03-20 - Calibrated\Solved Images'
ground_based = True

# ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Landolt Stars\2009_Landolt_Standard_Stars.txt'
# filepath = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Landolt Stars\SA108\NEOS_SCI_2020121233640_clean.fits'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Landolt Stars'
# ground_based = False

reference_stars, ref_star_positions = read_ref_stars(ref_stars_file)
large_table_columns = init_large_table_columns()

for dirpath, dirnames, filenames in os.walk(directory):
    for filename in filenames:
        if filename.endswith(".fits"):
        # if filename.endswith("_clean.fits"):
            filepath = os.path.join(dirpath, filename)
            hdr, imgdata = read_fits_file(filepath)
            exptime = hdr['EXPTIME']
            # exptime = hdr['AEXPTIME']
            bkg, bkg_std = calculate_img_bkg(imgdata)
            irafsources = detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
            if not irafsources:
                continue
            fwhm, fwhm_std = calculate_fwhm(irafsources)
            photometry_result = perform_photometry(irafsources, fwhm, imgdata, bkg=bkg)
            instr_mags = calculate_magnitudes(photometry_result, exptime)
            instr_mags_sigma = calculate_magnitudes_sigma(photometry_result, exptime)
            wcs = WCS(hdr)
            skypositions = convert_pixel_to_ra_dec(irafsources, wcs)
            # print(skypositions)
            altazpositions = None
            if ground_based:
                altazpositions = convert_ra_dec_to_alt_az(skypositions, hdr)
                # print(altazpositions)
            matched_stars = find_ref_stars(reference_stars, 
                                           ref_star_positions,
                                           skypositions,
                                           instr_mags,
                                           instr_mags_sigma,
                                           ground_based=ground_based,
                                           altazpositions=altazpositions)
            if not matched_stars:
                continue
            # photometry_result.pprint_all()
            # print(instr_mags)
            # print(instr_mags_sigma)
            # print(matched_stars)
            
            large_table_columns = update_large_table_columns(large_table_columns, matched_stars, hdr, exptime, name_key='HIP')

large_stars_table = create_large_stars_table(large_table_columns, ground_based=ground_based)
large_stars_table.pprint_all()
stars_table = group_each_star(large_stars_table, 
                              ground_based=ground_based, 
                              keys=['Name', 'X_rounded'])
stars_table.pprint_all()