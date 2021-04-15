# -*- coding: utf-8 -*-
"""
AstroFunctions.py.

This file holds all of the functions that will likely be implemented in a final version of the image processor.

Created on Thu Apr 15 10:14:43 2021

@author: Jack Wawrow
"""
from astropy.io import fits, ascii
import numpy as np


def read_ref_stars(ref_stars_file):
    """
    Read a file containing information regarding the reference stars to be used to calculate the transforms.
    
    Parameters
    ----------
        ref_stars_file : string
            Location of the reference stars file.
            
    Returns
    -------
        reference_stars : <astropy.table> object
            Table with the data extracted from ref_stars_file.
    """
    try:
        reference_stars = ascii.read(ref_stars_file, format='basic', delimiter='\t', guess=False, encoding='UTF-8')
    except Exception:
        reference_stars = ascii.read(ref_stars_file, encoding='UTF-8')
    reference_stars = reference_stars.filled(np.nan)
    return reference_stars


def read_fits_file(filepath):
    """
    Read a fits file and return its header and data.

    Parameters
    ----------
    filepath : string
        Location of the fits file.

    Returns
    -------
    hdr : <astropy.io.fits.header.Header> object
        Header from the fits file.
    imgdata : numpy.ndarray
        Data from the fits file.

    """
    with fits.open(filepath) as hdul:
        hdr = hdul[0].header
        imgdata = hdul[0].data
    return hdr, imgdata

