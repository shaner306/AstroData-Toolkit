# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 11:00:58 2021

@author: shane

Credit to NEOSSAT Dark Subtraction Algorithm. - https://github.com/jasonfrowe/neossat/tree/master/neossat

"""
from tkinter.filedialog import SaveFileDialog
import numpy as np
from photutils import CircularAperture
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import re
import sys
from os.path import dirname
src_path = dirname(dirname(__file__))
sys.path.append(os.path.join(src_path, 'general_tools'))
import AstroFunctions as astro


def plot_histogram(scidata, imstat, sigscalel, sigscaleh):
    """
    Plot the histogram of Image Counts (ADU) VS Frequency.
    Parameters
    ----------
    scidata : numpy.ndarray
        The image to be plotted.
    imstat : list
        The statistics of the image.
    sigscalel : float
        The lower limit of the sigma scale.
    sigscaleh : float
        The upper limit of the sigma scale.

    Returns
    -------
    None
    """

    matplotlib.rcParams.update({'font.size': 24})  # Adjust font.

    flat = scidata.flatten()
    vmin = np.min(flat[flat > imstat[2] - imstat[3]*sigscalel])
    vmax = np.max(flat[flat < imstat[2] + imstat[3]*sigscaleh])

    plt.figure(figsize=(12, 6))  # Adjust size of figure.
    image_hist = plt.hist(scidata.flatten(), 100, range=(vmin, vmax))
    plt.xlabel('Image Counts (ADU)')
    plt.ylabel('Number Count')
    plt.show()

    return


def plot_image_wsource(scidata, imstat, sigscalel, sigscaleh, sources=None, xy=None, figname=None, display=True):
    """
    Plot the image with the sources highlighted.

    **Part of NEOSSAT Dark Subtraction Algorithm.** - https://github.com/jasonfrowe/neossat/tree/master/neossat

    Parameters
    ----------
    scidata : numpy.ndarray
        The image to be plotted.
    imstat : list
        The statistics of the image.
    sigscalel : float
        The lower limit of the sigma scale.
    sigscaleh : float
        The upper limit of the sigma scale.
    sources : list
        The sources to be highlighted.
    xy : list
        The x and y coordinates of the sources.
    figname : str
        The name of the figure.
    display : bool
        Whether to display the figure.

    Returns
    -------
    None
    """

    eps = 1.0e-9
    sigscalel = -np.abs(sigscalel)  # Expected to be negative.
    sigscaleh = np.abs(sigscaleh)  # Expected to be positive.

    matplotlib.rcParams.update({'font.size': 24})  # Adjust font.

    flat = scidata.flatten()
    vmin = np.min(flat[flat > imstat[2] + imstat[3]*sigscalel]) - imstat[0] + eps
    vmax = np.max(flat[flat < imstat[2] + imstat[3]*sigscaleh]) - imstat[0] + eps

    if sources is not None:
        positions = np.column_stack([sources['xcentroid'], sources['ycentroid']])
    elif xy is not None:
        positions = np.column_stack(xy)
    else:
        raise ValueError('Either sources or xy must be give.')

    apertures = CircularAperture(positions, r=4.)

    plt.figure(figsize=(20, 20))  # Adjust size of figure.
    imgplot = plt.imshow(scidata - imstat[0], norm=LogNorm(), vmin=vmin, vmax=vmax)  # TODO scidata index?
    apertures.plot(color='red', lw=1.5, alpha=0.5)
    for i in range(len(positions)):
        plt.annotate('{}'.format(i), positions[i])
    plt.axis((-0.5, scidata.shape[1]-0.5, -0.5, scidata.shape[0]-0.5))
    plt.xlabel("Column (Pixels)")
    plt.ylabel("Row (Pixels)")

    if figname is not None:
        plt.savefig(figname)

    if display:
        plt.show()

    plt.close()

    return


def plot_image(scidata, imstat, sigscalel, sigscaleh):
    """"""

    eps = 1.0e-9
    sigscalel = -np.abs(sigscalel)  # Expected to be negative.
    sigscaleh = np.abs(sigscaleh)  # Expected to be positive.

    matplotlib.rcParams.update({'font.size': 24})  # Adjust font.

    flat = scidata.flatten()
    vmin = np.min(flat[flat > imstat[2] + imstat[3] * sigscalel]) - imstat[0] + eps
    vmax = np.max(flat[flat < imstat[2] + imstat[3] * sigscaleh]) - imstat[0] + eps

    plt.figure(figsize=(20, 20))  # Adjust size of figure.
    imgplot = plt.imshow(scidata[:, :] - imstat[0], norm=LogNorm(), vmin=vmin, vmax=vmax)
    plt.axis((0, scidata.shape[1], 0, scidata.shape[0]))
    plt.xlabel("Column (Pixels)")
    plt.ylabel("Row (Pixels)")
    plt.show()

    return


def plot_match_confirmation(wcs, imgdata, matched_stars, reference_stars, unique_id, save_loc, save_plots=False, name_key='Name'):
    viridis = cm.get_cmap('viridis', 2)
    colours = viridis(np.linspace(0, 1, 2))
    if save_plots:
        save_loc = os.path.join(save_loc, 'Annotated Images')
        if not os.path.exists(save_loc):
            os.mkdir(save_loc)
    img_star_x, img_star_y = wcs.world_to_pixel(matched_stars.img_star_loc)
    field_in_img = astro.get_field_name(matched_stars, name_key=name_key)
    ref_names = np.array(reference_stars[name_key])
    field_names = np.empty(np.shape(ref_names), dtype=object)
    for i, name in enumerate(ref_names):
            split_string = re.split('[^a-zA-Z0-9]', str(name))
            if len(split_string) > 1:
                field_names[i] = ' '.join(split_string[:-1])
            elif len(split_string) == 1:
                field_names[i] = split_string[0]
    mask = field_names == field_in_img
    field_stars = reference_stars[mask]
    field_star_loc = SkyCoord(ra=field_stars['RA'], dec=field_stars['Dec'], 
                                  unit=(u.hourangle, u.deg))
    ref_star_x, ref_star_y = wcs.world_to_pixel(field_star_loc)
    mask_y = ((ref_star_y >= 0) & (ref_star_y <= np.shape(imgdata)[0]))
    mask_x = ((ref_star_x >= 0) & (ref_star_x <= np.shape(imgdata)[1]))
    mask_outside_fov = (mask_x & mask_y)
    ref_star_x = ref_star_x[mask_outside_fov]
    ref_star_y = ref_star_y[mask_outside_fov]
    field_stars = field_stars[mask_outside_fov]
    try:
        num_img_stars = len(matched_stars.img_instr_mag)
    except TypeError:
        return
    num_field_stars = len(field_stars)
    fig = plt.figure(figsize=(12,8))
    ax = plt.subplot(projection=wcs)
    ax.imshow(imgdata, cmap='gray', norm=LogNorm(), interpolation='nearest')
    ax.scatter(ref_star_x, ref_star_y,
                s=100, edgecolor=colours[0], facecolor='none', label=f'Reference Star from File ({num_field_stars})')

    ax.scatter(img_star_x, img_star_y,
                s=100, edgecolor=colours[1], facecolor='none', label=f'Reference Star from Image ({num_img_stars})')
    ax.grid(color='gray', ls='solid')

    ref_star_names = np.array(field_stars[name_key])
    app_mags = np.array(field_stars['V_ref'])
    i=0
    for x,y in  zip(ref_star_x,ref_star_y):
        ref_star_name = ref_star_names[i]
        app_mag = app_mags[i]
        i = i + 1
        plt.annotate(f"{ref_star_name} ({app_mag})", (x, y), textcoords="offset points", xytext=(0, 10), ha='center')

    # HIP_title = reference_stars.colnames[0]
    # ref_name = reference_stars[HIP_title][possible_ref_star_index]
    ax.set_ylabel('Declination (J2000)')
    ax.set_xlabel('Right Ascension (J2000)')
    plt.title(f"{field_in_img}: {unique_id}")
    plt.legend()
    if save_plots:
        plt.savefig(os.path.join(save_loc, f"{unique_id}.png"))
    plt.show()
    plt.close()  
    return num_img_stars, num_field_stars