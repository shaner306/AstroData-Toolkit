# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 11:00:58 2021

@author: shane
"""

from tkinter.filedialog import SaveFileDialog
import numpy as np

from photutils import CircularAperture

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os


def plot_histogram(scidata, imstat, sigscalel, sigscaleh):
    """"""

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
    """"""

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


def plot_match_confirmation(wcs, imgdata, matched_stars, unique_id, save_loc, save_plots=False, name_key='Name'):
    if save_plots:
        save_loc = os.path.join(save_loc, 'Annotated Images')
        if not os.path.exists(save_loc):
            os.mkdir(save_loc)
    ref_star_x, ref_star_y = wcs.world_to_pixel(matched_stars.ref_star_loc)
    img_star_x, img_star_y = wcs.world_to_pixel(matched_stars.img_star_loc)
    fig = plt.figure(figsize=(12,8))
    ax = plt.subplot(projection=wcs)
    ax.imshow(imgdata, cmap='gray', norm=LogNorm(), interpolation='nearest')
    ax.scatter(ref_star_x, ref_star_y,
                s=100, edgecolor='red', facecolor='none', label='Reference Star from File')

    ax.scatter(img_star_x, img_star_y,
                s=100, edgecolor='green', facecolor='none', label='Reference Star from Image')
    ax.grid(color='gray', ls='solid')

    ref_star_names = np.array(matched_stars.ref_star[name_key])
    app_mags = np.array(matched_stars.ref_star['V_ref'])
    i=0
    try:
        for x,y in  zip(ref_star_x,ref_star_y):

            ref_star_name = ref_star_names[i]
            app_mag = app_mags[i]
            i = i + 1
            plt.annotate(f"{ref_star_name} ({app_mag})", (x, y), textcoords="offset points", xytext=(0, 10), ha='center')
    except TypeError:
        return
    # HIP_title = reference_stars.colnames[0]
    # ref_name = reference_stars[HIP_title][possible_ref_star_index]
    ax.set_ylabel('Declination (J2000)')
    ax.set_xlabel('Right Ascension (J2000)')
    plt.title(unique_id)
    plt.legend()
    if save_plots:
        plt.savefig(os.path.join(save_loc, f"{unique_id}.png"))
    plt.show()
    plt.close()  
    return