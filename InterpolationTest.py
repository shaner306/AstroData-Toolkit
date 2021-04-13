from astropy.io import fits, ascii
import astropy.units as u
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.table import Table, QTable, hstack, unique, vstack
from astropy.time import Time
import datetime
from photutils.detection import IRAFStarFinder
from photutils.psf import DAOGroup, BasicPSFPhotometry, IntegratedGaussianPRF
import numpy as np
import os
from math import sqrt, atan
import matplotlib
from matplotlib import patches
from matplotlib.colors import LogNorm
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import ctypes
from tkinter import *
from matplotlib.lines import Line2D
from shutil import copy2, rmtree
import cv2 as cv
from photutils.aperture import RectangularAperture
matplotlib.use('TkAgg')

# Just used to setup the debugging code.
sat_names = ['Intelsat 10-02', 'MEV 2', 'THOR 6', 'THOR 5']

b_sats_table = ascii.read('C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/b_instr_mags.csv')
g_sats_table = ascii.read('C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/g_instr_mags.csv')
r_sats_table = ascii.read('C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/r_instr_mags.csv')
b_uncertainty_table = ascii.read('C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/b_uncertainty.csv')
g_uncertainty_table = ascii.read('C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/g_uncertainty.csv')
r_uncertainty_table = ascii.read('C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/r_uncertainty.csv')
b_fwhm_table = ascii.read('C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/b_fwhm.csv')
g_fwhm_table = ascii.read('C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/g_fwhm.csv')
r_fwhm_table = ascii.read('C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/r_fwhm.csv')
large_sats_table = vstack([b_sats_table, g_sats_table, r_sats_table])
large_sats_table.sort('Time (JD)')

sat_fwhm_table = vstack([b_fwhm_table, g_fwhm_table, r_fwhm_table])
sat_fwhm_table.sort('Time (JD)')
avg_sat_fwhm = []
avg_sat_fwhm_std = []

arcsec_per_pix = 1.5599859295701453

for row in sat_fwhm_table:
    sat_full_row_numpy = np.array(list(row))
    sat_row_numpy = sat_full_row_numpy[2:].astype(float)
    avg_sat_fwhm.append(np.nanmean(sat_row_numpy) * arcsec_per_pix)
    avg_sat_fwhm_std.append(np.nanstd(sat_row_numpy) * arcsec_per_pix)
# avg_sat_fwhm_arcsec = avg_sat_fwhm * arcsec_per_pix
# New code. Should be able to be transplanted directly into SCInstrMagLightCurve.py
times_list = np.array(large_sats_table['Time (JD)'])
times_obj = Time(times_list, format='jd', scale='utc')
times_datetime = times_obj.to_value('datetime')

fig, ax = plt.subplots()
plt.errorbar(times_datetime, avg_sat_fwhm, yerr=avg_sat_fwhm_std, fmt='o', capsize=2)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.title('Seeing - 20 Mar')
plt.ylabel('FWHM (arcsec)')
plt.xlabel('Time (UTC)')
plt.show(block=True)
# plt.pause(3)
plt.close()

for sat in sat_names:
    b_interpolated = np.interp(times_list, b_sats_table['Time (JD)'][~np.isnan(b_sats_table[sat])],
                               b_sats_table[sat][~np.isnan(b_sats_table[sat])])
    g_interpolated = np.interp(times_list, g_sats_table['Time (JD)'][~np.isnan(g_sats_table[sat])],
                               g_sats_table[sat][~np.isnan(g_sats_table[sat])])
    r_interpolated = np.interp(times_list, r_sats_table['Time (JD)'][~np.isnan(r_sats_table[sat])],
                               r_sats_table[sat][~np.isnan(r_sats_table[sat])])
    b_interpolated[np.isnan(large_sats_table[sat])] = np.nan
    g_interpolated[np.isnan(large_sats_table[sat])] = np.nan
    r_interpolated[np.isnan(large_sats_table[sat])] = np.nan
    g_regular = np.full(len(g_interpolated), np.nan)
    g_uncertainty = np.full(len(g_interpolated), np.nan)
    mask = np.in1d(large_sats_table['Time (JD)'], g_sats_table['Time (JD)'])
    # print(mask)
    g_regular[mask] = g_sats_table[sat]
    g_uncertainty[mask] = g_uncertainty_table[sat]

    fig, ax = plt.subplots()
    ax.plot(times_datetime, g_interpolated, 'ko', markersize=3)
    ax.errorbar(times_datetime, g_regular, yerr=g_uncertainty, fmt='ko', markersize=3, capsize=1, label='g')
    ax.set_ylabel('Instrumental Magnitude')
    # ax.set_ylim([min(g_interpolated)*1.05, max(g_interpolated)*0.5])
    ax.set_ylim(-12.2, -3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.gca().invert_yaxis()
    ax2 = ax.twinx()
    ax2.plot(times_datetime, b_interpolated - g_interpolated, 'bo', label='b-g', markersize=3)
    ax2.plot(times_datetime, b_interpolated - r_interpolated, 'go', label='b-r', markersize=3)
    ax2.plot(times_datetime, g_interpolated - r_interpolated, 'ro', label='g-r', markersize=3)
    ax2.set_ylabel('Colour Index')
    ax2.set_ylim([-5, 2])
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.gca().invert_yaxis()
    ax.set_xlabel("Time (UTC)")
    fig.legend()
    plt.title(f'{sat} Light Curve - 20 Mar')
    plt.show(block=True)
    # plt.pause(5)
    plt.close()


