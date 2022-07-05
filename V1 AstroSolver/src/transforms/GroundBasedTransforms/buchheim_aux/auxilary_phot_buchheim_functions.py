# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:31:30 2022

@author: mstew
"""
import perform_photometry
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
from astropy import table
from astropy.coordinates import EarthLocation, AltAz, SkyCoord,\
    match_coordinates_sky
from astropy.io import fits, ascii
from astropy.modeling.fitting import LevMarLSQFitter, FittingWithOutlierRemoval, LinearLSQFitter
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
from photutils.aperture import RectangularAperture,CircularAperture,CircularAnnulus
from photutils.background import Background2D
from photutils.background import MeanBackground
from photutils.background import MedianBackground
from photutils.background import SExtractorBackground
from photutils.detection import IRAFStarFinder
from photutils.psf import DAOGroup, BasicPSFPhotometry, IntegratedGaussianPRF
from skimage import measure
from tqdm import tqdm
from random import shuffle, choice
from astropy.nddata import CCDData
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import SigmaClip
from photutils.aperture import aperture_photometry,ApertureStats
from photutils.utils import calc_total_error


import perform_photometry
import auxilary_phot_boyde_functions as boyde_aux


import auxilary_phot_warner_functions as warn_aux
import AstroFunctions as astro


def calculate_slopes_Buchheim(stars_table,
                              different_filter_list,
                              save_plots,
                              **kwargs):
    stars_for_second_order_extinction,\
        multiple_stars = astro.get_stars_with_multiple_observations(
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
    if num_filters == 1:
        multiple_filters = False
    else:
        multiple_filters = True
    all_indices, all_indices_formatted = astro.get_all_indicies_combinations(
        different_filter_list, num_filters, multiple_filters)
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
            fit, or_fit, line_init = astro.init_linear_fitting(sigma=2.5)
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
                fit, or_fit, line_init = astro.init_linear_fitting(
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
        colour_indices = astro.get_all_colour_indices(different_filter)
        colour_index, _ = astro.get_colour_index_lower(different_filter)
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
        fit, or_fit, line_init = astro.init_linear_fitting(
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
        colour_indices_list = astro.get_all_colour_indices(instr_filter)
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
            colour_indices = astro.get_all_colour_indices(different_filter)
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

