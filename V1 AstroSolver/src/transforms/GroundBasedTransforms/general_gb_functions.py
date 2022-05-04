import os
from collections import namedtuple
from itertools import groupby
from math import sqrt

import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt

from AstroFunctions import Delta, get_all_colour_indices, get_app_mag_and_index_AVG, init_linear_fitting


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