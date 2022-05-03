# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:31:56 2022

@author: mstew

This file contains specific auxillary funcitons which correspond to the warner transformations
"""
import os
from math import sqrt

import matplotlib.cm as cm
import numpy as np
from astropy.table import Table, hstack
from matplotlib import pyplot as plt

import AstroFunctions as astro


def apply_transforms_Warner(stars_table,
                            exoatmospheric_table_verify,
                            Warner_final_transform_table,
                            hidden_transform_table,
                            different_filter_list,
                            save_plots,
                            **kwargs):
    for different_filter in different_filter_list:
        app_mag, app_mag_sigma, app_filter, colour_index = astro.get_app_mag_and_index_AVG(
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
        fit, or_fit, line_init = astro.init_linear_fitting(
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


def exoatmospheric_mags_warner(stars_table,
                               extinction_table_warner,
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
                mask = ((extinction_table_warner['filter'] == different_filter) & (
                        extinction_table_warner['CI'] == colour_index))
                row_of_extinctions = extinction_table_warner[mask]
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
            app_mag, app_mag_sigma, app_mag_filter, _ = astro.get_app_mag_and_index(
                exoatmospheric_table, different_filter)
        except KeyError:
            exoatmospheric_table.rename_column('V_sigma', 'e_V')
            app_mag, app_mag_sigma, app_mag_filter, _ = astro.get_app_mag_and_index(
                exoatmospheric_table, different_filter)
        colour_indices = astro.get_all_colour_indices(different_filter)
        for colour_index in colour_indices:
            ci_plot = np.arange(min(exoatmospheric_table[colour_index]) - 0.1,
                                max(exoatmospheric_table[colour_index]) + 0.1,
                                step=0.01)
            exoatmospheric_mag = np.array(
                exoatmospheric_table[f"{different_filter} {colour_index}"])
            y = app_mag - exoatmospheric_mag
            x = exoatmospheric_table[colour_index]
            fit, or_fit, line_init = astro.init_linear_fitting(niter=100, sigma=2.5)
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
        colour_index, ci = astro.get_colour_index_lower(different_filter)
        try:
            ci0_index, _ = astro.get_colour_index_lower(ci[0])
            ci1_index, _ = astro.get_colour_index_lower(ci[1])
            positive_instr_mag = exoatmospheric_table[f"{ci[0]} {ci0_index}"]
            negative_instr_mag = exoatmospheric_table[f"{ci[1]} {ci1_index}"]
            table_ci = ci
        except KeyError:
            if 'v' in ci:
                table_ci = ci.replace('v', 'g')
            else:
                table_ci = ci
            try:
                ci0_index, _ = astro.get_colour_index_lower(ci[0])
                ci1_index, _ = astro.get_colour_index_lower(ci[1])
                positive_instr_mag =\
                    exoatmospheric_table[f"{table_ci[0]} {ci0_index}"]
                negative_instr_mag =\
                    exoatmospheric_table[f"{table_ci[1]} {ci1_index}"]
            except KeyError:
                if 'b' in ci:
                    table_ci = table_ci.replace('b', 'u')
                # table_ci = ci.replace('v', 'g')
                ci0_index, _ = astro.get_colour_index_lower(ci[0])
                ci1_index, _ = astro.get_colour_index_lower(ci[1])
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
        fit, or_fit, line_init = astro.init_linear_fitting(
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
    colour_index, ci = astro.get_colour_index_lower(different_filter)
    try:
        ci0_index, _ = astro.get_colour_index_lower(ci[0])
        ci1_index, _ = astro.get_colour_index_lower(ci[1])
        positive_instr_mag = stars_table[ci[0]]
        negative_instr_mag = stars_table[ci[1]]
        table_ci = ci
    except KeyError:
        if 'v' in ci:
            table_ci = ci.replace('v', 'g')
        else:
            table_ci = ci
        try:
            ci0_index, _ = astro.get_colour_index_lower(ci[0])
            ci1_index, _ = astro.get_colour_index_lower(ci[1])
            positive_instr_mag = stars_table[table_ci[0]]
            negative_instr_mag = stars_table[table_ci[1]]
        except KeyError:
            if 'b' in ci:
                table_ci = table_ci.replace('b', 'u')
            # table_ci = ci.replace('v', 'g')
            ci0_index, _ = astro.get_colour_index_lower(ci[0])
            ci1_index, _ = astro.get_colour_index_lower(ci[1])
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
            colour_indices = astro.get_all_colour_indices(different_filter)
            colour_index, _ = astro.get_colour_index_lower(different_filter)
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



def second_order_extinction_calc_warner(slopes_table,
                                        different_filter_list,
                                        save_plots, **kwargs):
    '''
    Calculate First and Second order Extinction Coeffieicients using Warners method 

    Parameters
    ----------
    slopes_table : astropy.table.table.Table
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
    different_filter_list : list
        different filter list
    save_plots : boolean
        Boolean to save plots
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    extinction_table_Warner : astropy.table.table.Table
        DESCRIPTION.

    '''
    filter_column = []
    CI_column = []
    k_primeprime_column = []
    k_prime_column = []
    k_primeprime_sigma_column = []
    k_prime_sigma_column = []
    for different_filter in different_filter_list:
        colour_indices = astro.get_all_colour_indices(different_filter)
        for colour_index in colour_indices:
            # print(different_filter)
            # print(colour_index)
            # print(slopes_table.pprint_all())
            filter_column.append(different_filter)
            CI_column.append(colour_index)
            # print(slopes_table[colour_index])

            #  FIXME: lenght of ci_plot might not be numerically stable. See Numpy Docs on arange
            ci_plot = np.arange(min(slopes_table[colour_index][~np.isnan(slopes_table[colour_index])]) - 0.1,
                                max(slopes_table[colour_index][~np.isnan(
                                    slopes_table[colour_index])]) + 0.1,
                                step=0.01)
            # fit, or_fit, line_init = init_linear_fitting(sigma=1.5)
            fit, or_fit, line_init = astro.init_linear_fitting(
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
def calculate_slopes_Warner(stars_table, different_filter_list, save_plots, **kwargs):
    '''
    Calculates the slope of Each Stars Magntitude with Time and plots it

    Parameters
    ----------
    stars_table :  astropy.table.table.Table
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

    different_filter_list : List
        Different Filter List
    save_plots : Boolean
        True= Save Plots
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    slopes_table : Table
        DESCRIPTION.

    '''

    stars_for_second_order_extinction, multiple_stars = astro.get_stars_with_multiple_observations(
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




