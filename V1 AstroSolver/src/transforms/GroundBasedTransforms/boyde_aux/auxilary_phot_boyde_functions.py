# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:04:14 2022

@author: mstew

"""


#%% Boyde Slopes 1 


import numpy
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from astropy.modeling.fitting import  FittingWithOutlierRemoval, LinearLSQFitter
from astropy.modeling.models import Linear1D

from astropy.stats import sigma_clip


from matplotlib import pyplot as plt

import os
from astropy.table import Table


import AstroFunctions as astro
from astropy.io import ascii
from astropy.wcs import WCS


def calculate_boyde_slopes(filepath,Boyde_Table, save_plots, sav_loc,stars_table,img_filter,img_filters):
    '''


    Parameters
    ----------
    matched_stars : Table
        matched stars is a subsect of the reference stars which have arrays corresponding to the image stars- iraf sources detected in the image
    img_filter : str
        Jouhnson-Cousins Filter for the particular calculation file
    reference_stars_colnames : Table
        Lists the reference stars
    save_plots : bool
        boolean which determines whether the plots will be saved
    filepath : str
        describes file path of the image
    Boyde_Table : astropy.Table
        Includes:
            Image Name, 
            C prime,
            z prime, 
            Index,
            Airmass, 
            Airmass of the image

            Colour Filter:
                filter used when the image was taken


    Returns
    -------
    Boyde_Table : Astropy.Table
        Includes:
            Image Name:
                    the basename of the inputted image file. Describes the Image Name
            C:
                    The  slope produced from linear regression of the Standard 
                    Magntiude-intrumental Magnitude vs Colour Index of the Landolt stars.
                    This slope is used determine the colour transform term T in Steps 2 and 3
            Z-prime: 
                    the y- intercept of the linear regressed data. 
                    Used to determine k prime and zero point in Step 2 and 3
            Index: 
                    The colour index used in linear regression to calculate C and z Prime
            Airmass:
                    The predicted center airmass of the image
            Colour Filter:
                    Colour filter used in the capture of the image
            Step 1 Standard Deviation:
                    The standard Deviation produced by the residuals
                    of the plotted points against the fitted line   
            Number of Valid Matched Stars:
                    The number of stars matched in the image
    Outputs
    ------
    If save_plot==True:
        Save the Linear Regression Plots created with each image
            
            

    '''
    # hdr, imgdata = astro.read_fits_file(filepath)
    #
    # #TODO: Add Keyword Editor to except multiple
    #
    # try:
    #     predicted_airmass = float(hdr['AIRMASS'])
    # except KeyError:
    #     # AIRMASS header not found, attempt to estimate from WCS data
    #     wcs = WCS(hdr)
    #     altazpositions = astro.convert_ra_dec_to_alt_az((
    #         wcs.pixel_to_world(wcs.pixel_shape[0] / 2, wcs.pixel_shape[1] / 2)), hdr)
    #     predicted_airmass = altazpositions.secz.value
    #
    # img_filter=str(hdr['FILTER'])
    # if len(img_filter)>1:
    #     img_filter=img_filter.split(' ')[0] #for filter keywords such as 'B BESSEL'
    # airmass_std=np.sqrt(sum((matched_stars.img_star_airmass-predicted_airmass)**2)/(np.count_nonzero(matched_stars.img_star_airmass)-1))
    #
    #
    # airmass=predicted_airmass
    #
    #
    # fit = LinearLSQFitter(calc_uncertainties=True)
    # line_init = Linear1D()
    # or_fit = FittingWithOutlierRemoval(fit, sigma_clip, niter=300, sigma=3)
    #
    # # Create Table with headers
    # # Create a dict object of column names that contain the colour indices
    # # Keys are reference table column index
    # # Values are column names
    #
    # colour_indices = {}
    # colour_indices_errors={}
    # reference_magntiude_column = matched_stars.ref_star.colnames.index("V_ref")
    # reference_magntiude_column_error = matched_stars.ref_star.colnames.index("e_V")
    #
    #
    # for column_num, colname in enumerate(matched_stars.ref_star.colnames):                                              #Dynamically create colour_indices and colour_indices_errors
    #     if (('-' in colname) and ('e' not in colname)):
    #         colour_indices[column_num] = colname
    #     if (('-' in colname) and ('e' in colname)):
    #         colour_indices_errors[column_num]=colname
    #
    #
    #
    # # Iterate through indices
    # for colour_index in colour_indices:
    #
    #     # Iterate through macthed reference stars
    #     y_data = []
    #     x_data = []
    #     y_data_e=[]
    #     labels=[]
    #
    #     for row, matched_ref_star_vals in enumerate(matched_stars.ref_star):
    #
    #         # Reference Magntiude- instrumental magntiude
    #
    #         # Find Reference Magntitude:
    #         ref_mag = matched_ref_star_vals[reference_magntiude_column]
    #         ref_mag_e=matched_ref_star_vals[reference_magntiude_column_error]
    #
    #         # Instrumental Magnitude
    #
    #         ins_mag_sigma=matched_stars.img_instr_mag_sigma[row]
    #         matched_star_stars_table_row = np.where(stars_table['Name']==matched_stars.ref_star['Name'][row])
    #         #mag_unc_column = [i for i,colname in enumerate(stars_table.colnames) if str() == colname]
    #         #mag_unc_list=(stars_table[matched_star_stars_table_row][str(img_filter.lower()+'_sigma')])
    #
    #         # If multiple airmasses are in the same file dataset, chose the data that has the closest airmass to the
    #         # image
    #         closest_airmass=min(stars_table[matched_star_stars_table_row]['X_'+str(img_filter.lower())].value,key=lambda x:abs(
    # x-np.mean(matched_stars.img_star_airmass)))
    #         closest_airmass_ind = np.where(stars_table[matched_star_stars_table_row]['X_' + str(img_filter.lower())] == closest_airmass)
    #         mag_unc= (stars_table[matched_star_stars_table_row][str(img_filter.lower() + '_sigma')])[closest_airmass_ind]
    #
    #         #  Mean value of the instrumental magnitude calcualted in stars table
    #         try:
    #             ins_mag=(stars_table[matched_star_stars_table_row][img_filter.lower()])[closest_airmass_ind].value[0]
    #         except IndexError:
    #             continue
    #         #mag_unc=matched_stars.img_instr_mag_sigma[row]
    #
    #         if mag_unc==0 or mag_unc==np.nan: #or mag_unc<np.nanmean((stars_table[str(img_filter.lower(
    #             # )+'_sigma')])*0.25):
    #             # If magnitude uncertainty is 0 or np.nan skip entry or less than 0.25 of the mean exclude the value
    #             continue
    #
    #         # TODO: Restructure so it's easier to read
    #         # Difference between the two
    #         if (img_filter == 'V') or (img_filter == 'G'):
    #             diff_mag = ref_mag-ins_mag
    #             diff_mag_error=ref_mag_e
    #         else:
    #             '''
    #             Using the filter name, find and match the filter to an index in the the reference star column name
    #             i.e. An image with a B filter with us B to find the B-V index.
    #
    #             '''
    #
    #             if str(img_filter)+'-V' in colour_indices[colour_index] and str(img_filter)+'-V' in\
    #                     matched_stars.ref_star.colnames:
    #
    #                 matched_index=[i for i in range(len(matched_stars.ref_star.colnames)) if str(img_filter)+'-V' in matched_stars.ref_star.colnames[i]]
    #                 ref_filter_mag = matched_ref_star_vals[matched_index[0]] + ref_mag
    #                 diff_mag = ref_filter_mag - ins_mag
    #                 diff_mag_error = (matched_ref_star_vals[matched_index[1]] + ref_mag_e)
    #             elif 'V-'+ str(img_filter) in matched_stars.ref_star.colnames and 'V-'+ str(img_filter) in \
    #                     colour_indices[colour_index]:
    #                 matched_index = [i for i in range(len(matched_stars.ref_star.colnames)) if
    #                                  ('V-' + str(img_filter) in matched_stars.ref_star.colnames[i])]
    #                 ref_filter_mag = -(matched_ref_star_vals[matched_index[0]] - ref_mag)
    #                 diff_mag = ref_filter_mag - ins_mag
    #                 diff_mag_error = (matched_ref_star_vals[matched_index[1]] + ref_mag_e)
    #             else:
    #                 # No img filter in colour index
    #                 break
    #         if ((str(diff_mag_error+mag_unc) != 'nan')
    #             and (str(diff_mag) != 'nan')
    #             and (str(x_data) != 'nan')
    #             and ((('V-'+ str(img_filter) in colour_indices[colour_index])
    #                 or (str(img_filter) + '-V' in colour_indices[colour_index]))
    #                 or ((img_filter == 'V') or (img_filter == 'G') and  'V' in colour_indices[colour_index] )
    #                 )
    #         ):
    #
    #
    #             star_name=matched_ref_star_vals['Name']
    #             y_data.append(diff_mag)
    #             y_data_e.append((diff_mag_error+mag_unc.value[0]+ins_mag_sigma))
    #             labels.append(star_name)
    #             x_data.append(matched_stars.ref_star[row][colour_index])
    #
    #         else:
    #             continue
    #
    #     # Implement Linear Regression Model
    #
    #     #fit, or_fit, line_init=init_linear_fitting(sigma=2.0)
    #     #fitted_line, mask = fit(line_init, x_data, y_data)
    #     #x_nan_indices = np.isnan(x_data)
    #     #y_nan_indices = np.isnan(y_data)
    #     # fitted_line,mask=or_fit(line_init,x_data[~x_nan_indices],y_data[~y_nan_indices])
    #
    #     fitted_line1, mask = or_fit(
    #         line_init, np.array(x_data), np.array(y_data),weights=1.0/(np.array(y_data_e)**2))
    #     filtered_data = np.ma.masked_array(y_data, mask=mask)
    #
    #
    #     # TODO: Confirm with Don that this is the right equations
    #     # Calculate the residuals and accuracy of the computed data
    #     residuals = filtered_data - (fitted_line1((x_data)))
    #     filtered_data_e = np.ma.masked_array(y_data_e,mask=mask)
    #     std_residuals=np.sqrt(sum(residuals.data[np.where(residuals.mask==False)]**2)/(np.count_nonzero(residuals.mask == False)-1))
    #     se_residuals=std_residuals/np.sqrt(np.count_nonzero(residuals.mask == False))
    #     ave_data_e=np.mean(filtered_data_e)
    #     ave_std=abs(std_residuals)+abs(ave_data_e)
    grouped_stars=stars_table.group_by('Name')
    # Create run class



    rows = []
    for group in grouped_stars.groups:
        for row, star in enumerate(group):
            rows.append(row)
    grouped_stars.add_column(rows, name='run')
    grouped_stars = grouped_stars.group_by('run')

    colour_indices = []

    if img_filter=='G' or img_filter=='V':
        for colname in grouped_stars.colnames:
            for img_filter_instance in img_filters:

                if ('V-'+img_filter_instance.upper()) in colname and 'e' not in colname:
                    colour_indices.append(colname)
                elif (img_filter_instance.upper()+'-V') in colname and 'e' not in colname:
                    colour_indices.append(colname)
        for group in grouped_stars.groups:
            for colour_index in colour_indices:
                y_data = []
                x_data = []
                y_data_e = []
                labels = []
                airmasses = []

                for row, instr_mag in enumerate(group[img_filter.lower()]):
                    if str(instr_mag) != 'nan' and \
                            str(group[colour_index][row]) != 'nan' and \
                            str(img_filter.lower() + '_sigma') != 'nan' and \
                            str(group[img_filter.lower() + '_sigma'][row]) != 'nan' and (
                            img_filter == 'G' or img_filter == 'V'):
                                diff_mag = group['V_ref'][row] - instr_mag
                                colour_index_value = group[colour_index][row]
                                diff_mag_e = group[img_filter.lower() + '_sigma'][row]

                                average_mag_e = np.mean(group[img_filter.lower() + '_sigma'])
                                std_mag_e = np.std(group[img_filter.lower() + '_sigma'])
                                if diff_mag_e < average_mag_e - std_mag_e or diff_mag_e > average_mag_e + std_mag_e or \
                                        diff_mag_e == 0:
                                    continue
                                star_ave_airmass = group['X_' + img_filter.lower()][row]
                                if (isinstance(diff_mag, float) or isinstance(diff_mag, int)) and \
                                        (isinstance(colour_index_value, float) or isinstance(colour_index_value,
                                                                                             int)) and \
                                        (isinstance(diff_mag_e, float) or isinstance(diff_mag_e, int)) and \
                                        (isinstance(star_ave_airmass, float) or isinstance(star_ave_airmass, int)):
                                    x_data.append(colour_index_value)
                                    y_data.append(diff_mag)
                                    y_data_e.append(diff_mag_e)
                                    labels.append(group['Name'][row])
                                    airmasses.append(star_ave_airmass)

                    else:
                        continue
                try:
                    Boyde_Table = plot_bodye_slope1(Boyde_Table, airmasses, x_data, y_data, y_data_e,
                                                    colour_index,
                                                    img_filter,
                                                    labels, sav_loc, filepath, save_plots)
                    ascii.write(Boyde_Table, os.path.join(
                        sav_loc, 'Boyde_Table1.csv'), format='csv',overwrite=True)
                except:
                    print('Could Not Save Boyde Table')


    else:
        for colname in grouped_stars.colnames:
            if ('V-' + img_filter.upper()) in colname and 'e' not in colname:
                colour_indices.append(colname)
            elif (img_filter.upper() + '-V') in colname and 'e' not in colname:
                colour_indices.append(colname)


        for group in grouped_stars.groups:

            for colour_index in colour_indices:
                y_data = []
                x_data = []
                y_data_e = []
                labels = []
                airmasses=[]

                # TODO: Fix All of these redundant conditions

                for row,instr_mag in enumerate(group[img_filter.lower()]):
                    if str(instr_mag) != 'nan' and \
                            str(group[colour_index][row]) != 'nan' and \
                            str(img_filter.lower() + '_sigma') != 'nan' and \
                            str(group[img_filter.lower() + '_sigma'][row]) != 'nan':
                        if str(img_filter) + '-V' in colour_index:
                            diff_mag = group[colour_index][row]+group['V_ref'][row]  - instr_mag
                            colour_index_value=group[colour_index][row]
                            diff_mag_e=group[img_filter.lower() + '_sigma'][row]

                            # TODO: Delete this when better correlator/detector method is created
                            average_mag_e=np.mean(group[img_filter.lower() + '_sigma'])
                            std_mag_e=np.std(group[img_filter.lower() + '_sigma'])
                            if diff_mag_e<average_mag_e-std_mag_e or diff_mag_e>average_mag_e+std_mag_e or \
                                    diff_mag_e==0:
                                continue


                            star_ave_airmass=group['X_' + img_filter.lower()][row]
                            if (isinstance(diff_mag,float) or isinstance(diff_mag,int)) and \
                                (isinstance(colour_index_value,float) or isinstance(colour_index_value,int)) and \
                                (isinstance(diff_mag_e,float) or isinstance(diff_mag_e,int)) and \
                                    (isinstance(star_ave_airmass,float) or isinstance(star_ave_airmass,int)):

                                x_data.append(colour_index_value)
                                y_data.append(diff_mag)
                                y_data_e.append(diff_mag_e)
                                labels.append(group['Name'][row])
                                airmasses.append(star_ave_airmass)

                            else:
                                continue

                        elif 'V-' +str(img_filter) in colour_index:
                            diff_mag=-(group[colour_index][row]-group['V_ref'][row]) - instr_mag
                            diff_mag_e = group[img_filter.lower() + '_sigma'][row]

                            # FIXME: Delete this when better correlator/detector is implemented
                            average_mag_e = np.mean(group[img_filter.lower() + '_sigma'])
                            std_mag_e = np.std(group[img_filter.lower() + '_sigma'])
                            if diff_mag_e < average_mag_e - std_mag_e or diff_mag_e > average_mag_e + std_mag_e or \
                                    diff_mag_e == 0:
                                # Filters out low variability stars associated with small sample size
                                continue
                            colour_index_value = group[colour_index][row]
                            star_ave_airmass = group['X_' + img_filter.lower()][row]

                            if (isinstance(diff_mag, float) or isinstance(diff_mag, int)) and \
                                    (isinstance(colour_index_value, float) or isinstance(colour_index_value, int)) and \
                                    (isinstance(diff_mag_e, float) or isinstance(diff_mag_e, int)) and \
                                    (isinstance(star_ave_airmass, float) or isinstance(star_ave_airmass, int)):

                                x_data.append(colour_index_value)
                                y_data.append(diff_mag)
                                y_data_e.append(diff_mag_e)
                                labels.append(group['Name'][row])
                                airmasses.append(star_ave_airmass)
                            else:
                                continue
                        else:
                            continue
                try:
                    Boyde_Table=plot_bodye_slope1(Boyde_Table, airmasses, x_data, y_data, y_data_e, colour_index,
                                                              img_filter,
                                                  labels,sav_loc,filepath,save_plots)
                    ascii.write(Boyde_Table, os.path.join(
                        sav_loc, 'Boyde_Table1.csv'), format='csv', overwrite=True)

                except:
                    print('Could Not Save Boyde Table')


    return Boyde_Table


def plot_bodye_slope1(Boyde_Table,airmasses,x_data,y_data,y_data_e,colour_index,img_filter,labels,
                      sav_loc,filepath,save_plots=True):
    airmass = np.mean(airmasses)
    airmass_std = np.std(airmasses)

    fit = LinearLSQFitter(calc_uncertainties=True)
    line_init = Linear1D()
    or_fit = FittingWithOutlierRemoval(fit, sigma_clip, niter=300, sigma=3)

    # fit, or_fit, line_init=astro.init_linear_fitting(sigma=2.0)
    # fitted_line, mask = fit(line_init, x_data, y_data)
    # fitted_line,mask=or_fit(line_init,x_data,y_data)
    fitted_line1, mask = or_fit(
        line_init, np.array(x_data), np.array(y_data), weights=1.0 / (np.array(y_data_e) ** 2))
    filtered_data = np.ma.masked_array(y_data, mask=mask)
    residuals = filtered_data - (fitted_line1((x_data)))

    filtered_data_e = np.ma.masked_array(y_data_e, mask=mask)
    std_residuals = np.sqrt(
        sum(residuals.data[np.where(residuals.mask == False)] ** 2) / (np.count_nonzero(residuals.mask == False) - 1))
    se_residuals = std_residuals / np.sqrt(np.count_nonzero(residuals.mask == False))
    ave_data_e = np.mean(filtered_data_e)
    ave_std = abs(std_residuals) + abs(ave_data_e)

    # Plotting #
    if save_plots and x_data != [] and y_data != [] and y_data_e != [] and labels != []:
        # print('Save Plots')

        plt.figure()
        plt.errorbar(x_data, y_data, yerr=y_data_e, fmt='ko', fillstyle='none', label='Clipped Data')
        # plt.plot(x_data, y_data, 'ro',
        #         fillstyle='none', label='Clipped Data')
        plt.plot(x_data, filtered_data, 'ro', label='Filtered Data')
        plt.plot(x_data, fitted_line1(x_data), 'k:',
                 label=f"Z': {fitted_line1.intercept.value:.3f}, C: {fitted_line1.slope.value:.3f}")
        # plt.plot(x_data, fitted_line1(x_data), 'k:', label=("Z':"+str(fitted_line1.intercept.value)+"C:"+str(fitted_line1.slope.value)))
        plt.xlabel(colour_index)
        plt.ylabel(img_filter.upper() + '_ref' + '-' + img_filter.lower() + '_inst')

        i = 0
        for x, y in zip(x_data, y_data):
            label = labels[i]
            i = i + 1
            plt.annotate(label, (x, y), textcoords="offset points", xytext=(0, 10), ha='center')

        plt.legend()

        plt.title('Step 1: V_ref-v_inst vs.' +
                  colour_index + '_in Filter ' + img_filter + " With image " + os.path.basename(filepath))

        plt.savefig(str(sav_loc) + 'Boyde_step_1_' +
                    str(colour_index) + '_' + str(img_filter) + os.path.basename(filepath) + '.png')

        plt.close('all')

    # Boyde_Table=Table(names=['Image Name','C prime','Z-prime','Index (i.e. B-V)','Z Prime','Airmass','Colour Filter'])

    Boyde_Table.add_row([os.path.basename(filepath), fitted_line1.slope.value,
                         fitted_line1.intercept.value, colour_index, airmass, airmass_std, img_filter, ave_std,
                         np.count_nonzero(residuals.mask == False)])
    return Boyde_Table

#%% Boyde Slopes 2
def calculate_boyde_slope_2(Boyde_Table,save_loc,match_stars_lim, save_plots=True):
    '''
    See calculate_boyde_slope_1
    
    Calculates Step 2 and 3 in the Boyde Method and adds them to the previously built Boyde Table in Step 1
    
    Parameters
    ---------
    save_plots
    Inputs:
        Boyde_Table
    '''

    Boyde_Table_grouped = Boyde_Table.group_by(
        ['Colour Filter', 'Index (i.e. B-V)'])

    # Initialize Linear Regression Model
    fit = LinearLSQFitter(calc_uncertainties=True)
    line_init = Linear1D()
    or_fit = FittingWithOutlierRemoval(fit, sigma_clip, niter=300, sigma=3)

    k_prime_prime_array = []
    colour_transform_array = []
    step2_uncertainty_array=[]
    k_prime_array = []
    zero_point_array = []
    step3_uncertainty_array=[]

    for i, index1 in enumerate(Boyde_Table_grouped.groups.indices[1:]):
        x_data = []
        y_data = []
        y_data_e=[]
        x_data_e=[]
        for image in np.arange(Boyde_Table_grouped.groups.indices[i], index1):
            if Boyde_Table_grouped[image]['Number of Valid Matched Stars'] >= match_stars_lim:

                y_data.append(Boyde_Table_grouped[image]['C'])
                x_data.append(Boyde_Table_grouped[image]['Airmass'])
                y_data_e.append(Boyde_Table_grouped[image]['Step1_Standard_Deviation'])
                x_data_e.append(Boyde_Table_grouped[image]['Air Mass Std'])
        # Fit the Data
        fitted_line1, mask = or_fit(
            line_init, np.array(x_data), np.array(y_data),weights=1/(np.array(y_data_e)**2))
        filtered_data = np.ma.masked_array(y_data, mask=mask)

        # Get K Prime Prime and Colour transform
        k_prime_prime = fitted_line1.slope.value
        colour_transform = fitted_line1.intercept.value

        # Get Uncertainty from Slope 2 fitted line
        residuals = filtered_data - (fitted_line1((x_data)))
        filtered_data_e = np.ma.masked_array(y_data_e,mask=mask)
        std_residuals=np.sqrt(sum(residuals.data[np.where(residuals.mask==False)]**2)/(np.count_nonzero(residuals.mask == False)-1))

        if save_plots and x_data!=[] and y_data!=[] and x_data_e !=[] and y_data_e !=[] :

            plt.figure()
            plt.errorbar(x_data, y_data, yerr=y_data_e, xerr=x_data_e, fmt='ko', fillstyle='none', label='Clipped Data')

            plt.plot(x_data, filtered_data, 'ro', label='Filtered Data')
            plt.plot(x_data, fitted_line1(x_data), 'k:', label=('k":'+str(k_prime_prime)+'T:'+str(colour_transform)))
            plt.xlabel('Mean Airmass')
            plt.ylabel('C')
            plt.legend()

            plt.title('Boydes Second Slope of ' +
                      Boyde_Table_grouped[Boyde_Table_grouped.groups.indices[i]]['Colour Filter'] + ' ' + Boyde_Table_grouped[Boyde_Table_grouped.groups.indices[i]]['Index (i.e. B-V)'])

            plt.savefig(str(save_loc)+'Boyde_step_2_' + str(Boyde_Table_grouped[Boyde_Table_grouped.groups.indices[i]]['Index (i.e. B-V)'])+'_'+str(
                Boyde_Table_grouped[Boyde_Table_grouped.groups.indices[i]]['Colour Filter'])+'.png')
            #plt.show(block=False)
            plt.clf()
            plt.close('all')

        for image in np.arange(Boyde_Table_grouped.groups.indices[i], index1):
            k_prime_prime_array.append(k_prime_prime)
            colour_transform_array.append(colour_transform)
            step2_uncertainty_array.append(std_residuals)

    # Step 3
    for i, index1 in enumerate(Boyde_Table_grouped.groups.indices[1:]):
        x_data = []
        y_data = []
        y_data_e2=[]
        x_data_e2=[]
        for image in np.arange(Boyde_Table_grouped.groups.indices[i], index1):
            if Boyde_Table_grouped[image]['Number of Valid Matched Stars'] >= match_stars_lim:
                if str(Boyde_Table_grouped[image]['Step1_Standard_Deviation'])=='inf' or str(Boyde_Table_grouped[image]['Step1_Standard_Deviation'])=='nan':
                    continue
                x_data.append(Boyde_Table_grouped[image]['Airmass'])
                y_data.append(Boyde_Table_grouped[image]['Z-prime'])
                y_data_e2.append(Boyde_Table_grouped[image]['Step1_Standard_Deviation'])
                x_data_e2.append(Boyde_Table_grouped[image]['Air Mass Std'])

        fitted_line2, mask = or_fit(
            line_init, np.array(x_data), np.array(y_data),weights=1/(np.array(y_data_e2)**2))
        filtered_data = np.ma.masked_array(y_data, mask=mask)


        # Calcualte Residuals
        residuals = filtered_data - (fitted_line2((x_data)))
        std_residuals=np.sqrt(sum(residuals.data[np.where(residuals.mask==False)]**2)/(np.count_nonzero(residuals.mask == False)-1))



        # Get Important Values
        k_prime = fitted_line2.slope.value
        zero_point = fitted_line2.intercept.value

        if save_plots and x_data!=[] and x_data_e2!=[] and y_data !=[] and y_data_e2 !=[]:

            plt.figure()
# =============================================================================
#             plt.plot(x_data, y_data, 'ro',
#                      fillstyle='none', label='Clipped Data')
#
# =============================================================================
            plt.errorbar(x_data, y_data, yerr=y_data_e2,xerr=x_data_e2, fmt='ko', fillstyle='none', label='Clipped Data')
            plt.plot(x_data, filtered_data, 'ro', label='Filtered Data')
            plt.plot(x_data, fitted_line2(x_data), 'k:', label=f"k': {k_prime:.3f}, Zp: {zero_point:.3f}")
            # plt.plot(x_data, fitted_line2(x_data), 'k:', label=("k':"+str(k_prime)+'Zp: '+ str(zero_point)))
# label=f'C_{unique_filter}{ci_plot} = {kprimeprime_fci:.3f} * X + {t_fci:.3f}')
            plt.xlabel('Mean Airmass')
            plt.ylabel('Z_prime')
            plt.legend()
            plt.title('Boydes Third Slope of ' +
                      Boyde_Table_grouped[Boyde_Table_grouped.groups.indices[i]]['Colour Filter'] + ' ' + Boyde_Table_grouped[Boyde_Table_grouped.groups.indices[i]]['Index (i.e. B-V)'])


            plt.savefig(str(save_loc)+'Boyde_step_3_' + str(Boyde_Table_grouped[Boyde_Table_grouped.groups.indices[i]]['Index (i.e. B-V)'])+'_'+str(
                Boyde_Table_grouped[Boyde_Table_grouped.groups.indices[i]]['Colour Filter'])+'.png')
            #plt.show()
            plt.close()

        for image in np.arange(Boyde_Table_grouped.groups.indices[i], index1):
            k_prime_array.append(k_prime)
            zero_point_array.append(zero_point)
            step3_uncertainty_array.append(std_residuals)

    Boyde_Table_grouped['k_prime_prime'] = k_prime_prime_array
    Boyde_Table_grouped['colour_transform'] = colour_transform_array
    Boyde_Table_grouped['k_prime'] = k_prime_array
    Boyde_Table_grouped['zero_point'] = zero_point_array
    Boyde_Table_grouped['step_2_error']=step2_uncertainty_array
    Boyde_Table_grouped['step_3_error']=step3_uncertainty_array

    return Boyde_Table_grouped


def create_coefficeint_output(Boyde_table_grouped,save_loc):
    '''
    Calculates the coefficients for the particular data
    Useful for comparing different dataset

    Parameters
    ----------
    Boyde_table_grouped:
       Image_Name:
        Image Name
       C:
        Slope value produced from Boyde Step 1
       Z-prime:
        The y-intercept of the slope produced from Boyde Step 1

       Index:
        The Colour index used in the calculations
       Airmass:
        The Predicted Airmass of the image stemmin
       Air Mass Standard Deviation
       Colour Filter
       Slope 1 Standard Deviaiton
       Number of Matched Stars
       k_prime_prime
       colour_transform
       k_prime
       zero_point


    Returns
    -------
    coefficient_data:
        k_prime_bbv
        k-prime_vbv
        k_prime_rvr
        k_prime_prime_bbv
        k_prime_prime_vbv
        k_prime_prime_rvr
        T_bbv
        T_vbv
        T_rvr
        Z_bbv
        Z_vbv
        Z_rvr
    '''
    coefficient_data=Table(
        names=['Coefficient','Value','Error'],
        dtype=['str','float64','float64'])

    # TODO: Add Dydnamic Version

    # Work around since Table.loc can only handle one index
    rows=[row for row in np.array(Boyde_table_grouped) if
                    ('B' in row and 'B-V' in row)]
    #Add bbv data

    # TODO: Check to see if this error is right.
    index_bbv=Table(np.array(rows))[0]
    coefficient_data.add_row(['k_prime_bbv',index_bbv['k_prime'],index_bbv[
        'step_3_error']])
    coefficient_data.add_row(['T_bbv', index_bbv['colour_transform'],index_bbv[
        'step_2_error']])
    coefficient_data.add_row(['k_prime_prime_bbv', index_bbv['k_prime_prime'],
                        index_bbv['step_2_error']])
    coefficient_data.add_row(['Z_bbv', index_bbv['zero_point'],index_bbv['step_3_error']])

    rows = [row for row in np.array(Boyde_table_grouped) if
            ('G' in row and 'B-V' in row)]
    index_vbv = Table(np.array(rows))[0]
    # Add k_prime_vbv
    coefficient_data.add_row(['k_prime_vbv',index_vbv['k_prime'],index_vbv[
        'step_3_error']])
    coefficient_data.add_row(['T_vbv', index_vbv['colour_transform'], index_vbv[
        'step_2_error']])
    coefficient_data.add_row(['k_prime_prime_vbv', index_vbv['k_prime_prime'], index_vbv[
        'step_2_error']])
    coefficient_data.add_row(['Z_vbv', index_vbv['zero_point'], index_vbv[
        'step_3_error']])

    rows = [row for row in np.array(Boyde_table_grouped) if
            ('R' in row and 'V-R' in row)]
    index_rvr = Table(np.array(rows))[0]
    # Add k_prime_rvr
    coefficient_data.add_row(['k_prime_rvr',index_rvr['k_prime'],index_rvr[
        'step_3_error']])
    coefficient_data.add_row(['T_rvr', index_rvr['colour_transform'], index_rvr[
        'step_2_error']])
    coefficient_data.add_row(['k_prime_prime_rvr', index_rvr['k_prime_prime'], index_rvr[
        'step_2_error']])
    coefficient_data.add_row(['Z_rvr', index_rvr['zero_point'], index_rvr[
        'step_3_error']])

    # Work around to easily sort data
    # FIXME: Get sorted_coefficient_data to output correctly
    #coefficient_df=coefficient_data.argsort(keys='Coefficients')
    #coefficient_df.sort('Cofficient')
    #coefficient_data=coefficient_df.from_pandas()
    coefficient_data.sort(['Coefficient'])
    ascii.write(coefficient_data, os.path.join(
        str(os.path.dirname(save_loc)), 'Coefficient_data.csv'), format='csv', overwrite=True)

    return coefficient_data





