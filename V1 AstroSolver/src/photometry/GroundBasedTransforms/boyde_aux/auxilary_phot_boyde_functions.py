# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:04:14 2022

@author: mstew

"""


#%% Boyde Slopes 1 


import numpy
import numpy as np

from astropy.modeling.fitting import  FittingWithOutlierRemoval, LinearLSQFitter
from astropy.modeling.models import Linear1D

from astropy.stats import sigma_clip

from matplotlib import pyplot as plt

import os


import AstroFunctions as astro

def calculate_boyde_slopes(matched_stars, filepath, Boyde_Table, save_plots, sav_loc,stars_table):
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
    hdr, imgdata = astro.read_fits_file(filepath)

    predicted_airmass = float(hdr['AIRMASS'])
    img_filter=str(hdr['FILTER'])
    
    airmass_std=np.sqrt(sum((matched_stars.img_star_airmass-predicted_airmass)**2)/(np.count_nonzero(matched_stars.img_star_airmass)-1))
    
    
    airmass=predicted_airmass
    
    # Create Table with headers

    # Create a dict object of column names that contain the colour indices
    # Keys are reference table column index
    # Values are column names

    fit = LinearLSQFitter(calc_uncertainties=True)
    line_init = Linear1D()
    or_fit = FittingWithOutlierRemoval(fit, sigma_clip, niter=300, sigma=3)

    colour_incides = {}
    colour_indices_errors={}
    reference_magntiude_column = matched_stars.ref_star.colnames.index("V_ref")
    reference_magntiude_column_error = matched_stars.ref_star.colnames.index("e_V")
    for column_num, colname in enumerate(matched_stars.ref_star.colnames):
        if (('-' in colname) and ('e' not in colname)):
            colour_incides[column_num] = colname
        if (('-' in colname) and ('e' in colname)):
            colour_indices_errors[column_num]=colname



    # Iterate through indices
    for colour_index in colour_incides:

        # Iterate through macthed reference stars
        y_data = []
        x_data = []
        y_data_e=[]

        for row, matched_ref_star_vals in enumerate(matched_stars.ref_star):

            # Reference Magntiude- instrumental magntiude

            # Find Reference Magntitude:
            ref_mag = matched_ref_star_vals[reference_magntiude_column]
            ref_mag_e=matched_ref_star_vals[reference_magntiude_column_error]

            # Instrumental Magnitude
            ins_mag = matched_stars.img_instr_mag[row]

            mag_unc_row = np.where(stars_table['Name']==matched_stars.ref_star['Name'][row])
            #mag_unc_column = [i for i,colname in enumerate(stars_table.colnames) if str() == colname]
            mag_unc=np.max(stars_table[mag_unc_row][str(img_filter.lower()+'_sigma')])
            if mag_unc==0 or mag_unc==np.nan or mag_unc<np.nanmean((stars_table[str(img_filter.lower()+'_sigma')])*0.25):
                # If magnitude uncertainty is 0 or np.nan skip entry or less than 0.25 of the mean exclude the value
                continue


            # Difference between the two
            if (img_filter == 'V') or (img_filter == 'G'):
                diff_mag = ref_mag-ins_mag
                diff_mag_error=ref_mag_e
            else:
                matched_index=[i for i in range(len(matched_stars.ref_star.colnames)) if str(img_filter)+'-V' in matched_stars.ref_star.colnames[i]]
                if len(matched_index) ==0:
                    matched_index=[i for i in range(len(matched_stars.ref_star.colnames)) if ('V-'+ str(img_filter) in matched_stars.ref_star.colnames[i])]
                    if len(matched_index)==0:
                        print('Could not find corresponding index')
                        break
                    else:
                        ref_filter_mag=-(matched_ref_star_vals[matched_index[0]]-ref_mag)
                        diff_mag=ref_filter_mag-ins_mag
                        diff_mag_error=(matched_ref_star_vals[matched_index[1]]+ref_mag_e)

                    
                    
                else:
                    ref_filter_mag=matched_ref_star_vals[matched_index[0]]+ref_mag
                    diff_mag=ref_filter_mag-ins_mag
                    diff_mag_error=(matched_ref_star_vals[matched_index[1]]+ref_mag_e)

                    
                    
               
                    
            y_data.append(diff_mag)
            y_data_e.append(diff_mag_error+mag_unc)
            x_data.append(matched_stars.ref_star[row][colour_index])

        # Implement Linear Regression Model

        #fit, or_fit, line_init=init_linear_fitting(sigma=2.0)
        #fitted_line, mask = fit(line_init, x_data, y_data)
        #x_nan_indices = np.isnan(x_data)
        #y_nan_indices = np.isnan(y_data)
        # fitted_line,mask=or_fit(line_init,x_data[~x_nan_indices],y_data[~y_nan_indices])

        fitted_line1, mask = or_fit(
            line_init, np.array(x_data), np.array(y_data),weights=1.0/(np.array(y_data_e)**2))
        filtered_data = np.ma.masked_array(y_data, mask=mask)


        # TODO: Confirm with Don that this is the right equations
        # Calculate the residuals and accuracy of the computed data
        residuals = filtered_data - (fitted_line1((x_data)))
        filtered_data_e = np.ma.masked_array(y_data_e,mask=mask)
        std_residuals=np.sqrt(sum(residuals.data[np.where(residuals.mask==False)]**2)/(np.count_nonzero(residuals.mask == False)-1))
        se_residuals=std_residuals/np.sqrt(np.count_nonzero(residuals.mask == False))
        ave_data_e=np.mean(filtered_data_e)
        ave_std=abs(std_residuals)+abs(ave_data_e)


        # Plotting #
        if save_plots:
            # print('Save Plots')

            plt.figure()
            plt.errorbar(x_data,y_data,yerr=y_data_e,fmt='ko',fillstyle='none',label='Clipped Data')
            #plt.plot(x_data, y_data, 'ro',
            #         fillstyle='none', label='Clipped Data')
            plt.plot(x_data, filtered_data, 'ro', label='Filtered Data')
            plt.plot(x_data, fitted_line1(x_data), 'k:', label=("Z':"+str(fitted_line1.intercept.value)+"C:"+str(fitted_line1.slope.value)))
            plt.xlabel(colour_incides[colour_index])
            plt.ylabel('V_ref-v_inst')
            plt.legend()
            
            plt.title('Step 1: V_ref-v_inst vs.' +
                      colour_incides[colour_index] + '_in Filter ' + img_filter +" With image " + os.path.basename(filepath))

            plt.savefig(str(sav_loc)+'Boyde_step_1_' +
                        str(colour_incides[colour_index])+'_'+str(img_filter)+ os.path.basename(filepath) +'.png')

            plt.close()

        #Boyde_Table=Table(names=['Image Name','C prime','Z-prime','Index (i.e. B-V)','Z Prime','Airmass','Colour Filter'])

        Boyde_Table.add_row([os.path.basename(filepath), fitted_line1.slope.value,
                            fitted_line1.intercept.value, colour_incides[colour_index], airmass,airmass_std, img_filter,ave_std,np.count_nonzero(residuals.mask == False)])

    return Boyde_Table

#%% Boyde Slopes 2

def calculate_boyde_slope_2(Boyde_Table, save_plots, save_loc,match_stars_lim):
    '''
    See calculate_boyde_slope_1
    
    Calculates Step 2 and 3 in the Boyde Method and adds them to the previously built Boyde Table in Step 1
    
    Parameters
    ---------
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
        
        if save_plots is True:

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
            #plt.show()
            plt.close()

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

        if save_plots is True:

            plt.figure()
# =============================================================================
#             plt.plot(x_data, y_data, 'ro',
#                      fillstyle='none', label='Clipped Data')
#             
# =============================================================================
            plt.errorbar(x_data, y_data, yerr=y_data_e2,xerr=x_data_e2, fmt='ko', fillstyle='none', label='Clipped Data')
            plt.plot(x_data, filtered_data, 'ro', label='Filtered Data')
            plt.plot(x_data, fitted_line2(x_data), 'k:', label=("k':"+str(k_prime)+'Zp: '+ str(zero_point)))
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

    return Boyde_Table_grouped

