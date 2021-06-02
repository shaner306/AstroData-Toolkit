# -*- coding: utf-8 -*-
"""
Run the Space-based transform calculation.
Created on Mon May 31 12:35:34 2021

@author: jmwawrow
"""
import AstroFunctions as astro
import os

subfolder_list = [
    '2020-04-30',
    '2020-05-15',
    '2021-05-17',
    '2021-05-19',
    '2021-05-21',
    '2021-05-23'
    ]


ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\FITS Tutorial\Reference_stars_Apr29.txt'

plot_results = True
save_plots = False
file_suffix = "_clean.fits"
exposure_key = 'AEXPTIME'
name_key = 'Name'
transform_index_list = ['(B-V)', '(V-R)', '(V-I)']

for subfolder in subfolder_list:
    unique_id = subfolder
    directory = f'C:\\Users\\jmwawrow\\Documents\\DRDC_Code\\NEOSSat Landolt Stars\\{subfolder}'
    save_loc = os.path.join(directory, 'Outputs')
    
    sb_final_transform_table = astro._main_sb_transform_calc(directory, 
                                                             ref_stars_file, 
                                                             plot_results=plot_results, 
                                                             save_plots=save_plots,
                                                             file_suffix=file_suffix, 
                                                             exposure_key=exposure_key,  
                                                             name_key=name_key,
                                                             transform_index_list=transform_index_list,
                                                             save_loc=save_loc,
                                                             unique_id=unique_id)
    sb_final_transform_table.pprint_all()