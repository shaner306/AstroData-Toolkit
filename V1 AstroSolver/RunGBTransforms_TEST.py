# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 10:49:48 2021

@author: jmwawrow
"""

import AstroFunctions as astro
import os
from astropy.io import ascii

subfolder_list = [
    '2021-03-20 - unprocessed'
    ]

subfolder_list = [
    '2020 10 23 - 2x2 - unprocessed', 
    '2020 10 31 - 2x2 - unprocessed', 
    '2020 11 04 - 2x2 - unprocessed', 
    '2020 11 10 - 2x2 - unprocessed', 
    '2020 11 21 - 2x2 - unprocessed', 
    '2020 11 30 - 2x2 - unprocessed', 
    '2020-12-12 - unprocessed', 
    '2020-12-24 - unprocessed', 
    '2020-12-25 - unprocessed', 
    '2020-12-28 - unprocessed', 
    '2021-02-07 - unprocessed', 
    '2021-03-10 - unprocessed', 
    '2021-03-20 - unprocessed', 
    '2021-03-21 - unprocessed', 
    '2021-03-22 - unprocessed', 
    '2021-03-23 - unprocessed', 
    '2021-03-31 - unprocessed', 
    '2021-04-21 - unprocessed', 
    '2021-04-24 - unprocessed', 
    '2021-04-25 - unprocessed'
    ]


ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\Combined_Ref_Stars.txt'

plot_results = True
save_plots = True
file_suffix = ".fits"
exposure_key = 'EXPTIME'
name_key = 'Name'
transform_index_list = ['(B-V)', '(V-R)', '(V-I)']

for subfolder in subfolder_list:
    unique_id = subfolder
    directory = f'F:\\Intelsat 10-02\\{subfolder}\\Stars\\corrected_lights'
    save_loc = os.path.join(directory, 'Outputs_TEST')
    
    sb_final_transform_table = astro._main_gb_transform_calc_TEST(directory, 
                                                             ref_stars_file, 
                                                             plot_results=plot_results, 
                                                             save_plots=save_plots,
                                                             file_suffix=file_suffix, 
                                                             exposure_key=exposure_key,  
                                                             name_key=name_key,
                                                             transform_index_list=transform_index_list,
                                                             save_loc=save_loc,
                                                             unique_id=unique_id)
    # sb_final_transform_table.pprint(max_lines=40)
    ascii.write(sb_final_transform_table, os.path.join(save_loc, 'large_stars_table.csv'), format='csv')