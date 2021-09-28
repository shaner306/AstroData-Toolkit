# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 11:26:24 2021

@author: jmwawrow
"""
import AstroFunctions as astro
import os
from astropy.io import ascii

ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\Combined_Ref_Stars.txt'
# ref_stars_file = "C:\\Users\\jmwawrow\\Documents\\GitHub\\Astro2\\Reference Star Files\\Reference_stars_Apr29.txt"

plot_results = True
save_plots = True
file_suffix = ".fits"
exposure_key = 'EXPTIME'
name_key = 'Name'

# for subfolder in subfolder_list:
subfolder = '2021-09-17'
print(subfolder)
unique_id = subfolder
# directory = f'F:\\Intelsat 10-02\\{subfolder}\\Stars\\corrected_lights'
directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-09-17\corrected_lights\Calculation'
save_loc = os.path.join(directory, 'Outputs_TESTING')

large_stars_table = astro._main_gb_transform_calc_Warner(directory, 
                                                         ref_stars_file, 
                                                         plot_results=plot_results, 
                                                         save_plots=save_plots,
                                                         file_suffix=file_suffix, 
                                                         exposure_key=exposure_key,  
                                                         name_key=name_key,
                                                         save_loc=save_loc,
                                                         unique_id=unique_id)