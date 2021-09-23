# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 10:25:21 2021

@author: jmwawrow
"""

import AstroFunctions as astro
from astropy.io import ascii
import os

ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\Combined_Ref_Stars.txt'

# directory = r'F:\Intelsat 10-02\2021-04-25 - unprocessed\Stars\corrected_lights'
# transforms_file = r'F:\Intelsat 10-02\2021-04-25 - unprocessed\Stars\corrected_lights\Outputs\_gb_final_transforms.csv'
# gb_final_transforms = ascii.read(transforms_file)

# subfolder_list = [
#     '2020 10 23 - 2x2 - unprocessed', 
#     '2020 10 31 - 2x2 - unprocessed',
#     '2020 11 04 - 2x2 - unprocessed',
#     '2020 11 10 - 2x2 - unprocessed',
#     '2020 11 21 - 2x2 - unprocessed',
#     '2020 11 30 - 2x2 - unprocessed',
#     '2020-12-12 - unprocessed', 
#     '2020-12-24 - unprocessed', 
#     '2020-12-25 - unprocessed', 
#     '2020-12-28 - unprocessed', 
#     '2021-02-07 - unprocessed',
#     '2021-03-10 - unprocessed', 
#     '2021-03-20 - unprocessed', 
#     '2021-03-21 - unprocessed', 
#     '2021-03-22 - unprocessed',
#     '2021-03-23 - unprocessed', 
#     '2021-03-31 - unprocessed', 
#     '2021-04-21 - unprocessed', 
#     '2021-04-24 - unprocessed', 
#     '2021-04-25 - unprocessed'
#     ]

# subfolder_list = [
#     '2021-04-21 - unprocessed', 
#     '2021-04-24 - unprocessed', 
#     '2021-04-25 - unprocessed'
#     ]

# Subfolders with bad/without transforms:
# subfolder_list = [
#     '2020 11 04 - 2x2 - unprocessed',
#     '2020 11 10 - 2x2 - unprocessed',
#     '2020 11 21 - 2x2 - unprocessed',
#     '2020 11 30 - 2x2 - unprocessed',
#     '2021-02-07 - unprocessed',
#     '2021-03-09 - unprocessed',
#     '2021-03-23 - unprocessed'
#     ]
# 2x2 binned with good transforms:
# subfolder_list = [
#     '2020 10 23 - 2x2 - unprocessed', 
#     '2020 10 31 - 2x2 - unprocessed'
#     ]

plot_results = True
save_plots = True
file_suffix = ".fits"
exposure_key = 'EXPTIME'
name_key = 'Name'

# for subfolder in subfolder_list:
#     print(subfolder)
# directory = f"F:\\Intelsat 10-02\\{subfolder}\\Stars\\corrected_lights"
directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-09-17\corrected_lights\Verification'
save_loc = os.path.join(directory, 'Outputs', 'VERIFICATION')
# transforms_file = f"F:\\Intelsat 10-02\\{subfolder}\\Stars\\corrected_lights\\Outputs\\_gb_final_transforms.csv"
transforms_file = "C:\\Users\\jmwawrow\\Documents\\DRDC_Code\\Intelsat 10-02\\2021-09-17\\corrected_lights\\Calculation\\Outputs\\_gb_final_transforms.csv"
gb_final_transforms = ascii.read(transforms_file)
app_mag_table = astro.verify_gb_transforms(directory, 
                                            ref_stars_file, 
                                            gb_final_transforms,
                                            plot_results=plot_results, 
                                            save_plots=save_plots,
                                            file_suffix=file_suffix, 
                                            exposure_key=exposure_key,  
                                            name_key=name_key,
                                            save_loc=save_loc)