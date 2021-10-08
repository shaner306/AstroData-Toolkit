# -*- coding: utf-8 -*-
"""
Run the light curve generation process.
Created on Thu Jun  3 11:16:13 2021

@author: jmwawrow
"""
import AstroFunctions as astro
from astropy.io import ascii
import os

# directory = 'D:\\Transfer to mac\\2021-03-10 - Calibrated\\Intelsat 10-02 Post Eclipse\\LIGHT\\B_lim'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-03-20 - Calibrated\Intelsat 10-02\Test'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-04-21\Intelsat 10-02 ALL'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Observations\2016-072'
directory = r'F:\Intelsat 10-02\2020 11 30 - 2x2 - unprocessed\Intelsat 10-02\corrected_lights'

# transforms_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-04-21\Solved Stars\Outputs\gb_final_transforms.csv'
# transforms_file = r'F:\Intelsat 10-02\2020 11 30 - 2x2 - unprocessed\Stars\corrected_lights\Outputs\_gb_final_transforms.csv'
# gb_final_transforms = ascii.read(transforms_file)
# gb_final_transforms.pprint_all()
gb_final_transforms = None

subfolder_list = [
    '2020-12-28 - unprocessed', 
    '2021-03-10 - unprocessed', 
    '2021-03-20 - unprocessed', 
    '2021-03-31 - unprocessed'#, 
    # '2021-04-21 - unprocessed', 
    # '2021-04-25 - unprocessed'
    ]

temp_dir = 'tmp'
save_loc = 'Outputs_Warner'
file_suffix = '.fits'
ecct_cut=0.5
max_distance_from_sat = 20
size = 20
max_num_nan = 5
plot_results = 0

for subfolder in reversed(subfolder_list):
    print(subfolder)
    transforms_file = f'F:\\Intelsat 10-02\\{subfolder}\\Stars\\corrected_lights\\Outputs_Warner\\_gb_final_transforms.csv'
    gb_final_transforms = ascii.read(transforms_file)
    directory = f"F:\\Intelsat 10-02\\{subfolder}\\Intelsat 10-02\\corrected_lights"
    save_loc = os.path.join(directory, 'Outputs_Warner')
    sat_dict, app_sat_dict, sats_table, uncertainty_table, sat_auxiliary_table = astro._main_sc_lightcurve(directory, 
                                                                              gb_final_transforms=gb_final_transforms,
                                                                              temp_dir=temp_dir, 
                                                                              save_loc=save_loc,
                                                                              file_suffix=file_suffix,
                                                                              ecct_cut=ecct_cut,
                                                                              max_distance_from_sat=max_distance_from_sat, 
                                                                              size=size, 
                                                                              max_num_nan=max_num_nan, 
                                                                              plot_results=plot_results)

# print(sat_dict)
# sat_auxiliary_table.pprint_all()
# uncertainty_table.pprint_all()
# sats_table.pprint_all()