# -*- coding: utf-8 -*-
"""
Run the light curve generation process.
Created on Thu Jun  3 11:16:13 2021

@author: jmwawrow
"""
import AstroFunctions as astro
from astropy.io import ascii

# directory = 'D:\\Transfer to mac\\2021-03-10 - Calibrated\\Intelsat 10-02 Post Eclipse\\LIGHT\\B_lim'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-03-20 - Calibrated\Intelsat 10-02\Test'
directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-04-21\Intelsat 10-02 ALL'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Observations\2016-072'

transforms_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-04-21\Solved Stars\Outputs\gb_final_transforms.csv'
gb_final_transforms = ascii.read(transforms_file)
# gb_final_transforms.pprint_all()

temp_dir = 'tmp'
file_suffix = '.fits'
ecct_cut=0.5
max_distance_from_sat = 20
size = 20
max_num_nan = 5
plot_results = 0

sat_dict, sats_table, uncertainty_table, sat_auxiliary_table = astro._main_sc_lightcurve(directory, 
                                                                          gb_final_transforms=gb_final_transforms,
                                                                          temp_dir=temp_dir, 
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