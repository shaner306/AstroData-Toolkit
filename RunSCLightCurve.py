# -*- coding: utf-8 -*-
"""
Run the light curve generation process..
Created on Thu Jun  3 11:16:13 2021

@author: jmwawrow
"""
import AstroFunctions as astro

directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-03-20 - Calibrated\Intelsat 10-02\LIGHT'
temp_dir = 'tmp'
max_distance_from_sat = 20
size = 20
max_num_nan = 5
plot_results = 0

sats_table, uncertainty_table, sat_fwhm_table = astro._main_sc_lightcurve(directory, 
                                                                          temp_dir=temp_dir, 
                                                                          max_distance_from_sat=max_distance_from_sat, 
                                                                          size=size, 
                                                                          max_num_nan=max_num_nan, 
                                                                          plot_results=plot_results)

sat_fwhm_table.pprint_all()
uncertainty_table.pprint_all()
sats_table.pprint_all()