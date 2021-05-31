# -*- coding: utf-8 -*-
"""
Run the Ground-based transform calculation.

Created on Mon May 31 11:51:04 2021

@author: jmwawrow
"""
import AstroFunctions as astro
import os

directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021_J132_46927_DESCENT\May 18 2021\Landolt Fields\Solved'
ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\FITS Tutorial\Reference_stars_Apr29.txt'

plot_results = True
save_plots = False
remove_large_airmass = False
file_suffix=".fit"
exposure_key='EXPTIME'
lat_key='OBSGEO-B'
lon_key='OBSGEO-L'
elev_key='OBSGEO-H'
# file_suffix=".fits"
# exposure_key='EXPTIME'
# lat_key='SITELAT'
# lon_key='SITELONG'
# elev_key='SITEELEV'
name_key='Name'

save_loc = os.path.join(directory, 'Outputs')
gb_final_transforms = astro._main_gb_transform_calc(directory, 
                                                    ref_stars_file, 
                                                    plot_results=plot_results, 
                                                    save_plots=save_plots, 
                                                    remove_large_airmass=remove_large_airmass, 
                                                    file_suffix=file_suffix, 
                                                    exposure_key=exposure_key, 
                                                    lat_key=lat_key, 
                                                    lon_key=lon_key, 
                                                    elev_key=elev_key, 
                                                    name_key=name_key,
                                                    save_loc=save_loc)

gb_final_transforms.pprint_all()