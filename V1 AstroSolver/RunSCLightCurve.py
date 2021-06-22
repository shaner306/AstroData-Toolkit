# -*- coding: utf-8 -*-
"""
Run the light curve generation process.
Created on Thu Jun  3 11:16:13 2021

@author: jmwawrow
"""
import AstroFunctions as astro

# directory = 'D:\\Transfer to mac\\2021-03-10 - Calibrated\\Intelsat 10-02 Post Eclipse\\LIGHT\\B_lim'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-03-20 - Calibrated\Intelsat 10-02\Test'
directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Observations\2016-072'

temp_dir = 'tmp'
file_suffix = '.FIT'
ecct_cut=0.75
max_distance_from_sat = 20
size = 20
max_num_nan = 1
plot_results = 0

sat_dict, sats_table, uncertainty_table, sat_fwhm_table = astro._main_sc_lightcurve(directory, 
                                                                          temp_dir=temp_dir, 
                                                                          file_suffix=file_suffix,
                                                                          ecct_cut=ecct_cut,
                                                                          max_distance_from_sat=max_distance_from_sat, 
                                                                          size=size, 
                                                                          max_num_nan=max_num_nan, 
                                                                          plot_results=plot_results)

print(sat_dict)
sat_fwhm_table.pprint_all()
uncertainty_table.pprint_all()
sats_table.pprint_all()

import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
import matplotlib.dates as mdates
times_list = np.array(sats_table['Time (JD)'])
times_obj = Time(times_list, format='jd', scale='utc')
times_datetime = times_obj.to_value('datetime')
fig, ax = plt.subplots()
plt.errorbar(times_datetime, sats_table['NEOSSat'], yerr=uncertainty_table['NEOSSat'], fmt='o', markersize=3, capsize=2)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.gca().invert_yaxis()
plt.ylabel('Instrumental Magnitude')
plt.xlabel('Time (UTC)')
plt.show(block=True)
# plt.pause(3)
plt.close()
# fig, ax = plt.subplots()
# plt.plot(times_datetime, sat_fwhm_table['FWHM (arcsec)'], 'o')
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# plt.ylabel('FWHM (arcsec)')
# plt.xlabel('Time (UTC)')
# plt.show(block=True)
# # plt.pause(3)
# plt.close()