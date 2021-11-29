# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 12:32:58 2021

@author: jack.wawrow
"""

import AstroFunctions as astro
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings('default', category=UserWarning)

directory = r'C:\Users\jack.wawrow\Documents\2021 10 21 - ZWO with C8\2021 10 21 - Pointing Run\corrected_lights'
plot_results = True
save_plots = True
file_suffix = ".fits"
exposure_key = 'EXPTIME'
lat_key = 'OBSGEO-B', 
lon_key = 'OBSGEO-L', 
elev_key = 'OBSGEO-H',
# transforms_file = r"C:\Users\jack.wawrow\Documents\Suffield\2021 10 19\Standards\Test\corrected_lights\Outputs_Warner_FewerStars\_gb_final_transforms.csv"
# gb_final_transforms = ascii.read(transforms_file)
gb_final_transforms = None
save_loc = os.path.join(directory, 'Outputs_SkySurvey')


star_aux_table = astro._sky_survey_calc(directory, 
                                        plot_results=plot_results, 
                                        save_plots=save_plots,
                                        file_suffix=file_suffix, 
                                        exposure_key=exposure_key,
                                        gb_final_transforms=gb_final_transforms,
                                        lat_key=lat_key, 
                                        lon_key=lon_key, 
                                        elev_key=elev_key,
                                        save_loc=save_loc)

theta = star_aux_table['Azimuth']
r = 90 - star_aux_table['Elevation']
z = star_aux_table['BSB']
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(7,7))
ax.set_theta_zero_location("N")
# norm = matplotlib.colors.Normalize(vmin=14, vmax=22.13)
norm = matplotlib.colors.Normalize(vmin=np.percentile(z, 10), vmax=max(z))
m = cm.ScalarMappable(cmap=plt.get_cmap('Greys'), norm=norm)
m.set_array([])
plt.colorbar(m)
# Change contourf in the line below to scatter if you have only 1D theta, r and brightness values
ax.scatter(theta, r, c=z, cmap=plt.get_cmap('Greys'), norm=norm)
rlabels = ax.get_ymajorticklabels()
for label in rlabels:
	label.set_color('black')
plt.show()
plt.close()