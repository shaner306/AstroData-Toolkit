# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 09:39:59 2021

@author: jack.wawrow
"""

import AstroFunctions as astro
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
import os
import numpy as np

survey_file = r"C:\Users\jack.wawrow\Documents\Suffield\2021 10 19\Standards\No Flats\Plate Solved\corrected_lights\Outputs_SkySurvey\auxiliary_table.csv"
# survey_file = r"C:\Users\jack.wawrow\Documents\Suffield\2021 10 19\Sky Survey\Automated Pointing Run 004\Outputs_SkySurvey\auxiliary_table.csv"
star_aux_table = ascii.read(survey_file)

theta = star_aux_table['Azimuth']
r = 90 - star_aux_table['Elevation']
z = star_aux_table['BSB']
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(7,7))
ax.set_theta_zero_location("N")
custom_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("",
            ["#44EFF2", "#E0ED42", "#EA8730", "#DB3585", "#3182D7", "#2F79C7", "#2C6EB4", "#2966A6", "#265D97",
             "#235387", "#204C7A", "#1D446E", "#1A3D62", "#173555", "#132C47", "#10253B", "#0C1C2C", "#08121C",
             "#050B11", "#010203"])

norm = matplotlib.colors.Normalize(vmin=14, vmax=22.13)
m = cm.ScalarMappable(cmap=custom_cmap, norm=norm)
m.set_array([])
plt.colorbar(m)
# norm = matplotlib.colors.Normalize(vmin=np.percentile(z[~np.isnan(z)], 7), vmax=max(z[~np.isnan(z)]))
# m = cm.ScalarMappable(cmap=plt.get_cmap('Greys'), norm=norm)
# m.set_array([])
# plt.colorbar(m)
# Change contourf in the line below to scatter if you have only 1D theta, r and brightness values
# ax.scatter(theta[~np.isnan(z)], r[~np.isnan(z)], c=z[~np.isnan(z)], cmap=plt.get_cmap('Greys'), norm=norm)
ax.scatter(theta[~np.isnan(z)], r[~np.isnan(z)], c=z[~np.isnan(z)], cmap=custom_cmap, norm=norm)
# ax.tricontourf(theta[~np.isnan(z)], r[~np.isnan(z)], z[~np.isnan(z)], cmap=custom_cmap, norm=norm)
rlabels = ax.get_ymajorticklabels()
for label in rlabels:
	label.set_color('black')
# plt.savefig(r'C:\Users\jack.wawrow\Documents\Suffield\2021 10 19\Standards\No Flats\Plate Solved\corrected_lights\Outputs_SkySurvey\BSB_plot')
plt.show()
plt.close()