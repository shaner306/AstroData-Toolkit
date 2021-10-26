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
from astropy.time import Time
import matplotlib.dates as mdates

survey_file = r"C:\Users\jack.wawrow\Documents\Suffield\2021 10 19\Standards\Test\corrected_lights\Outputs_SkySurvey\auxiliary_table.csv"
# survey_file = r"C:\Users\jack.wawrow\Documents\Suffield\2021 10 19\Sky Survey\Automated Pointing Run 004\Outputs_SkySurvey\auxiliary_table.csv"
star_aux_table = ascii.read(survey_file)

theta = star_aux_table['Azimuth'][star_aux_table['BSB'] > 5]
r = 90 - star_aux_table['Elevation'][star_aux_table['BSB'] > 5]
z = star_aux_table['BSB'][star_aux_table['BSB'] > 5]
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
# plt.savefig(r'C:\Users\jack.wawrow\Documents\Suffield\2021 10 19\Standards\Test\corrected_lights\Outputs_SkySurvey\BSB_plot')
plt.show()
plt.close()

fwhm = star_aux_table['FWHM']
fwhm_sigma = star_aux_table['FWHM_sigma']
times_list = np.array(star_aux_table['Time (JD)'])
times_obj = Time(times_list, format='jd', scale='utc')
times_datetime = times_obj.to_value('datetime')

fig, ax = plt.subplots()
_, _, bars = ax.errorbar(times_datetime, fwhm, yerr=fwhm_sigma, fmt='o', markersize=2, capsize=0, elinewidth=0.75)
[bar.set_alpha(0.3) for bar in bars]
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax.set_ylabel('FWHM')
ax.set_xlabel('Time (UTC)')
plt.title('FWHM v. Time')
plt.show()
plt.close()
