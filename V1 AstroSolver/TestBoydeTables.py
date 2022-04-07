# %% Import Section

#-*- coding: utf-8 -*-
"""
Created on Thu Apr  7 12:14:23 2022

@author: stewe
"""
from ccdproc import ImageFileCollection
import csv
import os
from astropy.nddata import CCDData
import astropy.units as u
import numpy as np
from astropy import table
from astropy.table import Table, QTable, hstack
import matplotlib.pyplot as plt


general_folder_location = r'D:\School\Work - Winter 2022\Work\2022-03-16\2022-03-16\SA23\LIGHT'

# %% Get Boyd Table
Boyde_Table_location = os.path.join(
    general_folder_location, 'Outputs\\Boyde_Table.csv')

file = open(Boyde_Table_location)
csvreader = csv.reader(file)
header = []

header = next(csvreader)
Boyde_Table = Table(names=header, dtype=('str', 'float64', 'float64',
                    'str', 'float64', 'str', 'float64', 'float64', 'float64', 'float64'))
boyde_rows = []
for row in csvreader:

    Boyde_Table.add_row(row)
file.close()


# %% Get File Collection of the images
file_paths = []
for dirpath, dirnames, files in os.walk(general_folder_location):
    for file in files:
        if file.endswith('.fits'):
            file_paths.append(os.path.join(str(dirpath), str(file)))
all_fits = ImageFileCollection(filenames=file_paths)

# %% Load Star Aux Table

aux_file_path = r'D:/School/Work - Winter 2022/Work/2022-03-16/2022-03-16/SA23/Outputs/auxiliary_table.csv'
file = open(aux_file_path)
csvreader = csv.reader(file)
header = []
header = next(csvreader)
Aux_Table = Table(names=header, dtype=['str', 'float64', 'str', 'float64', 'float64', 'float64',
                  'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64' ])
for row in csvreader:
    Aux_Table.add_row(row)
file.close()

# %% Plot Airmass


plot_el = []
plot_airmass = []
plot_k_prime = []
plot_fits_predicted_airmass = []
plot_el_predicted = []
for image in Boyde_Table:
    plot_el.append(float(Aux_Table['Elevation'][[np.where(Aux_Table['filename'] == image['Image Name'])]]))
    plot_airmass.append(float(Aux_Table['X'][[np.where(Aux_Table['filename'] == image['Image Name'])]]))
# =============================================================================
#     plot_airmass.append(float(image['Average Airmass']))
# =============================================================================
    for file_path in file_paths:
        if os.path.basename(file_path)==image['Image Name']:
            plot_fits_predicted_airmass.append(1/np.cos(np.deg2rad(90-(CCDData.read(file_path, unit='adu').header['CENTALT']))))
            plot_el_predicted.append(
                             float(CCDData.read(file_path, unit='adu').header['CENTALT']))
            
# =============================================================================
#     plot_el.append(float((star_aux_table['Elevation'][np.where(
#         star_aux_table['filename'] == image['Image Name'])]).value))
#     plot_airmass.append(float(image['Average Airmass']))
#     plot_k_prime.append(float(image['k_prime']))
#     for file_path in file_paths:
#         if os.path.basename(file_path) == image["Image Name"]:
#             plot_fits_predicted_airmass.append(
#                 1/np.cos(np.deg2rad(90-(CCDData.read(file_path, unit='adu').header['CENTALT']))))
#             plot_el_predicted.append(
#                 float(CCDData.read(file_path, unit='adu').header['CENTALT']))
# plt.figure()
# =============================================================================
plt.plot(plot_el, plot_airmass, 'bo', fillstyle='none',
         label='Boyde Table ')
plt.ylabel('Airmass')
plt.xlim([42.7,43.4])
plt.xlabel('Elevation Angle')
plt.title('Comparing Predicted Centrepoint vs. Averaged Matched Star Airmass ')
plt.plot(plot_el_predicted, plot_fits_predicted_airmass, 'go',
         fillstyle='none', label='Predicted Centrepoint Airmass from FITS Elevation')

# plt.plot(plot_el,plot_k_prime,'ro',label='k_prime')
plt.legend()
plt.show()


