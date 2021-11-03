# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 10:40:03 2021

@author: jack.wawrow
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

file = fits.open(r"C:\Users\jack.wawrow\Documents\Suffield\2021 10 21\2021 10 21 - Raw Landolts\Plate Solved\1021_standards_corrected_lights\2021_10_21_2x2_Light_G_0.10C_7.000secs_SA38_00004856.fits")
print(repr(file[0].header))
plt.imshow(file[0].data, cmap='gray', norm=LogNorm())
plt.show(block=True)
plt.close()
file.close()