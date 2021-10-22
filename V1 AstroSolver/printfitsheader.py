# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 10:40:03 2021

@author: jack.wawrow
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

file = fits.open(r"C:\Users\jack.wawrow\Documents\Suffield\2021 10 19\Standards\No Flats\Plate Solved\corrected_lights\20211019_00003389_Light_G_3x3_SA113_SF1_9.60C.fits")
print(repr(file[0].header))
plt.imshow(file[0].data, cmap='gray', norm=LogNorm())
plt.show(block=True)
plt.close()
file.close()