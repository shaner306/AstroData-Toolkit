# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 10:40:03 2021

@author: jack.wawrow
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

file = fits.open(r"C:\Users\jack.wawrow\Documents\2021 10 26 - Automated Pointing Run\corrected_lights\2021_10_26_2x2_Light_V_0.00C_5.000secs_Target_4_00005368.fits")
print(repr(file[0].header))
# plt.imshow(file[0].data, cmap='gray', norm=LogNorm())
# plt.show(block=True)
# plt.close()
file.close()