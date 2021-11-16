# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 10:40:03 2021

@author: jack.wawrow
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

file = fits.open(r"G:\Suffield\2021 10 26 - ZWO with C14\Sky Survey\2021 10 26 - Pointing Run 1\Plate Solved\corrected_lights\2021_10_26_2x2_0.00C_Light_5.000secs_G_Target_36_00000349.fits")
print(repr(file[0].header))
# plt.imshow(file[0].data, cmap='gray', norm=LogNorm())
# plt.show(block=True)
# plt.close()
file.close()