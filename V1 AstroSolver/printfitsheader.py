# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 10:40:03 2021

@author: jack.wawrow
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

file = fits.open(r"C:\Users\jack.wawrow\Documents\Suffield\2021 10 21\2021 10 21 - Automated Pointing Run\2x2\2021 10 21 - FLAT\2021_10_21_2x2_FlatField_G_-0.30C_1.000secs_00004426.fits")
print(repr(file[0].header))
plt.imshow(file[0].data, cmap='gray')#, norm=LogNorm())
plt.show(block=True)
plt.close()
file.close()