# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 10:40:03 2021

@author: jack.wawrow
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

file = fits.open(r"D:\Intelsat 10-02\2021-09-17\master_frame_data\master_flat_filter_G.fits")
print(repr(file[0].header))
plt.imshow(file[0].data, cmap='gray', norm=LogNorm())
plt.show(block=True)
plt.close()
file.close()