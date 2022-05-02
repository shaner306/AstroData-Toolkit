# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 09:10:52 2021

@author: jack.wawrow
"""

import numpy as np
import os
import AstroFunctions as astro

flatpath = r'G:\Suffield\2021 10 26 - ZWO with C14\2021 10 26 - Flats - 2x2'
file_suffix = ".fits"
for dirpath, dirnames, filenames in os.walk(flatpath):
    for filename in filenames:
        if filename.endswith(file_suffix):
            filepath = os.path.join(dirpath, filename)
            hdr, imgdata = astro.read_fits_file(filepath)
            bkg, _ = astro.calculate_img_bkg(imgdata, 3.0)
            percent = (bkg / (2**hdr['BITPIX'])) * 100
            print(f"{filename}: {bkg} ({percent:.1f} %)")
            # print(f"mean: {mean_value:.3f}, median: {median_value:.3f}")