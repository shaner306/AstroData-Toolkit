# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 11:47:19 2021

@author: shane
"""
from numpy import mean
import numpy
import scipy
from scipy import ndimage
import numpy
import pandas as pd
#import win32com.client as win32
#import win32com
import os
#import pywin32_system32
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime
import astropy
from astropy.io import fits
import PIL
import cv2
from astropy.stats import SigmaClip
from photutils.background import SExtractorBackground
from photutils.background import MeanBackground
from photutils.background import Background2D
import skimage
from scipy import ndimage
from skimage import measure
from skimage import filters

def BackgroundIteration(image, tolerance):       
    old_mean = 1e9;
    old_rms   = 1e9;
    
    new_mean = 2e9;
    new_rms = 2e9;
    
    while abs(new_rms - old_rms) > (tolerance * old_rms):
        old_mean = float(new_mean)
        old_rms = float(new_rms)
        image = myclip(image, (old_mean - 2 * old_rms), (old_mean + 2 * old_rms))
    	
        if (np.size(image) == 0):
            new_mean = 0;
            new_rms = 2e9;
            break;
      
            
       	new_mean = mean(image);
       	new_rms   = numpy.std(image);
        retval = [new_mean, new_rms];
       	return new_mean, new_rms
    
    
def myclip(x1,lo,hi):
        vector = np.vectorize(np.float)
        x= vector(x1)
        
        float(hi)
        float(lo)
        print(x)
        print(hi)
        print(lo)

       	y = (x * np.any(x<= hi)) + ((hi) * np.any(x> hi))
       	y = (y * np.any(y>lo)) + ((lo) * np.any(y<=lo))
        return y



streak = 'D:\\Wawrow\\2. Observational Data\\2021-02-07 - Calibrated\\Intelsat 10-02\\LIGHT\\G\\0257_3x3_-10.00_5.00_G_20-42-55.fits'
imagehdularray = fits.open(streak)

streak_array = [];         
sigma_clip = 2.5;           
edge_protect = 10;          
min_obj_pixels = 5;
SNRLimit = 2.5;

date=imagehdularray[0].header['DATE-OBS']
exposuretime=imagehdularray[0].header['EXPTIME']
imagesizeX=imagehdularray[0].header['NAXIS1']
imagesizeY=imagehdularray[0].header['NAXIS2']
fitsdata =  imagehdularray[0].data
sigma_clip = SigmaClip(sigma=2.5)
bkg = SExtractorBackground(sigma_clip)
bkg_value = bkg.calc_background(fitsdata)

print(bkg_value)
bkg = MeanBackground(sigma_clip)
bkg_value = bkg.calc_background(fitsdata)
bkg_estimator1 = SExtractorBackground()
bkg_estimator2 = SExtractorBackground()
#bkg = Background2D(fitsdata, (2, 2), filter_size=(3,3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2) Closest Approximate to Matlab Result
bkg = Background2D(fitsdata, (10, 10), filter_size=(3,3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2)
bg_rem = fitsdata - bkg.background
print(bkg_value)
print(mean(bg_rem))

bg_rem[1:edge_protect,1:edge_protect] = 0;
bg_rem[imagesizeX - edge_protect:imagesizeX, :] = 0
bg_rem[:, 1:edge_protect] = 0
bg_rem[:, imagesizeY - edge_protect:imagesizeY] = 0
im_mean = mean(bg_rem)

im_rms=np.std(fitsdata)
im_mean, im_rms = BackgroundIteration(bg_rem, 0.1);
low_clip = im_mean + 2.5 * im_rms;
high_clip = 161

binary_image = np.zeros((1224,1832))

bg_rem[bg_rem<= low_clip]
binary_image = (binary_image * bg_rem[bg_rem<= low_clip]) + (1 * bg_rem[bg_rem> low_clip])
th, im_th = cv2.threshold(bg_rem, low_clip, 1, cv2.THRESH_BINARY)
print(im_mean)
blobs_labels = measure.label(im_th, background=0)

# im = cv2.imread('data/src/lena_square_half.png')
# th, im_th = cv2.threshold(im, 128, 255, cv2.THRESH_BINARY)

num_labels, labels_im = cv2.connectedComponents(im_th)
#num_sourcepix = cv2.connectedComponentsWithStats(binary_image, np.array(connected_image), np.array(stats), np.array(centroids), 4, np.int(CV_32S))
