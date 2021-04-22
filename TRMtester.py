# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 10:33:21 2021
@author: shane

-----
INSTRUCTIONS
1. MUST RUN ON 32 BIT Python - Pinpoint will not run on 64 bit python code
2. Reference Data
3. Output Transforms + Standard Magnitudes
    - Errors
4. Output Star Data
    - FWHM
    - Magnitudes
    - Instrumental Mag
    - Flux, Stellar Flux, Visual Magnitude, gaussian data, sigmas
5. TRM mode
    - Martin Levesque Method
    - Brad Wallace Method
    - Photutils combination
    - SExtractor (background Extraction method)- Currently using
5a. Background Extraction
    -Sextractor
    -Polynomial Form Fit
    -Mean Background
    *Filter and Box Size TBD

6. Creating Stars file and outputing solutions
7. Creating Light Curves

"""
from numpy import mean
import numpy
import scipy
from scipy import ndimage
import numpy
import pandas as pd

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
        #image = myclip(image, (old_mean - 2 * old_rms), (old_mean + 2 * old_rms))
    	
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
        #print(x)
        #print(hi)
        #print(lo)

       	y = (x * np.any(x<= hi)) + ((hi) * np.any(x> hi))
       	y = (y * np.any(y>lo)) + ((lo) * np.any(y<=lo))
        return y
def PointSourceFluxExtraction(mask_x, mask_y, flux_image):
      
    num_elem_x = mask_x.size;
    num_elem_y = mask_y.size;
    sum1 = 0;
    pix_flux = np.zeros((num_elem_x))    
    for i in range(num_elem_x):
        x_pix = mask_x[i];
        y_pix = mask_y[i];
        
        sum1 = sum1 + flux_image[x_pix, y_pix]
        pix_flux[i] = flux_image[x_pix, y_pix]
    
    object_flux = sum1;
    max_pixel_flux = max(pix_flux);
    return object_flux, max_pixel_flux
def MomentCalculation(xmask, ymask, xc, yc, p, q):
        num_pix = xmask.size;
        mom = sum((xmask - xc)**p * (ymask - yc)**q) / num_pix;
        moment=mom
        return moment

def EccentricityCalculation(m11, m02, m20):
    eccent = numpy.sqrt((m20 - m02)**2 + (4*m11**2))/ (m20+m02);
    return eccent

def Compact(num_pix, m02, m20):
    compact = (num_pix/(m02 + m20))
    return compact

def WeightedCentroid(mask_x, mask_y, flux_image):
  
    num_elem_x = mask_x.size;
    num_elem_y = mask_y.size;
    x_wt_sum = 0;
    y_wt_sum = 0;
    flux_sum = 0;
    #print("2")
    if num_elem_x != num_elem_y:
        object_flux = -999;
        #print("3")
        return;
    else:
      for i in range(num_elem_x):
                       
                x_pix = mask_x[i];
                y_pix = mask_y[i];
                
                x_wt_sum =x_wt_sum + (x_pix * flux_image[x_pix, y_pix]);
                y_wt_sum = y_wt_sum + (y_pix * flux_image[x_pix, y_pix]);
                flux_sum = flux_sum + flux_image[x_pix, y_pix];
            
    x_centroid = x_wt_sum / flux_sum;
    y_centroid = y_wt_sum / flux_sum;
    
    
    x_var_sum = 0;
    y_var_sum = 0;
    flux_sum = 0;
    #print("2")
    for i in range(num_elem_x):
                   
            x_pix = mask_x[i];
            y_pix = mask_y[i];
            
            x_var_sum =x_var_sum + ((x_pix-x_centroid)**2 * flux_image[x_pix, y_pix]);
            y_var_sum = y_var_sum + ((y_pix-y_centroid)**2 * flux_image[x_pix, y_pix]);
            flux_sum = flux_sum + flux_image[x_pix, y_pix];
            
    x_rms = numpy.sqrt(x_var_sum/flux_sum);
    y_rms = numpy.sqrt(y_var_sum/flux_sum);
    return x_centroid, x_rms, y_centroid, y_rms

streak = 'D:\\Wawrow\\2. Observational Data\\2021-02-07 - Calibrated\\Intelsat 10-02\\LIGHT\\G\\0150_3x3_-10.00_5.00_G_19-43-13.fits'
imagehdularray = fits.open(streak)

streak_array = [];         
sigma_clip = 3.5;           
edge_protect = 10;          
min_obj_pixels = 5;
SNRLimit = 0;

date=imagehdularray[0].header['DATE-OBS']
exposuretime=imagehdularray[0].header['EXPTIME']
imagesizeX=imagehdularray[0].header['NAXIS1']
imagesizeY=imagehdularray[0].header['NAXIS2']
fitsdata =  imagehdularray[0].data
sigma_clip = SigmaClip(sigma=2.5)
bkg = SExtractorBackground(sigma_clip)
bkg_value = bkg.calc_background(fitsdata)

#print(bkg_value)
bkg = MeanBackground(sigma_clip)
bkg_value = bkg.calc_background(fitsdata)
bkg_estimator1 = SExtractorBackground()
bkg_estimator2 = SExtractorBackground()
#bkg = Background2D(fitsdata, (2, 2), filter_size=(3,3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2) Closest Approximate to Matlab Result
bkg = Background2D(fitsdata, (50,50), filter_size=(3,3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2)
bg_rem = fitsdata - bkg.background
#print(bkg_value)
#print(mean(bg_rem))

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
#binary_image = (binary_image * bg_rem[bg_rem<= low_clip]) + (1 * bg_rem[bg_rem> low_clip])
th, im_th = cv2.threshold(bg_rem, low_clip, 1, cv2.THRESH_BINARY)
#print(im_mean)
connected_image = measure.label(im_th, background=0)
# plt.subplot(133)
# plt.imshow(connected_image, cmap='nipy_spectral')
# plt.axis('off')
# plt.tight_layout()
# plt.show()
#im = cv2.imread(bg_rem)
# th, im_th = cv2.threshold(im, 128, 255, cv2.THRESH_BINARY)

#num_labels, labels_im = cv2.connectedComponents(im_th)
#num_sourcepix = cv2.connectedComponentsWithStats(binary_image, np.array(connected_image), np.array(stats), np.array(centroids), 4, np.int(CV_32S))
num_sourcepix =numpy.zeros(shape=(100000,1))
[size_x, size_y] = 1224,1832
            
for x in range(size_x):
    for y in range(size_y):
        pixval = connected_image[x,y]
        
        if (pixval != 0):
            num_sourcepix[pixval, 0] = num_sourcepix[pixval, 0] + 1;
 
[valid_sources, temp] = numpy.nonzero(num_sourcepix > min_obj_pixels)
num_valid_sources = valid_sources.size

centroid_x = np.zeros((num_valid_sources,1));
centroid_y = np.zeros((num_valid_sources,1));
rms_x_pos = np.zeros((num_valid_sources,1));
rms_y_pos = np.zeros((num_valid_sources,1));
m11        = np.zeros((num_valid_sources,1));
m02        = np.zeros((num_valid_sources,1));
m20        = np.zeros((num_valid_sources,1));
ecct       = np.zeros((num_valid_sources,1));
compact    = np.zeros((num_valid_sources,1));
obj_flux   = np.zeros((num_valid_sources,1));
obj_max1   = np.zeros((num_valid_sources,1));
length     = np.zeros((num_valid_sources,1));


for j in range(num_valid_sources):

    vsj = valid_sources[j];  

    [mask_x, mask_y] = numpy.nonzero(connected_image == vsj);
    obj_flux[j], obj_max1[j] = PointSourceFluxExtraction(mask_x, mask_y, bg_rem);
    
  
    centroid_x[j] = mean(mask_x);
    rms_x_pos[j] = numpy.std(mask_x);
    centroid_y[j] = mean(mask_y);
    rms_y_pos[j] = numpy.std(mask_y);

    m11[j] = MomentCalculation(mask_x, mask_y, centroid_x[j], centroid_y[j], 1,1);
    m02[j] = MomentCalculation(mask_x, mask_y, centroid_x[j], centroid_y[j], 0,2);
    m20[j] = MomentCalculation(mask_x, mask_y, centroid_x[j], centroid_y[j], 2,0);
    compact[j] = Compact(num_sourcepix[vsj], m02[j], m20[j]);
    ecct[j] = EccentricityCalculation(m11[j], m02[j], m20[j]);

    x_length = (max(mask_x) - min(mask_x));
    y_length = (max(mask_y) - min(mask_y));
    length[j] = numpy.sqrt(x_length**2 + y_length**2);

            

[compact_mean, compact_rms] = BackgroundIteration(compact,0.1);
[ecct_mean, ecct_rms] = BackgroundIteration(ecct,0.1)
compact_cut = compact_mean  + 1 * compact_rms;  
ecct_cut = 0.7; 
stars = numpy.nonzero(ecct < ecct_cut);
streaks = numpy.nonzero(ecct > ecct_cut);
stars= np.delete(stars, 1,0)
streaks= np.delete(streaks, 1,0)
sda = valid_sources[stars]
num_pix_in_stars = num_sourcepix[sda]
[mean_starpix, rms_starpix] = BackgroundIteration(num_pix_in_stars, 0.1);

pix_cutoff = mean_starpix + 10 * rms_starpix;

num_stars = stars.size;

stellar_flux_SNR = np.zeros((num_valid_sources,1));
            
xmin = edge_protect;
xmax = imagesizeX - edge_protect;
ymin = edge_protect;
ymax = imagesizeY - edge_protect;
streaksize = streaks.size

for k in range(streaksize):
    
    real_star_num = streaks[0,k]  
    vsj = valid_sources[real_star_num]   
    [mask_x, mask_y] = numpy.nonzero(connected_image == vsj)

    [cen_x, rms_x, cen_y, rms_y] = WeightedCentroid(mask_x, mask_y, bg_rem)
    
    temp_SNR = obj_max1[real_star_num]/im_rms;
    stellar_flux_SNR[k] = temp_SNR;
    #print(temp_SNR)
    if temp_SNR > SNRLimit:
        if ((cen_x > xmin) & (cen_x < xmax) & (cen_y > ymin) & (cen_y < ymax)):
            stellar_flux_SNR[k] = temp_SNR;

            if streak_array== []:
                streak_arrayelement = [cen_x, rms_x, cen_y, rms_y, obj_flux[real_star_num], stellar_flux_SNR[k], exposuretime];
                streak_array.append(streak_arrayelement)
                #print(STARS, '%5.4f %5.4f 10 10 100 %5.0f 0 0.00\n',cen_y, cen_x, obj_flux(rsn));
               # print(Streaks_Detected, [num2str(cen_x) ',' num2str(rms_x) ',' num2str(cen_y) ',' num2str(rms_y) ',' num2str(obj_max1(1,rsn)) ',' num2str(temp_SNR) ',' fpath1(i).name '\r\n']);
            else:
                new_element = [cen_x, rms_x, cen_y, rms_y, obj_flux[real_star_num], stellar_flux_SNR[k], exposuretime];
                streak_array.append(new_element)
               
