# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 11:06:28 2021

@author: shane
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
from photutils.background import ModeEstimatorBackground
from photutils.background import MedianBackground
import skimage
from scipy import ndimage
from skimage import measure
from skimage import filters


def BackgroundEstimationMulti(fitsdata, sigma_clip, bkgmethod, printval):
   
    sigma_clip = SigmaClip(sigma=2.5)
    bkg = SExtractorBackground(sigma_clip)
    bkg_value1 = bkg.calc_background(fitsdata)
    
    #print(bkg_value)
    bkg = MeanBackground(sigma_clip)
    bkg_value2 = bkg.calc_background(fitsdata)
    
    bkg = MedianBackground(sigma_clip)
    bkg_value3 = bkg.calc_background(fitsdata)
    
    bkg = ModeEstimatorBackground(sigma_clip)
    bkg_value4 = bkg.calc_background(fitsdata)
    
    
    bkg_estimator1 = SExtractorBackground()
    bkg_estimator2 = SExtractorBackground()
    #bkg = Background2D(fitsdata, (2, 2), filter_size=(3,3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2) Closest Approximate to Matlab Result
    bkg = Background2D(fitsdata, (50,50), filter_size=(3,3),sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2)
    bg_rem = fitsdata - bkg.background

    if printval == 1:
        print("Background Solutions")
        print("Sigma Clip: " + str(sigma_clip))
        print("Using: "+ bkgmethod)
        print("---------------------------")
        print("SExtractor Background: " + str(mean(bkg.background)))
        print("SExtractor Background(Filtered): " + str(mean(bkg.background))+"\n "+ "     " + "Box Size: " +"50x50" +"\n "+ "     " + "Filter Size: " +"3x3")
        print("Mean Background: " + str(bkg_value2))
        print("Median Background: " + str(bkg_value3))
        print("Mode Estimator Background: " + str(bkg_value4))
        print("Remaining Background (subtracted): " + str(bg_rem))
        print("Polyfit Background: Not Implemented Yet")
        
    else:
        return
    
    
        
    
    
    
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


def SaturatedStars(imgin):

# =============================================================================
# % extracts information about saturated stars
# %
# %  Saturated-star geometry:
# %
# %               A, B, C, D and E : saturated pixels with the value: 65535
# %
# %       D              Scan direction: top-down first and then left to right
# %       D
# %       D       A       : fisrt detedted saturated pixel and following pixel
# %       B                 on the same column.
# %      BBB      B and D : pixels are found upward from the center.
# %     ABBBBB    C and E : pixels are found downward from the center.
# %     ACCCC
# %      CCC      B and C : normal bixels belonging to the star core
# %       C       D and E : anormal pixels created by sensor bleeding
# %       E
# %       E
# %       E
# % _________________________________________________________________
# %
# % output structure   sat_stars:
# %
# %  sat_stars.nb_stars:   scalar. This is the number of stars found, and the
# %                        number of lines in the structure members
# %  sat_stars.nb_sat_pixels:  scalar. Total number of saturated pixels found
# %
# %  sat_stars.center:   matrix. Each line is a star center. First sample is
# %                      the line number, second sample is the column number.
# %  sat_stars.width:    vector. Each line is a star width. Only one column.
# %  sat_stars.height:   vector. Each line is a star height. Only one column.
# %  sat_stars.star_pixels:  vector. Each line is a star's number of
# %                          saturated pixels. Only one column.
# %
# %
# 
# =============================================================================

    #error(nargchk(1,1,nargin));
    dx=imgin[0].header['NAXIS1']
    dy=imgin[0].header['NAXIS2']
    MAX_PIXEL_VALUE = 65535;

    idds=imgin[0].data
    imgbin = idds >= MAX_PIXEL_VALUE

    
    
    
   
    N = 0; #% number of stars
    s_nbp = []; #% star's number of saturated pixels
    s_centers = []; #% star's centroid
    s_sizes   = []; #% star's height and width
    s_pixels  = []; #% one star's saturated pixels, re-used for each star
    

    for i in range(dx):
        
        for j in range(dy):
            
            if imgbin[j,i]:
                N = N + 1; # new star found!
                jj = j;
                ii = i;
                
               # extract saturated star
                
                # get all A pixels
                imgbin[jj,ii] = 0;
                s_nbp[N,1] = 1
                s_pixels[s_nbp[N,1], s_nbp[N,1]] = [jj, ii]
                if jj+1 <= dy:
                    keepgoing = imgbin(jj+1, ii); 
                else:
                    keepgoing = 0;
                
                    
                while keepgoing:

                    jj = jj + 1;
                    imgbin[jj,ii] = 0;
                    s_nbp[N,1] = s_nbp[N,1] + 1;
                    s_pixels[s_nbp[N,1], s_nbp[N,1]] = [jj, ii]
                    if jj+1 <= dy:
                        keepgoing = imgbin(jj+1, ii); 
                    else:
                        keepgoing = 0;
                   
                j_center = j;
    
              
                nb_col_pixels = s_nbp(N,1);
                while nb_col_pixels>0 & ii<dx:
                    ii = ii + 1;
                    nb_col_pixels = 0;
           
                    jj = j_center;
                    while( (jj>0) & (imgbin(jj,ii)) ):
                        imgbin[jj,ii] = 0;
                        nb_col_pixels = nb_col_pixels + 1;
                        s_nbp[N,1] = s_nbp[N,1] + 1;
                        s_pixels[s_nbp[N,1], s_nbp[N,1]] = [jj, ii]
                        jj = jj - 1;
                    
                
                    jj = j_center + 1;
                    while( (jj<=dy) & (imgbin(jj,ii)) ):
                        imgbin[jj,ii] = 0;
                        nb_col_pixels = nb_col_pixels + 1;
                        s_nbp[N,1] = s_nbp[N,1] + 1;
                        s_pixels[s_nbp[N,1], s_nbp[N,1]] = [jj, ii]
                        jj = jj + 1;
                    
                    
          
                
                s_sizes[N,2] = max(s_pixels[1:s_nbp[N,1],2]) - min(s_pixels[1:s_nbp[N,1],2]) + 1;
               
                s_sizes[N,1] = max(s_pixels[1:s_nbp[N,1],1]) - min(s_pixels[1:s_nbp[N,1],1]) + 1;
                
                s_centers[N,1] = sum(s_pixels[1:s_nbp[N,1],1]/s_nbp[N,1]);
                s_centers[N,2] = sum(s_pixels[1:s_nbp[N,1],2]/s_nbp[N,1]);
      
    return  N, s_nbp, s_centers, s_sizes, s_pixels