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
import pandas as pd
import win32com
import os
import pywin32_system32
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
from photutils.datasets import make_100gaussians_image
from photutils.segmentation import detect_threshold
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft
from astropy.stats import gaussian_fwhm_to_sigma
from photutils.segmentation import detect_sources
import skimage
from scipy import ndimage
from skimage import measure
from skimage import filters
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.segmentation import SourceCatalog
from photutils.segmentation import deblend_sources
from astropy.visualization import simple_norm
import sep



def BackgroundIteration(image, tolerance):       
    old_mean = 1e9
    old_rms   = 1e9
    
    new_mean = 2e9
    new_rms = 2e9
    
    while abs(new_rms - old_rms) > (tolerance * old_rms):
        old_mean = float(new_mean)
        old_rms = float(new_rms)
        #image = myclip(image, (old_mean - 2 * old_rms), (old_mean + 2 * old_rms))
    	
        if (np.size(image) == 0):
            new_mean = 0
            new_rms = 2e9
            break
      
            
       	new_mean = mean(image)
       	new_rms   = numpy.std(image)
        retval = [new_mean, new_rms]
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
      
    num_elem_x = mask_x.size
    num_elem_y = mask_y.size
    sum1 = 0
    pix_flux = np.zeros((num_elem_x))    
    for i in range(num_elem_x):
        x_pix = mask_x[i]
        y_pix = mask_y[i]
        
        sum1 = sum1 + flux_image[x_pix, y_pix]
        pix_flux[i] = flux_image[x_pix, y_pix]
    
    object_flux = sum1
    max_pixel_flux = max(pix_flux)
    return object_flux, max_pixel_flux
def MomentCalculation(xmask, ymask, xc, yc, p, q):
        num_pix = xmask.size
        mom = sum((xmask - xc)**p * (ymask - yc)**q) / num_pix
        moment=mom
        return moment

def EccentricityCalculation(m11, m02, m20):
    eccent = numpy.sqrt((m20 - m02)**2 + (4*m11**2))/ (m20+m02)
    return eccent

def Compact(num_pix, m02, m20):
    compact = (num_pix/(m02 + m20))
    return compact

def WeightedCentroid(mask_x, mask_y, flux_image):
  
    num_elem_x = mask_x.size
    num_elem_y = mask_y.size
    x_wt_sum = 0
    y_wt_sum = 0
    flux_sum = 0
    #print("2")
    if num_elem_x != num_elem_y:
        object_flux = -999
        #print("3")
        return
    else:
      for i in range(num_elem_x):
                       
                x_pix = mask_x[i]
                y_pix = mask_y[i]
                
                x_wt_sum =x_wt_sum + (x_pix * flux_image[x_pix, y_pix])
                y_wt_sum = y_wt_sum + (y_pix * flux_image[x_pix, y_pix])
                flux_sum = flux_sum + flux_image[x_pix, y_pix]
            
    x_centroid = x_wt_sum / flux_sum
    y_centroid = y_wt_sum / flux_sum
    
    
    x_var_sum = 0
    y_var_sum = 0
    flux_sum = 0
    #print("2")
    for i in range(num_elem_x):
                   
            x_pix = mask_x[i]
            y_pix = mask_y[i]
            
            x_var_sum =x_var_sum + ((x_pix-x_centroid)**2 * flux_image[x_pix, y_pix])
            y_var_sum = y_var_sum + ((y_pix-y_centroid)**2 * flux_image[x_pix, y_pix])
            flux_sum = flux_sum + flux_image[x_pix, y_pix]
            
    x_rms = numpy.sqrt(x_var_sum/flux_sum)
    y_rms = numpy.sqrt(y_var_sum/flux_sum)
    return x_centroid, x_rms, y_centroid, y_rms

streak122 = r'D:\2021_J172\June 22 2021'
streak1 = 'D:\\Breeze-M_R_B_38746U'
streak= r'D:\trm-stars-images\new folder'
streak13 = r'D:\Solved Stars\Tycho 3023_1724\LIGHT\B\0000_3x3_-10.00_5.00_B_21-22-59.fits'

file_suffix = (".fits", ".fit", ".fts",".FIT")

for dirpath, dirnames, filenames in os.walk(streak):
    for filename in filenames:
        print(filename)
        if (filename.endswith(file_suffix)):
            #print(file)
            filepath = os.path.join(dirpath, filename)
            
            STARS = open(filepath+'.stars', "w")
            imagehdularray = fits.open(filepath)
            
            streak_array = []         
            sigma_clip = 2.5           
            edge_protect = 10          
            min_obj_pixels = 5
            SNRLimit = 0
            pix_frac = 0;
            moffat_avg = 0;
            gauss_avg = 0;
            star_count = 0;
            mstar_count = 0;
            count =0
            
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
            threshold = detect_threshold(fitsdata, nsigma=2)
            sigma_clip = SigmaClip(sigma=2.5)
            bkg = SExtractorBackground(sigma_clip)
            bkg_value = bkg.calc_background(fitsdata)
            #fitsdata=fitsdata-bkg_value
            #fitsdata=bg_rem
            
            sigma = 2.5 * gaussian_fwhm_to_sigma  # FWHM = 3.
            kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
            convolved_data = convolve_fft(fitsdata, kernel, normalize_kernel=True)
            #segm = detect_sources(convolved_data, threshold, npixels=2)
            segm = detect_sources(fitsdata, threshold, npixels=5)
            segm_deblend = deblend_sources(convolved_data, segm, npixels=5,
                                           nlevels=32, contrast=0.001)
            
            # norm = ImageNormalize(stretch=SqrtStretch())
            
            
            # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
            # ax1.imshow(fitsdata, origin='lower', cmap='Greys_r', norm=norm)
            # ax1.set_title('Data')
            # cmap = segm.make_cmap(seed=123)
            # ax2.imshow(segm, origin='lower', cmap=cmap, interpolation='nearest')
            # ax2.set_title('Segmentation Image')
            
            cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=convolved_data)
            tbl = cat.to_table()
            tbl['xcentroid'].info.format = '.2f'  # optional format
            tbl['ycentroid'].info.format = '.2f'
            tbl['kron_flux'].info.format = '.2f'
            #print(tbl)
            
            # norm = simple_norm(fitsdata, 'sqrt')
            # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
            # ax1.imshow(fitsdata, origin='lower', cmap='Greys_r', norm=norm)
            # ax1.set_title('FitsData')
            # cmap = segm_deblend.make_cmap(seed=123)
            # ax2.imshow(segm_deblend, origin='lower', cmap=cmap,
            #            interpolation='nearest')
            # ax2.set_title('Segmentation Image')
            # cat.plot_kron_apertures((2.5, 1.0), axes=ax1, color='white', lw=1.5)
            # cat.plot_kron_apertures((2.5, 1.0), axes=ax2, color='white', lw=1.5)
            
            newCenY=list(tbl['ycentroid'])
            newCenX=list(tbl['xcentroid'])
            newflux=list(tbl['segment_flux'])

            for i in range(len(newCenY)):
                
                streak_line='{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(newCenX[i]), float(newCenY[i]),  newflux[i])
                STARS.write(streak_line+"\n")

            STARS.close()
            
            bg_rem[1:edge_protect,1:edge_protect] = 0
            bg_rem[imagesizeX - edge_protect:imagesizeX, :] = 0
            bg_rem[:, 1:edge_protect] = 0
            bg_rem[:, imagesizeY - edge_protect:imagesizeY] = 0
            im_mean = mean(bg_rem)
            
            im_rms=np.std(fitsdata)
            im_mean, im_rms = BackgroundIteration(bg_rem, 0.1)
            low_clip = 110
            high_clip = 161
            
            binary_image = np.zeros((imagesizeX,imagesizeY))
            
            bg_rem[bg_rem<= low_clip]
            #binary_image = (binary_image * bg_rem[bg_rem<= low_clip]) + (1 * bg_rem[bg_rem> low_clip])
            
            
            
            th, im_th = cv2.threshold(bg_rem, low_clip, 1, cv2.THRESH_BINARY)
            #print(im_mean)
            connected_image = measure.label(im_th, background=0)
            plt.subplot(133)
            plt.imshow(connected_image, cmap='nipy_spectral')
            plt.axis('off')
            plt.tight_layout()
            #plt.show()
            #im = cv2.imread(bg_rem)
            # th, im_th = cv2.threshold(im, 128, 255, cv2.THRESH_BINARY)
            
            #num_labels, labels_im = cv2.connectedComponents(im_th)
            #num_sourcepix = cv2.connectedComponentsWithStats(binary_image, np.array(connected_image), np.array(stats), np.array(centroids), 4, np.int(CV_32S))
            num_sourcepix =numpy.zeros(shape=(100000,1))
            [size_x, size_y] = imagesizeX,imagesizeY
                        
            for x in range(0,size_y):
                for y in range(0,size_x):
                    pixval = connected_image[x,y]
                    
                    if (pixval != 0):
                        num_sourcepix[pixval, 0] = num_sourcepix[pixval, 0] + 1
             
            [valid_sources, temp] = numpy.nonzero(num_sourcepix > min_obj_pixels)
            num_valid_sources = valid_sources.size
            
            centroid_x = np.zeros((num_valid_sources,1))
            centroid_y = np.zeros((num_valid_sources,1))
            rms_x_pos = np.zeros((num_valid_sources,1))
            rms_y_pos = np.zeros((num_valid_sources,1))
            m11        = np.zeros((num_valid_sources,1))
            m02        = np.zeros((num_valid_sources,1))
            m20        = np.zeros((num_valid_sources,1))
            ecct       = np.zeros((num_valid_sources,1))
            compact    = np.zeros((num_valid_sources,1))
            obj_flux   = np.zeros((num_valid_sources,1))
            obj_max1   = np.zeros((num_valid_sources,1))
            length     = np.zeros((num_valid_sources,1))
            
            
            for j in range(num_valid_sources):
                a = [0,0]
                vsj = valid_sources[j]  
            
                [mask_x, mask_y] = numpy.nonzero(connected_image == vsj)
                obj_flux[j], obj_max1[j] = PointSourceFluxExtraction(mask_x, mask_y, bg_rem)
                
              
                centroid_x[j] = mean(mask_x)
                rms_x_pos[j] = numpy.std(mask_x)
                centroid_y[j] = mean(mask_y)
                rms_y_pos[j] = numpy.std(mask_y)
            
                m11[j] = MomentCalculation(mask_x, mask_y, centroid_x[j], centroid_y[j], 1,1)
                m02[j] = MomentCalculation(mask_x, mask_y, centroid_x[j], centroid_y[j], 0,2)
                m20[j] = MomentCalculation(mask_x, mask_y, centroid_x[j], centroid_y[j], 2,0)
                compact[j] = Compact(num_sourcepix[vsj], m02[j], m20[j])
                ecct[j] = EccentricityCalculation(m11[j], m02[j], m20[j])
            
                x_length = (max(mask_x) - min(mask_x))
                y_length = (max(mask_y) - min(mask_y))
                length[j] = numpy.sqrt(x_length**2 + y_length**2)
                
                
               
            
                Zp = 21
                vmag= Zp - 2.5*np.log10(obj_flux[j]/exposuretime)
            
               
                if obj_max1[j] < 60000:  #%fit unsaturated stars only
                      if (centroid_x[j] > 10) and (centroid_x[j] < (imagesizeY-10)) and (centroid_y[j] > 10) and ( centroid_y[j] < (imagesizeX-10)):
                       #%Find middle pixel value
                        
                       [cen_x, rms_x, cen_y, rms_y] = WeightedCentroid(mask_x, mask_y, 0*bg_rem+1)
            
                       if (centroid_x[j] > 10) and (centroid_x[j] < (imagesizeX-10)) and (centroid_y[j] > 10) and (centroid_y[j] < (imagesizeY-10)):
                           
                            #mid_pix_val = bg_rem(round(cen_x),round(cen_y))
                            cenx=int(centroid_y[j])
                            ceny=int(centroid_x[j])
                            mid_pix_valPP = bg_rem[ceny,cenx]
                            
                            if vmag < 13:
                                #Fit a moffat profile
                                r = np.zeros(len(mask_x))  #holds radial distance from centroid
                                S = np.zeros(len(mask_x))  #holds intensity
                                np.delete(S, -1)
                                for q in range(0,len(mask_x)):
                                   r[q] = np.sqrt((mask_x[q]+0.5-(ceny+1))**2 + (mask_y[q]+0.5-(cenx+1))**2)
                                   S[q] = bg_rem[mask_x[q],mask_y[q]]
                                
               
                                C_index = np.argmin(r)
                                r[C_index] = 0; #%centroid radial value
                                C = S[C_index]
                                #%a holds [alpha Beta] moffat parameters
                               # %Fix a(2) Beta parameter to 1.5
                                #a = [0,0]
                                #print(C)
                                fun = lambda a: sum((S - (C/((1+(r**2)/(a[0]**2))**1.5)))**2)
                                aguess = 1
                                a = scipy.optimize.fmin(func=fun, x0=aguess, disp=0)
                                
                                
                                #%b holds [alpha Beta] moffat parameters
                                
                                fung = lambda b: sum((S - (C*np.exp(-(r**2)/(2*(b**2)))))**2)
                                bguess = 2;
                                b = scipy.optimize.fmin(func=fung, x0=bguess, disp=0)
                                #print(b)
                                #%Optional plot the fits:
                                
                                #plt.scatter(r,S);
                                #E = lambda a,r: (C/((1+(r**2)/(a[0]^2))**1.5))
                                #F = lambda b,r:(C*np.exp(-(r**2)/(2*(b**2))))
                                #plot=plt(E,[0,max(r)])
                                
                                #h = plt.gca().get_children()
                                
                                #plot.set(h(1),'color','red')
                                
                                #plot= plt(F,[0,max(r)])
                                #plt.axis([0,max(r),0,60000])
                                
                                #h = plt.gca().get_children()
                                #plot.set(h(1),'color','green')
                                # Output results
                               
                            else: 
                                a = [0,0]
                                b = 0
                                
                     
                            
                            pix_frac = pix_frac + mid_pix_valPP/obj_flux[j];
                            
                            if vmag < 13 and a[0]<4:
                                #mstar_count = mstar_count +1;
                                #print(a[0])
                                count = count+1
                                moffat_avg = moffat_avg + a[0];
                                gauss_avg = gauss_avg + b;
            
                        
            
            [compact_mean, compact_rms] = BackgroundIteration(compact,0.1)
            [ecct_mean, ecct_rms] = BackgroundIteration(ecct,0.1)
            compact_cut = compact_mean  + 1 * compact_rms  
            ecct_cut = 0.5
            
            stars = numpy.nonzero(ecct < ecct_cut)
            streaks = numpy.nonzero(ecct > ecct_cut)
            stars= np.delete(stars, 1,0)
            streaks= np.delete(streaks, 1,0)
            
            sda = valid_sources[stars]
            num_pix_in_stars = num_sourcepix[sda]
            [mean_starpix, rms_starpix] = BackgroundIteration(num_pix_in_stars, 0.1)
            
            pix_cutoff = mean_starpix + 10 * rms_starpix
            
            num_stars = stars.size
            
            stellar_flux_SNR = np.zeros((num_valid_sources,1))
                        
            xmin = edge_protect
            xmax = imagesizeX - edge_protect
            ymin = edge_protect
            ymax = imagesizeY - edge_protect
            streaksize = streaks.size
            
            
            
            
            
            for k in range(streaksize):
                
                real_star_num = streaks[0,k]  
                vsj = valid_sources[real_star_num]   
                [mask_x, mask_y] = numpy.nonzero(connected_image == vsj)
            
                [cen_x, rms_x, cen_y, rms_y] = WeightedCentroid(mask_x, mask_y, bg_rem)
                
                temp_SNR = obj_max1[real_star_num]/im_rms
                stellar_flux_SNR[k] = temp_SNR
                #print(temp_SNR)
                if temp_SNR > SNRLimit:
                    if ((cen_x > xmin) & (cen_x < xmax) & (cen_y > ymin) & (cen_y < ymax)):
                        stellar_flux_SNR[k] = temp_SNR
            
                        if streak_array== []:
                            streak_arrayelement = [cen_x, rms_x, cen_y, rms_y, obj_flux[real_star_num], stellar_flux_SNR[k], exposuretime]
                            streak_array.append(streak_arrayelement)
                            flux=float(obj_flux[real_star_num,0])
                            streak_line='{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(cen_y), float(cen_x),  flux)
                           
                            #STARS.write(streak_line+"\n")
                           # print(Streaks_Detected, [num2str(cen_x) ',' num2str(rms_x) ',' num2str(cen_y) ',' num2str(rms_y) ',' num2str(obj_max1(1,rsn)) ',' num2str(temp_SNR) ',' fpath1(i).name '\r\n'])
                        else:
                            new_element = [cen_x, rms_x, cen_y, rms_y, obj_flux[real_star_num], stellar_flux_SNR[k], exposuretime]
                            flux=float(obj_flux[real_star_num,0])
                            streak_line='{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(cen_y), float(cen_x),  flux)
                            #STARS.write(streak_line+"\n")
                            streak_array.append(new_element)
            
            #avg_pix_frac = pix_frac/star_count
            #moffat_avg = moffat_avg/count
            #gauss_avg = gauss_avg/count
            #FWHM= 2*gauss_avg*0.7664
            #print("FWMM: " + str(FWHM))
           # STARS.close()







# [bg_mean, bg_rms] = determine_bg_iteratively(enlarged_bg_image, 0.1);
# pix_frac = 0;
# moffat_avg = 0;
# gauss_avg = 0;
# star_count = 0;
# mstar_count = 0;
# #Get the matched star intesity and max pixel     

# nmstars = p.MatchedStars.Count
# mstars = p.MatchedStars;

# for j in nmstars:
#     mstar = mstars.Item(j)
#     rawflux = mstar.RawFlux
#     Zp = p.MagZeroPoint
#     vmag= Zp - 2.5*log10(rawflux/exptime)
#     #StarXY = [mstar.X mstar.Y]
#     InstrumentalMag= 1;
#     #%Get the bg avg and std
#     ppbgsigma = p.ImageBackgroundSigma;
#     ppbgmean = p.ImageBackgroundmean;
#     SQmean = bg_mean;
#     SQsigma = im_rms;
# #           Get the flux and max_ppixel of each star
#     X = round(mstar.X)+1
#     if X > image_size_y:
#         X = image_size_y;
    
#     Y = round(mstar.Y)+1
#     if Y > image_size_x:
#         Y = image_size_x
    
#     st_index = connected_image(Y,X);
#     [mask_x, mask_y] = numpy.nonzero(connected_image == vsj)
#     mobj_flux[j], mobj_max1[j] = PointSourceFluxExtraction(mask_x, mask_y, bg_rem)
   
#     if mobj_max1[j] < 60000:  #%fit unsaturated stars only
# #                Border checks
#           if (X > 10) and (X < (image_size_y-10)) and (Y > 10) and (Y < (image_size_x-10)):
#            #%Find middle pixel value
            
#            [cen_x, rms_x, cen_y, rms_y] = WeightedCentroid(mask_x, mask_y, 0*bg_rem+1)

#            if (cen_x > 10) and (cen_x < (image_size_x-10)) and (cen_y > 10) and (cen_y < (image_size_y-10)):
#                 mid_pix_val = bg_rem(round(cen_x),round(cen_y))
                
#                 mid_pix_valPP = bg_rem(Y,X)
                
#                 if vmag < 13:
#                     #Fit a moffat profile
#                     r = np.zeros(1,mask_x);  #holds radial distance from centroid
#                     S = np.zeros(1,mask_x);  #holds intensity
#                     for q in mask_x:
#                        r[q] = np.sqrt((mask_x(q)+0.5-(mstar.Y+1))^2 + (mask_y(q)+0.5-(mstar.X+1))^2)
#                        S[q] = bg_rem(mask_x(q),mask_y(q))
                    
   
#                     C_index = numpy.nonzero(r==min(r),1)
#                     r[C_index] = 0; #%centroid radial value
#                     C = S[C_index]
#                     #%a holds [alpha Beta] moffat parameters
#                    # %Fix a(2) Beta parameter to 1.5
#                     fun = lambda a: sum((S - (C/((1+(r**2)/(a(1)^2))**1.5)))**2)
#                     aguess = 1
#                     a = scipy.optimize.fmin(func=fun, x0=aguess)
                    
                    
#                     #%b holds [alpha Beta] moffat parameters
                    
#                     fung = lambda b: sum((S - (C*math.exp(-(r**2)/(2*(b^2)))))**2)
#                     bguess = 2;
#                     b = scipy.optimize.fmin(func=fung, x0=bguess)
                    
#                     #%Optional plot the fits:
                    
#                     plt.scatter(r,S);
#                     E = lambda a,r: (C/((1+(r**2)/(a(1)^2))**1.5))
#                     F = lambda b,r:(C*math.exp(-(r**2)/(2*(b^2))))
#                     plot=plt(E,[0,max(r)])
                    
#                     h = plt.gca().get_children()
                    
#                     plot.set(h(1),'color','red')
                    
#                     plot= plt(F,[0,max(r)])
#                     plt.axis([0,max(r),0,60000])
                    
#                     h = plt.gca().get_children()
#                     plot.set(h(1),'color','green')
#                     # Output results
                   
#                 else: 
#                     a = [0,0]
#                     b = 0
                    
         
#                 star_count = star_count +1;
#                 pix_frac = pix_frac + mid_pix_valPP/rawflux;
#                 if vmag < 13:
#                     mstar_count = mstar_count +1;
#                     moffat_avg = moffat_avg + a(1);
#                     gauss_avg = gauss_avg + b;
