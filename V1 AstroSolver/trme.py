# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 10:19:36 2021

@author: shane
"""

streak1 = r'D:\Transfer to mac\2021-03-10 - Calibrated\Intelsat 10-02 Post Eclipse\LIGHT\B_lim\0066_3x3_-10.00_5.00_B_21-23-04.fits'
streak = 'D:\\Breeze-M_R_B_38746U\\CAN_OTT.00018670.BREEZE-M_R_B_#38746U.FIT'
STARS = open("CAN_OTT.00018670.BREEZE-M_R_B_#38746U.FIT.stars", "w")


streak_array = []         
sigma_clip = 2.5           
edge_protect = 10          
min_obj_pixels = 5
SNRLimit = 0


imagehdularray = fits.open(streak1)
date=imagehdularray[0].header['DATE-OBS']
exposuretime=imagehdularray[0].header['EXPTIME']
imagesizeX=imagehdularray[0].header['NAXIS1']
imagesizeY=imagehdularray[0].header['NAXIS2']
fitsdata =  imagehdularray[0].data


'''Background Modification'''
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

'''Image with Background noise removed'''
bg_rem = fitsdata - bkg.background

print(mean(bkg.background))

bg_rem[1:edge_protect,1:edge_protect] = 0
bg_rem[imagesizeX - edge_protect:imagesizeX, :] = 0
bg_rem[:, 1:edge_protect] = 0
bg_rem[:, imagesizeY - edge_protect:imagesizeY] = 0
im_mean = mean(bg_rem)




im_rms=np.std(fitsdata)
im_mean, im_rms = BackgroundIteration(bg_rem, 0.1)
low_clip = im_mean + 2.5 * im_rms
high_clip = 161

binary_image = np.zeros((imagesizeX,imagesizeY))

bg_rem[bg_rem<= low_clip]
th, im_th = cv2.threshold(bg_rem, low_clip, 1, cv2.THRESH_BINARY)
#print(im_mean)
connected_image = measure.label(im_th, background=0)

num_sourcepix =numpy.zeros(shape=(100000,1))
[size_x, size_y] = imagesizeX-1,imagesizeY-1
            
for x in range(1,size_y-1):
    for y in range(1,size_x-1):
        pixval = connected_image[x-1,y-1]
        
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

            

[compact_mean, compact_rms] = BackgroundIteration(compact,0.1)
[ecct_mean, ecct_rms] = BackgroundIteration(ecct,0.1)
compact_cut = compact_mean  + 1 * compact_rms  
ecct_cut = 0.7



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