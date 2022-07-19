# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 16:03:42 2022

@author: mstew
"""
import ctypes
import os
import tkinter as tk
from collections import namedtuple

import cv2 as cv
import numpy
import numpy as np
import scipy
from astropy import units as u
from astropy.io import ascii
from astropy.stats import SigmaClip
from astropy.table import Table, hstack
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
    NavigationToolbar2Tk
from matplotlib.colors import LogNorm
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from photutils.aperture import RectangularAperture
from photutils.background import Background2D
from photutils.background import SExtractorBackground
from skimage import measure

import sys
from os.path import dirname
src_path = dirname(dirname(dirname(__file__)))
sys.path.append(os.path.join(src_path, 'general_tools'))
sys.path.append(os.path.join(src_path, 'photometry'))
sys.path.append(os.path.join(src_path, 'Streakdetection'))
import AstroFunctions as astro
import perform_photometry
import StreakDetection as sd
# from AstroFunctions import copy_and_rename, read_fits_file, add_new_time_and_filter, change_sat_positions, \
#     calculate_background_sky_brightness, convert_fwhm_to_arcsec_trm, get_image_airmass, calculate_magnitudes, \
#     calculate_magnitudes_sigma, check_if_sat, determine_if_change_sat_positions, remove_temp_dir, determine_num_filters, \
#     interpolate_sats, save_interpolated_light_curve, get_all_indicies_combinations, calculate_timeseries_colour_indices, \
#     choose_indices_to_plot, axis_limits_multiband_gui, apply_gb_timeseries_transforms, choose_aux_data_to_plot, \
#     axis_limits_singleband_gui


def TRM_sat_detection(filepath,
                      ecct_cut=0.5,
                      sigma_clip=2.5,
                      edge_protect=10,
                      SNRLimit=0,
                      min_obj_pixels=5,
                      pix_frac=0,
                      moffat_avg=0,
                      gauss_avg=0,
                      star_count=0,
                      mstar_count=0,
                      count=0):
    def BackgroundIteration(image, tolerance):
        # old_mean = 1e9
        old_rms = 1e9

        new_mean = 2e9
        new_rms = 2e9

        while abs(new_rms - old_rms) > (tolerance * old_rms):
            # old_mean = float(new_mean)
            old_rms = float(new_rms)
            # image = myclip(image, (old_mean - 2 * old_rms),
            # (old_mean + 2 * old_rms))

            if (np.size(image) == 0):
                new_mean = 0
                new_rms = 2e9
                break

            new_mean = np.mean(image)
            new_rms = np.std(image)
            # retval = [new_mean, new_rms]
            return new_mean, new_rms

    def myclip(x1, lo, hi):
        vector = np.vectorize(np.float)
        x = vector(x1)

        float(hi)
        float(lo)
        # print(x)
        # print(hi)
        # print(lo)

        y = (x * np.any(x <= hi)) + ((hi) * np.any(x > hi))
        y = (y * np.any(y > lo)) + ((lo) * np.any(y <= lo))
        return y

    def PointSourceFluxExtraction(mask_x, mask_y, flux_image):

        num_elem_x = mask_x.size
        # num_elem_y = mask_y.size
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
        mom = sum((xmask - xc) ** p * (ymask - yc) ** q) / num_pix
        moment = mom
        return moment

    def EccentricityCalculation(m11, m02, m20):
        eccent = np.sqrt((m20 - m02) ** 2 + (4 * m11 ** 2)) / (m20 + m02)
        return eccent

    def Compact(num_pix, m02, m20):
        compact = (num_pix / (m02 + m20))
        return compact

    def WeightedCentroid(mask_x, mask_y, flux_image):

        num_elem_x = mask_x.size
        num_elem_y = mask_y.size
        x_wt_sum = 0
        y_wt_sum = 0
        flux_sum = 0
        # print("2")
        if num_elem_x != num_elem_y:
            # object_flux = -999
            # print("3")
            return
        else:
            for i in range(num_elem_x):
                x_pix = mask_x[i]
                y_pix = mask_y[i]

                x_wt_sum = x_wt_sum + (x_pix * flux_image[x_pix, y_pix])
                y_wt_sum = y_wt_sum + (y_pix * flux_image[x_pix, y_pix])
                flux_sum = flux_sum + flux_image[x_pix, y_pix]

        x_centroid = x_wt_sum / flux_sum
        y_centroid = y_wt_sum / flux_sum

        x_var_sum = 0
        y_var_sum = 0
        flux_sum = 0
        # print("2")
        for i in range(num_elem_x):
            x_pix = mask_x[i]
            y_pix = mask_y[i]

            x_var_sum = x_var_sum + \
                ((x_pix - x_centroid) ** 2 * flux_image[x_pix, y_pix])
            y_var_sum = y_var_sum + \
                ((y_pix - y_centroid) ** 2 * flux_image[x_pix, y_pix])
            flux_sum = flux_sum + flux_image[x_pix, y_pix]

        x_rms = np.sqrt(x_var_sum / flux_sum)
        y_rms = np.sqrt(y_var_sum / flux_sum)
        return x_centroid, x_rms, y_centroid, y_rms

    # streak1 = r'D:\Transfer to mac\2021-03-10 - Calibrated\Intelsat 10-02
    # Post Eclipse\LIGHT\B_lim\0066_3x3_-10.00_5.00_B_21-23-04.fits'
    # streak = \
    # 'D:\\Breeze-M_R_B_38746U\\CAN_OTT.00018674.BREEZE-M_R_B_#38746U.FIT'
    # streak12 =\
    #  r'D:\Transfer to mac\trm-stars-images\NEOS_SCI_2021099173229frame.fits'
    hdr, fitsdata = astro.read_fits_file(filepath)
    # STARS = open(filepath, "w")
    # imagehdularray = fits.open(filepath)

    streak_array = []
    # sigma_clip = 2.5
    # edge_protect = 10
    # min_obj_pixels = 5
    # SNRLimit = 0
    # pix_frac = 0;
    # moffat_avg = 0;
    # gauss_avg = 0;
    # star_count = 0;
    # mstar_count = 0;
    # count =0

    # date=imagehdularray[0].header['DATE-OBS']
    # exposuretime=imagehdularray[0].header['EXPTIME']
    # imagesizeX=imagehdularray[0].header['NAXIS1']
    # imagesizeY=imagehdularray[0].header['NAXIS2']
    imagesizeX = hdr['NAXIS1']
    imagesizeY = hdr['NAXIS2']
    exposuretime = hdr['EXPTIME']
    # fitsdata =  imagehdularray[0].data
    sigma_clip = SigmaClip(sigma=2.5)
    # bkg = SExtractorBackground(sigma_clip)
    # bkg_value = bkg.calc_background(fitsdata)

    # print(bkg_value)
    # bkg = MeanBackground(sigma_clip)
    # bkg_value = bkg.calc_background(fitsdata)
    # bkg_estimator1 = SExtractorBackground()
    bkg_estimator2 = SExtractorBackground()
    # bkg = Background2D(fitsdata, (2, 2),
    # filter_size=(3,3),sigma_clip=sigma_clip,
    # bkg_estimator=bkg_estimator2) Closest Approximate to Matlab Result
    bkg = Background2D(fitsdata, (50, 50), filter_size=(
        3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator2)
    bg_rem = fitsdata - bkg.background

    # print(np.mean(bkg.background))

    bg_rem[1:edge_protect, 1:edge_protect] = 0
    bg_rem[imagesizeX - edge_protect:imagesizeX, :] = 0
    bg_rem[:, 1:edge_protect] = 0
    bg_rem[:, imagesizeY - edge_protect:imagesizeY] = 0
    im_mean = np.mean(bg_rem)

    im_rms = np.std(fitsdata)
    im_mean, im_rms = BackgroundIteration(bg_rem, 0.1)
    low_clip = im_mean + 2.5 * im_rms
    # high_clip = 161

    # binary_image = np.zeros((imagesizeX,imagesizeY))

    bg_rem[bg_rem <= low_clip]
    # binary_image = (binary_image * bg_rem[bg_rem<= low_clip]) +
    # (1 * bg_rem[bg_rem> low_clip])

    th, im_th = cv.threshold(bg_rem, low_clip, 1, cv.THRESH_BINARY)
    # print(im_mean)
    connected_image = measure.label(im_th, background=0)
    # plt.subplot(133)
    # plt.imshow(connected_image, cmap='nipy_spectral')
    # plt.axis('off')
    # plt.tight_layout()
    # plt.show()
    # im = cv2.imread(bg_rem)
    # th, im_th = cv2.threshold(im, 128, 255, cv2.THRESH_BINARY)

    # num_labels, labels_im = cv2.connectedComponents(im_th)
    # num_sourcepix = cv2.connectedComponentsWithStats(binary_image,
    # np.array(connected_image), np.array(stats), np.array(centroids),
    # 4, np.int(CV_32S))
    num_sourcepix = np.zeros(shape=(100000, 1))
    [size_x, size_y] = imagesizeX, imagesizeY

    for x in range(0, size_y):
        for y in range(0, size_x):
            pixval = connected_image[x, y]

            if (pixval != 0):
                num_sourcepix[pixval, 0] = num_sourcepix[pixval, 0] + 1

    [valid_sources, temp] = np.nonzero(num_sourcepix > min_obj_pixels)
    num_valid_sources = valid_sources.size
    if num_valid_sources == 0:
        return

    centroid_x = np.zeros((num_valid_sources, 1))
    centroid_y = np.zeros((num_valid_sources, 1))
    rms_x_pos = np.zeros((num_valid_sources, 1))
    rms_y_pos = np.zeros((num_valid_sources, 1))
    m11 = np.zeros((num_valid_sources, 1))
    m02 = np.zeros((num_valid_sources, 1))
    m20 = np.zeros((num_valid_sources, 1))
    ecct = np.zeros((num_valid_sources, 1))
    compact = np.zeros((num_valid_sources, 1))
    obj_flux = np.zeros((num_valid_sources, 1))
    obj_max1 = np.zeros((num_valid_sources, 1))
    length = np.zeros((num_valid_sources, 1))

    for j in range(num_valid_sources):

        vsj = valid_sources[j]
        a = [0, 0]

        [mask_x, mask_y] = np.nonzero(connected_image == vsj)
        obj_flux[j], obj_max1[j] = PointSourceFluxExtraction(
            mask_x, mask_y, bg_rem)

        centroid_x[j] = np.mean(mask_x)
        rms_x_pos[j] = np.std(mask_x)
        centroid_y[j] = np.mean(mask_y)
        rms_y_pos[j] = np.std(mask_y)

        m11[j] = MomentCalculation(
            mask_x, mask_y, centroid_x[j], centroid_y[j], 1, 1)
        m02[j] = MomentCalculation(
            mask_x, mask_y, centroid_x[j], centroid_y[j], 0, 2)
        m20[j] = MomentCalculation(
            mask_x, mask_y, centroid_x[j], centroid_y[j], 2, 0)
        compact[j] = Compact(num_sourcepix[vsj], m02[j], m20[j])
        ecct[j] = EccentricityCalculation(m11[j], m02[j], m20[j])

        x_length = (max(mask_x) - min(mask_x))
        y_length = (max(mask_y) - min(mask_y))
        length[j] = np.sqrt(x_length ** 2 + y_length ** 2)

        Zp = 21
        vmag = Zp - 2.5 * np.log10(obj_flux[j] / exposuretime)

        if obj_max1[j] < 60000:  # %fit unsaturated stars only
            if (centroid_x[j] > 10) and (centroid_x[j] < (imagesizeY - 10)) \
                and (centroid_y[j] > 10) and (
                    centroid_y[j] < (imagesizeX - 10)):
                # %Find middle pixel value

                [cen_x, rms_x, cen_y, rms_y] = WeightedCentroid(
                    mask_x, mask_y, 0 * bg_rem + 1)

                if (centroid_x[j] > 10) and\
                    (centroid_x[j] < (imagesizeX - 10)) and\
                        (centroid_y[j] > 10) and\
                        (centroid_y[j] < (imagesizeY - 10)):
                    # mid_pix_val = bg_rem(round(cen_x),round(cen_y))
                    cenx = int(centroid_y[j])
                    ceny = int(centroid_x[j])
                    mid_pix_valPP = bg_rem[ceny, cenx]

                    if vmag < 13:
                        # Fit a moffat profile
                        # holds radial distance from centroid
                        r = np.zeros(len(mask_x))
                        S = np.zeros(len(mask_x))  # holds intensity
                        np.delete(S, -1)
                        for q in range(0, len(mask_x)):
                            r[q] = np.sqrt(
                                (mask_x[q] + 0.5 - (ceny + 1)) ** 2 +
                                (mask_y[q] + 0.5 - (cenx + 1)) ** 2)
                            S[q] = bg_rem[mask_x[q], mask_y[q]]

                        C_index = np.argmin(r)
                        r[C_index] = 0  # %centroid radial value
                        C = S[C_index]
                        # %a holds [alpha Beta] moffat parameters
                        # %Fix a(2) Beta parameter to 1.5
                        # a = [0,0]
                        # print(C)

                        def fun(a): return sum(
                            (S - (C / ((1 + (r ** 2) / (a[0] ** 2)) ** 1.5)))
                            ** 2)
                        aguess = 1
                        a = scipy.optimize.fmin(func=fun, x0=aguess, disp=0)
                        # print(a)

                        # %b holds [alpha Beta] moffat parameters

                        def fung(b): return sum(
                            (S - (C * np.exp(-(r ** 2) / (2 * (b ** 2))))) **
                            2)
                        bguess = 2
                        b = scipy.optimize.fmin(func=fung, x0=bguess, disp=0)
                        # print(b)
                        # %Optional plot the fits:

                        # plt.scatter(r,S);
                        # E = lambda a,r: (C/((1+(r**2)/(a[0]^2))**1.5))
                        # F = lambda b,r:(C*np.exp(-(r**2)/(2*(b**2))))
                        # plot=plt(E,[0,max(r)])

                        # h = plt.gca().get_children()

                        # plot.set(h(1),'color','red')

                        # plot= plt(F,[0,max(r)])
                        # plt.axis([0,max(r),0,60000])

                        # h = plt.gca().get_children()
                        # plot.set(h(1),'color','green')
                        # Output results

                    else:
                        a = [0, 0]
                        b = 0

                    pix_frac = pix_frac + mid_pix_valPP / obj_flux[j]

                    if vmag < 13 and a[0] < 4:
                        # mstar_count = mstar_count +1;
                        # print(a[0])
                        count = count + 1
                        moffat_avg = moffat_avg + a[0]
                        gauss_avg = gauss_avg + b

    # [compact_mean, compact_rms] = BackgroundIteration(compact,0.1)
    [ecct_mean, ecct_rms] = BackgroundIteration(ecct, 0.1)
    # compact_cut = compact_mean  + 1 * compact_rms
    # ecct_cut = 0.5

    stars = np.nonzero(ecct < ecct_cut)
    streaks = np.nonzero(ecct > ecct_cut)
    stars = np.delete(stars, 1, 0)
    streaks = np.delete(streaks, 1, 0)
    sda = valid_sources[stars]
    if len(sda) == 0 or count == 0:
        return
    num_pix_in_stars = num_sourcepix[sda]
    [mean_starpix, rms_starpix] = BackgroundIteration(num_pix_in_stars, 0.1)

    sat_x = centroid_y[stars].flatten()
    sat_y = centroid_x[stars].flatten()
    # sat_array = np.empty((len(stars[0]), 2))
    # for sat_index in range(len(stars[0])):
    #     sat_array[sat_index][0] = float(sat_x[0][sat_index])
    #     sat_array[sat_index][1] = float(sat_y[0][sat_index])

    pix_cutoff = mean_starpix + 10 * rms_starpix

    num_stars = stars.size

    stellar_flux_SNR = np.zeros((num_valid_sources, 1))

    xmin = edge_protect
    xmax = imagesizeX - edge_protect
    ymin = edge_protect
    ymax = imagesizeY - edge_protect
    streaksize = streaks.size

    for k in range(streaksize):

        real_star_num = streaks[0, k]
        vsj = valid_sources[real_star_num]
        [mask_x, mask_y] = np.nonzero(connected_image == vsj)

        [cen_x, rms_x, cen_y, rms_y] = WeightedCentroid(mask_x, mask_y, bg_rem)

        temp_SNR = obj_max1[real_star_num] / im_rms
        stellar_flux_SNR[k] = temp_SNR
        # print(temp_SNR)
        if temp_SNR > SNRLimit:
            if ((cen_x > xmin) & (cen_x < xmax) & (cen_y > ymin) &
                    (cen_y < ymax)):
                stellar_flux_SNR[k] = temp_SNR

                if streak_array == []:
                    streak_arrayelement = [cen_x, rms_x, cen_y, rms_y,
                                           obj_flux[real_star_num],
                                           stellar_flux_SNR[k],
                                           exposuretime]
                    streak_array.append(streak_arrayelement)
                    # flux=float(obj_flux[real_star_num,0])
                    # streak_line='{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(cen_y), float(cen_x),  flux)

                    # STARS.write(streak_line+"\n")
                    # print(Streaks_Detected,
                    # [num2str(cen_x) ','
                    # num2str(rms_x) ','
                    # num2str(cen_y) ','
                    # num2str(rms_y) ','
                    # num2str(obj_max1(1,rsn))
                    # ',' num2str(temp_SNR)
                    # ',' fpath1(i).name '\r\n'])
                else:
                    new_element = [cen_x,
                                   rms_x,
                                   cen_y,
                                   rms_y,
                                   obj_flux[real_star_num],
                                   stellar_flux_SNR[k],
                                   exposuretime]
                    # flux=float(obj_flux[real_star_num,0])
                    # streak_line='{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(cen_y), float(cen_x),  flux)
                    # STARS.write(streak_line+"\n")
                    streak_array.append(new_element)

    # avg_pix_frac = pix_frac/star_count
    moffat_avg = moffat_avg / count
    gauss_avg = gauss_avg / count
    FWHM = float(2 * gauss_avg * 0.7664)

    # STARS.close()
    return sat_x, sat_y, bkg.background, FWHM


def set_sat_positions(imgdata,
                      filecount,
                      set_sat_positions_bool,
                      max_distance_from_sat=25,
                      norm=LogNorm(),
                      cmap_set='Set1'):
    """
    Initialize the names and locations of the satellites to create a light
    curve for.

    Parameters
    ----------
    imgdata : numpy.ndarray
        Data from the fits file.
    filecount : float
        Number of files that will be used to create the light curve.
    set_sat_positions_bool : bool
        Decides whether or not to initialize the satellite positions.
    max_distance_from_sat : int, optional
        Maximimum number of pixels that a source can be away from the defined
        sat position to be considered the sat. 
        The default is 25.
    norm : TYPE, optional
        DESCRIPTION. The default is LogNorm(). #TODO
    cmap_set : string, optional
        CMAP set to use for plotting the satellite positions. The default is
        'Set1'.

    Returns
    -------
    TYPE
        DESCRIPTION.
        TODO

    """

    def mbox(title, text, style):
        return ctypes.windll.user32.MessageBoxW(0, text, title, style)

    def set_sat_position(event, x, y, flags, params):
        global sat_locs
        if event == cv.EVENT_LBUTTONDOWN:
            sat_locs.append([x, y])

    def return_entry(event=None):
        """Get and print the content of the entry."""
        # global entry
        global content
        content = entry.get()
        root.destroy()

    def onclick(event, sat_locs):
        # global num_sat_event
        sat_locs.append([event.xdata, event.ydata])

    # global set_sat_positions_bool
    while set_sat_positions_bool:
        # global sat_locs
        sat_locs = []
        # global num_sat_event
        # num_sat_event = 0
        # return satloc
        # mbox('Information',
        #       'Please select the positions of the satellites on the
        # following image. Press any key when finished.',
        #       0)
        # cv.namedWindow('TestImage')
        # cv.setMouseCallback('TestImage', set_sat_position)
        # logdata = cv.normalize(imgdata, None, alpha=0, beta=5,
        # norm_type=cv.NORM_MINMAX, dtype=cv.CV_32F)
        # cv.imshow('TestImage', logdata)
        # cv.waitKey(0)
        # cv.destroyAllWindows()
        try:
            fig = Figure()
            ax = fig.add_subplot()
            root = tk.Tk()
            ax.imshow(imgdata, cmap='gray', norm=norm, interpolation='nearest')
            canvas = FigureCanvasTkAgg(fig, master=root)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
            canvas.mpl_connect('button_press_event',
                               lambda event: onclick(event, sat_locs))
            root.update()
            root.mainloop()
        except ValueError:
            fig = Figure()
            ax = fig.add_subplot()
            root = tk.Tk()
            ax.imshow(abs(imgdata), cmap='gray', norm=LogNorm(
                vmin=1), interpolation='nearest')
            canvas = FigureCanvasTkAgg(fig, master=root)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, root)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
            canvas.mpl_connect('button_press_event',
                               lambda event: onclick(event, sat_locs))
            root.update()
            root.mainloop()
        sat_locs = np.array(sat_locs)
        print(sat_locs)
        num_sats = len(sat_locs)
        num_nans = np.zeros(num_sats, dtype=int)
        names = np.empty(num_sats + 2, dtype=object)
        names[0] = 'Time (JD)'
        names[1] = 'Filter'
        date_col = np.empty((filecount, 1))
        date_col.fill(np.nan)
        filter_col = np.empty((filecount, 1), dtype=object)
        data = np.empty((filecount, num_sats))
        data.fill(np.nan)
        auxiliary_data = np.empty((filecount, 3))
        auxiliary_data.fill(np.nan)

        for i, name in enumerate(names[2:]):
            try:
                fig = Figure()
                ax = fig.add_subplot()
                ax.imshow(imgdata, cmap='gray', norm=norm,
                          interpolation='nearest')
                sat_aperture = RectangularAperture(sat_locs[i],
                                                   w=max_distance_from_sat * 2,
                                                   h=max_distance_from_sat * 2)
                sat_aperture.plot(axes=ax, color='r', lw=1.5, alpha=0.5)
                root = tk.Tk()
                root.title("Set Satellite Position")
                img_frame = tk.Frame(root)
                img_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
                input_frame = tk.Frame(root)
                input_frame.pack(side=tk.RIGHT)
                label = tk.Label(input_frame, text='Enter Satellite')
                label.pack()
                entry = tk.Entry(input_frame)
                entry.bind("<Return>", return_entry)
                entry.pack(padx=5)
                button = tk.Button(input_frame, text="OK",
                                   command=return_entry)
                button.pack()
                canvas = FigureCanvasTkAgg(fig, master=img_frame)
                canvas.draw()
                toolbar = NavigationToolbar2Tk(canvas, img_frame)
                toolbar.update()
                canvas.get_tk_widget().pack(side=tk.LEFT,
                                            fill=tk.BOTH,
                                            expand=1)
                root.update()
                root.focus_force()
                entry.focus_set()
                root.mainloop()
            except ValueError:
                fig = Figure()
                ax = fig.add_subplot()
                ax.imshow(abs(imgdata), cmap='gray', norm=LogNorm(
                    vmin=1), interpolation='nearest')
                sat_aperture =\
                    RectangularAperture(sat_locs[i],
                                        w=max_distance_from_sat * 2,
                                        h=max_distance_from_sat * 2)
                sat_aperture.plot(axes=ax, color='r', lw=1.5, alpha=0.5)
                root = tk.Tk()
                root.title("Set Satellite Position")
                img_frame = tk.Frame(root)
                img_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
                input_frame = tk.Frame(root)
                input_frame.pack(side=tk.RIGHT)
                label = tk.Label(input_frame, text='Enter Satellite')
                label.pack()
                entry = tk.Entry(input_frame)
                entry.bind("<Return>", return_entry)
                entry.pack(padx=5)
                button = tk.Button(input_frame, text="OK",
                                   command=return_entry)
                button.pack()
                canvas = FigureCanvasTkAgg(fig, master=img_frame)
                canvas.draw()
                toolbar = NavigationToolbar2Tk(canvas, img_frame)
                toolbar.update()
                canvas.get_tk_widget().pack(side=tk.LEFT,
                                            fill=tk.BOTH,
                                            expand=1)
                root.update()
                root.focus_force()
                entry.focus_set()
                root.mainloop()
            names[i + 2] = content
            print(
                f"Satellite {names[i + 2]} at location ({sat_locs[i, 0]}, {sat_locs[i, 1]})")
        print(names)

        cmap = plt.get_cmap(cmap_set)
        colours = [cmap(i) for i in range(0, num_sats)]
        legend_elements = []
        window = tk.Tk()
        window.title('Plotting in Tkinter Test')
        img_frame = tk.Frame(window)
        img_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        input_frame = tk.Frame(window)
        input_frame.pack(side=tk.RIGHT)
        label = tk.Label(
            input_frame, text='Are the satellite positions correct?')
        label.pack()
        yes_no = tk.IntVar()
        yes_btn = tk.Radiobutton(
            input_frame, text='Yes', variable=yes_no, value=1)
        yes_btn.pack(anchor=tk.W, padx=5)
        no_btn = tk.Radiobutton(input_frame, text='No',
                                variable=yes_no, value=2)
        no_btn.pack(anchor=tk.W, padx=5)
        closebutton = tk.Button(input_frame, text='OK', command=window.destroy)
        closebutton.pack()
        try:
            fig = Figure()
            ax = fig.add_subplot()
            ax.imshow(imgdata, cmap='gray', norm=LogNorm(),
                      interpolation='nearest')
            for i in range(0, num_sats):
                sat_aperture = RectangularAperture(sat_locs[i],
                                                   w=max_distance_from_sat * 2,
                                                   h=max_distance_from_sat * 2)
                sat_aperture.plot(axes=ax, color=colours[i], lw=1.5, alpha=0.5)
                legend_elements.append(Line2D([0],
                                              [0],
                                              color='w',
                                              marker='s',
                                              markerfacecolor=colours[i],
                                              markersize=7,
                                              label=names[i + 2]))
            fig.legend(handles=legend_elements, framealpha=1)
            canvas = FigureCanvasTkAgg(fig, master=img_frame)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, img_frame)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            window.mainloop()
        except ValueError:
            fig = Figure()
            ax = fig.add_subplot()
            ax.imshow(abs(imgdata), cmap='gray', norm=LogNorm(
                vmin=1), interpolation='nearest')
            for i in range(0, num_sats):
                sat_aperture = RectangularAperture(sat_locs[i],
                                                   w=max_distance_from_sat * 2,
                                                   h=max_distance_from_sat * 2)
                sat_aperture.plot(axes=ax, color=colours[i], lw=1.5, alpha=0.5)
                legend_elements.append(Line2D([0],
                                              [0],
                                              color='w',
                                              marker='s',
                                              markerfacecolor=colours[i],
                                              markersize=7,
                                              label=names[i + 2]))
            fig.legend(handles=legend_elements, framealpha=1)
            canvas = FigureCanvasTkAgg(fig, master=img_frame)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, img_frame)
            toolbar.update()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            window.mainloop()
        # Have a way for the user to confirm the satellite locations.
        # If it is wrong, then decide whether to change
        # set_sat_positions to True/False or change_sat_positions
        if yes_no.get() == 1:
            set_sat_positions_bool = False
        else:
            continue
        # elif yes_no.get() == 2:
        #     sat_locs = []
        #     names = []
        # print(names[0])
        # print(date_col)
        sat_names = names[2:]
        date_table = Table(names=[names[0]], data=date_col)
        filter_table = Table(names=[names[1]], data=filter_col)
        data_table = Table(names=names[2:], data=data)
        auxiliary_column_names = ["FWHM", "Airmass", "BSB"]
        ausiliary_columns = Table(
            names=auxiliary_column_names, data=auxiliary_data)
        sats_table = hstack(
            [date_table, filter_table, data_table], join_type='exact')
        uncertainty_table = hstack(
            [date_table, filter_table, data_table], join_type='exact')
        sat_auxiliary_table = hstack(
            [date_table, filter_table, ausiliary_columns], join_type='exact')
        sats_table.pprint_all()
    sat_information = namedtuple('sat_information',
                                 ['sats_table',
                                  'uncertainty_table',
                                  'sat_auxiliary_table',
                                  'sat_locs',
                                  'num_sats',
                                  'num_nans',
                                  'sat_names'])
    return set_sat_positions_bool, sat_information(sats_table,
                                                   uncertainty_table,
                                                   sat_auxiliary_table,
                                                   sat_locs,
                                                   num_sats,
                                                   num_nans,
                                                   sat_names)
    # return sats_table, uncertainty_table, sat_fwhm_table, sat_locs,
    # num_sats, num_nans, sat_names


def plot_detected_sats(filename,
                       plot_results,
                       imgdata,
                       irafsources,
                       sat_information,
                       max_distance_from_sat=20,
                       norm=LogNorm()):
    """
    TODO.

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    plot_results : TYPE
        DESCRIPTION.
    imgdata : TYPE
        DESCRIPTION.
    irafsources : TYPE
        DESCRIPTION.
    sat_information : TYPE
        DESCRIPTION.
    max_distance_from_sat : TYPE, optional
        DESCRIPTION. The default is 20.
    norm : TYPE, optional
        DESCRIPTION. The default is LogNorm().

    Returns
    -------
    None.

    """
    if plot_results != 0:
        fig, ax = plt.subplots()
        ax.imshow(imgdata, cmap='gray', norm=norm, interpolation='nearest')
        ax.scatter(irafsources['xcentroid'], irafsources['ycentroid'],
                   s=100, edgecolor='red', facecolor='none')
        for i in range(0, sat_information.num_sats):
            rect = patches.Rectangle((sat_information.sat_locs[i, 0] -
                                      max_distance_from_sat,
                                      sat_information.sat_locs[i, 1] -
                                      max_distance_from_sat),
                                     width=max_distance_from_sat * 2,
                                     height=max_distance_from_sat * 2,
                                     edgecolor='green', facecolor='none')
            ax.add_patch(rect)
        plt.title(filename)
        if plot_results == 1:
            plt.show(block=False)
            plt.pause(2)
        elif plot_results == 2:
            plt.show()
        plt.close()


def _main_sc_lightcurve(directory,
                        gb_final_transforms=None,
                        temp_dir='tmp',
                        save_loc='Outputs',
                        file_suffix=(".fits", ".fit", ".fts"),
                        ecct_cut=0.5,
                        max_distance_from_sat=20,
                        size=25,
                        max_num_nan=5,
                        plot_results=0):
    filecount, filenames = astro.copy_and_rename(directory=directory,
                                           file_suffix=file_suffix,
                                           temp_dir=temp_dir,
                                           debugging=True)
    # remove_temp_dir(temp_dir=temp_dir)
    # return
    set_sat_positions_bool = True
    change_sat_positions_bool = False
    num_nan = 0
    for filenum, file in enumerate(filenames):
        filepath = f"{temp_dir}/{file}"
        hdr, imgdata = astro.read_fits_file(filepath)
        if set_sat_positions_bool:
            set_sat_positions_bool, sat_information = set_sat_positions(
                imgdata, filecount, set_sat_positions_bool)
        sat_information = astro.add_new_time_and_filter(
            hdr, sat_information, filenum)
        if change_sat_positions_bool:
            change_sat_positions_bool,\
                sat_information\
                = astro.change_sat_positions(filenames,
                                       filenum,
                                       num_nan,
                                       sat_information,
                                       change_sat_positions_bool,
                                       gb_final_transforms=gb_final_transforms,
                                       max_distance_from_sat=max_distance_from_sat,
                                       size=size,
                                       temp_dir=temp_dir,
                                       cmap_set='Set1',
                                       plot_results=plot_results)
        exptime = hdr['EXPTIME'] * u.s
        # bkg, bkg_std = calculate_img_bkg(imgdata)

        # try:
        # sat_x, sat_y, bkg_trm, fwhm = TRM_sat_detection(
        #     filepath, ecct_cut=ecct_cut)
        tbl, bkg_trm, fwhms = sd.streak_detection_single(filepath, sigma=3.0)
        fwhm = np.mean(fwhms)
        fwhm_std = np.mean(fwhms)
        sat_x, sat_y = sd.filter_sats_stars(tbl, ecct_cut=ecct_cut)
        # except TypeError:
        if len(sat_x) == 0 or len(sat_y) == 0:
            print("No satellites detected.")
            continue
        if fwhm < 0:
            continue
        # bkg, bkg_std = calculate_img_bkg(imgdata)
        # irafsources = detecting_stars(imgdata, bkg, bkg_std)
        # if not irafsources:
        #     sat_information.num_nans[:] = 0
        #     continue
        # plot_detected_sats(filenames[filenum],
        #                    plot_results,
        #                    imgdata,
        #                    irafsources,
        #                    sat_information,
        #                    max_distance_from_sat=max_distance_from_sat,
        #                    norm=LogNorm())
        # fwhms, fwhm, fwhm_std = calculate_fwhm(irafsources)
        bkg = np.median(bkg_trm)
        bsb = astro.calculate_background_sky_brightness(bkg,
                                                  hdr,
                                                  exptime,
                                                  gb_final_transforms,
                                                  focal_length_key='FOCALLEN',
                                                  xpixsz_key='XPIXSZ',
                                                  ypixsz_key='YPIXSZ')
        # TODO: Incorporate the standard deviation.
        _, fwhm_arcsec, _ = astro.convert_fwhm_to_arcsec(hdr, fwhms, fwhm, fwhm_std)
        airmass = astro.get_image_airmass(hdr)
        photometry_result = perform_photometry.perform_PSF_photometry_sat(
            sat_x, sat_y, fwhm, imgdata, bkg_trm)

        fluxes=photometry_result['flux_fit']
        fluxes_unc=photometry_result['flux_unc']
        instr_mags = astro.calculate_magnitudes(fluxes, exptime)
        instr_mags_sigma, snr = astro.calculate_magnitudes_sigma(
            fluxes,fluxes_unc, exptime)
        sat_information = astro.check_if_sat(sat_information,
                                       filenum,
                                       sat_x,
                                       sat_y,
                                       instr_mags,
                                       instr_mags_sigma,
                                       fwhm_arcsec,
                                       airmass,
                                       bsb,
                                       max_distance_from_sat=max_distance_from_sat)
        change_sat_positions_bool, num_nan = \
            astro.determine_if_change_sat_positions(sat_information,
                                                                               filenum,
                                                                               change_sat_positions_bool,
                                                                               max_num_nan=max_num_nan)
        # del hdr
        # del imgdata
    # try:
    astro.remove_temp_dir(temp_dir=temp_dir)
    # except PermissionError as e:
    #     print(e)
    if not os.path.exists(save_loc):
        os.mkdir(save_loc)
    sats_table = sat_information.sats_table
    ascii.write(
        sats_table, output=f"{save_loc}/Measured_Magnitudes.csv", format='csv')
    uncertainty_table = sat_information.uncertainty_table
    ascii.write(uncertainty_table,
                output=f"{save_loc}/Measured_Magnitude_Uncertainties.csv", format='csv')
    sat_auxiliary_table = sat_information.sat_auxiliary_table
    for aux_data in sat_auxiliary_table.columns[2:]:
        if all(np.isnan(sat_auxiliary_table[aux_data])):
            sat_auxiliary_table.remove_column(aux_data)
    ascii.write(sat_auxiliary_table,
                output=f"{save_loc}/Auxiliary_Information.csv", format='csv')

    # Place rows without nan or masked values into a new table which can be processed easier later
    # TODO: Make this code more elegant
    sats_table2 = Table(sats_table[0])
    uncertainty_table2 = Table(uncertainty_table[0])
    sat_auxiliary_table2 = Table(sat_auxiliary_table[0])
    for row in range(1, len(sats_table)):
        switch = True
        for element in range(0, len(sats_table.columns)):
            if ((str(sats_table[row][element]) == 'nan') or
                    (bool(sats_table[row][element]) is False)):
                switch = False
        for element in range(0, len(uncertainty_table.columns)):
            if ((str(uncertainty_table[row][element]) == 'nan') or
                    (bool(uncertainty_table[row][element]) is False)):
                switch = False
        for element in range(0, len(sat_auxiliary_table.columns)):
            if ((str(sat_auxiliary_table[row][element]) == 'nan') or
                    (bool(sat_auxiliary_table[row][element]) is False)):
                switch = False
        if switch is True:
            sats_table2.add_row(sats_table[row])
            uncertainty_table2.add_row(uncertainty_table[row])
            sat_auxiliary_table2.add_row(sat_auxiliary_table[row])
        switch = True
    sats_table = sats_table2
    uncertainty_table = uncertainty_table2
    sat_auxiliary_table = sat_auxiliary_table2

    unique_filters, num_filters, multiple_filters = \
        astro.determine_num_filters(
        sats_table)
    if multiple_filters:
        sat_dict = astro.interpolate_sats(
            sats_table, uncertainty_table, unique_filters)
        if not gb_final_transforms:
            app_sat_dict = None
            astro.save_interpolated_light_curve(sat_dict, save_loc)
            all_indices, all_indices_formatted = \
                astro.get_all_indicies_combinations(unique_filters,
                                                                               num_filters,
                                                                               multiple_filters)
            colour_indices_dict = astro.calculate_timeseries_colour_indices(
                sat_dict, all_indices)
            astro.save_interpolated_light_curve(
                colour_indices_dict, save_loc, suffix="Colour Indices")
            filters_to_plot, indices_to_plot, aux_data_to_plot\
                = astro.choose_indices_to_plot(unique_filters,
                                         num_filters,
                                         all_indices_formatted,
                                         sat_auxiliary_table)
            fig = astro.axis_limits_multiband_gui(sat_dict,
                                            colour_indices_dict,
                                            sat_auxiliary_table,
                                            filters_to_plot,
                                            indices_to_plot,
                                            aux_data_to_plot,
                                            save_loc)
        else:
            app_sat_dict = astro.apply_gb_timeseries_transforms(
                gb_final_transforms,
                                                          sat_dict,
                                                          sat_auxiliary_table,
                                                          unique_filters)
            astro.save_interpolated_light_curve(app_sat_dict, save_loc)
            all_indices, all_indices_formatted = \
                astro.get_all_indicies_combinations(unique_filters,
                                                                               num_filters,
                                                                               multiple_filters)
            colour_indices_dict = astro.calculate_timeseries_colour_indices(
                app_sat_dict, all_indices)
            astro.save_interpolated_light_curve(
                colour_indices_dict, save_loc, suffix="Colour Indices")
            filters_to_plot, indices_to_plot, aux_data_to_plot = \
                astro.choose_indices_to_plot(unique_filters,
                                                                                        num_filters,
                                                                                        all_indices_formatted,
                                                                                        sat_auxiliary_table)
            fig = astro.axis_limits_multiband_gui(app_sat_dict,
                                            colour_indices_dict,
                                            sat_auxiliary_table,
                                            filters_to_plot,
                                            indices_to_plot,
                                            aux_data_to_plot,
                                            save_loc)
    else:
        sat_dict = None
        app_sat_dict = None
        aux_data_to_plot = astro.choose_aux_data_to_plot(sat_auxiliary_table)
        astro.axis_limits_singleband_gui(sats_table,
                                   uncertainty_table,
                                   sat_auxiliary_table,
                                   aux_data_to_plot,
                                   save_loc)
        # plot_light_curve_singleband(sats_table, uncertainty_table,
        # sat_auxiliary_table, aux_data_to_plot, save_loc)
    return sat_dict, app_sat_dict, sats_table, uncertainty_table, sat_auxiliary_table