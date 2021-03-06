
import os
import sys
import multiprocessing as mp

import tqdm
import numpy as np
from scipy import optimize
from scipy import fftpack

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval
from photutils import DAOStarFinder, CircularAperture, aperture_photometry

import sys
from os.path import dirname
src_path = dirname(dirname(__file__))
sys.path.append(os.path.join(src_path, 'general_tools'))
sys.path.append(os.path.join(src_path, 'data_visualization'))

import utils
import Visualize


def columncor(scidata, bpix):
    """"""

    scidata_masked = np.ma.array(scidata, mask=scidata < bpix)
    med = np.ma.median(scidata_masked, axis=0)
    scidata_colcor = scidata - med

    return scidata_colcor


def combine(imagefiles, ilow, ihigh, bpix):
    """Usage: masterimage = combine(imagefiles)"""

    # Read in first image.
    image1 = utils.read_fitsdata(imagefiles[0])

    # Get the number of images and the image dimensions.
    nfiles = len(imagefiles)
    n1, n2 = image1.shape

    # Allocate the necessary arrays.
    masterimage = np.zeros((n1, n2))
    allfitsdata = np.zeros((nfiles, n1, n2))

    # Store first image in array.
    allfitsdata[0, :, :] = image1

    # Loop over the remaining images.
    for icount, filename in enumerate(imagefiles):

        # Skip first image (already in array).
        if icount > 0:
            allfitsdata[icount, :, :] = utils.read_fitsdata(filename)

    for i in range(n1):
        for j in range(n2):
            pixels = []
            for k in range(nfiles):
                if allfitsdata[k, i, j] > bpix:  # Exclude bad-pixels.
                    pixels.append(allfitsdata[k, i, j])
            pixels = np.array(pixels)
            npixels = len(pixels)
            if npixels < 1:
                masterimage[i, j] = bpix
            elif npixels == 1:
                masterimage[i, j] = pixels[0]
            else:
                pixels = np.sort(pixels)
                i1 = 0 + ilow
                i2 = npixels - ihigh
                if i1 > i2:
                    i1 = npixels/2
                    i2 = npixels/2
                masterimage[i, j] = np.sum(pixels[i1:i2])/float(i2-i1)

                # print(pixels)
                # print(i1, i2)
                # print(pixels[i1:i2])
                # print(masterimage[i, j])
                # input()

    return masterimage


def fourierd2d(a, xn, yn, xoff, yoff):
    """"""

    n = len(a)  # Number of parameters in model.
    tpi = 2.0*np.pi
    m = a[0]*np.ones([int(xn), int(yn)])  # Zero point.
    xx, yy = np.indices((xn, yn))
    for k in range(1, n, 4):
        m += a[k]*np.sin(tpi*(a[k+1]*(xx - xoff) + a[k+2]*(yy - yoff)) + a[k+3])

    return m


def func(a, xn, yn, xoff, yoff, overscan):  # TODO more descriptive function name.
    """"""

    model = fourierd2d(a, xn, yn, xoff, yoff)
    sqmeanabs = np.sqrt(np.mean(np.abs(overscan)))

    diff = (overscan - model)/sqmeanabs
    diffflat = diff.flatten()

    return diffflat


def funcphase(aoff, a, xn, yn, scidata_in, stdcut=None):
    """Determine phase offset for science image."""

    xoff = aoff[0]
    yoff = aoff[1]
    model = fourierd2d(a, xn, yn, xoff, yoff)
    sqmeanabs = np.sqrt(np.mean(np.abs(scidata_in)))

    if sqmeanabs > 0:  # TODO how could this not be >0? Maybe not finite?
        diff = (scidata_in - model)/sqmeanabs
    else:
        diff = (scidata_in - model)

    if stdcut is not None:

        diffflat = diff.flatten()
        diffflat[np.abs(diffflat) > stdcut] = 0.0

        return diffflat

    return diff


def fouriercor(scidata_in, a):
    """Apply Fourier correction from overscan."""

    xn, yn = scidata_in.shape
    scidata_z = scidata_in - np.median(scidata_in)

    # Perform an initial fit.
    aoff = np.array([0.5, 0.5])
    stdcut = 1.0e30
    aoff, ier = optimize.leastsq(funcphase, aoff, args=(a, xn, yn, scidata_z, stdcut), factor=1)

    # Apply a sigma cut, to reduce the effect of stars in the image.
    diff = funcphase(aoff, a, xn, yn, scidata_z)
    stdcut = 3.0*np.std(diff)
    aoff, ier = optimize.leastsq(funcphase, aoff, args=(a, xn, yn, scidata_z, stdcut), factor=1)

    # Remove the final model from the data.
    xoff, yoff = aoff
    model = fourierd2d(a, xn, yn, xoff, yoff)
    scidata_cor = scidata_in - model

    return scidata_cor


def overscan_cor(scidata_c, overscan, a, bpix):
    """"""

    scidata_co = fouriercor(scidata_c, a)

    # General Overscan correction.
    xn = overscan.shape[0]
    yn = overscan.shape[1]
    model = fourierd2d(a, xn, yn, 0.0, 0.0)
    overscan_cor1 = overscan - model

#    row_cor = np.mean(overscan_cor1, axis=1, keepdims=True)
#    scidata_cor = scidata_co - row_cor

    row_cor, _, _ = sigma_clipped_stats(overscan_cor1, axis=1)
    scidata_cor = scidata_co - row_cor[:, np.newaxis]

    return scidata_cor


def darkprocess(workdir, darkfile, xsc, ysc, xov, yov, snrcut, fmax, xoff, yoff, T, bpix):
    """"""

    info = 0

    filename = os.path.join(workdir, darkfile)
    scidata = utils.read_fitsdata(filename)

    # Crop Science Image.
    sh = scidata.shape
    strim = np.array([sh[0] - xsc, sh[0], sh[1] - ysc, sh[1]])
    scidata_c = np.copy(scidata[strim[0]:strim[1], strim[2]:strim[3]])

    # Crop Overscan.
    sh = scidata.shape
    otrim = np.array([sh[0] - xov, sh[0], 0, yov])
    overscan = np.copy(scidata[otrim[0]:otrim[1], otrim[2]:otrim[3]])
    med = np.median(overscan, axis=0)
    overscan = overscan - med
    mean = np.mean(med)
    overscan = overscan + mean  # Add mean back to overscan (this is the BIAS).

    # Fourier Decomp of overscan.
    a = fourierdecomp(overscan, snrcut, fmax, xoff, yoff, T, bpix, info=info)

    # Apply overscan correction to science raster.
    scidata_cor = overscan_cor(scidata_c, overscan, a, bpix)

    return scidata_cor


def scale_image(image, ref_image, mind=0, maxd=8000, b1=100, m1=0.3, m2=1.3, tp=2000):
    """"""

    data = image.flatten()
    ref_data = ref_image.flatten()

    mask = (data > mind) & (data < maxd) & (ref_data > mind) & (ref_data < maxd)

    if np.sum(mask) > 10:
        data_bin, ref_data_bin, err_bin = utils.bindata(data[mask], ref_data[mask], 50)

        x0 = [b1, m1, m2, tp]
        ans = optimize.least_squares(ls_seg_func, x0, args=[data_bin, ref_data_bin, err_bin])

        scaled_data = seg_func(ans.x, data)
        scaled_image = scaled_data.reshape(image.shape)

    return scaled_image


def scale_image_zscale(image, ref_image, b1=100, m1=0.3, m2=1.3, tp=2000):
    """"""

    data = image.flatten()
    ref_data = ref_image.flatten()

    min1, max1 = ZScaleInterval().get_limits(data)
    min2, max2 = ZScaleInterval().get_limits(ref_data)

    mask = (data > min1) & (data < max1) & (ref_data > min2) & (ref_data < max2)

    if np.sum(mask) > 10:
        data_bin, ref_data_bin, err_bin = utils.bindata(data[mask], ref_data[mask], 50)

        x0 = [b1, m1, m2, tp]
        ans = optimize.least_squares(ls_seg_func, x0, args=[data_bin, ref_data_bin, err_bin])

        scaled_data = seg_func(ans.x, data)
        scaled_image = scaled_data.reshape(image.shape)

    return scaled_image


def combinedarks(alldarkdata, mind=0, maxd=8000, b1=100, m1=0.3, m2=1.3, tp=2000):
    """
    mind,maxd : range of data to consider when matching frames.  Keeping maxd relatively low avoids stars
    [b1,m1,m2,tp] - initial guess for solution.
    b1=y-intercept for first segment
    m1=slope for first segment
    m2=slope for second segment
    tp=division point from first to second segment
    """

    darkscaled = []
    ndark = len(alldarkdata)
    for i in range(1, ndark):

        image = alldarkdata[i]
        ref_image = alldarkdata[0]

        # newdark = scale_image(image, ref_image, mind=mind, maxd=maxd, b1=b1, m1=m1, m2=m2, tp=tp)
        newdark = scale_image_zscale(image, ref_image, b1=b1, m1=m1, m2=m2, tp=tp)
        darkscaled.append(newdark)

    darkscaled = np.array(darkscaled)
    darkavg = np.median(darkscaled, axis=0)

    return darkavg


def find_line_model(points):
    """ find a line model for the given points
    :param points selected points for model fitting
    :return line model
    """

    # [WARNING] vertical and horizontal lines should be treated differently
    #           here we just add some noise to avoid division by zero

    # find a line model for these points
    m = (points[1, 1] - points[0, 1])/(points[1, 0] - points[0, 0] + sys.float_info.epsilon)  # Slope (gradient) of the line.
    c = points[1, 1] - m*points[1, 0]  # y-intercept of the line.

    return m, c


def find_intercept_point(m, c, x0, y0):
    """ find an intercept point of the line model with
        a normal from point (x0,y0) to it
    :param m slope of the line model
    :param c y-intercept of the line model
    :param x0 point's x coordinate
    :param y0 point's y coordinate
    :return intercept point
    """

    # Intersection point with the model.
    x = (x0 + m*y0 - m*c)/(1 + m**2)
    y = (m*x0 + (m**2)*y0 - (m**2)*c)/(1 + m**2) + c

    return x, y


def darkcorrect(scidata, masterdark, bpix):
    """Usage: m, c = darkcorrect(scidata, masterdark, pbix)"""

    masked_dark = np.ma.array(masterdark, mask=masterdark < bpix)
    maxd = masked_dark.max()
    mind = masked_dark.min()

    mask = (scidata > bpix) & (masterdark > bpix) & (scidata > mind) & (scidata < maxd)

    x = masterdark[mask]
    y = scidata[mask]

    n_samples = len(x)

    # Ransac parameters.
    ransac_iterations = 20  # Number of iterations.
    ransac_threshold = 3  # Threshold.
    ransac_ratio = 0.6  # Ratio of inliers required to assert that a model fits well to data.

    # data = np.hstack((x, y))
    data = np.vstack((x, y)).T
    ratio = 0.
    model_m = 0.
    model_c = 0.

    for it in range(ransac_iterations):
        # Pick up two random points.
        n = 2

        all_indices = np.arange(x.shape[0])
        np.random.shuffle(all_indices)

        indices_1 = all_indices[:n]
        indices_2 = all_indices[n:]

        maybe_points = data[indices_1, :]
        test_points = data[indices_2, :]

        # Find a line model for these points.
        m, c = find_line_model(maybe_points)

        x_list = []
        y_list = []
        num = 0

        # Find orthogonal lines to the model for all testing points.
        for ind in range(test_points.shape[0]):

            x0 = test_points[ind, 0]
            y0 = test_points[ind, 1]

            # Find an intercept point of the model with a normal from point (x0, y0).
            x1, y1 = find_intercept_point(m, c, x0, y0)

            # Distance from point to the model.
            dist = np.sqrt((x1 - x0)**2 + (y1 - y0)**2)

            # Check whether it's an inlier or not.
            if dist < ransac_threshold:
                x_list.append(x0)
                y_list.append(y0)
                num += 1

        x_inliers = np.array(x_list)
        y_inliers = np.array(y_list)

        # In case a new model is better - cache it.
        if num/float(n_samples) > ratio:
            ratio = num/float(n_samples)
            model_m = m
            model_c = c

        # print ('  inlier ratio = ', num/float(n_samples))
        # print ('  model_m = ', model_m)
        # print ('  model_c = ', model_c)

        # Plot the current step.
        # ransac_plot(it, x_noise,y_noise, m, c, False, x_inliers, y_inliers, maybe_points)

        # We are done in case we have enough inliers.
        if num > n_samples*ransac_ratio:
            # print ('The model is found !')
            break

    # print ('\nFinal model:\n')
    # print ('  ratio = ', ratio)
    # print ('  model_m = ', model_m)
    # print ('  model_c = ', model_c)

    return model_m, model_c


def seg_func(x0, data):
    """"""

    b1 = x0[0]
    m1 = x0[1]
    m2 = x0[2]
    tp = x0[3]
    # print(x0)

    b2 = m1*tp + b1 - m2*tp

    ans = np.where(data < tp, m1 * data + b1, m2 * data + b2)
    ans = ans.flatten()

    return ans


def ls_seg_func(x0, data1, data2, derr):
    """"""

    ans = seg_func(x0, data1)

    diff = (data2 - ans)/(derr + 1.0e-20)

    return diff


def fourierdecomp(overscan, snrcut, fmax, xoff, yoff, T, bpix, info=0):
    """"""

    # Count number of frequencies.
    nfreqs = 0

    # Calculate Median of overscan region.
    med_overscan = np.median(overscan)
    std_overscan = np.std(overscan - med_overscan)

    # Size of Overscan.
    xn = overscan.shape[0]
    yn = overscan.shape[1]

    # Oversampled overscan.
    overscan_os = np.zeros([xn*T, yn*T])
    overscan_os[:xn, :yn] = np.copy(overscan)

    # Initialize model.
    a = np.zeros(1)
    model = np.zeros([xn, yn])

    # Frequency Grid.
    # xf = np.linspace(0.0, 1.0/(2.0), T*xn//2)
    xf = np.append(np.linspace(0.0, 1.0/2.0, T*xn//2), -np.linspace(1.0/2.0, 0.0, T*xn//2))
    yf = np.linspace(0.0, 1.0/2.0, T*yn//2)

    if fmax > 0:
        loop = True
    else:
        loop = False

    while loop:

        # Remove median, model and then calculate FFT.
        overscan_os[:xn, :yn] = overscan - med_overscan - model
        ftoverscan = fftpack.fft2(overscan_os)
        ftoverscan_abs = np.abs(ftoverscan)/(xn*yn)  # Amplitude.

        if info >= 2:
            # Plot the FFT.
            imstat = utils.imagestat(ftoverscan_abs, bpix)
            Visualize.plot_image(np.transpose(np.abs(ftoverscan_abs[:T*xn, :T*yn//2])), imstat, 0.0, 10.0)

        mean_ftoverscan_abs = np.mean(ftoverscan_abs[T:T*(xn - 1), T:T*yn//2])
        std_ftoverscan_abs = np.std(ftoverscan_abs[T:T*(xn - 1), T:T*yn//2])
        if info >= 1:
            print('mean, std:', mean_ftoverscan_abs, std_ftoverscan_abs)

        # Locate Frequency with largest amplitude.
        max_array = ftoverscan_abs[T:T*(xn - 1), T:T*yn//2]
        maxi, maxj = np.unravel_index(np.argmax(max_array), max_array.shape)
        maxi += T
        maxj += T

        snr = (ftoverscan_abs[maxi, maxj] - mean_ftoverscan_abs)/std_ftoverscan_abs
        if info >= 1:
            print('SNR,i,j,amp: ', snr, maxi, maxj, ftoverscan_abs[maxi, maxj])
            # print(ftoverscan_abs[maxi, maxj], xf[maxi], yf[maxj], np.angle(ftoverscan[maxi, maxj]))

        if (snr > snrcut) and (nfreqs < fmax):
            nfreqs += 1
            a = np.append(a, [ftoverscan_abs[maxi, maxj], xf[maxi], yf[maxj], np.angle(ftoverscan[maxi, maxj])+np.pi/2])
            # a1 = np.array([ftoverscan_abs[maxi, maxj], xf[maxi], yf[maxj], np.angle(ftoverscan[maxi, maxj])-np.pi/2])
            # a1 = np.append(0.0, a1)
            if info >= 1:
                print('Next Mode: (amp, xf, yf, phase)')
                print(ftoverscan_abs[maxi, maxj], xf[maxi], yf[maxj], np.angle(ftoverscan[maxi, maxj]) + np.pi/2)
            ans = optimize.leastsq(func, a, args=(xn, yn, xoff, yoff, overscan - med_overscan), factor=1)
            # a = np.append(a, ans[0][1:])
            a = ans[0]
            model = fourierd2d(a, xn, yn, xoff, yoff)
            n = len(a)
            if info >= 1:
                print("---Solution---")
                print("zpt: ", a[0])
                for k in range(1, n, 4):
                    print(a[k], a[k+1], a[k+2], a[k+3])
                print("--------------")

        else:
            loop = False
            # Remove median, model and then calculate FFT.
            # overscan_os[:xn, :yn] = overscan - med_overscan - model
            # ftoverscan = fft2(overscan_os)
            # ftoverscan_abs = np.abs(ftoverscan)/(xn*yn)  # Amplitude.

            # Plot the FFT.
            # if info >= 2:
            #     imstat = imagestat(ftoverscan_abs, bpix)
            #     plot_image(np.abs(np.transpose(ftoverscan_abs[:T, :T*yn//2])), imstat, 0.0, 10.0)

            if info >= 1:
                print('Done')

    return a


def clean_sciimage(filename, darkavg, xsc, ysc, xov, yov, snrcut, fmax, xoff, yoff, T, info, bpix):
    """"""

    cor = False  # Updates from Hamza.
    dark = False
    scidata = utils.read_fitsdata(filename)
    scidata_cor = None
    scidata_cord = None

    # Updates from Hamza.
    # This part examines if overscan+cropping is needed and check if dark is valid.
    hdul = fits.open(filename, mode="update")
    hdr = hdul[0].header
    NAXIS1 = hdr['NAXIS1']
    NAXIS2 = hdr['NAXIS2']
    hdul.close()
    # Set flags of what needs to be performed.
    if not NAXIS2 == xsc or not NAXIS1 == ysc:
        cor = True
    if len(darkavg) != 0:
        dark = True

    # If we only need to perform dark correction then set the scidata_cor to scidata.
    if not cor and dark:
        scidata_cor = scidata

    if cor:

        # Crop Science Image.
        sh = scidata.shape
        strim = np.array([sh[0]-xsc, sh[0], sh[1]-ysc, sh[1]])
        scidata_c = np.copy(scidata[strim[0]:strim[1], strim[2]:strim[3]])

        # Crop Overscan.
        sh = scidata.shape

        otrim = np.array([sh[0]-xov, sh[0], 0, yov])
        overscan = np.copy(scidata[otrim[0]:otrim[1], otrim[2]:otrim[3]])
        mean = 0.0
        for i in range(yov):
            med = np.median(overscan[:, i])
            overscan[:, i] = overscan[:, i]-med
            mean = mean+med
        mean = mean/yov
        overscan = overscan+mean  # Add mean back to overscan (this is the BIAS).

        if info >= 2:
            imstat = utils.imagestat(overscan, bpix)
            Visualize.plot_image(np.transpose(overscan), imstat, 0.3, 3.0)

        # Fourier Decomp of overscan.
        a = fourierdecomp(overscan, snrcut, fmax, xoff, yoff, T, bpix, info=info)

        if info >= 2:
            xn = overscan.shape[0]
            yn = overscan.shape[1]
            model = fourierd2d(a, xn, yn, xoff, yoff)
            imstat = utils.imagestat(overscan-model, bpix)
            Visualize.plot_image(np.transpose(overscan-model), imstat, 0.3, 3.0)

        # Apply overscan correction to science raster
        scidata_cor = overscan_cor(scidata_c, overscan, a, bpix)

    if dark:
        # Apply Dark correction

        # OLD Dark correction REQUIRES meddif FORTRAN external.
        # image1 = darkavg
        # mind = darkavg.min()
        # maxd = darkavg.max()
        # image2 = scidata_cor
        # data1 = image1.flatten()
        # data2 = image2.flatten()
        # data1t = data1[(data1 > mind) & (data1 < maxd) & (data2 > mind) & (data2 < maxd)]
        # data2t = data2[(data1 > mind) & (data1 < maxd) & (data2 > mind) & (data2 < maxd)]
        # data1 = np.copy(data1t)
        # data2 = np.copy(data2t)
        # ndata = len(data1)
        # abdev = 1.0
        # if ndata > 3:
        #    a, b = medfit.medfit(data1, data2, ndata, abdev)
        # else:
        #    a = 0.0
        #    b = 1.0
        # scidata_cord = scidata_cor-(a+b*darkavg)

        # New Dark-correction. Not extensively tested. No Fortran dependence.
        # newdark = scale_image(darkavg, scidata_cor)
        newdark = scale_image_zscale(darkavg, scidata_cor)
        scidata_cord = scidata_cor - newdark

    # Return if only clipping and overscan is performed.
    if cor and not dark:
        # print ("Only performed clipping and overscan")
        return scidata_cor
    # Return if only dark correction is performed.
    elif not cor and dark:
        # print ("Only performed dark correction")
        return scidata_cord
    # Return if both clipping, overscan and dark correction is performed.
    elif cor and dark:
        # print ("Performed both clipping and dark correction")
        return scidata_cord
    # Return original scidata if nothing was performed.
    else:
        # print ("No clipping, overscan and dark correction requested")
        return scidata


def lightprocess(filename, date, darkavg, xsc, ysc, xov, yov, snrcut, fmax, xoff, yoff, T, photap, bpix):
    """"""

    info = 0

    scidata_cord = clean_sciimage(filename, darkavg, xsc, ysc, xov, yov, snrcut, fmax, xoff, yoff, T, info, bpix)

    mean, median, std = sigma_clipped_stats(scidata_cord, sigma=3.0, maxiters=5)

    daofind = DAOStarFinder(fwhm=2.0, threshold=5.*std)
    sources = daofind(scidata_cord - median)

    positions = np.column_stack([sources['xcentroid'], sources['ycentroid']])
    apertures = CircularAperture(positions, r=photap)
    phot_table = aperture_photometry(scidata_cord-median, apertures)

    # photall.append(phot_table)
    # photstat.append([mean, median, std])

    return [phot_table, date, mean, median, std, scidata_cord]


def lightprocess_save(filename, savedir, darkavg, xsc, ysc, xov, yov, snrcut, fmax, xoff, yoff, T, bpix):
    """"""

    info = 0

    scidata_cord = clean_sciimage(filename, darkavg, xsc, ysc, xov, yov, snrcut, fmax, xoff, yoff, T, info, bpix)

    # Set up new file name for cleaned image.
    base = os.path.basename(filename)
    x = os.path.splitext(base)
    newfile = os.path.join(savedir, x[0] + "_cord.fits")

    # Write the file.
    # header = fits.getheader(filename)
    # fits.writeto(newfile, scidata_cord, header, overwrite=True)
    header = fits.getheader(filename)  # Make copy of original header to insert.
    header['BZERO'] = 0  # Make sure BZERO and BSCALE are set.
    header['BSCALE'] = 1.0
    hdu = fits.PrimaryHDU(scidata_cord)
    hdu.scale('int32')  # Scaling to 32-bit integers.

    for key in header:
        if key not in ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'EXTEND']:
            hdu.header.append((key, header[key]))

    hdu.writeto(newfile, overwrite=True)

    return info


def Darksub(target, obspath, savedir, **kwargs):
    """Process all observations of a specific target in a specific directory."""

    # Make sure output directory exists.
    utils.ensure_dir(savedir)

    # Unpack parameters.
    T = kwargs.pop('T', 8)
    bpix = kwargs.pop('bpix', -1.0e10)
    snrcut = kwargs.pop('snrcut', 10.0)
    fmax = kwargs.pop('fmax', 2)
    xoff = kwargs.pop('xoff', 0)
    yoff = kwargs.pop('yoff', 0)
    nproc = kwargs.pop('nproc', 4)

    print('Processing observations in directory {}'.format(obspath))

    # Create a table of all fits files in the specified diretcory.
    obs_table = utils.observation_table(obspath)

    # Parse the observation to select the desired target and appropriate darks.
    light_table, dark_table = utils.parse_observation_table(obs_table, target)
    nlight, ndark = len(light_table), len(dark_table)

    # Process the dark images.
    print('Processing {} dark images to create the masterdark.'.format(ndark))

    # Use multiprocessing to process all dark images.
    results = []
    pbar = tqdm.tqdm(total=ndark)
    with mp.Pool(nproc) as p:

        for i in range(ndark):

            darkfile = dark_table['FILENAME'][i]
            xsc, ysc = dark_table['xsc'][i], dark_table['ysc'][i]
            xov, yov = dark_table['xov'][i], dark_table['yov'][i]

            args = (obspath, darkfile, xsc, ysc, xov, yov, snrcut, fmax, xoff, yoff, T, bpix)

            results.append(p.apply_async(darkprocess, args=args, callback=lambda x: pbar.update()))

        p.close()
        p.join()

        alldarkdata = [result.get() for result in results]

    pbar.close()

    # Combine the processed darks to obtain a master dark.
    masterdark = combinedarks(alldarkdata)

    # Display the master dark.
#    imstat = utils.imagestat(masterdark, bpix)
#    visualize.plot_image(masterdark, imstat, 0.3, 10.0)

    # Save the masterdark.
    head, tail = os.path.split(obspath)
    darkname = 'masterdark_{}_{}.fits'.format(tail, target)
    darkname = os.path.join(savedir, darkname)
    hdu = fits.PrimaryHDU(masterdark)
    hdu.writeto(darkname, overwrite=True)

    # Clear memory.
    del alldarkdata

    # Process the light images.
    print('Processing {} light images.'.format(nlight))

    # Use multiproessing to process all light images.
    pbar = tqdm.tqdm(total=nlight)
    with mp.Pool(nproc) as p:

        for i in range(nlight):

            filename = os.path.join(obspath, light_table['FILENAME'][i])
            xsc, ysc = light_table['xsc'][i], light_table['ysc'][i]
            xov, yov = light_table['xov'][i], light_table['yov'][i]

            args = (filename, savedir, masterdark, xsc, ysc, xov, yov, snrcut, fmax, xoff, yoff, T, bpix)

            p.apply_async(lightprocess_save, args=args, callback=lambda x: pbar.update())

        p.close()
        p.join()

    pbar.close()

    return

# if __name__ == '__Darksub__':
#     Darksub('SA-111', 'D:\\NEOSSat-SA-111','D:\\NEOSSat-SA-111\\clean')

