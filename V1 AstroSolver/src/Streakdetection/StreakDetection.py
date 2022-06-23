
import os
import numpy as np
import numpy
from scipy import fft
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.table import QTable
from astropy.stats import gaussian_fwhm_to_sigma
from photutils.background import Background2D
from photutils.background import SExtractorBackground
from photutils.segmentation import SourceCatalog
from photutils.segmentation import deblend_sources
from photutils.segmentation import detect_sources
from photutils.segmentation import detect_threshold
from photutils.segmentation import make_source_mask

streak1 = r'/Users/home/Sync/'
streak = r'/Users/home/Downloads/2022-108-neossat-crosstalk-example'
streak = r'/Users/home/Downloads/2020_J107_Ottawa_IS901'
file_suffix = (".fits", ".fit", ".FIT", '.fts')
trm = True
numFits = 0  # Number of images run through the Streak Detection

useMatchedFilter=False # NOT WORKING - Experimental

# TODO Count fits not total files
# TODO Iterative matching based on results
# TODO Identify critical inputs


def make_kernal_line(angle, length, option=None):

    #Generate a kernal line based on angle and length
    # angle=numpy.arctan(numpy.tan(angle))
    angle = numpy.deg2rad(angle)
    xsize=
    dx = int(numpy.ceil(length * max(np.abs(np.sin(angle)),
                                     np.cos(angle))))

    if dx % 2 == 0:
        dx += 1
    dy = dx
    linekernal = np.zeros((dx, dy), dtype=float)
    cx = (dx + 1) / 2
    cy = cx
    angle = angle.value
    if (angle > (0.7854) or angle < (-0.7854)):
        vertIncrement = np.cos(angle) / np.sin(angle)
        for j in range(dy):
            x = int(np.round(cx + (j - cy) * vertIncrement))

            linekernal[x][j] = 1
    else:
        horizIncrement = np.sin(angle) / np.cos(angle)
        for i in range(dx):
            y = int(np.round(cy + (i - cx) * horizIncrement))
            linekernal[i][y] = 1

    linekernal = linekernal / np.sum(np.sum(linekernal))
    return linekernal




def make_matched_filter(kernal, xsize, ysize):
    #Generate a matched filter based on angle and length
    #kernal = make_kernal_line(angle, length)   # Generate the kernal

    col = len(kernal)
    row = len(kernal[0])
    print(col)
    print(row)
    midrow = int(np.floor(row / 2))
    midcol = int(np.floor(col / 2))
    xsize = int(xsize)
    ysize = int(ysize)
    if xsize > col and ysize > row:
        filter = np.zeros((xsize, ysize), dtype=float)
        print(row - midrow)
        filter[int((xsize - row) / 2):int(((xsize + row) / 2) - 1)][
        int((ysize - col) / 2): int((ysize + col) / 2 - 1)] = kernal
        #filter[0:(row-midrow)][0:(col-midcol)] = kernal[(midrow+1):row][(midcol+1):col]
        # filter[1:(row - midrow)][(ysize-midcol+1):ysize] = kernal[midrow + 1:row][1:col]
        # filter[xsize - midrow+1:xsize][ysize-midcol+1:ysize] = kernal[1:midrow][1:midcol]
        # filter[xsize - midrow+1:xsize][1:col - midcol] = kernal[1:midrow][midcol + 1:col]
    ft = fft.fft2(filter)
    return ft
def write_pinpoint_file(STARS, tbl, astrometryNetFile):
    # Write a file with the results of the Streak Detection
    # filename = 'pinpoint_results.txt'

def


def streak_detection(imageDir, sigma=5.0, streakLength=5, TRM=True, useMask=True):
    fluxSat = []
    xlocSat = []
    ylocSat = []
    imageSats = []
    #satelliteDetections = QTable([imageSats, xlocSat, ylocSat, fluxSat],
    #                             names=('Image', 'X Location (pixel)', 'Y Location (pixel)', 'Instrumental Magnitude'),
    #                             meta={'name': 'Satellite Detections'})
    numFits = 0
    for dirpath, dirnames, filenames in os.walk(imageDir):
        print(f'{len(filenames)} Files Detected')
        for filename in filenames:
            if (filename.endswith(file_suffix)):
                numFits += 1
                filepath = os.path.join(dirpath,
                                        filename)

                # Create .STARS (pinpoint) and .txt(astrometry) detection files
                STARS = open(filepath + '.stars', "w")
                AstrometryNetFile = open(filepath + '.txt', 'w')

                fitsdata = fits.open(filepath)[0].data
                sigma_clip = SigmaClip(sigma, maxiters=10)
                mask = make_source_mask(fitsdata, nsigma=sigma, npixels=streakLength,
                                        dilate_size=11)
                bkg_estimator = SExtractorBackground()

                if trm:
                    try:
                        bkg = Background2D(fitsdata, (30, 30),filter_size=(3,3),mask=mask,
                                           bkg_estimator=bkg_estimator)
                    except:
                        bkg = SExtractorBackground(fitsdata)
                    threshold = bkg.background + (2.0 * bkg.background_rms)
                    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
                    convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
                    segm = detect_sources(fitsdata, threshold, npixels=streakLength)
                    segm_deblend = deblend_sources(fitsdata, segm, npixels=streakLength,
                                                   nlevels=50, contrast=0.0001)
                    cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=None)

                else:
                    mask = make_source_mask(fitsdata, nsigma=sigma, npixels=5, dilate_size=11)
                    bkg = Background2D(fitsdata, (30, 30), filter_size=(3, 3),
                                       mask=mask, bkg_estimator=bkg_estimator)
                    threshold = bkg.background + (2.5 * bkg.background_rms)
                    #sigma = 5.5 * gaussian_fwhm_to_sigma  # FWHM = 3.
                    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
                    convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
                    segm = detect_sources(convolved_data, threshold, npixels=10)
                    segm_deblend = deblend_sources(convolved_data, segm, npixels=10,
                                                   nlevels=32, contrast=0.001)
                    cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=convolved_data)

                tbl = cat.to_table()
                #print(str(len(tbl)) + "Sources Detected")

                if useMatchedFilter==True:
                    xmax = cat.bbox_xmax
                    xmin = cat.bbox_xmin
                    ymax = cat.bbox_ymax
                    ymin = cat.bbox_ymin
                    cenX = cat.xcentroid
                    cenY = cat.ycentroid
                    length = numpy.zeros((len(cenY), 1))
                    angles = cat.orientation
                    medangle = numpy.median(angles)
                    xsiz = len(fitsdata)
                    ysiz = len(fitsdata[0])
                    #print(angles)
                    for i in range(len(xmax)): #Calculate Length of streaks
                        x_length = xmax[i] - xmin[i]
                        y_length = ymax[i] - ymin[i]
                        length[i] = numpy.sqrt(x_length ** 2 + y_length ** 2)
                    #print(length)
                    medlen=numpy.median(length)

                    kernel=make_kernal_line(medangle,medlen) #Create a matched filter
                    #filter=make_matched_filter(kernel,xsiz,ysiz)
                    convolved_data = convolve_fft(fitsdata, kernel, normalize_kernel=True)
                    segm = detect_sources(convolved_data, threshold, npixels=10)
                    segm_deblend = deblend_sources(convolved_data, segm, npixels=10,
                                                      nlevels=32, contrast=0.001)
                    cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=convolved_data)
                    tbl=cat.to_table()
                    #print(str(len(tbl)) + "Sources Detected - with Matched Filter")
                    #print(cat.orientation)

                for i in range(len(tbl['segment_flux'])):
                    streak_line = '{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(tbl['xcentroid'][i]),
                                                                                  float(tbl['ycentroid'][i]),
                                                                                  tbl['segment_flux'][i])
                    STARS.write(streak_line + "\n")
                    if cat.eccentricity[i] < 0.5:
                        #satelliteDetections.append(filename, float(tbl['xcentroid'][i]), float(tbl['ycentroid'][i]),
                        #                           float(tbl['segment_flux'][i]))
                        print(
                            f"Potential satellite detected at {cat.xcentroid[i]:.3f}, "
                            f"{cat.ycentroid[i]:.3f} with a flux of "
                            f"{cat.segment_flux[i]} +/-{cat.segment_fluxerr[i]} ")

                tbl.sort('segment_flux', reverse=True)
                count = 0
                for i in range(len(tbl['xcentroid'])):
                    streak_line2 = '{:.4f} {:.4f}'.format(float(tbl['xcentroid'][i]),
                                                          float(tbl['ycentroid'][i]))
                    AstrometryNetFile.write(streak_line2 + "\n")
                    count += 1
                print(f'{filename} index created.')
                print(f" {len(tbl['xcentroid'])} streaks detected.")

                STARS.close()
                AstrometryNetFile.close()

    if numFits == 0:
        print("No Valid Images Detected")
        return None

    elif trm == True:
        #if len(satelliteDetections) == 0:
            #return "No valid satellite detections - Check images or adjust parameters. "
        #else:
            #return satelliteDetections
        return 0
    else:
        #return satelliteDetections
        return 0
