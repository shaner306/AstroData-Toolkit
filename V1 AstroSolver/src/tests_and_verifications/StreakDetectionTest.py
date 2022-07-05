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
from photutils.segmentation import make_source_mask
from scipy.ndimage import rotate
from matplotlib import pyplot as plt
import os
streak1 = r'/Users/home/Sync/'
streak2 = r'/Users/home/Downloads/2022-108-neossat-crosstalk-example'
streak = r'/Users/home/Downloads/2020_J107_Ottawa_IS901'

class DetectStreaks():

    def __init__(self):
        #for dirpath, dirnames, filenames in os.walk(imagedir):
        #    print(f'{len(filenames)} Files Detected')
        self.imageset = ""
        self.fitsdata=None
        self.solved_images=[]
        self.failed_images=[]
        self.detectedsatellites={}
        self.sigma_clip=SigmaClip(sigma=3.0,maxiters=10)
        self.bkg_estimator = SExtractorBackground()
        self.keepfiles=True

    def read_fitsdata(self, filename):
        try:
            fitsdata = fits.open(filename)[0].data
        except:
            print("File Error - No Image Data Detected.")
        return fitsdata

    def create_stars_file(self, filepath):
        STARS = open(filepath + '.stars', "w")
        AstrometryNetFile = open(filepath + '.txt', 'w')
        return STARS, AstrometryNetFile

    def delete_stars_file(self, filepath):
        os.remove(filepath + '.stars')
        os.remove(filepath + '.txt')
        return self

    def create_mask(self, fitsdata, sigma, streak_length, **kwargs):
        mask = make_source_mask(fitsdata, nsigma=sigma, npixels=streakLength,
                                dilate_size=11, **kwargs)
        return self, mask

    def background_estimation(self, fitsdata, mask, **kwargs):
        try:
            bkg = Background2D(fitsdata, (30, 30), filter_size=(3, 3), mask=mask,
                               bkg_estimator=self.bkg_estimator)
        except:
            bkg = SExtractorBackground(fitsdata)
        return self, bkg

    def generate_threshold(self, bkg, magnitude):
        threshold = bkg.background + (magnitude * bkg.background_rms)
        return self, threshold

    def detect_sources(self, fitsdata, threshold, mask, **kwargs):

        if convolve_data:
            kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
            convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
            segm = detect_sources(fitsdata, threshold, npixels=streakLength)
            segm_deblend = deblend_sources(fitsdata, segm, npixels=streakLength,
                                           nlevels=50, contrast=0.0001)
            cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=None)
        else:
            kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
            convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
            segm = detect_sources(fitsdata, threshold, npixels=streakLength)
            segm_deblend = deblend_sources(fitsdata, segm, npixels=streakLength,
                                           nlevels=50, contrast=0.0001)
            cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=None)
            tbl=cat.to_table()
        return self, tbl

    def write_to_file(self, tbl, STARS, AstrometryNetFile):
        for i in range(len(tbl['segment_flux'])):
            streak_line = '{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(tbl['xcentroid'][i]),
                                                                          float(tbl['ycentroid'][i]),
                                                                          tbl['segment_flux'][i])
            STARS.write(streak_line + "\n")
            if cat.eccentricity[i] < 0.5:
                # satelliteDetections.append(filename, float(tbl['xcentroid'][i]), float(tbl['ycentroid'][i]),
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

        return self

    def separate_streaks_and_stars(self, tbl, **kwargs):
        """
        Separate the detected sources into stars and streaks based on
        eccentricity of the detected source, if the source is a circle it is a satellite,
        otherwise the streaks represent stars.
        Parameters
        ----------
        tbl - SourceCatalog object
            Combined table of all sources detected in the image.

        kwargs

        Returns
        -------
        tbl_stars - SourceCatalog object
            Table of stars detected in the image.
        tbl_streaks - SourceCatalog object
            Table of streaks detected in the image.

        """
        tbl_stars = tbl[tbl['eccentricity'] < 0.5]
        tbl_streaks = tbl[tbl['eccentricity'] >= 0.5]

        return self, tbl_stars, tbl_streaks

    def predict_params(self, **kwargs):
        return self

    def run(self):
        for filename in self.imageset:
            self.fitsdata = self.read_fitsdata(filename)
            self.mask = self.create_mask(self.fitsdata, sigma, streakLength)
            self.bkg = self.background_estimation(self.fitsdata, self.mask)
            self.threshold = self.generate_threshold(self.bkg, self.bkg.background_rms, 2.5)
            self.cat = self.detect_sources(self.fitsdata, self.threshold, self.mask)
            while (len(self.cat) > 100 or len(self.cat) < 10) and solveattempts <=5:

                if len(self.cat) > 100:
                    magnitude += 0.25
                elif len(self.cat) < 10:
                    magnitude -= 0.25
                solveattempts += 1
                self.threshold = self.generate_threshold(self.bkg, self.bkg.background_rms, magnitude)
                self.cat = self.detect_sources(self.fitsdata, self.threshold, self.mask)

            self.write_to_file(self.cat, self.STARS, self.AstrometryNetFile)
            self.solved_images.append(filename)
        return self.solved_images, self.failed_images

    def make_kernel_line(self, angle, length, **kwargs):
        """
            Generate a matched filter kernel based on angle and length of a star streak
            Parameters
            ----------
            angle: float radians
                Median angle of the image streaks in radians
            length : float
                Median length of the image streaks


            Returns
            -------
            linekernel : numpy.ndarray
                Convolution kernel
        """
        angle = numpy.deg2rad(angle)
        angle = np.arctan(np.tan(angle))

        xsize = int(numpy.ceil(length * max(np.abs(np.sin(angle)), np.cos(angle))))
        xsize = xsize if xsize % 2 == 1 else xsize + 1
        ysize = xsize
        linekernel = np.zeros((xsize, ysize))
        cx = (xsize + 1) / 2
        cy = cx

        print(angle)
        if angle > (np.pi / 4) or angle < (-np.pi / 4):
            vertIncrement = np.cos(angle) / np.sin(angle)
            for j in range(1, ysize - 1):
                x = int(np.round(cx + (j - cy) * vertIncrement))

                linekernel[x][j] = 1
        else:
            horizIncrement = np.sin(angle) / np.cos(angle)
            for i in range(1, xsize - 1):
                y = int(np.round(cy + (i - cx) * horizIncrement))
                linekernel[i][y] = 1

        linekernel = linekernel / (np.sum(np.sum(linekernel)))
        return linekernel

    def expand_matched_filter(self,image, linekernel):
        """
        Enlarge the kernal shape (matched filter) to the same size of the image to be filtered
        (this is an 'enlarge' resize function) and calculate its Fourier transform.

        The end result is a ready to use matched filter (ready to be multiplied by the
        Fourier transform of the image to be filtered).

        Parameters
        ----------
        image : numpy.ndarray
            Image to be filtered - Calculate size
        linekernel : numpy.ndarray
            Convolution kernel to be used for the matched filter

        Returns
        -------
        matched_filter : numpy.ndarray
            Matched filter ready to be multiplied by the Fourier transform of the image to be filtered
        """

        return matched_filter


def make_kernel_line(angle, length, **kwargs):
    """
        Generate a matched filter kernel based on angle and length of a star streak
        Parameters
        ----------
        angle: float radians
            Median angle of the image streaks in radians
        length : float
            Median length of the image streaks


        Returns
        -------
        linekernel : numpy.ndarray
            Convolution kernel
    """
    angle = numpy.deg2rad(angle)
    # angle=angle.value
    angle = np.arctan(np.tan(angle))

    xsize = int(numpy.ceil(length * max(np.abs(np.sin(angle)), np.cos(angle))))
    xsize = xsize if xsize % 2 == 1 else xsize + 1
    ysize = xsize
    linekernel = np.zeros((xsize, ysize))
    cx = (xsize + 1) / 2
    cy = cx

    print(angle)
    if angle > (np.pi / 4) or angle < (-np.pi / 4):
        vertIncrement = np.cos(angle) / np.sin(angle)
        for j in range(1,ysize-1):
            x = int(np.round(cx + (j - cy) * vertIncrement))
            print(x)
            linekernel[x][j] = 1
    else:
        horizIncrement = np.sin(angle) / np.cos(angle)
        for i in range(1,xsize-1):
            y = int(np.round(cy + (i - cx) * horizIncrement))
            print(y)
            linekernel[i][y] = 1

    linekernel = linekernel / (np.sum(np.sum(linekernel)))
    return linekernel


def make_matched_filter(kernal, xsize, ysize):
    # Generate a matched filter based on angle and length
    # kernal = make_kernal_line(angle, length)   # Generate the kernal

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

        # filter[0:(row-midrow)][0:(col-midcol)] = kernal[(midrow+1):row][(midcol+1):col]
        # filter[1:(row - midrow)][(ysize-midcol+1):ysize] = kernal[midrow + 1:row][1:col]
        # filter[xsize - midrow+1:xsize][ysize-midcol+1:ysize] = kernal[1:midrow][1:midcol]
        # filter[xsize - midrow+1:xsize][1:col - midcol] = kernal[1:midrow][midcol + 1:col]
    ft = fft.fft2(filter)
    return ft


# def write_pinpoint_file(STARS, tbl, astrometryNetFile):
# Write a file with the results of the Streak Detection
# filename = 'pinpoint_results.txt'

def streak_detection(imageDir, sigma=5.0, streakLength=5, TRM=True, useMask=True, **kwargs):
    """
    Parameters
    ----------
    imageDir
    sigma
    streakLength
    TRM
    useMask

    Returns
    -------

    """
    fluxSat = []
    xlocSat = []
    ylocSat = []
    imageSats = []
    # satelliteDetections = QTable([imageSats, xlocSat, ylocSat, fluxSat],
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
                sigma_clip = SigmaClip(sigma, maxiters=5)
                mask = make_source_mask(fitsdata, nsigma=sigma, npixels=streakLength,
                                        dilate_size=11)
                bkg_estimator = SExtractorBackground()

                if trm:
                    try:
                        bkg = Background2D(fitsdata, (30, 30), filter_size=(3, 3), mask=mask,
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
                    # sigma = 5.5 * gaussian_fwhm_to_sigma  # FWHM = 3.
                    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
                    convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
                    segm = detect_sources(convolved_data, threshold, npixels=10)
                    segm_deblend = deblend_sources(convolved_data, segm, npixels=10,
                                                   nlevels=32, contrast=0.001)
                    cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=convolved_data)

                tbl = cat.to_table()
                # print(str(len(tbl)) + "Sources Detected")

                if useMatchedFilter == True:
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
                    # print(angles)
                    for i in range(len(xmax)):  # Calculate Length of streaks
                        x_length = xmax[i] - xmin[i]
                        y_length = ymax[i] - ymin[i]
                        length[i] = numpy.sqrt(x_length ** 2 + y_length ** 2)
                    # print(length)
                    medlen = numpy.median(length)

                    kernel = make_kernal_line(medangle, medlen)  # Create a matched filter
                    # filter=make_matched_filter(kernel,xsiz,ysiz)
                    convolved_data = convolve_fft(fitsdata, kernel, normalize_kernel=True)
                    segm = detect_sources(convolved_data, threshold, npixels=10)
                    segm_deblend = deblend_sources(convolved_data, segm, npixels=10,
                                                   nlevels=32, contrast=0.001)
                    cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=convolved_data)
                    tbl = cat.to_table()
                    # print(str(len(tbl)) + "Sources Detected - with Matched Filter")
                    # print(cat.orientation)

                for i in range(len(tbl['segment_flux'])):
                    streak_line = '{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(tbl['xcentroid'][i]),
                                                                                  float(tbl['ycentroid'][i]),
                                                                                  tbl['segment_flux'][i])
                    STARS.write(streak_line + "\n")
                    if cat.eccentricity[i] < 0.5:
                        # satelliteDetections.append(filename, float(tbl['xcentroid'][i]), float(tbl['ycentroid'][i]),
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
        # if len(satelliteDetections) == 0:
        # return "No valid satellite detections - Check images or adjust parameters. "
        # else:
        # return satelliteDetections
        return 0
    else:
        # return satelliteDetections
        return 0


streak1 = r'/Users/home/Sync/'
streak2 = r'/Users/home/Downloads/2022-108-neossat-crosstalk-example'
streak = r'/Users/home/Downloads/2020_J107_Ottawa_IS901'

file_suffix = {".fits", ".fit", ".FIT", '.fts'} # file suffixes
trm = True  # True if the file is a Track Rate Mode image
numFits = 0  # Number of images run through the Streak Detection

useMatchedFilter = False   # NOT WORKING - Experimental

# TODO Count fits not total files
# TODO Iterative matching based on results
# TODO Identify critical inputs
# TODO Convert to Class

#streak_detection(streak, file_suffix, trm, numFits, useMatchedFilter)