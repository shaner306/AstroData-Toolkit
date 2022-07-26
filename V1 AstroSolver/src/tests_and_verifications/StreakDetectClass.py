import os

import matplotlib.pyplot as plt
import numpy
import numpy as np
import scipy
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.io import fits
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.table import QTable
from photutils.utils import circular_footprint
from photutils import detect_threshold
from photutils.background import Background2D
from photutils.background import SExtractorBackground
from photutils.segmentation import SourceCatalog
from photutils.segmentation import deblend_sources
from photutils.segmentation import detect_sources
from photutils.segmentation import make_source_mask
from scipy import signal
from photutils.detection import DAOStarFinder
import pandas as pd
from astropy.io import ascii
from astropy.table import Table
from scipy import misc

streak22 = r'/Users/home/Downloads/Landolt Fields/New Folder With Items'
streak = r'/Users/home/Downloads/KEPLER2022-166to177/166'
streak2 = r"/Users/home/Downloads/2022-108-neossat-crosstalk-example"
streak1 = r'/Users/home/Downloads/2020_J107_Ottawa_IS901'


class DetectStreaks:

    def __init__(self, imageset, sigma, trm=True,
                 streak_length=5.0, **kwargs):
        """
        Initialize the class with the following parameters:
        """
        self.imageset = imageset
        self.fitsdata = None
        self.solved_images = []
        self.failed_images = []
        self.detectedsatellites = QTable()  # Table of detected satellites
        self.sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
        self.bkg_estimator = SExtractorBackground()
        self.keepfiles = True
        self.magnitude = 2.0  # Tested "Best Guess" starting point for threshold magnitude
        self.file_suffix = (".fits", ".fit", ".FIT", ".fts")
        self.streakLength = streak_length
        self.trm = trm
        self.sigma = sigma

    def __str__(self):
        """
        Return a string representation of the object.
        Returns
        -------

        """

        # TODO Add this functionality
        return 'StreakDetect(Solved Images=' + str(len(self.solved_images)) + \
               ', Failed Images=' + str(len(self.failed_images)) + ')'

    def __repr__(self):
        """
        Return a string representation of the object,
        for solve parameters, and solve results for each image, and relevant statistics.
        Returns
        -------
        string:
            List of Solved Images
                Number of streaks detected
                Number of stars detected
                Solve Params
                Image Statistics

            List of Failed Images
                Last used Params
                Image Statistics
        """
        # TODO Add this functionality
        return 'StreakDetect(Solved Images=' + str(len(self.solved_images)) + \
               ', Failed Images=' + str(len(self.failed_images)) + ')'

    def create_stars_file(self, filepath):
        STARS = open(filepath + '.stars', "w")
        AstrometryNetFile = open(filepath + '.txt', 'w')
        return STARS, AstrometryNetFile

    def read_fitsdata(self, filename):
        """
        Read in the FITS file and return the data.

        Parameters
        ----------
        filename - string:
            The name of the FITS file to read.

        Returns
        -------
        fitsdata - numpy.ndarray:
            The data from the FITS file.
        """
        fitsdata = 0
        try:
            fitsdata = fits.open(filename)[0].data
        except:
            print("File Error - No Image Data Detected.")
        return fitsdata

    def delete_stars_file(self, filepath):
        os.remove(filepath + '.stars')
        os.remove(filepath + '.txt')
        return self

    def create_mask(self, fitsdata, **kwargs):
        '''
        Create masks for both the streaks and the circular point sources on the image.
        '''


        #mask = make_source_mask(fitsdata, nsigma=self.sigma, npixels=self.streakLength,
        #                        dilate_size=11, **kwargs)
        sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
        threshold = detect_threshold(fitsdata, nsigma=2.0, sigma_clip=sigma_clip)

        segment_img = detect_sources(fitsdata, threshold, npixels=5)
        footprint = circular_footprint(radius=5)
        mask = segment_img.make_source_mask(footprint=footprint)
        #mean, median, std = sigma_clipped_stats(fitsdata, sigma=3.0, mask=mask)


        return mask

    def background_estimation(self, fitsdata, mask, **kwargs):
        try:
            bkg = Background2D(fitsdata, (30, 30), filter_size=(3, 3), mask=mask,
                               bkg_estimator=self.bkg_estimator)
        except:
            bkg = SExtractorBackground(fitsdata)
        return bkg

    def generate_threshold(self, bkg, **kwargs):
        threshold = bkg.background + (self.magnitude * bkg.background_rms)
        return threshold

    def predict_angle_and_streak_length(self, fitsdata, **kwargs):
        '''
        Predict the angle and streak length of the streak.
        Parameters
        ----------
        fitsdata : numpy.ndarray
            Image data to be processed.
        kwargs
        Returns
        -------
        angle : float
            Angle of the streak.
        streak_length : float
            Length of the streak.
        '''
        sigma = 5

        sigma_clip = SigmaClip(sigma, maxiters=5)
        mask = make_source_mask(fitsdata, nsigma=sigma, npixels=5.0,
                                dilate_size=11)
        bkg_estimator = SExtractorBackground()
        try:
            bkg = Background2D(fitsdata, (30, 30), filter_size=(3, 3), mask=mask,
                               bkg_estimator=bkg_estimator)
        except:
            bkg = SExtractorBackground(fitsdata)
        threshold = bkg.background + (2.0 * bkg.background_rms)
        kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
        convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
        segm = detect_sources(fitsdata, threshold, npixels=5.0)
        segm_deblend = deblend_sources(fitsdata, segm, npixels=5.0,
                                       nlevels=50, contrast=0.001)
        cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=None)
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
        self.streakLength = int(np.round(medlen, 0))

    def remove_saturated_pixels(self, fitsdata, **kwargs):
        """
        Remove saturated pixels from the image, if the pixel value is above 65535.
        Parameters
        ----------
        fitsdata - numpy.ndarray:
            The data to be cleaned.
        Returns
        -------
        fitsdata - numpy.ndarray:
            The image data with saturated pixels removed.

        """

        fitsdata[fitsdata >= 65535] = 0
        return fitsdata

    def detect_sources(self, fitsdata, threshold, **kwargs):

        if not self.trm:
            kernel = Gaussian2DKernel(self.sigma, x_size=3, y_size=3)
            convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
            segm = detect_sources(fitsdata, threshold, npixels=self.streakLength)
            segm_deblend = deblend_sources(fitsdata, segm, npixels=self.streakLength,
                                           nlevels=50, contrast=0.0001)
            cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=None)
        else:
            kernel = Gaussian2DKernel(self.sigma, x_size=3, y_size=3)
            convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
            segmstars = detect_sources(fitsdata, threshold, npixels=self.streakLength)
            segmstreaks= detect_sources(fitsdata, threshold, npixels=15)
            segm_deblend = deblend_sources(fitsdata, segmstreaks, npixels=15,
                                           nlevels=50, contrast=0.001)
            plt.imshow(segm_deblend, cmap='gray')
            plt.show()
            cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=None)

        tbl = cat.to_table()
        return cat, tbl

    def write_to_file(self, filename, filepath, cat, tbl, STARS, AstrometryNetFile):

        tbl.write(filepath+'streakvalues.ecsv', overwrite=True)
        for i in range(len(tbl['segment_flux'])):
            streak_line = '{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'. \
                format(float(tbl['xcentroid'][i]),
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
        print(f"{filename} index created.")
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
        if len(tbl_stars) == 0:
            print("No stars detected.")
        else:
            for i in range(len(tbl_stars)):
                print(f"Star detected at {tbl_stars['xcentroid'][i]:.3f}, "
                      f"{tbl_stars['ycentroid'][i]:.3f} with a flux of "
                      f"{tbl_stars['segment_flux'][i]} +/-{tbl_stars['segment_fluxerr'][i]} ")

        return tbl_stars, tbl_streaks

    def killBrightSpots(self, fitsdata, sigma=3.0): #NOT USED AT THE MOMENT
        '''
        Remove all pixels with a brightness greater than sigma*sigma*mean(data)
        Parameters
        ----------
        data
        sigma

        Returns
        -------

        '''

        x, y = fitsdata.shape
        spots = numpy.zeros((x, y))
        imgout = fitsdata
        outwindowindex = [-7, -4, -1, 1, 4, 7]
        k3 = numpy.ones((3, 3))
        k3 = k3 / 9
        noise = sigma * 10
        means = scipy.signal.convolve2d(fitsdata, k3)
        for i in range(9, x - 8):
            for j in range(9, y - 8):
                if means[i, j] > noise:
                    out = means[i + outwindowindex[0]:i + outwindowindex[5],
                          j + outwindowindex[0]:j + outwindowindex[5]]
                    out[2:5, 2:5] = 0
                    maxmeanout = numpy.mean(out)
                    if (5 * maxmeanout) < means[i, j]:
                        spots[i - 3:i + 3, j - 3:j + 3] = fitsdata[i - 3:i + 3, j - 3:j + 3]
                        imgout[i - 3:i + 3, j - 3:j + 3] = 0
        return imgout, spots

    def remove_dark_column(self, fitsdata, **kwargs):
        """
        Remove dark columns from the image. Dark columns are vertical strips on the image
        where the average pixel value is zero prior to background subtraction.
        They only occur at the edge of the image, and once detected, clip them from the image.
        Parameters
        ----------
        fitsdata : numpy.ndarray
            Image data to be processed.
        kwargs
        Returns
        -------
        fitsdata : numpy.ndarray
            Image data with dark columns removed.
        """

        # clip 75 pixel width column from left of the image
        fitsdata = fitsdata[:, 75:]
        return fitsdata

    def run(self):
        """
        Run the Streak Detection pipeline.

        Parameters
        ----------
        self

        Returns
        -------
        Array of Failed and Solved Images

        """
        numFits = 0
        solveattempts = 0
        for dirpath, dirnames, filenames in os.walk(self.imageset):
            print(f'{len(filenames)} Files Detected')
            for filename in filenames:
                if filename.endswith((".fits", ".fit", ".FIT", ".fts")):
                    print(filename)
                    numFits += 1
                    filepath = os.path.join(dirpath, filename)
                    STARS, AstrometryNetFile = self.create_stars_file(filepath)
                    fitsdata = self.read_fitsdata(filepath)
                    # fitsdata, spots = self.killBrightSpots(fitsdata)
                    fitsdata = self.remove_dark_column(fitsdata)
                    # self.predict_angle_and_streak_length(fitsdata)

                    mask = self.create_mask(fitsdata)
                    bkg = self.background_estimation(fitsdata, mask)
                    threshold = self.generate_threshold(bkg)
                    cat, tbl = self.detect_sources(fitsdata, threshold)
                    while ((len(cat) >= 100) or (len(cat) <= 10)) and (solveattempts <= 5):
                        if len(cat) >= 100:
                            self.magnitude += 0.25
                        elif len(cat) <= 10:
                            self.magnitude -= 0.25
                        solveattempts += 1
                        threshold = self.generate_threshold(bkg)
                        cat, tbl = self.detect_sources(fitsdata, threshold)

                    self.write_to_file(filename,filepath, cat, tbl, STARS, AstrometryNetFile)
                    self.solved_images.append(filename)

        if len(self.solved_images) > 0:
            print(f"{len(self.solved_images)} images solved.")
        else:
            print("No images solved.")

        return self.solved_images, self.failed_images


# TODO: Print Full table of detected streaks and info
# TODO: Add in option to detect stars and streaks separately
# TODO: Remove Saturated pixels
# TODO: Remove faint stars
# TODO: Remove anything that doesnt have similar moments to the average streak
# Todo Add in matched filtering

d = DetectStreaks(streak, 5.0)
d.run()
