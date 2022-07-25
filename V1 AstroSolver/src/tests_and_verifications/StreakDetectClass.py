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

    def __init__(self, imageset, sigma, trm=True,
                 streaklength=10.0, **kwargs):
        """
        Initialize the class with the following parameters:
        """
        # for dirpath, dirnames, filenames in os.walk(imagedir):
        #    print(f'{len(filenames)} Files Detected')
        self.imageset = ""
        self.fitsdata = None
        self.solved_images = []
        self.failed_images = []
        self.detectedsatellites = {}
        self.sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
        self.bkg_estimator = SExtractorBackground()
        self.keepfiles = True
        self.magnitude = 2.5 #Tested "Best Guess" starting point for threshold magnitude
        self.file_suffix = {".fits", ".fit", ".FIT", '.fts'}
        self.sigma = sigma
        self.streakLength = 10

    def __str__(self):
        """
        Return a string representation of the object.
        Returns
        -------

        """

        #TODO Add this functionality
        return 'StreakDetect(Solved Images=' + str(len(self.solved_images)) +\
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
        #TODO Add this functionality


        return 'StreakDetect(Solved Images=' + str(len(self.solved_images)) +\
               ', Failed Images=' + str(len(self.failed_images)) + ')'


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
        threshold = bkg.background + (self.magnitude * bkg.background_rms)
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
            tbl = cat.to_table()
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
        if len(tbl_stars) == 0:
            print("No stars detected.")
        else:
            for i in range(len(tbl_stars)):
                print(f"Star detected at {tbl_stars['xcentroid'][i]:.3f}, "
                      f"{tbl_stars['ycentroid'][i]:.3f} with a flux of "
                      f"{tbl_stars['segment_flux'][i]} +/-{tbl_stars['segment_fluxerr'][i]} ")

        return self, tbl_stars, tbl_streaks

    def predict_params(self, **kwargs):
        return self

    def run(self):
        for filename in self.imageset:

            self.fitsdata = self.read_fitsdata(filename)
            self.mask = self.create_mask(self.fitsdata, sigma, streakLength)
            self.bkg = self.background_estimation(self.fitsdata, self.mask)
            self.threshold = self.generate_threshold(self.bkg, self.bkg.background_rms)
            self.cat = self.detect_sources(self.fitsdata, self.threshold, self.mask)
            while (len(self.cat) > 100 or len(self.cat) < 10) and solveattempts <= 5:

                if len(self.cat) > 100:
                    magnitude += 0.25
                elif len(self.cat) < 10:
                    magnitude -= 0.25
                solveattempts += 1
                self.threshold = self.generate_threshold(self.bkg, self.bkg.background_rms, magnitude)
                self.cat = self.detect_sources(self.fitsdata, self.threshold, self.mask)

            self.write_to_file(self.cat, self.STARS, self.AstrometryNetFile)
            self.solved_images.append(filename)
        if len(self.solved_images) > 0:
            print(f"{len(self.solved_images)} images solved.")
        else:
            print("No images solved.")

        return self.solved_images, self.failed_images




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

file_suffix = {".fits", ".fit", ".FIT", '.fts'}  # file suffixes
trm = True  # True if the file is a Track Rate Mode image
numFits = 0  # Number of images run through the Streak Detection

useMatchedFilter = False  # NOT WORKING - Experimental

# TODO Count fits not total files
# TODO Iterative matching based on results
# TODO Identify critical inputs
# TODO Convert to Class

# streak_detection(streak, file_suffix, trm, numFits, useMatchedFilter)