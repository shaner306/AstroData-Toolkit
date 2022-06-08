
import matplotlib.pyplot as plt
#import pywin32_system32
#import win32com
import os
import numpy
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import SqrtStretch
from astropy.visualization import simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.background import Background2D
from photutils.background import SExtractorBackground
from photutils.segmentation import SourceCatalog
from photutils.segmentation import deblend_sources
from photutils.segmentation import detect_sources
from photutils.segmentation import detect_threshold
from photutils.segmentation import make_source_mask
streak =r'/Users/home/Sync/'
streak1= r'/Users/home/Downloads/2021 10 21 - ZWO with C8/Test'

file_suffix = (".fits", ".fit", ".FIT")
trm=1


def StreakDetection(streak, TRM):

    for dirpath, dirnames, filenames in os.walk(streak):
        print(f'{len(filenames)} Images Detected')
        for filename in filenames:
            if (filename.endswith(file_suffix)):
                #sum1+=1
                filepath = os.path.join(dirpath, filename)
                STARS = open(filepath + '.stars', "w")
                AstrometryNetFile=open(filepath+'.txt', 'w')

                imagehdularray = fits.open(filepath)
                fitsdata = imagehdularray[0].data

                sigma_clip = SigmaClip(sigma=3)
                mask = make_source_mask(fitsdata, nsigma=2, npixels=5, dilate_size=11)
                bkg_estimator = SExtractorBackground()
                if trm==1:
                    bkg = Background2D(fitsdata, (30, 30), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
                    # fitsdata = fitsdata - bkg.background
                    # threshold = detect_threshold(fitsdata, nsigma=3.0)

                    threshold = bkg.background + (1.5 * bkg.background_rms)
                    sigma = 2.5 * gaussian_fwhm_to_sigma  # FWHM = 3.
                    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
                    convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
                    segm = detect_sources(fitsdata, threshold, npixels=20)
                    segm_deblend = deblend_sources(fitsdata, segm, npixels=20,
                                                   nlevels=50, contrast=0.0001)
                else:
                    bkg = Background2D(fitsdata, (30, 30), filter_size=(3,3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
                    #fitsdata = fitsdata - bkg.background
                    #threshold = detect_threshold(fitsdata, nsigma=3.0)

                    threshold = bkg.background + (2.5 * bkg.background_rms)
                    sigma = 2.5 * gaussian_fwhm_to_sigma  # FWHM = 3.
                    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
                    convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
                    segm = detect_sources(convolved_data, threshold, npixels=5)
                    segm_deblend = deblend_sources(convolved_data, segm, npixels=5,
                                               nlevels=32, contrast=0.001)

                cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=None)

                tbl = cat.to_table()


                newCenY = list(tbl['ycentroid'])
                newCenX = list(tbl['xcentroid'])
                newflux = list(tbl['segment_flux'])


                for i in range(len(newCenY)):
                    streak_line = '{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(newCenX[i]), float(newCenY[i]), newflux[i])
                    STARS.write(streak_line + "\n")

                tbl.sort('segment_flux', reverse=True)
                newCenY = list(tbl['ycentroid'])
                newCenX = list(tbl['xcentroid'])
                newflux = list(tbl['segment_flux'])

                for i in range(len(newCenY)):
                    streak_line2 = '{:.4f} {:.4f}'.format(float(newCenX[i]), float(newCenY[i]))
                    AstrometryNetFile.write(streak_line2 + "\n")
                print(f'{filename} index created.')
                print(f' {len(newCenY)} streaks detected.')

                STARS.close()
                AstrometryNetFile.close()
    #print(f'{sum1} streak detection files generated')

