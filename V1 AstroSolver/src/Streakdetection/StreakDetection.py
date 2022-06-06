
import matplotlib.pyplot as plt
#import pywin32_system32
#import win32com
import os
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

streak = 'D:\\Breeze-M_R_B_38746U'

file_suffix = (".fits", ".fit", ".fts")

for dirpath, dirnames, filenames in os.walk(streak):
    for filename in filenames:
        if (filename.endswith(file_suffix)):
            filepath = os.path.join(dirpath, filename)
            STARS = open(filepath + '.stars', "w")

            imagehdularray = fits.open(filepath)
            fitsdata = imagehdularray[0].data

            sigma_clip = SigmaClip(sigma=3)
            bkg_estimator = SExtractorBackground()
            bkg = Background2D(fitsdata, (30, 30), filter_size=(3, 3),
                               sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
            bg_rem = fitsdata - bkg.background
            threshold = detect_threshold(fitsdata, nsigma=3.0)

            sigma = 2.5 * gaussian_fwhm_to_sigma  # FWHM = 3.
            kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
            convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
            segm = detect_sources(fitsdata, threshold, npixels=5)
            segm_deblend = deblend_sources(fitsdata, segm, npixels=5,
                                       nlevels=32, contrast=0.001)

            cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=convolved_data)
            tbl = cat.to_table()

            newCenY = list(tbl['ycentroid'])
            newCenX = list(tbl['xcentroid'])
            newflux = list(tbl['segment_flux'])

            for i in range(len(newCenY)):
                streak_line = '{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'\
                    .format(float(newCenX[i]), float(newCenY[i]), newflux[i])
                STARS.write(streak_line + "\n")



            STARS.close()