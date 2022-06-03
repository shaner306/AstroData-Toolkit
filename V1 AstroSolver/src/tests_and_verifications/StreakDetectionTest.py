import matplotlib.pyplot as plt
# import pywin32_system32
# import win32com
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

streak122 = r'D:\Transfer to mac\2021-03-10 - Calibrated\Intelsat 10-02 Post Eclipse\LIGHT\B_lim\0066_3x3_-10.00_5.00_B_21-23-04.fits'
streak = 'D:\\Breeze-M_R_B_38746U'
streak2 = r'/Users/home/Downloads/CAN_OTT.00041100.22108.FIT'
streak1 = r'D:\Transfer to mac\trm-stars-images\NEOS_SCI_2021099173229frame.fits'
streak13 = r'D:\Solved Stars\Tycho 3023_1724\LIGHT\B\0000_3x3_-10.00_5.00_B_21-22-59.fits'
file_suffix = (".fits", ".fit", ".fts")

for dirpath, dirnames, filenames in os.walk(streak):
    for filename in filenames:
        if (filename.endswith(file_suffix)):
            filepath = os.path.join(dirpath, filename)
            STARS = open(filepath + '.stars', "w")

            imagehdularray = fits.open(streak2)
            date = imagehdularray[0].header['DATE-OBS']
            exposuretime = imagehdularray[0].header['EXPTIME']
            imagesizeX = imagehdularray[0].header['NAXIS1']
            imagesizeY = imagehdularray[0].header['NAXIS2']
            fitsdata = imagehdularray[0].data

            sigma_clip = SigmaClip(sigma=3)
            bkg_estimator = SExtractorBackground()
            bkg = Background2D(fitsdata, (30, 30), filter_size=(3, 3), sigma_clip=sigma_clip,
                               bkg_estimator=bkg_estimator)
            bg_rem = fitsdata - bkg.background
            threshold = detect_threshold(fitsdata, nsigma=3.0)

            plt.imshow(fitsdata - bkg.background, origin='lower',
                       cmap='Greys_r', interpolation='nearest')

            sigma = 2.5 * gaussian_fwhm_to_sigma  # FWHM = 3.
            kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
            convolved_data = convolve(fitsdata, kernel, normalize_kernel=True)
            segm = detect_sources(fitsdata, threshold, npixels=5)
            segm_deblend = deblend_sources(fitsdata, segm, npixels=5,
                                           nlevels=32, contrast=0.001)

            norm = ImageNormalize(stretch=SqrtStretch())

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
            ax1.imshow(fitsdata, origin='lower', cmap='Greys_r', norm=norm)
            ax1.set_title('Data')
            cmap = segm.make_cmap(seed=123)
            ax2.imshow(segm, origin='lower', cmap=cmap, interpolation='nearest')
            ax2.set_title('Segmentation Image')

            cat = SourceCatalog(fitsdata, segm_deblend, convolved_data=convolved_data)
            tbl = cat.to_table()

            tbl['xcentroid'].info.format = '.2f'  # optional format
            tbl['ycentroid'].info.format = '.2f'
            tbl['kron_flux'].info.format = '.2f'
            print(tbl)

            # norm = simple_norm(fitsdata, segm, 'sqrt')
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
            ax1.imshow(fitsdata, origin='lower', cmap='Greys_r', norm=norm)
            ax1.set_title('FitsData')
            cmap = segm_deblend.make_cmap(seed=123)
            ax2.imshow(segm_deblend, origin='lower', cmap=cmap,
                       interpolation='nearest')
            ax2.set_title('Segmentation Image')
            cat.plot_kron_apertures((2.5, 1.0), axes=ax1, color='white', lw=1.5)
            cat.plot_kron_apertures((2.5, 1.0), axes=ax2, color='white', lw=1.5)

            newCenY = list(tbl['ycentroid'])
            newCenX = list(tbl['xcentroid'])
            newflux = list(tbl['segment_flux'])

            for i in range(len(newCenY)):
                streak_line = '{:.4f} {:.4f} 10 10 100 {:5.0f} 0 0.00'.format(float(newCenX[i]), float(newCenY[i]),
                                                                              newflux[i])
                STARS.write(streak_line + "\n")

            STARS.close()