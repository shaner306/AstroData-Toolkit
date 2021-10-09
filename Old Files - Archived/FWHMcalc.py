from astropy.io import fits
import astropy.units as u
from astropy.stats import sigma_clipped_stats, SigmaClip, gaussian_fwhm_to_sigma, gaussian_sigma_to_fwhm
from astropy.nddata import Cutout2D
from astropy.modeling.functional_models import Gaussian1D
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.table import Table
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle
import datetime
from photutils.detection import IRAFStarFinder
from photutils import CircularAperture, make_source_mask
from photutils.psf import DAOGroup, BasicPSFPhotometry, IntegratedGaussianPRF, IterativelySubtractedPSFPhotometry
from photutils.background import MMMBackground
import numpy as np
import os
from math import sqrt, atan
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from time import process_time

class UnitError(Exception):
    pass

# The directory where the HIP 88427 files are stored.
# I used HIP 88427 because it has already been plate solved with Pinpoint, so I can use that to verify my results.
directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2020 12 28\2020 12 28 - HIP 88427 - G Band - PinPoint'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2020 12 28\2020 12 28 - Intelsat 10-02 - G Band'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-02-07 - Calibrated\Intelsat 10-02\LIGHT\G'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2020 12 28\2020 12 28 - HIP 88427 - G Band - Unsolved'
# Controls if the plots will be shown. Useful to be able to control when debugging.
plot_results = False
size = 25                                                                                                               # Size of the cutout to plot and fit the gaussian to (pixels).
hsize = int((size - 1) / 2)                                                                                             # Half of the size of the cutout.
gaussian1D = Gaussian1D()                                                                                               # Initialize the 1D Gaussian that will be used to fit the data.
fitter = LevMarLSQFitter()                                                                                              # Initialize the fitter that will be used to fit the 1D Gaussian.
df = pd.DataFrame()
brightest_obj = []

for entry in os.scandir(directory):                                                                                     # Loop over every fits file in the directory.
    if entry.path.endswith(".fits"):                                                                                    # Only open the file if it is a .fits file.
        with fits.open(entry.path) as image:                                                                            # Open the image in the 'with' statement. Best practice.
            # Store the header and image data as variables that can be easily accessed.                                 #####
            hdr = image[0].header                                                                                       # Store the fits header as a variable.
            wcs = WCS(hdr)                                                                                              # Read World Coordinate data from the image header.
            # if wcs.wcs.ctype[0] == '':
            #     wcs.wcs.ctype[0] = 'RA---TAN'
            #     wcs.wcs.ctype[1] = 'DEC--TAN'
            #     if hdr.comments['RA'][1:4] == "deg" or hdr.comments['RA'][1:4] == "dms":
            #         RA = hdr['RA']  # Store the right ascensions in an array with the same size and position as dates
            #     elif hdr.comments['RA'][1:4] == "hms":
            #         # Convert to degrees.
            #         # print("Converting RA from hms to degrees.")
            #         RA_hms = Angle(hdr['RA'], unit=u.hourangle)
            #         RA = RA_hms.to_string(unit=u.deg)
            #         # print(RA_hms)
            #         # print(RA[i])
            #     else:
            #         print("Unknown RA units.")
            #         raise UnitError
            #     if hdr.comments['DEC'][1:4] == "deg" or hdr.comments['DEC'][1:4] == "dms":
            #         dec = hdr['DEC']  # Store the declinations in an array with the same size and position as dates
            #     elif hdr.comments['DEC'][1:4] == "hms":
            #         # Convert to degrees
            #         # print("Converting DEC from hms to degrees.")
            #         dec_hms = Angle(hdr['DEC'], unit=u.hourangle)
            #         dec = dec_hms.to_string(unit=u.deg)
            #     else:
            #         print("Unknown RA units.")
            #         raise UnitError
            #     wcs.wcs.crval[0] = RA
            #     wcs.wcs.crval[1] = dec
            #     wcs.wcs.crpix[0] = hdr['NAXIS1'] / 2
            #     wcs.wcs.crpix[1] = hdr['NAXIS2'] / 2
            #     # wcs.wcs.cd = [[0.0002629, 0.00034456],
            #     #               [-0.00034422, 0.00026317]]
            #     # wcs.wcs.cd = [[0.00003, 0.00003],
            #     #               [-0.00003, 0.00003]]
            # print(wcs)
            imgdata = image[0].data                                                                                     # Store the image as a variable.
            exptime = hdr['EXPTIME'] * u.s                                                                              # Store the exposure time with unit seconds.
            # mask = make_source_mask(imgdata, nsigma=2, npixels=5, dilate_size=11, sigclip_iters=None)                 # I tried this and found it made no difference on the small dataset tested.
            mean_val, median_val, std_val = sigma_clipped_stats(imgdata)                                                # Calculate background stats.
            iraffind = IRAFStarFinder(threshold=median_val+3*std_val, fwhm=2)                                           # Find stars using IRAF.
            irafsources = iraffind(imgdata - median_val)                                                                # Subtract background median value.
            irafpositions = np.transpose((irafsources['xcentroid'], irafsources['ycentroid']))                          # Store source positions as a numpy array.
            skypositions = wcs.pixel_to_world(irafsources['xcentroid'],irafsources['ycentroid'])                        # Convert star positions from pixels to RA and dec.
            print(skypositions)
            # Convert the RA/dec of the stars into Az/El
            ### This will only work for images that have been plate solved with PinPoint a priori. I am working on a way
            # to fix that. #####
            date = np.datetime64(datetime.datetime.strptime(hdr['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f'))                   # Store the dates in an array that can be used later
            lat = hdr['SITELAT']
            long = hdr['SITELONG']
            elev = hdr['SITEELEV']
            current_location = EarthLocation.from_geodetic(lon=long, lat=lat, height=elev)
            current_aa = AltAz(location=current_location, obstime=date)
            altazpositions = skypositions.transform_to(current_aa)
            print(altazpositions)
            iraf_fwhms = irafsources['fwhm']                                                                            # Save FWHM in a list.
            iraf_fwhms = np.array(iraf_fwhms)                                                                           # Convert FWHM list to numpy array.
            # Print information about the file                                                                          #####
            print(f"\n{os.path.basename(entry.path)}")                                                                  # Print filename.
            # Calculate statistics of the FWHMs given by the loop over each source.                                     #####
            iraf_sdom = iraf_fwhms.std() / sqrt(len(iraf_fwhms))                                                        # Standard deviation of the mean. Not used.
            num_IRAF_sources = len(irafsources)                                                                         # Number of stars IRAF found.
            print(f"No. of IRAF sources: {num_IRAF_sources}")                                                           # Print number of IRAF stars found.
            iraf_fwhm = iraf_fwhms.mean()                                                                               # Calculate IRAF FWHM mean.
            iraf_std = iraf_fwhms.std()                                                                                 # Calculate IRAF standard deviation.
            print(f"IRAF Calculated FWHM (pixels): {iraf_fwhm:.3f} +/- {iraf_std:.3f}")                                 # Print IRAF FWHM.
            if 'FWHM' in hdr.keys():                                                                                    # Only if it has been plate solved.
                pinpointFWHM = hdr['FWHM']                                                                              # Store the FWHM from Pinpoint as a variable.
                print(f"From PinPoint:{hdr['HISTORY'][-2]}")                                                            # Print the number of stars that PinPoint was able to match.
                print(f"PinPoint FWHM (pixels): {pinpointFWHM:.3f}")                                                    # Print Pinpoint's FWHM result.
                iraf_error = 100 * (iraf_fwhm - pinpointFWHM) / pinpointFWHM                                            # Percent Error of IRAF FWHM.
                print(f"IRAF Error: {iraf_error:.1f} %")                                                                # Print IRAF percent error.
            # Calculate the flux of each star using PSF photometry. This uses the star positions calculated by
            # IRAFStarFinder earlier.
            daogroup = DAOGroup(2*iraf_fwhm)                                                                            # Groups overlapping stars together.
            # sigma_clip = SigmaClip()
            # mmm_bkg = MMMBackground(sigma_clip=sigma_clip)
            # bkg_value = mmm_bkg(imgdata)
            psf_model = IntegratedGaussianPRF(sigma=iraf_fwhm*gaussian_fwhm_to_sigma)                                   # Defime the PSF model to be used for photometry.
            psf_model.x_0.fixed = True                                                                                  # Don't change the initial 'guess' of the star x positions to be provided.
            psf_model.y_0.fixed = True                                                                                  # Don't change the initial 'guess' of the star y positions to be provided.
            # Provide the initial guesses for the x-y positions and the flux. The flux will be fit using the psf_model,
            # so that will change, but the star positions will remain the same.
            pos = Table(names=['x_0', 'y_0', 'flux_0'],
                        data=[irafsources['xcentroid'], irafsources['ycentroid'], irafsources['flux']])
            # Initialize the photometry to be performed. Do not estimate the background, as it will be subtracted from
            # the image when fitting the PSF.
            photometry = BasicPSFPhotometry(group_maker=daogroup,
                                            bkg_estimator=None,
                                            psf_model=psf_model,
                                            fitter=LevMarLSQFitter(),
                                            fitshape=size)
            # No Improvement over BasicPSFPhotometry
            # iterative_photometry = IterativelySubtractedPSFPhotometry(group_maker=daogroup,
            #                                                           bkg_estimator=mmm_bkg,
            #                                                           psf_model=psf_model,
            #                                                           fitter=LevMarLSQFitter(),
            #                                                           finder=iraffind,
            #                                                           niters=None,
            #                                                           fitshape=size)
            # Perform the photometry on the background subtracted image. Also pass the fixed x-y positions and the
            # initial guess for the flux.
            result_tab = photometry(image=imgdata - median_val, init_guesses=pos)
            # iterative_result_tab = iterative_photometry(image=imgdata)
            # iterative_result_tab.sort('flux_fit', reverse=True)
            fluxes = result_tab['flux_fit']                                                                             # Store the fluxes as a list.
            fluxes = np.array(fluxes) * u.ct                                                                            # Convert the fluxes to a numpy array and add the unit of count to it.
            fluxes = fluxes / exptime                                                                                   # Normalize the fluxes by exposure time (unit is now counts / second)
            instr_mags_units = u.Magnitude(fluxes)                                                                      # Convert the fluxes to an instrumental magnitude.
            instr_mags = instr_mags_units.value                                                                         # Store the magnitudes without the unit attached.
            # iterative_fluxes = iterative_result_tab['flux_fit']
            # iterative_fluxes = np.array(iterative_fluxes) * u.ct
            # iterative_instr_mag_units = u.Magnitude(iterative_fluxes)
            # iterative_instr_mag = iterative_instr_mag_units.value
            # print(instr_mags)                                                                                           # Print the instrumental magnitudes.
            # print(iterative_instr_mag)
            # PSF_sigmas = iterative_result_tab['sigma_fit']
            # PSF_sigmas = np.array(PSF_sigmas)
            # PSF_fwhms = PSF_sigmas * gaussian_sigma_to_fwhm
            # PSF_fwhm = PSF_fwhms.mean()
            # PSF_fwhm_std = PSF_fwhms.std()
            # print(f"PSF fwhm (pixels): {PSF_fwhm:.3f} +/- {PSF_fwhm_std:.3f}")
            # Fit a 1D gaussian in both x and y (note that variables 'y_vals' and 'z_vals' correspond to the image data #####
            # for the x and y positions, respectively). This was used to verify the results obtained by IRAFStarFinder. #####
            stddev = iraf_fwhm * gaussian_fwhm_to_sigma                                                                 # Initial guess used for the Gaussian1D stddev.
            gaussianFWHMs = []                                                                                          # Initialize array to store FWHMs from the gaussian fits. Difficult to preallocate.
            for i in range(0, num_IRAF_sources):                                                                        # Loop over every star detected by IRAFStarFinder.
                current_pos = irafpositions[i, :]                                                                       # Get the current star's position from the irafpositions table.
                x_pos = int(current_pos[0])                                                                             # Store the star's x coordinate as its own variable.
                y_pos = int(current_pos[1])                                                                             # Store the star's y coordinate as its own variable.
                if (x_pos > hsize and x_pos < (imgdata.shape[1] - 1 - hsize)) and (
                        y_pos > hsize) and (y_pos < (imgdata.shape[0] - 1 - hsize)):                                    # Only consider stars that are at least hsize away from the image border.
                    gaussian1D = Gaussian1D(amplitude=irafsources['flux'][i], mean=x_pos, stddev=stddev)                # Initialize the gaussian1D function to be fit. Initial values given to hopefully improve the quality of the fit.
                    x_vals = np.arange(start=x_pos - hsize - 1, stop=x_pos + hsize, step=1)                             # Initialze the array to be used as 'x' in the gaussian fit.
                    y_cutout = Cutout2D(imgdata - median_val, (x_pos, y_pos), (1, size))                                # Cutout the bacground subtracted image data in a 1 x size array centred on the star's position.
                    y_vals = np.array(y_cutout.data)                                                                    # Convert the cutout's data to a numpy array.
                    y_vals = np.squeeze(y_vals)                                                                         # Change the array's shape to be the same as 'x_vals' to be used as 'y' in the gaussian fit.
                    fitted_model = fitter(gaussian1D, x_vals, y_vals)                                                   # Fit the gaussian in the x direction.
                    gaussianFWHMs.append(fitted_model.stddev * gaussian_sigma_to_fwhm)                                  # Store the stddev of the gaussian as a FWHM.
                    z_cutout = Cutout2D(imgdata - median_val, (x_pos, y_pos), (size, 1))                                # Cutout the bacground subtracted image data in a size x 1 array centred on the star's position.
                    z_vals = np.array(z_cutout.data)                                                                    # Convert the cutout's data to a numpy array.
                    z_vals = np.squeeze(z_vals)                                                                         # Change the array's shape to be the same as 'x_vals' to be used as 'y' in the gaussian fit.
                    fitted_model = fitter(gaussian1D, x_vals, y_vals)                                                   # Fit the gaussian in the y direction.
                    gaussianFWHMs.append(fitted_model.stddev * gaussian_sigma_to_fwhm)                                  # Store the stddev of the gaussian as a FWHM.
            gaussianFWHMs = np.array(gaussianFWHMs)                                                                     # Convert the FWHM list to a numpy array.
            gaussianFWHM = gaussianFWHMs.mean()                                                                         # Calculate the mean of the FWHMs from the gaussian fits.
            gaussianFWHM_std = gaussianFWHMs.std()                                                                      # Calculate the standard deviation of the FWHMs from the gaussian fits.
            print(f"Gaussian Calculated FWHM: {gaussianFWHM:.3f} +/- {gaussianFWHM_std:.3f}")                           # Print the FWHM calculated by the gaussian fits.
            percent_error_gaussian_v_iraf = 100 * (gaussianFWHM - iraf_fwhm) / iraf_fwhm                                # Percent difference between the gaussian fit and IRAFStarFinder FWHMs.
            print(f"Percent Error (gaussian v. IRAF): {percent_error_gaussian_v_iraf:.1f} %")
            if 'FWHM' in hdr.keys():# Print the gaussian fit vs. IRAFStarFinder FWHMs percent difference.
                percent_error_gaussian_v_pinpoint = 100 * (gaussianFWHM - pinpointFWHM) / pinpointFWHM                  # Percent difference between the gaussian fit and PinPoint FWHMs.
                print(f"Percent Error (gaussian v. PinPoint): {percent_error_gaussian_v_pinpoint:.1f} %")               # Print the gaussian fit vs. PinPoint FWHMs percent difference.

            # Calculate the FWHM in units of arcseconds as opposed to pixels. This allows the user to more easily
            # interpret the differences because it is more universal to interpret FWHM as arcsec instead of pixels,
            # as the number of arcsec / pixel varies depending on the CCD and telescope. This can also enable the
            # capability to image the same site with different hardware and still be able to compare the seeing as a
            # function of AltAz without needing to modify anything.

            focal_length = hdr['FOCALLEN'] * u.mm                                                                       # Store the telescope's focal length with unit millimetres.
            xpixsz = hdr['XPIXSZ']                                                                                      # Store the size of the x pixels.
            ypixsz = hdr['XPIXSZ']                                                                                      # Store the size of the y pixels.
            if xpixsz == ypixsz:                                                                                        # If the pixels are square.
                pixsz = xpixsz * u.um                                                                                   # Store the pixel size with unit micrometre.
                # Can find FOV by finding deg/pix and then multiplying by the x and y number of pix (NAXIS).
                rad_per_pix = atan(pixsz / focal_length) * u.rad                                                        # Calculate the angular resolution of each pixel. Store with unit radians.
                arcsec_per_pix = rad_per_pix.to(u.arcsec)                                                               # Convert the per pixel angular resultion to arcseconds.
                iraf_FWHM_arcsec = iraf_fwhm * arcsec_per_pix.value                                                     # Convert the IRAFStarFinder FWHM from pixels to arcsec.
                iraf_std_arcsec = iraf_std * arcsec_per_pix                                                             # Convert the IRAFStarFinder FWHM standard deviation from pixels to arcsec.
                gaussianFWHM_arcsec = gaussianFWHM * arcsec_per_pix.value                                               # Convert the gaussian fit FWHM from pixels to arcsec.
                gaussianFWHM_std_arcsec = gaussianFWHM_std * arcsec_per_pix                                             # Convert the gaussian fit FWHM standard deviation from pixels to arcsec.
                if 'FWHM' in hdr.keys():
                    pinpointFWHM_arcsec = pinpointFWHM * arcsec_per_pix                                                 # Convert the PinPoint FWHM from pixels to arcsec.
                print("Converting FWHM from pixels to arcsec...")                                                       # Print this line to more easily differentiate pixel values from arcsec.
                print(f"IRAF Calculated FWHM (arcsec): {iraf_FWHM_arcsec:.3f} +/- {iraf_std_arcsec:.3f}")               # Print the IRAFStarFinder FWHM in arcsec.
                print(f"FWHM from Gaussian fit (arcsec): {gaussianFWHM_arcsec:.3f} +/- {gaussianFWHM_std_arcsec:.3f}")  # Print the gaussian fit FWHM in arcsec.
                if 'FWHM' in hdr.keys():                                                                                # Only if the image has been plate solved by PinPoint.
                    print(f"PinPoint FWHM (arcsec): {pinpointFWHM_arcsec:.3f}")                                         # Print the PinPoint FWHM in arcsec.

            # Plot the 5 brightest and 5 dimmest stars in each image. This is to get a better understanding of what a
            # good estimate of the FWHM should be (in pixels).

            fig, ax = plt.subplots(nrows=2, ncols=5, squeeze=True)                                                      # Initialize a 2 x 5 matplotlib subplot environment.
            rows = [f"5 {row} Stars" for row in ["Brightest", "Dimmest"]]                                               # The array that will be printed beside each row.
            pad = 5                                                                                                     # Padding to add beside the row labels.
            # Print the labels for the rows.
            for axis, row in zip(ax[:, 0], rows):
                axis.annotate(row, xy=(0, 0.5), xytext=(-axis.yaxis.labelpad - pad, 0),
                              xycoords=axis.yaxis.label, textcoords='offset points',
                              size='large', ha='right', va='center', rotation=90)
            fig.tight_layout(w_pad=-1)                                                                                  # Size of the subplots. Makes sure they aren't overlapping.
            fig.subplots_adjust(top=0.9)                                                                                # Brings the subplots down so that they aren't overlapping the title.
            j = 0                                                                                                       # Initialize variable for the second loop.
            for i in range(0, 5):                                                                                       # Loop over the 5 brightest stars.
                # The below repeats the x direction gaussian fit from earlier. The differences are that it plots the
                # image/fit and that it only calculates it for 5 images.
                current_pos = irafpositions[i, :]
                x_pos = int(current_pos[0])
                y_pos = int(current_pos[1])
                gaussian1D = Gaussian1D(amplitude=irafsources['flux'][i], mean=x_pos, stddev=stddev)
                x_vals = np.arange(start=x_pos - hsize - 1, stop=x_pos + hsize, step=1)
                x_plots = np.arange(start=x_pos - hsize - 1, stop=x_pos + hsize, step=0.01)
                y_cutout = Cutout2D(imgdata - median_val, (x_pos, y_pos), (1, size))
                y_vals = np.array(y_cutout.data)
                y_vals = np.squeeze(y_vals)
                fitted_model = fitter(gaussian1D, x_vals, y_vals)
                # Below plots the image data as well as the fit calculated above.
                ax[0, i].plot(x_vals, y_vals + median_val, 'o', label='Data')                                           # Plot the image data. The background has been added back.
                ax[0, i].plot(x_plots, fitted_model(x_plots) + median_val, label='Gaussian Fit')                        # Plot the gaussian fit. The background has been added to it.
                # ax[0, i].set_ylabel("Flux")
                ax[0, i].set_xlabel("X Location (pixels)")                                                              # Display the x axis title.
                ax[0, i].xaxis.set_major_locator(plt.MaxNLocator(6))                                                    # Show a maximum of 6 tick marks along the x axis.
            for i in range(-6, -1):                                                                                     # Loop over the 5 dimmest stars in the image.
                # This loop is the exact same as the one above, except that it is over the 5 dimmest stars, rather than
                # the 5 brightest images.
                current_pos = irafpositions[i, :]
                x_pos = int(current_pos[0])
                y_pos = int(current_pos[1])
                gaussian1D = Gaussian1D(amplitude=irafsources['flux'][i], mean=x_pos, stddev=stddev)
                x_vals = np.arange(start=x_pos - hsize - 1, stop=x_pos + hsize, step=1)
                x_plots = np.arange(start=x_pos - hsize - 1, stop=x_pos + hsize, step=0.01)
                y_cutout = Cutout2D(imgdata - median_val, (x_pos, y_pos), (1, size))
                y_vals = np.array(y_cutout.data)
                y_vals = np.squeeze(y_vals)
                fitted_model = fitter(gaussian1D, x_vals, y_vals)
                ax[1, j].plot(x_vals, y_vals + median_val, 'o', label='Data')
                ax[1, j].plot(x_plots, fitted_model(x_plots) + median_val, label='Gaussian Fit')
                # ax[1, j].set_ylabel("Flux")
                ax[1, j].set_xlabel("X Location (pixels)")
                ax[1, j].xaxis.set_major_locator(plt.MaxNLocator(6))
                j += 1                                                                                                  # Increase the iterator.
            # Set the plot title. Include the number of stars, filename, PinPoint FWHM, IRAFStarFinder FWHM, and the
            # gaussian fit FWHM in the title.
            if 'FWHM' in hdr.keys():
                plt.suptitle(
                    f"Plot of Gaussian Fit of Brightest and Dimmest Stars (of {num_IRAF_sources} total stars) in "
                    f"{os.path.basename(entry.path)}\n"
                    f"PinPoint FWHM = {pinpointFWHM:.3f} pixels ({pinpointFWHM_arcsec:.3f})\n"
                    f"IRAFStarFinder FWHM = {iraf_fwhm:.3f} +/- {iraf_std:.3f} pixels "
                    f"({iraf_FWHM_arcsec:.3f} +/- {iraf_std_arcsec:.3f})\n"
                    f"FWHM from gaussian fit = {gaussianFWHM:.3f} +/- {gaussianFWHM_std:.3f} pixels "
                    f"({gaussianFWHM_arcsec:.3f} +/- {gaussianFWHM_std_arcsec:.3f})",
                    size='medium')
            else:
                plt.suptitle(
                    f"Plot of Gaussian Fit of Brightest and Dimmest Stars (of {num_IRAF_sources} total stars) in "
                    f"{os.path.basename(entry.path)}\n"
                    f"IRAFStarFinder FWHM = {iraf_fwhm:.3f} +/- {iraf_std:.3f} pixels "
                    f"({iraf_FWHM_arcsec:.3f} +/- {iraf_std_arcsec:.3f})\n"
                    f"FWHM from gaussian fit = {gaussianFWHM:.3f} +/- {gaussianFWHM_std:.3f} pixels "
                    f"({gaussianFWHM_arcsec:.3f} +/- {gaussianFWHM_std_arcsec:.3f})",
                    size='medium')
            handles, labels = ax[0,0].get_legend_handles_labels()                                                       # Get the legend information from the last subplot.
            fig.legend(handles, labels, loc='upper right')                                                              # Display the legend in the upper right corner.
            if plot_results:                                                                                            # Only if plot_results is True.
                plt.show()                                                                                              # Show the plot.
            plt.close()                                                                                                 # Close the plot. Clears it from RAM.

            # Plot the visualization of the stars that were detected by the algorithm used in this file and (if
            # applicable) the stars that PinPoint found.
            if plot_results and os.path.exists(f"{entry.path}.stars"):                                                  # If it has been plate solved. And plot_results is True.
                f = open(f"{entry.path}.stars", "r")                                                                    # Open PinPoint star locations.f
                pinpoint_x = []                                                                                         # Initialize x locations.
                pinpoint_y = []                                                                                         # Initialize y locations.
                for line in f:                                                                                          # Iterate over .stars file.
                    fields = line.split(" ")                                                                            # File is delimited by a space.
                    pinpoint_x.append(float(fields[0]))                                                                 # First value is the x position.
                    pinpoint_y.append(float(fields[1]))                                                                 # Second value is the y position.
                pinpoint_positions = np.transpose((pinpoint_x, pinpoint_y))                                             # Save x and y in a way that can be plotted.
                # Convert positions to circles with radius 'r' for plotting.
                irafapertures = CircularAperture(irafpositions, r=10.)
                pinpointapertures = CircularAperture(pinpoint_positions, r=10.)
                norm = ImageNormalize(stretch=SqrtStretch())                                                            # Normalize image using a Square Root Stretch.
                plt.imshow(imgdata, cmap='gray', origin='lower', norm=norm, interpolation='nearest')                    # Plot image.
                irafapertures.plot(color='b', lw=1.5, alpha=0.5)                                                        # Plot IRAF stars.
                pinpointapertures.plot(color='g', lw=1.5, alpha=0.5)                                                    # Plot PinPoint stars.
                # Define the colours to be used for the legend. Unfortunately this cannot be automated as the legend
                # would include an entry for each individual star. This does not have to be automated, however, because
                # the plotting is only used to verify results.
                legend_elements = [
                    Line2D([0], [0], color='w', marker='o', markerfacecolor='b', markersize=7, label='IRAF'),
                    Line2D([0], [0], color='w', marker='o', markerfacecolor='g', markersize=7, label='PinPoint'),
                    Line2D([0], [0], color='w', marker='o', markerfacecolor='c', markersize=7, label='PinPoint + IRAF'),
                ]
                plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')                         # Include legend in the plot.
                # plt.colorbar()
                plt.title(label=os.path.basename(entry.path))                                                           # Title is the filename.
                plt.show()                                                                                              # Show the image and overlays.
                plt.close()                                                                                             # Close plot. Clears it from RAM.
                f.close()                                                                                               # Close .stars file. Best practice.
            elif plot_results and not os.path.exists(f"{entry.path}.stars"):                                            # If it has not been plate solved. And plot_results is True.
                irafapertures = CircularAperture(irafpositions, r=10.)                                                  # Plot IRAF stars.
                norm = ImageNormalize(stretch=SqrtStretch())                                                            # Normalize image using a Square Root Stretch.
                plt.imshow(imgdata, cmap='gray', origin='lower', norm=norm, interpolation='nearest')                    # Plot image.
                irafapertures.plot(color='b', lw=1.5, alpha=0.5)                                                        # Plot IRAF stars.
                plt.title(label=os.path.basename(entry.path))                                                           # Title is the filename.
                plt.show()                                                                                              # Show the image and overlays.
                plt.close()                                                                                             # Close plot. Clears it from RAM.
            # print(irafsources['peak'] + median_val)                                                                     # Akin to 'max_pixel' from Shane's spreadsheet.
            print(result_tab['x_0', 'y_0', 'flux_fit', 'flux_unc'])                                                     # Print the fluxes and their uncertainty for the current image.
            repeat_bkg = np.full(iraf_fwhms.shape, median_val)
            repeat_std = np.full(iraf_fwhms.shape, std_val)
            repeat_exptime = np.full(iraf_fwhms.shape, exptime.value)
            repeat_fwhm = np.full(iraf_fwhms.shape, iraf_fwhm)
            repeat_fwhm_std = np.full(iraf_fwhms.shape, iraf_std)
            repeat_filename = np.full(iraf_fwhms.shape, os.path.basename(entry.path), dtype=object)
            df2 = pd.DataFrame({'X': irafsources['xcentroid'],
                       'Y': irafsources['ycentroid'],
                       'flux': result_tab['flux_fit'],
                       'max_pixel': irafsources['peak'] + median_val,
                       'median_bkg': repeat_bkg,
                       'bkg_sigma': repeat_std,
                       'Exposure': repeat_exptime,
                       'FWHM': repeat_fwhm,
                       'FWHM sigma': repeat_fwhm_std,
                       'Instrumental Magnitude': instr_mags,
                       'FileName': repeat_filename})
            df = df.append(df2)
            # print(df)
            # print(iterative_result_tab['x_fit', 'y_fit', 'flux_fit'])
            # print(iterative_result_tab['flux_fit'])
            # print(f"Sig_clip background: {median_val}")
            # print(f"MMMBackground: {bkg_value}")
            # num_iters = 1000
            # start = process_time()
            # for i in range(0, num_iters):
            #     sigma_clipped_stats(imgdata)
            # end = process_time()
            # avg = (end - start) / num_iters
            # print(f"Sig_clip_time: {avg}")
            # start = process_time()
            # for i in range(0, num_iters):
            #     mmm_bkg(imgdata)
            # end = process_time()
            # avg = (end - start) / num_iters
            # print(f"MMMBackground time: {avg}")
            print("")                                                                                                   # Print a new line between each file iteration.
# df.to_csv('detected_stars.csv', index=False)