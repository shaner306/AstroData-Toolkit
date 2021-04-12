from astropy import table
from astropy.coordinates import EarthLocation, AltAz, Angle, SkyCoord, match_coordinates_sky
from astropy.io import fits, ascii
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.table import Table, QTable
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
import datetime
from math import floor, atan
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import os
from photutils import CircularAperture
from photutils.detection import IRAFStarFinder
from photutils.psf import DAOGroup, BasicPSFPhotometry, IntegratedGaussianPRF
import win32com.client as win32

# The directory that holds all of the plate solved images for the night.
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2020 12 28\Solved Images'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-03-10 - Calibrated\Solved Images'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-03-20 - Calibrated\Solved Images'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-03-21 - Calibrated\Solved Images'
directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Landolt Stars'

ground_based = False

catloc = r'C:\Program Files (x86)\PinPoint\UCAC4'

f = win32.Dispatch("Pinpoint.plate")
f.DetachFITS

size = 25                                                                                                               # Size of the cutout to plot and fit the gaussian to (pixels).
hsize = int((size - 1) / 2)                                                                                             # Half of the size of the cutout.
fitter = LevMarLSQFitter()                                                                                              # Initialize the fitter that will be used to fit the PSF.
# ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\FITS Tutorial\Reference_stars.csv'                             # Location of the file containing the reference stars.
# ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Landolt Stars\Landolt_reference.csv'
ref_stars_file = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Landolt Stars\2009_Landolt_Standard_Stars.txt'

# Initialize all of the arrays that will be added to the large AstroPy table containing all of the information on the
# reference stars that were detected in each image.
ref_star_name = []
times = []
flux_table = []
exposure = []
ref_star_RA = []
ref_star_dec = []
img_star_RA = []
img_star_dec = []
angular_separation = []
ref_star_x = []
ref_star_y = []
img_star_x = []
img_star_y = []
ref_star_mag = []
img_star_mag = []
filters = []
B_V_apparents = []
U_B_apparents = []
V_R_apparents = []
V_I_apparents = []
V_sigma_apparents = []
img_star_secz = []
X_rounded = []

try:
    reference_stars = ascii.read(ref_stars_file, format='basic', delimiter='\t', guess=False, encoding='UTF-8')         # Read the reference stars file to an AstroPy Table as a csv with a Tab delimiter using UTF-8 encoding.
except Exception as e:                                                                                                  # Usually if the file is not in the expected format. If that's not why the exception was thrown, then the next line will cause one as well.
    reference_stars = ascii.read(ref_stars_file, encoding='UTF-8')                                                      # Read the reference stars file in any way that ascii.read can, but still with UTF-8 encoding.
# reference_stars = ascii.read(ref_stars_file, encoding='UTF-8')
reference_stars = reference_stars.filled(np.nan)
print(reference_stars)                                                                                                  # Print the AstroPy table containing all of the information from the reference stars file.
ref_star_positions = SkyCoord(ra=reference_stars['RA'], dec=reference_stars['Dec'], unit=(u.hourangle, u.deg))          # Convert the stars' locations to a SkyCoord object. Assumes RA has units hms and dec has units dms or deg.
print(ref_star_positions)                                                                                               # Print the SkyCoord object of the reference stars' locations.
# The for loops below loop over all files in all subfolders of 'directory' and then opens each fits file one by one.
for dirpath, dirnames, filenames in os.walk(directory):
    for filename in filenames:
        # if filename.endswith(".fits"):
        if filename.endswith("_clean.fits"):
            with fits.open(os.path.join(dirpath, filename)) as image:                                                   # Open the fits file. Closes once the with statement ends to avoid keeping it in RAM.
                hdr = image[0].header                                                                                   # Store the fits header as a variable.
                imgdata = image[0].data                                                                                 # Store the image as a variable.
            print(filename)
            # if 'PLTSOLVD' in hdr.keys():
            #     if hdr['PLTSOLVD'] == True:
            #         print("Plate had already been solved by PinPoint.")
            #     else:
            #         f.AttachFITS(os.path.join(dirpath, filename))
            #         f.Declination = f.targetDeclination
            #         f.RightAscension = f.targetRightAscension
            #         # yBin = 4.33562092816E-004 * 3600
            #         # xBin = 4.33131246330E-004 * 3600
            #         focal_length = hdr['FOCALLEN'] * u.mm                                                               # Store the telescope's focal length with unit millimetres.
            #         xpixsz = hdr['XPIXSZ'] * u.um                                                                       # Store the size of the x pixels.
            #         ypixsz = hdr['YPIXSZ'] * u.um                                                                       # Store the size of the y pixels.
            #         rad_per_xpix = atan(xpixsz / focal_length) * u.rad                                                  # Calculate the angular resolution of each pixel. Store with unit radians.
            #         arcsec_per_xpix = rad_per_xpix.to(u.arcsec)                                                         # Convert the per pixel angular resultion to arcseconds.
            #         rad_per_ypix = atan(ypixsz / focal_length) * u.rad                                                  # Calculate the angular resolution of each pixel. Store with unit radians.
            #         arcsec_per_ypix = rad_per_ypix.to(u.arcsec)                                                         # Convert the per pixel angular resultion to arcseconds.
            #         xBin = arcsec_per_xpix.value
            #         yBin = arcsec_per_ypix.value
            #         f.ArcsecperPixelHoriz = xBin
            #         f.ArcsecperPixelVert = yBin
            #         f.Catalog = 11
            #         f.CatalogPath = catloc
            #         f.CatalogMaximumMagnitude = 13
            #         f.CatalogExpansion = 0.8
            #         f.SigmaAboveMean = 3.0
            #         f.FindImageStars
            #         f.FindCatalogStars
            #         f.MaxSolveTime = 60
            #         f.MaxMatchResidual = 1.5
            #         flag = 0
            #         f.FindCatalogStars()
            #         f.Solve()
            #         f.MatchedStars.count
            #         f.FindImageStars()
            #         f.updateFITS()
            #         f.DetachFITS()
            #         print("Plate was solved by PinPoint.")
            #         with fits.open(os.path.join(dirpath, filename)) as image:                                           # Open the fits file. Closes once the with statement ends to avoid keeping it in RAM.
            #             hdr = image[0].header                                                                           # Store the fits header as a variable.
            #             imgdata = image[0].data                                                                         # Store the image as a variable.
            # else:
            #     f.AttachFITS(os.path.join(dirpath, filename))
            #     f.Declination = f.targetDeclination
            #     f.RightAscension = f.targetRightAscension
            #     focal_length = hdr['FOCALLEN'] * u.mm                                                                   # Store the telescope's focal length with unit millimetres.
            #     xpixsz = hdr['XPIXSZ'] * u.um                                                                           # Store the size of the x pixels.
            #     ypixsz = hdr['YPIXSZ'] * u.um                                                                           # Store the size of the y pixels.
            #     rad_per_xpix = atan(xpixsz / focal_length) * u.rad                                                      # Calculate the angular resolution of each pixel. Store with unit radians.
            #     arcsec_per_xpix = rad_per_xpix.to(u.arcsec)                                                             # Convert the per pixel angular resultion to arcseconds.
            #     rad_per_ypix = atan(ypixsz / focal_length) * u.rad                                                      # Calculate the angular resolution of each pixel. Store with unit radians.
            #     arcsec_per_ypix = rad_per_ypix.to(u.arcsec)                                                             # Convert the per pixel angular resultion to arcseconds.
            #     xBin = arcsec_per_xpix.value
            #     yBin = arcsec_per_ypix.value
            #     f.ArcsecperPixelHoriz = xBin
            #     f.ArcsecperPixelVert = yBin
            #     f.Catalog = 11
            #     f.CatalogPath = catloc
            #     f.CatalogMaximumMagnitude = 13
            #     f.CatalogExpansion = 0.8
            #     f.SigmaAboveMean = 3.0
            #     f.FindImageStars
            #     f.FindCatalogStars
            #     f.MaxSolveTime = 60
            #     f.MaxMatchResidual = 1.5
            #     flag = 0
            #     f.FindCatalogStars()
            #     f.Solve()
            #     f.MatchedStars.count
            #     f.FindImageStars()
            #     f.updateFITS()
            #     f.DetachFITS()
            #     print("Plate was solved by PinPoint.")
            #     with fits.open(os.path.join(dirpath, filename)) as image:                                               # Open the fits file. Closes once the with statement ends to avoid keeping it in RAM.
            #         hdr = image[0].header                                                                               # Store the fits header as a variable.
            #         imgdata = image[0].data                                                                             # Store the image as a variable.
            wcs = WCS(hdr)                                                                                              # Read World Coordinate data from the image header.
            if 'EXPTIME' in hdr.keys():
                exptime = hdr['EXPTIME'] * u.s                                                                          # Store the exposure time with unit seconds.
            elif 'AEXPTIME' in hdr.keys():
                exptime = hdr['AEXPTIME'] * u.s
            mean_val, median_val, std_val = sigma_clipped_stats(imgdata, maxiters=None)                                 # Calculate background stats.
            iraffind = IRAFStarFinder(threshold=median_val+3*std_val, fwhm=2)                                           # Find stars using IRAF.
            irafsources = iraffind(imgdata - median_val)                                                                # Subtract background median value.
            # irafsources = iraffind(imgdata)
            # irafsources.sort('flux', reverse=True)                                                                      # Sort the stars by flux, from greatest to least.
            try:
                irafpositions = np.transpose((irafsources['xcentroid'], irafsources['ycentroid']))                          # Store source positions as a numpy array.
            except TypeError as e:
                print(e)
                print("Moving onto the next file.")
                continue
            skypositions = wcs.pixel_to_world(irafsources['xcentroid'], irafsources['ycentroid'])                       # Convert star positions from pixels to RA and dec.
            # date = np.datetime64(datetime.datetime.strptime(hdr['DATE-OBS'], '%Y-%m-%dT%H:%M:%S.%f'))                   # Store the dates in an array that can be used later
            date = Time(hdr['DATE-OBS'], format='fits')
            if ground_based:
                lat = hdr['SITELAT']                                                                                        # Store the latitude that was given in the FITS header.
                long = hdr['SITELONG']                                                                                      # Store the longitude that was given in the FITS header.
                elev = hdr['SITEELEV']                                                                                      # Store the elevation that was given in the FITS header.
                current_location = EarthLocation.from_geodetic(lon=long, lat=lat, height=elev)                              # Create an EarthLocation object of the observation site read from the FITS header.
                current_aa = AltAz(location=current_location, obstime=date)                                                 # Create an Altitude/Azimuth reference frame with the observation location and time.
                altazpositions = skypositions.transform_to(current_aa)                                                      # Convert the stars' RA/dec to Altitude/Azimuth (AltAz or Az/El).
            iraf_fwhms = irafsources['fwhm']                                                                            # Save FWHM in a list.
            iraf_fwhms = np.array(iraf_fwhms)                                                                           # Convert FWHM list to numpy array.
            num_IRAF_sources = len(irafsources)                                                                         # Number of stars IRAF found.
            iraf_fwhm = iraf_fwhms.mean()                                                                               # Calculate IRAF FWHM mean.
            # print(iraf_fwhm)
            iraf_std = iraf_fwhms.std()                                                                                 # Calculate IRAF standard deviation.
            daogroup = DAOGroup(2 * iraf_fwhm)                                                                          # Groups overlapping stars together.
            psf_model = IntegratedGaussianPRF(sigma=iraf_fwhm * gaussian_fwhm_to_sigma)                                 # Defime the PSF model to be used for photometry.
            psf_model.x_0.fixed = True                                                                                  # Don't change the initial 'guess' of the star x positions to be provided.
            psf_model.y_0.fixed = True                                                                                  # Don't change the initial 'guess' of the star y positions to be provided.
            # Provide the initial guesses for the x-y positions and the flux. The flux will be fit using the psf_model,
            # so that will change, but the star positions will remain the same.
            pos = Table(names=['x_0', 'y_0', 'flux_0'],
                        data=[irafsources['xcentroid'], irafsources['ycentroid'], irafsources['flux']])
            # Initialize the type of photometry to be performed. Do not estimate the background as it will be subtracted
            # from the image when fitting the PSF.
            photometry = BasicPSFPhotometry(group_maker=daogroup,
                                            bkg_estimator=None,
                                            psf_model=psf_model,
                                            fitter=LevMarLSQFitter(),
                                            fitshape=size)
            result_tab = photometry(image=imgdata - median_val, init_guesses=pos)                                       # Perform the PSF photometry on the background subtracted image using stars found by IRAFStarFinder earlier.
            fluxes = result_tab['flux_fit']                                                                             # Store the fluxes as a list.
            flux_uncs = result_tab['flux_unc']
            flux_uncs = np.array(flux_uncs) * u.ct
            flux_uncs = flux_uncs / exptime
            fluxes = np.array(fluxes) * u.ct                                                                            # Convert the fluxes to a numpy array and add the unit of count to it.
            fluxes = fluxes / exptime                                                                                   # Normalize the fluxes by exposure time (unit is now counts / second)
            instr_mags_units = u.Magnitude(fluxes)                                                                      # Convert the fluxes to an instrumental magnitude.
            instr_mags = instr_mags_units.value                                                                         # Store the magnitudes without the unit attached.
            snr = (fluxes / flux_uncs).value
            # print(snr)
            instr_mags_sigma = 1.0857 / np.sqrt(snr)
            # print(instr_mags_sigma)
            # print(skypositions)                                                                                         # Print the RA/dec of the stars detected by IRAFStarFinder.
            idx, sep2d, dist3d = match_coordinates_sky(ref_star_positions, skypositions)                                # Match the reference stars with the stars detected in the image.
            # print(idx)                                                                                                  # Print the index of the closest image star to each star in the reference stars file.
            # print(sep2d)                                                                                                # Print the 2D angular separation between the closest image star from each star in the reference stars file.
            # print(dist3d)
            # max_ref_sep = 3 * u.deg                                                                                     # The threshold for deciding if any of the reference stars are close to the image stars.
            max_ref_sep = 30 * u.arcsec
            possible_ref_star_index = np.where(sep2d < max_ref_sep)                                                     # Finds the indicies where the separation is less than max_ref_sep.
            if len(possible_ref_star_index[0]) == 0:                                                                    # If no stars are within max_ref_sep of the reference stars.
                print("No reference star detected in the image.")
                continue
                # raise Exception(f"No star in the image was within {max_ref_sep} of any reference star.")                # This will stop the program with an error. Haven't decided if this or 'continue' is the best approach.
            elif len(possible_ref_star_index[0]) == 1:                                                                  # If there is only 1 reference star within max_ref_sep of the image stars.
                possible_ref_star_index = int(possible_ref_star_index[0])                                               # Convert the array of length 1 to an int.
            else:                                                                                                       # If there are more than 1 reference stars possibly in the image.
                possible_ref_star_index = possible_ref_star_index[0]                                                    # Changes the variable from a tuple to an array.
            # print(possible_ref_star_index)                                                                              # Print the index or indicies of the stars from the reference stars file.
            possible_image_star_index = idx[possible_ref_star_index]                                                    # Store the index or indicies of which star(s) in the image correspond to the possible reference stars.
            # print(possible_image_star_index)                                                                            # Print the index or indicies of the image star(s) that are likely the same as the reference star(s).
            possible_ref_star = reference_stars[possible_ref_star_index]                                                # Store the row(s) of the table containing the reference star(s) that are likely found in the image.
            print(possible_ref_star)                                                                                    # Print the row(s) of the table containing the reference star(s) that are likely found in the image.
            possible_ref_star_loc = ref_star_positions[possible_ref_star_index]                                         # Store the RA/dec of the possible reference star(s) from the file.
            print("Given reference star location:")                                                                     # Print a little title describing what the next line will be.
            print(possible_ref_star_loc)                                                                                # Print the SkyCoord object of the possible reference star(s) from the file.
            possible_image_star_loc = skypositions[possible_image_star_index]                                           # Store the RA/dec of the possible image star(s) matching the possible reference star(s).
            print("Detected reference star loation:")                                                                   # Print a little title describing what the next line will be.
            print(possible_image_star_loc)                                                                              # Print the SkyCoord object of the possible image star(s) matching the possible reference star(s).
            if ground_based:
                possible_image_star_altaz = altazpositions[possible_image_star_index]                                       # Store the Alt/Az of the possible image star(s) matching the possible reference star(s).
                print("Detected reference star AltAz:")                                                                     # Print a little title describing what the next line will be.
                print(possible_image_star_altaz)                                                                            # Print the Alt/Az of the possible image star(s) matching the possible reference star(s).
            possible_ref_star_x, possible_ref_star_y = wcs.world_to_pixel(possible_ref_star_loc)                        # Convert the reference star(s) from the file from RA/dec to pixels.
            print(f"Given reference star pixel location: ({possible_ref_star_x}, {possible_ref_star_y})")               # Print the pixel location of the possible reference star(s) from the file.
            possible_image_star_x, possible_image_star_y = wcs.world_to_pixel(possible_image_star_loc)                  # Convert the possible image star(s) location from RA/dec to pixels.
            print(f"Detected reference star pixel location: ({possible_image_star_x}, {possible_image_star_y})")        # Print the pixel location of the possible image star(s) matching the possible reference star(s).
            print("")                                                                                                   # Print a line break.
            # print(f"Angular separation between reference and image star: "
            #       f"{sep2d[possible_ref_star_index].to(u.arcsec):.3f}")                                                 # Print the angular separation between the reference star and possible image star.
            # Plot the locations of the image star(s) matching the possible reference star(s) and the reference star(s)
            # as a confirmation that they are the same star(s).
            ax = plt.subplot(projection=wcs)                                                                            # Create matplotlib axes with a projection onto the RA/dec from the WCS.
            ax.imshow(imgdata, cmap='gray', norm=LogNorm(), interpolation='nearest')                                    # Plot the image onto the projected axes.
            ax.scatter(possible_ref_star_x, possible_ref_star_y,
                       s=199, edgecolor='red', facecolor='none', label='Reference Star from File')                      # Plot a circle where the possible reference star(s) from the file are located in the image.

            ax.scatter(possible_image_star_x, possible_image_star_y,
                       s=199, edgecolor='green', facecolor='none', label='Reference Star from Image')                   # Plot a circle where the possible image star(s) matching the possible reference star(s) are located in the image.
            ax.grid(color='gray', ls='solid')                                                                           # Plot the grid of RA/dec onto the image.
            HIP_title = reference_stars.colnames[0]                                                                     # This stores the name of the firs column to be used to index for the title.
            ref_name = reference_stars[HIP_title][possible_ref_star_index]                                              # The name of the reference star from the reference stars file.
            ax.set_ylabel('Declination (J2000)')                                                                        # Y (dec) axis title.
            ax.set_xlabel('Right Ascension (J2000)')                                                                    # X (RA) axis title.
            # plt.title(f"Confirmation of Detection of Reference Star HIP {ref_name}")                                    # Title includes the name of the possible reference star(s) from the file.
            plt.legend()                                                                                                # Show the legend for which is the image star and which is the reference star from a file.
            # plt.show(block=False)                                                                                       # Show the plot.
            # plt.pause(2)
            plt.close()                                                                                                 # Close the plot. Clears it from RAM.
            print("")                                                                                                   # Print a blank line.
            # This is how it is currently set up to determine if the image star and the star from the file are the same.
            # If the angular separation is less than 10 arcsec, it continues assuming they are the same.
            if max(sep2d[possible_ref_star_index]) < 10 * u.arcsec:
                correct_star_bool = True
            else:
                correct_star_bool = False
            # This is the other option to determine if the image and reference stars match up. This involves user input
            # for each image. It is currently commented out to make debugging faster, as the code can iterate over each
            # image without pausing for user input.

            # y_or_n_input = False                                                                                        # Whether or not the user has input 'y' or 'n' inside the following while loop.
            # while not y_or_n_input:                                                                                     # This loop continues to ask the user until the response is correct.
            #     correct_star = input("Was the correct reference star found in the image? (y/n): ")                      # Ask the user whether or not the image star(s) is the same as the reference star(s).
            #     if correct_star == 'y' or correct_star == 'Y':                                                          # If the correct star was found.
            #         correct_star_bool = True                                                                            # The correct star was found.
            #         y_or_n_input = True                                                                                 # Exits the while loop.
            #     elif correct_star == 'n' or correct_star == 'N':                                                        # If the image star(s) does not match the reference star(s).
            #         correct_star_bool = False                                                                           # The wrong star(s) was found.
            #         y_or_n_input = True                                                                                 # Exits the while loop.
            #     else:                                                                                                   # The user did not enter an acceptable value (y or n).
            #         print("Please enter 'y' or 'n'")                                                                    # Reminds the user to please enter either 'y' or 'n'.
            #         y_or_n_input = False                                                                                # Does not exit the while loop.
            # correct_star_bool = True
            if not correct_star_bool:                                                                                   # If the wrong star was found.
                print("Wrong star. Moving onto next file.")                                                             # Confirms to the user that the wrong star was found.
                continue                                                                                                # Goes to the next iteration of the loop (the next image file).
            print("Hopefully only shows if it is the right star!")                                                      # Confirms that the code works. Will delete later.
            filter_name = hdr['FILTER']                                                                                 # Read the filter name from the FITS header.
            if ground_based:
                image_star_secz = possible_image_star_altaz.secz                                                            # Calculate sec(z). This is used for the airmass.
            ang_separation = sep2d[possible_ref_star_index]                                                             # Store the angular separation in a variable.
            flux = (fluxes[possible_image_star_index] * exptime).value
            print(flux)
            instr_mag = instr_mags[possible_image_star_index]                                                           # Store the instrumental magnitude of the image star(s).
            instr_mag_sigma = instr_mags_sigma[possible_image_star_index]
            filter_name_repeat = np.full(len(instr_mag), filter_name)
            time_repeat = np.full(len(instr_mag), date.jd)
            exposure_repeat = np.full(len(instr_mag), exptime.value)
            # print(f"m = {instr_mag:.3f} +/- {instr_mag_sigma:.3f}")
            img_star_icrs = possible_image_star_loc.icrs                                                                # Store the RA/dec of the image star in ICRS.
            V_apparent = possible_ref_star['V']                                                                         # Store the apparent V band magnitude of the reference star(s).
            try:
                B_V_apparent = possible_ref_star['(B-V)']                                                               # Store the (B-V) index of the reference star(s).
                U_B_apparent = possible_ref_star['(U-B)']                                                               # Store the (U-B) index of the reference star(s).
                V_R_apparent = possible_ref_star['(V-R)']                                                               # Store the (V-R) index of the reference star(s).
                V_I_apparent = possible_ref_star['(V-I)']                                                               # Store the (V-R) index of the reference star(s).
                V_sigma_apparent = possible_ref_star['V_sigma']
            except KeyError as e:
                B_V_apparent = possible_ref_star['B-V']
                U_B_apparent = possible_ref_star['U-B']
                V_R_apparent = possible_ref_star['V-R']
                V_I_apparent = possible_ref_star['V-I']
                V_sigma_apparent = possible_ref_star['e_V']
            # print(instr_mag)                                                                                            # Print the instrumental magnitude of the image star(s).
            # print(V_apparent)                                                                                           # Print the apparent V band magnitude of the reference star(s).
            # print(B_V_apparent)                                                                                         # Print the (B-V) index of the reference star(s).
            # print(U_B_apparent)                                                                                         # Print the (U-B) index of the reference star(s).
            # print(V_R_apparent)                                                                                         # Print the (V-R) index of the reference star(s).
            # print(V_I_apparent)                                                                                         # Print the (V-I) index of the reference star(s).
            # The following append to the lists that were created before the for loop began. These lists will be used to
            # create the large table containing all of the stars' information. Using lists was recommended instead of
            # appending to a table each iteration for performance reasons.
            ref_star_name.extend(list(ref_name))
            times.extend(list(time_repeat))
            ref_star_RA.extend(list(possible_ref_star_loc.ra.to(u.hourangle)))
            ref_star_dec.extend(list(possible_ref_star_loc.dec))
            img_star_RA.extend(list(img_star_icrs.ra.to(u.hourangle)))
            img_star_dec.extend(list(img_star_icrs.dec))
            angular_separation.extend(list(ang_separation.to(u.arcsec)))
            ref_star_x.extend(list(possible_ref_star_x))
            ref_star_y.extend(list(possible_ref_star_y))
            img_star_x.extend(list(possible_image_star_x))
            img_star_y.extend(list(possible_image_star_y))
            ref_star_mag.extend(list(V_apparent))
            img_star_mag.extend(list(instr_mag))
            filters.extend(list(filter_name_repeat))
            flux_table.extend(list(flux))
            exposure.extend(list(exposure_repeat))
            B_V_apparents.extend(list(B_V_apparent))
            U_B_apparents.extend(list(U_B_apparent))
            V_R_apparents.extend(list(V_R_apparent))
            V_I_apparents.extend(list(V_I_apparent))
            V_sigma_apparents.extend(list(V_sigma_apparent))
            if ground_based:
                img_star_secz.append(image_star_secz)                                                                   # The airmass that will be used for calculations.
                # Perhaps change this to check if the difference to this value is less than the difference to one already
                # included in the 'X_rounded' array. This may not work, but is a step that can be tried.
                X_rounded.append(round(image_star_secz.value, 1))
                # X_rounded.append(floor(image_star_secz.value * 10) / 10.0)                                              # Rounds down the airmass to use as an identifier. Not a perfect solution.
# Create the table containing most of the desired star information. This is the table that will be read later to convert
# the data to a format similar to Appendices B and C from A Practical Guide to Photometry.
large_stars_table = QTable(
    data=[
        ref_star_name,
        times,
        ref_star_RA,
        ref_star_dec,
        img_star_RA,
        img_star_dec,
        angular_separation,
        ref_star_x,
        ref_star_y,
        img_star_x,
        img_star_y,
        flux_table,
        exposure,
        img_star_mag,
        filters,
        ref_star_mag,
        B_V_apparents,
        U_B_apparents,
        V_R_apparents,
        V_I_apparents,
        V_sigma_apparents
    ],
    names=[
        'Name',
        'Time (JD)',
        'RA_ref',
        'dec_ref',
        'RA_img',
        'dec_img',
        'angular_separation',
        'x_ref',
        'y_ref',
        'x_img',
        'y_img',
        'flux',
        'exposure',
        'mag_instrumental',
        'filter',
        'V',
        '(B-V)',
        '(U-B)',
        '(V-R)',
        '(V-I)',
        'V_sigma'
    ]
)
# large_stars_table = QTable(
#     data=[
#         ref_star_name,
#         ref_star_RA,
#         ref_star_dec,
#         img_star_RA,
#         img_star_dec,
#         angular_separation,
#         ref_star_x,
#         ref_star_y,
#         img_star_x,
#         img_star_y,
#         img_star_mag,
#         filters,
#         ref_star_mag,
#         B_V_apparents,
#         U_B_apparents,
#         V_R_apparents,
#         V_I_apparents,
#         img_star_secz,
#         X_rounded
#     ],
#     names=[
#         'HIP',
#         'RA_ref',
#         'dec_ref',
#         'RA_img',
#         'dec_img',
#         'angular_separation',
#         'x_ref',
#         'y_ref',
#         'x_img',
#         'y_img',
#         'mag_instrumental',
#         'filter',
#         'V',
#         '(B-V)',
#         '(U-B)',
#         '(V-R)',
#         '(V-I)',
#         'X',
#         'X_rounded'
#     ]
# )
large_stars_table.pprint_all()                                                                                          # Print the large stars table.
# unique_stars = table.unique(large_stars_table, keys=['HIP', 'X_rounded'])                                               # Get all unique combinations of star names and X_rounded in the large stars table.
unique_stars = table.unique(large_stars_table, keys=['Name'])
unique_stars.pprint_all()                                                                                               # Print the information of the unique observations.
N = len(unique_stars)                                                                                                   # Get the number of unique observations. Used to preallocate the table.
# Preallocate the table containing the averages of each observation. The format is similar to how the data was presented
# in Appendices B and C of A Practical Guide. This currently assumes observations in the BGR bands.
nan_array = np.empty(N)
nan_array.fill(np.nan)
if ground_based:
    stars_table = Table(
        names=[
            'Name',
            'b',
            'g',
            'r',
            'b_sigma',
            'g_sigma',
            'r_sigma',
            'X_b',
            'X_g',
            'X_r',
            'X_b_sigma',
            'X_g_sigma',
            'X_r_sigma',
            'V',
            '(B-V)',
            '(U-B)',
            '(V-R)',
            '(V-I)',
            'X_rounded'
        ],
        data=[
            np.empty(N, dtype=object),
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array
        ]
    )
else:
    stars_table = Table(
        names=[
            'Name',
            'flux',
            'exposure',
            'clear',
            'clear_sigma',
            'V',
            '(B-V)',
            '(U-B)',
            '(V-R)',
            '(V-I)',
            'V_sigma'
        ],
        data=[
            np.empty(N, dtype=object),
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array,
            nan_array
        ]
    )
# Write the apparent magnitudes and colour indices to the stars table as those will not change.
stars_table['V'] = unique_stars['V']
stars_table['(B-V)'] = unique_stars['(B-V)']
stars_table['(U-B)'] = unique_stars['(U-B)']
stars_table['(V-R)'] = unique_stars['(V-R)']
stars_table['(V-I)'] = unique_stars['(V-I)']
stars_table['V_sigma'] = unique_stars['V_sigma']
stars_table['exposure'] = unique_stars['exposure']
if ground_based:
    stars_table['X_rounded'] = unique_stars['X_rounded']
i = 0                                                                                                                   # Initialize the iterator that will be used for indexing.
for star in unique_stars['Name']:                                                                                       # Loop over every unique observation.
    stars_table['Name'][i] = star                                                                                       # Write the star name to the stars table.
    mask = large_stars_table['Name'] == star                                                                            # Indices of all entries in the large stars table that match the current star name.
    current_star_table = large_stars_table[mask]                                                                        # Create a new table using the entries that satisfy the mask created previously.
    current_star_table.pprint_all()                                                                                     # Print the entries that match the current star name.
    unique_filters = table.unique(current_star_table, keys='filter')                                                    # Get all combinations of filters for the current star.
    unique_filters.pprint_all()                                                                                         # Print the entries containing the first instance of each unique filter.
    if ground_based:
        unique_X_rounded = unique_stars['X_rounded'][i]                                                                 # Get the value of X_rounded for the current position in the iteration over the unique stars table.
    for unique_filter in unique_filters['filter']:                                                                      # Loop over every different filter.
        if ground_based:
            mask = ((current_star_table['filter'] == unique_filter) & (current_star_table['X_rounded'] == unique_X_rounded))# Indicies of all entries in the current star table that match the current filter AND X_rounded.
        else:
            mask = ((current_star_table['filter'] == unique_filter))
        current_star_filter_table = current_star_table[mask]                                                            # Create a new table using the entries that satisfy the mask created previously.
        current_star_filter_table.pprint_all()                                                                          # Print the entries that match the current filter and X_rounded.
        # print(unique_filter)                                                                                            # Print the name of the filter.
        # print(unique_X_rounded)                                                                                         # Print the current X_rounded.
        fluxes_numpy = np.array(current_star_filter_table['flux'])
        mean_flux = fluxes_numpy.mean()
        std_flux = fluxes_numpy.std()
        exposure_numpy = np.array(current_star_filter_table['exposure'])
        mean_exposure = exposure_numpy.mean()
        std_exposure = exposure_numpy.std()
        mags_numpy = np.array(current_star_filter_table['mag_instrumental'])                                            # Store the instrumental magnitudes for the current filter at the current airmass as a numpy array.
        mean_mag = mags_numpy.mean()                                                                                    # Calculate the mean of the instrumental magnitude numpy array.
        std_mag = mags_numpy.std()                                                                                      # Calculate the standard deviation of the instrumental magnitude numpy array.
        filter_column = unique_filter.lower()                                                                           # Convert the current filter to a lowercase string.
        sigma_column = f'{filter_column}_sigma'                                                                         # Column name where the standard deviation of the instrumental magnitude will be stored.
        stars_table[filter_column][i] = mean_mag                                                                        # Store the mean instrumental magnitude.
        stars_table[sigma_column][i] = std_mag                                                                          # Store the standard deviation of the instrumental magnitude.
        stars_table['flux'][i] = mean_flux
        # stars_table['exposure'][i] = mean_exposure
        if ground_based:
            X_numpy = np.array(current_star_filter_table['X'])                                                          # Store the airmasses for the current filter of the current observation as a numpy array.
            mean_X = X_numpy.mean()                                                                                     # Calculate the mean of the airmass numpy array.
            std_X = X_numpy.std()                                                                                       # Calculate the standard deviation of the airmass numpy array.
            X_column = f'X_{filter_column}'                                                                             # Column name where the mean airmass of the current filter/observation combination will be stored.
            X_std_column = f'X_{filter_column}_sigma'                                                                   # Column name where the standard deviation airmass of the current filter/observation combination will be stored.
            stars_table[X_column][i] = mean_X                                                                           # Store the mean airmass.
            stars_table[X_std_column][i] = std_X                                                                        # Store the standard deviation of the airmass.
    i += 1                                                                                                              # Increase the iterator.
stars_table.pprint_all()                                                                                                # Print the formatted stars table.
# ascii.write(stars_table, r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Landolt Stars\Detected_stars.csv',
#             format='csv', delimiter=',')

# ascii.write(stars_table, 'stars_table_test.csv', format='csv')

# I think I need a way to not include duplicate stars.
if ground_based:
    minimum_x_rounded = min(stars_table['X_rounded'])
    min_airmass_cutoff = 25
    multiply_min_airmass_cutoff = 1 + (min_airmass_cutoff / 100.0)
    mask = stars_table['X_rounded'] <= multiply_min_airmass_cutoff * minimum_x_rounded
    transform_table = stars_table[mask]
    if len(stars_table) > 1:
        while len(transform_table) < 2:
            print(f"No 2 stars within {min_airmass_cutoff}% of the lowest airmass.")
            min_airmass_cutoff += 5
            print(f"Trying again with {min_airmass_cutoff}%...")
            multiply_min_airmass_cutoff = 1 + (min_airmass_cutoff / 100.0)
            mask = stars_table['X_rounded'] <= multiply_min_airmass_cutoff * minimum_x_rounded
            transform_table = stars_table[mask]
else:
    transform_table = stars_table
# transform_table.pprint_all()
# Calculate the transforms.
max_V_sigma = max(stars_table['V_sigma'])
max_clear_sigma = max(stars_table['clear_sigma'])
bv_plot = np.arange(start=min(stars_table['(B-V)'])-0.2, stop=max(stars_table['(B-V)'])+0.2, step=0.1)
err_sum = np.nan_to_num(stars_table['V_sigma'], nan=max_V_sigma) + \
          np.nan_to_num(stars_table['clear_sigma'], nan=max_clear_sigma)
err_sum = np.array(err_sum)
# print(err_sum)
# c_vbv, zprime_v = np.polyfit(stars_table['(B-V)'], stars_table['V'] - stars_table['g'], 1)
# plt.plot(stars_table['(B-V)'], stars_table['V'] - stars_table['g'], 'o')
[c_vbv, zprime_vbv], residual, _, _, _ = np.polyfit(stars_table['(B-V)'], stars_table['V'] - stars_table['clear'], 1,
                                                  full=True, w=1/err_sum)
# [c_vbv, zprime_v], residual, _, _, _ = np.polyfit(stars_table['(B-V)'], stars_table['V'] - stars_table['clear'], 1,
#                                                   full=True)
# c_vbv, zprime_v = np.polyfit(stars_table['(B-V)'], stars_table['V'] - stars_table['clear'], 1)
# plt.plot(stars_table['(B-V)'], stars_table['V'] - stars_table['clear'], 'o')
plt.errorbar(stars_table['(B-V)'], stars_table['V'] - stars_table['clear'], yerr=err_sum, fmt='o', capsize=2)
plt.plot(bv_plot, c_vbv * bv_plot + zprime_vbv)
plt.ylabel("V-c")
plt.xlabel("(B-V)")
plt.title(f"(V-c) = {c_vbv:.3f} * (B-V) + {zprime_vbv:.3f}")
plt.show()
# plt.show(block=False)
# plt.pause(3)
plt.close()
# print(f"T_v = {tv}")
# print(f"Z_p (assuming k=0) = {zp_k0}")
print(residual)

vi_plot = np.arange(start=min(stars_table['(V-I)'])-0.2, stop=max(stars_table['(V-I)'])+0.2, step=0.1)
# print(err_sum)
# c_vvi, zprime_vvi = np.polyfit(stars_table['(B-V)'], stars_table['V'] - stars_table['g'], 1)
# plt.plot(stars_table['(V-I)'], stars_table['V'] - stars_table['g'], 'o')
[c_vvi, zprime_vvi], residual, _, _, _ = np.polyfit(stars_table['(V-I)'], stars_table['V'] - stars_table['clear'], 1,
                                                  full=True, w=1/err_sum)
# [c_vvi, zprime_vvi], residual, _, _, _ = np.polyfit(stars_table['(V-I)'], stars_table['V'] - stars_table['clear'], 1,
#                                                   full=True)
# c_vvi, zprime_vvi = np.polyfit(stars_table['(V-I)'], stars_table['V'] - stars_table['clear'], 1)
# plt.plot(stars_table['(V-I)'], stars_table['V'] - stars_table['clear'], 'o')
plt.errorbar(stars_table['(V-I)'], stars_table['V'] - stars_table['clear'], yerr=err_sum, fmt='o', capsize=2)
plt.plot(vi_plot, c_vvi * vi_plot + zprime_vvi)
plt.ylabel("V-c")
plt.xlabel("(V-I)")
plt.title(f"(V-c) = {c_vvi:.3f} * (V-I) + {zprime_vvi:.3f}")
plt.show()
# plt.show(block=False)
# plt.pause(3)
plt.close()

vr_plot = np.arange(start=min(stars_table['(V-R)'])-0.2, stop=max(stars_table['(V-R)'])+0.2, step=0.1)
# print(err_sum)
# c_vvr, zprime_vvr = np.polyfit(stars_table['(V-R)'], stars_table['V'] - stars_table['g'], 1)
# plt.plot(stars_table['(V-R)'], stars_table['V'] - stars_table['g'], 'o')
[c_vvr, zprime_vvr], residual, _, _, _ = np.polyfit(stars_table['(V-R)'], stars_table['V'] - stars_table['clear'], 1,
                                                  full=True, w=1/err_sum)
# [c_vvr, zprime_vvr], residual, _, _, _ = np.polyfit(stars_table['(V-R)'], stars_table['V'] - stars_table['clear'], 1,
#                                                   full=True)
# c_vvr, zprime_vvr = np.polyfit(stars_table['(V-R)'], stars_table['V'] - stars_table['clear'], 1)
# plt.plot(stars_table['(V-R)'], stars_table['V'] - stars_table['clear'], 'o')
plt.errorbar(stars_table['(V-R)'], stars_table['V'] - stars_table['clear'], yerr=err_sum, fmt='o', capsize=2)
plt.plot(vr_plot, c_vvr * vr_plot + zprime_vvr)
plt.ylabel("V-c")
plt.xlabel("(V-R)")
plt.title(f"(V-c) = {c_vvr:.3f} * (V-R) + {zprime_vvr:.3f}")
plt.show()
# plt.show(block=False)
# plt.pause(3)
plt.close()

if ground_based:
    x_plot = np.arange(start=1, stop=3.1, step=0.1)
    # t_vbv, k_doubleprime_vbv = np.polyfit(stars_table['X_g'], )

    # Calculate the extinction (without the transforms).

    x_plot = np.arange(start=1, stop=4, step=0.1)
    # k_g, offset_g = np.polyfit(stars_table['X_g'], stars_table['V'] - stars_table['g'], 1)
    # k_b, offset_b = np.polyfit(stars_table['X_b'], stars_table['V'] + stars_table['(B-V)'] - stars_table['b'], 1)
    # k_r, offset_r = np.polyfit(stars_table['X_r'], stars_table['V'] - stars_table['(V-R)'] - stars_table['r'], 1)
    # plt.plot(stars_table['X_g'], stars_table['V'] - stars_table['g'], 'go')
    # plt.plot(stars_table['X_b'], stars_table['V'] + stars_table['(B-V)'] - stars_table['b'], 'bo')
    # plt.plot(stars_table['X_r'], stars_table['V'] - stars_table['(V-R)'] - stars_table['r'], 'ro')
    k_g, offset_g = np.polyfit(stars_table['X_g'][~np.isnan(stars_table['X_g'])], stars_table['g'][~np.isnan(stars_table['g'])], 1)
    k_b, offset_b = np.polyfit(stars_table['X_b'][~np.isnan(stars_table['X_b'])], stars_table['b'][~np.isnan(stars_table['b'])], 1)
    k_r, offset_r = np.polyfit(stars_table['X_r'][~np.isnan(stars_table['X_r'])], stars_table['r'][~np.isnan(stars_table['r'])], 1)
    plt.plot(stars_table['X_g'], stars_table['g'], 'go')
    plt.plot(stars_table['X_b'], stars_table['b'], 'bo')
    plt.plot(stars_table['X_r'], stars_table['r'], 'ro')
    plt.plot(x_plot, k_b*x_plot+offset_b, 'b', label=f"b={k_b:.3f}X + {offset_b:.3f}")
    plt.plot(x_plot, k_r*x_plot+offset_r, 'r', label=f"r={k_r:.3f}X + {offset_r:.3f}")
    plt.plot(x_plot, k_g*x_plot+offset_g, 'g', label=f"g={k_g:.3f}X + {offset_g:.3f}")
    plt.title("Extinction Plot 21 Mar")
    plt.gca().invert_yaxis()
    plt.ylabel("Instrumental Magnitude")
    plt.xlabel("Airmass")
    plt.legend()
    # plt.show()
    plt.close()

# Calculate the extinction for each star.
# diff_stars = table.unique(stars_table, keys='HIP')
# k_gs = []
# k_bs = []
# k_rs = []
# k_g_repeat = np.full(len(diff_stars), abs(k_g))
# k_b_repeat = np.full(len(diff_stars), abs(k_b))
# k_r_repeat = np.full(len(diff_stars), abs(k_r))
# for star in diff_stars['HIP']:
#     mask = stars_table['HIP'] == star
#     current_star = stars_table[mask]
#     k_g_star, offset_g = np.polyfit(current_star['X_g'], current_star['V'] - current_star['g'], 1)
#     k_b_star, offset_b = np.polyfit(current_star['X_b'], current_star['V'] + current_star['(B-V)'] - current_star['b'], 1)
#     k_r_star, offset_r = np.polyfit(current_star['X_r'], current_star['V'] - current_star['(V-R)'] - current_star['r'], 1)
#     k_gs.append(abs(k_g_star))
#     k_bs.append(abs(k_b_star))
#     k_rs.append(abs(k_r_star))
# k_gs_avg = np.full(len(diff_stars), np.mean(k_gs))
# k_bs_avg = np.full(len(diff_stars), np.mean(k_bs))
# k_rs_avg = np.full(len(diff_stars), np.mean(k_rs))
# k_gs_std = np.full(len(diff_stars), np.std(k_gs))
# k_bs_std = np.full(len(diff_stars), np.std(k_bs))
# k_rs_std = np.full(len(diff_stars), np.std(k_rs))
# extinction_table = Table(
#     data=[
#         diff_stars['HIP'],
#         k_gs,
#         k_g_repeat,
#         k_gs_avg,
#         k_gs_std,
#         k_bs,
#         k_b_repeat,
#         k_bs_avg,
#         k_bs_std,
#         k_rs,
#         k_r_repeat,
#         k_rs_avg,
#         k_rs_std
#     ],
#     names=[
#         'HIP',
#         'k_g_per_star',
#         'k_g_per_night',
#         'k_g_avg_per_star',
#         'k_g_sigma_per_star',
#         'k_b_per_star',
#         'k_b_per_night',
#         'k_b_avg_per_star',
#         'k_b_sigma_per_star',
#         'k_r_per_star',
#         'k_r_per_night',
#         'k_r_avg_per_star',
#         'k_r_sigma_per_star'
#     ]
# )
# extinction_table.pprint_all()

# ascii.write(extinction_table, 'extinction_table_10Mar.csv', format='csv')