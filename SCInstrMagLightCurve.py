from astropy.io import fits, ascii
import astropy.units as u
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma, SigmaClip
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.table import Table, QTable, hstack, unique
from astropy.time import Time
import datetime
from photutils.background import Background2D, SExtractorBackground
from photutils.detection import IRAFStarFinder
from photutils.psf import DAOGroup, BasicPSFPhotometry, IntegratedGaussianPRF
import numpy as np
import os
from math import sqrt, atan
import matplotlib
from matplotlib import patches
from matplotlib.colors import LogNorm
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import ctypes
import win32com.client as win32
from tkinter import *
from matplotlib.lines import Line2D
from shutil import copy2, rmtree
import cv2 as cv
from photutils.aperture import RectangularAperture
matplotlib.use('TkAgg')

# The directory where the Intelsat 10-02 files are stored.
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-03-20 - Calibrated\Intelsat 10-02'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021_J132_46927_DESCENT\2021_J132_46927_DESCENT\May 11 2021'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021_J132_46927_DESCENT\May 18 2021\46927'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021 02 07\2021 02 07 - Intelsat 10-02 - G Band'
directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\Intelsat 10-02\2021-04-25 - Calibrated - Intelsat 10-02'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\NEOSSat Observations\2016-111\2016-111'
# stars_directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-03-20 - Calibrated\Zpoint Test'
# directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-04-21\Intelsat 10-02 ALL'
# stars_directory = r'C:\Users\jmwawrow\Documents\DRDC_Code\2021-04-21\Zpoint Test'
# catloc = r'C:\Program Files (x86)\PinPoint\UCAC4'

date_string = '07 Feb'

# b_zpoints = []
# g_zpoints = []
# r_zpoints = []

# for dirpath, dirnames, filenames in os.walk(stars_directory):
#     for filename in filenames:
#         if filename.endswith(".fits"):
#             filepath = os.path.join(dirpath, filename)
#             # print(filepath)d
#             with fits.open(filepath) as image:
#                 hdr = image[0].header
            
#             focal_length = hdr['FOCALLEN'] * u.mm
#             xpixsz = hdr['XPIXSZ'] * u.um
#             x_rad_per_pix = atan(xpixsz / focal_length) * u.rad
#             x_arcsec_per_pix_units = x_rad_per_pix.to(u.arcsec)
#             x_arcsec_per_pix = x_arcsec_per_pix_units.value
#             ypixsz = hdr['XPIXSZ'] * u.um
#             y_rad_per_pix = atan(ypixsz / focal_length) * u.rad
#             y_arcsec_per_pix_units = y_rad_per_pix.to(u.arcsec)
#             y_arcsec_per_pix = y_arcsec_per_pix_units.value
            
#             instr_filter = hdr['FILTER']
            
#             f = win32.Dispatch("Pinpoint.plate")
#             f.DetachFITS
            
#             f.AttachFITS(filepath)
#             #print(c[o])
#             f.Declination = f.targetDeclination
#             f.RightAscension = f.targetRightAscension
#             f.ArcsecperPixelHoriz  = x_arcsec_per_pix
#             f.ArcsecperPixelVert = y_arcsec_per_pix
#             if instr_filter == 'B':
#                 f.ColorBand = 1
#             elif instr_filter == 'V' or instr_filter == 'G':
#                 f.ColorBand = 2
#             elif instr_filter == 'R':
#                 f.ColorBand = 3
#             f.Catalog = 11
#             f.CatalogPath = catloc
#             f.CatalogMaximumMagnitude = 13
#             f.CatalogExpansion = 0.8
#             f.SigmaAboveMean = 3.0
#             f.FindImageStars
#             f.FindCatalogStars
#             f.MaxSolveTime = 60
#             f.MaxMatchResidual = 1.5
#             flag = 0
#             f.FindCatalogStars()
#             try:
#                 f.Solve()
#             except Exception as e:
#                 print(e)
#                 continue

#             zpoint = f.MagZeroPoint
            
#             f.DetachFITS()
#             if instr_filter == 'B':
#                 b_zpoints.append(zpoint)
#             elif instr_filter == 'V' or instr_filter == 'G':
#                 g_zpoints.append(zpoint)
#             elif instr_filter == 'R':
#                 r_zpoints.append(zpoint)
            
# b_zpoints = np.array(b_zpoints)
# g_zpoints = np.array(g_zpoints)
# r_zpoints = np.array(r_zpoints)

# b_zpoint = b_zpoints.mean()
# b_zpoint_std = b_zpoints.std()
# print(f"B band ZMag = {b_zpoint:.3f} +/- {b_zpoint_std:.3f}")
# g_zpoint = g_zpoints.mean()
# g_zpoint_std = g_zpoints.std()
# print(f"G band ZMag = {g_zpoint:.3f} +/- {g_zpoint_std:.3f}")
# r_zpoint = r_zpoints.mean()
# r_zpoint_std = r_zpoints.std()
# print(f"R band ZMag = {r_zpoint:.3f} +/- {r_zpoint_std:.3f}")

b_zpoint = 0
b_zpoint_std = 0
print(f"B band ZMag = {b_zpoint:.3f} +/- {b_zpoint_std:.3f}")
g_zpoint = 0
g_zpoint_std = 0
print(f"G band ZMag = {g_zpoint:.3f} +/- {g_zpoint_std:.3f}")
r_zpoint = 0
r_zpoint_std = 0
print(f"R band ZMag = {r_zpoint:.3f} +/- {r_zpoint_std:.3f}")

"""
#### For the mouse click: ####
# Initial satellite position.
while defining_sats_initially == True:
    Have opencv display the image and get the mouse click position.
    Display the box that was just created.
    Have the user enter the name of the sat specified.
    Confirm if the user wants to add more sats.
    if user doesn't want more sats:
        defining_sats_initially = False
if X nans in a row (when objects were detected):
    change_sats_position = True
while change_sats_position:
    Ask the user which sat to change the position for (A drop-down menu?)
    Use opencv to set the new position like earlier.
    Confirm if the user wants to change the position of more sats.
    if user doesn't want to edit more sats:
        change_sats_position = False
        Go back X images (will need to find out how to do that)
"""
set_sat_positions = True
cmap = plt.get_cmap('Set1')

def return_entry(event=None):
    """Gets and prints the content of the entry"""
    # global entry
    global content
    content = entry.get()
    root.destroy()


def set_sat_position(event, x, y, flags, params):
    global sat_locs
    if event == cv.EVENT_LBUTTONDOWN:
        sat_locs.append([x, y])


def change_sat_position(event, x, y, flags, params):
    global sat_locs
    if event == cv.EVENT_LBUTTONDOWN:
        sat_locs[params[0]] = [x, y]
        cv.destroyAllWindows()


def mbox(title, text, style):
    return ctypes.windll.user32.MessageBoxW(0, text, title, style)


# Controls if the plots will be shown. Useful to be able to control when debugging.
# 0 - Don't show the plots of each image.
# 1 - Show the plots for 3 seconds.
# 2 - Show the plots until closed manually.
plot_results = 0
size = 25                                                                                                               # Size of the cutout to plot and fit the gaussian to (pixels).
hsize = int((size - 1) / 2)                                                                                             # Half of the size of the cutout.
fitter = LevMarLSQFitter()                                                                                              # Initialize the fitter that will be used to fit the 1D Gaussian.
max_distance_from_sat = 20
max_num_nan = 5
num_nan = 0
change_sat_positions = False
none_sats = False
sigma_clip = SigmaClip(sigma=3.)
bkg_estimator = SExtractorBackground()
debugging = True
temp_dir = 'tmp'
if not debugging:
    try:
        os.mkdir(temp_dir)
    except FileExistsError as e:
        print(e)
else:
    if os.path.exists(temp_dir):
        rmtree(temp_dir)
    os.mkdir(temp_dir)

filecount = 0
for dirpth, _, files in os.walk(directory):
    for file in files:
        if file.endswith(".fits"):
            with fits.open(os.path.join(dirpth, file)) as image:
                hdr = image[0].header
            t = Time(hdr['DATE-OBS'], format='fits', scale='utc')
            try:
                filter_name = hdr['FILTER']
            except KeyError:
                filter_name = 'C'
            t_datetime = t.to_datetime()
            new_filename = f'{t_datetime.strftime("%Y%m%d%H%M%S")}{filter_name}.fits'
            copy2(os.path.join(dirpth, file), f'{temp_dir}/{new_filename}')
            filecount += 1

filenames = sorted(os.listdir(temp_dir))
for filename in filenames:
    print(filename)

for filenum, file in enumerate(filenames):
    with fits.open(f"{temp_dir}/{file}") as image:
        hdr = image[0].header                                                                                       # Store the fits header as a variable.
        imgdata = image[0].data                                                                                     # Store the image as a variable.
    print(file)
    # box_size_x = int(hdr['NAXIS1'] / 31)
    # box_size_y = int(hdr['NAXIS2'] / 31)
    # bkg = Background2D(imgdata, (box_size_x, box_size_y), filter_size=(3, 3),
    #                    sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    # imgdata = imgdata - bkg.background
    while set_sat_positions:
        sat_locs = []
        mbox('Information',
             'Please select the positions of the satellites on the following image. Press any key when finished.',
             0)
        cv.namedWindow('TestImage')
        cv.setMouseCallback('TestImage', set_sat_position)
        logdata = cv.normalize(imgdata, None, alpha=0, beta=1, norm_type=cv.NORM_MINMAX, dtype=cv.CV_32F)
        cv.imshow('TestImage', logdata)
        cv.waitKey(0)
        cv.destroyAllWindows()
        sat_locs = np.array(sat_locs)
        print(sat_locs)
        num_sats = len(sat_locs)
        num_nans = np.zeros(num_sats, dtype=int)
        names = np.empty(num_sats + 2, dtype=object)
        names[0] = 'Time (JD)'
        names[1] = 'Filter'
        date_col = np.empty((filecount, 1))
        date_col.fill(np.nan)
        filter_col = np.empty((filecount, 1), dtype=object)
        data = np.empty((filecount, num_sats))
        data.fill(np.nan)

        for i, name in enumerate(names[2:]):
            fig = Figure()
            ax = fig.add_subplot()
            ax.imshow(imgdata, cmap='gray', norm=LogNorm(), interpolation='nearest')
            sat_aperture = RectangularAperture(sat_locs[i], w=max_distance_from_sat * 2,
                                               h=max_distance_from_sat * 2)
            sat_aperture.plot(axes=ax, color='r', lw=1.5, alpha=0.5)
            root = Tk()
            root.title("Set Satellite Position")
            img_frame = Frame(root)
            img_frame.pack(side=LEFT, fill=BOTH, expand=1)
            input_frame = Frame(root)
            input_frame.pack(side=RIGHT)
            label = Label(input_frame, text='Enter Satellite')
            label.pack()
            entry = Entry(input_frame)
            entry.bind("<Return>", return_entry)
            entry.pack(padx=5)
            button = Button(input_frame, text="OK", command=return_entry)
            button.pack()
            canvas = FigureCanvasTkAgg(fig, master=img_frame)
            canvas.draw()
            toolbar = NavigationToolbar2Tk(canvas, img_frame)
            toolbar.update()
            canvas.get_tk_widget().pack(side=LEFT, fill=BOTH, expand=1)
            root.update()
            root.focus_force()
            entry.focus_set()
            root.mainloop()
            names[i + 2] = content
            print(f"Satellite {names[i + 2]} at location ({sat_locs[i, 0]}, {sat_locs[i, 1]})")
        print(names)

        colours = [cmap(i) for i in range(0, num_sats)]
        legend_elements = []
        window = Tk()
        window.title('Plotting in Tkinter Test')
        img_frame = Frame(window)
        img_frame.pack(side=LEFT, fill=BOTH, expand=1)
        input_frame = Frame(window)
        input_frame.pack(side=RIGHT)
        label = Label(input_frame, text='Are the satellite positions correct?')
        label.pack()
        yes_no = IntVar()
        yes_btn = Radiobutton(input_frame, text='Yes', variable=yes_no, value=1)
        yes_btn.pack(anchor=W, padx=5)
        no_btn = Radiobutton(input_frame, text='No', variable=yes_no, value=2)
        no_btn.pack(anchor=W, padx=5)
        closebutton = Button(input_frame, text='OK', command=window.destroy)
        closebutton.pack()
        fig = Figure()
        ax = fig.add_subplot()
        ax.imshow(imgdata, cmap='gray', norm=LogNorm(), interpolation='nearest')
        for i in range(0, num_sats):
            sat_aperture = RectangularAperture(sat_locs[i], w=max_distance_from_sat * 2,
                                               h=max_distance_from_sat * 2)
            sat_aperture.plot(axes=ax, color=colours[i], lw=1.5, alpha=0.5)
            legend_elements.append(Line2D([0], [0], color='w', marker='s', markerfacecolor=colours[i], markersize=7,
                                          label=names[i + 2]))
        fig.legend(handles=legend_elements, framealpha=1)
        canvas = FigureCanvasTkAgg(fig, master=img_frame)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, img_frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        window.mainloop()
        # Have a way for the user to confirm the satellite locations. If it is wrong, then decide whether to change
        # set_sat_positions to True/False or change_sat_positions
        if yes_no.get() == 1:
            set_sat_positions = False
        else:
            continue
        # elif yes_no.get() == 2:
        #     sat_locs = []
        #     names = []
        # print(names[0])
        # print(date_col)
        date_table = Table(names=[names[0]], data=date_col)
        filter_table = Table(names=[names[1]], data=filter_col)
        data_table = Table(names=names[2:], data=data)
        sats_table = hstack([date_table, filter_table, data_table], join_type='exact')
        uncertainty_table = hstack([date_table, filter_table, data_table], join_type='exact')
        sat_fwhm_table = hstack([date_table, filter_table, data_table], join_type='exact')
        sats_table.pprint_all()

    sat_names = names[2:]
    t = Time(hdr['DATE-OBS'], format='fits', scale='utc')
    sats_table['Time (JD)'][filenum] = t.jd
    try:
        sats_table['Filter'][filenum] = hdr['FILTER']
    except KeyError:
        sats_table['Filter'][filenum] = 'C'
    uncertainty_table['Time (JD)'][filenum] = t.jd
    try:
        uncertainty_table['Filter'][filenum] = hdr['FILTER']
    except KeyError:
        uncertainty_table['Filter'][filenum] = 'C'
    sat_fwhm_table['Time (JD)'][filenum] = t.jd
    try:
        sat_fwhm_table['Filter'][filenum] = hdr['FILTER']
    except KeyError:
        sat_fwhm_table['Filter'][filenum] = 'C'
    if change_sat_positions:
        # Change the position
        # Display filenames[filenum - num_nan]
        print(filenames[filenum - num_nan])
        # sat_checked = np.zeros(num_sats, dtype=IntVar())
        with fits.open(f"{temp_dir}/{filenames[filenum - num_nan]}") as image:
            hdr = image[0].header
            imgdata = image[0].data
        root = Tk()
        root.title("Current satellite positions")
        input_frame = Frame(root)
        input_frame.pack(side=RIGHT)
        img_frame = Frame(root)
        img_frame.pack(side=LEFT, fill=BOTH, expand=1)
        label = Label(input_frame, text='Select the satellite(s) whose position you would like to change.')
        label.grid(row=0)
        sat_checked = []
        for sat_num, sat in enumerate(sat_names):
            sat_checked.append(IntVar())
            checkbutton = Checkbutton(input_frame, text=sat, variable=sat_checked[sat_num])
            checkbutton.grid(row=sat_num+1, sticky=W, padx=5)
        none_select = IntVar()
        checkbutton = Checkbutton(input_frame, text="None", variable=none_select)
        checkbutton.grid(row=num_sats+2, sticky=W, padx=5)
        closebutton = Button(input_frame, text='OK', command=root.destroy)
        closebutton.grid(row=num_sats+3)
        legend_elements = []
        fig = Figure()
        ax = fig.add_subplot()
        ax.imshow(imgdata, cmap='gray', norm=LogNorm(), interpolation='nearest')
        for i in range(0, num_sats):
            sat_aperture = RectangularAperture(sat_locs[i], w=max_distance_from_sat * 2,
                                               h=max_distance_from_sat * 2)
            sat_aperture.plot(axes=ax, color=colours[i], lw=1.5, alpha=0.5)
            legend_elements.append(Line2D([0], [0], color='w', marker='s', markerfacecolor=colours[i], markersize=7,
                                          label=names[i + 2]))
        fig.legend(handles=legend_elements, framealpha=1)
        canvas = FigureCanvasTkAgg(fig, master=img_frame)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, img_frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=LEFT, fill=BOTH, expand=1)
        root.mainloop()
        if none_select.get() == 1:
            none_sats = True
        elif none_select.get() == 0:
            none_sats = False
        if not none_sats:
            sat_checked_int = np.empty(len(sat_checked), dtype=int)
            for sat_num, sat in enumerate(sat_checked):
                sat_checked_int[sat_num] = sat.get()
            print(sat_checked_int)
            sat_checked_mask = sat_checked_int == 1
            print(sat_checked_mask)
            for sat in sat_names[sat_checked_mask]:
                print(sat)
                index = np.where(sat_names == sat)
                mbox('Information',
                     f'Please select the new position of {sat} on the following image.',
                     0)
                cv.namedWindow('TestImage')
                cv.setMouseCallback('TestImage', change_sat_position, index)
                logdata = cv.normalize(imgdata, None, alpha=0, beta=1, norm_type=cv.NORM_MINMAX, dtype=cv.CV_32F)
                cv.imshow('TestImage', logdata)
                cv.waitKey(0)
            print(sat_locs)
            # If some selection is none
            # none_sats = True
            for reversing_index in range(1, num_nan+1):
                with fits.open(f"{temp_dir}/{filenames[filenum - reversing_index]}") as image:
                    hdr = image[0].header
                    imgdata = image[0].data
                print(filenames[filenum - reversing_index])
                # Do the photometry on the images. Should change this to a bunch of function calls, rather than copy
                # and pasting the code.
                exptime = hdr['EXPTIME'] * u.s  # Store the exposure time with unit seconds.
                mean_val, median_val, std_val = sigma_clipped_stats(imgdata)  # Calculate background stats.
                iraffind = IRAFStarFinder(threshold=4*std_val, fwhm=2)  # Find stars using IRAF.
                irafsources = iraffind(imgdata - median_val)  # Subtract background median value.
                try:
                    irafsources.sort('flux', reverse=True)  # Sort the stars by flux, from greatest to least.
                except Exception as e:
                    # Reset the NaN boolean as it doesn't count.
                    num_nans[:] = 0
                    continue
                irafpositions = np.transpose((irafsources['xcentroid'],
                                              irafsources['ycentroid']))  # Store source positions as a numpy array.
                if plot_results != 0:
                    fig, ax = plt.subplots()
                    ax.imshow(imgdata, cmap='gray', norm=LogNorm(), interpolation='nearest')
                    ax.scatter(irafsources['xcentroid'], irafsources['ycentroid'],
                               s=100, edgecolor='red', facecolor='none')
                    for i in range(0, num_sats):
                        rect = patches.Rectangle(
                            (sat_locs[i, 0] - max_distance_from_sat, sat_locs[i, 1] - max_distance_from_sat),
                            width=max_distance_from_sat * 2, height=max_distance_from_sat * 2,
                            edgecolor='green', facecolor='none')
                        ax.add_patch(rect)
                    plt.title(os.path.basename(entry.path))
                    if plot_results == 1:
                        plt.show(block=False)
                        plt.pause(2)
                    elif plot_results == 2:
                        plt.show()
                    plt.close()
                iraf_fwhms = irafsources['fwhm']  # Save FWHM in a list.
                iraf_fwhms = np.array(iraf_fwhms)  # Convert FWHM list to numpy array.
                # Print information about the file                                                                              #####
                # Calculate statistics of the FWHMs given by the loop over each source.                                         #####
                iraf_sdom = iraf_fwhms.std() / sqrt(len(iraf_fwhms))  # Standard deviation of the mean. Not used.
                num_IRAF_sources = len(irafsources)  # Number of stars IRAF found.
                print(f"No. of IRAF sources: {num_IRAF_sources}")  # Print number of IRAF stars found.
                iraf_fwhm = iraf_fwhms.mean()  # Calculate IRAF FWHM mean.
                iraf_std = iraf_fwhms.std()  # Calculate IRAF standard deviation.
                print(f"IRAF Calculated FWHM (pixels): {iraf_fwhm:.3f} +/- {iraf_std:.3f}")  # Print IRAF FWHM.
                # Calculate the flux of each star using PSF photometry. This uses the star positions calculated by
                # IRAFStarFinder earlier.
                daogroup = DAOGroup(2 * iraf_fwhm)  # Groups overlapping stars together.
                # sigma_clip = SigmaClip()
                # mmm_bkg = MMMBackground(sigma_clip=sigma_clip)
                # bkg_value = mmm_bkg(imgdata)
                psf_model = IntegratedGaussianPRF(
                    sigma=iraf_fwhm * gaussian_fwhm_to_sigma)  # Defime the PSF model to be used for photometry.
                psf_model.x_0.fixed = True  # Don't change the initial 'guess' of the star x positions to be provided.
                psf_model.y_0.fixed = True  # Don't change the initial 'guess' of the star y positions to be provided.
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
                # Perform the photometry on the background subtracted image. Also pass the fixed x-y positions and the
                # initial guess for the flux.
                result_tab = photometry(image=imgdata - median_val, init_guesses=pos)
                fluxes = result_tab['flux_fit']  # Store the fluxes as a list.
                fluxes = np.array(
                    fluxes) * u.ct  # Convert the fluxes to a numpy array and add the unit of count to it.
                fluxes = fluxes / exptime  # Normalize the fluxes by exposure time (unit is now counts / second)
                flux_uncs = result_tab['flux_unc']
                flux_uncs = np.array(flux_uncs) * u.ct
                flux_uncs = flux_uncs / exptime
                snr = (fluxes / flux_uncs).value
                instr_mags_sigma = 1.0857 / np.sqrt(snr)
                instr_mags_units = u.Magnitude(fluxes)  # Convert the fluxes to an instrumental magnitude.
                instr_mags = instr_mags_units.value  # Store the magnitudes without the unit attached.
                try:
                    if hdr['FILTER'] == 'B':
                        instr_mags_sigma += b_zpoint_std
                        instr_mags += b_zpoint
                    elif hdr['FILTER'] == 'V' or hdr['FILTER'] == 'G':
                        instr_mags_sigma += g_zpoint_std
                        instr_mags += g_zpoint
                    elif hdr['FILTER'] == 'R':
                        instr_mags_sigma += r_zpoint_std
                        instr_mags += r_zpoint
                except KeyError:
                    pass
                # Calculate the FWHM in units of arcseconds as opposed to pixels.
                try:
                    focal_length = hdr['FOCALLEN'] * u.mm  # Store the telescope's focal length with unit millimetres.
                    xpixsz = hdr['XPIXSZ']  # Store the size of the x pixels.
                    ypixsz = hdr['XPIXSZ']  # Store the size of the y pixels.
                    if xpixsz == ypixsz:  # If the pixels are square.
                        pixsz = xpixsz * u.um  # Store the pixel size with unit micrometre.
                        # Can find FOV by finding deg/pix and then multiplying by the x and y number of pix (NAXIS).
                        rad_per_pix = atan(
                            pixsz / focal_length) * u.rad  # Calculate the angular resolution of each pixel. Store with unit radians.
                        arcsec_per_pix = rad_per_pix.to(
                            u.arcsec)  # Convert the per pixel angular resultion to arcseconds.
                        iraf_FWHM_arcsec = iraf_fwhm * arcsec_per_pix.value  # Convert the IRAFStarFinder FWHM from pixels to arcsec.
                        iraf_std_arcsec = iraf_std * arcsec_per_pix  # Convert the IRAFStarFinder FWHM standard deviation from pixels to arcsec.
                        iraf_FWHMs_arcsec = iraf_fwhms * arcsec_per_pix.value
                        print(
                            f"IRAF Calculated FWHM (arcsec): {iraf_FWHM_arcsec:.3f} +/- {iraf_std_arcsec:.3f}")  # Print the IRAFStarFinder FWHM in arcsec.
                except KeyError:
                    iraf_FWHMs_arcsec = iraf_fwhms
                # print(irafsources['peak'] + median_val)                                                                         # Akin to 'max_pixel' from Shane's spreadsheet.
                # print(result_tab['x_0', 'y_0', 'flux_fit', 'flux_unc'])                                                         # Print the fluxes and their uncertainty for the current image.
                for obj_index, obj in enumerate(irafsources):
                    obj_x = obj['xcentroid']
                    obj_y = obj['ycentroid']
                    for sat_num, sat in enumerate(sat_locs, start=2):
                        sat_x = sat[0]
                        sat_y = sat[1]
                        if abs(sat_x - obj_x) < max_distance_from_sat and abs(
                                sat_y - obj_y) < max_distance_from_sat:
                            sats_table[filenum - reversing_index][sat_num] = instr_mags[obj_index]
                            uncertainty_table[filenum - reversing_index][sat_num] = instr_mags_sigma[obj_index]
                            sat_fwhm_table[filenum - reversing_index][sat_num] = iraf_FWHMs_arcsec[obj_index]
                print(sats_table[filenum - reversing_index])
            num_nans[sat_checked_mask] = 0
        change_sat_positions = False
    exptime = hdr['EXPTIME'] * u.s                                                                                  # Store the exposure time with unit seconds.
    mean_val, median_val, std_val = sigma_clipped_stats(imgdata)                                                    # Calculate background stats.
    iraffind = IRAFStarFinder(threshold=4*std_val, fwhm=2)                                               # Find stars using IRAF.
    irafsources = iraffind(imgdata - median_val)                                                                    # Subtract background median value.
    try:
        irafsources.sort('flux', reverse=True)                                                                      # Sort the stars by flux, from greatest to least.
    except Exception as e:
        # Reset the NaN boolean as it doesn't count.
        num_nans[:] = 0
        # change_sat_positions = True
        continue
    irafpositions = np.transpose((irafsources['xcentroid'], irafsources['ycentroid']))                              # Store source positions as a numpy array.
    if plot_results != 0:
        fig, ax = plt.subplots()
        ax.imshow(imgdata, cmap='gray', norm=LogNorm(), interpolation='nearest')
        ax.scatter(irafsources['xcentroid'], irafsources['ycentroid'],
                    s=100, edgecolor='red', facecolor='none')
        for i in range(0, num_sats):
            rect = patches.Rectangle((sat_locs[i,0] - max_distance_from_sat, sat_locs[i,1] - max_distance_from_sat),
                                     width=max_distance_from_sat * 2, height=max_distance_from_sat * 2,
                                     edgecolor='green', facecolor='none')
            ax.add_patch(rect)
        plt.title(os.path.basename(entry.path))
        if plot_results == 1:
            plt.show(block=False)
            plt.pause(2)
        elif plot_results == 2:
            plt.show()
        plt.close()
    iraf_fwhms = irafsources['fwhm']                                                                                # Save FWHM in a list.
    iraf_fwhms = np.array(iraf_fwhms)                                                                               # Convert FWHM list to numpy array.
    # Print information about the file                                                                              #####
    # Calculate statistics of the FWHMs given by the loop over each source.                                         #####
    iraf_sdom = iraf_fwhms.std() / sqrt(len(iraf_fwhms))                                                            # Standard deviation of the mean. Not used.
    num_IRAF_sources = len(irafsources)                                                                             # Number of stars IRAF found.
    print(f"No. of IRAF sources: {num_IRAF_sources}")                                                               # Print number of IRAF stars found.
    iraf_fwhm = iraf_fwhms.mean()                                                                                   # Calculate IRAF FWHM mean.
    iraf_std = iraf_fwhms.std()                                                                                     # Calculate IRAF standard deviation.
    print(f"IRAF Calculated FWHM (pixels): {iraf_fwhm:.3f} +/- {iraf_std:.3f}")                                     # Print IRAF FWHM.
    # Calculate the flux of each star using PSF photometry. This uses the star positions calculated by
    # IRAFStarFinder earlier.
    daogroup = DAOGroup(2*iraf_fwhm)                                                                                # Groups overlapping stars together.
    # sigma_clip = SigmaClip()
    # mmm_bkg = MMMBackground(sigma_clip=sigma_clip)
    # bkg_value = mmm_bkg(imgdata)
    psf_model = IntegratedGaussianPRF(sigma=iraf_fwhm*gaussian_fwhm_to_sigma)                                       # Defime the PSF model to be used for photometry.
    psf_model.x_0.fixed = True                                                                                      # Don't change the initial 'guess' of the star x positions to be provided.
    psf_model.y_0.fixed = True                                                                                      # Don't change the initial 'guess' of the star y positions to be provided.
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
    # Perform the photometry on the background subtracted image. Also pass the fixed x-y positions and the
    # initial guess for the flux.
    result_tab = photometry(image=imgdata - median_val, init_guesses=pos)
    fluxes = result_tab['flux_fit']                                                                                 # Store the fluxes as a list.
    fluxes = np.array(fluxes) * u.ct                                                                                # Convert the fluxes to a numpy array and add the unit of count to it.
    fluxes = fluxes / exptime                                                                                       # Normalize the fluxes by exposure time (unit is now counts / second)
    flux_uncs = result_tab['flux_unc']
    flux_uncs = np.array(flux_uncs) * u.ct
    flux_uncs = flux_uncs / exptime
    snr = (fluxes / flux_uncs).value
    instr_mags_sigma = 1.0857 / np.sqrt(snr)
    instr_mags_units = u.Magnitude(fluxes)                                                                          # Convert the fluxes to an instrumental magnitude.
    instr_mags = instr_mags_units.value                                                                             # Store the magnitudes without the unit attached.
    try:
        if hdr['FILTER'] == 'B':
            instr_mags_sigma += b_zpoint_std
            instr_mags += b_zpoint
        elif hdr['FILTER'] == 'V' or hdr['FILTER'] == 'G':
            instr_mags_sigma += g_zpoint_std
            instr_mags += g_zpoint
        elif hdr['FILTER'] == 'R':
            instr_mags_sigma += r_zpoint_std
            instr_mags += r_zpoint
    except KeyError:
        pass
    # Calculate the FWHM in units of arcseconds as opposed to pixels.
    try:
        focal_length = hdr['FOCALLEN'] * u.mm                                                                           # Store the telescope's focal length with unit millimetres.
        xpixsz = hdr['XPIXSZ']                                                                                          # Store the size of the x pixels.
        ypixsz = hdr['XPIXSZ']                                                                                          # Store the size of the y pixels.
        if xpixsz == ypixsz:                                                                                            # If the pixels are square.
            pixsz = xpixsz * u.um                                                                                       # Store the pixel size with unit micrometre.
            # Can find FOV by finding deg/pix and then multiplying by the x and y number of pix (NAXIS).
            rad_per_pix = atan(pixsz / focal_length) * u.rad                                                            # Calculate the angular resolution of each pixel. Store with unit radians.
            arcsec_per_pix = rad_per_pix.to(u.arcsec)                                                                   # Convert the per pixel angular resultion to arcseconds.
            iraf_FWHM_arcsec = iraf_fwhm * arcsec_per_pix.value                                                         # Convert the IRAFStarFinder FWHM from pixels to arcsec.
            iraf_std_arcsec = iraf_std * arcsec_per_pix                                                                 # Convert the IRAFStarFinder FWHM standard deviation from pixels to arcsec.
            iraf_FWHMs_arcsec = iraf_fwhms * arcsec_per_pix.value
            print(f"IRAF Calculated FWHM (arcsec): {iraf_FWHM_arcsec:.3f} +/- {iraf_std_arcsec:.3f}")                   # Print the IRAFStarFinder FWHM in arcsec.
    except KeyError:
        iraf_FWHMs_arcsec = iraf_fwhms
    # print(irafsources['peak'] + median_val)                                                                         # Akin to 'max_pixel' from Shane's spreadsheet.
    # print(result_tab['x_0', 'y_0', 'flux_fit', 'flux_unc'])                                                         # Print the fluxes and their uncertainty for the current image.
    for obj_index, obj in enumerate(irafsources):
        obj_x = obj['xcentroid']
        obj_y = obj['ycentroid']
        for sat_num, sat in enumerate(sat_locs, start=2):
            sat_x = sat[0]
            sat_y = sat[1]
            if abs(sat_x - obj_x) < max_distance_from_sat and abs(sat_y - obj_y) < max_distance_from_sat:
                sats_table[filenum][sat_num] = instr_mags[obj_index]
                uncertainty_table[filenum][sat_num] = instr_mags_sigma[obj_index]
                sat_fwhm_table[filenum][sat_num] = iraf_FWHMs_arcsec[obj_index]
                sat[0] = obj_x
                sat[1] = obj_y
    print(sats_table[filenum])
    sat_mags = np.array(list(sats_table[filenum]))
    mask = np.isnan(sat_mags[2:].astype(float))
    num_nans[mask] += 1
    num_nans[~mask] = 0
    print(num_nans)
    num_nan = max(num_nans)
    # print(num_nan)
    # if num_nan >= max_num_nan:
    if num_nan != 0 and (num_nan % max_num_nan) == 0:
        change_sat_positions = True
rmtree(temp_dir)
sats_table.pprint_all()
# plt.plot(sats_table['Time (JD)'], sats_table['NEOSSat'], 'o')
# plt.show(block=True)
# plt.close()
unique_filters = unique(sats_table, keys='Filter')
for filter_ in unique_filters['Filter']:
    mask = sats_table['Filter'] == filter_
    if filter_ == 'B':
        b_sats_table = sats_table[mask]
        b_uncertainty_table = uncertainty_table[mask]
        b_fwhm_table = sat_fwhm_table[mask]
    elif filter_ == 'G':
        g_sats_table = sats_table[mask]
        g_uncertainty_table = uncertainty_table[mask]
        g_fwhm_table = sat_fwhm_table[mask]
    elif filter_ == 'R':
        r_sats_table = sats_table[mask]
        r_uncertainty_table = uncertainty_table[mask]
        r_fwhm_table = sat_fwhm_table[mask]
# ascii.write(b_sats_table, 'C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/b_instr_mags.csv', format='csv', overwrite=True)
# ascii.write(g_sats_table, f'C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Light curves/{date_string} g_instr_mags.csv', format='csv', overwrite=True)
# ascii.write(r_sats_table, 'C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/r_instr_mags.csv', format='csv', overwrite=True)
# ascii.write(b_uncertainty_table, 'C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/b_uncertainty.csv', format='csv', overwrite=True)
# ascii.write(g_uncertainty_table, f'C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Light curves/{date_string} g_uncertainty.csv', format='csv', overwrite=True)
# ascii.write(r_uncertainty_table, 'C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/r_uncertainty.csv', format='csv', overwrite=True)
# ascii.write(b_fwhm_table, 'C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/b_fwhm.csv', format='csv', overwrite=True)
# ascii.write(g_fwhm_table, 'C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/g_fwhm.csv', format='csv', overwrite=True)
# ascii.write(r_fwhm_table, 'C:/Users/jmwawrow/Documents/DRDC_Code/FITS Tutorial/CSV files/Mar 20 Light Curve/r_fwhm.csv', format='csv', overwrite=True)
avg_sat_fwhm = []
avg_sat_fwhm_std = []
for row in g_fwhm_table:
    sat_full_row_numpy = np.array(list(row))
    sat_row_numpy = sat_full_row_numpy[2:].astype(float)
    avg_sat_fwhm.append(np.nanmean(sat_row_numpy))
    avg_sat_fwhm_std.append(np.nanstd(sat_row_numpy))
times_list = np.array(sats_table['Time (JD)'])
times_obj = Time(times_list, format='jd', scale='utc')
times_datetime = times_obj.to_value('datetime')
g_times_list = np.array(g_fwhm_table['Time (JD)'])
g_times_obj = Time(g_times_list, format='jd', scale='utc')
g_times_datetime = g_times_obj.to_value('datetime')
fig, ax = plt.subplots()
plt.errorbar(g_times_datetime, avg_sat_fwhm, yerr=avg_sat_fwhm_std, fmt='o', markersize=3, capsize=2)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.title(f'G Band Seeing - {date_string}')
plt.ylabel('FWHM (arcsec)')
plt.xlabel('Time (UTC)')
plt.show(block=True)
# plt.pause(3)
plt.close()
date_string = '25 Apr'
for sat in sat_names:
    b_interpolated = np.interp(times_list, b_sats_table['Time (JD)'][~np.isnan(b_sats_table[sat])],
                               b_sats_table[sat][~np.isnan(b_sats_table[sat])])
    g_interpolated = np.interp(times_list, g_sats_table['Time (JD)'][~np.isnan(g_sats_table[sat])],
                               g_sats_table[sat][~np.isnan(g_sats_table[sat])])
    r_interpolated = np.interp(times_list, r_sats_table['Time (JD)'][~np.isnan(r_sats_table[sat])],
                               r_sats_table[sat][~np.isnan(r_sats_table[sat])])
    b_interpolated[np.isnan(sats_table[sat])] = np.nan
    g_interpolated[np.isnan(sats_table[sat])] = np.nan
    r_interpolated[np.isnan(sats_table[sat])] = np.nan
    g_regular = np.full(len(g_interpolated), np.nan)
    g_uncertainty = np.full(len(g_interpolated), np.nan)
    mask = np.in1d(sats_table['Time (JD)'], g_sats_table['Time (JD)'])
    # print(mask)
    g_regular[mask] = g_sats_table[sat]
    g_uncertainty[mask] = g_uncertainty_table[sat]

    fig, ax = plt.subplots()
    ax.plot(times_datetime, g_interpolated, 'ko', markersize=3)
    ax.errorbar(times_datetime, g_regular, yerr=g_uncertainty, fmt='ko', markersize=3, capsize=1, label='g')
    ax.set_ylabel('Magnitude')
    ax.set_ylim([min(g_interpolated)*1.05, max(g_interpolated)*0.75])
    # ax.set_ylim(-12.2, -3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.gca().invert_yaxis()
    ax2 = ax.twinx()
    ax2.plot(times_datetime, b_interpolated - g_interpolated, 'bo', label='b-g', markersize=3)
    ax2.plot(times_datetime, b_interpolated - r_interpolated, 'go', label='b-r', markersize=3)
    ax2.plot(times_datetime, g_interpolated - r_interpolated, 'ro', label='g-r', markersize=3)
    ax2.set_ylabel('Colour Index')
    ax2.set_ylim([-6.5, 3.75])
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.gca().invert_yaxis()
    ax.set_xlabel("Time (UTC)")
    fig.legend()
    plt.title(f'{sat} Light Curve - {date_string}')
    plt.show(block=True)
    # plt.pause(5)
    plt.close()


# times_obj = Time(b_sats_table['Time (JD)'], format='jd', scale='utc')
# fig, ax = plt.subplots()
# for sat_num, sat in enumerate(sat_names, start=1):
#     ax.plot(times_obj.to_value('datetime'), b_sats_table[sat], 'o', label=sat)
# plt.ylabel("Instrumental Magnitude (b)")
# plt.gca().invert_yaxis()
# plt.xlabel("Time (UTC)")
# plt.title("B Band Light Curve - 20 Mar")
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# plt.legend()
# plt.show()
# # plt.show(block=False)
# # plt.pause(3)
# plt.close()

# times_obj = Time(g_sats_table['Time (JD)'], format='jd', scale='utc')
# fig, ax = plt.subplots()
# for sat_num, sat in enumerate(sat_names, start=1):
#     ax.plot(times_obj.to_value('datetime'), g_sats_table[sat], 'o', label=sat)
# plt.ylabel("Instrumental Magnitude (g)")
# plt.gca().invert_yaxis()
# plt.xlabel("Time (UTC)")
# plt.title("G Band Light Curve - 20 Mar")
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# plt.legend()
# plt.show()
# # plt.show(block=False)
# # plt.pause(3)
# plt.close()

# times_obj = Time(r_sats_table['Time (JD)'], format='jd', scale='utc')
# fig, ax = plt.subplots()
# for sat_num, sat in enumerate(sat_names, start=1):
#     ax.plot(times_obj.to_value('datetime'), r_sats_table[sat], 'o', label=sat)
# plt.ylabel("Instrumental Magnitude (r)")
# plt.gca().invert_yaxis()
# plt.xlabel("Time (UTC)")
# plt.title("R Band Light Curve - 20 Mar")
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# plt.legend()
# plt.show()
# # plt.show(block=False)
# # plt.pause(3)
# plt.close()
