import os
from warnings import warn

import matplotlib
import matplotlib.cm as cm
import matplotlib.dates as mdates
import numpy as np
from astropy.io import ascii
import astropy.units as u
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from math import atan
from tqdm import tqdm

import sys
from os.path import dirname
src_path = dirname(dirname(__file__))
sys.path.append(os.path.join(src_path, 'general_tools'))

import AstroFunctions as astro


def _sky_survey_calc(directory,
                     plot_results=False,
                     save_plots=False,
                     file_suffix=(".fits", ".fit", ".fts"),
                     exposure_key='EXPTIME',
                     gb_final_transforms=None,
                     lat_key='SITELAT',
                     lon_key='SITELONG',
                     elev_key='SITEELEV',
                     **kwargs):
    # TODO: Docstring.
    # TODO: Fix errors when save_plots = False.
    # large_table_columns = init_large_table_columns()
    star_aux_table_columns = astro.init_star_aux_table_columns()

    filecount = 0
    file_paths = []
    file_names = []
    for dirpth, _, files in os.walk(directory):
        for file in files:
            if file.endswith(file_suffix):

                file_paths.append(os.path.join(dirpth, file))
                file_names.append(file)
    #             filecount += 1
    # filenames = sorted(os.listdir(temp_dir))

    save_loc = kwargs.get('save_loc')
    if not os.path.exists(save_loc):
        os.mkdir(save_loc)

    # for dirpath, dirnames, filenames in os.walk(directory):
    #     for filename in tqdm(filenames):
    #         if filename.endswith(file_suffix):
    # TODO: Edit this line to properly read the checkpoint
    # [100:]), start=100):
    for file_num, filepath in enumerate(tqdm(file_paths)):
        # filepath = os.path.join(dirpath, filename)
        hdr, imgdata = astro.read_fits_file(filepath)
        exptime = hdr[exposure_key]
        bkg, bkg_std = astro.calculate_img_bkg(imgdata)
        # print(bkg)
        irafsources = astro.detecting_stars(imgdata, bkg=bkg, bkg_std=bkg_std)
        if not irafsources:
            warn("Could not detect sources.")
            fwhm = np.nan
            fwhm_std = np.nan
            fwhm_arcsec = np.nan
            fwhm_arcsec_std = np.nan
            num_sources = 0
            t = Time(hdr['DATE-OBS'], format='fits', scale='utc')
            time = t.jd
            # time = hdr['DATE-OBS']
            img_filter = hdr['FILTER']
            # background_sky_brightness = astro.calculate_background_sky_brightness(
            #     bkg, hdr, exptime, gb_final_transforms)
            wcs = WCS(hdr)
            pixel_scale_deg = proj_plane_pixel_scales(wcs) * u.deg
            pixel_scale_arcsec = pixel_scale_deg.to(u.arcsec)
            x_arcsec_per_pix = pixel_scale_arcsec[0]
            y_arcsec_per_pix = pixel_scale_arcsec[1]
            square_arcsec_per_pix = x_arcsec_per_pix * y_arcsec_per_pix
            background_sky_brightness = astro.calculate_bsb_TEMP(bkg, hdr, exptime, square_arcsec_per_pix, gb_final_transforms)
            background_sky_brightness_sigma = astro.calculate_BSB_sigma(
                bkg, bkg_std, exptime)
            # try:
            azimuth = hdr['CENTAZ']
            elevation = hdr['CENTALT']
            airmass = hdr['AIRMASS']
            # except KeyError:
            #     continue
            star_aux_table_columns = astro.update_star_aux_columns(star_aux_table_columns,
                                                             file_names[file_num],
                                                             time,
                                                             img_filter,
                                                             fwhm,
                                                             fwhm_std,
                                                             fwhm_arcsec,
                                                             fwhm_arcsec_std,
                                                             num_sources,
                                                             background_sky_brightness,
                                                             background_sky_brightness_sigma,
                                                             azimuth,
                                                             elevation,
                                                             airmass)
            continue
        num_sources = len(irafsources)
        fwhms, fwhm, fwhm_std = astro.calculate_fwhm(irafsources)
        # photometry_result = perform_photometry.perform_PSF_photometry(irafsources, fwhm, imgdata, bkg=bkg)
        # fluxes = np.array(photometry_result['flux_fit'])
        # instr_mags = calculate_magnitudes(photometry_result, exptime)
        # instr_mags_sigma, snr = calculate_magnitudes_sigma(photometry_result, exptime)
        # try:
        wcs = WCS(hdr)
        skypositions = astro.convert_pixel_to_ra_dec(irafsources, wcs)
        try:
            altazpositions = astro.convert_ra_dec_to_alt_az(skypositions,
                                                        hdr,
                                                        lat_key=lat_key,
                                                        lon_key=lon_key,
                                                        elev_key=elev_key)
            azimuth = np.mean(altazpositions.az).degree
            elevation = np.mean(altazpositions.alt).degree
            airmass = np.mean(altazpositions.secz)
        except AttributeError as e:
            print(filepath)
            azimuth = hdr['CENTAZ']
            elevation = hdr['CENTALT']
            airmass = hdr['AIRMASS']
            # continue
        # except Exception:
        #     # TODO: change to an if/else statement.
        #     try:
        #         azimuth = hdr['CENTAZ']
        #         elevation = hdr['CENTALT']
        #         airmass = hdr['AIRMASS']
        #     except KeyError:
        #         continue
        wcs = WCS(hdr)
        pixel_scale_deg = proj_plane_pixel_scales(wcs) * u.deg
        pixel_scale_arcsec = pixel_scale_deg.to(u.arcsec)
        x_arcsec_per_pix = pixel_scale_arcsec[0]
        y_arcsec_per_pix = pixel_scale_arcsec[1]
        if x_arcsec_per_pix.value == 3600 and y_arcsec_per_pix.value == 3600:
            try:
                focal_length = hdr['FOCALLEN'] * u.mm
                xpixsz = hdr['XPIXSZ'] * u.um
                ypixsz = hdr['YPIXSZ'] * u.um
                x_rad_per_pix = atan(xpixsz / focal_length) * u.rad
                x_arcsec_per_pix = x_rad_per_pix.to(u.arcsec)
                y_rad_per_pix = atan(ypixsz / focal_length) * u.rad
                y_arcsec_per_pix = y_rad_per_pix.to(u.arcsec)
            except KeyError:
                print("Could not determine arcsec^2 / pix.")
                continue
        arcsec_per_pix = x_arcsec_per_pix
        square_arcsec_per_pix = x_arcsec_per_pix * y_arcsec_per_pix
        fwhms_arcsec, fwhm_arcsec, fwhm_arcsec_std = astro.convert_fwhm_to_arcsec_TEMP(
            arcsec_per_pix, fwhms, fwhm, fwhm_std)
        t = Time(hdr['DATE-OBS'], format='fits', scale='utc')
        time = t.jd
        # time = hdr['DATE-OBS']
        img_filter = hdr['FILTER']
        # background_sky_brightness = astro.calculate_background_sky_brightness(
        #     bkg, hdr, exptime, gb_final_transforms)
        background_sky_brightness = astro.calculate_bsb_TEMP(bkg, hdr, exptime, square_arcsec_per_pix, gb_final_transforms)
        background_sky_brightness_sigma = astro.calculate_BSB_sigma(
            bkg, bkg_std, exptime)
        # azimuth = hdr['CENTAZ']
        # elevation = hdr['CENTALT']
        # airmass = hdr['AIRMASS']
        star_aux_table_columns = astro.update_star_aux_columns(star_aux_table_columns,
                                                         file_names[file_num],
                                                         time,
                                                         img_filter,
                                                         fwhm,
                                                         fwhm_std,
                                                         fwhm_arcsec,
                                                         fwhm_arcsec_std,
                                                         num_sources,
                                                         background_sky_brightness,
                                                         background_sky_brightness_sigma,
                                                         azimuth,
                                                         elevation,
                                                         airmass)
        # if (file_num % 10) == 0:
        #     star_aux_table = create_star_aux_table(star_aux_table_columns)
        #     ascii.write(star_aux_table, os.path.join(save_loc, 'auxiliary_table.csv'), format='csv')
        #     with open(os.path.join(save_loc, 'checkpoint.txt'), 'a') as f:
        #         f.write(str(file_num))
        #         f.write('\n')
        #         f.write(filepath)
        #         f.write('\n')
    star_aux_table = astro.create_star_aux_table(star_aux_table_columns)
    ascii.write(star_aux_table, os.path.join(
        save_loc, 'auxiliary_table.csv'), format='csv')
    with open(os.path.join(save_loc, 'NightlyStats.txt'), 'a') as f:
        f.write('Parameter')
        f.write('\t')
        f.write('Value')
        f.write('\n')
    # TODO: Make these weighted means
    min_bsb = max(star_aux_table['BSB'][star_aux_table['BSB'] > 5])
    max_bsb = min(star_aux_table['BSB'][star_aux_table['BSB'] > 5])
    mean_bsb = np.mean(star_aux_table['BSB'][star_aux_table['BSB'] > 5])
    std_bsb = np.std(star_aux_table['BSB'][star_aux_table['BSB'] > 5])

    min_fwhm_arcsec = min(
        star_aux_table['FWHM_arcsec'][~np.isnan(star_aux_table['FWHM_arcsec'])])
    max_fwhm_arcsec = max(
        star_aux_table['FWHM_arcsec'][~np.isnan(star_aux_table['FWHM_arcsec'])])
    mean_fwhm_arcsec = np.mean(
        star_aux_table['FWHM_arcsec'][~np.isnan(star_aux_table['FWHM_arcsec'])])
    std_fwhm_arcsec = np.std(
        star_aux_table['FWHM_arcsec'][~np.isnan(star_aux_table['FWHM_arcsec'])])

    min_fwhm_pixel = min(
        star_aux_table['FWHM_pixel'][~np.isnan(star_aux_table['FWHM_arcsec'])])
    max_fwhm_pixel = max(
        star_aux_table['FWHM_pixel'][~np.isnan(star_aux_table['FWHM_arcsec'])])
    mean_fwhm_pixel = np.mean(
        star_aux_table['FWHM_pixel'][~np.isnan(star_aux_table['FWHM_arcsec'])])
    std_fwhm_pixel = np.std(
        star_aux_table['FWHM_pixel'][~np.isnan(star_aux_table['FWHM_arcsec'])])

    with open(os.path.join(save_loc, 'NightlyStats.txt'), 'a') as f:
        f.write(f'Minimum BSB:\t{min_bsb}')
        f.write('\n')
        f.write(f'Maximum BSB:\t{max_bsb}')
        f.write('\n')
        f.write(f'Mean BSB:\t{mean_bsb}')
        f.write('\n')
        f.write(f'Standard Deviation of BSB:\t{std_bsb}')
        f.write('\n')
        f.write(f'Minimum FWHM (arcsec):\t{min_fwhm_arcsec}')
        f.write('\n')
        f.write(f'Maximum FWHM (arcsec):\t{max_fwhm_arcsec}')
        f.write('\n')
        f.write(f'Mean FWHM (arcsec):\t{mean_fwhm_arcsec}')
        f.write('\n')
        f.write(f'Standard Deviation of FWHM (arcsec):\t{std_fwhm_arcsec}')
        f.write('\n')
        f.write(f'Minimum FWHM (pixels):\t{min_fwhm_pixel}')
        f.write('\n')
        f.write(f'Maximum FWHM (pixels):\t{max_fwhm_pixel}')
        f.write('\n')
        f.write(f'Mean FWHM (pixels):\t{mean_fwhm_pixel}')
        f.write('\n')
        f.write(f'Standard Deviation of FWHM (pixels):\t{std_fwhm_pixel}')
        f.write('\n')

    theta = star_aux_table['Azimuth'][star_aux_table['BSB'] > 5]
    r = star_aux_table['Elevation'][star_aux_table['BSB'] > 5]
    z = star_aux_table['BSB'][star_aux_table['BSB'] > 5]
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(7, 7))
    ax.set_theta_zero_location("N")
    norm = matplotlib.colors.Normalize(vmin=18, vmax=22)
    m = cm.ScalarMappable(cmap=plt.get_cmap('viridis_r'), norm=norm)
    m.set_array([])
    plt.colorbar(m)
    # Change contourf in the line below to scatter if you have only 1D theta,
    # r and brightness values
    ax.scatter(theta[~np.isnan(z)], r[~np.isnan(z)], c=z[~np.isnan(z)],
               cmap=plt.get_cmap('viridis_r'), norm=norm)
    rlabels = ax.get_ymajorticklabels()
    ax.set_rlim(bottom=90, top=15)
    for label in rlabels:
        label.set_color('black')
    plt.title('BSB (mag/arcsec$^2$)')
    plt.savefig(os.path.join(save_loc, 'BSB_polar.png'))
    plt.show()
    plt.close()

    times_list = np.array(
        star_aux_table['Time (JD)'][star_aux_table['BSB'] > 5])
    times_obj = Time(times_list, format='jd', scale='utc')
    times_datetime = times_obj.to_value('datetime')

    bsb = star_aux_table['BSB'][star_aux_table['BSB'] > 5]
    bsb_sigma = star_aux_table['BSB_sigma'][star_aux_table['BSB'] > 5]
    fig, ax = plt.subplots()
    _, _, bars = ax.errorbar(times_datetime,
                             bsb,
                             yerr=bsb_sigma,
                             fmt='o',
                             markersize=2,
                             capsize=0,
                             elinewidth=0.75)
    [bar.set_alpha(0.3) for bar in bars]
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.set_ylabel('BSB (mag/arcsec^2)')
    ax.set_xlabel('Time (UTC)')
    plt.title('BSB v. Time')
    plt.savefig(os.path.join(save_loc, 'BSB_v_time.png'))
    plt.show()
    plt.close()

    bsb = star_aux_table['BSB'][star_aux_table['BSB'] > 5]
    bsb_sigma = star_aux_table['BSB_sigma'][star_aux_table['BSB'] > 5]
    elevation = star_aux_table['Elevation'][star_aux_table['BSB'] > 5]
    fig, ax = plt.subplots()
    _, _, bars = ax.errorbar(bsb,
                             elevation,
                             xerr=bsb_sigma,
                             fmt='o',
                             markersize=2,
                             capsize=0,
                             elinewidth=0.75)
    [bar.set_alpha(0.3) for bar in bars]
    ax.set_xlabel('BSB (mag/arcsec$^2$)')
    ax.set_ylabel('Elevation')
    plt.title('BSB v. Elevation')
    plt.savefig(os.path.join(save_loc, 'BSB_v_elevation.png'))
    plt.show()
    plt.close()

    times_list = np.array(star_aux_table['Time (JD)'])
    times_obj = Time(times_list, format='jd', scale='utc')
    times_datetime = times_obj.to_value('datetime')

    fwhm_arcsec = star_aux_table['FWHM_arcsec']
    fwhm_arcsec_sigma = star_aux_table['FWHM_arcsec_sigma']
    fig, ax = plt.subplots()
    _, _, bars = ax.errorbar(times_datetime,
                             fwhm_arcsec,
                             yerr=fwhm_arcsec_sigma,
                             fmt='o',
                             markersize=2,
                             capsize=0,
                             elinewidth=0.75)
    [bar.set_alpha(0.3) for bar in bars]
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.set_ylabel('FWHM (arcsec)')
    ax.set_xlabel('Time (UTC)')
    plt.title('FWHM (arcsec) v. Time')
    plt.savefig(os.path.join(save_loc, 'FWHM_arcsec.png'))
    plt.show()
    plt.close()

    fwhm = star_aux_table['FWHM_pixel']
    fwhm_sigma = star_aux_table['FWHM_pixel_sigma']
    fig, ax = plt.subplots()
    _, _, bars = ax.errorbar(times_datetime,
                             fwhm,
                             yerr=fwhm_sigma,
                             fmt='o',
                             markersize=2,
                             capsize=0,
                             elinewidth=0.75)
    [bar.set_alpha(0.3) for bar in bars]
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.set_ylabel('FWHM (pixels)')
    ax.set_xlabel('Time (UTC)')
    plt.title('FWHM (pixels) v. Time')
    plt.savefig(os.path.join(save_loc, 'FWHM.png'))
    plt.show()
    plt.close()

    fwhm_arcsec = star_aux_table['FWHM_arcsec']
    fwhm_arcsec_sigma = star_aux_table['FWHM_arcsec_sigma']
    elevation = star_aux_table['Elevation']
    fig, ax = plt.subplots()
    _, _, bars = ax.errorbar(fwhm_arcsec,
                             elevation,
                             xerr=fwhm_arcsec_sigma,
                             fmt='o',
                             markersize=2,
                             capsize=0,
                             elinewidth=0.75)
    [bar.set_alpha(0.3) for bar in bars]
    ax.set_xlabel('FWHM (arcsec)')
    ax.set_ylabel('Elevation')
    plt.title('FWHM (arcsec) v. Elevation')
    plt.savefig(os.path.join(save_loc, 'FWHM_arcsec_v_elevation.png'))
    plt.show()
    plt.close()

    fwhm = star_aux_table['FWHM_pixel']
    fwhm_sigma = star_aux_table['FWHM_pixel_sigma']
    elevation = star_aux_table['Elevation']
    fig, ax = plt.subplots()
    _, _, bars = ax.errorbar(fwhm,
                             elevation,
                             xerr=fwhm_sigma,
                             fmt='o',
                             markersize=2,
                             capsize=0,
                             elinewidth=0.75)
    [bar.set_alpha(0.3) for bar in bars]
    ax.set_xlabel('FWHM (pixel)')
    ax.set_ylabel('Elevation')
    plt.title('FWHM (pixel) v. Elevation')
    plt.savefig(os.path.join(save_loc, 'FWHM_pixel_v_elevation.png'))
    plt.show()
    plt.close()

    theta = star_aux_table['Azimuth']
    r = star_aux_table['Elevation']
    z = star_aux_table['FWHM_arcsec']
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(7, 7))
    ax.set_theta_zero_location("N")
    norm = matplotlib.colors.Normalize(vmin=np.percentile(z[~np.isnan(z)], 7),
                                       vmax=max(z[~np.isnan(z)]))
    m = cm.ScalarMappable(cmap=plt.get_cmap('viridis_r'), norm=norm)
    m.set_array([])
    plt.colorbar(m)
    # Change contourf in the line below to scatter if you have only 1D theta,
    # r and brightness values
    ax.scatter(theta[~np.isnan(z)], r[~np.isnan(z)], c=z[~np.isnan(z)],
               cmap=plt.get_cmap('viridis_r'), norm=norm)
    rlabels = ax.get_ymajorticklabels()
    ax.set_rlim(bottom=90, top=15)
    for label in rlabels:
        label.set_color('black')
    plt.title('FWHM (arcsec)')
    plt.savefig(os.path.join(save_loc, 'FWHM_arcsec_polar.png'))
    plt.show()
    plt.close()

    theta = star_aux_table['Azimuth']
    r = star_aux_table['Elevation']
    z = star_aux_table['FWHM_pixel']
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(7, 7))
    ax.set_theta_zero_location("N")
    norm = matplotlib.colors.Normalize(vmin=np.percentile(z[~np.isnan(z)], 7),
                                       vmax=max(z[~np.isnan(z)]))
    m = cm.ScalarMappable(cmap=plt.get_cmap('viridis_r'), norm=norm)
    m.set_array([])
    plt.colorbar(m)
    # Change contourf in the line below to scatter if you have only 1D theta,
    # r and brightness values
    ax.scatter(theta[~np.isnan(z)], r[~np.isnan(z)], c=z[~np.isnan(z)],
               cmap=plt.get_cmap('viridis_r'), norm=norm)
    rlabels = ax.get_ymajorticklabels()
    ax.set_rlim(bottom=90, top=15)
    for label in rlabels:
        label.set_color('black')
    plt.title('FWHM (pixels)')
    plt.savefig(os.path.join(save_loc, 'FWHM_pixel_polar.png'))
    plt.show()
    plt.close()

    return star_aux_table

