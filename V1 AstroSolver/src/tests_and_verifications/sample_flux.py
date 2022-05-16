'''
This script was built to sample the flux of matched stars and track the flux over an observation period

Steps:
1. Use Pinpoint to match the images with catalgoue stars
2. Calculate Flux of matched stars
3. Compare flux over time

'''
##
from astropy.table import Table,QTable
import astropy.units as u
import os

import win32com
import win32com.client

import AstroFunctions as astro
import numpy as np
import matplotlib.pyplot as plt


##
def get_matched_stars(image_dir,catloc,max_mag,sigma, catexp, match_residual,max_solve_time,catalog,use_sextractor,inner_ap=12,outer_ap=24):
    '''
    A similar function to that of pinpoint_solve. This version is built specifically for sampling the flux of matched
    stars detected by Pinpoint

    Parameters
    ----------
    image_dir: str
        describes the path to the images
    catloc: str
        str describing the path to the catalogue
    max_mag:
    sigma
    match_reisudal
    max_solve_time
    catalog
    use_sextractor
    inner_ap:
        inner aperture in arcseconds
    outer_ap:
        outer aperture in arcseconds
    Returns
    -------
    flux_table: AstroPy.QTable
        QTable which describes the fluxes of all the matched stars in the images

    '''
    matched_star_dictionary={}
    file_suffix=(".fits",".fit",".fts")
    for dirpath,dirnames,filenames in os.walk(image_dir):
        for filename in filenames:
            if (filename.endswith(file_suffix)):
                filepath=os.path.join(dirpath,filename)
                f=win32com.client.Dispatch("Pinpoint.plate")

                "Import Data from FITS Image"
                header = 0
                imagehdularray, date, exposure_Time, imagesizeX, imagesizeY, \
                fitsdata, filt, header, XPIXSZ, YPIXSZ, wcs = \
                    astro.fits_header_import(filepath)


                """Pinpoint Solve"""
                try:
                    # Attaches FITS image to the created Plate
                    f.AttachFITS(filepath)
                    # We set the declination of the image to the target object's declination
                    f.Declination = f.targetDeclination
                    f.RightAscension = f.targetRightAscension

                    x_arcsecperpixel, y_arcsecperpixel = \
                        astro.calc_ArcsecPerPixel(header)
                    # yBin = 4.33562092816E-004*3600;
                    # xBin =  4.33131246330E-004*3600;
                    # f.ArcsecperPixelHoriz  = 4.556
                    # f.ArcsecperPixelVert =  4.556
                    # if space_based_bool==1:
                    #     yBin = 4.33562092816E-004*3600*2
                    # %Image Specific Pixel Size in arcsec /
                    # Obtained from FITS Header and converted from deg
                    #     xBin =\
                    # 4.33131246330E-004*3600*2#%Image Specific Pixel Size in arcsec /
                    # Obtained from FITS Header and converted from deg
                    # else:
                    #     yBin = 4.33562092816E-004*3600
                    # %Image Specific Pixel Size in arcsec /
                    # Obtained from FITS Header and converted from deg
                    #     xBin =  4.33131246330E-004*3600
                    # %Image Specific Pixel Size in arcsec /
                    # Obtained from FITS Header and converted from deg
                    # f.ArcsecperPixelHoriz  = xBin
                    # %CCD Pixel scale on CD
                    # f.ArcsecperPixelVert = yBin
                    f.ArcsecperPixelHoriz = x_arcsecperpixel
                    # %CCD Pixel scale on CD
                    f.ArcsecperPixelVert = y_arcsecperpixel

                    "Pinpoint Solve Inputs"
                    # TODO Add Inputs for pinpoint solving to GUI
                    f.Catalog = 11
                    f.CatalogPath = catloc
                    f.UseSExtractor = use_sextractor
                    f.CatalogMaximumMagnitude = max_mag
                    # print(f.CatalogMaximumMagnitude)
                    f.CatalogExpansion = catexp
                    f.MaxSolveTime = max_solve_time
                    f.MaxMatchResidual = match_residual
                    f.SigmaAboveMean = sigma
                    f.RemoveHotPixels()

                    "Pinpoint Solving"

                    # FIXME f.Solve intronsicly calls Find Catalog Stars and Image Stars
                    f.FindCatalogStars()
                    # print(f.CatalogStars.Count)
                    f.FindImageStars()
                    # print(f.ImageStars.Count)
                    f.Solve()
                    matched_stars=f.MatchedStars


                    iterative_exclusion_list=[] # only allows one instance of the star/frame
                    for i,star in enumerate(matched_stars):

                        try:

                            if star.MatchedCatalogStar.Identification in iterative_exclusion_list:
                                continue
                            else:

                                # FIXME: Greedy Algorithm


                                matched_star_dictionary[star.MatchedCatalogStar.Identification].append(
                                    f.MeasureFlux(star.X, star.Y, inner_ap, outer_ap) * u.count)
                                iterative_exclusion_list.append(star.MatchedCatalogStar.Identification)

                        except:
                            # First Instance of the star
                            matched_star_dictionary[star.MatchedCatalogStar.Identification] = ([
                                f.MeasureFlux(star.X, star.Y, inner_ap, outer_ap) * u.count])
                            iterative_exclusion_list.append(star.MatchedCatalogStar.Identification)
                    f=None
                except:
                    print('Could Not Identify ')
    return matched_star_dictionary


##
def statistics_of_matched_stars(matched_star_dictionary,save_loc):
    '''
    Calculates the Mean and Standard Deviations of the
    Parameters
    ----------
    flux_table
    save_loc: str
        Respresents the save location of the pass/fail file

    Returns
    -------

    '''
    if os.path.exists(save_loc):
        with open(save_loc,'r+') as f:
            f.truncate(0) # clear prexisitng file

    for i,matched_star in enumerate(matched_star_dictionary):
        array_of_match_star_items=[matched_star_item.value for matched_star_item in matched_star_dictionary[matched_star]]
        standard_dev=np.std(array_of_match_star_items)
        mean=np.mean(array_of_match_star_items)
        pass_fail_str='PASS'
        if standard_dev > 0.10*mean:
            pass_fail_str='FAIL'
        with open(save_loc,'a') as file:
            file.write(f"{matched_star} : {mean} +/- {standard_dev} -- {pass_fail_str} \n ")

##
def plot_measure_flux(matched_star_dictionary,save_loc):
    #maxcount=max(len(v) for v in matched_star_dictionary.values())
    #x_data=range()
    for i,matched_star in enumerate(matched_star_dictionary):

        array_of_match_star_items=[matched_star_item.value for matched_star_item in matched_star_dictionary[matched_star]]
        y_data=array_of_match_star_items
        x_data=[*range(0,len(y_data))]
        plt.plot(x_data,y_data,'--')

        if i==100:
            break

    plt.title('Flux Variations/Frame in the Matched Stars as Detected by Pinpoint')
    plt.xlabel('Frame')
    plt.ylabel('Flux (in counts)')
    plt.savefig(save_loc)



## Example of GC 279 B Filter
image_dir=r"D:\School\Work - Winter 2022\Work\2022-03-16\2022-03-16\GD 279\LIGHT\B"
catalogue_dir=r"D:\School\StarCatalogues\USNO UCAC4"
matched_star_dictionary=get_matched_stars(image_dir, catalogue_dir, 13, 3, 0.8, 1.5, 60, 11, False, 12, 24)
statistics_of_matched_stars(matched_star_dictionary,r"D:\School\Work - Winter 2022\Work\2022-03-16\2022-03-16\GD 279\LIGHT\B\flux_calculations.txt")
plot_measure_flux(matched_star_dictionary,r"D:\School\Work - Winter 2022\Work\2022-03-16\2022-03-16\GD 279\LIGHT\B\fluxplots.png")