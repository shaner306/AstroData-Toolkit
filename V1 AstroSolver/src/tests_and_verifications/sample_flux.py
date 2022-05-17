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
def get_matched_stars(image_dir: str, catloc: str, max_mag: float,
                      sigma: float, catexp: float,
                      match_residual: float,
                      max_solve_time: float,
                      catalog: int,
                      use_sextractor: bool,
                      inner_ap: object = 12,
                      outer_ap: object = 24) -> object:
    '''
    A similar function to that of pinpoint_solve. This version is built specifically for sampling the flux of matched
    stars detected by Pinpoint

    Parameters
    ----------
    image_dir: str
        describes the path to the images
    catloc: str
        str describing the path to the catalogue
    max_mag: float
        float denoting the maximum magntitude for selection of stars from
        the reference catalogue. Default is 20
    sigma
    match_reisudal: float
        max match residual. The maximum positional error (in arc-seconds)
        for matching image and catalog stars
    max_solve_time:float
        the maximum solve time Pinpoint will spend on a given image

    catalog: int
        catalog number deifining which catalog to use. Default is 11 (UCAC4)
        - UCAC4: 11
        - ATLAS: 12
        - USNOA2.0: 5
        - USNO_B: 7

    use_sextractor: bool
        option to call and use source extractor
    inner_ap:
        inner aperture in arcseconds
    outer_ap:
        outer aperture in arcseconds
    Returns
    -------
    flux_table: AstroPy.QTable
        QTable which describes the fluxes of all the matched stars in the images

    '''
    matched_star_dictionary = {}
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
                    f.DetachFITS()

                    f = None

                except:
                    print('Could Not Identify ')
    return matched_star_dictionary


##
def statistics_of_matched_stars(matched_star_dictionary,save_loc,
                                pass_threshold=0.1):
    '''
    Calculates the Mean and Standard Deviations of the matched stars.


    Parameters
    ----------
    flux_table
    save_loc: str
        Respresents the save location of the pass/fail file
    pass_threshold: float
        the percentage that is multiplied by the mean of the star's flux
        value to trigger the FAIL condition. FAIL conditions show that the
        flux deviates largely from it's average value



    Returns
    -------

    Outputs
    -------


    '''
    if os.path.exists(save_loc):
        with open(save_loc,'r+') as f:
            f.truncate(0) # clear prexisitng file
    with open(save_loc, 'a') as file:
        file.write(f"Matched Stars: {len(matched_star_dictionary)},"
                   f"Pass Threshold: {pass_threshold} * star mean  \n")
        for i,matched_star in enumerate(matched_star_dictionary):
            array_of_match_star_items=[matched_star_item.value for matched_star_item in matched_star_dictionary[matched_star]]
            standard_dev=np.std(array_of_match_star_items)
            mean=np.mean(array_of_match_star_items)
            pass_fail_str='PASS'
            if standard_dev > pass_threshold*mean:
                pass_fail_str='FAIL'


            file.write(f"{matched_star:18} : {mean:.2f} +/-"
                       f" {standard_dev:.2f} --"
                       f" {pass_fail_str} \n ")

##
def plot_measure_flux(matched_star_dictionary,save_loc):
    #maxcount=max(len(v) for v in matched_star_dictionary.values())
    #x_data=range()
    for i,matched_star in enumerate(matched_star_dictionary):

        array_of_match_star_items=[matched_star_item.value for matched_star_item in matched_star_dictionary[matched_star]]
        y_data=array_of_match_star_items
        x_data=[*range(0,len(y_data))]
        plt.plot(x_data,y_data,'--')

        if i==30:
            break

    plt.title('Flux Variations/Frame in the Matched Stars as Detected by Pinpoint')
    plt.xlabel('Frame')
    plt.ylabel('Flux (in counts)')
    plt.savefig(save_loc)
    plt.clf()
    plt.close()



## Singular Request

#
# image_dir=r"C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-10-26\Automated Pointing Runs 009-014\Automated Pointing Run 009"
# catalogue_dir=r"C:\Users\mstew\Documents\School and Work\Winter 2022\Work\StarCatalogues\USNO UCAC4"
# matched_star_dictionary=get_matched_stars(image_dir, catalogue_dir, 10 ,3, \
#                         0.8, 1.5, 60, 11, False, 12, 24)
# statistics_of_matched_stars(matched_star_dictionary,
#                             r"C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-10-26\Automated Pointing Runs 009-014\Automated Pointing Run 009\flux_calculations.txt")
# plot_measure_flux(matched_star_dictionary,r"C:\Users\mstew\Documents\School "
#                                           r"and Work\Winter "
#                                           r"2022\Work\2021-10-26\Automated Pointing Runs 009-014\Automated Pointing Run 009\fluxplots.png")
#
# # Reset the Variables when complete
# matched_star_dictionary={}

#%%
## Mutliple Requests
dataset_folder = r"C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-04-25\Siderial Stare Mode -Reduced"
catalog_dir = r"C:\Users\mstew\Documents\School and Work\Winter 2022\Work\StarCatalogues\USNO UCAC4"
refstar_dir = r"C:\Users\mstew\Documents\GitHub\Astro2\Reference Star Files\Reference_stars_2022_02_17_d.txt"

# Pinpoint Solve Parameters
max_mag = 13
sigma = 3
max_solve_time = 60  # Seconds
match_residual = 1.5
catalog = 11
catalog_exp = 0.8
use_sextractor = False
all_sky_solve = False
space_based_bool = False
photometry_method = 'aperture'
aperture_estimation_mode = 'mean'

image_path = dataset_folder

list_subfolders_with_paths = [folder.path for folder in os.scandir(image_path) if folder.is_dir()]

target_dirs=[]

#TODO : Clean up Code
for subfolder in list_subfolders_with_paths:
    try:
        for subfolder2 in (os.scandir(subfolder)):
            try:
                for subfolder3 in (os.scandir(subfolder2)):
                    if subfolder3.is_dir():
                        target_dirs.append(subfolder3.path)
                    else:
                        continue
            except:
                continue
    except:
        continue








for dirs in target_dirs:



    matched_star_dictionary_result=get_matched_stars(dirs,
                                            catalog_dir,
                                            max_mag,
                                            sigma,
                                            catalog_exp,
                                            match_residual,
                                            max_solve_time,
                                            catalog,
                                            space_based_bool,
                                            use_sextractor,
                                            )
    statistics_of_matched_stars(matched_star_dictionary_result,
                                dirs+"/flux_calculations.txt")
    plot_measure_flux(matched_star_dictionary_result,
                     dirs+"/fluxplots.png")

    matched_star_dictionary_result={}




# for dirpath,dirname,filename in os.walk()
# image_dir=r"D:\School\Work - Winter 2022\Work\2022-03-16\2022-03-16"
# catalogue_dir=r"D:\School\StarCatalogues\USNO UCAC4"
# matched_star_dictionary=get_matched_stars(image_dir, catalogue_dir, 13, 3, 0.8, 1.5, 60, 11, False, 12, 24)
#statistics_of_matched_stars(matched_star_dictionary,)



