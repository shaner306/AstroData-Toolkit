'''
This script was built to sample the flux of matched stars and track the flux over an observation period

Steps:
1. Use Pinpoint to match the images with catalgoue stars
2. Calculate Flux of matched stars
3. Compare flux over time
4. Query SIMBAD for Star Name

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

                    #focal_length = header['FOCALLEN']
                    #x_arcsecperpixel, y_arcsecperpixel = \
                    #    astro.calc_ArcsecPerPixel(XPIXSZ,YPIXSZ,focal_length)
                    x_arcsecperpixel=2.5
                    y_arcsecperpixel=2.5

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



                except Exception as e:
                    print(e)
                    continue

    return matched_star_dictionary


##
def statistics_of_matched_stars(matched_star_dictionary, save_loc,
                                 inner_rad, outer_rad,std_pass_threshold=0.1):
    '''
    Calculates the Mean and Standard Deviations of the matched stars.


    Parameters
    ----------
    matched_star_dictionary: the stars matched to the image in the previous
    step

    save_loc: str
        Respresents the save location of the pass/fail file
    std_pass_threshold: float
        The Stadard deviation threshold. If data has greater value then the
        threshold* mean then the std_pass_threshold then the star is
        considered to be too variable and thus does not get passed on to the
        passed_matched_stars dict and a fail is written beside their name in the outputted
        statistic text file.



    Returns
    -------
    pass_matched_stars: dict
        stars that pass the threshold test (std < mean* threshold_std)

    Outputs
    -------


    '''

    passed_matched_stars={} # Create dict for passed matched stars
    if os.path.exists(save_loc):
        with open(save_loc,'r+') as f:
            f.truncate(0) # clear prexisitng file if one exists
    with open(save_loc, 'a') as file:

        # Write Header
        file.write(f"Matched Stars: {len(matched_star_dictionary)},"
                   f"Pass Threshold: {std_pass_threshold} ,"
                   f"Inner Radius (arcsec): {inner_rad},"
                   f"Outer Radius(arcsecs): {outer_rad}\n")

        #Write Matched Stars and Their statistics (mean +/- 1 std)
        for i,matched_star in enumerate(matched_star_dictionary):
            array_of_match_star_items=[matched_star_item.value for matched_star_item in matched_star_dictionary[matched_star]]
            standard_dev=np.std(array_of_match_star_items)
            mean=np.mean(array_of_match_star_items)



            if standard_dev > np.mean(array_of_match_star_items)*std_pass_threshold:
                pass_fail_str='FAIL'
            else:
                passed_matched_stars[matched_star]=matched_star_dictionary[
                    matched_star]
                pass_fail_str='PASS'


            file.write(f"{matched_star:18} : {mean:12.2f} +/-"
                       f" {standard_dev:12.2f} --"
                       f" {pass_fail_str} \n")
    return passed_matched_stars

##
def plot_measure_flux(matched_star_dictionary,save_loc,
                      passed_matched_star_dictionary,passed_stars_sav_loc):
    ''''
    '''
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
    for i, matched_passed_star in enumerate(passed_matched_star_dictionary):
        array_of_passed_match_star_items = [matched_passed_star_time.value for
                                     matched_passed_star_time in
                                     matched_star_dictionary[matched_passed_star]]
        y_data = array_of_passed_match_star_items
        x_data = [*range(0, len(y_data))]
        plt.plot(x_data, y_data, '--')

        if i == 30: # Only print 30 fluxes
            break

    plt.title(
        ' Passed Stars Flux Variations/Frame in the Matched Stars as Detected '
        'by Pinpoint')
    plt.xlabel('Frame')
    plt.ylabel('Flux (in counts)')
    plt.savefig(passed_stars_sav_loc)
    plt.clf()
    plt.close()



#%% Query the Stars in SIMBAD
##
from astroquery.simbad import Simbad
def query_stars(passed_matched_stars):
    passed_matched_stars_list=passed_matched_stars.keys()
    repacked_passed_matched_star_list = [
        (passed_matched_stars_list_item.split('-')[
             0] + ' '
                  '' + passed_matched_stars_list_item.split('-')[1] + '-' +
         passed_matched_stars_list_item.split('-')[2]) for
        passed_matched_stars_list_item in passed_matched_stars_list]
    result_table=Simbad.query_objectids(repacked_passed_matched_star_list)







####################### Testing the Functions ################################

## Singular Request

inner_rad=32
outer_rad=48

image_dir=r"C:\Users\mstew\Documents\School and Work\Summer 2022\ImageProcessor\Landolt Fields\images\SA 26"
catalogue_dir=r"C:\Users\mstew\Documents\School and Work\Winter 2022\Work\StarCatalogues\USNO UCAC4"
matched_star_dictionary=get_matched_stars(image_dir, catalogue_dir, 18 ,2, \
                        0.8, 0, 60, 11, False, 32, 48)
passed_matched_stars=statistics_of_matched_stars(matched_star_dictionary,
                            r"C:\Users\mstew\Documents\School and "
                            r"Work\Summer 2022\ImageProcessor\Landolt Fields\images\SA 26\flux_calculations.txt",
                            inner_rad,outer_rad)
plot_measure_flux(matched_star_dictionary,r"C:\Users\mstew\Documents\School "
                                          r"and Work\Summer 2022\ImageProcessor\Landolt Fields\images\SA 26\fluxplots.png",
                  passed_matched_stars,r"C:\Users\mstew\Documents\School and Work\Summer 2022\ImageProcessor\Landolt Fields\images\SA 26\passed_stars_fluxplots.png"

                                        )

# Reset the Variables when complete
matched_star_dictionary={}

#%%
## Mutliple Requests

image_path = r"C:\Users\mstew\Documents\School and Work\Summer 2022\ImageProcessor\Landolt Fields\images\SA 26"
catalog_dir = r"C:\Users\mstew\Documents\School and Work\Winter 2022\Work\StarCatalogues\USNO UCAC4"
refstar_dir = r"C:\Users\mstew\Documents\GitHub\Astro2\Reference Star Files\Reference_stars_2022_02_17_d.txt"



inner_rad=32
outer_rad=48
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

    matched_star_dictionary_result=get_matched_stars(dirs, catalog_dir, 15 ,
                                                     2, 0.8, 0, 60, 11,
                                                     False, inner_rad, outer_rad)

    passed_matched_stars=statistics_of_matched_stars(
        matched_star_dictionary_result,
                                dirs+r"\flux_calculations.txt",inner_rad,
        outer_rad)
    plot_measure_flux(matched_star_dictionary_result,
                     dirs+r"\allfluxplots.png",passed_matched_stars,
                      dirs+r"\passed_stars_fluxplots.png")

    matched_star_dictionary_result={}

# for dirpath,dirname,filename in os.walk()
# image_dir=r"D:\School\Work - Winter 2022\Work\2022-03-16\2022-03-16"
# catalogue_dir=r"D:\School\StarCatalogues\USNO UCAC4"
# matched_star_dictionary=get_matched_stars(image_dir, catalogue_dir, 13, 3, 0.8, 1.5, 60, 11, False, 12, 24)
#statistics_of_matched_stars(matched_star_dictionary,)



