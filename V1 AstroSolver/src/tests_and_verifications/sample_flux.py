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
class stars:
    '''
    list_of_stars: list of strings
    A list of star string names

    '''
    list_of_stars=[]

    def __init__(self,name):
        self.name=name
        self.list_of_stars.append(name)

        self.counts=[]
        self.xs=[]
        self.ys=[]
        self.ras=[]
        self.decs=[]


    def update(self,flux,x,y,ra,dec):
        '''

        Parameters
        ----------
        flux: value of the fluxes
        x
        y
        ra
        dec

        Returns
        -------

        '''
        self.counts.append(flux)
        self.xs.append(x)
        self.ys.append(y)
        self.ras.append(ra)
        self.decs.append(dec)





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
    matched_star_dictionary={}
    matched_star_collection={}
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
                    x_arcsecperpixel, y_arcsecperpixel = \
                        astro.calc_ArcsecPerPixel(header)
                    #x_arcsecperpixel=2.5
                    #y_arcsecperpixel=2.5


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
                            if star.MatchedCatalogStar.Identification in list(matched_star_dictionary.keys()):

                                if star.MatchedCatalogStar.Identification in iterative_exclusion_list:
                                    continue
                                else:
                                    # FIXME: Greedy Algorithm
                                    measured_flux = f.MeasureFlux(star.X, star.Y, inner_ap, outer_ap) * u.count
                                    #simple_matched_star_dictionary[star.MatchedCatalogStar.Identification].append(
                                    #    measured_flux)


                                    matched_star_dictionary[star.MatchedCatalogStar.Identification].append(
                                        measured_flux)
                                    iterative_exclusion_list.append(star.MatchedCatalogStar.Identification)

                                    # Extract matched class from dictionary object
                                    star_instance=[]
                                    for i in matched_star_collection:
                                        if i == star.MatchedCatalogStar.Identification:
                                            star_instance.append(matched_star_collection[i])
                                    star_instance=star_instance[0]
                                    star_instance.update(measured_flux,star.X,star.Y,f.RightAscension,f.Declination)
                            else:
                                # First Instance of the star
                                measured_flux = f.MeasureFlux(star.X, star.Y, inner_ap, outer_ap) * u.count
                                matched_star_dictionary[star.MatchedCatalogStar.Identification] = ([
                                    measured_flux])
                                iterative_exclusion_list.append(star.MatchedCatalogStar.Identification)
                                star_instance=stars(star.MatchedCatalogStar.Identification)
                                star_instance.update(measured_flux,star.X,star.Y,f.RightAscension,f.Declination)
                                matched_star_collection[star_instance.name]=star_instance # Add Star to Collection
                        except Exception as e:
                            raise (e)
                            break


                    f.DetachFITS()

                    f = None



                except Exception as e:
                    print(e)
                    continue

    return matched_star_dictionary,matched_star_collection


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



            if standard_dev > np.mean(array_of_match_star_items)*std_pass_threshold or standard_dev==0:
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
    for i,passed_matched_star in enumerate(repacked_passed_matched_star_list):
        try:
            if 'result_table' in locals():
                result_table.add_row(Simbad.query_object(passed_matched_star))
            else:
                result_table = Simbad.query_object(passed_matched_star)
        except Exception as e:
            print(e)
            continue










####################### Testing the Functions ################################

## Singular Request

inner_rad=32
outer_rad=48

image_dir=r"D:\School\Work - Winter 2022\Work\2021-03-21\HIP 2894\LIGHT\B"
catalogue_dir=r"D:\School\StarCatalogues\USNO UCAC4"
matched_star_dictionary,matched_star_collection=get_matched_stars(image_dir, catalogue_dir, 15 ,3, \
                        0.8, 2, 60, 11, False, 32, 48)
passed_matched_stars=statistics_of_matched_stars(matched_star_dictionary,
                            r"D:\School\Work - Winter 2022\Work\2021-03-21\HIP 2894\LIGHT\B\flux_calculations.txt",
                            inner_rad,outer_rad)
plot_measure_flux(matched_star_dictionary,r"D:\School\Work - Winter 2022\Work\2021-03-21\HIP 2894\LIGHT\B\fluxplots.png",
                  passed_matched_stars,r"D:\School\Work - Winter 2022\Work\2021-03-21\HIP 2894\LIGHT\B\passed_stars_fluxplots.png"

                                        )

# Reset the Variables when complete
matched_star_dictionary={}

#%%
## Mutliple Requests

image_path = r"D:\School\Work - Winter 2022\Work\2021-03-21"
catalog_dir = r"D:\School\StarCatalogues\USNO UCAC4"
refstar_dir = r"C:\Users\stewe\Documents\GitHub\Astro2\Reference Star Files\Reference_stars_2022_02_17_d.txt"



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

    matched_star_dictionary,matched_star_collection=get_matched_stars(dirs, catalog_dir, 15,
                                                     3, 0.8, 2, 60, 11,
                                                     False, inner_rad, outer_rad,)

    passed_matched_stars=statistics_of_matched_stars(
        matched_star_dictionary,
                                dirs+r"\flux_calculations.txt",inner_rad,
        outer_rad)
    plot_measure_flux(matched_star_dictionary,
                     dirs+r"\allfluxplots.png",passed_matched_stars,
                      dirs+r"\passed_stars_fluxplots.png")

    #matched_star_dictionary_result={}

# for dirpath,dirname,filename in os.walk()
# image_dir=r"D:\School\Work - Winter 2022\Work\2022-03-16\2022-03-16"
# catalogue_dir=r"D:\School\StarCatalogues\USNO UCAC4"
# matched_star_dictionary=get_matched_stars(image_dir, catalogue_dir, 13, 3, 0.8, 1.5, 60, 11, False, 12, 24)
#statistics_of_matched_stars(matched_star_dictionary,)


##
