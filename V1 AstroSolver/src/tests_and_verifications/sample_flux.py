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
    def set_pass_bool(self, pass_bool):
        self.pass_bool=pass_bool
    def add_other_ids(self,other_ids):
        self.other_ids=other_ids






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
                                    star_instance.update(measured_flux,
                                                         star.X,star.Y,
                                                         star.RightAscension,
                                                         star.Declination)
                            else:
                                # First Instance of the star
                                measured_flux = f.MeasureFlux(star.X, star.Y, inner_ap, outer_ap) * u.count
                                matched_star_dictionary[star.MatchedCatalogStar.Identification] = ([
                                    measured_flux])
                                iterative_exclusion_list.append(star.MatchedCatalogStar.Identification)
                                star_instance=stars(star.MatchedCatalogStar.Identification)
                                star_instance.update(measured_flux,star.X,
                                                     star.Y,
                                                     star.RightAscension,
                                                     star.Declination)
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
def statistics_of_matched_stars(matched_star_collection, save_loc,
                                 inner_rad, outer_rad,std_pass_threshold=0.1):
    '''
    Calculates the Mean and Standard Deviations of the matched stars.


    Parameters
    ----------
    matched_star_collection: the stars matched to the image in the previous
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
        file.write(f"Matched Stars: {len(matched_star_collection)},"
                   f"Pass Threshold: {std_pass_threshold} ,"
                   f"Inner Radius (arcsec): {inner_rad},"
                   f"Outer Radius(arcsecs): {outer_rad}\n")

        #Write Matched Stars and Their statistics (mean +/- 1 std)
        for i,matched_star in enumerate(matched_star_collection):

            array_of_match_star_items=[flux.value for flux in matched_star_collection[matched_star].counts]
            # standard_dev=np.std(array_of_match_star_items)
            # mean=np.mean(array_of_match_star_items)
            standard_dev=np.std(array_of_match_star_items)
            mean=np.mean(array_of_match_star_items)


            if standard_dev > np.mean(array_of_match_star_items)*std_pass_threshold or standard_dev==0:
                pass_fail_str='FAIL'
                matched_star_collection[matched_star].set_pass_bool(False)
            else:
                passed_matched_stars[matched_star]=matched_star_collection[
                    matched_star]
                pass_fail_str='PASS'
                matched_star_collection[matched_star].set_pass_bool(True)


            file.write(f"{matched_star:18} : {mean:12.2f} +/-"
                       f" {standard_dev:12.2f} --"
                       f" {pass_fail_str} \n")
    return passed_matched_stars

##
def plot_measure_flux(matched_star_dictionary, save_loc,
                      passed_matched_star_collection, passed_stars_sav_loc):
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
    for i, matched_passed_star in enumerate(passed_matched_star_collection):
        array_of_passed_match_star_items = [matched_passed_star_time.value for
                                     matched_passed_star_time in
                                     passed_matched_star_collection[
                                                matched_passed_star].counts]
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


def query_passed_matched_stars(passed_matched_stars):
    if 'ids' not in Simbad.get_votable_fields():
        Simbad.add_votable_fields('ids')
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
                if result_table != None:
                    result_table.add_row(Simbad.query_object(passed_matched_star))
                else:
                    result_table = Simbad.query_object(passed_matched_star)
            else:
                result_table = Simbad.query_object(passed_matched_star)
        except Exception as e:
            print(e)
            continue
    return result_table

## Find all IDS corresponding to the stars in the reference star file
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
def find_ids_for_ref_stars(ref_stars_file):
    '''
    Utilizing the power of astroquery with region look up, ids for reference stars are looked up



    Parameters
    ----------
    ref_stars_file: str
    string directory which the reference star data will be read from

    Returns
    -------
    resultant_table: Astropy.Table
    table which contains the values extracted from the reference star file
    and the various id's queried.
    '''

    if 'ids' not in Simbad.get_votable_fields():
        Simbad.add_votable_fields('ids')
    other_ids = []
    with open(ref_stars_file,'r') as file:
        reference_star_file_lines=file.readlines()
        formatted_ras=[]
        formated_decs=[]

        resultant_table = Table(names=((reference_star_file_lines[0]).split('\t'))[0:17],
                              dtype=['str','str','str','float64','float64',
                                     'float64','float64','float64','float64','float64',
                                     'float64','float64','float64','float64','float64',
                                     'float64', 'float64'
                                     ])
        for i,reference_star_file_line in enumerate(reference_star_file_lines):
            if i==0:
                continue # title
            else:

                reference_star_table_line=(reference_star_file_line.split('\t'))[0:17]
                resultant_table.add_row(reference_star_table_line)
                # FIXME: Add proper string spliting
                formatted_ra=f"{(reference_star_table_line[1])[0:3]}h" \
                             f" {(reference_star_table_line[1])[3:6]}m" \
                             f"{(reference_star_table_line[1])[6:14]}s"
                formatted_ras.append(formatted_ra)
                formatted_dec=f"{(reference_star_table_line[2])[0:3]}deg"\
                           f"{(reference_star_table_line[2])[3:6]}m"\
                            f"{reference_star_table_line[2][6:14]}s"
                formated_decs.append(formatted_dec)


                # Will need to query region instead of querying name due to
                # non-standardized naming scheme


                search_name=reference_star_table_line[0]

    # Vectorize then query
    query_result=Simbad.query_region(SkyCoord(formatted_ras,formated_decs))
    for i,reference_star_table_line in enumerate(resultant_table):
        search_name=reference_star_table_line[0]

        for j,query_line in enumerate(query_result):
            if query_line['SCRIPT_NUMBER_ID']==(i+1) :



                #FIXME: Greedy Algorithm
                other_ids.append(query_line['IDS'])
                break

    resultant_table.add_column(other_ids,name='IDS')




    return resultant_table,formatted_ras,formated_decs, other_ids
# TODO: Find Offline Version of this


#%%
##
def  find_ids_for_matched_stars(passed_matched_stars):
    average_ras=[]
    average_decs=[]
    if 'ids' not in Simbad.get_votable_fields():
        Simbad.add_votable_fields('ids')

    # extract all ras, decs average them out then query using averages
    for passed_matched_star in passed_matched_stars:
        average_ra=f"{np.mean(passed_matched_stars[passed_matched_star].ras)}h"
        average_ras.append(average_ra)
        average_dec=f"{np.mean(passed_matched_stars[passed_matched_star].decs)}deg"
        average_decs.append(average_dec)
    try:
        query_results=Simbad.query_region(SkyCoord(ra=average_ras,
                                                   dec=average_decs))
    except:
        individual_query=True
        print('Vectorized Query did not work, Trying individual Query')
    if 'individual_query' in locals():
        if individual_query==True:
            # perform individual queries
            for i in range(len(average_ras)):
                try:
                    if query_results in locals() and query_results != None:
                        query_results.add_row(Simbad.query_region(SkyCoord(
                            ra=average_ras[i], dec=average_decs[i])))
                    else:
                        query_results=Simbad.query_region(SkyCoord(
                            ra=average_ras[i], dec=average_decs[i]))
                except Exception as e:
                    print(e)
                    continue


    # Set the other ids to the passed_stars class
    for i,passed_matched_star in enumerate(passed_matched_stars):
        for query_result in query_results:

            if query_result['SCRIPT_NUMBER_ID']==i+1:
                passed_matched_stars[passed_matched_star].add_other_ids(
                    query_result['IDS'])
                break
    return query_results

#%%
##
def compare_ref_ids_with_matched_stars(passed_matched_stars,
                                       reference_star_table):
    passed_reference_star_table=Table(names=reference_star_table,dtype=reference_star_table)

    # Unpack all passed_matched_stars
    passed_id_list=[]
    for passed_matched_star in passed_matched_stars:
        if hasattr(passed_matched_stars[passed_matched_star],'other_ids'):
            for other_id in passed_matched_stars[
                passed_matched_star].other_ids.split('|'):
                    passed_id_list.append(other_id)



    #TODO: Find better way to match IDs to star.
    #for reference_star in reference_star_table:
        # Search reference Star ID in

    return passed_id_list



####################### Testing the Functions ################################

## Singular Request

# inner_rad=32
# outer_rad=48
#
# image_dir=r"D:\School\Work - Winter 2022\Work\2021-03-21\HIP 2894\LIGHT\B"
# catalogue_dir=r"D:\School\StarCatalogues\USNO UCAC4"
# matched_star_dictionary,matched_star_collection=get_matched_stars(image_dir, catalogue_dir, 15 ,3, \
#                         0.8, 2, 60, 11, False, 32, 48)
# passed_matched_stars=statistics_of_matched_stars(matched_star_collection,
#                             r"D:\School\Work - Winter 2022\Work\2021-03-21\HIP 2894\LIGHT\B\flux_calculations.txt",
#                             inner_rad,outer_rad)
# plot_measure_flux(matched_star_dictionary,r"D:\School\Work - Winter 2022\Work\2021-03-21\HIP 2894\LIGHT\B\fluxplots.png",
#                   passed_matched_stars,r"D:\School\Work - Winter 2022\Work\2021-03-21\HIP 2894\LIGHT\B\passed_stars_fluxplots.png"
#
#                                         )
#
# # Reset the Variables when complete
# matched_star_dictionary={}

#%%
## Mutliple Requests

image_path = r"C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-03-21"
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

    matched_star_dictionary,matched_star_collection=get_matched_stars(dirs, catalog_dir, 15,
                                                     3, 0.8, 2, 60, 11,
                                                     False, inner_rad, outer_rad,)

    passed_matched_stars=statistics_of_matched_stars(
        matched_star_collection,
        dirs+r"\flux_calculations.txt",inner_rad,
        outer_rad)
    plot_measure_flux(matched_star_dictionary,
        dirs+r"\allfluxplots.png",passed_matched_stars,
        dirs+r"\passed_stars_fluxplots.png")


    # Get Other IDS of passed_matched stars using ra and dec
    find_ids_for_matched_stars(passed_matched_stars)


    # Get Other IDS of reference stars
    reference_star_table, reference_formatted_ras, reference_formated_decs, \
    reference_other_ids = find_ids_for_ref_stars(refstar_dir)

    # See if matched stars ids are in reference star and create new
    # reference star
    compare_ref_ids_with_matched_stars(passed_matched_stars,reference_star_table)


##


