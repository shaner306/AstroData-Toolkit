# -*- coding: utf-8 -*-
"""
Batch AstroSolves Images

Created on Mon Mar 21 11:48:29 2022

@author: mstew
"""
import os
import Main
import AstroFunctions as astro
from astropy.nddata import CCDData
from pathlib import Path


#%% Batch Reduce Images

# Manually set the image Reduction Parameters

#bias_frames=r'C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-09-17 - processed\2021-09-17 - unprocessed\2022 01 17 - Bias - 3x3 - 0 sec'
#dark_frames=r'C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-09-17 - processed\2021-09-17 - unprocessed\2021 09 17 - Dark - 3x3 - 10 sec'
#flat_frames=r'C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-09-17 - processed\2021-09-17 - unprocessed\2021 09 17 - Flats - 3x3'
path=r'C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-09-17 - uncorrected\2021-09-17 - non-reduced'

create_master_dark=False
create_master_flat=False
create_master_bias=False
correct_outliers_params = {'Outlier Boolean': False,

                           'Hot Pixel': False,
                           'Dark Frame Threshold Bool': False,
                           'Dark Frame Threshold Min':  -50,
                           'Dark Frame Threshold Max': 100,
                           'ccdmask': False,
                           'Cosmic Rays Bool': False,
                           'Replace Bool': False,
                           'Replace Mode':  'Interpolate',
                           'Multiple Flat Combination':False,
                           'Save Corrected Flats': False,
                           'Radius of local Averaging': 1,
                           }

create_master_dir=False
use_existing_masters=True
exisiting_masters_dir=r'C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-09-17 - uncorrected\master_frame_data'

scalable_dark_bool=True


list_subfolders_with_paths= [f.path for f in os.scandir(path) if f.is_dir()]



for dirs in list_subfolders_with_paths:
    reduced_dirs=[dirs,exisiting_masters_dir]
    
    sav_loc=Path(str(dirs)+'_Corrected')
    sav_loc.mkdir(exist_ok=True)
    
    Main.Image_reduce(reduced_dirs,
                      create_master_dark,
                      create_master_flat,
                      create_master_bias,
                      correct_outliers_params,
                      create_master_dir,
                      use_existing_masters,
                      exisiting_masters_dir,
                      scalable_dark_bool,
                      sav_loc
                      )


        
#%% Batch Solve           
dataset_folder=r'C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-09-17 - processed\2021-09-17 - non-reduced'
catalog_dir=r'C:\Users\mstew\Documents\School and Work\Winter 2022\Work\StarCatalogues\USNO UCAC4'
refstar_dir=r'C:/Users/mstew/Documents/GitHub/Astro2/Reference Star Files/Reference_stars_2022_02_17_d.txt'


# Pinpoint Solve Parameters
max_mag=13
sigma=3
max_solve_time=60
match_residual=1.5
catalog=11
catalog_exp=0.8
use_sextractor=True
all_sky_solve=False
space_based_bool=False




path=dataset_folder

list_subfolders_with_paths= [f.path for f in os.scandir(path) if f.is_dir()]

for dirs in list_subfolders_with_paths:
    
    
    for dirpath,dirnames,files in os.walk(dirs):
        for name in files:
            if name.lower().endswith(('.fits','.fit','.fts')):
                sample_image=os.path.join(dirpath,name)
                break
    Sample_image = CCDData.read(sample_image,unit='adu')
    
    Main.pinpoint_solve(dirs,
                        catalog_dir,
                        max_mag,
                        sigma,
                        catalog_exp,
                        match_residual,
                        max_solve_time,
                        catalog,
                        space_based_bool,
                        use_sextractor,
                        all_sky_solve)
    
    try:
        plot_results = True
        save_plots = True
        exposure_key = 'EXPTIME'
        name_key = 'Name'
        unique_id = 'GBO'
        # For St. John's
        # # Uncomment these lines if you are processing
        # data from St. John's.
        file_suffix = (".fits", ".fit", ".fts")
        lat_key = 'SITELAT'
        lon_key = 'SITELONG'
        elev_key = 'SITEELEV'
        # For ORC GBO
        # # Uncomment these lines if you are processing data
        # from the ORC GBO.
        # file_suffix = ".fit"
        # lat_key = 'OBSGEO-B'
        # lon_key = 'OBSGEO-L'
        # elev_key = 'OBSGEO-H'
        # elev_key = 'OBSGEO-H'
        try: 
            if Sample_image.meta['Correctd'] == True:
                save_loc = os.path.join(dirs, 'Corrected_Outputs')
            else:
                save_loc = os.path.join(dirs, 'Outputs')
        except KeyError:
            save_loc = os.path.join(dirs, 'Outputs')
                
        Warner_final_transform_table =\
            astro._main_gb_transform_calc_Warner(
                dirs,
                refstar_dir,
                plot_results=plot_results,
                save_plots=save_plots,
                file_suffix=file_suffix,
                exposure_key=exposure_key,
                name_key=name_key, lat_key=lat_key,
                lon_key=lon_key, elev_key=elev_key,
                save_loc=save_loc, unique_id=unique_id)
    except:
        print('Error')
        