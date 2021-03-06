# -*- coding: utf-8 -*-
"""
Batch AstroSolves Images

Created on Mon Mar 21 11:48:29 2022

@author: mstes

"""
'''
Batch Reduced Image Script

Batch_solves purpose is to enable the user to solve multiple catalogs of data.
Each module is broken down into sections which may be operated individual or
 the script can be run as a whole if the whole image processing pipeline is 
 needed.

'''

## Import Section
import os
from pathlib import Path

from astropy.nddata import CCDData

import AstroFunctions as astro
import Main
import main_transforms as transforms
#import pinpoint

# %% Batch Reduce Images
##

# Manually set the image Reduction Parameters
# User must already have prebuilt master files in order to avaoid repeatititon
# of master_frame creation. Developement of master frames can be done using
# the GUI and only inputting Bias, Dark and Flat Frames

'''
 Steps: Create Master frame data manually using the GUI
____
image_path: string
 must contain only light images

create_master_dark: Boolean
    default = False
create_master_flat=False
create_master_bias=False

'''
def batch_reduced_images(create_master_dark,
                         create_master_flat,
                         create_master_bias,
                         create_master_dir,
                         correct_outliers_params,

                         use_existing_masters,
                         exisiting_masters_dir,
                         scalable_dark_bool,
                         ):

    list_subfolders_with_paths = [f.path for f in os.scandir(image_path) if f.is_dir()]

    for dirs in list_subfolders_with_paths:
        sav_loc = Path(str(dirs) + '_Outlier_Corrected')
        sav_loc.mkdir(exist_ok=True)
        reduced_dirs = [dirs, exisiting_masters_dir]
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

#%% Call the Image Reduciton Module
##

# Define the Science image Path and image Path

#bias_frames=r'C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-09-17 - processed\2021-09-17 - unprocessed\2022 01 17 - Bias - 3x3 - 0 sec'
#dark_frames=r'C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-09-17 - processed\2021-09-17 - unprocessed\2021 09 17 - Dark - 3x3 - 10 sec'
#flat_frames=r'C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-09-17 - processed\2021-09-17 - unprocessed\2021 09 17 - Flats - 3x3'
image_path= r"D:\School\Work - Winter 2022\Work\2022-03-16\2022-01-16- Raw\SiderialStareMode"

# Set Parameters

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
                           'Force Offset':False
                           }

create_master_dir=False

use_existing_masters=True
exisiting_masters_dir=r"D:\School\Work - Winter 2022\Work\2022-03-16\2022-01-16- Raw\master_frame_data"


scalable_dark_bool=True

batch_reduced_images(create_master_dark,
                     create_master_flat,
                     create_master_bias,
                     create_master_dir,
                     correct_outliers_params,

                     use_existing_masters,
                     exisiting_masters_dir,
                     scalable_dark_bool)

list_subfolders_with_paths= [f.path for f in os.scandir(image_path) if f.is_dir()]



for dirs in list_subfolders_with_paths:
    sav_loc = Path(str(dirs) + '_Corrected')
    sav_loc.mkdir(exist_ok=True)
    if use_existing_masters:
        reduced_dirs=[dirs,exisiting_masters_dir]
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

    else:
        reduced_dirs=[dirs,exisiting_masters_dir]
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

    



        
# %% Batch Solve
## Batch Solve
#
'''
Batch Solving is built to iterate through the Starfields in a certain dataset 
instead of only iterating through a singular starfield

Current Example Structure:

Mar 16 2022
    |
    ---> SA23
           |
           ---> B
                |
                ---> *.fits
           |
           ---> G
           ````
    
'''

dataset_folder=r"D:\School\Work - Winter 2022\Work\2022-03-16\2022-03-16 - Copy"
catalog_dir=r"D:\School\StarCatalogues\USNO UCAC4"
refstar_dir=r"C:\Users\stewe\Documents\GitHub\Astro2\Reference Star Files\Reference_stars_2022_02_17_d.txt"


# Pinpoint Solve Parameters
max_mag=13
sigma=3
max_solve_time=60 # Seconds
match_residual=1.5
catalog=11
catalog_exp=0.8
use_sextractor=False
all_sky_solve=False
space_based_bool=False
photometry_method='aperture'
aperture_estimation_mode='mean'




image_path=dataset_folder

list_subfolders_with_paths= [f.path for f in os.scandir(image_path) if f.is_dir()]


for dirs in list_subfolders_with_paths:
    
    
    for dirpath,dirnames,files in os.walk(dirs):
        for name in files:
            if name.lower().endswith(('.fits','.fit','.fts')):
                sample_image=os.path.join(dirpath,name)
                break
    Sample_image = CCDData.read(sample_image,unit='adu')
    
# TODO: Create script that determines if WCS data is in the image
    # pinpoint.pinpoint_solve(dirs,
    #                     catalog_dir,
    #                     max_mag,
    #                     sigma,
    #                     catalog_exp,
    #                     match_residual,
    #                     max_solve_time,
    #                     catalog,
    #                     space_based_bool,
    #                     use_sextractor,
    #                     all_sky_solve)


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
                
        NewBoydMethod=transforms._main_gb_new_boyd_method(
          dirs,
          refstar_dir,
          plot_results=plot_results,
          save_plots=save_plots,
          file_suffix=file_suffix,
          exposure_key=exposure_key,
          name_key=name_key, lat_key=lat_key,
          lon_key=lon_key, elev_key=elev_key,
          save_loc=save_loc, unique_id=unique_id,photometry_method=photometry_method,
            aperture_estimation_mode=aperture_estimation_mode)
    except Exception as e:
        print(e)
        continue
     
# %% Create Combined Large Star Table
##
import csv
dataset_folder=r'C:\Users\mstew\Documents\School and Work\Winter 2022\Work\2021-10-30\Siderial Stare Mode Reduced No Dark Scale'
first_switch=True
file=dataset_folder+"\\" + dataset_folder.split('\\')[-1] +"_Combined_Large_Star_Table.csv"
with open(file,"a+",newline='\n') as f:
    writer=csv.writer(f,delimiter=',')
    for dirpath,dirname,file in os.walk(dataset_folder):
        for name in file:
            if name=='large_stars_table.csv':
                f2=open(os.path.join(dirpath,name),'r',newline='\n')
                csvreader=csv.reader(f2,delimiter=',')
                if first_switch==True:
                    print('first line')
                    first_switch=False
                else:
                    next(csvreader)
                for row in csvreader:
                    writer.writerow(row)
    f.close()
    
# %% Create Combined Boyd Tables
##
import csv
import os
dataset_folder=r"D:\School\Work - Winter 2022\Work\2022-03-16\2022-03-16 Results\UncappedDynamicSources"
first_switch=True
file=dataset_folder+"\\" + dataset_folder.split('\\')[-1] +"Boyde_Table1_Combined.csv"
# =============================================================================
# file=dataset_folder+"\\" + dataset_folder.split('\\')[-1] +"\\Boyde_Table1_Combined.csv"
# if os.path.isdir(file) == True:
#     with open(file,"w+",newline='\n') as f:
#         writer=csv.writer(f,delimiter=',')
#         for dirpath,dirname,file in os.walk(dataset_folder):
#             for name in file:
#                 if name=='Boyde_Table1.csv':
#                     f2=open(os.path.join(dirpath,name),'r',newline='\n')
#                     csvreader=csv.reader(f2,delimiter=',')
#                     if first_switch==True:
#                         print('first line')
#                         first_switch=False
#                     else:
#                         next(csvreader)
#                     for row in csvreader:
#                         writer.writerow(row)
# =============================================================================
with open(file,"a+",newline='\n') as f:
    writer=csv.writer(f,delimiter=',')
    for dirpath,dirname,file in os.walk(dataset_folder):
        for name in file:
            if name=='Boyde_Table1.csv':
                f2=open(os.path.join(dirpath,name),'r',newline='\n')
                csvreader=csv.reader(f2,delimiter=',')
                if first_switch==True:
                    print('first line')
                    first_switch=False
                else:
                    next(csvreader)
                for row in csvreader:
                    writer.writerow(row)

# %% Perform steps 2 and 3 of Boyd Method
##
from astropy.table import Table
from astropy.io import ascii
import csv
import os
import AstroFunctions as astro
import auxilary_phot_boyde_functions as aux_boyde

#dataset_folder=r'D:\School\Work - Winter 2022\Work\2021-04-21\Siderial Stare Mode\Post'

#file=dataset_folder+"\\" + dataset_folder.split('\\')[-1] +"Boyde_Table1_Combined.csv"
file=r"D:\School\Work - Winter 2022\Work\2022-03-16\2022-03-16 Results\UncappedDynamicSources\UncappedDynamicSourcesBoyde_Table1_Combined.csv"
dataset_folder=os.path.dirname(file)
header=[]

with open(file,'r',newline='\n') as f:
    csvreader=csv.reader(f)
    header=next(csvreader)
    Big_Boyde_Table=Table(names=header,dtype=['str','float64','float64','str','float64','float64','str','float64','int'])
    for row in csvreader:
        Big_Boyde_Table.add_row(row)
        
Boyde_Table2=aux_boyde.calculate_boyde_slope_2(Big_Boyde_Table,str(dataset_folder+'\\Output'),4,save_plots=True)
ascii.write(Boyde_Table2, os.path.join(
    str(os.path.dirname(file)), 'Combined_Boyde_Table2.csv'), format='csv')

#% SAve Coefficeint Data
###
# FIXME: Get sorted_coefficient_data to output correctly

sorted_coefficient_data=aux_boyde.create_coefficeint_output(Boyde_Table2,file)
