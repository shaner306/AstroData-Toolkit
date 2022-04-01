# -*- coding: utf-8 -*-
"""
Created on Fri May 28 11:56:48 2021

@author: shane

Find:
    Image Folder
    Catalog Folder
        What Catalog
    Reference Star Sheet
    Ground Based or Space Based
    TRM or SSM
    Add Moffat and Gaussian Distributions
    Background Subtraction Method
    Pinpoint Astrometry Solve
    Print Docs
    Astro Reducer Method
"""

import PySimpleGUI as sg
import os
import os.path
import Main
import AstroFunctions as astro
from astropy.nddata import CCDData
from pathlib import Path
imagefolder = 0
catalogfolder = 0
refdoc = 0


def Gui():
    windowopen = False
    sg.theme("Default1")

    # Pinpoint Solve Parameters
    column2 = [[sg.Text(
        'Pinpoint Solve Parameters',
        background_color='#F7F3EC',
        justification='center',
        size=(30, 1))],
        [sg.T("Maximum Mag.:"),
         sg.T("     "),
         sg.InputText('13', size=(5, 5),
                      key="-IN51-"),
         sg.T("Sigma Above Mean:  "),
         sg.InputText('3.0',
                      size=(5, 5),
                      key="-IN52-")],
        [sg.T("Max Solve Time (sec):"),
         sg.InputText('60',
                      size=(5, 5),
                      key="-IN56-"),
         sg.T("Max Match Residual:"),
         sg.InputText('1.5',
                      size=(5, 5),
                      key="-IN53-")],
        [sg.T("Catalog Expansion:"),
         sg.T(""),
         sg.InputText('0.8',
                      size=(5, 5),
                      key="-IN55-"),
         sg.T("Catalog:"),
         sg.T("            "),
         sg.InputText('UCAC4',
                      size=(7, 5),
                      disabled=True,
                      key="-IN54-")],
        [sg.T("Use SExtractor?       "),
         sg.Radio('Yes',
                  "RADIO2",
                  default=False,
                  key="-IN57-"),
         sg.T(" "),
         sg.Radio('No',
                  "RADIO2",
                  default=True,
                  key="-IN58-")],
        [sg.T("All Sky Solve?          "),
         sg.Radio('Yes',
                  "RADIO4",
                  default=False,
                  key="-IN59-"),
         sg.T(" "),
         sg.Radio('No',
                  "RADIO4",
                  default=True,
                  key="-IN60-")]]
    # Analysis Parameters
    column1 = [[sg.Text('Analysis Parameters',
                        background_color='#F7F3EC',
                        justification='center',
                        size=(30, 1))],
               [sg.Frame('Source Capture Mode',
                         [[sg.T(""),
                           sg.Radio('Star Stare Mode',
                                    "RADIO3",
                                    default=True,
                                    key="-IN91-"),
                           sg.T("            "),
                           sg.Radio('Track Rate Mode',
                                    "RADIO3",
                                    default=False,
                                    key="-IN92-")]],
                         size=(200, 100))],
               [sg.Frame('Image Source',
                         [[sg.T(""),
                           sg.Radio('Ground Based',
                                    "RADIO1",
                                    default=True,
                                    key="-IN82-"),
                           sg.T("                "),
                           sg.Radio('Space Based',
                                    "RADIO1",
                                    default=False,
                                    key="-IN83-")]])], ]
    # column2.update(disabled=True)
    tab1_layout = [[sg.T("Input Folders")],
                   [sg.T("   ")],
                   [sg.Text("Image Folder:      "),
                    sg.Input("C:/Users/mstew/Documents/School and Work/Winter 2022/Work/Suffield Data/2022-01-17 - Amazonas 2 and SA/SA26/LIGHT",
                             key="-IN2-",
                             change_submits=True),
                    sg.FolderBrowse(key="-IN1-"), sg.Text('')],
                   [sg.Text("Catalog Folder:    "),
                    sg.Input(r'C:/Users/mstew/Documents/School and Work/Winter 2022/Work/StarCatalogues/USNO UCAC4',
                             key="-IN3-",
                             change_submits=True),
                    sg.FolderBrowse(key="-IN4-")],
                   [sg.Text("Reference Stars:  "),
                    sg.Input
                    (r'C:/Users/mstew/Documents/GitHub/Astro2/Reference Star Files/Reference_stars_2022_02_17_d.txt',
                     key="-IN5-", change_submits=True),
                    sg.FileBrowse(key="-IN6-")],
                   [sg.T(""), sg.Checkbox('Save Data to Folder',
                                          default=True,
                                          key="-IN7-")],
                   [sg.T(""), sg.Checkbox('Pinpoint Solve',
                                          default=False,
                                          key="-IN100-")],
                   [sg.T(""),sg.Checkbox('New Boyd Method',
                                         default=False,
                                         key='NewBoydMethod')],
                   # 1N100- PinPoint Solve
                   [sg.Column(column1),
                    sg.Column(column2)],
                   [sg.T("   ")],
                   [sg.T("   "), sg.Button("Solve"), sg.Cancel()]]

    tab2_column1 = [[sg.Text('Correct Outliear Parameters',
                             background_color='#F7F3EC',
                             justification='center',
                             size=(30, 1))],
                    [sg.Checkbox('Hot Pixel Removal',
                                 default=False,
                                 key="-IN1012-1-")],
                    [sg.Checkbox('Dark Frame Threshold Bool',
                                 default=False,
                                 key="-IN1012-2-")],

                    [sg.T("Dark Frame Threshold Min"),
                     sg.T("     "),
                     sg.InputText('-50', size=(5, 5),
                                  key="-IN1012-3-")],
                    [sg.T("Dark Frame Threshold Max"),
                     sg.InputText('100',
                                  size=(5, 5),
                                  key="-IN1012-4-")],
                    [sg.Checkbox("ccdmask",
                     default=False,
                     key="-IN1012-5-")],
                    [sg.Checkbox("Cosmic Rays Removal",
                     default=False,
                     key="IN1012-6")]
                    ]
    tab2_column2 = [[sg.Text('Replace Outliers Options',
                             background_color='#F7F3EC',
                             justification='center',
                             size=(30, 1))],
                    [sg.Checkbox('Replace Mask Values',
                                 default=True,
                                 key="1IN10121")],
                    [sg.Radio('Average Background',
                              "RADIO2",
                              default=True,
                              key="1IN10121-1"),

                    sg.Radio('Local Averaging',
                             "RADIO2",
                             default=False,
                             key="1IN10121-2")],
                    [sg.T("Radius of Local Averaging"),
                     sg.T("     "),
                     sg.InputText('1', size=(5, 5),
                                  key="1IN10121-2-1")],
                    [sg.Checkbox('Multiple Flat Combination',
                                 default=False,
                                 key='1IN10122')],
                    [sg.Checkbox('Save Corrected Flats',
                                default=False,
                                key='1IN10123')],
                    [sg.Checkbox('Force Offset',
                                 default=False,
                                 key='ForceOffset')],
                    ]

    tab2_layout = [
        [sg.T('Image Reduction Script')],
        [sg.T("   ")],
        [sg.Text("Image Folder:    "),
         sg.Input(key="-IN200-",
                  change_submits=True),
         sg.FolderBrowse(key="-IN120-"), sg.Text('')],
        [sg.Text("Bias Images:    "),
          sg.Input(key="-IN20-" ,change_submits=True),
          sg.FolderBrowse(key="-IN12-"), sg.Text('')],
        [sg.Text("Dark Images:    "),
          sg.Input(key="-IN30-" ,change_submits=True),
          sg.FolderBrowse(key="-IN40-")],
        [sg.Text("Flat Images:     "),
          sg.Input(key="-IN50-" ,change_submits=True),
          sg.FolderBrowse(key="-IN60-")],
        [sg.T(""), sg.Checkbox('Create Master Flats',
                               default=True,
                               key="-IN71-"),
         sg.T(""), sg.Checkbox('Create Master Darks',
                               default=True,
                               key="-IN1010-"),
         sg.T(""), sg.Checkbox('Create Master Bias',
                               default=True,
                               key="-IN1011-"),
         sg.T(""), sg.Checkbox('Correct for Outliers',
                               default=True,
                               key="-IN1012-")],
        [sg.T(""),sg.Checkbox('Use Exisitng Masters',
                              default=True,
                              key='-1N109-'),
         sg.Text("Master Images:    "),
           sg.Input(key="-IN 109-1" ,change_submits=True),
           sg.FolderBrowse(key="-1N109-2")],
                                              

        [sg.T(""), sg.Checkbox('Space Based',
                               default=True,
                               key="-IN1013-"),
         sg.T("Target:"),
         sg.InputText('SA-111',
                      size=(10, 5),
                      key="-IN1014-")],
        [sg.Column(tab2_column1), sg.Column(tab2_column2)],
        [sg.T("   ")],
        [sg.T(" "), sg.Button("Reduce"), sg.Cancel()]
    ]

# Layout
    layout = [[sg.TabGroup([[sg.Tab('AstroSolver', tab1_layout),
               sg.Tab('Image Reduction', tab2_layout)
    ]])]]
    if windowopen is False:
        window = sg.Window('AstroSolver', layout)
        windowopen is True
    while True:  # Read Events
        window.Refresh()
        event, values = window.read()
        # window["-IN56-"].update("hello")
        # if values["-IN100-"]:
        #     window["-IN56-"].update(disabled=True)
        #     print("test")
        # print(values["-IN2-"])
        if event == sg.WIN_CLOSED or event == "Exit":
            window.close()
            break
        elif event == "Reduce":
             
            use_existing_masters = values['-1N109-']
            exisiting_masters_dir = values['-1N109-2']
            
            reduce_dir = values["-IN200-"]
            if use_existing_masters is True:
                reduce_dirs = [values["-IN200-"],exisiting_masters_dir] 
            else:
                reduce_dirs = [values["-IN200-"],values["-IN30-"],values["-IN50-"], values["-IN20-"]] 
            
            
            
            
            use_existing_masters = values['-1N109-']
            exisiting_masters_dir = values['-1N109-2']
            
            #### FLAGS #####
            
            # TODO: Create Functions for this
            
            if os.path.isdir(reduce_dir):
                for dirpath,dirnames,files in os.walk(reduce_dir):
                    for name in files:
                        if name.lower().endswith(('.fits','.fit','.fts')):
                            sample_image_path=os.path.join(dirpath,name)
                            break
                sample_science_image = CCDData.read(sample_image_path,unit='adu')
                try:
                    if sample_science_image.header['Correctd'] is True:
                        Popup_string=sg.popup_yes_no("Images Are Already Reduced by this program, Continue?")
                        if Popup_string=='No':
                            window.close()
                            quit()
                except KeyError:
                    print('Could not find Correctd keyword')
            
            # TODO: Create Function for this        
            
            # Find Sample Dark
            
            if use_existing_masters:
                dark_look_dir = exisiting_masters_dir
            else:
                dark_look_dir= values["-IN30-"]
            try:    
                for dirpath,dirnames,files in os.walk(dark_look_dir):
                    for name in files:
                        if name.lower().endswith(('.fits','.fit','.fts')):
                            CCDData_sample=CCDData.read(os.path.join(dirpath,name),unit='adu')
                            
                            try: 
                                if "dark" in CCDData_sample.header["IMAGETYP"].lower():
                                    sample_dark_image = CCDData_sample
                                    break
                                else:
                                    continue 
                            except KeyError: 
                                raise KeyError("WARNING -- Cound not find Keyword when looking for Dark Frame")
            except:
                raise KeyError("WARNING -- Could not find Dark Sample image")
            try: 
                if sample_science_image.header["EXPTIME"]!=sample_dark_image.header["EXPTIME"]:
                    popup_string=sg.popup_yes_no("Science images might not have the same exposure time as the dark frames, Continue?")
                    if popup_string=='No':
                        break
                        window.close()
                        
                
            except UnboundLocalError:
                print('No Sample Prediction')
                
                                
            # TODO : Come up with better methods for improving this
            
            if use_existing_masters is True:
                
                    create_master_dir=False
            else:
                    create_master_dir=True
            
            # TODO: Come up with a better method for this 
            scalable_dark_bool=True
            if (values["-IN20-"] == '') and (values["-IN30-"]!='') and (values["-IN40-"]!=''):
                      
                if (len(reduce_dirs)==4) & (use_existing_masters is False):
                    print('None Scalabale Dark Detected')
                    scalable_dark_bool=False
                    del reduce_dirs[reduce_dirs.index('')]
                    print('Deleted Bias in redcued dir')
     
            create_master_flat = values["-IN71-"]
            create_master_dark = values["-IN1010-"]
            create_master_bias = values["-IN1011-"]
            
            if values["1IN10121-1"] is True:
                replace_mode = 'Ave'
            elif values["1IN10121-2"] is True:
                replace_mode = 'Interpolate'
                
            correct_outliers_params = {'Outlier Boolean': values["-IN1012-"],

                                       'Hot Pixel': values["-IN1012-1-"],
                                       'Dark Frame Threshold Bool': values["-IN1012-2-"],
                                       'Dark Frame Threshold Min':  values["-IN1012-3-"],
                                       'Dark Frame Threshold Max': values["-IN1012-4-"],
                                       'ccdmask': values["-IN1012-5-"],
                                       'Cosmic Rays Bool': values["IN1012-6"],
                                       'Replace Bool': values["1IN10121"],
                                       'Replace Mode': replace_mode,
                                       'Multiple Flat Combination':values["1IN10122"],
                                       'Save Corrected Flats': values ["1IN10123"],
                                       'Radius of local Averaging': values["1IN10121-2-1"],
                                       'Force Offset': values['ForceOffset']
                                       }

            
            correct_light_dir = Path(reduce_dirs[0], 'corrected_lights')
            correct_light_dir.mkdir(exist_ok=True)
            correct_light_directory = reduce_dirs[0] + '\\corrected_lights'
            sav_loc = correct_light_directory
            
            target = values["-IN1014-"]
            if values["-IN1013-"]:  # Space Based Observations is True
                try:
                    Main.DarkSub(target, reduce_dir,
                                 'D:\\NEOSSat-SA-111\\test')
                    print("Reduce Space-Based Images ---- Started")
                    break
                    window.close()
                except:
                    print("Input Error")
                    break
                    window.update()
            else:
                
                
                try:
                    
                    
                    
                    
                    Main.Image_reduce(reduce_dirs,
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
                    print("Reduce Space-Based Images ---- Started")
                    break
                    window.close()
                except Exception as ex:
                    print(ex)

        elif event == "Solve":
            
            image_dir = values["-IN2-"]
            catalog_dir = values["-IN3-"]
            refstar_dir = values["-IN5-"]
            save_data = values["-IN7-"]
            plot_data = values["-IN1014-"]
            
            # Check to See if image has been corrected already
            for dirpath,dirnames,files in os.walk(image_dir):
                for name in files:
                    
                    
                    if name.lower().endswith(('.fits','.fit','.fts')):
                        CCData_sample=CCDData.read(os.path.join(dirpath,name),unit='adu')
                        sample_science_image=CCData_sample
                        break
                        
            
            try: 
                if sample_science_image.header['Correctd'] is False:
                    Popup_string=sg.popup_yes_no("Images Aren't Reduced, Continue?")
                    if Popup_string=='No':
                        break
                        window.close()
                        
                    
                
                # Do Nothing if Image is already Reduced
                
            except:
                sample_science_image.meta['Correctd'] = False
                
                # Prompt User about Confirming to Solve the image despite it not being redcued
                Popup_string=sg.popup_yes_no("Images Aren't Reduced, Continue?")
                if Popup_string=='No':
                    break
                    window.close()
                    
                
                        
                
            
            
            
            if values["-IN82-"] is True:  # Ground Based Observation
                space_based_bool = 0   # Space Based Boolean=0
            else:
                space_based_bool = 1
            window.Refresh()
            if values["-IN100-"]:  # PinPoint Solve Key is True
                max_mag = values["-IN51-"]
                sigma = values["-IN52-"]
                match_residual = values["-IN53-"]
                catalog = values["-IN54-"]
                catalog_exp = values["-IN55-"]
                max_solve_time = values["-IN56-"]
                if values["-IN57-"] is True:
                    use_sextractor = True
                else:
                    use_sextractor = False
                if values["-IN59-"] is True:
                    all_sky_solve = True
                else:
                    all_sky_solve = False
                try:
                    print("Pinpoint Solve Images ---- Started")
                    Main.pinpoint_solve(image_dir,
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
                    
                    window.close()
                except:
                    print("Pinpoint Error. Please See Instructions")
                    window.refresh()
            else:
                # print("yes")
                pass
                # Main.Ground_based_transforms(image_dir,refstar_dir)
            # print("yes")
            if values["-IN91-"] is True:  # Star Stare Mode is True
                if values["-IN82-"] is True:  # Ground Based is True
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
                            
                        try: 
                            if sample_science_image.meta['Correctd'] == True:
                                save_loc = os.path.join(image_dir, 'Corrected_Outputs')
                            else:
                                save_loc = os.path.join(image_dir, 'Outputs')
                        except KeyError:
                            save_loc = os.path.join(image_dir, 'Outputs')
                        
                        ## New Boyd Method ##
                        if values['NewBoydMethod'] is False:
                            
                            Warner_final_transform_table =\
                                astro._main_gb_transform_calc_Warner(
                                    image_dir,
                                    refstar_dir,
                                    plot_results=plot_results,
                                    save_plots=save_plots,
                                    file_suffix=file_suffix,
                                    exposure_key=exposure_key,
                                    name_key=name_key, lat_key=lat_key,
                                    lon_key=lon_key, elev_key=elev_key,
                                    save_loc=save_loc, unique_id=unique_id)
                        else:    
                            NewBoydMethod=astro._main_gb_new_boyd_method(
                              image_dir,
                              refstar_dir,
                              plot_results=plot_results,
                              save_plots=save_plots,
                              file_suffix=file_suffix,
                              exposure_key=exposure_key,
                              name_key=name_key, lat_key=lat_key,
                              lon_key=lon_key, elev_key=elev_key,
                              save_loc=save_loc, unique_id=unique_id)
                        
                          
                            
                        # Main.Ground_based_transforms(image_dir,refstar_dir)
                        # print (image_dir,refstar_dir)
                        break
                        window.close()
                    except Exception as e:
                        print(e)
                        print("Input Error. Please See Instructions1")
                        # window.update()
                else:
                    try:
                        Main.space_based_transform(image_dir, refstar_dir)
                        print("Reducing Images ---- Started")
                        window.close()
                    except:
                        print(" space based Input Error.\
                              Please See Instructions")
                        # window.update()
            else:
                try:
                    temp_dir = 'tmp'
                    max_distance_from_sat = 20
                    size = 20
                    max_num_nan = 5
                    plot_results = 0
                    sats_table,\
                        uncertainty_table,\
                        sat_fwhm_table = \
                        astro._main_sc_lightcurve(image_dir,
                                                  temp_dir=temp_dir,
                                                  max_distance_from_sat=max_distance_from_sat,
                                                  size=size,
                                                  max_num_nan=max_num_nan,
                                                  plot_results=plot_results)
                    sat_fwhm_table.pprint_all()
                    uncertainty_table.pprint_all()
                    sats_table.pprint_all()
                    window.Refresh()
                    print(image_dir)
                    window.close()
                except:
                    print("Input Error. Please See Instructions2")
                    # window.update()
                    continue
    window.close()


GUI = Gui()
