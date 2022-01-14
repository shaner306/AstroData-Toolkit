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
imagefolder = 0
catalogfolder = 0
refdoc = 0


def Gui():
    windowopen = False
    sg.theme("Default1")
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
                    sg.Input("D:\Solved Stars\Tycho 3023_1724",
                             key="-IN2-",
                             change_submits=True),
                    sg.FolderBrowse(key="-IN1-"), sg.Text('')],
                   [sg.Text("Catalog Folder:    "),
                    sg.Input("D:\\squid\\UCAC4",
                             key="-IN3-",
                             change_submits=True),
                    sg.FolderBrowse(key="-IN4-")],
                   [sg.Text("Reference Stars:  "),
                    sg.Input
                    (r'D:\Astro2\Reference Star Files\
                     Reference_stars_Apr29.txt',
                     key="-IN5-", change_submits=True),
                    sg.FileBrowse(key="-IN6-")],
                   [sg.T(""), sg.Checkbox('Save Data to Folder',
                                          default=True,
                                          key="-IN7-")],
                   [sg.T(""), sg.Checkbox('Pinpoint Solve',
                                          default=False,
                                          key="-IN100-")],
                   [sg.Column(column1),
                    sg.Column(column2)],
                   [sg.T("   ")],
                   [sg.T("   "), sg.Button("Solve"), sg.Cancel()]]

    tab2_layout = [
                   [sg.T('Image Reduction Script')],
                   [sg.T("   ")],
                   [sg.Text("Image Folder:    "),
                    sg.Input(key="-IN200-",
                             change_submits=True),
                    sg.FolderBrowse(key="-IN120-"), sg.Text('')],
              # [sg.Text("Bias Images:    "), 
              #  sg.Input(key="-IN20-" ,change_submits=True), 
              #  sg.FolderBrowse(key="-IN12-"), sg.Text('')],
              # [sg.Text("Dark Images:    "), 
              #  sg.Input(key="-IN30-" ,change_submits=True), 
              #  sg.FolderBrowse(key="-IN40-")],
              # [sg.Text("Flat Images:     "), 
              #  sg.Input(key="-IN50-" ,change_submits=True), 
              #  sg.FileBrowse(key="-IN60-")],
                  [sg.T(""), sg.Checkbox('Use Flats', 
                                         default=True, 
                                         key="-IN71-"),
                   sg.T(""), sg.Checkbox('Use Darks', 
                                         default=True,
                                         key="-IN1010-"),
                     sg.T(""), sg.Checkbox('Use Bias',
                                           default=True,
                                           key="-IN1011-")],
                  [sg.T(""), sg.Checkbox('Space Based',
                                         default=True,
                                         key="-IN1012-"),
                       sg.T("Target:"),
                       sg.InputText('SA-111',
                                    size=(10, 5),
                                    key="-IN1013-")],
                  [sg.T("   ")],
                  [sg.T(" "), sg.Button("Reduce"), sg.Cancel()]]

# Layout
    layout = [[sg.TabGroup([[sg.Tab('AstroSolver', tab1_layout),
               sg.Tab('Image Reduction', tab2_layout)
                   ]])]]
    if windowopen is False:
        window = sg.Window('AstroSolver', layout)
        windowopen is True
    while True:
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
            reduce_dir = values["-IN200-"]
            create_master_flat = values["-IN71-"]
            create_master_dark = values["-IN1010-"]
            create_master_bias = values["-IN1011-"]
            target = values["-IN1013-"]
            if values["-IN1012-"]:
                try:
                    Main.DarkSub(target, reduce_dir,
                                 'D:\\NEOSSat-SA-111\\test')
                    print("Reduce Space-Based Images ---- Started")
                    window.close()
                except:
                    print("Input Error")
                    window.update()
            else:
                try:
                    Main.Image_reduce(reduce_dir,
                                      create_master_dark,
                                      create_master_flat,
                                      create_master_bias)
                    print("Reduce Space-Based Images ---- Started")
                    window.close()
                except:
                    print("Input Error")
                    window.update()

        elif event == "Solve":
            image_dir = values["-IN2-"]
            catalog_dir = values["-IN3-"]
            refstar_dir = values["-IN5-"]
            save_data = values["-IN7-"]
            plot_data = values["-IN1013-"]
            
            if values["-IN82-"] is True:
                sb = 0
            else:
                sb = 1
            window.Refresh()
            if values["-IN100-"]:
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
                # try:    
                Main.pinpoint_solve(image_dir,
                                    catalog_dir,
                                    max_mag,
                                    sigma,
                                    catalog_exp,
                                    match_residual,
                                    max_solve_time,
                                    catalog,
                                    sb,
                                    use_sextractor,
                                    all_sky_solve)
                # print ("Reducing Images ---- Started")
                window.close()
                # except:
                #     print("Pinpoint Error. Please See Instructions")
                #     window.refresh()
            else:
                # print("yes")
                pass
                # Main.Ground_based_transforms(image_dir,refstar_dir)
            # print("yes")
            if values["-IN91-"] is True:
                if values["-IN82-"] is True:
                    try:
                        plot_results = True
                        save_plots = True
                        exposure_key = 'EXPTIME'
                        name_key = 'Name'
                        unique_id = 'GBO'
                        # For St. John's
                        # # Uncomment these lines if you are processing
                        # data from St. John's.
                        file_suffix = ".fits"
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
                        save_loc = os.path.join(image_dir, 'Outputs')
                        Warner_final_transform_table =\
                            astro._main_gb_transform_calc_Warner\
                            (image_dir,
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
                        window.close()
                    except Exception as e:
                        print(e)
                        print("Input Error. Please See Instructions")
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
                                                  max_distance_from_sat=
                                                  max_distance_from_sat,
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
                    print("Input Error. Please See Instructions")
                    # window.update()
                    continue
    window.close()
GUI = Gui()
