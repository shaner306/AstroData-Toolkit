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
import os, os.path
import Main
imagefolder=0
catalogfolder=0
refdoc = 0



def Gui ():
    windowopen=False
    sg.theme("Default1")
   
    column2 = [[sg.Text('Pinpoint Solve Parameters', background_color='#F7F3EC', justification='center', size=(30, 1))],
               
            [sg.T("Maximum Mag.:"),sg.T("     "), sg.InputText('13', size=(5,5), key="-IN51-"), sg.T("Sigma Above Mean:  "), sg.InputText('2.5', size=(5,5), key="-IN52-")],      
            [sg.T("Max Solve Time (sec):"), sg.InputText('60', size=(5,5),key="-IN56-"), sg.T("Max Match Residual:"), sg.InputText('1.5', size=(5,5), key="-IN53-")],
            [sg.T("Catalog Expansion:"),sg.T(""), sg.InputText('0.8', size=(5,5),key="-IN55-"), sg.T("Catalog:"), sg.T("            "), sg.InputText('UCAC4', size=(7,5), disabled=True,  key="-IN54-")]
            ]   
    
    column1 = [[sg.Text('Analysis Parameters', background_color='#F7F3EC', justification='center', size=(30, 1))],
               
            [sg.Frame('Source Capture Mode',[[
              sg.T(""), sg.Radio('Star Stare Mode', "RADIO3", default=True),
               sg.T("            "), sg.Radio('Track Rate Mode', "RADIO3", default=False, key="-IN92-",disabled=True)]], size=(200,100))],
                  
              [sg.Frame('Image Source',[[      
              sg.T(""), sg.Radio('Ground Based', "RADIO1", default=True, key="-IN82-"),
              sg.T("                "), sg.Radio('Space Based', "RADIO1", default=False, disabled=True, key="-IN83-")]])], 
              ] 
    #column2.update(disabled=True)
    
    tab1_layout = [
    [sg.T("AstroSolver Processor")],[sg.T("Version 0.1")], [sg.T("   ")],
              [sg.Text("Image Folder:      "), 
               sg.Input("D:\Solved Stars\Tycho 3023_1724", key="-IN2-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN1-"), sg.Text('')],
              [sg.Text("Catalog Folder:    "), 
               sg.Input(key="-IN3-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN4-")],
              [sg.Text("Reference Stars:  "), 
               sg.Input ("D:\Astro2\Reference Star Files\Reference_stars.csv", key="-IN5-" ,change_submits=True), 
               sg.FileBrowse(key="-IN6-")],
              
              
              [sg.T(""), sg.Checkbox('Save Data to Folder', default=True, key="-IN7-")],
              [sg.T(""), sg.Checkbox('Pinpoint Solve', default=False, key="-IN100-")],
              [sg.Column(column1),sg.Column(column2)],
              [sg.T("   ")],
              [sg.T("   "), sg.Button("Solve"), sg.Cancel()]
              ]


    tab2_layout = [
    [sg.T('Image Reduction Script')],
    [sg.T("   ")],
    [sg.Text("Image Folder:    "), 
               sg.Input(key="-IN200-" ,change_submits=True), 
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
              
              [sg.T(""), sg.Checkbox('Use Flats', default=True, key="-IN71-"),
              sg.T(""), sg.Checkbox('Use Darks', default=True, key="-IN1010-"),
              sg.T(""), sg.Checkbox('Use Bias', default=True, key="-IN1011-")],
              [sg.T(""), sg.Checkbox('Space Based', default=True, key="-IN1012-"), sg.T("Target:"), sg.InputText('SA-111', size=(10,5),key="-IN1013-")],      
            
              [sg.T("   ")],
              [sg.T(" "), sg.Button("Reduce"), sg.Cancel()]]
    

    
    layout = [[sg.TabGroup([[sg.Tab('AstroSolver', tab1_layout),
                   sg.Tab('Image Reduction', tab2_layout)
                   ]])]]
    
    
    
        # f.Catalog = 5;
        #         f.CatalogPath = catloc;
        #         f.CatalogMaximumMagnitude = 13;
        #         f.CatalogExpansion = 0.8;
        #         f.SigmaAboveMean = 2.5; 
        #         f.FindImageStars; 
        #         f.FindCatalogStars; 
        #         f.MaxSolveTime = 60; 
        #         f.MaxMatchResidual = 1.5; 
    #window['-OUTPUT-'].update(values['-IN-'])
        # simple version for working with CWD
        
        
       

        # path joining version for other paths
    #DIR = 'D:\\Wawrow\\2. Observational Data\\2021-03-10 - Calibrated\\HIP 46066\\LIGHT\\B\\'
    #path, dirs, files = next(os.walk(DIR))
    #file_count = len(files)
    
    # layout = [[sg.T(" ")], 
    #           [sg.Text("Image Folder: "), 
    #            sg.Input(key="-IN2-" ,change_submits=True), 
    #            sg.FolderBrowse(key="-IN1-")],
    #           [sg.Text("Catalog Folder: "), 
    #            sg.Input(key="-IN3-" ,change_submits=True), 
    #            sg.FolderBrowse(key="-IN4-")],
    #           [sg.Text("Reference Stars: "), 
    #            sg.Input(key="-IN5-" ,change_submits=True), 
    #            sg.FileBrowse(key="-IN6-")],
    #          [sg.T("                   "), sg.Checkbox('Print On:', default=True, key="-IN7-")],
    #           [sg.T("         "), sg.Radio('Permission Granted', "RADIO1", default=False, key="-IN8-")],
    #           [sg.T("         "), sg.Radio('Permission not Granted', "RADIO1", default=True)],
    #           [sg.Radio('Track Rate Mode', "RADIO3", default=False, key="-IN9-"),
    #            sg.Radio('Star Stare Mode', "RADIO3", default=True)],
    #           [sg.Button("Submit")]]
    
    ###Building Window
    if windowopen==False:
        window = sg.Window('AstroSolver', layout)
        windowopen==True
    
        
    while True:
        event, values = window.read()
        # window["-IN56-"].update("hello")
        
        
        # if values["-IN100-"]:
        #     window["-IN56-"].update(disabled=True)
        #     print("test")
           
                        
                        
        #print(values["-IN2-"])
        if event == sg.WIN_CLOSED or event=="Exit":
            window.close()
            break
        
        elif event == "Reduce":
            
            reduce_dir =values["-IN200-"]
            create_master_flat = values["-IN71-"]
            create_master_dark = values["-IN1010-"]
            create_master_bias = values["-IN1011-"]
            target = values["-IN1013-"]
            
            window.close()
            
            if values["-IN1012-"]:
               Main.DarkSub(target, reduce_dir, 'D:\\NEOSSat-SA-111\\clean2')

              
            else:
                Main.Image_reduce(reduce_dir, create_master_dark, create_master_flat, create_master_bias)
          
                
          
        elif event == "Solve":
            image_dir =values["-IN200-"]
            catalog_dir = values["-IN71-"]
            refstar_dir = values["-IN1010-"]
            save_data = values["-IN1011-"]
            plot_data = values["-IN1013-"]
            
            Main.Ground_based_transforms(image_dir,refstar_dir)
            
            
            if values["-IN100-"]:
                max_mag = values["-IN51-"]
                sigma= values["-IN52-"]
                match_residual = values["-IN53-"]
                catalog = values["-IN54-"]
                catalog_exp = values["-IN55-"]
                max_solve_time = values["-IN56-"]
                Main.pinpoint_solve(image_dir, catalog_dir, max_mag, sigma, catalog_exp, match_residual, max_solve_time, catalog)
            else:
                continue
            if values["-IN82-"]:
                Main.Ground_based_transforms(image_dir,refstar_dir)
                
            else:
                
                Main.space_based_transforms(image_dir,refstar_dir)
                
                
            window.close()
            
            
            
            
            
    window.close()           
            
GUI = Gui()
