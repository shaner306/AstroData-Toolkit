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

def Gui ():

    sg.theme("Default1")
   
    column2 = [[sg.Text('Pinpoint Solve Parameters', background_color='#F7F3EC', justification='center', size=(30, 1))],
               
            [sg.T("Maximum Mag.:"), sg.InputText('13', size=(5,5)), sg.T("Sigma Above Mean:"), sg.InputText('2.5', size=(5,5))],      
            [sg.T("Max Solve Time (sec):"), sg.InputText('60', size=(5,5)), sg.T("Max Match Residual:"), sg.InputText('1.5', size=(5,5))],
            [sg.T("Catalog Expansion:"), sg.InputText('0.8', size=(5,5)), sg.T("Catalog:"), sg.InputText('UCAC3', size=(7,5))]
            ]   
    
    column1 = [[sg.Text('Analysis Parameters', background_color='#F7F3EC', justification='center', size=(30, 1))],
               
            [sg.Frame('Source Capture Mode',[[
              sg.T(""), sg.Radio('Track Rate Mode', "RADIO3", default=False, key="-IN92-"),
               sg.T("            "), sg.Radio('Star Stare Mode', "RADIO3", default=True)]], size=(200,100))],
                  
              [sg.Frame('Image Source',[[      
              sg.T(""), sg.Radio('Ground Based', "RADIO1", default=False, key="-IN82-"),
              sg.T("                "), sg.Radio('Space Based', "RADIO1", default=True)]])], 
              ] 
    
    
    tab1_layout = [
    [sg.T("AstroSolver Processor")],[sg.T("Version 0.1")], [sg.T("   ")],
              [sg.Text("Image Folder:      "), 
               sg.Input(key="-IN2-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN1-"), sg.Text('')],
              [sg.Text("Catalog Folder:    "), 
               sg.Input(key="-IN3-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN4-")],
              [sg.Text("Reference Stars:  "), 
               sg.Input(key="-IN5-" ,change_submits=True), 
               sg.FileBrowse(key="-IN6-")],
              
              
              [sg.T(""), sg.Checkbox('Save Data to Folder', default=True, key="-IN7-")],
              [sg.T(""), sg.Checkbox('Pinpoint Solve', default=True, key="-IN100-")],
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
              [sg.Text("Bias Images:    "), 
               sg.Input(key="-IN20-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN12-"), sg.Text('')],
              [sg.Text("Dark Images:    "), 
               sg.Input(key="-IN30-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN40-")],
              [sg.Text("Flat Images:     "), 
               sg.Input(key="-IN50-" ,change_submits=True), 
               sg.FileBrowse(key="-IN60-")],
              [sg.T("   ")],
              [sg.T("Clean Images Folder:   ")],
              [sg.T("   ")],
              [sg.T(" "), sg.Button("Reduce"), sg.Cancel()]]
    
    column2 = [[sg.Text('Pinpoint Solve Parameters', background_color='#F7F3EC', justification='center', size=(30, 1))],
               
            [sg.T("Maximum Mag.:"), sg.InputText('13', size=(5,5)), sg.T("Sigma Above Mean:"), sg.InputText('10', size=(5,5))],      
            [sg.T("Max Solve Time:"), sg.InputText('10', size=(5,5)), sg.T("Max Match Residual:"), sg.InputText('10', size=(5,5))],
            [sg.T("Catalog Expansion:"), sg.InputText('10', size=(5,5)), sg.T("Max Match Residual:"), sg.InputText('10', size=(5,5))]
            ]   
    
    column1 = [[sg.Text('Pinpoint Solve Parameters', background_color='#F7F3EC', justification='center', size=(30, 1))],
               
            [sg.Frame('Source Capture Mode',[[
              sg.T(""), sg.Radio('Track Rate Mode', "RADIO3", default=False, key="-IN9-"),
               sg.T("            "), sg.Radio('Star Stare Mode', "RADIO3", default=True)]])],
                  
              [sg.Frame('Image Source',[[      
              sg.T(""), sg.Radio('Ground Based', "RADIO1", default=False, key="-IN8-"),
              sg.T("                "), sg.Radio('Space Based', "RADIO1", default=True)]])], 
              ] 
    
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
    window = sg.Window('AstroSolver', layout)
        
    while True:
        event, values = window.read()
        #print(values["-IN2-"])
        if event == sg.WIN_CLOSED or event=="Exit":
            window.close()
            break
        
        elif event == "Submit":
            
            imagefolder =values["-IN2-"]
            catalogfolder = values["-IN3-"]
            refdoc = values["-IN5-"]
            #print(values["-IN2-"])
            window.close()
            return imagefolder, catalogfolder, refdoc
    window.close()           
            
Gui()
#imagefolder, catalogfolder, refdoc = Gui()
#print(imagefolder, catalogfolder, refdoc)