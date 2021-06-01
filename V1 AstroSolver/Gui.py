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

def Gui ():

    sg.theme("Default1")
    
    layout = [[sg.T("AstroSolver Processor V0.1")], [sg.T("   ")],
              [sg.Text("Image Folder: "), 
               sg.Input(key="-IN2-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN1-")],
              [sg.Text("Catalog Folder: "), 
               sg.Input(key="-IN3-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN4-")],
              [sg.Text("Reference Stars: "), 
               sg.Input(key="-IN5-" ,change_submits=True), 
               sg.FileBrowse(key="-IN6-")],
              [sg.Button("Submit")]]
    
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
    window = sg.Window('My File Browser', layout, size=(600,250))
        
    while True:
        event, values = window.read()
        #print(values["-IN2-"])
        if event == sg.WIN_CLOSED or event=="Exit":
            break
        elif event == "Submit":
            
            imagefolder =values["-IN2-"]
            catalogfolder = values["-IN3-"]
            refdoc = values["-IN5-"]
            #print(values["-IN2-"])
            window.close()
            return imagefolder, catalogfolder, refdoc
            
            
        
imagefolder, catalogfolder, refdoc = Gui()
print(imagefolder, catalogfolder, refdoc)