# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 11:32:45 2021

@author: shane
"""
import pandas as pd
import time
import os
import datetime as dt
import numpy as np
import PySimpleGUI as sg

tab1_layout = [
    [sg.T("AstroSolver Processor")],[sg.T("Version 0.1")], [sg.T("   ")],
              [sg.Text("Image Folder:      "), 
               sg.Input(key="-IN2-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN1-"), sg.Text('s')],
              [sg.Text("Catalog Folder:    "), 
               sg.Input(key="-IN3-" ,change_submits=True), 
               sg.FolderBrowse(key="-IN4-")],
              [sg.Text("Reference Stars:  "), 
               sg.Input(key="-IN5-" ,change_submits=True), 
               sg.FileBrowse(key="-IN6-")],
              [sg.T("   ")],
              [sg.Button("Solve"), sg.Cancel()]]


tab2_layout = [
    [sg.T('Tab 2')],
    [sg.InputText('data 2', key='file 2'), sg.FileBrowse('data 2 download')],
               ]


layout = [[sg.TabGroup([[sg.Tab('AstroSolver', tab1_layout),
                   sg.Tab('Image Reduction', tab2_layout)
                   ]])]]

window = sg.Window('Multiple tabs one window').Layout(layout)

def sum_tab1(file_1):
    df=pd.read_excel(file_1)
    df["sum"]=df["column 1"]+df["column 2"]
    folder = os.path.dirname(file_1)
    writer = pd.ExcelWriter(folder + "\\sum_output.xlsx")
    df.to_excel(writer, index=False)
    writer.save()

def multiply_tab2(file_2):
    df=pd.read_excel(file_2)
    df["multiply"]=df["column 1"]*df["column 2"]
    folder = os.path.dirname(file_2)
    writer = pd.ExcelWriter(folder + "\\multiply_output.xlsx")
    df.to_excel(writer, index=False)
    writer.save()

while True:
        button, values = window.Read()
        if button in (None, "Quit", "Cancel"):
            window.Close()
            break
        elif button == 'Run':
            if values["file 1"] != "data 1":
                print("Running first tab")
                sum_tab1(values["file 1"])
                print("Output generation completed. Please Close Window")
     
            if values["file 2"] != "data 2":
                print("Running second tab")
                multiply_tab2(values["file 2"])
                print("Completed. Please Close Window")
            
window.Close()