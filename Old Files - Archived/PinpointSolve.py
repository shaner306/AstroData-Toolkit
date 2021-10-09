# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 14:14:59 2021

@author: shane
"""

import pandas as pd
import win32com.client as win32
import win32com
import os
import pywin32_system32
import math
import numpy as np
inbox = 'D:\\Wawrow\\2. Observational Data\\2021-03-10 - Calibrated\\HIP 46066\\LIGHT\\B'
inbox2 = 'D:\\trm-stars-images\\NEOS_SCI_2021099173159frame.fits'
streak2 = 'D:\\Wawrow\\2. Observational Data\\2021-02-07 - Calibrated\\Intelsat 10-02\\LIGHT\\G\\0013_3x3_-10.00_5.00_G_18-32-45.fits'

streak = 'D:\\Breeze-M_R_B_38746U\\CAN_OTT.00018670.BREEZE-M_R_B_#38746U.FIT'
f =win32com.client.Dispatch("Pinpoint.plate")
catloc = 'D:\squid\\USNOA20-All';
refstars_doc = 'D:\\Reference_stars.xlsx'
refstars = pd.read_excel(refstars_doc)
refstars.head()
HIP= refstars["HIP"]
erad = refstars["erad"]
edec= refstars["edec"]
vref= refstars["V"]
bvindex=refstars["(B-V)"]
vrindex=refstars["(V-R)"]
refstarsfin= np.column_stack((HIP, erad,edec,vref))
#print(refstarsfin)



         
     


def getFileList(directory = os.path.dirname(inbox)):
    list = os.listdir(inbox)
    listSize = len(list)
    return [list, listSize]

print(getFileList())
fpath1 = getFileList();
c=fpath1[0]
o=0;
filepathall=[];
for i in c:
    filepath2 = inbox+"\\"+c[o]
    filepathall.append(filepath2)
    o=o+1;

o=0;
f.DetachFITS
for i in filepathall:
    f = win32com.client.Dispatch("Pinpoint.plate")
   #try:
    f.AttachFITS(streak)
    #print(c[o])
    f.TargetDeclination = -6.66777925531915
    f.TargetRightAscension = 139.239438112044/15
    f.Declination = -6.66777925531915
    f.RightAscension =  139.239438112044/15
    #yBin = 3.0;
    #xBin = 3.0; 
    yBin = 3.0
    xBin =  3.0; 
    f.ArcsecperPixelHoriz  = 4.33562092816E-004*3600*3
    f.ArcsecperPixelVert = 4.33562092816E-004*3600*3
    f.UseFaintStars=1
    f.Catalog = 5;
    f.CatalogPath = catloc;
    f.CatalogMaximumMagnitude = 12;
    f.CatalogExpansion = 0.6
    f.SigmaAboveMean = 1.25; 
    f.MaxSolveTime = 60; 
    f.MaxMatchResidual = 8; 
    f.CacheImageStars = 1;
    flag = 0;
    f.FindCatalogStars()
    f.Solve()
    f.MatchedStars.count
    f.FindImageStars()
    print(str(f.RightAscension) +" "+ str(f.Declination))
    
    #except:None
       