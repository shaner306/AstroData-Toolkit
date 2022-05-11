# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 11:22:19 2021

@author: shane
"""

import pandas as pd
import win32com.client as win32
import win32com
import os
import pywin32_system32
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime

#inbox = 'D:\\Wawrow\\2. Observational Data\\2021-03-10 - Calibrated\\HIP 46066\\LIGHT\\B'
inbox = 'D:\\2021-03-10 - Calibrated\\Intelsat 10-02 Post Eclipse\\LIGHT\\G'

f = win32com.client.Dispatch("Pinpoint.plate")
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
p=1;
f.DetachFITS
vintelsat=[]
vobsdates=[]
for i in filepathall:
    f = win32com.client.Dispatch("Pinpoint.plate")
    try:
        f.AttachFITS(filepathall[o])
        #print(c[o])
        dateofimage=f.ExposureStartTime;
        dates = matplotlib.dates.date2num(dateofimage)
        
        #print(newdate)
        f.Declination = f.targetDeclination;
        f.RightAscension = f.targetRightAscension; 
        yBin = 4.33562092816E-004*3600;
        xBin =  4.33131246330E-004*3600; 
        f.ArcsecperPixelHoriz  = xBin;
        f.ArcsecperPixelVert = yBin;
         
        f.Catalog = 5;
        f.CatalogPath = catloc;
        f.CatalogMaximumMagnitude = 13;
        f.CatalogExpansion = 0.8;
        f.SigmaAboveMean = 3.0;
        f.MaxSolveTime = 60; 
        f.MaxMatchResidual = 1.5; 
        flag = 0;
        f.solve()
        nmstars = f.MatchedStars.Count
        mstars = f.MatchedStars;
        print(":"+ str(nmstars))
        print("Reference Stars Located:")
        print("")
        
        MatchedNew=[]
        for j in range(1,nmstars):
            mstar = f.ImageStars.Item(j)
            for i in range(len(erad)):
                if round(mstar.Declination,2)==round(edec[i],2) and round(mstar.RightAscension,2)==round(erad[i],2):
                    MatchedNew.append(HIP[i])
        print(MatchedNew)

        f.DetachFITS()
        f=None
    except:
        print("Pinpoint failed to solve")
