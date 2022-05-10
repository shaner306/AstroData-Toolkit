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
        f.FindImageStars; 
        f.FindCatalogStars; 
        f.MaxSolveTime = 60; 
        f.MaxMatchResidual = 1.5; 
        flag = 0;
        f.FindCatalogStars()
        #f.MatchedStars.count
        f.FindImageStars()
        refx=[]
        refy=[]
        vref2=[]
        HIP2=[]
        vrindexdet=[]
        bvindexdet=[]
        vmagsats=[]
        
        nmstars = f.ImageStars.Count
        mstars = f.ImageStars;
        print(":"+ str(nmstars))
        print("Reference Stars Located:")
        print("")
        
        
        for j in range(1,nmstars):
            mstar = f.ImageStars.Item(j)
            exptime = f.ExposureInterval  
            rawflux = mstar.RawFlux;            
            Zp = 22
              
            vmag= Zp - 2.5*(math.log10(rawflux/exptime))
            instrsmag= -2.5*(math.log10(rawflux/exptime))
            vmagsats.append(instrsmag)
            #print(mstar.X)
            if mstar.X>950 and mstar.X<980 and mstar.Y>580 and mstar.Y<600:
                intelsatmag=instrsmag
                #print(mstar.X)
                print(intelsatmag)
            else:
                None; 
                #print("no")
            #print(mstar.X)
            #print(mstar.Y)
            #print(vmag2)
            #print(Zp)
            #print("")
            
            
        if p<66:
            #print(p)
            #print(intelsatmag)
            if intelsatmag>-100:
                #print("yes")
                vintelsat.append(intelsatmag) 
                vobsdates.append(dates)
                #print("yes")
            else:
                    None;
        else:
        
            if min(vmagsats)<0 or min(vmagsats)>0:
                
                vintelsat.append(min(vmagsats)) 
                vobsdates.append(dates)
                print("non")
            else:
                None;

        f.DetachFITS()
        p=p+1
        f=None
    except:
        None;
plot=plt.scatter(vobsdates,vintelsat)
ax = plot.axes
ax.invert_yaxis()
n=0;
flag = 1;