# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 10:33:21 2021
@author: shane

-----
INSTRUCTIONS
1. MUST RUN ON 32 BIT Python - Pinpoint will not run on 64 bit python code
2. Reference Data
3. Output Transforms + Standard Magnitudes
    - Errors
4. Output Star Data
    - FWHM
    - Magnitudes
    - Instrumental Mag
    - Flux, Stellar Flux, Visual Magnitude, gaussian data, sigmas
5. TRM mode
    - Martin Levesque Method
    - Brad Wallace Method
    - Photutils combination
    - SExtractor (background Extraction method)- Currently using
5a. Background Extraction
    -Sextractor
    -Polynomial Form Fit
    -Mean Background
    *Filter and Box Size TBD

6. Creating Stars file and outputing solutions
7. Creating Light Curves

"""
import pandas as pd
import win32com.client as win32
import win32com
import os
import pywin32_system32
import math
import numpy as np
inbox = 'D:\\Wawrow\\2. Observational Data\\2021-03-10 - Calibrated\\HIP 46066\\LIGHT\\B'
#inbox = 'D:\\2021-03-10 - Calibrated\\Intelsat 10-02 Post Eclipse\\LIGHT'

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
f.DetachFITS
for i in filepathall:
    f = win32com.client.Dispatch("Pinpoint.plate")
    try:
        f.AttachFITS(filepathall[o])
        #print(c[o])
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
        f.Solve()
        f.MatchedStars.count
        f.FindImageStars()
        #print(f.ImageStars)
        o=o+1;
        s=1
        q=0
        refx=[]
        refy=[]
        vref2=[]
        HIP2=[]
        vrindexdet=[]
        bvindexdet=[]

        for s in range(89):
            
            try:
                f.SkyToXy(erad[s],edec[s]);
                refSTAR_X  = f.ScratchX   #the x result of skytoxy
                refSTAR_Y = f.ScratchY
                if refSTAR_X>0 and refSTAR_X<1900:
                    
                    refx.append(refSTAR_X)
                    refy.append(refSTAR_Y)
                    vref2.append(vref[s])
                    HIP2.append(HIP[s])
                    vrindexdet.append(vrindex[s])
                    bvindexdet.append(bvindex[s])
                # else:
                #   print("Star Found outside bounds of image")
                #s=s+1
            except:
                if s > 89:
                    print("Stop")
                # else:
                #     print('Pinpoint can''t process coords')
                                                          
        nmstars = f.ImageStars.Count
        mstars = f.ImageStars;
        print("Matched Stars:"+ str(nmstars))
        print("Reference Stars Located:")
        print("")
        for i in range(1,nmstars):
            #print(i)
            mstar = f.ImageStars.Item(i)
            X_min = mstar.X - 0.5*mstar.Width
            X_max = mstar.X + 0.5*mstar.Width
            Y_min = mstar.Y - 0.5*mstar.Height
            Y_max = mstar.Y + 0.5*mstar.Height
            #print(i)
            length = len(refx)
            
            exptime = f.ExposureInterval  
            rawflux = mstar.RawFlux;            
            Zp = f.MagZeroPoint;
            
            vmag= Zp - 2.5*(math.log10(rawflux/exptime))
            #print(vmag)

            for j in range(length):
                #print("ref" +str(refx[j]))
                #print(X_max)
                if (refx[j] > X_min) and (refx[j] < X_max):
                    if (refy[j] > Y_min) and (refy[j] < Y_max):
                          if abs(vmag - vref2[j]) < 0.5:
                                  
                                  print("HIP: " +str(HIP2[j]))
                                  print("Located at: X: " +str(mstar.X) + " Y: " + str(mstar.Y))
                                  #print("matched X:" + str(X_max))
                                  #print(str(vref2[j]))
                                  #print(mstar.ColorMagnitude)
                                  print("Reference Mag: "+str(vref2[j]) + " vs " + "Detected Mag: " + str(vmag))
                                  print("")
                                  Bvtransform=(vref2[j]-vmag)/bvindexdet[j]
                                  print("B-V Transform: " + str(Bvtransform))
                                  Vrtransform=(vref2[j]-vmag)/vrindexdet[j]
                                  print("V-R Transform: " + str(Vrtransform))
                 
            #StarXY = [mstar.X mstar.Y]
            #InstrumentalMag= -2.5*log10(mid_pix_valPP*1.00857579708099)
            #ppbgsigma = f.ImageBackgroundSigma;
            #ppbgmean = f.ImageBackgroundmean;
        
        f.DetachFITS()
        f=None
    except:
        None;
n=0;
flag = 1;
