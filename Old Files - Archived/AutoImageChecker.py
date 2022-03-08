import os
from astropy.io import fits, ascii
from dateutil import parser
from datetime import timedelta
from matplotlib import pyplot as plt

def read_fits_file(filepath):
    with fits.open(filepath, memmap=False) as hdul:
        hdr = hdul[0].header
    return hdr

def time_since_Epoch(timeList):
    convertedList=[]
    for i in timeList:  # Convert timestamps of residuals to Time since Epoch
        dt = parser.parse(i)
        timestamp = int(dt.timestamp())
        timestamp = str(float(timestamp))
        convertedList.append(timestamp)
    return convertedList

def create_Plot(x,y,xlabel,ylabel):
    a = plt.scatter(x, y, label='Residuals')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f'{xlabel} Vs {ylabel}')
    plt.legend()
    plt.savefig(f'{xlabel} Vs {ylabel}.png')
    plt.show()
    return None


imageDirectory = "/Users/home/Downloads/processedLEOcal" #Image Directory
g = open("/Users/home/Downloads/data-product-saral-fixed.txt", "r") # Saral Residuals
f = open("/Users/home/Downloads/data-product-H2C-fixed.txt","r") #H2C Residuals
p = open("/Users/home/Downloads/data-product-S6A-fixed.txt","r") #S6A Residuals

objects,exptime,ela_min,sun_min,crval1,crval2,crota2,julian,date=[],[],[],[],[],[],[],[],[] #Header Data taken from images
imageDatesCode,residualEpoch=[],[] #Timestamps converted to Unix seconds since Epoch

#Matched data between Image headers and Residuals
newobjects,newexptime,newela_min,newsun_min,newcrval1,newcrval2,newcrota2,newjulian,newdate,matchedRA,matchedDec=[],[],[],[],[],[],[],[],[],[],[]

#Data taken from residual document
result = []
result2 = []
result3 = []
residualsRA = []
residualsDec = []
datesres=[]

for dirpath, dirnames, filenames in os.walk(imageDirectory):
    for filename in filenames:
        if filename.endswith('.fits'):
            filepath = os.path.join(dirpath, filename)
            hdr = read_fits_file(filepath)
            objects.append(hdr['OBJECT'])
            exptime.append(hdr['AEXPTIME'])
            ela_min.append(hdr['ELA_MIN'])
            sun_min.append(hdr['SUN_MIN'])
            crval1.append(hdr['CRVAL1'])
            crval2.append(hdr['CRVAL2'])
            crota2.append(hdr['CROTA2'])
            julian.append(hdr['JD-OBS'])
            date.append(hdr['DATE-OBS'])
            #print(f'{filepath} Read')

for x in f:
    result.append(x.split())
for y in g:
    result2.append(y.split())
for z in p:
    result3.append(y.split())

for x in range(0, len(result) - 1):
    if x >= 8:
        try:
            if result[x][7] == "Right":
                residualsRA.append(result[x][9])
                st = result[x][0] + " " + result[x][1] + " " + result[x][2] + " " + result[x][3]
                datesres.append(st)
            else:
                residualsDec.append(result[x][8])
        except:
            break

for x in range(0, len(result2) - 1):
    if x >= 8:
        try:
            if result2[x][7] == "Right":
                residualsRA.append(result2[x][9])
                st=result2[x][0]+" "+result2[x][1]+" "+result2[x][2]+" "+result2[x][3]
                datesres.append(st)
            else:
                residualsDec.append(result2[x][8])
        except:
            break

for x in range(0, len(result3) - 1):
    if x >= 8:
        try:
            if result3[x][7] == "Right":
                residualsRA.append(result3[x][9])
                st=result3[x][0]+" "+result3[x][1]+" "+result3[x][2]+" "+result3[x][3]
                datesres.append(st)
            else:
                residualsDec.append(result3[x][8])
        except:
            break



imageDatesCode=time_since_Epoch(date)
residualEpoch=time_since_Epoch(datesres)

for i in range(len(residualEpoch)):
    for j in range(len(imageDatesCode)):
        if float(residualEpoch[i])<=float(imageDatesCode[j])+5 and float(residualEpoch[i])>=float(imageDatesCode[j])-5: # Checks if the dates are within a 10 second window of each other
            newobjects.append(objects[j])
            newexptime.append(exptime[j])
            newela_min.append(ela_min[j])
            newsun_min.append(sun_min[j])
            newcrval1.append(crval1[j])
            newcrval2.append(crval2[j])
            newcrota2.append(crota2[j])
            newjulian.append(julian[j])
            matchedRA.append(residualsRA[i])
            matchedDec.append(residualsDec[i])
        else:
            continue

RA=[float(i) for i in matchedRA]
Dec=[float(i) for i in matchedDec]

#RA Plots
create_Plot(RA,newexptime,"Residual (RA)","Exposure time (s)")
create_Plot(RA,newcrval1,"Residual (RA)","Ref. Pixel RA")
create_Plot(RA,newcrval2,"Residual (RA)","Ref. Pixel Declination")
create_Plot(RA,newcrota2,"Residual (RA)","Rotation Angle Between Axes")
create_Plot(RA,newela_min,"Residual (RA)","Minimum Earth Limb Angle (deg)")
create_Plot(RA,newsun_min,"Residual (RA)","Minimum Sun Angle (deg)")
#Dec Plots
create_Plot(Dec,newexptime,"Residual (Dec)","Exposure time (s)")
create_Plot(Dec,newcrval1,"Residual (Dec)","Ref. Pixel RA")
create_Plot(Dec,newcrval2,"Residual (Dec)","Ref. Pixel Declination")
create_Plot(Dec,newcrota2,"Residual (Dec)","Rotation Angle Between Axes")
create_Plot(Dec,newela_min,"Residual (Dec)","Minimum Earth Limb Angle (deg)")
create_Plot(Dec,newsun_min,"Residual (Dec)","Minimum Sun Angle (deg)")
#RA vs Dec Bullet Plot
create_Plot(RA,Dec,"Residual (RA)","Residual (Dec)")
