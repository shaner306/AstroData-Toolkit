# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from scipy.stats import rayleigh
import os
from astropy.io import fits, ascii
from datetime import datetime
import time



f=open("/Users/home/Downloads/data-product-H2C-fixed.txt","r")
g=open("/Users/home/Downloads/data-product-saral-fixed.txt","r")

result=[]
result2=[]
residualsRA=[]
residualsDec=[]


for x in f:
    result.append(x.split())
for y in g:
    result2.append(y.split())

for x in range(0,len(result)-1):
    if x >= 8:
     try:
        if result[x][7]=="Right":
            residualsRA.append(result[x][9])


        else:
            residualsDec.append(result[x][8])
     except:
        break
    
for x in range(0,len(result2)-1):
    if x >= 8:
     try:
        if result2[x][7]=="Right":
            residualsRA.append(result2[x][9])

        else:
            residualsDec.append(result2[x][8])
     except:
        break

RA=np.array(residualsRA).astype(np.float64)
DEC=np.array(residualsDec).astype(np.float64)



rms = np.sqrt(np.mean(RA**2))
mean, std= scipy.stats.norm.fit(RA)
mean2, std2 = scipy.stats.norm.fit(DEC)

#bins=[-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12]
bins=[-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12]
n1, bins1, patches1 = plt.hist(x=RA, bins=bins, range= (-10,10), color='#0504aa',
                            alpha=0.5, rwidth=0.95)
l= plt.plot(bins, 100*scipy.stats.norm.pdf(bins, mean, std), 'r--', linewidth=2)
mu = RA.mean()
median = np.median(RA)
sigma = RA.std()
textstr = '\n'.join((
    r'$\mu=%.2f$ arcsec' % (mu, ),
    r'$\mathrm{median}=%.2f$ arcsec' % (median, ),
    r'$\sigma=%.2f$ arcsec' % (sigma, )))
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# place a text box in upper left in axes coords
plt.text(-12.6, 27.5, textstr, fontsize=14,
        verticalalignment='top', bbox=props)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Residuals (arcsec)')
plt.ylabel('Frequency')
plt.xticks(np.arange(-12, 14, 2.0))
plt.title('Right Ascension Residual Histogram')
#plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n1.max()
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
plt.savefig("Calibration Residual Frequency - Right Ascension.png")
plt.show()


n2, bins2, patches2 = plt.hist(x=DEC, bins=bins, range= (-10,10), color='#0504aa',
                            alpha=0.5, rwidth=0.95)
l = plt.plot(bins, 100*scipy.stats.norm.pdf(bins, mean2, std2), 'r--', linewidth=2)
mu = DEC.mean()
median = np.median(DEC)
sigma = DEC.std()
rms = np.sqrt(np.mean(DEC**2))
textstr = '\n'.join((
    r'$\mu=%.2f$ arcsec' % (mu, ),
    r'$\mathrm{median}=%.2f$ arcsec' % (median, ),
    r'$\sigma=%.2f$ arcsec' % (sigma, )
    #r'$\rms=%.2f$ arcsec' % (rms, )
))
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# place a text box in upper left in axes coords
plt.text(-12.6, 27.5, textstr, fontsize=14,
        verticalalignment='top', bbox=props)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Residuals (arcsec)')
plt.ylabel('Frequency')
plt.xticks(np.arange(-12, 14, 2.0))
plt.title('Declination Residual Histogram')
#plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n2.max()
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
plt.savefig("Calibration Residual Frequency - Declination.png")
plt.show()

n3, bins3, patches3 = plt.hist(x=Datapoints, bins=bins, range= (-10,10), color='#0504aa',
                            alpha=0.5, rwidth=0.95)
l= plt.plot(bins3, 100*scipy.stats.norm.pdf(bins3, mean3, std3), 'r--', linewidth=2)
mu = Datapoints.mean()
median = np.median(Datapoints)
sigma = Datapoints.std()
textstr = '\n'.join((
    r'$\mu=%.2f$ arcsec' % (mu, ),
    r'$\mathrm{median}=%.2f$ arcsec' % (median, ),
    r'$\sigma=%.2f$ arcsec' % (sigma, )))
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# place a text box in upper left in axes coords
plt.text(-12.6, 27.5, textstr, fontsize=14,
        verticalalignment='top', bbox=props)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Residuals (arcsec)')
plt.ylabel('Frequency')
plt.xticks(np.arange(-12, 14, 2.0))
plt.title('Combined Residuals (arcsec)')
#plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n3.max()
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
plt.savefig("Calibration Residual Frequency - Combined.png")
plt.show()


f.close()
