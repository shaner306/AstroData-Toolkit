# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

f=open("/Users/home/Downloads/data-product-saral-2.txt","r")

result=[]
residualsRA=[]
residualsDec=[]
for x in f:
    result.append(x.split())

for x in range(0,len(result)-1):
    if x >= 8:
     try:
        if result[x][7]=="Right":
            residualsRA.append(result[x][9])
        else:
            residualsDec.append(result[x][8])
     except:
        break;


RA=np.array(residualsRA).astype(np.float64)
DEC=np.array(residualsDec).astype(np.float64)

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
    r'$\mu=%.2f$' % (mu, ),
    r'$\mathrm{median}=%.2f$' % (median, ),
    r'$\sigma=%.2f$' % (sigma, )))
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# place a text box in upper left in axes coords
plt.text(-10, 18.0, textstr, fontsize=14,
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

textstr = '\n'.join((
    r'$\mu=%.2f$' % (mu, ),
    r'$\mathrm{median}=%.2f$' % (median, ),
    r'$\sigma=%.2f$' % (sigma, )))
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# place a text box in upper left in axes coords
plt.text(-10, 18.0, textstr, fontsize=14,
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



f.close()