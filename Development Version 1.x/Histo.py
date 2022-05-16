# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import matplotlib.pyplot as plt
import numpy as np

f=open("/Users/home/Downloads/Data Product.txt","r")

result=[]
residualsRA=[]
residualsDec=[]
for x in f:
    result.append(x.split())

for x in range(0,len(result)):
    if x >= 7:
        if result[x][7]=="Right":
            residualsRA.append(result[x][9])
        else:
            residualsDec.append(result[x][8])



# An "interface" to matplotlib.axes.Axes.hist() method
n, bins, patches = plt.hist(x=data_new, bins=[-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12], range= (-10,10), color='#0504aa',
                            alpha=0.5, rwidth=0.95)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Residuals (arcsec)')
plt.ylabel('Frequency')
plt.xticks(np.arange(-12, 14, 2.0))
plt.title('Pre-Fit Residual Histogram')
plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

f.close()